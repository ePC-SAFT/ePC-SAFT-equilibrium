#include "held2.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>

namespace epcsaft_equilibrium {
namespace {

constexpr double kCandidateRadius = 1.0e-3;
constexpr double kVolumeLower = 0.5;
constexpr double kVolumeUpper = 1.5;
constexpr double kNumericalTolerance = 1.0e-8;
constexpr double kKktTolerance = 1.0e-7;
constexpr double kMergeTolerance = 1.0e-6;
constexpr double kDistinctPhaseTolerance = 1.0e-3;
constexpr double kActivePhaseTolerance = 1.0e-10;

double maximum_abs(const std::vector<double>& values) {
    double result = 0.0;
    for (double value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

double maximum_abs_difference(
    const std::vector<double>& left,
    const std::vector<double>& right
) {
    if (left.size() != right.size()) {
        throw std::invalid_argument("HELD2 Stage III vectors have different sizes");
    }
    double result = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        result = std::max(result, std::abs(left[index] - right[index]));
    }
    return result;
}

double charge_residual(
    const std::vector<double>& charges,
    const std::vector<double>& fractions
) {
    double result = 0.0;
    for (std::size_t index = 0; index < charges.size(); ++index) {
        result += charges[index] * fractions[index];
    }
    return result;
}

std::string ipopt_status_name(Ipopt::ApplicationReturnStatus status) {
    switch (status) {
        case Ipopt::Solve_Succeeded:
            return "solve_succeeded";
        case Ipopt::Solved_To_Acceptable_Level:
            return "solved_to_acceptable_level";
        case Ipopt::Infeasible_Problem_Detected:
            return "infeasible_problem_detected";
        case Ipopt::Search_Direction_Becomes_Too_Small:
            return "search_direction_too_small";
        case Ipopt::Diverging_Iterates:
            return "diverging_iterates";
        case Ipopt::Maximum_Iterations_Exceeded:
            return "maximum_iterations_exceeded";
        case Ipopt::Restoration_Failed:
            return "restoration_failed";
        case Ipopt::Error_In_Step_Computation:
            return "error_in_step_computation";
        case Ipopt::Invalid_Number_Detected:
            return "invalid_number_detected";
        case Ipopt::NonIpopt_Exception_Thrown:
            return "non_ipopt_exception_thrown";
        default:
            return "ipopt_status_" + std::to_string(static_cast<int>(status));
    }
}

std::size_t independent_retained_index(const Held2Coordinates& coordinates) {
    if (coordinates.independent_indices.size() != 1) {
        throw std::invalid_argument(
            "manufactured HELD2 Stage III requires one independent modified composition"
        );
    }
    return static_cast<std::size_t>(
        std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            coordinates.independent_indices.front()
        ) - coordinates.retained_indices.begin()
    );
}

Held2StageIIINlpEvaluation evaluate_stage_iii(
    const Held2Coordinates& coordinates,
    double feed,
    std::size_t phase_count,
    const std::vector<double>& variables,
    const std::vector<double>& equality_multipliers,
    double objective_factor
) {
    constexpr std::size_t block_size = 3;
    const std::size_t variable_count = phase_count * block_size;
    constexpr std::size_t constraint_count = 2;
    if (phase_count < 2 || variables.size() != variable_count
        || equality_multipliers.size() != constraint_count) {
        throw std::invalid_argument("HELD2 Stage III dimensions do not match the candidate set");
    }
    Held2StageIIINlpEvaluation result;
    result.objective_gradient.assign(variable_count, 0.0);
    result.constraints.assign(constraint_count, 0.0);
    result.constraint_jacobian.assign(constraint_count * variable_count, 0.0);
    result.lagrangian_hessian.assign(variable_count * variable_count, 0.0);
    result.constraints[0] = -1.0;
    result.constraints[1] = -feed;

    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const std::size_t offset = phase * block_size;
        const double fraction = variables[offset];
        const double composition = variables[offset + 1];
        const double volume = variables[offset + 2];
        if (!std::isfinite(fraction) || !std::isfinite(composition)
            || !std::isfinite(volume) || volume <= 0.0) {
            throw std::invalid_argument("HELD2 Stage III variables must be finite and physical");
        }
        const Held2StateEvaluation state = evaluate_held2_manufactured_state(
            coordinates,
            {composition},
            std::log(volume)
        );
        const double composition_gradient = state.gradient[0];
        const double volume_gradient = state.gradient[1] / volume;
        const double composition_hessian = state.hessian[0];
        const double composition_volume_hessian = state.hessian[1] / volume;
        const double volume_hessian =
            (state.hessian[3] - state.gradient[1]) / (volume * volume);

        result.objective += fraction * state.objective;
        result.objective_gradient[offset] = state.objective;
        result.objective_gradient[offset + 1] = fraction * composition_gradient;
        result.objective_gradient[offset + 2] = fraction * volume_gradient;
        result.constraints[0] += fraction;
        result.constraints[1] += fraction * composition;
        result.constraint_jacobian[offset] = 1.0;
        result.constraint_jacobian[variable_count + offset] = composition;
        result.constraint_jacobian[variable_count + offset + 1] = fraction;

        auto add_symmetric = [&result, variable_count](
                                 std::size_t row,
                                 std::size_t column,
                                 double value
                             ) {
            result.lagrangian_hessian[row * variable_count + column] += value;
            if (row != column) {
                result.lagrangian_hessian[column * variable_count + row] += value;
            }
        };
        add_symmetric(
            offset + 1,
            offset,
            objective_factor * composition_gradient + equality_multipliers[1]
        );
        add_symmetric(offset + 2, offset, objective_factor * volume_gradient);
        add_symmetric(
            offset + 1,
            offset + 1,
            objective_factor * fraction * composition_hessian
        );
        add_symmetric(
            offset + 2,
            offset + 1,
            objective_factor * fraction * composition_volume_hessian
        );
        add_symmetric(
            offset + 2,
            offset + 2,
            objective_factor * fraction * volume_hessian
        );
    }
    result.lagrangian_gradient = result.objective_gradient;
    for (std::size_t variable = 0; variable < variable_count; ++variable) {
        result.lagrangian_gradient[variable] *= objective_factor;
        for (std::size_t constraint = 0; constraint < constraint_count; ++constraint) {
            result.lagrangian_gradient[variable] += equality_multipliers[constraint]
                * result.constraint_jacobian[constraint * variable_count + variable];
        }
    }
    return result;
}

struct StageIIIRun {
    bool solver_converged = false;
    std::string solver_status;
    std::string callback_error;
    std::vector<double> variables;
    std::vector<double> equality_multipliers;
    std::vector<double> lower_bound_multipliers;
    std::vector<double> upper_bound_multipliers;
};

class Held2StageIIITnlp final : public Ipopt::TNLP {
public:
    Held2StageIIITnlp(
        Held2Coordinates coordinates,
        double feed,
        std::vector<std::array<double, 2>> candidates,
        std::vector<double> initial
    )
        : coordinates_(std::move(coordinates)),
          feed_(feed),
          candidates_(std::move(candidates)),
          initial_(std::move(initial)) {}

    bool get_nlp_info(
        Ipopt::Index& n,
        Ipopt::Index& m,
        Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag,
        IndexStyleEnum& index_style
    ) override {
        n = static_cast<Ipopt::Index>(initial_.size());
        m = 2;
        nnz_jac_g = static_cast<Ipopt::Index>(3 * candidates_.size());
        nnz_h_lag = n * (n + 1) / 2;
        index_style = TNLP::C_STYLE;
        return true;
    }

    bool get_bounds_info(
        Ipopt::Index n,
        Ipopt::Number* x_l,
        Ipopt::Number* x_u,
        Ipopt::Index m,
        Ipopt::Number* g_l,
        Ipopt::Number* g_u
    ) override {
        if (n != static_cast<Ipopt::Index>(initial_.size()) || m != 2) {
            return false;
        }
        for (std::size_t phase = 0; phase < candidates_.size(); ++phase) {
            const std::size_t offset = 3 * phase;
            x_l[offset] = 0.0;
            x_u[offset] = 1.0;
            x_l[offset + 1] = std::max(
                coordinates_.independent_lower_bounds.front(),
                candidates_[phase][0] - kCandidateRadius
            );
            x_u[offset + 1] = std::min(
                coordinates_.independent_upper_bounds.front(),
                candidates_[phase][0] + kCandidateRadius
            );
            x_l[offset + 2] = kVolumeLower;
            x_u[offset + 2] = kVolumeUpper;
        }
        g_l[0] = 0.0;
        g_u[0] = 0.0;
        g_l[1] = 0.0;
        g_u[1] = 0.0;
        return true;
    }

    bool get_starting_point(
        Ipopt::Index n,
        bool init_x,
        Ipopt::Number* x,
        bool init_z,
        Ipopt::Number*,
        Ipopt::Number*,
        Ipopt::Index m,
        bool init_lambda,
        Ipopt::Number*
    ) override {
        if (n != static_cast<Ipopt::Index>(initial_.size()) || m != 2 || !init_x
            || init_z || init_lambda) {
            return false;
        }
        std::copy(initial_.begin(), initial_.end(), x);
        return true;
    }

    bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool, Ipopt::Number& value) override {
        try {
            value = evaluate(n, x, {0.0, 0.0}, 1.0).objective;
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
    }

    bool eval_grad_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number* gradient
    ) override {
        try {
            const auto result = evaluate(n, x, {0.0, 0.0}, 1.0);
            std::copy(result.objective_gradient.begin(), result.objective_gradient.end(), gradient);
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
    }

    bool eval_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Index m,
        Ipopt::Number* constraints
    ) override {
        try {
            if (m != 2) {
                return false;
            }
            const auto result = evaluate(n, x, {0.0, 0.0}, 1.0);
            std::copy(result.constraints.begin(), result.constraints.end(), constraints);
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
    }

    bool eval_jac_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Index m,
        Ipopt::Index nonzero_count,
        Ipopt::Index* rows,
        Ipopt::Index* columns,
        Ipopt::Number* values
    ) override {
        if (m != 2 || nonzero_count != static_cast<Ipopt::Index>(3 * candidates_.size())) {
            return false;
        }
        Ipopt::Index position = 0;
        if (values == nullptr) {
            for (std::size_t phase = 0; phase < candidates_.size(); ++phase) {
                rows[position] = 0;
                columns[position++] = static_cast<Ipopt::Index>(3 * phase);
            }
            for (std::size_t phase = 0; phase < candidates_.size(); ++phase) {
                rows[position] = 1;
                columns[position++] = static_cast<Ipopt::Index>(3 * phase);
                rows[position] = 1;
                columns[position++] = static_cast<Ipopt::Index>(3 * phase + 1);
            }
            return true;
        }
        try {
            const auto result = evaluate(n, x, {0.0, 0.0}, 1.0);
            for (std::size_t phase = 0; phase < candidates_.size(); ++phase) {
                values[position++] = result.constraint_jacobian[3 * phase];
            }
            for (std::size_t phase = 0; phase < candidates_.size(); ++phase) {
                values[position++] = result.constraint_jacobian[n + 3 * phase];
                values[position++] = result.constraint_jacobian[n + 3 * phase + 1];
            }
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
    }

    bool eval_h(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number objective_factor,
        Ipopt::Index m,
        const Ipopt::Number* lambda,
        bool,
        Ipopt::Index nonzero_count,
        Ipopt::Index* rows,
        Ipopt::Index* columns,
        Ipopt::Number* values
    ) override {
        if (m != 2 || nonzero_count != n * (n + 1) / 2) {
            return false;
        }
        Ipopt::Index position = 0;
        if (values == nullptr) {
            for (Ipopt::Index row = 0; row < n; ++row) {
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    rows[position] = row;
                    columns[position++] = column;
                }
            }
            return true;
        }
        try {
            if (lambda == nullptr) {
                return false;
            }
            const std::vector<double> multipliers(lambda, lambda + m);
            const auto result = evaluate(n, x, multipliers, objective_factor);
            for (Ipopt::Index row = 0; row < n; ++row) {
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    values[position++] = result.lagrangian_hessian[
                        static_cast<std::size_t>(row * n + column)
                    ];
                }
            }
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
    }

    void finalize_solution(
        Ipopt::SolverReturn status,
        Ipopt::Index n,
        const Ipopt::Number* x,
        const Ipopt::Number* z_l,
        const Ipopt::Number* z_u,
        Ipopt::Index m,
        const Ipopt::Number*,
        const Ipopt::Number* lambda,
        Ipopt::Number,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*
    ) override {
        solver_converged_ = status == Ipopt::SUCCESS
            || status == Ipopt::STOP_AT_ACCEPTABLE_POINT;
        if (n == static_cast<Ipopt::Index>(initial_.size()) && m == 2 && x != nullptr
            && z_l != nullptr && z_u != nullptr && lambda != nullptr) {
            variables_.assign(x, x + n);
            lower_bound_multipliers_.assign(z_l, z_l + n);
            upper_bound_multipliers_.assign(z_u, z_u + n);
            equality_multipliers_.assign(lambda, lambda + m);
        }
    }

    [[nodiscard]] StageIIIRun result(const std::string& solver_status) const {
        return {
            solver_converged_,
            solver_status,
            callback_error_,
            variables_,
            equality_multipliers_,
            lower_bound_multipliers_,
            upper_bound_multipliers_,
        };
    }

private:
    [[nodiscard]] Held2StageIIINlpEvaluation evaluate(
        Ipopt::Index n,
        const Ipopt::Number* x,
        const std::vector<double>& multipliers,
        double objective_factor
    ) const {
        if (n != static_cast<Ipopt::Index>(initial_.size())) {
            throw std::invalid_argument("HELD2 Stage III variable count changed");
        }
        return evaluate_stage_iii(
            coordinates_,
            feed_,
            candidates_.size(),
            std::vector<double>(x, x + n),
            multipliers,
            objective_factor
        );
    }

    Held2Coordinates coordinates_;
    double feed_;
    std::vector<std::array<double, 2>> candidates_;
    std::vector<double> initial_;
    bool solver_converged_ = false;
    std::string callback_error_;
    std::vector<double> variables_;
    std::vector<double> equality_multipliers_;
    std::vector<double> lower_bound_multipliers_;
    std::vector<double> upper_bound_multipliers_;
};

StageIIIRun run_stage_iii(
    const Held2Coordinates& coordinates,
    double feed,
    const std::vector<std::array<double, 2>>& candidates,
    const std::vector<double>& initial
) {
    auto* raw_problem = new Held2StageIIITnlp(coordinates, feed, candidates, initial);
    Ipopt::SmartPtr<Ipopt::TNLP> problem = raw_problem;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application = IpoptApplicationFactory();
    application->Options()->SetStringValue("option_file_name", "");
    application->Options()->SetIntegerValue("print_level", 0);
    application->Options()->SetStringValue("sb", "yes");
    application->Options()->SetIntegerValue("max_iter", 300);
    application->Options()->SetNumericValue("tol", 1.0e-10);
    application->Options()->SetNumericValue("acceptable_tol", 1.0e-9);
    application->Options()->SetIntegerValue("acceptable_iter", 0);
    application->Options()->SetNumericValue("constr_viol_tol", 1.0e-10);
    application->Options()->SetStringValue("jacobian_approximation", "exact");
    application->Options()->SetStringValue("hessian_approximation", "exact");
    application->Options()->SetStringValue("nlp_scaling_method", "none");
    application->Options()->SetNumericValue("bound_relax_factor", 0.0);
    application->Options()->SetStringValue("honor_original_bounds", "yes");
    application->Options()->SetStringValue("check_derivatives_for_naninf", "yes");
    const Ipopt::ApplicationReturnStatus initialized = application->Initialize();
    if (initialized != Ipopt::Solve_Succeeded) {
        return {false, ipopt_status_name(initialized)};
    }
    const Ipopt::ApplicationReturnStatus status = application->OptimizeTNLP(problem);
    return raw_problem->result(ipopt_status_name(status));
}

}  // namespace

Held2StageIIINlpEvaluation evaluate_held2_manufactured_stage_iii_nlp(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::vector<std::array<double, 2>>& candidates,
    const std::vector<double>& variables,
    const std::vector<double>& equality_multipliers
) {
    const Held2Coordinates coordinates = make_held2_coordinates(charges);
    const std::vector<double> modified_feed = held2_transform_physical_fractions(
        coordinates,
        physical_feed
    );
    return evaluate_stage_iii(
        coordinates,
        modified_feed[independent_retained_index(coordinates)],
        candidates.size(),
        variables,
        equality_multipliers,
        1.0
    );
}

Held2StageIIIResult solve_held2_manufactured_stage_iii(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::vector<std::array<double, 2>>& candidates
) {
    Held2StageIIIResult result;
    result.input_candidate_count = static_cast<int>(candidates.size());
    if (candidates.size() < 2) {
        result.failure_reason = "candidate_set_incomplete";
        return result;
    }
    const Held2Coordinates coordinates = make_held2_coordinates(charges);
    const std::vector<double> modified_feed = held2_transform_physical_fractions(
        coordinates,
        physical_feed
    );
    const std::size_t independent_retained = independent_retained_index(coordinates);
    const double feed = modified_feed[independent_retained];
    const auto minimum = std::min_element(
        candidates.begin(),
        candidates.end(),
        [](const auto& left, const auto& right) { return left[0] < right[0]; }
    );
    const auto maximum = std::max_element(
        candidates.begin(),
        candidates.end(),
        [](const auto& left, const auto& right) { return left[0] < right[0]; }
    );
    if (minimum == candidates.end() || maximum == candidates.end()
        || minimum->at(0) >= feed || maximum->at(0) <= feed) {
        result.failure_reason = "candidate_set_does_not_bracket_feed";
        return result;
    }
    const std::size_t left_index = static_cast<std::size_t>(minimum - candidates.begin());
    const std::size_t right_index = static_cast<std::size_t>(maximum - candidates.begin());
    const double trace_fraction = 1.0e-6;
    std::vector<double> initial(3 * candidates.size(), 0.0);
    double reserved_fraction = 0.0;
    double reserved_balance = 0.0;
    for (std::size_t phase = 0; phase < candidates.size(); ++phase) {
        initial[3 * phase + 1] = candidates[phase][0];
        initial[3 * phase + 2] = candidates[phase][1];
        if (phase != left_index && phase != right_index) {
            initial[3 * phase] = trace_fraction;
            reserved_fraction += trace_fraction;
            reserved_balance += trace_fraction * candidates[phase][0];
        }
    }
    const double remaining_fraction = 1.0 - reserved_fraction;
    const double remaining_balance = feed - reserved_balance;
    const double separation = candidates[right_index][0] - candidates[left_index][0];
    const double right_fraction =
        (remaining_balance - remaining_fraction * candidates[left_index][0]) / separation;
    const double left_fraction = remaining_fraction - right_fraction;
    if (left_fraction < 0.0 || right_fraction < 0.0) {
        result.failure_reason = "candidate_set_initialization_infeasible";
        return result;
    }
    initial[3 * left_index] = left_fraction;
    initial[3 * right_index] = right_fraction;

    const StageIIIRun run = run_stage_iii(coordinates, feed, candidates, initial);
    result.solver_status = run.solver_status;
    if (!run.solver_converged || !run.callback_error.empty()
        || run.variables.size() != initial.size()) {
        result.numerical_status = "not_converged";
        result.failure_reason = run.callback_error.empty()
            ? "stage_iii_solver_not_converged"
            : run.callback_error;
        return result;
    }
    const Held2StageIIINlpEvaluation nlp = evaluate_stage_iii(
        coordinates,
        feed,
        candidates.size(),
        run.variables,
        run.equality_multipliers,
        1.0
    );
    std::vector<double> kkt = nlp.lagrangian_gradient;
    if (run.lower_bound_multipliers.size() != kkt.size()
        || run.upper_bound_multipliers.size() != kkt.size()) {
        result.numerical_status = "not_converged";
        result.failure_reason = "stage_iii_multiplier_evidence_missing";
        return result;
    }
    for (std::size_t index = 0; index < kkt.size(); ++index) {
        kkt[index] -= run.lower_bound_multipliers[index];
        kkt[index] += run.upper_bound_multipliers[index];
    }
    result.kkt_stationarity_inf_norm = maximum_abs(kkt);
    if (maximum_abs(nlp.constraints) > kNumericalTolerance
        || result.kkt_stationarity_inf_norm > kKktTolerance) {
        result.numerical_status = "not_converged";
        result.failure_reason = "stage_iii_numerical_certificate_failed";
        return result;
    }
    result.numerical_status = "converged";

    std::vector<Held2StageIIIPhase> refined;
    refined.reserve(candidates.size());
    for (std::size_t phase = 0; phase < candidates.size(); ++phase) {
        const double fraction = run.variables[3 * phase];
        const double composition = run.variables[3 * phase + 1];
        const double volume = run.variables[3 * phase + 2];
        const Held2StateEvaluation state = evaluate_held2_manufactured_state(
            coordinates,
            {composition},
            std::log(volume)
        );
        bool merged = false;
        for (Held2StageIIIPhase& existing : refined) {
            if (maximum_abs_difference(existing.modified_fractions, state.modified_fractions)
                    <= kMergeTolerance
                && std::abs(existing.volume - volume) <= kMergeTolerance) {
                const double combined = existing.phase_fraction + fraction;
                if (combined > 0.0) {
                    for (std::size_t index = 0; index < existing.modified_fractions.size(); ++index) {
                        existing.modified_fractions[index] =
                            (existing.phase_fraction * existing.modified_fractions[index]
                             + fraction * state.modified_fractions[index])
                            / combined;
                    }
                    for (std::size_t index = 0; index < existing.physical_fractions.size(); ++index) {
                        existing.physical_fractions[index] =
                            (existing.phase_fraction * existing.physical_fractions[index]
                             + fraction * state.physical_amounts[index])
                            / combined;
                    }
                    existing.volume =
                        (existing.phase_fraction * existing.volume + fraction * volume) / combined;
                }
                existing.phase_fraction = combined;
                ++result.retired_duplicate_count;
                merged = true;
                break;
            }
        }
        if (!merged) {
            refined.push_back({
                fraction,
                state.modified_fractions,
                state.physical_amounts,
                volume,
            });
        }
    }
    refined.erase(
        std::remove_if(
            refined.begin(),
            refined.end(),
            [](const Held2StageIIIPhase& phase) {
                return phase.phase_fraction <= kActivePhaseTolerance;
            }
        ),
        refined.end()
    );
    if (refined.size() < 2) {
        result.physical_status = "rejected";
        result.failure_reason = "collapsed_phase_set";
        return result;
    }
    result.phases = refined;
    result.objective = nlp.objective;
    std::vector<double> modified_balance(modified_feed.size(), 0.0);
    std::vector<double> ordinary_balance(physical_feed.size(), 0.0);
    result.phase_charge_inf_norm = 0.0;
    for (const Held2StageIIIPhase& phase : refined) {
        for (std::size_t index = 0; index < modified_balance.size(); ++index) {
            modified_balance[index] += phase.phase_fraction * phase.modified_fractions[index];
        }
        for (std::size_t index = 0; index < ordinary_balance.size(); ++index) {
            ordinary_balance[index] += phase.phase_fraction * phase.physical_fractions[index];
        }
        result.phase_charge_inf_norm = std::max(
            result.phase_charge_inf_norm,
            std::abs(charge_residual(charges, phase.physical_fractions))
        );
    }
    result.modified_balance_inf_norm = maximum_abs_difference(
        modified_balance,
        modified_feed
    );
    result.ordinary_balance_inf_norm = maximum_abs_difference(
        ordinary_balance,
        physical_feed
    );
    result.pressure_stationarity_inf_norm = 0.0;
    result.modified_potential_mixed_gap = 0.0;
    result.minimum_phase_distance = std::numeric_limits<double>::infinity();
    std::vector<std::vector<double>> potentials;
    for (const Held2StageIIIPhase& phase : refined) {
        const double composition = phase.modified_fractions[independent_retained];
        const Held2StateEvaluation state = evaluate_held2_manufactured_state(
            coordinates,
            {composition},
            std::log(phase.volume)
        );
        result.pressure_stationarity_inf_norm = std::max(
            result.pressure_stationarity_inf_norm,
            std::abs(state.pressure_stationarity_relative)
        );
        potentials.push_back(state.modified_potentials);
    }
    for (std::size_t left = 0; left < refined.size(); ++left) {
        for (std::size_t right = left + 1; right < refined.size(); ++right) {
            result.minimum_phase_distance = std::min(
                result.minimum_phase_distance,
                maximum_abs_difference(
                    refined[left].modified_fractions,
                    refined[right].modified_fractions
                )
            );
        }
    }
    for (std::size_t retained = 0; retained < modified_feed.size(); ++retained) {
        double lower_bound = 1.0e-10 * coordinates.modified_factors[retained];
        bool trace = false;
        for (const Held2StageIIIPhase& phase : refined) {
            trace = trace || phase.modified_fractions[retained] <= 10.0 * lower_bound;
        }
        if (trace) {
            ++result.trace_component_count;
            continue;
        }
        ++result.certified_modified_potential_count;
        for (std::size_t left = 0; left < potentials.size(); ++left) {
            for (std::size_t right = left + 1; right < potentials.size(); ++right) {
                const double absolute_gap = std::abs(
                    potentials[left][retained] - potentials[right][retained]
                );
                const double mixed_scale = 1.0 + std::max(
                    std::abs(potentials[left][retained]),
                    std::abs(potentials[right][retained])
                );
                result.modified_potential_mixed_gap = std::max(
                    result.modified_potential_mixed_gap,
                    absolute_gap / mixed_scale
                );
            }
        }
    }
    if (result.trace_component_count > 0) {
        result.trace_refinement_status = "complementarity_refinement_required";
        result.physical_status = "not_adjudicated";
        result.failure_reason = "trace_component_requires_log_refinement";
        return result;
    }
    result.trace_refinement_status = "not_required";
    result.enumeration_objective_gap = result.objective
        - held2_manufactured_enumerated_objective(feed);
    const bool physical = result.modified_balance_inf_norm < kNumericalTolerance
        && result.ordinary_balance_inf_norm < kNumericalTolerance
        && result.phase_charge_inf_norm < kNumericalTolerance
        && result.pressure_stationarity_inf_norm < kNumericalTolerance
        && result.modified_potential_mixed_gap < kNumericalTolerance
        && result.minimum_phase_distance > kDistinctPhaseTolerance
        && std::abs(result.enumeration_objective_gap) < kNumericalTolerance;
    if (!physical) {
        result.physical_status = "rejected";
        result.failure_reason = "stage_iii_physical_certificate_failed";
        return result;
    }
    result.physical_status = "accepted";
    result.feedback = "none";
    return result;
}

}  // namespace epcsaft_equilibrium
