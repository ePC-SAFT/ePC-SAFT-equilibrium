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
constexpr double kConstraintLowerInfinity = -2.0e19;

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
    int recoverable_domain_rejection_count = 0;
    std::string last_domain_rejection;
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
            0,
            {},
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
        StageIIIRun failure;
        failure.solver_status = ipopt_status_name(initialized);
        return failure;
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
        double lower_bound =
            kHeld2ModifiedLowerScale * coordinates.modified_factors[retained];
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

namespace {

std::vector<std::size_t> independent_retained_positions(
    const Held2Coordinates& coordinates
) {
    std::vector<std::size_t> positions;
    positions.reserve(coordinates.independent_indices.size());
    for (std::size_t component : coordinates.independent_indices) {
        const auto position = std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            component
        );
        if (position == coordinates.retained_indices.end()) {
            throw std::invalid_argument("HELD2 independent coordinate is not retained");
        }
        positions.push_back(static_cast<std::size_t>(
            position - coordinates.retained_indices.begin()
        ));
    }
    return positions;
}

double dependent_modified_fraction_lower(const Held2Coordinates& coordinates) {
    const auto position = std::find(
        coordinates.retained_indices.begin(),
        coordinates.retained_indices.end(),
        coordinates.dependent_index
    );
    if (position == coordinates.retained_indices.end()) {
        throw std::invalid_argument("HELD2 dependent coordinate is not retained");
    }
    return kHeld2ModifiedLowerScale
        * coordinates.modified_factors[static_cast<std::size_t>(
            position - coordinates.retained_indices.begin()
        )];
}

std::size_t general_stage_iii_constraint_count(
    std::size_t dimension,
    std::size_t phase_count
) {
    return dimension + 1 + phase_count;
}

double general_stage_iii_composition_sum_upper(
    const Held2Coordinates& coordinates
) {
    return 1.0 - dependent_modified_fraction_lower(coordinates);
}

bool general_stage_iii_simplex_domain_valid(
    const Held2Coordinates& coordinates,
    std::size_t phase_count,
    const std::vector<double>& variables
) {
    const std::size_t dimension = coordinates.independent_indices.size();
    const std::size_t block_size = dimension + 2;
    if (variables.size() != phase_count * block_size) {
        return false;
    }
    const double composition_sum_upper =
        general_stage_iii_composition_sum_upper(coordinates);
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const std::size_t offset = phase * block_size;
        double sum = 0.0;
        for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
            sum += variables[offset + 1 + coordinate];
        }
        if (!std::isfinite(sum) || sum > composition_sum_upper) {
            return false;
        }
    }
    return true;
}

Held2StageIIINlpEvaluation evaluate_general_stage_iii_algebraic(
    const std::vector<double>& feed,
    std::size_t phase_count,
    const std::vector<double>& variables
) {
    const std::size_t dimension = feed.size();
    const std::size_t block_size = dimension + 2;
    const std::size_t variable_count = phase_count * block_size;
    const std::size_t constraint_count =
        general_stage_iii_constraint_count(dimension, phase_count);
    if (phase_count < 2 || variables.size() != variable_count) {
        throw std::invalid_argument("HELD2 Stage III dimensions do not match the candidate set");
    }
    Held2StageIIINlpEvaluation result;
    result.objective_gradient.assign(variable_count, 0.0);
    result.constraints.assign(constraint_count, 0.0);
    result.constraint_jacobian.assign(constraint_count * variable_count, 0.0);
    result.lagrangian_hessian.assign(variable_count * variable_count, 0.0);
    result.constraints[0] = -1.0;
    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
        result.constraints[coordinate + 1] = -feed[coordinate];
    }
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const std::size_t offset = phase * block_size;
        const double fraction = variables[offset];
        result.constraints[0] += fraction;
        result.constraint_jacobian[offset] = 1.0;
        double independent_sum = 0.0;
        for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
            const double independent = variables[offset + 1 + coordinate];
            independent_sum += independent;
            result.constraints[coordinate + 1] += fraction * independent;
            result.constraint_jacobian[
                (coordinate + 1) * variable_count + offset
            ] = independent;
            result.constraint_jacobian[
                (coordinate + 1) * variable_count + offset + 1 + coordinate
            ] = fraction;
            result.constraint_jacobian[
                (dimension + 1 + phase) * variable_count + offset + 1 + coordinate
            ] = 1.0;
        }
        result.constraints[dimension + 1 + phase] = independent_sum;
    }
    return result;
}

Held2StageIIINlpEvaluation evaluate_general_stage_iii(
    const Held2Coordinates& coordinates,
    const std::vector<double>& feed,
    const Held2StateEvaluator& evaluator,
    std::size_t phase_count,
    const std::vector<double>& variables,
    const std::vector<double>& equality_multipliers,
    double objective_factor
) {
    const std::size_t dimension = coordinates.independent_indices.size();
    const std::size_t block_size = dimension + 2;
    const std::size_t variable_count = phase_count * block_size;
    const std::size_t constraint_count =
        general_stage_iii_constraint_count(dimension, phase_count);
    if (phase_count < 2 || feed.size() != dimension
        || variables.size() != variable_count
        || equality_multipliers.size() != constraint_count) {
        throw std::invalid_argument("HELD2 Stage III dimensions do not match the candidate set");
    }
    Held2StageIIINlpEvaluation result = evaluate_general_stage_iii_algebraic(
        feed, phase_count, variables
    );
    if (!general_stage_iii_simplex_domain_valid(
            coordinates, phase_count, variables
        )) {
        throw std::domain_error(
            "independent modified fractions do not leave the declared dependent lower bound"
        );
    }
    const auto add_symmetric = [&result, variable_count](
                                   std::size_t row,
                                   std::size_t column,
                                   double value
                               ) {
        result.lagrangian_hessian[row * variable_count + column] += value;
        if (row != column) {
            result.lagrangian_hessian[column * variable_count + row] += value;
        }
    };
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const std::size_t offset = phase * block_size;
        const double fraction = variables[offset];
        const std::vector<double> independent(
            variables.begin() + static_cast<std::ptrdiff_t>(offset + 1),
            variables.begin() + static_cast<std::ptrdiff_t>(offset + 1 + dimension)
        );
        const double phase_coordinate = variables[offset + 1 + dimension];
        const Held2StateEvaluation state = evaluator(independent, phase_coordinate);
        if (state.gradient.size() != dimension + 1
            || state.hessian.size() != (dimension + 1) * (dimension + 1)) {
            throw std::invalid_argument("HELD2 Stage III evaluator derivative dimensions changed");
        }
        result.objective += fraction * state.objective;
        result.objective_gradient[offset] = state.objective;
        for (std::size_t local = 0; local < dimension + 1; ++local) {
            const std::size_t variable = offset + 1 + local;
            result.objective_gradient[variable] = fraction * state.gradient[local];
            const double balance_multiplier = local < dimension
                ? equality_multipliers[local + 1]
                : 0.0;
            add_symmetric(
                variable,
                offset,
                objective_factor * state.gradient[local] + balance_multiplier
            );
            for (std::size_t other = 0; other < dimension + 1; ++other) {
                add_symmetric(
                    variable,
                    offset + 1 + other,
                    objective_factor * fraction
                        * state.hessian[local * (dimension + 1) + other]
                        * (local == other ? 1.0 : 0.5)
                );
            }
        }
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

class Held2GeneralStageIIITnlp final : public Ipopt::TNLP {
public:
    Held2GeneralStageIIITnlp(
        Held2Coordinates coordinates,
        std::vector<double> feed,
        Held2StateEvaluator evaluator,
        std::vector<Held2StageIICandidate> candidates,
        std::array<double, 2> phase_coordinate_bounds,
        std::vector<double> initial
    )
        : coordinates_(std::move(coordinates)),
          feed_(std::move(feed)),
          evaluator_(std::move(evaluator)),
          candidates_(std::move(candidates)),
          phase_coordinate_bounds_(phase_coordinate_bounds),
          initial_(std::move(initial)) {}

    bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                      Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style) override {
        const std::size_t dimension = feed_.size();
        n = static_cast<Ipopt::Index>(initial_.size());
        m = static_cast<Ipopt::Index>(
            general_stage_iii_constraint_count(dimension, candidates_.size())
        );
        nnz_jac_g = static_cast<Ipopt::Index>(
            candidates_.size() * (1 + 3 * dimension)
        );
        nnz_h_lag = n * (n + 1) / 2;
        index_style = TNLP::C_STYLE;
        return true;
    }

    bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                         Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) override {
        const std::size_t dimension = feed_.size();
        const std::size_t block_size = dimension + 2;
        const std::size_t equality_count = dimension + 1;
        const std::size_t constraint_count =
            general_stage_iii_constraint_count(dimension, candidates_.size());
        if (n != static_cast<Ipopt::Index>(initial_.size())
            || m != static_cast<Ipopt::Index>(constraint_count)) {
            return false;
        }
        for (std::size_t phase = 0; phase < candidates_.size(); ++phase) {
            const std::size_t offset = phase * block_size;
            x_l[offset] = 0.0;
            x_u[offset] = 1.0;
            for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                x_l[offset + 1 + coordinate] = std::max(
                    coordinates_.independent_lower_bounds[coordinate],
                    candidates_[phase].independent_modified_fractions[coordinate]
                        - kCandidateRadius
                );
                x_u[offset + 1 + coordinate] = std::min(
                    coordinates_.independent_upper_bounds[coordinate],
                    candidates_[phase].independent_modified_fractions[coordinate]
                        + kCandidateRadius
                );
            }
            x_l[offset + 1 + dimension] = phase_coordinate_bounds_[0];
            x_u[offset + 1 + dimension] = phase_coordinate_bounds_[1];
        }
        std::fill(g_l, g_l + static_cast<std::ptrdiff_t>(equality_count), 0.0);
        std::fill(g_u, g_u + static_cast<std::ptrdiff_t>(equality_count), 0.0);
        const double composition_sum_upper =
            general_stage_iii_composition_sum_upper(coordinates_);
        for (std::size_t phase = 0; phase < candidates_.size(); ++phase) {
            g_l[equality_count + phase] = kConstraintLowerInfinity;
            g_u[equality_count + phase] = composition_sum_upper;
        }
        return true;
    }

    bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                            bool init_z, Ipopt::Number*, Ipopt::Number*,
                            Ipopt::Index m, bool init_lambda, Ipopt::Number*) override {
        if (n != static_cast<Ipopt::Index>(initial_.size())
            || m != static_cast<Ipopt::Index>(general_stage_iii_constraint_count(
                feed_.size(), candidates_.size()
            )) || !init_x
            || init_z || init_lambda) {
            return false;
        }
        std::copy(initial_.begin(), initial_.end(), x);
        return true;
    }

    bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool, Ipopt::Number& value) override {
        if (!precheck_objective_domain(n, x)) return false;
        try { value = evaluate(n, x, {}, 1.0).objective; return true; }
        catch (const std::exception& error) { callback_error_ = error.what(); return false; }
    }
    bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool,
                     Ipopt::Number* gradient) override {
        if (!precheck_objective_domain(n, x)) return false;
        try { const auto e = evaluate(n, x, {}, 1.0); std::copy(e.objective_gradient.begin(), e.objective_gradient.end(), gradient); return true; }
        catch (const std::exception& error) { callback_error_ = error.what(); return false; }
    }
    bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool, Ipopt::Index m,
                Ipopt::Number* constraints) override {
        try { const auto e = evaluate_algebraic(n, x); if (m != static_cast<Ipopt::Index>(e.constraints.size())) return false; std::copy(e.constraints.begin(), e.constraints.end(), constraints); return true; }
        catch (const std::exception& error) { callback_error_ = error.what(); return false; }
    }
    bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool, Ipopt::Index m,
                    Ipopt::Index nnz, Ipopt::Index* rows, Ipopt::Index* columns,
                    Ipopt::Number* values) override {
        const std::size_t dimension = feed_.size();
        const std::size_t block_size = dimension + 2;
        const Ipopt::Index expected = static_cast<Ipopt::Index>(candidates_.size() * (1 + 3 * dimension));
        if (m != static_cast<Ipopt::Index>(general_stage_iii_constraint_count(
                dimension, candidates_.size()
            )) || nnz != expected) return false;
        Ipopt::Index position = 0;
        if (values == nullptr) {
            for (std::size_t phase = 0; phase < candidates_.size(); ++phase) { rows[position] = 0; columns[position++] = static_cast<Ipopt::Index>(phase * block_size); }
            for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                for (std::size_t phase = 0; phase < candidates_.size(); ++phase) {
                    const std::size_t offset = phase * block_size;
                    rows[position] = static_cast<Ipopt::Index>(coordinate + 1); columns[position++] = static_cast<Ipopt::Index>(offset);
                    rows[position] = static_cast<Ipopt::Index>(coordinate + 1); columns[position++] = static_cast<Ipopt::Index>(offset + 1 + coordinate);
                }
            }
            for (std::size_t phase = 0; phase < candidates_.size(); ++phase) {
                const std::size_t offset = phase * block_size;
                for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                    rows[position] = static_cast<Ipopt::Index>(dimension + 1 + phase);
                    columns[position++] = static_cast<Ipopt::Index>(
                        offset + 1 + coordinate
                    );
                }
            }
            return true;
        }
        try {
            const auto e = evaluate_algebraic(n, x);
            for (std::size_t phase = 0; phase < candidates_.size(); ++phase) values[position++] = e.constraint_jacobian[phase * block_size];
            for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                for (std::size_t phase = 0; phase < candidates_.size(); ++phase) {
                    const std::size_t offset = phase * block_size;
                    values[position++] = e.constraint_jacobian[(coordinate + 1) * static_cast<std::size_t>(n) + offset];
                    values[position++] = e.constraint_jacobian[(coordinate + 1) * static_cast<std::size_t>(n) + offset + 1 + coordinate];
                }
            }
            for (std::size_t phase = 0; phase < candidates_.size(); ++phase) {
                for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                    values[position++] = 1.0;
                }
            }
            return true;
        } catch (const std::exception& error) { callback_error_ = error.what(); return false; }
    }
    bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool,
                Ipopt::Number objective_factor, Ipopt::Index m,
                const Ipopt::Number* lambda, bool, Ipopt::Index nnz,
                Ipopt::Index* rows, Ipopt::Index* columns, Ipopt::Number* values) override {
        if (m != static_cast<Ipopt::Index>(general_stage_iii_constraint_count(
                feed_.size(), candidates_.size()
            )) || nnz != n * (n + 1) / 2) return false;
        Ipopt::Index position = 0;
        if (values == nullptr) {
            for (Ipopt::Index row = 0; row < n; ++row) for (Ipopt::Index column = 0; column <= row; ++column) { rows[position] = row; columns[position++] = column; }
            return true;
        }
        if (!precheck_objective_domain(n, x)) return false;
        try {
            if (lambda == nullptr) return false;
            const std::vector<double> multipliers(lambda, lambda + m);
            const auto e = evaluate(n, x, multipliers, objective_factor);
            for (Ipopt::Index row = 0; row < n; ++row) for (Ipopt::Index column = 0; column <= row; ++column) values[position++] = e.lagrangian_hessian[static_cast<std::size_t>(row * n + column)];
            return true;
        } catch (const std::exception& error) { callback_error_ = error.what(); return false; }
    }

    void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n,
                           const Ipopt::Number* x, const Ipopt::Number* z_l,
                           const Ipopt::Number* z_u, Ipopt::Index m,
                           const Ipopt::Number*, const Ipopt::Number* lambda,
                           Ipopt::Number, const Ipopt::IpoptData*,
                           Ipopt::IpoptCalculatedQuantities*) override {
        solver_converged_ = status == Ipopt::SUCCESS || status == Ipopt::STOP_AT_ACCEPTABLE_POINT;
        if (x != nullptr && z_l != nullptr && z_u != nullptr && lambda != nullptr) {
            variables_.assign(x, x + n); lower_.assign(z_l, z_l + n); upper_.assign(z_u, z_u + n); multipliers_.assign(lambda, lambda + m);
        }
    }

    StageIIIRun result(const std::string& status) const {
        return {
            solver_converged_,
            status,
            callback_error_,
            variables_,
            multipliers_,
            lower_,
            upper_,
            recoverable_domain_rejection_count_,
            last_domain_rejection_,
        };
    }

private:
    bool precheck_objective_domain(Ipopt::Index n, const Ipopt::Number* x) {
        if (x == nullptr || n != static_cast<Ipopt::Index>(initial_.size())) {
            callback_error_ = "HELD2 Stage III callback dimensions changed";
            return false;
        }
        const std::vector<double> variables(x, x + n);
        if (general_stage_iii_simplex_domain_valid(
                coordinates_, candidates_.size(), variables
            )) {
            return true;
        }
        ++recoverable_domain_rejection_count_;
        last_domain_rejection_ =
            "independent modified fractions do not leave the declared dependent lower bound";
        return false;
    }

    Held2StageIIINlpEvaluation evaluate_algebraic(
        Ipopt::Index n,
        const Ipopt::Number* x
    ) const {
        if (x == nullptr || n != static_cast<Ipopt::Index>(initial_.size())) {
            throw std::invalid_argument("HELD2 Stage III callback dimensions changed");
        }
        return evaluate_general_stage_iii_algebraic(
            feed_, candidates_.size(), std::vector<double>(x, x + n)
        );
    }

    Held2StageIIINlpEvaluation evaluate(Ipopt::Index n, const Ipopt::Number* x,
                                        std::vector<double> multipliers,
                                        double factor) const {
        if (multipliers.empty()) {
            multipliers.assign(
                general_stage_iii_constraint_count(feed_.size(), candidates_.size()),
                0.0
            );
        }
        return evaluate_general_stage_iii(coordinates_, feed_, evaluator_, candidates_.size(), std::vector<double>(x, x + n), multipliers, factor);
    }
    Held2Coordinates coordinates_; std::vector<double> feed_; Held2StateEvaluator evaluator_;
    std::vector<Held2StageIICandidate> candidates_; std::array<double, 2> phase_coordinate_bounds_{};
    std::vector<double> initial_; bool solver_converged_ = false; std::string callback_error_;
    int recoverable_domain_rejection_count_ = 0;
    std::string last_domain_rejection_;
    std::vector<double> variables_, multipliers_, lower_, upper_;
};

StageIIIRun run_general_stage_iii(
    const Held2Coordinates& coordinates,
    const std::vector<double>& feed,
    const Held2StateEvaluator& evaluator,
    const std::vector<Held2StageIICandidate>& candidates,
    const std::array<double, 2>& phase_coordinate_bounds,
    const std::vector<double>& initial
) {
    auto* raw = new Held2GeneralStageIIITnlp(coordinates, feed, evaluator, candidates, phase_coordinate_bounds, initial);
    Ipopt::SmartPtr<Ipopt::TNLP> problem = raw;
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
    if (application->Initialize() != Ipopt::Solve_Succeeded) return {};
    const auto status = application->OptimizeTNLP(problem);
    return raw->result(ipopt_status_name(status));
}

}  // namespace

std::tuple<
    int,
    int,
    int,
    std::vector<double>,
    std::vector<double>,
    std::vector<int>,
    std::vector<int>,
    std::vector<double>> inspect_held2_stage_iii_tnlp(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const std::vector<Held2StageIICandidate>& candidates,
    const std::array<double, 2>& phase_coordinate_bounds,
    const std::vector<double>& variables
) {
    const std::vector<double> modified_feed =
        held2_transform_physical_fractions(coordinates, physical_feed);
    const std::vector<std::size_t> positions =
        independent_retained_positions(coordinates);
    std::vector<double> feed;
    for (std::size_t position : positions) {
        feed.push_back(modified_feed[position]);
    }
    const Held2StateEvaluator unused_evaluator = [](
        const std::vector<double>&,
        double
    ) -> Held2StateEvaluation {
        throw std::logic_error("HELD2 Stage III schema inspection invoked the evaluator");
    };
    Held2GeneralStageIIITnlp problem(
        coordinates,
        feed,
        unused_evaluator,
        candidates,
        phase_coordinate_bounds,
        variables
    );
    Ipopt::Index variable_count = 0;
    Ipopt::Index constraint_count = 0;
    Ipopt::Index jacobian_nonzero_count = 0;
    Ipopt::Index hessian_nonzero_count = 0;
    Ipopt::TNLP::IndexStyleEnum index_style = Ipopt::TNLP::FORTRAN_STYLE;
    if (!problem.get_nlp_info(
            variable_count,
            constraint_count,
            jacobian_nonzero_count,
            hessian_nonzero_count,
            index_style
        ) || index_style != Ipopt::TNLP::C_STYLE) {
        throw std::logic_error("HELD2 Stage III schema inspection failed");
    }
    std::vector<double> variable_lower(static_cast<std::size_t>(variable_count));
    std::vector<double> variable_upper(static_cast<std::size_t>(variable_count));
    std::vector<double> constraint_lower(static_cast<std::size_t>(constraint_count));
    std::vector<double> constraint_upper(static_cast<std::size_t>(constraint_count));
    if (!problem.get_bounds_info(
            variable_count,
            variable_lower.data(),
            variable_upper.data(),
            constraint_count,
            constraint_lower.data(),
            constraint_upper.data()
        )) {
        throw std::logic_error("HELD2 Stage III bound inspection failed");
    }
    std::vector<Ipopt::Index> raw_rows(
        static_cast<std::size_t>(jacobian_nonzero_count)
    );
    std::vector<Ipopt::Index> raw_columns(
        static_cast<std::size_t>(jacobian_nonzero_count)
    );
    std::vector<double> jacobian_values(
        static_cast<std::size_t>(jacobian_nonzero_count)
    );
    if (!problem.eval_jac_g(
            variable_count,
            variables.data(),
            true,
            constraint_count,
            jacobian_nonzero_count,
            raw_rows.data(),
            raw_columns.data(),
            nullptr
        ) || !problem.eval_jac_g(
            variable_count,
            variables.data(),
            false,
            constraint_count,
            jacobian_nonzero_count,
            nullptr,
            nullptr,
            jacobian_values.data()
        )) {
        throw std::logic_error("HELD2 Stage III Jacobian inspection failed");
    }
    std::vector<int> rows(raw_rows.begin(), raw_rows.end());
    std::vector<int> columns(raw_columns.begin(), raw_columns.end());
    return {
        static_cast<int>(variable_count),
        static_cast<int>(constraint_count),
        static_cast<int>(jacobian_nonzero_count),
        std::move(constraint_lower),
        std::move(constraint_upper),
        std::move(rows),
        std::move(columns),
        std::move(jacobian_values),
    };
}

std::tuple<bool, int, std::string, std::string>
probe_held2_stage_iii_objective_trial(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const std::vector<Held2StageIICandidate>& candidates,
    const Held2StateEvaluator& evaluator,
    const std::array<double, 2>& phase_coordinate_bounds,
    const std::vector<double>& variables
) {
    const std::vector<double> modified_feed =
        held2_transform_physical_fractions(coordinates, physical_feed);
    const std::vector<std::size_t> positions =
        independent_retained_positions(coordinates);
    std::vector<double> feed;
    for (std::size_t position : positions) {
        feed.push_back(modified_feed[position]);
    }
    Held2GeneralStageIIITnlp problem(
        coordinates,
        feed,
        evaluator,
        candidates,
        phase_coordinate_bounds,
        variables
    );
    double objective = 0.0;
    const bool accepted = problem.eval_f(
        static_cast<Ipopt::Index>(variables.size()),
        variables.data(),
        true,
        objective
    );
    const StageIIIRun evidence = problem.result("not_run");
    return {
        accepted,
        evidence.recoverable_domain_rejection_count,
        evidence.last_domain_rejection,
        evidence.callback_error,
    };
}

Held2StageIIINlpEvaluation evaluate_held2_stage_iii_nlp(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const Held2StateEvaluator& evaluator,
    std::size_t phase_count,
    const std::vector<double>& variables,
    const std::vector<double>& equality_multipliers
) {
    const std::vector<double> modified_feed = held2_transform_physical_fractions(coordinates, physical_feed);
    const std::vector<std::size_t> positions = independent_retained_positions(coordinates);
    std::vector<double> feed;
    for (std::size_t position : positions) feed.push_back(modified_feed[position]);
    return evaluate_general_stage_iii(coordinates, feed, evaluator, phase_count, variables, equality_multipliers, 1.0);
}

Held2StageIIIResult solve_held2_stage_iii(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const std::vector<Held2StageIICandidate>& candidates,
    const Held2StateEvaluator& evaluator,
    const std::array<double, 2>& phase_coordinate_bounds
) {
    Held2StageIIIResult result;
    result.input_candidate_count = static_cast<int>(candidates.size());
    const std::size_t dimension = coordinates.independent_indices.size();
    const std::size_t block_size = dimension + 2;
    if (candidates.size() < 2) { result.failure_reason = "candidate_set_incomplete"; return result; }
    const std::vector<double> modified_feed = held2_transform_physical_fractions(coordinates, physical_feed);
    const std::vector<std::size_t> positions = independent_retained_positions(coordinates);
    std::vector<double> feed; for (std::size_t position : positions) feed.push_back(modified_feed[position]);
    std::vector<double> initial(block_size * candidates.size(), 0.0);
    for (std::size_t phase = 0; phase < candidates.size(); ++phase) {
        if (candidates[phase].independent_modified_fractions.size() != dimension) { result.failure_reason = "candidate_dimension_changed"; return result; }
        const std::size_t offset = phase * block_size;
        initial[offset] = 1.0 / static_cast<double>(candidates.size());
        std::copy(candidates[phase].independent_modified_fractions.begin(), candidates[phase].independent_modified_fractions.end(), initial.begin() + static_cast<std::ptrdiff_t>(offset + 1));
        initial[offset + 1 + dimension] = candidates[phase].phase_coordinate;
    }
    const StageIIIRun run = run_general_stage_iii(coordinates, feed, evaluator, candidates, phase_coordinate_bounds, initial);
    result.solver_status = run.solver_status;
    if (!run.solver_converged || !run.callback_error.empty()
        || run.variables.size() != initial.size()) {
        result.numerical_status = "not_converged";
        if (!run.callback_error.empty()) {
            result.failure_reason = run.callback_error;
        } else if (run.recoverable_domain_rejection_count > 0) {
            result.failure_reason = "stage_iii_simplex_domain_exhausted";
        } else {
            result.failure_reason = "stage_iii_solver_not_converged";
        }
        return result;
    }
    Held2StageIIINlpEvaluation nlp;
    try {
        nlp = evaluate_general_stage_iii(
            coordinates,
            feed,
            evaluator,
            candidates.size(),
            run.variables,
            run.equality_multipliers,
            1.0
        );
    } catch (const std::exception& error) {
        result.numerical_status = "not_converged";
        result.failure_reason = error.what();
        return result;
    }
    std::vector<double> kkt = nlp.lagrangian_gradient;
    if (run.lower_bound_multipliers.size() != kkt.size() || run.upper_bound_multipliers.size() != kkt.size()) { result.numerical_status = "not_converged"; result.failure_reason = "stage_iii_multiplier_evidence_missing"; return result; }
    for (std::size_t index = 0; index < kkt.size(); ++index) { kkt[index] -= run.lower_bound_multipliers[index]; kkt[index] += run.upper_bound_multipliers[index]; }
    result.kkt_stationarity_inf_norm = maximum_abs(kkt);
    double constraint_violation = 0.0;
    for (std::size_t row = 0; row < dimension + 1; ++row) {
        constraint_violation = std::max(
            constraint_violation, std::abs(nlp.constraints[row])
        );
    }
    const double composition_sum_upper =
        general_stage_iii_composition_sum_upper(coordinates);
    for (std::size_t phase = 0; phase < candidates.size(); ++phase) {
        constraint_violation = std::max(
            constraint_violation,
            std::max(
                0.0,
                nlp.constraints[dimension + 1 + phase] - composition_sum_upper
            )
        );
    }
    if (constraint_violation > kNumericalTolerance
        || result.kkt_stationarity_inf_norm > kKktTolerance) {
        result.numerical_status = "not_converged";
        result.failure_reason = "stage_iii_numerical_certificate_failed";
        return result;
    }
    result.numerical_status = "converged"; result.objective = nlp.objective;
    std::vector<Held2StateEvaluation> states;
    for (std::size_t phase = 0; phase < candidates.size(); ++phase) {
        const std::size_t offset = phase * block_size;
        const std::vector<double> independent(run.variables.begin() + static_cast<std::ptrdiff_t>(offset + 1), run.variables.begin() + static_cast<std::ptrdiff_t>(offset + 1 + dimension));
        const auto state = evaluator(independent, run.variables[offset + 1 + dimension]);
        const double fraction = run.variables[offset];
        if (fraction <= kActivePhaseTolerance) continue;
        bool merged = false;
        for (std::size_t existing = 0; existing < result.phases.size(); ++existing) {
            if (maximum_abs_difference(result.phases[existing].modified_fractions, state.modified_fractions) <= kMergeTolerance && std::abs(result.phases[existing].volume - state.volume) <= kMergeTolerance) {
                result.phases[existing].phase_fraction += fraction; ++result.retired_duplicate_count; merged = true; break;
            }
        }
        if (!merged) { result.phases.push_back({fraction, state.modified_fractions, state.physical_amounts, state.volume}); states.push_back(state); }
    }
    if (result.phases.size() < 2) { result.physical_status = "rejected"; result.failure_reason = "collapsed_phase_set"; return result; }
    std::vector<double> modified_balance(modified_feed.size(), 0.0), ordinary_balance(physical_feed.size(), 0.0);
    result.minimum_phase_distance = std::numeric_limits<double>::infinity();
    for (const auto& phase : result.phases) {
        for (std::size_t i = 0; i < modified_balance.size(); ++i) modified_balance[i] += phase.phase_fraction * phase.modified_fractions[i];
        for (std::size_t i = 0; i < ordinary_balance.size(); ++i) ordinary_balance[i] += phase.phase_fraction * phase.physical_fractions[i];
        result.phase_charge_inf_norm = std::max(result.phase_charge_inf_norm, std::abs(charge_residual(coordinates.charges, phase.physical_fractions)));
    }
    result.modified_balance_inf_norm = maximum_abs_difference(modified_balance, modified_feed);
    result.ordinary_balance_inf_norm = maximum_abs_difference(ordinary_balance, physical_feed);
    for (const auto& state : states) result.pressure_stationarity_inf_norm = std::max(result.pressure_stationarity_inf_norm, std::abs(state.pressure_stationarity_relative));
    for (std::size_t left = 0; left < result.phases.size(); ++left) for (std::size_t right = left + 1; right < result.phases.size(); ++right) result.minimum_phase_distance = std::min(result.minimum_phase_distance, maximum_abs_difference(result.phases[left].modified_fractions, result.phases[right].modified_fractions));
    for (std::size_t retained = 0; retained < modified_feed.size(); ++retained) {
        const double lower_bound =
            kHeld2ModifiedLowerScale * coordinates.modified_factors[retained];
        bool trace = false; for (const auto& phase : result.phases) trace = trace || phase.modified_fractions[retained] <= 10.0 * lower_bound;
        if (trace) { ++result.trace_component_count; continue; }
        ++result.certified_modified_potential_count;
        for (std::size_t left = 0; left < states.size(); ++left) for (std::size_t right = left + 1; right < states.size(); ++right) {
            const double gap = std::abs(states[left].modified_potentials[retained] - states[right].modified_potentials[retained]);
            const double scale = 1.0 + std::max(std::abs(states[left].modified_potentials[retained]), std::abs(states[right].modified_potentials[retained]));
            result.modified_potential_mixed_gap = std::max(result.modified_potential_mixed_gap, gap / scale);
        }
    }
    if (result.trace_component_count > 0) { result.trace_refinement_status = "complementarity_refinement_required"; result.failure_reason = "trace_component_requires_log_refinement"; return result; }
    result.trace_refinement_status = "not_required";
    const bool physical = result.modified_balance_inf_norm < kNumericalTolerance && result.ordinary_balance_inf_norm < kNumericalTolerance && result.phase_charge_inf_norm < kNumericalTolerance && result.pressure_stationarity_inf_norm < kNumericalTolerance && result.modified_potential_mixed_gap < kNumericalTolerance && result.minimum_phase_distance > kDistinctPhaseTolerance;
    if (!physical) { result.physical_status = "rejected"; result.failure_reason = "stage_iii_physical_certificate_failed"; return result; }
    result.physical_status = "accepted"; result.feedback = "none"; return result;
}

}  // namespace epcsaft_equilibrium
