#include "held2.hpp"
#include "held2_progress.hpp"
#include "held2_tolerances.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>
#include <Highs.h>

namespace epcsaft_equilibrium {
namespace {

constexpr double kCandidateRadius = 1.0e-3;
constexpr double kConstraintLowerInfinity = -1.0e19;
constexpr double kHeld2ModifiedLowerScale = 1.0e-10;

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

double maximum_abs(const std::vector<double>& values) {
    double result = 0.0;
    for (double value : values) result = std::max(result, std::abs(value));
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
        case Ipopt::Solve_Succeeded: return "solve_succeeded";
        case Ipopt::Solved_To_Acceptable_Level: return "solved_to_acceptable_level";
        case Ipopt::Infeasible_Problem_Detected: return "infeasible_problem_detected";
        case Ipopt::Search_Direction_Becomes_Too_Small: return "search_direction_too_small";
        case Ipopt::Maximum_Iterations_Exceeded: return "maximum_iterations_exceeded";
        case Ipopt::Invalid_Number_Detected: return "invalid_number_detected";
        default: return "ipopt_status_" + std::to_string(static_cast<int>(status));
    }
}

}  // namespace

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

bool remaining_candidate_balance_feasible(
    const std::vector<double>& feed,
    const std::vector<Held2StageIICandidate>& candidates,
    std::size_t removed
) {
    if (removed >= candidates.size() || candidates.size() <= 2) {
        return false;
    }
    if (feed.size() == 1) {
        double lower = std::numeric_limits<double>::infinity();
        double upper = -std::numeric_limits<double>::infinity();
        for (std::size_t phase = 0; phase < candidates.size(); ++phase) {
            if (phase == removed
                || candidates[phase].independent_modified_fractions.size() != 1) {
                continue;
            }
            lower = std::min(
                lower, candidates[phase].independent_modified_fractions[0]
            );
            upper = std::max(
                upper, candidates[phase].independent_modified_fractions[0]
            );
        }
        return feed[0] >= lower - kHeld2Stage3ModifiedBalance.atol
            && feed[0] <= upper + kHeld2Stage3ModifiedBalance.atol;
    }
    const std::size_t column_count = candidates.size() - 1;
    HighsModel model;
    model.lp_.num_col_ = static_cast<HighsInt>(column_count);
    model.lp_.num_row_ = static_cast<HighsInt>(feed.size() + 1);
    model.lp_.col_cost_.assign(column_count, 0.0);
    model.lp_.col_lower_.assign(column_count, 0.0);
    model.lp_.col_upper_.assign(column_count, 1.0);
    model.lp_.row_lower_.reserve(feed.size() + 1);
    model.lp_.row_upper_.reserve(feed.size() + 1);
    model.lp_.row_lower_.push_back(1.0);
    model.lp_.row_upper_.push_back(1.0);
    for (double value : feed) {
        model.lp_.row_lower_.push_back(value);
        model.lp_.row_upper_.push_back(value);
    }
    model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
    model.lp_.a_matrix_.start_.push_back(0);
    for (std::size_t phase = 0; phase < candidates.size(); ++phase) {
        if (phase == removed) {
            continue;
        }
        if (candidates[phase].independent_modified_fractions.size() != feed.size()) {
            return false;
        }
        model.lp_.a_matrix_.index_.push_back(0);
        model.lp_.a_matrix_.value_.push_back(1.0);
        for (std::size_t coordinate = 0; coordinate < feed.size(); ++coordinate) {
            model.lp_.a_matrix_.index_.push_back(
                static_cast<HighsInt>(coordinate + 1)
            );
            model.lp_.a_matrix_.value_.push_back(
                candidates[phase].independent_modified_fractions[coordinate]
            );
        }
        model.lp_.a_matrix_.start_.push_back(
            static_cast<HighsInt>(model.lp_.a_matrix_.index_.size())
        );
    }
    Highs highs;
    if (highs.setOptionValue("output_flag", false) == HighsStatus::kError
        || highs.setOptionValue("threads", 1) == HighsStatus::kError
        || highs.passModel(model) == HighsStatus::kError
        || highs.run() == HighsStatus::kError) {
        return false;
    }
    return highs.getModelStatus() == HighsModelStatus::kOptimal
        && highs.getInfo().primal_solution_status == kSolutionStatusFeasible;
}

Held2StageIIILifecycleStep make_general_lifecycle_step(
    int solve_index,
    std::size_t active_candidate_count,
    std::string action,
    std::string solver_status,
    std::string decision_reason,
    int candidate_index = -1,
    double phase_fraction = 0.0,
    double lower_bound_multiplier = 0.0,
    double reduced_derivative = 0.0,
    double complementarity = 0.0,
    std::vector<double> candidate = {},
    double candidate_volume = 0.0
) {
    return {
        solve_index,
        static_cast<int>(active_candidate_count),
        candidate_index,
        std::move(action),
        phase_fraction,
        lower_bound_multiplier,
        reduced_derivative,
        complementarity,
        std::move(candidate),
        candidate_volume,
        std::move(solver_status),
        std::move(decision_reason),
    };
}

class Held2GeneralStageIIITnlp final : public Ipopt::TNLP {
public:
    Held2GeneralStageIIITnlp(
        Held2Coordinates coordinates,
        std::vector<double> feed,
        Held2StateEvaluator evaluator,
        std::vector<Held2StageIICandidate> candidates,
        std::vector<std::array<double, 2>> phase_coordinate_bounds,
        std::vector<double> initial,
        Held2ProgressObserver* observer,
        int solve_index
    )
        : coordinates_(std::move(coordinates)),
          feed_(std::move(feed)),
          evaluator_(std::move(evaluator)),
          candidates_(std::move(candidates)),
          phase_coordinate_bounds_(phase_coordinate_bounds),
          initial_(std::move(initial)),
          observer_(observer),
          solve_index_(solve_index) {}

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
            x_l[offset + 1 + dimension] = phase_coordinate_bounds_[phase][0];
            x_u[offset + 1 + dimension] = phase_coordinate_bounds_[phase][1];
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

    bool intermediate_callback(
        Ipopt::AlgorithmMode,
        Ipopt::Index iteration,
        Ipopt::Number objective,
        Ipopt::Number primal_residual,
        Ipopt::Number dual_residual,
        Ipopt::Number barrier,
        Ipopt::Number step_norm,
        Ipopt::Number,
        Ipopt::Number dual_step,
        Ipopt::Number primal_step,
        Ipopt::Index line_search_steps,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*
    ) override {
        Held2ProgressEvent progress;
        progress.kind = Held2ProgressKind::LocalIteration;
        progress.stage = "STAGE III IPOPT";
        progress.iteration = static_cast<int>(iteration);
        progress.attempt = solve_index_;
        progress.objective = objective;
        progress.primal_residual = primal_residual;
        progress.dual_residual = dual_residual;
        progress.complementarity = barrier;
        progress.step_norm = step_norm;
        progress.dual_step = dual_step;
        progress.primal_step = primal_step;
        progress.line_search_steps = static_cast<int>(line_search_steps);
        observe_held2(observer_, progress);
        return true;
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
    std::vector<Held2StageIICandidate> candidates_;
    std::vector<std::array<double, 2>> phase_coordinate_bounds_;
    std::vector<double> initial_; bool solver_converged_ = false; std::string callback_error_;
    int recoverable_domain_rejection_count_ = 0;
    std::string last_domain_rejection_;
    std::vector<double> variables_, multipliers_, lower_, upper_;
    Held2ProgressObserver* observer_ = nullptr;
    int solve_index_ = -1;
};

StageIIIRun run_general_stage_iii(
    const Held2Coordinates& coordinates,
    const std::vector<double>& feed,
    const Held2StateEvaluator& evaluator,
    const std::vector<Held2StageIICandidate>& candidates,
    const std::vector<std::array<double, 2>>& phase_coordinate_bounds,
    const std::vector<double>& initial,
    Held2ProgressObserver* observer,
    int solve_index
) {
    auto* raw = new Held2GeneralStageIIITnlp(
        coordinates,
        feed,
        evaluator,
        candidates,
        phase_coordinate_bounds,
        initial,
        observer,
        solve_index
    );
    Ipopt::SmartPtr<Ipopt::TNLP> problem = raw;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application = IpoptApplicationFactory();
    application->Options()->SetStringValue("option_file_name", "");
    application->Options()->SetIntegerValue("print_level", 0);
    application->Options()->SetStringValue("sb", "yes");
    application->Options()->SetIntegerValue("max_iter", 300);
    application->Options()->SetNumericValue("tol", kHeld2IpoptTarget.atol);
    application->Options()->SetNumericValue(
        "acceptable_tol", kHeld2IpoptAcceptable.atol
    );
    application->Options()->SetIntegerValue("acceptable_iter", 0);
    application->Options()->SetNumericValue(
        "constr_viol_tol", kHeld2IpoptConstraint.atol
    );
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
    const std::vector<std::array<double, 2>>& phase_coordinate_bounds,
    Held2ProgressObserver* observer
) {
    Held2StageIIIResult result;
    result.input_candidate_count = static_cast<int>(candidates.size());
    const std::size_t dimension = coordinates.independent_indices.size();
    const std::size_t block_size = dimension + 2;
    if (candidates.size() < 2) {
        result.failure_reason = "candidate_set_incomplete";
        return result;
    }
    if (phase_coordinate_bounds.size() != candidates.size()) {
        result.failure_reason = "candidate_bound_count_changed";
        return result;
    }

    const std::vector<double> modified_feed =
        held2_transform_physical_fractions(coordinates, physical_feed);
    const std::vector<std::size_t> positions =
        independent_retained_positions(coordinates);
    std::vector<double> feed;
    for (std::size_t position : positions) {
        feed.push_back(modified_feed[position]);
    }

    std::vector<Held2StageIICandidate> active_candidates = candidates;
    std::vector<std::array<double, 2>> active_bounds = phase_coordinate_bounds;
    StageIIIRun accepted_run;
    Held2StageIIINlpEvaluation accepted_nlp;
    std::vector<Held2StateEvaluation> accepted_states;
    bool active_set_accepted = false;

    for (std::size_t lifecycle = 0; lifecycle <= candidates.size(); ++lifecycle) {
        std::vector<double> initial(block_size * active_candidates.size(), 0.0);
        for (std::size_t phase = 0; phase < active_candidates.size(); ++phase) {
            const auto& candidate = active_candidates[phase];
            if (candidate.independent_modified_fractions.size() != dimension) {
                result.failure_reason = "candidate_dimension_changed";
                return result;
            }
            const std::size_t offset = phase * block_size;
            initial[offset] = 1.0 / static_cast<double>(active_candidates.size());
            std::copy(
                candidate.independent_modified_fractions.begin(),
                candidate.independent_modified_fractions.end(),
                initial.begin() + static_cast<std::ptrdiff_t>(offset + 1)
            );
            initial[offset + 1 + dimension] = candidate.phase_coordinate;
        }

        StageIIIRun run = run_general_stage_iii(
            coordinates,
            feed,
            evaluator,
            active_candidates,
            active_bounds,
            initial,
            observer,
            result.stage_iii_solve_count + 1
        );
        ++result.stage_iii_solve_count;
        result.solver_status = run.solver_status;
        if (!run.solver_converged || !run.callback_error.empty()
            || run.variables.size() != initial.size()) {
            result.numerical_status = "not_converged";
            result.failure_reason = !run.callback_error.empty()
                ? run.callback_error
                : run.recoverable_domain_rejection_count > 0
                    ? "stage_iii_simplex_domain_exhausted"
                    : "stage_iii_solver_not_converged";
            result.lifecycle.push_back(make_general_lifecycle_step(
                result.stage_iii_solve_count,
                active_candidates.size(),
                "solve_failed",
                run.solver_status,
                result.failure_reason
            ));
            return result;
        }

        Held2StageIIINlpEvaluation nlp;
        try {
            nlp = evaluate_general_stage_iii(
                coordinates,
                feed,
                evaluator,
                active_candidates.size(),
                run.variables,
                run.equality_multipliers,
                1.0
            );
        } catch (const std::exception& error) {
            result.numerical_status = "not_converged";
            result.failure_reason = error.what();
            result.lifecycle.push_back(make_general_lifecycle_step(
                result.stage_iii_solve_count,
                active_candidates.size(),
                "numerical_certificate_failed",
                run.solver_status,
                result.failure_reason
            ));
            return result;
        }

        std::vector<double> kkt = nlp.lagrangian_gradient;
        if (run.lower_bound_multipliers.size() != kkt.size()
            || run.upper_bound_multipliers.size() != kkt.size()) {
            result.numerical_status = "not_converged";
            result.failure_reason = "stage_iii_multiplier_evidence_missing";
            result.lifecycle.push_back(make_general_lifecycle_step(
                result.stage_iii_solve_count,
                active_candidates.size(),
                "numerical_certificate_failed",
                run.solver_status,
                result.failure_reason
            ));
            return result;
        }

        double constraint_violation = 0.0;
        for (std::size_t row = 0; row < dimension + 1; ++row) {
            constraint_violation = std::max(
                constraint_violation, std::abs(nlp.constraints[row])
            );
        }
        const double composition_sum_upper =
            general_stage_iii_composition_sum_upper(coordinates);
        result.dual_sign_violation_inf_norm = 0.0;
        result.bound_complementarity_inf_norm = 0.0;
        result.minimum_phase_fraction = std::numeric_limits<double>::infinity();
        for (std::size_t phase = 0; phase < active_candidates.size(); ++phase) {
            const std::size_t offset = phase * block_size;
            constraint_violation = std::max(
                constraint_violation,
                std::max(
                    0.0,
                    nlp.constraints[dimension + 1 + phase]
                        - composition_sum_upper
                )
            );
            for (std::size_t local = 0; local < block_size; ++local) {
                const std::size_t index = offset + local;
                double lower = 0.0;
                double upper = 1.0;
                if (local > 0 && local <= dimension) {
                    const std::size_t coordinate = local - 1;
                    lower = std::max(
                        coordinates.independent_lower_bounds[coordinate],
                        active_candidates[phase]
                            .independent_modified_fractions[coordinate]
                            - kCandidateRadius
                    );
                    upper = std::min(
                        coordinates.independent_upper_bounds[coordinate],
                        active_candidates[phase]
                            .independent_modified_fractions[coordinate]
                            + kCandidateRadius
                    );
                } else if (local == block_size - 1) {
                    lower = active_bounds[phase][0];
                    upper = active_bounds[phase][1];
                }
                kkt[index] -= run.lower_bound_multipliers[index];
                kkt[index] += run.upper_bound_multipliers[index];
                result.dual_sign_violation_inf_norm = std::max({
                    result.dual_sign_violation_inf_norm,
                    -run.lower_bound_multipliers[index],
                    -run.upper_bound_multipliers[index],
                });
                result.bound_complementarity_inf_norm = std::max({
                    result.bound_complementarity_inf_norm,
                    std::abs((run.variables[index] - lower)
                        * run.lower_bound_multipliers[index]),
                    std::abs((upper - run.variables[index])
                        * run.upper_bound_multipliers[index]),
                });
            }
            result.minimum_phase_fraction = std::min(
                result.minimum_phase_fraction, run.variables[offset]
            );
        }
        result.kkt_stationarity_inf_norm = maximum_abs(kkt);
        result.kkt_evidence_available = true;
        if (!audit_held2_tolerance(
                kHeld2Stage3ModifiedBalance, constraint_violation
            ).passed
            || !audit_held2_tolerance(
                kHeld2Stage3Stationarity,
                result.kkt_stationarity_inf_norm
            ).passed
            || !audit_held2_tolerance(
                kHeld2Stage3DualSign,
                result.dual_sign_violation_inf_norm
            ).passed
            || !audit_held2_tolerance(
                kHeld2Stage3Complementarity,
                result.bound_complementarity_inf_norm
            ).passed) {
            result.numerical_status = "not_converged";
            result.failure_reason = "stage_iii_numerical_certificate_failed";
            result.lifecycle.push_back(make_general_lifecycle_step(
                result.stage_iii_solve_count,
                active_candidates.size(),
                "numerical_certificate_failed",
                run.solver_status,
                result.failure_reason,
                -1,
                0.0,
                0.0,
                0.0,
                result.bound_complementarity_inf_norm
            ));
            return result;
        }
        result.numerical_status = "converged";
        Held2ProgressEvent numerical_progress;
        numerical_progress.kind = Held2ProgressKind::Certificate;
        numerical_progress.stage = "STAGE III NUMERICAL";
        numerical_progress.status = "passed";
        numerical_progress.primal_residual = constraint_violation;
        numerical_progress.dual_residual = result.kkt_stationarity_inf_norm;
        numerical_progress.complementarity =
            result.bound_complementarity_inf_norm;
        observe_held2(observer, numerical_progress);

        bool active_set_changed = false;
        for (std::size_t phase = 0; phase < active_candidates.size(); ++phase) {
            const std::size_t offset = phase * block_size;
            const bool remaining_feasible = remaining_candidate_balance_feasible(
                feed, active_candidates, phase
            );
            const Held2StageIIIRetirementDecision decision =
                held2_stage_iii_retirement_decision(
                    run.variables[offset],
                    run.lower_bound_multipliers[offset],
                    run.upper_bound_multipliers[offset],
                    nlp.lagrangian_gradient[offset],
                    remaining_feasible
                );
            result.lifecycle.push_back(make_general_lifecycle_step(
                result.stage_iii_solve_count,
                active_candidates.size(),
                decision.retire ? "retire_kkt_inactive" : "retain_phase",
                run.solver_status,
                decision.reason,
                static_cast<int>(phase),
                run.variables[offset],
                run.lower_bound_multipliers[offset],
                nlp.lagrangian_gradient[offset],
                decision.complementarity_inf_norm,
                active_candidates[phase].independent_modified_fractions,
                active_candidates[phase].volume
            ));
            if (!decision.retire) {
                continue;
            }
            active_candidates.erase(
                active_candidates.begin() + static_cast<std::ptrdiff_t>(phase)
            );
            active_bounds.erase(
                active_bounds.begin() + static_cast<std::ptrdiff_t>(phase)
            );
            ++result.retired_inactive_count;
            ++result.active_set_resolve_count;
            active_set_changed = true;
            break;
        }
        if (active_set_changed) {
            continue;
        }

        std::vector<Held2StateEvaluation> states;
        states.reserve(active_candidates.size());
        for (std::size_t phase = 0; phase < active_candidates.size(); ++phase) {
            const std::size_t offset = phase * block_size;
            const std::vector<double> independent(
                run.variables.begin() + static_cast<std::ptrdiff_t>(offset + 1),
                run.variables.begin()
                    + static_cast<std::ptrdiff_t>(offset + 1 + dimension)
            );
            states.push_back(evaluator(
                independent, run.variables[offset + 1 + dimension]
            ));
        }
        for (std::size_t left = 0;
             left < active_candidates.size() && !active_set_changed;
             ++left) {
            for (std::size_t right = left + 1;
                 right < active_candidates.size();
                 ++right) {
                const double merge_distance = std::max(
                    maximum_abs_difference(
                        states[left].modified_fractions,
                        states[right].modified_fractions
                    ),
                    std::abs(
                        std::log(states[left].volume)
                            - std::log(states[right].volume)
                    )
                );
                if (!audit_held2_tolerance(
                        kHeld2PhaseMerge,
                        merge_distance
                    ).passed) {
                    continue;
                }
                const std::size_t removed =
                    run.variables[left * block_size]
                        <= run.variables[right * block_size]
                    ? left : right;
                if (!remaining_candidate_balance_feasible(
                        feed, active_candidates, removed
                    )) {
                    continue;
                }
                result.lifecycle.push_back(make_general_lifecycle_step(
                    result.stage_iii_solve_count,
                    active_candidates.size(),
                    "merge_duplicate",
                    run.solver_status,
                    "duplicate_state_certified",
                    static_cast<int>(removed),
                    run.variables[removed * block_size],
                    run.lower_bound_multipliers[removed * block_size],
                    nlp.lagrangian_gradient[removed * block_size],
                    std::max(
                        std::abs(run.variables[removed * block_size]
                            * run.lower_bound_multipliers[removed * block_size]),
                        std::abs((1.0 - run.variables[removed * block_size])
                            * run.upper_bound_multipliers[removed * block_size])
                    ),
                    active_candidates[removed]
                        .independent_modified_fractions,
                    active_candidates[removed].volume
                ));
                active_candidates.erase(
                    active_candidates.begin()
                        + static_cast<std::ptrdiff_t>(removed)
                );
                active_bounds.erase(
                    active_bounds.begin() + static_cast<std::ptrdiff_t>(removed)
                );
                ++result.retired_duplicate_count;
                ++result.active_set_resolve_count;
                active_set_changed = true;
                break;
            }
        }
        if (active_set_changed) {
            continue;
        }

        result.lifecycle.push_back(make_general_lifecycle_step(
            result.stage_iii_solve_count,
            active_candidates.size(),
            "accept_active_set",
            run.solver_status,
            "active_set_certified"
        ));
        accepted_run = std::move(run);
        accepted_nlp = std::move(nlp);
        accepted_states = std::move(states);
        active_set_accepted = true;
        break;
    }

    if (!active_set_accepted) {
        result.numerical_status = "not_converged";
        result.failure_reason = "stage_iii_active_set_lifecycle_exhausted";
        return result;
    }

    result.objective = accepted_nlp.objective;
    for (std::size_t phase = 0; phase < active_candidates.size(); ++phase) {
        const std::size_t offset = phase * block_size;
        result.phases.push_back({
            accepted_run.variables[offset],
            accepted_states[phase].modified_fractions,
            accepted_states[phase].physical_amounts,
            accepted_states[phase].volume,
        });
    }
    if (result.phases.size() < 2) {
        result.physical_status = "rejected";
        result.failure_reason = "collapsed_phase_set";
        return result;
    }

    std::vector<double> modified_balance(modified_feed.size(), 0.0);
    std::vector<double> ordinary_balance(physical_feed.size(), 0.0);
    result.minimum_phase_distance = std::numeric_limits<double>::infinity();
    for (const Held2StageIIIPhase& phase : result.phases) {
        for (std::size_t index = 0; index < modified_balance.size(); ++index) {
            modified_balance[index] +=
                phase.phase_fraction * phase.modified_fractions[index];
        }
        for (std::size_t index = 0; index < ordinary_balance.size(); ++index) {
            ordinary_balance[index] +=
                phase.phase_fraction * phase.physical_fractions[index];
        }
        result.phase_charge_inf_norm = std::max(
            result.phase_charge_inf_norm,
            std::abs(charge_residual(
                coordinates.charges, phase.physical_fractions
            ))
        );
        double phase_scale = 0.0;
        for (std::size_t index = 0; index < coordinates.charges.size(); ++index) {
            phase_scale += std::abs(
                coordinates.charges[index] * phase.physical_fractions[index]
            );
        }
        result.phase_charge_scale = std::max(
            result.phase_charge_scale, phase_scale
        );
    }
    result.modified_balance_inf_norm =
        maximum_abs_difference(modified_balance, modified_feed);
    result.ordinary_balance_inf_norm =
        maximum_abs_difference(ordinary_balance, physical_feed);
    for (const Held2StateEvaluation& state : accepted_states) {
        result.pressure_stationarity_inf_norm = std::max(
            result.pressure_stationarity_inf_norm,
            std::abs(state.pressure_stationarity_relative)
        );
    }
    for (std::size_t left = 0; left < result.phases.size(); ++left) {
        for (std::size_t right = left + 1;
             right < result.phases.size();
             ++right) {
            result.minimum_phase_distance = std::min(
                result.minimum_phase_distance,
                std::max(
                    maximum_abs_difference(
                        result.phases[left].modified_fractions,
                        result.phases[right].modified_fractions
                    ),
                    std::abs(
                        std::log(result.phases[left].volume)
                            - std::log(result.phases[right].volume)
                    )
                )
            );
        }
    }
    for (std::size_t retained = 0; retained < modified_feed.size(); ++retained) {
        const double lower_bound =
            kHeld2ModifiedLowerScale * coordinates.modified_factors[retained];
        bool trace = false;
        for (const Held2StageIIIPhase& phase : result.phases) {
            trace = trace
                || phase.modified_fractions[retained] <= 10.0 * lower_bound;
        }
        if (trace) {
            ++result.trace_component_count;
            continue;
        }
        ++result.certified_modified_potential_count;
        for (std::size_t left = 0; left < accepted_states.size(); ++left) {
            for (std::size_t right = left + 1;
                 right < accepted_states.size();
                 ++right) {
                const double gap = std::abs(
                    accepted_states[left].modified_potentials[retained]
                    - accepted_states[right].modified_potentials[retained]
                );
                const double scale = std::max(
                    std::abs(accepted_states[left].modified_potentials[retained]),
                    std::abs(accepted_states[right].modified_potentials[retained])
                );
                if (gap > result.modified_potential_mixed_gap) {
                    result.modified_potential_mixed_gap = gap;
                    result.modified_potential_scale = scale;
                }
            }
        }
    }
    result.physical_evidence_available = true;
    result.phase_identity_evidence_available = true;
    if (result.trace_component_count > 0) {
        result.trace_refinement_status =
            "complementarity_refinement_required";
        result.failure_reason = "trace_component_requires_log_refinement";
        return result;
    }
    result.trace_refinement_status = "not_required";
    const bool duplicate_identity = audit_held2_tolerance(
        kHeld2PhaseMerge, result.minimum_phase_distance
    ).passed;
    const bool distinct_identity = audit_held2_tolerance(
        kHeld2PhaseDistinct, result.minimum_phase_distance
    ).passed;
    if (duplicate_identity) {
        result.phase_identity_status = "duplicate";
    } else if (distinct_identity) {
        result.phase_identity_status = "confidently_distinct";
    } else {
        result.phase_identity_status = "unresolved";
        result.physical_status = "rejected";
        result.failure_reason = "stage_iii_phase_identity_unresolved";
        return result;
    }
    const bool physical =
        audit_held2_tolerance(
            kHeld2Stage3ModifiedBalance,
            result.modified_balance_inf_norm
        ).passed
        && audit_held2_tolerance(
            kHeld2Stage3ExplicitBalance,
            result.ordinary_balance_inf_norm
        ).passed
        && audit_held2_tolerance(
            kHeld2Stage3Charge,
            result.phase_charge_inf_norm,
            result.phase_charge_scale
        ).passed
        && audit_held2_tolerance(
            kHeld2Stage3Pressure,
            result.pressure_stationarity_inf_norm
        ).passed
        && audit_held2_tolerance(
            kHeld2Stage3Potential,
            result.modified_potential_mixed_gap,
            result.modified_potential_scale
        ).passed
        && distinct_identity
        && audit_held2_tolerance(
            kHeld2PhaseActivity,
            result.minimum_phase_fraction
        ).passed
        && audit_held2_tolerance(
            kHeld2Stage3Stationarity,
            result.kkt_stationarity_inf_norm
        ).passed
        && audit_held2_tolerance(
            kHeld2Stage3DualSign,
            result.dual_sign_violation_inf_norm
        ).passed
        && audit_held2_tolerance(
            kHeld2Stage3Complementarity,
            result.bound_complementarity_inf_norm
        ).passed;
    if (!physical) {
        result.physical_status = "rejected";
        result.failure_reason = "stage_iii_physical_certificate_failed";
        Held2ProgressEvent physical_progress;
        physical_progress.kind = Held2ProgressKind::Certificate;
        physical_progress.stage = "STAGE III PHYSICAL";
        physical_progress.status = "failed";
        physical_progress.reason = result.failure_reason;
        physical_progress.primal_residual = std::max(
            result.modified_balance_inf_norm,
            result.ordinary_balance_inf_norm
        );
        physical_progress.dual_residual = result.modified_potential_mixed_gap;
        physical_progress.complementarity =
            result.bound_complementarity_inf_norm;
        observe_held2(observer, physical_progress);
        return result;
    }
    result.physical_status = "accepted";
    result.feedback = "none";
    Held2ProgressEvent physical_progress;
    physical_progress.kind = Held2ProgressKind::Certificate;
    physical_progress.stage = "STAGE III PHYSICAL";
    physical_progress.status = "passed";
    physical_progress.primal_residual = std::max(
        result.modified_balance_inf_norm,
        result.ordinary_balance_inf_norm
    );
    physical_progress.dual_residual = result.modified_potential_mixed_gap;
    physical_progress.complementarity = result.bound_complementarity_inf_norm;
    observe_held2(observer, physical_progress);
    return result;
}
}  // namespace epcsaft_equilibrium
