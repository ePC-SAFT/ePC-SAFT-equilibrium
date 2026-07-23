#include "held2_stage_ii_upper.hpp"
#include "held2_tolerances.hpp"

#include <Highs.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace epcsaft_equilibrium {
namespace {

void validate_problem(const Held2StageIIUpperProblem& problem) {
    const std::size_t dimension = problem.multiplier_lower_bounds.size();
    if (problem.multiplier_upper_bounds.size() != dimension) {
        throw std::invalid_argument("HELD2 upper LP multiplier bounds have different sizes");
    }
    for (std::size_t index = 0; index < dimension; ++index) {
        const double lower = problem.multiplier_lower_bounds[index];
        const double upper = problem.multiplier_upper_bounds[index];
        if (std::isnan(lower) || std::isnan(upper) || lower > upper) {
            throw std::invalid_argument("HELD2 upper LP multiplier bounds are invalid");
        }
    }
    if (std::isnan(problem.value_lower_bound)
        || std::isnan(problem.value_upper_bound)
        || problem.value_lower_bound > problem.value_upper_bound) {
        throw std::invalid_argument("HELD2 upper LP value bounds are invalid");
    }
    for (const Held2StageIIUpperCut& cut : problem.cuts) {
        if (!std::isfinite(cut.intercept) || cut.slopes.size() != dimension
            || !std::all_of(cut.slopes.begin(), cut.slopes.end(), [](double value) {
                return std::isfinite(value);
            })) {
            throw std::invalid_argument("HELD2 upper LP cut data are invalid");
        }
    }
}

double to_highs_bound(double value) {
    if (value == std::numeric_limits<double>::infinity()) {
        return kHighsInf;
    }
    if (value == -std::numeric_limits<double>::infinity()) {
        return -kHighsInf;
    }
    return value;
}

}  // namespace

Held2StageIIUpperResult solve_held2_stage_ii_upper_highs(
    const Held2StageIIUpperProblem& problem
) {
    validate_problem(problem);
    const HighsInt dimension = static_cast<HighsInt>(
        problem.multiplier_lower_bounds.size()
    );
    const HighsInt value_column = dimension;

    HighsModel model;
    model.lp_.num_col_ = dimension + 1;
    model.lp_.num_row_ = static_cast<HighsInt>(problem.cuts.size());
    model.lp_.sense_ = ObjSense::kMaximize;
    model.lp_.col_cost_.assign(static_cast<std::size_t>(dimension + 1), 0.0);
    model.lp_.col_cost_[static_cast<std::size_t>(value_column)] = 1.0;
    model.lp_.col_lower_.reserve(static_cast<std::size_t>(dimension + 1));
    model.lp_.col_upper_.reserve(static_cast<std::size_t>(dimension + 1));
    for (HighsInt index = 0; index < dimension; ++index) {
        model.lp_.col_lower_.push_back(
            to_highs_bound(problem.multiplier_lower_bounds[static_cast<std::size_t>(index)])
        );
        model.lp_.col_upper_.push_back(
            to_highs_bound(problem.multiplier_upper_bounds[static_cast<std::size_t>(index)])
        );
    }
    model.lp_.col_lower_.push_back(to_highs_bound(problem.value_lower_bound));
    model.lp_.col_upper_.push_back(to_highs_bound(problem.value_upper_bound));
    model.lp_.row_lower_.assign(problem.cuts.size(), -kHighsInf);
    model.lp_.row_upper_.reserve(problem.cuts.size());
    model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
    model.lp_.a_matrix_.start_.reserve(static_cast<std::size_t>(dimension + 2));
    model.lp_.a_matrix_.start_ = {0};
    for (const Held2StageIIUpperCut& cut : problem.cuts) {
        model.lp_.row_upper_.push_back(cut.intercept);
    }
    for (HighsInt index = 0; index < dimension; ++index) {
        for (std::size_t row = 0; row < problem.cuts.size(); ++row) {
            const double coefficient =
                -problem.cuts[row].slopes[static_cast<std::size_t>(index)];
            if (coefficient != 0.0) {
                model.lp_.a_matrix_.index_.push_back(static_cast<HighsInt>(row));
                model.lp_.a_matrix_.value_.push_back(coefficient);
            }
        }
        model.lp_.a_matrix_.start_.push_back(
            static_cast<HighsInt>(model.lp_.a_matrix_.index_.size())
        );
    }
    for (std::size_t row = 0; row < problem.cuts.size(); ++row) {
        model.lp_.a_matrix_.index_.push_back(static_cast<HighsInt>(row));
        model.lp_.a_matrix_.value_.push_back(1.0);
    }
    model.lp_.a_matrix_.start_.push_back(
        static_cast<HighsInt>(model.lp_.a_matrix_.index_.size())
    );

    Held2StageIIUpperResult result;
    Highs highs;
    result.solver_version = highs.version();
    if (highs.setOptionValue("output_flag", false) == HighsStatus::kError
        || highs.setOptionValue("threads", 1) == HighsStatus::kError
        || highs.setOptionValue(
               "primal_feasibility_tolerance",
               kHeld2IpoptTarget.atol
           ) == HighsStatus::kError
        || highs.setOptionValue(
               "dual_feasibility_tolerance",
               kHeld2IpoptTarget.atol
           ) == HighsStatus::kError
        || highs.passModel(model) == HighsStatus::kError) {
        result.solver_status = "setup_failed";
        return result;
    }
    const HighsStatus run_status = highs.run();
    const HighsModelStatus model_status = highs.getModelStatus();
    result.solver_status = highs.modelStatusToString(model_status);
    result.solver_finished = run_status != HighsStatus::kError;
    if (model_status == HighsModelStatus::kInfeasible) {
        result.outcome = "infeasible";
        return result;
    }
    if (model_status == HighsModelStatus::kUnbounded) {
        result.outcome = "unbounded";
        return result;
    }
    if (model_status == HighsModelStatus::kUnboundedOrInfeasible) {
        result.solver_status = "UnboundedOrInfeasible";
        return result;
    }
    if (run_status == HighsStatus::kError || model_status != HighsModelStatus::kOptimal) {
        return result;
    }

    const HighsInfo& info = highs.getInfo();
    const HighsSolution& solution = highs.getSolution();
    if (info.primal_solution_status != kSolutionStatusFeasible
        || info.dual_solution_status != kSolutionStatusFeasible
        || solution.col_value.size() != static_cast<std::size_t>(dimension + 1)
        || solution.col_dual.size() != static_cast<std::size_t>(dimension + 1)
        || solution.row_dual.size() != problem.cuts.size()) {
        result.solver_status = "optimal_without_complete_solution";
        return result;
    }

    result.upper_bound = solution.col_value[static_cast<std::size_t>(value_column)];
    result.multipliers.assign(solution.col_value.begin(), solution.col_value.begin() + dimension);
    for (double& multiplier : result.multipliers) {
        if (multiplier == 0.0) {
            multiplier = 0.0;
        }
    }
    result.cut_slacks.reserve(problem.cuts.size());
    result.cut_duals.reserve(problem.cuts.size());
    double primal_residual = 0.0;
    double primal_scale = std::abs(result.upper_bound);
    for (std::size_t row = 0; row < problem.cuts.size(); ++row) {
        const Held2StageIIUpperCut& cut = problem.cuts[row];
        double cut_value = cut.intercept;
        for (std::size_t index = 0; index < result.multipliers.size(); ++index) {
            cut_value += cut.slopes[index] * result.multipliers[index];
        }
        const double slack = cut_value - result.upper_bound;
        primal_scale = std::max({
            primal_scale, std::abs(cut_value), std::abs(cut.intercept)
        });
        result.cut_slacks.push_back(slack);
        result.cut_duals.push_back(solution.row_dual[row]);
        primal_residual = std::max(primal_residual, std::max(0.0, -slack));
        if (audit_held2_tolerance(kHeld2LpActiveCut, slack).passed) {
            result.active_cut_ids.push_back(cut.id);
        }
    }
    for (std::size_t index = 0; index < result.multipliers.size(); ++index) {
        const double value = result.multipliers[index];
        primal_residual = std::max(
            primal_residual,
            std::max(0.0, problem.multiplier_lower_bounds[index] - value)
        );
        primal_residual = std::max(
            primal_residual,
            std::max(0.0, value - problem.multiplier_upper_bounds[index])
        );
    }
    primal_residual = std::max(
        primal_residual,
        std::max(0.0, problem.value_lower_bound - result.upper_bound)
    );
    primal_residual = std::max(
        primal_residual,
        std::max(0.0, result.upper_bound - problem.value_upper_bound)
    );
    result.primal_residual_inf = primal_residual;
    result.primal_scale = primal_scale;

    double dual_residual = 0.0;
    double dual_scale = 1.0;
    double complementarity = 0.0;
    for (std::size_t row = 0; row < problem.cuts.size(); ++row) {
        dual_residual = std::max(
            dual_residual,
            std::max(0.0, -solution.row_dual[row])
        );
        dual_scale = std::max(dual_scale, std::abs(solution.row_dual[row]));
        complementarity = std::max(
            complementarity,
            std::abs(solution.row_dual[row] * result.cut_slacks[row])
        );
    }
    for (HighsInt column = 0; column < dimension + 1; ++column) {
        double residual = model.lp_.col_cost_[static_cast<std::size_t>(column)]
            - solution.col_dual[static_cast<std::size_t>(column)];
        for (std::size_t row = 0; row < problem.cuts.size(); ++row) {
            const double coefficient = column == value_column
                ? 1.0
                : -problem.cuts[row].slopes[static_cast<std::size_t>(column)];
            residual -= coefficient * solution.row_dual[row];
        }
        dual_residual = std::max(dual_residual, std::abs(residual));

        const std::size_t index = static_cast<std::size_t>(column);
        const double value = solution.col_value[index];
        const double lower = model.lp_.col_lower_[index];
        const double upper = model.lp_.col_upper_[index];
        const double reduced_cost = solution.col_dual[index];
        dual_scale = std::max({dual_scale, std::abs(residual), std::abs(reduced_cost)});
        const bool at_lower = std::isfinite(lower)
            && audit_held2_tolerance(kHeld2BoundActivity, value - lower).passed;
        const bool at_upper = std::isfinite(upper)
            && audit_held2_tolerance(kHeld2BoundActivity, upper - value).passed;
        if (at_lower && !at_upper) {
            dual_residual = std::max(dual_residual, std::max(0.0, reduced_cost));
            complementarity = std::max(
                complementarity,
                std::abs((value - lower) * reduced_cost)
            );
        } else if (at_upper && !at_lower) {
            dual_residual = std::max(dual_residual, std::max(0.0, -reduced_cost));
            complementarity = std::max(
                complementarity,
                std::abs((upper - value) * reduced_cost)
            );
        } else if (!at_lower && !at_upper) {
            dual_residual = std::max(dual_residual, std::abs(reduced_cost));
        }
    }
    result.dual_residual_inf = dual_residual;
    result.dual_scale = dual_scale;
    result.complementarity_inf = complementarity;
    result.primal_feasible = audit_held2_tolerance(
        kHeld2LpPrimal, primal_residual, primal_scale
    ).passed;
    result.dual_feasible = audit_held2_tolerance(
        kHeld2LpDual, dual_residual, dual_scale
    ).passed && audit_held2_tolerance(
        kHeld2LpComplementarity, complementarity
    ).passed;
    if (!result.primal_feasible || !result.dual_feasible) {
        result.solver_status = "optimal_but_residual_audit_failed";
        result.multipliers.clear();
        return result;
    }
    result.outcome = "optimal";
    return result;
}

Held2StageIIUpperResult solve_held2_stage_ii_upper_analytic_1d(
    const Held2StageIIUpperProblem& problem
) {
    validate_problem(problem);
    if (problem.multiplier_lower_bounds.size() != 1 || problem.cuts.empty()) {
        throw std::invalid_argument("analytic HELD2 upper oracle requires one dimension and cuts");
    }
    std::vector<double> candidates;
    const double lower = problem.multiplier_lower_bounds.front();
    const double upper = problem.multiplier_upper_bounds.front();
    if (std::isfinite(lower)) {
        candidates.push_back(lower);
    }
    if (std::isfinite(upper)) {
        candidates.push_back(upper);
    }
    for (std::size_t left = 0; left < problem.cuts.size(); ++left) {
        for (std::size_t right = left + 1; right < problem.cuts.size(); ++right) {
            const double slope_delta = problem.cuts[left].slopes.front()
                - problem.cuts[right].slopes.front();
            if (std::abs(slope_delta) <= std::numeric_limits<double>::epsilon()) {
                continue;
            }
            const double candidate = (
                problem.cuts[right].intercept - problem.cuts[left].intercept
            ) / slope_delta;
            if (candidate >= lower && candidate <= upper) {
                candidates.push_back(candidate);
            }
        }
    }
    Held2StageIIUpperResult result;
    result.solver = "analytic_1d_test_oracle";
    result.solver_status = "evaluated";
    result.solver_finished = true;
    double best = -std::numeric_limits<double>::infinity();
    double multiplier = 0.0;
    for (double candidate : candidates) {
        double envelope = std::numeric_limits<double>::infinity();
        for (const Held2StageIIUpperCut& cut : problem.cuts) {
            envelope = std::min(
                envelope,
                cut.intercept + cut.slopes.front() * candidate
            );
        }
        if (envelope > best) {
            best = envelope;
            multiplier = candidate;
        }
    }
    if (!std::isfinite(best)) {
        result.outcome = "indeterminate";
        return result;
    }
    result.outcome = "optimal";
    result.upper_bound = best;
    result.multipliers = {multiplier};
    result.primal_feasible = true;
    result.dual_feasible = true;
    result.primal_residual_inf = 0.0;
    result.primal_scale = std::max(std::abs(best), std::abs(multiplier));
    result.dual_residual_inf = 0.0;
    result.dual_scale = 1.0;
    result.complementarity_inf = 0.0;
    return result;
}

}  // namespace epcsaft_equilibrium
