#include "held2_step4.hpp"
#include "held2_tolerances.hpp"

#include <Highs.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace epcsaft_equilibrium {
namespace {

struct UpperCut {
    int id = -1;
    double intercept = 0.0;
    std::vector<double> slopes;
};

struct UpperResult {
    bool primal_feasible = false;
    bool dual_feasible = false;
    double upper_bound = 0.0;
    std::vector<double> multipliers;
    std::vector<double> cut_slacks;
    std::vector<int> active_cut_ids;
    double primal_residual_inf = std::numeric_limits<double>::infinity();
    double dual_residual_inf = std::numeric_limits<double>::infinity();
    double complementarity_inf = std::numeric_limits<double>::infinity();
};

void validate_cuts(const std::vector<UpperCut>& cuts, std::size_t dimension) {
    for (const UpperCut& cut : cuts) {
        if (!std::isfinite(cut.intercept) || cut.slopes.size() != dimension
            || !std::all_of(cut.slopes.begin(), cut.slopes.end(), [](double value) {
                return std::isfinite(value);
            })) {
            throw std::invalid_argument("HELD2 upper LP cut data are invalid");
        }
    }
}

UpperResult solve_upper_highs(
    const std::vector<UpperCut>& cuts,
    std::size_t size
) {
    validate_cuts(cuts, size);
    const HighsInt dimension = static_cast<HighsInt>(size);
    const HighsInt value_column = dimension;

    HighsModel model;
    model.lp_.num_col_ = dimension + 1;
    model.lp_.num_row_ = static_cast<HighsInt>(cuts.size());
    model.lp_.sense_ = ObjSense::kMaximize;
    model.lp_.col_cost_.assign(static_cast<std::size_t>(dimension + 1), 0.0);
    model.lp_.col_cost_[static_cast<std::size_t>(value_column)] = 1.0;
    model.lp_.col_lower_.assign(
        static_cast<std::size_t>(dimension + 1), -kHighsInf
    );
    model.lp_.col_upper_.assign(
        static_cast<std::size_t>(dimension + 1), kHighsInf
    );
    model.lp_.row_lower_.assign(cuts.size(), -kHighsInf);
    model.lp_.row_upper_.reserve(cuts.size());
    model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
    model.lp_.a_matrix_.start_.reserve(static_cast<std::size_t>(dimension + 2));
    model.lp_.a_matrix_.start_ = {0};
    for (const UpperCut& cut : cuts) {
        model.lp_.row_upper_.push_back(cut.intercept);
    }
    for (HighsInt index = 0; index < dimension; ++index) {
        for (std::size_t row = 0; row < cuts.size(); ++row) {
            const double coefficient =
                -cuts[row].slopes[static_cast<std::size_t>(index)];
            if (coefficient != 0.0) {
                model.lp_.a_matrix_.index_.push_back(static_cast<HighsInt>(row));
                model.lp_.a_matrix_.value_.push_back(coefficient);
            }
        }
        model.lp_.a_matrix_.start_.push_back(
            static_cast<HighsInt>(model.lp_.a_matrix_.index_.size())
        );
    }
    for (std::size_t row = 0; row < cuts.size(); ++row) {
        model.lp_.a_matrix_.index_.push_back(static_cast<HighsInt>(row));
        model.lp_.a_matrix_.value_.push_back(1.0);
    }
    model.lp_.a_matrix_.start_.push_back(
        static_cast<HighsInt>(model.lp_.a_matrix_.index_.size())
    );

    UpperResult result;
    Highs highs;
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
        return result;
    }
    const HighsStatus run_status = highs.run();
    const HighsModelStatus model_status = highs.getModelStatus();
    if (run_status == HighsStatus::kError || model_status != HighsModelStatus::kOptimal) {
        return result;
    }

    const HighsInfo& info = highs.getInfo();
    const HighsSolution& solution = highs.getSolution();
    if (info.primal_solution_status != kSolutionStatusFeasible
        || info.dual_solution_status != kSolutionStatusFeasible
        || solution.col_value.size() != static_cast<std::size_t>(dimension + 1)
        || solution.col_dual.size() != static_cast<std::size_t>(dimension + 1)
        || solution.row_dual.size() != cuts.size()) {
        return result;
    }

    result.upper_bound = solution.col_value[static_cast<std::size_t>(value_column)];
    result.multipliers.assign(solution.col_value.begin(), solution.col_value.begin() + dimension);
    for (double& multiplier : result.multipliers) {
        if (multiplier == 0.0) {
            multiplier = 0.0;
        }
    }
    result.cut_slacks.reserve(cuts.size());
    double primal_residual = 0.0;
    double primal_scale = std::abs(result.upper_bound);
    for (std::size_t row = 0; row < cuts.size(); ++row) {
        const UpperCut& cut = cuts[row];
        double cut_value = cut.intercept;
        for (std::size_t index = 0; index < result.multipliers.size(); ++index) {
            cut_value += cut.slopes[index] * result.multipliers[index];
        }
        const double slack = cut_value - result.upper_bound;
        primal_scale = std::max({
            primal_scale, std::abs(cut_value), std::abs(cut.intercept)
        });
        result.cut_slacks.push_back(slack);
        primal_residual = std::max(primal_residual, std::max(0.0, -slack));
        if (audit_held2_tolerance(kHeld2LpActiveCut, slack).passed) {
            result.active_cut_ids.push_back(cut.id);
        }
    }
    result.primal_residual_inf = primal_residual;

    double dual_residual = 0.0;
    double dual_scale = 1.0;
    double complementarity = 0.0;
    for (std::size_t row = 0; row < cuts.size(); ++row) {
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
        for (std::size_t row = 0; row < cuts.size(); ++row) {
            const double coefficient = column == value_column
                ? 1.0
                : -cuts[row].slopes[static_cast<std::size_t>(column)];
            residual -= coefficient * solution.row_dual[row];
        }
        dual_residual = std::max(dual_residual, std::abs(residual));

        const double reduced_cost =
            solution.col_dual[static_cast<std::size_t>(column)];
        dual_scale = std::max({dual_scale, std::abs(residual), std::abs(reduced_cost)});
        dual_residual = std::max(dual_residual, std::abs(reduced_cost));
    }
    result.dual_residual_inf = dual_residual;
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
        result.multipliers.clear();
    }
    return result;
}

}  // namespace

Held2Step4Result run_held2_step4(
    Held2PersistentState& state,
    Held2ProgressObserver* observer
) {
    Held2Step4Result result;
    result.timing.invocation_count = 1;
    const std::size_t dimension = state.feed.size();
    std::vector<UpperCut> cuts;
    for (const Held2MPoint& point : state.M) {
        std::vector<double> slopes(dimension);
        for (std::size_t index = 0; index < dimension; ++index) {
            slopes[index] =
                state.feed[index] - point.independent_modified_fractions[index];
        }
        cuts.push_back({
            static_cast<int>(point.insertion_id),
            point.reduced_gibbs,
            std::move(slopes),
        });
    }
    cuts.push_back({
        -1, state.feed_reduced_gibbs, std::vector<double>(dimension, 0.0),
    });
    const UpperResult upper = solve_upper_highs(cuts, dimension);
    result.certificate = Held2LpCertificate{
        upper.primal_feasible,
        upper.dual_feasible,
        upper.primal_residual_inf,
        upper.dual_residual_inf,
        upper.complementarity_inf,
    };
    if (!upper.primal_feasible || !upper.dual_feasible) {
        result.reason = "upper_lp_uncertified";
        result.timing.terminal_status = "indeterminate";
        result.timing.terminal_reason = result.reason;
        return result;
    }
    state.upper_bound = upper.upper_bound;
    state.multipliers = upper.multipliers;
    ++state.upper_solve_count;
    result.status = "complete";
    result.reason = "step4_complete";
    result.upper_bound = upper.upper_bound;
    result.multipliers = upper.multipliers;
    result.active_cut_ids = upper.active_cut_ids;
    result.timing.optimizer_solves = 1;
    result.timing.optimizer_iterations = 1;
    result.timing.terminal_status = result.status;
    result.timing.terminal_reason = result.reason;
    result.timing.next_step = 5;
    Held2ProgressEvent progress;
    progress.kind = Held2ProgressKind::StageIIUpper;
    progress.major_iteration = state.major_iteration;
    progress.upper_bound = upper.upper_bound;
    progress.primal_residual = upper.primal_residual_inf;
    progress.dual_residual = upper.dual_residual_inf;
    progress.count = static_cast<int>(state.M.size());
    progress.status = "certified";
    observe_held2(observer, progress);
    return result;
}

}  // namespace epcsaft_equilibrium
