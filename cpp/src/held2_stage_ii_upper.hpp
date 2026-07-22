#pragma once

#include <limits>
#include <string>
#include <vector>

namespace epcsaft_equilibrium {

struct Held2StageIIUpperCut {
    int id = -1;
    double intercept = 0.0;
    std::vector<double> slopes;
};

struct Held2StageIIUpperProblem {
    std::vector<Held2StageIIUpperCut> cuts;
    std::vector<double> multiplier_lower_bounds;
    std::vector<double> multiplier_upper_bounds;
    double value_lower_bound = -std::numeric_limits<double>::infinity();
    double value_upper_bound = std::numeric_limits<double>::infinity();
};

struct Held2StageIIUpperResult {
    std::string outcome = "indeterminate";
    std::string solver = "highs_lp";
    std::string solver_status;
    std::string solver_version;
    bool solver_finished = false;
    bool primal_feasible = false;
    bool dual_feasible = false;
    double upper_bound = 0.0;
    std::vector<double> multipliers;
    std::vector<double> cut_slacks;
    std::vector<double> cut_duals;
    std::vector<int> active_cut_ids;
    double primal_residual_inf = std::numeric_limits<double>::infinity();
    double dual_residual_inf = std::numeric_limits<double>::infinity();
};

[[nodiscard]] Held2StageIIUpperResult solve_held2_stage_ii_upper_highs(
    const Held2StageIIUpperProblem& problem
);

[[nodiscard]] Held2StageIIUpperResult solve_held2_stage_ii_upper_analytic_1d(
    const Held2StageIIUpperProblem& problem
);

}  // namespace epcsaft_equilibrium
