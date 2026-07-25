#pragma once

#include "held2_step3.hpp"

#include <limits>
#include <optional>
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
    double primal_scale = 0.0;
    double dual_residual_inf = std::numeric_limits<double>::infinity();
    double dual_scale = 0.0;
    double complementarity_inf = std::numeric_limits<double>::infinity();
};

struct Held2LpCertificate {
    bool primal_feasible = false;
    bool dual_feasible = false;
    double primal_residual_inf = std::numeric_limits<double>::infinity();
    double dual_residual_inf = std::numeric_limits<double>::infinity();
    double complementarity_inf = std::numeric_limits<double>::infinity();
};

struct Held2Step4Result {
    std::string status = "indeterminate";
    std::string reason = "not_run";
    std::optional<double> upper_bound;
    std::optional<std::vector<double>> multipliers;
    std::vector<int> active_cut_ids;
    std::optional<Held2LpCertificate> certificate;
    Held2StepTiming timing;
};

[[nodiscard]] Held2StageIIUpperResult solve_held2_stage_ii_upper_highs(
    const Held2StageIIUpperProblem& problem
);

[[nodiscard]] Held2StageIIUpperResult solve_held2_stage_ii_upper_analytic_1d(
    const Held2StageIIUpperProblem& problem
);

[[nodiscard]] Held2Step4Result run_held2_step4(
    Held2PersistentState& state,
    Held2ProgressObserver* observer = nullptr
);

}  // namespace epcsaft_equilibrium
