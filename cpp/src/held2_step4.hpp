#pragma once

#include "held2_step3.hpp"

#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace epcsaft_equilibrium {

struct Held2LpCertificate {
    bool primal_feasible = false;
    bool dual_feasible = false;
    double primal_residual_inf = std::numeric_limits<double>::infinity();
    double dual_residual_inf = std::numeric_limits<double>::infinity();
    double complementarity_inf = std::numeric_limits<double>::infinity();
};

struct Held2Step4CutEvidence {
    int id = -1;
    double intercept = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> slopes;
};

struct Held2Step4Result {
    std::string status = "indeterminate";
    std::string reason = "not_run";
    int major_iteration = -1;
    int upper_solve_count = -1;
    std::optional<double> upper_bound;
    std::optional<std::vector<double>> multipliers;
    std::vector<int> active_cut_ids;
    std::vector<Held2Step4CutEvidence> cut_snapshot;
    std::optional<Held2Coordinates> coordinate_snapshot;
    std::optional<Held2LpCertificate> certificate;
    Held2StepTiming timing;
};

[[nodiscard]] Held2Step4Result run_held2_step4(
    Held2PersistentState& state,
    Held2ProgressObserver* observer = nullptr
);

}  // namespace epcsaft_equilibrium
