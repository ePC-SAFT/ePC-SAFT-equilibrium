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

struct Held2Step4Result {
    std::string status = "indeterminate";
    std::string reason = "not_run";
    std::optional<double> upper_bound;
    std::optional<std::vector<double>> multipliers;
    std::vector<int> active_cut_ids;
    std::optional<Held2LpCertificate> certificate;
    Held2StepTiming timing;
};

[[nodiscard]] Held2Step4Result run_held2_step4(
    Held2PersistentState& state,
    Held2ProgressObserver* observer = nullptr
);

}  // namespace epcsaft_equilibrium
