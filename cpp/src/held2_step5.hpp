#pragma once

#include <functional>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "held2_step4.hpp"

namespace epcsaft_equilibrium {

using Held2PackingFractionEvaluator = std::function<double(
    const std::vector<double>&,
    double
)>;

struct Held2ResourceProfile {
    int step2_search_budget = 50;
    int step5_start_cap = 128;
    int step7_major_iteration_cap = 80;
};

struct Held2LocalCertificate {
    std::uint64_t start_ordinal = 0;
    std::string solver_status;
    bool finite_and_in_domain = false;
    bool accepted = false;
};

struct Held2Step5Result {
    std::string status = "indeterminate";
    std::string reason = "not_run";
    std::optional<double> lower_value;
    std::optional<Held2MPoint> terminal;
    std::uint64_t starts_consumed = 0;
    std::vector<Held2LocalCertificate> attempts;
    Held2StepTiming timing;
};

[[nodiscard]] Held2Step5Result run_held2_step5(
    const Held2Step1Result& step1,
    const Held2Step4Result& step4,
    Held2PersistentState& state,
    const Held2StateEvaluator& evaluator,
    const Held2PackingFractionEvaluator& packing_fraction,
    const Held2ResourceProfile& resources,
    Held2ProgressObserver* observer = nullptr
);

}  // namespace epcsaft_equilibrium
