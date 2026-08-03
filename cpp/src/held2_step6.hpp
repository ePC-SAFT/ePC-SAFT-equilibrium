#pragma once

#include "held2_step5.hpp"

namespace epcsaft_equilibrium {

struct Held2ThermodynamicAccessPolicy {
    bool packing_fraction_uses_provider = true;
};

struct Held2CandidateDecision {
    std::uint64_t insertion_id = 0;
    bool gap_passed = false;
    bool derivative_passed = false;
    bool pairwise_distinct = false;
    bool retained = false;
    std::string reason;
};

struct Held2Step6Result {
    std::string status = "indeterminate";
    std::string reason = "not_run";
    std::vector<Held2MPoint> candidates;
    std::vector<Held2CandidateDecision> decisions;
    Held2StepTiming timing;
};

[[nodiscard]] Held2Step6Result run_held2_step6(
    const Held2Step1Result& step1,
    const Held2Step4Result& step4,
    const Held2PersistentState& state,
    const Held2PackingFractionEvaluator& packing_fraction,
    Held2ProgressObserver* observer = nullptr,
    Held2ThermodynamicAccessPolicy access_policy = {}
);

}  // namespace epcsaft_equilibrium
