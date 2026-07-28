#pragma once

#include "held2_step6.hpp"

namespace epcsaft_equilibrium {

struct Held2Step7Result {
    std::string status = "indeterminate";
    std::string reason = "not_run";
    std::optional<int> next_step;
    Held2StepTiming timing;
};

[[nodiscard]] Held2Step7Result run_held2_step7(
    Held2PersistentState& state,
    const Held2Step5Result& step5,
    const Held2Step6Result& step6,
    const Held2ResourceProfile& resources,
    Held2ProgressObserver* observer = nullptr
);

}  // namespace epcsaft_equilibrium
