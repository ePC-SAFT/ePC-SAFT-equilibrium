#include "held2_step7.hpp"

namespace epcsaft_equilibrium {

Held2Step7Result run_held2_step7(
    Held2PersistentState& state,
    const Held2Step5Result& step5,
    const Held2Step6Result& step6,
    const Held2ResourceProfile& resources,
    Held2ProgressObserver*
) {
    Held2Step7Result result;
    result.timing.invocation_count = 1;
    if (step5.status != "complete" || step6.status != "complete") {
        result.reason = "invalid_step7_input";
        return result;
    }
    if (step6.candidates.size() >= 2) {
        result.status = "complete";
        result.reason = "stage_iii_ready";
        result.next_step = 8;
    } else {
        ++state.major_iteration;
        if (state.major_iteration >= resources.step7_major_iteration_cap) {
            result.reason = "stage_ii_major_iteration_limit";
        } else if (
            state.starts_consumed_in_epoch
                >= resources.step5_start_epoch_size
            && !state.start_epoch_added_member
        ) {
            result.reason = "stage_ii_stagnation";
        } else {
            result.status = "complete";
            result.reason = "continue_stage_ii";
            result.next_step = 4;
        }
    }
    result.timing.terminal_status = result.status;
    result.timing.terminal_reason = result.reason;
    result.timing.next_step = result.next_step.value_or(0);
    return result;
}

}  // namespace epcsaft_equilibrium
