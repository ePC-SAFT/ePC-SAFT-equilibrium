#pragma once

#include <string>
#include <vector>

#include "held2.hpp"
#include "held2_stage_i_direct.hpp"

namespace epcsaft_equilibrium {

struct Held2WorkflowTransition {
    int completed_step = 0;
    Held2NextAction next_action = Held2NextAction::TerminateIndeterminate;
    std::string outcome;
    std::string reason;
    int feedback_owner_step = 0;
};

struct Held2WorkflowState {
    std::string outcome = "running";
    std::string failure_stage;
    std::string failure_reason;
    int completed_step = 0;
    int feedback_owner_step = 0;
    Held2NextAction next_action = Held2NextAction::TerminateIndeterminate;
    std::vector<Held2WorkflowTransition> transitions;
};

class Held2WorkflowController {
public:
    [[nodiscard]] Held2NextAction complete_reference(
        const Held2PressureEnvelopeResult& reference
    );

    [[nodiscard]] Held2NextAction complete_stage_i(
        const Held2StageIDirectResult& stage_i
    );

    [[nodiscard]] Held2NextAction complete_stage_ii(
        const Held2StageIIResult& stage_ii
    );

    [[nodiscard]] Held2NextAction complete_stage_iii(
        const Held2StageIIIResult& stage_iii
    );

    [[nodiscard]] const Held2WorkflowState& state() const noexcept {
        return state_;
    }

private:
    void require_step(int completed_step) const;
    Held2NextAction transition(
        int completed_step,
        Held2NextAction next_action,
        std::string outcome,
        std::string reason = {},
        int feedback_owner_step = 0
    );

    Held2WorkflowState state_;
};

}  // namespace epcsaft_equilibrium
