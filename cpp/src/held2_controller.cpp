#include "held2_controller.hpp"

#include <stdexcept>
#include <utility>

namespace epcsaft_equilibrium {

void Held2WorkflowController::require_step(int completed_step) const {
    if (state_.completed_step != completed_step) {
        throw std::logic_error("HELD2 workflow transition is out of order");
    }
}

Held2NextAction Held2WorkflowController::transition(
    int completed_step,
    Held2NextAction next_action,
    std::string outcome,
    std::string reason,
    int feedback_owner_step
) {
    state_.completed_step = completed_step;
    state_.next_action = next_action;
    state_.outcome = std::move(outcome);
    state_.failure_reason = std::move(reason);
    state_.feedback_owner_step = feedback_owner_step;
    state_.transitions.push_back({
        completed_step,
        next_action,
        state_.outcome,
        state_.failure_reason,
        feedback_owner_step,
    });
    return next_action;
}

Held2NextAction Held2WorkflowController::complete_reference(
    const Held2PressureEnvelopeResult& reference
) {
    require_step(0);
    if (reference.outcome != "selected" || reference.selected_root_index < 0) {
        state_.failure_stage = "reference";
        return transition(
            1,
            Held2NextAction::TerminateIndeterminate,
            "indeterminate_reference",
            reference.failure_reason.empty()
                ? reference.outcome
                : reference.failure_reason,
            1
        );
    }
    return transition(
        1,
        Held2NextAction::ContinueStageII,
        "reference_selected"
    );
}

Held2NextAction Held2WorkflowController::complete_stage_i(
    const Held2StageIDirectResult& stage_i
) {
    require_step(1);
    const bool witness_valid =
        stage_i.negative_witness_index >= 0
        && static_cast<std::size_t>(stage_i.negative_witness_index)
            < stage_i.evaluations.size()
        && stage_i.evaluations[static_cast<std::size_t>(
               stage_i.negative_witness_index
           )].pressure_envelope.selected_root_index >= 0;
    if (stage_i.outcome == "negative_witness_found"
        && witness_valid) {
        return transition(
            3,
            Held2NextAction::ContinueStageII,
            "negative_witness_found"
        );
    }
    state_.failure_stage = "stage_i";
    const bool indeterminate = stage_i.outcome == "indeterminate";
    return transition(
        3,
        Held2NextAction::TerminateIndeterminate,
        stage_i.outcome == "negative_witness_found"
            ? "indeterminate_stage_i_witness"
            : indeterminate
            ? "indeterminate_stage_i"
            : "stage_i_finite_search_without_negative_witness",
        stage_i.outcome == "negative_witness_found"
            ? "stage_i_witness_has_no_selected_pressure_root"
            : stage_i.termination_reason,
        2
    );
}

Held2NextAction Held2WorkflowController::complete_stage_ii(
    const Held2StageIIResult& stage_ii
) {
    require_step(3);
    if (stage_ii.next_action == Held2NextAction::EnterStageIII
        && stage_ii.outcome == "candidate_set"
        && stage_ii.candidates.size() >= 2) {
        return transition(
            7,
            Held2NextAction::EnterStageIII,
            "stage_ii_candidate_set"
        );
    }
    state_.failure_stage = "stage_ii";
    return transition(
        7,
        Held2NextAction::TerminateIndeterminate,
        "indeterminate_stage_ii",
        stage_ii.outcome,
        stage_ii.next_action == Held2NextAction::RunGlobalEscalation ? 5 : 7
    );
}

Held2NextAction Held2WorkflowController::complete_stage_iii(
    const Held2StageIIIResult& stage_iii
) {
    require_step(7);
    if (stage_iii.physical_status == "accepted") {
        return transition(
            10,
            Held2NextAction::AcceptMultiphase,
            "physical_equilibrium_accepted"
        );
    }
    state_.failure_stage = "stage_iii";
    if (stage_iii.feedback == "return_to_stage_ii") {
        return transition(
            10,
            Held2NextAction::ReturnStageIIIFeedback,
            "stage_iii_feedback",
            stage_iii.failure_reason,
            5
        );
    }
    return transition(
        10,
        Held2NextAction::TerminateIndeterminate,
        "indeterminate_stage_iii",
        stage_iii.failure_reason,
        8
    );
}

}  // namespace epcsaft_equilibrium
