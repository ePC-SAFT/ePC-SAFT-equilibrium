#include "held2_step9.hpp"
#include "held2_tolerances.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace epcsaft_equilibrium {

Held2Step9Result run_held2_step9(
    const Held2Step4Result& step4,
    const Held2Step8Result& step8,
    const Held2StateEvaluator& evaluator,
    Held2ProgressObserver* observer
) {
    Held2Step9Result result;
    result.timing.invocation_count = 1;
    if (step4.status != "complete" || !step4.upper_bound
        || step8.outcome != Held2Step8Outcome::CertifiedFeasible
        || !step8.total_reduced_gibbs || !step8.nlp
        || step8.active_phases.size() < 2 || !evaluator) {
        result.reason = "invalid_step9_input";
        return result;
    }
    Held2PhysicalCertificate physical{
        step8.nlp->primal_residual_inf,
        step8.ordinary_balance_inf,
        step8.electroneutrality_inf,
        step8.pressure_residual_inf,
        std::max({
            step8.nlp->stationarity_residual_inf,
            step8.nlp->dual_sign_violation_inf,
            step8.nlp->complementarity_inf,
        }),
        false,
    };
    physical.accepted = step8.nlp->accepted;
    result.physical = physical;

    std::vector<std::pair<std::uint64_t, Held2StateEvaluation>> states;
    for (const Held2Phase& phase : step8.active_phases) {
        states.emplace_back(
            phase.stable_id,
            evaluator(
                phase.independent_modified_fractions,
                std::log(phase.volume)
            )
        );
        ++result.timing.provider_evaluations;
    }
    std::sort(
        states.begin(),
        states.end(),
        [](const auto& left, const auto& right) {
            return left.first < right.first;
        }
    );
    bool potentials_passed = true;
    const std::size_t component_count =
        states.front().second.modified_potentials.size();
    if (component_count == 0
        || !std::all_of(states.begin(), states.end(), [&](const auto& state) {
            return state.second.modified_potentials.size() == component_count;
        })) {
        result.reason = "modified_potential_evidence_missing";
        return result;
    }
    for (std::size_t phase = 0; phase + 1 < states.size(); ++phase) {
        for (std::size_t component = 0; component < component_count; ++component) {
            const double left =
                states[phase].second.modified_potentials[component];
            const double numerator = left
                - states[phase + 1].second.modified_potentials[component];
            const double ratio = left == 0.0
                ? numerator == 0.0
                    ? 0.0
                    : std::copysign(
                        std::numeric_limits<double>::infinity(), numerator
                    )
                : std::abs(numerator / left);
            const bool passed = std::isfinite(ratio)
                && audit_held2_tolerance(
                    kHeld2PaperPotentialRatio, ratio
                ).passed;
            potentials_passed = potentials_passed && passed;
            result.potential_comparisons.push_back({
                component,
                states[phase].first,
                states[phase + 1].first,
                numerator,
                left,
                ratio,
                passed,
            });
        }
    }
    const double gap = *step4.upper_bound - *step8.total_reduced_gibbs;
    result.free_energy_gap = gap;
    const bool gap_passed = gap >= -kHeld2PaperFreeEnergyGap.atol
        && audit_held2_tolerance(kHeld2PaperFreeEnergyGap, gap).passed;
    Held2ProgressEvent progress;
    progress.kind = Held2ProgressKind::Certificate;
    progress.stage = "STEP 9";
    progress.status = gap_passed && potentials_passed ? "passed" : "failed";
    progress.objective = *step8.total_reduced_gibbs;
    progress.upper_bound = *step4.upper_bound;
    progress.gap = gap;
    observe_held2(observer, progress);
    if (!gap_passed || !potentials_passed) {
        result.outcome = Held2Step9Outcome::PaperConvergenceFailed;
        result.next_action = !gap_passed
            ? Held2Step9Action::ReturnStageII
            : Held2Step9Action::RunStep10;
        result.reason = !gap_passed
            ? "paper_free_energy_convergence_failed"
            : "paper_potential_convergence_failed";
        result.timing.terminal_status = "complete";
        result.timing.terminal_reason = result.reason;
        result.timing.next_step = 4;
        return result;
    }
    result.outcome = Held2Step9Outcome::Converged;
    result.next_action = Held2Step9Action::RunStep10;
    result.reason = "step9_complete";
    result.timing.terminal_status = "complete";
    result.timing.terminal_reason = result.reason;
    result.timing.next_step = 10;
    return result;
}

}  // namespace epcsaft_equilibrium
