#include "held2_algorithm.hpp"

#include <chrono>
#include <cmath>
#include <ctime>

namespace epcsaft_equilibrium {
namespace {

template <typename Callable>
auto run_step(
    int step,
    Callable&& callable,
    Held2ProgressObserver* observer
) {
    Held2ProgressEvent progress;
    progress.kind = Held2ProgressKind::StageStart;
    progress.stage = "STEP " + std::to_string(step);
    observe_held2(observer, progress);
    const auto wall_start = std::chrono::steady_clock::now();
    const std::clock_t cpu_start = std::clock();
    auto result = callable();
    result.timing.step = step;
    result.timing.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start
    ).count();
    result.timing.cpu_seconds = static_cast<double>(
        std::clock() - cpu_start
    ) / static_cast<double>(CLOCKS_PER_SEC);
    progress = {};
    progress.kind = Held2ProgressKind::Certificate;
    progress.stage = "STEP " + std::to_string(step);
    progress.status = result.timing.terminal_status;
    progress.reason = result.timing.terminal_reason;
    progress.wall_seconds = result.timing.wall_seconds;
    progress.cpu_seconds = result.timing.cpu_seconds;
    progress.provider_evaluations = result.timing.provider_evaluations;
    progress.optimizer_solves = result.timing.optimizer_solves;
    progress.optimizer_iterations = result.timing.optimizer_iterations;
    observe_held2(observer, progress);
    return result;
}

}  // namespace

Held2AlgorithmResult run_held2_algorithm(
    const Held2ThermodynamicAccess& thermodynamics,
    const Held2Input& input,
    const Held2ResourceProfile& resources,
    Held2ProgressObserver* observer
) {
    Held2AlgorithmResult result;
    Held2ProgressEvent progress;
    progress.kind = Held2ProgressKind::CaseStart;
    progress.case_id = "installed-held2-paper-rewrite";
    progress.temperature_k = input.temperature_k;
    progress.pressure_pa = input.pressure_pa;
    observe_held2(observer, progress);
    const auto fail = [&](const char* stage, const std::string& reason) {
        result.failure_stage = stage;
        result.failure_reason = reason;
        return result;
    };
    result.step1 = run_step(1, [&] {
        return run_held2_step1(
        thermodynamics.component_ids,
        thermodynamics.charges,
        input.temperature_k,
        input.pressure_pa,
        input.overall_mole_fractions_provider_order,
        thermodynamics.volume_bounds_physical,
        thermodynamics.total_ion_mole_fraction_max
        );
    }, observer);
    result.step_timings.push_back(result.step1.timing);
    if (result.step1.status != "complete") {
        return fail("step1", result.step1.reason);
    }
    result.step2 = run_step(2, [&] {
        return run_held2_step2(
            result.step1,
            thermodynamics.evaluate,
            resources.step2_search_budget,
            observer
        );
    }, observer);
    result.step_timings.push_back(result.step2->timing);
    if (result.step2->outcome == Held2Step2Outcome::Indeterminate) {
        return fail("step2", result.step2->reason);
    }
    if (result.step2->outcome
        == Held2Step2Outcome::NoNegativeWitnessDetected) {
        const Held2StateEvaluation& reference = *result.step2->reference;
        result.phases = {{
            0,
            1.0,
            *result.step1.independent_feed,
            input.overall_mole_fractions_provider_order,
            reference.volume,
            thermodynamics.packing_fraction(
                *result.step1.independent_feed, reference.volume
            ),
        }};
        result.outcome = "one_phase_no_negative_witness_detected";
        return result;
    }
    const Held2PressureRootEvaluator pressure_roots = [&](const auto& composition) {
        return evaluate_held2_pressure_envelope(
            composition,
            (*result.step1.volume_bounds)(composition),
            thermodynamics.evaluate,
            64,
            8
        );
    };
    result.step3 = run_step(3, [&] {
        return run_held2_step3(
            result.step1, *result.step2, pressure_roots, observer
        );
    }, observer);
    result.step_timings.push_back(result.step3->timing);
    if (result.step3->status != "complete" || !result.step3->state) {
        return fail("step3", result.step3->reason);
    }
    Held2PersistentState state = *result.step3->state;
    while (true) {
        result.step4_history.push_back(run_step(4, [&] {
            return run_held2_step4(state, observer);
        }, observer));
        result.step_timings.push_back(result.step4_history.back().timing);
        const Held2Step4Result& step4 = result.step4_history.back();
        if (step4.status != "complete") {
            result.final_state = state;
            return fail("step4", step4.reason);
        }
        result.step5_history.push_back(run_step(5, [&] {
            return run_held2_step5(
                result.step1,
                step4,
                state,
                thermodynamics.evaluate,
                thermodynamics.packing_fraction,
                resources,
                observer
            );
        }, observer));
        result.step_timings.push_back(result.step5_history.back().timing);
        const Held2Step5Result& step5 = result.step5_history.back();
        if (step5.status != "complete") {
            result.final_state = state;
            return fail("step5", step5.reason);
        }
        result.step6_history.push_back(run_step(6, [&] {
            return run_held2_step6(
                result.step1,
                step4,
                state,
                thermodynamics.evaluate,
                thermodynamics.packing_fraction,
                observer
            );
        }, observer));
        result.step_timings.push_back(result.step6_history.back().timing);
        const Held2Step6Result& step6 = result.step6_history.back();
        if (step6.status != "complete") {
            result.final_state = state;
            return fail("step6", step6.reason);
        }
        const auto continue_stage_ii = [&](bool stage_iii_feedback) {
            result.step7_history.push_back(run_step(7, [&] {
                return run_held2_step7(
                    state,
                    step5,
                    step6,
                    resources,
                    stage_iii_feedback,
                    observer
                );
            }, observer));
            result.step_timings.push_back(
                result.step7_history.back().timing
            );
            return result.step7_history.back().status == "complete";
        };
        if (step6.candidates.size() < 2) {
            if (!continue_stage_ii(false)) {
                result.final_state = state;
                return fail("step7", result.step7_history.back().reason);
            }
            continue;
        }
        result.step8_history.push_back(run_step(8, [&] {
            return run_held2_step8(
                result.step1,
                step6,
                thermodynamics.evaluate,
                thermodynamics.packing_fraction,
                observer
            );
        }, observer));
        result.step_timings.push_back(result.step8_history.back().timing);
        const Held2Step8Result& step8 = result.step8_history.back();
        if (step8.outcome == Held2Step8Outcome::Indeterminate) {
            result.final_state = state;
            return fail("step8", step8.reason);
        }
        if (step8.outcome == Held2Step8Outcome::CertifiedInfeasible
            || step8.outcome == Held2Step8Outcome::InsufficientCandidates) {
            if (!continue_stage_ii(true)) {
                result.final_state = state;
                return fail(
                    "step7", result.step7_history.back().reason
                );
            }
            continue;
        }
        result.step9_history.push_back(run_step(9, [&] {
            return run_held2_step9(
                step4, step8, thermodynamics.evaluate, observer
            );
        }, observer));
        result.step_timings.push_back(result.step9_history.back().timing);
        const Held2Step9Result& step9 = result.step9_history.back();
        if (step9.outcome == Held2Step9Outcome::Indeterminate) {
            result.final_state = state;
            return fail("step9", step9.reason);
        }
        if (step9.outcome == Held2Step9Outcome::PaperConvergenceFailed) {
            if (!continue_stage_ii(true)) {
                result.final_state = state;
                return fail(
                    "step7", result.step7_history.back().reason
                );
            }
            continue;
        }
        result.step10 = run_step(10, [&] {
            return run_held2_step10(
                result.step1,
                step8,
                step9,
                thermodynamics.evaluate_trace
                    ? thermodynamics.evaluate_trace
                    : thermodynamics.evaluate,
                observer
            );
        }, observer);
        result.step_timings.push_back(result.step10->timing);
        result.final_state = state;
        result.upper_solve_count = state.upper_solve_count;
        if (result.step10->status != "complete") {
            return fail("step10", result.step10->reason);
        }
        result.phases = result.step10->phases;
        result.outcome = "physical_equilibrium_accepted";
        return result;
    }
}

}  // namespace epcsaft_equilibrium
