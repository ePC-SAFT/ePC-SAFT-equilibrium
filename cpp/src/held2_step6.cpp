#include "held2_step6.hpp"
#include "held2_tolerances.hpp"

#include <algorithm>
#include <cmath>

namespace epcsaft_equilibrium {

Held2Step6Result run_held2_step6(
    const Held2Step1Result& step1,
    const Held2Step4Result& step4,
    const Held2PersistentState& state,
    const Held2StateEvaluator& evaluator,
    const Held2PackingFractionEvaluator& packing_fraction,
    Held2ProgressObserver*
) {
    constexpr double kBoundGap = 1.0e-2;
    constexpr double kDerivativeRelative = 0.5;
    constexpr double kDerivativeFloor = 1.0e-8;
    constexpr double kPackingDistinct = 1.0e-3;
    constexpr double kCompositionDistinct = 1.0e-3;
    Held2Step6Result result;
    result.timing.invocation_count = 1;
    if (step4.status != "complete" || !step1.coordinates
        || !evaluator || !packing_fraction) {
        result.reason = "invalid_step6_input";
        return result;
    }
    for (const Held2MPoint& point : state.M) {
        Held2CandidateDecision decision;
        decision.insertion_id = point.insertion_id;
        const Held2StateEvaluation evaluation = evaluator(
            point.independent_modified_fractions, std::log(point.volume)
        );
        ++result.timing.provider_evaluations;
        double lower = evaluation.objective;
        for (std::size_t index = 0; index < state.feed.size(); ++index) {
            lower += state.multipliers[index] * (
                state.feed[index]
                - point.independent_modified_fractions[index]
            );
        }
        decision.gap_passed =
            std::abs(state.upper_bound - lower) <= kBoundGap;
        decision.derivative_passed = true;
        for (std::size_t index = 0; index < state.feed.size(); ++index) {
            if (point.independent_modified_fractions[index]
                <= step1.coordinates->independent_lower_bounds[index]
                    + kHeld2BoundActivity.atol) {
                continue;
            }
            decision.derivative_passed =
                decision.derivative_passed
                && std::abs(
                    evaluation.gradient[index] - state.multipliers[index]
                ) <= std::max(
                    kDerivativeFloor,
                    kDerivativeRelative * std::abs(state.multipliers[index])
                );
        }
        if (!decision.gap_passed || !decision.derivative_passed) {
            decision.reason = decision.gap_passed
                ? "derivative_failed"
                : "bound_gap_failed";
            result.decisions.push_back(std::move(decision));
            continue;
        }
        Held2MPoint candidate = point;
        candidate.reduced_gibbs = evaluation.objective;
        candidate.packing_fraction = packing_fraction(
            point.independent_modified_fractions, point.volume
        );
        decision.pairwise_distinct = std::all_of(
            result.candidates.begin(),
            result.candidates.end(),
            [&](const Held2MPoint& retained) {
                if (std::abs(
                        candidate.packing_fraction
                        - retained.packing_fraction
                    ) >= kPackingDistinct) {
                    return true;
                }
                for (std::size_t index = 0; index < state.feed.size(); ++index) {
                    if (std::abs(
                            candidate.independent_modified_fractions[index]
                            - retained.independent_modified_fractions[index]
                        ) >= kCompositionDistinct) {
                        return true;
                    }
                }
                return false;
            }
        );
        decision.retained = decision.pairwise_distinct;
        decision.reason = decision.retained ? "retained" : "not_distinct";
        if (decision.retained) {
            result.candidates.push_back(std::move(candidate));
        }
        result.decisions.push_back(std::move(decision));
    }
    result.status = "complete";
    result.reason = result.candidates.size() >= 2
        ? "candidate_set"
        : "insufficient_candidates";
    result.timing.terminal_status = result.status;
    result.timing.terminal_reason = result.reason;
    result.timing.next_step = result.candidates.size() >= 2 ? 8 : 7;
    return result;
}

}  // namespace epcsaft_equilibrium
