#include "held2_step3.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace epcsaft_equilibrium {
namespace {

std::vector<double> independent_from_physical(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical
) {
    const std::vector<double> modified =
        held2_transform_physical_fractions(coordinates, physical);
    std::vector<double> independent;
    for (std::size_t provider : coordinates.independent_indices) {
        const auto retained = std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            provider
        );
        independent.push_back(modified[static_cast<std::size_t>(
            retained - coordinates.retained_indices.begin()
        )]);
    }
    static_cast<void>(held2_lift_independent_fractions(
        coordinates, independent
    ));
    return independent;
}

}  // namespace

Held2Step3Result run_held2_step3(
    const Held2Step1Result& step1,
    const Held2Step2Result& step2,
    const Held2PressureRootEvaluator& pressure_roots,
    Held2ProgressObserver*
) {
    Held2Step3Result result;
    result.timing.invocation_count = 1;
    const auto fail = [&](const char* reason) {
        result.reason = reason;
        result.timing.terminal_status = "indeterminate";
        result.timing.terminal_reason = reason;
        return result;
    };
    if (step1.status != "complete" || !step1.coordinates
        || !step1.independent_feed || !pressure_roots
        || step2.outcome != Held2Step2Outcome::NegativeWitness
        || !step2.reference) {
        return fail("invalid_step3_input");
    }
    const Held2Coordinates& coordinates = *step1.coordinates;
    const std::size_t dimension = coordinates.independent_indices.size();
    const std::vector<double> feed_physical =
        held2_lift_independent_fractions(
            coordinates, *step1.independent_feed
        );
    Held2PersistentState state;
    state.feed = *step1.independent_feed;
    state.feed_reduced_gibbs = step2.reference->objective;
    state.upper_bound = state.feed_reduced_gibbs;
    state.M.push_back({
        0, state.feed, step2.reference->volume,
        std::numeric_limits<double>::quiet_NaN(),
        step2.reference->objective, "homogeneous_feed",
    });
    const double eliminated_charge =
        coordinates.charges[coordinates.eliminated_index];
    for (std::size_t distinguished = 0;
         distinguished < dimension;
         ++distinguished) {
        const std::size_t provider =
            coordinates.independent_indices[distinguished];
        for (int side = 0; side < 2; ++side) {
            std::vector<double> physical(coordinates.charges.size(), 0.0);
            physical[provider] = side == 0
                ? 0.5 * feed_physical[provider]
                : 0.5 * (
                    1.0 / (
                        1.0 - coordinates.charges[provider]
                            / eliminated_charge
                    ) + feed_physical[provider]
                );
            for (std::size_t compact = 0; compact < dimension; ++compact) {
                if (compact == distinguished) {
                    continue;
                }
                const std::size_t other =
                    coordinates.independent_indices[compact];
                physical[other] = (
                    1.0 / (
                        1.0 - coordinates.charges[other]
                            / eliminated_charge
                    ) - physical[provider]
                ) / static_cast<double>(coordinates.charges.size() - 1);
            }
            for (std::size_t other : coordinates.independent_indices) {
                physical[coordinates.eliminated_index] -=
                    coordinates.charges[other] / eliminated_charge
                    * physical[other];
            }
            physical[coordinates.dependent_index] = 1.0
                - std::accumulate(
                    physical.begin(), physical.end(), 0.0
                );
            std::vector<double> independent;
            try {
                independent = independent_from_physical(
                    coordinates, physical
                );
            } catch (...) {
                return fail("appendix_c_state_outside_domain");
            }
            Held2PressureEnvelopeResult envelope;
            try {
                envelope = pressure_roots(independent);
            } catch (...) {
                return fail("appendix_c_pressure_root_failed");
            }
            ++result.timing.provider_evaluations;
            if (envelope.outcome != "selected"
                || envelope.selected_root_index < 0) {
                return fail(
                    envelope.failure_reason.empty()
                    ? "appendix_c_pressure_root_indeterminate"
                    : envelope.failure_reason.c_str()
                );
            }
            const Held2PressureRoot& root = envelope.roots[
                static_cast<std::size_t>(envelope.selected_root_index)
            ];
            state.M.push_back({
                static_cast<std::uint64_t>(state.M.size()),
                std::move(independent),
                root.volume,
                std::numeric_limits<double>::quiet_NaN(),
                root.objective,
                side == 0 ? "appendix_c_lower" : "appendix_c_upper",
            });
        }
    }
    result.status = "complete";
    result.reason = "step3_complete";
    result.timing.terminal_status = result.status;
    result.timing.terminal_reason = result.reason;
    result.timing.next_step = 4;
    result.state = std::move(state);
    return result;
}

}  // namespace epcsaft_equilibrium
