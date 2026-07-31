#include "held2_step3.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace epcsaft_equilibrium {
namespace {

std::vector<double> bounding_point(
    const Held2Coordinates& coordinates,
    const std::vector<double>& feed,
    std::size_t coordinate,
    double direction
) {
    double distance = std::numeric_limits<double>::infinity();
    for (const Held2PolytopeConstraint& constraint :
         coordinates.polytope_constraints) {
        const double slope =
            direction * constraint.coefficients[coordinate];
        if (slope <= 0.0) {
            continue;
        }
        distance = std::min(
            distance,
            (
                constraint.upper_bound
                - std::inner_product(
                    constraint.coefficients.begin(),
                    constraint.coefficients.end(),
                    feed.begin(),
                    0.0
                )
            ) / slope
        );
    }
    if (!std::isfinite(distance) || distance <= 0.0) {
        throw std::invalid_argument("feed does not have two-sided bounds");
    }
    std::vector<double> point = feed;
    point[coordinate] += 0.5 * direction * distance;
    static_cast<void>(held2_lift_independent_fractions(
        coordinates, point
    ));
    return point;
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
    Held2PersistentState state;
    state.coordinates = coordinates;
    state.feed = *step1.independent_feed;
    state.feed_reduced_gibbs = step2.reference->objective;
    state.upper_bound = state.feed_reduced_gibbs;
    state.M.push_back({
        0, state.feed, step2.reference->volume,
        std::numeric_limits<double>::quiet_NaN(),
        step2.reference->objective, step2.reference->gradient,
        "homogeneous_feed",
    });
    for (std::size_t distinguished = 0;
         distinguished < dimension;
         ++distinguished) {
        for (int side = 0; side < 2; ++side) {
            std::vector<double> independent;
            try {
                independent = bounding_point(
                    coordinates,
                    state.feed,
                    distinguished,
                    side == 0 ? -1.0 : 1.0
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
                root.state.gradient,
                side == 0 ? "appendix_c_lower" : "appendix_c_upper",
            });
        }
    }
    if (dimension > 1) {
        if (!step2.negative_witness
            || step2.negative_witness->modified_fractions.size()
                != coordinates.retained_indices.size()
            || step2.negative_witness->reduced_gibbs_gradient.size()
                < dimension) {
            return fail("negative_witness_state_incomplete");
        }
        std::vector<double> witness_independent;
        witness_independent.reserve(dimension);
        for (std::size_t provider : coordinates.independent_indices) {
            const auto position = std::find(
                coordinates.retained_indices.begin(),
                coordinates.retained_indices.end(),
                provider
            );
            if (position == coordinates.retained_indices.end()) {
                return fail("negative_witness_coordinate_missing");
            }
            witness_independent.push_back(
                step2.negative_witness->modified_fractions[
                    static_cast<std::size_t>(
                        position - coordinates.retained_indices.begin()
                    )
                ]
            );
        }
        state.M.push_back({
            static_cast<std::uint64_t>(state.M.size()),
            std::move(witness_independent),
            step2.negative_witness->volume,
            std::numeric_limits<double>::quiet_NaN(),
            step2.negative_witness->reduced_gibbs,
            step2.negative_witness->reduced_gibbs_gradient,
            "stage_i_negative_witness",
        });
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
