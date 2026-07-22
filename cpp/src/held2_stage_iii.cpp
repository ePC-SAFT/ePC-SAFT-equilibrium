#include "held2.hpp"
#include "held2_tolerances.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace epcsaft_equilibrium {
namespace {

constexpr double kVolumeLower = 0.5;
constexpr double kVolumeUpper = 1.5;

Held2StateEvaluator manufactured_evaluator(const Held2Coordinates& coordinates) {
    return [coordinates](
        const std::vector<double>& independent,
        double log_volume
    ) {
        return evaluate_held2_manufactured_state(
            coordinates, independent, log_volume
        );
    };
}

std::vector<Held2StageIICandidate> manufactured_candidates(
    const Held2Coordinates& coordinates,
    const std::vector<std::array<double, 2>>& points
) {
    std::vector<Held2StageIICandidate> result;
    result.reserve(points.size());
    for (const auto& point : points) {
        const Held2StateEvaluation state = evaluate_held2_manufactured_state(
            coordinates, {point[0]}, std::log(point[1])
        );
        result.push_back({
            state.modified_fractions,
            {point[0]},
            point[1],
            std::log(point[1]),
            0.0,
        });
    }
    return result;
}

}  // namespace

Held2StageIIIRetirementDecision held2_stage_iii_retirement_decision(
    double phase_fraction,
    double lower_bound_multiplier,
    double upper_bound_multiplier,
    double reduced_derivative,
    bool remaining_balance_feasible
) {
    Held2StageIIIRetirementDecision result;
    result.complementarity_inf_norm = std::max(
        std::abs(phase_fraction * lower_bound_multiplier),
        std::abs((1.0 - phase_fraction) * upper_bound_multiplier)
    );
    result.stationarity_residual = std::abs(
        reduced_derivative - lower_bound_multiplier + upper_bound_multiplier
    );
    if (audit_held2_tolerance(kHeld2PhaseActivity, phase_fraction).passed) {
        result.reason = "phase_amount_active";
        return result;
    }
    if (!remaining_balance_feasible) {
        result.reason = "remaining_balance_infeasible";
        return result;
    }
    const double retirement_margin = std::min(
        lower_bound_multiplier, reduced_derivative
    );
    if (!audit_held2_tolerance(
            kHeld2PhaseRetirementMargin, retirement_margin
        ).passed) {
        result.reason = "descent_or_marginal_phase";
        return result;
    }
    const double sign_violation = std::max({
        0.0, -lower_bound_multiplier, -upper_bound_multiplier
    });
    if (!audit_held2_tolerance(kHeld2Stage3DualSign, sign_violation).passed) {
        result.reason = "invalid_multiplier_sign";
        return result;
    }
    if (!audit_held2_tolerance(
            kHeld2Stage3Complementarity,
            result.complementarity_inf_norm
        ).passed) {
        result.reason = "complementarity_failed";
        return result;
    }
    if (!audit_held2_tolerance(
            kHeld2Stage3Stationarity,
            result.stationarity_residual
        ).passed) {
        result.reason = "reduced_derivative_inconsistent";
        return result;
    }
    result.retire = true;
    result.reason = "kkt_inactive";
    return result;
}

Held2StageIIINlpEvaluation evaluate_held2_manufactured_stage_iii_nlp(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::vector<std::array<double, 2>>& candidates,
    const std::vector<double>& variables,
    const std::vector<double>& equality_multipliers
) {
    const Held2Coordinates coordinates = make_held2_coordinates(charges);
    return evaluate_held2_stage_iii_nlp(
        coordinates,
        physical_feed,
        manufactured_evaluator(coordinates),
        candidates.size(),
        variables,
        equality_multipliers
    );
}

Held2StageIIIResult solve_held2_manufactured_stage_iii(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::vector<std::array<double, 2>>& candidates
) {
    const Held2Coordinates coordinates = make_held2_coordinates(charges);
    return solve_held2_stage_iii(
        coordinates,
        physical_feed,
        manufactured_candidates(coordinates, candidates),
        manufactured_evaluator(coordinates),
        std::vector<std::array<double, 2>>(
            candidates.size(),
            {std::log(kVolumeLower), std::log(kVolumeUpper)}
        )
    );
}

}  // namespace epcsaft_equilibrium
