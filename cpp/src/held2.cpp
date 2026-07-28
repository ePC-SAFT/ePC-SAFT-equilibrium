#include "held2.hpp"
#include "held2_tolerances.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace epcsaft_equilibrium {
namespace {

void require_finite_vector(const std::vector<double>& values, const char* name) {
    if (!std::all_of(values.begin(), values.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
}

}  // namespace

Held2PressureEnvelopeResult evaluate_held2_pressure_envelope(
    const std::vector<double>& independent_modified_fractions,
    const std::array<double, 2>& molar_volume_bounds,
    const Held2StateEvaluator& evaluator,
    int initial_interval_count,
    int maximum_subdivision_depth
) {
    if (independent_modified_fractions.empty()) {
        throw std::invalid_argument(
            "HELD2 pressure envelope requires independent composition coordinates"
        );
    }
    require_finite_vector(
        independent_modified_fractions,
        "HELD2 pressure-envelope composition"
    );
    if (!std::isfinite(molar_volume_bounds[0])
        || !std::isfinite(molar_volume_bounds[1])
        || molar_volume_bounds[0] <= 0.0
        || molar_volume_bounds[1] <= molar_volume_bounds[0]) {
        throw std::invalid_argument(
            "HELD2 pressure-envelope volume bounds must be finite, positive, and ordered"
        );
    }
    if (initial_interval_count < 1 || maximum_subdivision_depth < 0) {
        throw std::invalid_argument(
            "HELD2 pressure-envelope scan controls must be nonnegative"
        );
    }

    struct EvaluatedPoint {
        Held2PressureScanPoint diagnostic;
        Held2StateEvaluation state;
    };
    struct WorkingInterval {
        double lower = 0.0;
        double upper = 0.0;
        int depth = 0;
    };
    struct Bracket {
        double lower = 0.0;
        double upper = 0.0;
    };

    Held2PressureEnvelopeResult result;
    std::map<double, EvaluatedPoint> cache;
    const auto evaluate_point = [
        &cache,
        &evaluator,
        &independent_modified_fractions,
        &result
    ](double log_volume) -> const EvaluatedPoint& {
        const auto found = cache.find(log_volume);
        if (found != cache.end()) {
            return found->second;
        }
        EvaluatedPoint point;
        point.diagnostic.log_volume = log_volume;
        point.diagnostic.volume = std::exp(log_volume);
        try {
            point.state = evaluator(
                independent_modified_fractions,
                log_volume
            );
            if (!std::isfinite(point.state.volume)
                || point.state.volume <= 0.0
                || !std::isfinite(point.state.objective)
                || !std::isfinite(point.state.pressure_stationarity_relative)
                || !std::isfinite(
                    point.state.pressure_stationarity_derivative_log_volume
                )
                || point.state.hessian.empty()
                || !std::isfinite(point.state.hessian.back())) {
                throw std::invalid_argument(
                    "HELD2 pressure-envelope evaluation returned incomplete evidence"
                );
            }
            point.diagnostic.volume = point.state.volume;
            point.diagnostic.pressure_residual =
                point.state.pressure_stationarity_relative;
            point.diagnostic.pressure_derivative_log_volume =
                point.state.pressure_stationarity_derivative_log_volume;
            point.diagnostic.objective = point.state.objective;
            point.diagnostic.valid = true;
        } catch (const std::exception& error) {
            point.diagnostic.failure = error.what();
            ++result.evaluation_failure_count;
        }
        return cache.emplace(log_volume, std::move(point)).first->second;
    };
    const auto brackets_zero = [](double left, double right) {
        return (left <= 0.0 && right >= 0.0)
            || (left >= 0.0 && right <= 0.0);
    };
    const auto pressure_bracketed = [&brackets_zero](
        const EvaluatedPoint& left,
        const EvaluatedPoint& right
    ) {
        return left.diagnostic.valid && right.diagnostic.valid
            && brackets_zero(
                left.diagnostic.pressure_residual,
                right.diagnostic.pressure_residual
            );
    };
    const auto stationary_bracketed = [&brackets_zero](
        const EvaluatedPoint& left,
        const EvaluatedPoint& right
    ) {
        return left.diagnostic.valid && right.diagnostic.valid
            && brackets_zero(
                left.diagnostic.pressure_derivative_log_volume,
                right.diagnostic.pressure_derivative_log_volume
            );
    };

    const double lower_log_volume = std::log(molar_volume_bounds[0]);
    const double upper_log_volume = std::log(molar_volume_bounds[1]);
    result.lower_log_volume = lower_log_volume;
    result.upper_log_volume = upper_log_volume;
    const double span = upper_log_volume - lower_log_volume;
    std::vector<WorkingInterval> pending;
    pending.reserve(static_cast<std::size_t>(initial_interval_count));
    for (int interval = 0; interval < initial_interval_count; ++interval) {
        const double lower = interval == 0
            ? lower_log_volume
            : lower_log_volume + span * static_cast<double>(interval)
                / static_cast<double>(initial_interval_count);
        const double upper = interval + 1 == initial_interval_count
            ? upper_log_volume
            : lower_log_volume + span * static_cast<double>(interval + 1)
                / static_cast<double>(initial_interval_count);
        pending.push_back({lower, upper, 0});
    }

    std::vector<Bracket> pressure_brackets;
    std::vector<Bracket> stationary_brackets;
    std::size_t cursor = 0;
    while (cursor < pending.size()) {
        const WorkingInterval interval = pending[cursor++];
        const double midpoint = 0.5 * (interval.lower + interval.upper);
        const EvaluatedPoint& left = evaluate_point(interval.lower);
        const EvaluatedPoint& middle = evaluate_point(midpoint);
        const EvaluatedPoint& right = evaluate_point(interval.upper);
        if (!left.diagnostic.valid || !middle.diagnostic.valid
            || !right.diagnostic.valid) {
            result.intervals.push_back({
                interval.lower,
                interval.upper,
                interval.depth,
                "invalid_evaluation",
            });
            continue;
        }

        const bool left_pressure = pressure_bracketed(left, middle);
        const bool right_pressure = pressure_bracketed(middle, right);
        const bool left_stationary = stationary_bracketed(left, middle);
        const bool right_stationary = stationary_bracketed(middle, right);
        const auto joint_root = [](const EvaluatedPoint& point) {
            return audit_held2_tolerance(
                       kHeld2RootPressure,
                       point.diagnostic.pressure_residual
                   ).passed
                && audit_held2_tolerance(
                       kHeld2RootStationary,
                       point.diagnostic.pressure_derivative_log_volume
                   ).passed;
        };
        const bool left_joint_root = joint_root(left) || joint_root(middle);
        const bool right_joint_root = joint_root(middle) || joint_root(right);
        const bool ambiguous =
            (left_pressure && left_stationary && !left_joint_root)
            || (right_pressure && right_stationary && !right_joint_root);
        if (ambiguous && interval.depth < maximum_subdivision_depth) {
            pending.push_back({interval.lower, midpoint, interval.depth + 1});
            pending.push_back({midpoint, interval.upper, interval.depth + 1});
            continue;
        }
        if (ambiguous) {
            ++result.refinement_failure_count;
            result.intervals.push_back({
                interval.lower,
                interval.upper,
                interval.depth,
                "unresolved_multiple_topology",
            });
            continue;
        }
        if (left_pressure) {
            pressure_brackets.push_back({interval.lower, midpoint});
        }
        if (right_pressure) {
            pressure_brackets.push_back({midpoint, interval.upper});
        }
        if (left_stationary) {
            stationary_brackets.push_back({interval.lower, midpoint});
        }
        if (right_stationary) {
            stationary_brackets.push_back({midpoint, interval.upper});
        }
        std::string status = "root_free_finite_scan";
        if ((left_pressure || right_pressure)
            && (left_stationary || right_stationary)) {
            status = "pressure_and_stationary_brackets";
        } else if (left_pressure || right_pressure) {
            status = "pressure_bracket";
        } else if (left_stationary || right_stationary) {
            status = "stationary_bracket";
        }
        result.intervals.push_back({
            interval.lower,
            interval.upper,
            interval.depth,
            std::move(status),
        });
    }

    const auto refine = [
        &evaluate_point,
        &brackets_zero,
        &result
    ](Bracket bracket, bool stationary) -> std::optional<EvaluatedPoint> {
        const auto residual = [stationary](const EvaluatedPoint& point) {
            return stationary
                ? point.diagnostic.pressure_derivative_log_volume
                : point.diagnostic.pressure_residual;
        };
        EvaluatedPoint left = evaluate_point(bracket.lower);
        EvaluatedPoint right = evaluate_point(bracket.upper);
        if (!left.diagnostic.valid || !right.diagnostic.valid
            || !brackets_zero(residual(left), residual(right))) {
            ++result.refinement_failure_count;
            return std::nullopt;
        }
        const double tolerance = stationary
            ? kHeld2RootStationary.atol
            : kHeld2RootPressure.atol;
        EvaluatedPoint best = std::abs(residual(left)) <= std::abs(residual(right))
            ? left
            : right;
        for (int iteration = 0; iteration < 120; ++iteration) {
            if (audit_held2_tolerance(
                    kHeld2RootLogVolumeWidth,
                    bracket.upper - bracket.lower
                ).passed
                && (stationary || std::abs(residual(best)) <= tolerance)) {
                return best;
            }
            const double midpoint = 0.5 * (bracket.lower + bracket.upper);
            if (!(midpoint > bracket.lower) || !(midpoint < bracket.upper)) {
                if (std::abs(residual(best)) <= tolerance) {
                    return best;
                }
                ++result.refinement_failure_count;
                return std::nullopt;
            }
            const EvaluatedPoint middle = evaluate_point(midpoint);
            if (!middle.diagnostic.valid) {
                ++result.refinement_failure_count;
                return std::nullopt;
            }
            if (std::abs(residual(middle)) < std::abs(residual(best))) {
                best = middle;
            }
            if (brackets_zero(residual(left), residual(middle))) {
                bracket.upper = midpoint;
                right = middle;
            } else {
                bracket.lower = midpoint;
                left = middle;
            }
        }
        ++result.refinement_failure_count;
        return std::nullopt;
    };

    std::vector<Held2PressureRoot> candidates;
    const auto make_root = [lower_log_volume, upper_log_volume](
        const EvaluatedPoint& point,
        std::string origin
    ) {
        Held2PressureRoot root;
        root.log_volume = point.diagnostic.log_volume;
        root.volume = point.state.volume;
        root.objective = point.state.objective;
        root.pressure_residual = point.state.pressure_stationarity_relative;
        root.pressure_derivative_log_volume =
            point.state.pressure_stationarity_derivative_log_volume;
        root.objective_curvature_log_volume = point.state.hessian.back();
        root.origin = std::move(origin);
        root.boundary = audit_held2_tolerance(
                            kHeld2RootBoundary,
                            root.log_volume - lower_log_volume
                        ).passed
            || audit_held2_tolerance(
                   kHeld2RootBoundary,
                   upper_log_volume - root.log_volume
               ).passed;
        root.state = point.state;
        return root;
    };
    for (const Bracket& bracket : pressure_brackets) {
        const std::optional<EvaluatedPoint> root = refine(bracket, false);
        if (root.has_value()) {
            candidates.push_back(make_root(*root, "sign_change"));
        }
    }

    std::vector<double> stationary_points;
    for (const Bracket& bracket : stationary_brackets) {
        const std::optional<EvaluatedPoint> stationary = refine(bracket, true);
        if (!stationary.has_value()) {
            continue;
        }
        const bool duplicate = std::any_of(
            stationary_points.begin(),
            stationary_points.end(),
            [&stationary](double retained) {
                return std::abs(
                    retained - stationary->diagnostic.log_volume
                ) <= kHeld2RootDuplicate.atol;
            }
        );
        if (duplicate) {
            continue;
        }
        stationary_points.push_back(stationary->diagnostic.log_volume);
        ++result.stationary_point_count;
        if (audit_held2_tolerance(
                kHeld2RootPressure,
                stationary->diagnostic.pressure_residual
            ).passed) {
            ++result.tangential_root_count;
            candidates.push_back(make_root(*stationary, "tangential"));
        }
    }

    std::sort(
        candidates.begin(),
        candidates.end(),
        [](const Held2PressureRoot& left, const Held2PressureRoot& right) {
            return left.log_volume < right.log_volume;
        }
    );
    for (Held2PressureRoot& candidate : candidates) {
        if (!result.roots.empty()
            && std::abs(
                result.roots.back().log_volume - candidate.log_volume
            ) <= kHeld2RootDuplicate.atol) {
            result.roots.back().boundary = result.roots.back().boundary
                || candidate.boundary;
            if (candidate.origin == "tangential") {
                result.roots.back().origin = candidate.origin;
            }
            ++result.deduplicated_root_count;
            continue;
        }
        result.roots.push_back(std::move(candidate));
    }

    std::vector<int> stable_indices;
    bool unclassified = false;
    for (std::size_t index = 0; index < result.roots.size(); ++index) {
        Held2PressureRoot& root = result.roots[index];
        if (root.boundary) {
            ++result.boundary_root_count;
        }
        const double mechanical_margin = std::min(
            std::abs(root.pressure_derivative_log_volume),
            std::abs(root.objective_curvature_log_volume)
        );
        const bool mechanical_margin_passed = audit_held2_tolerance(
            kHeld2MechanicalMargin,
            mechanical_margin
        ).passed;
        const bool marginal = !mechanical_margin_passed;
        if (marginal) {
            root.mechanical_class = "marginal";
            ++result.marginal_root_count;
        } else if (root.pressure_derivative_log_volume < 0.0
            && root.objective_curvature_log_volume > 0.0) {
            root.mechanical_class = "strict_stable";
            stable_indices.push_back(static_cast<int>(index));
        } else if (root.pressure_derivative_log_volume > 0.0
            && root.objective_curvature_log_volume < 0.0) {
            root.mechanical_class = "unstable";
        } else {
            unclassified = true;
        }
    }

    if (result.evaluation_failure_count != 0) {
        result.failure_reason = "evaluation_failed";
    } else if (result.refinement_failure_count != 0) {
        result.failure_reason = "refinement_failed";
    } else if (result.boundary_root_count != 0) {
        result.failure_reason = "boundary_root";
    } else if (result.marginal_root_count != 0) {
        result.failure_reason = "marginal_root";
    } else if (unclassified) {
        result.failure_reason = "mechanical_classification_failed";
    } else if (result.roots.empty()) {
        result.failure_reason = "root_not_found";
    } else if (stable_indices.empty()) {
        result.failure_reason = "strict_stable_root_not_found";
    } else {
        double lowest_objective = std::numeric_limits<double>::infinity();
        for (int index : stable_indices) {
            lowest_objective = std::min(
                lowest_objective,
                result.roots[static_cast<std::size_t>(index)].objective
            );
        }
        std::vector<int> lowest;
        for (int index : stable_indices) {
            const double objective =
                result.roots[static_cast<std::size_t>(index)].objective;
            const double scale = std::max(std::abs(objective), std::abs(lowest_objective));
            if (audit_held2_tolerance(
                    kHeld2StableObjectiveTie,
                    objective - lowest_objective,
                    scale
                ).passed) {
                lowest.push_back(index);
            }
        }
        if (lowest.size() != 1) {
            result.objective_tie_count = static_cast<int>(lowest.size()) - 1;
            result.failure_reason = "stable_objective_tie";
        } else {
            result.selected_root_index = lowest.front();
            result.outcome = "selected";
        }
    }

    std::sort(
        result.intervals.begin(),
        result.intervals.end(),
        [](const Held2PressureScanInterval& left,
           const Held2PressureScanInterval& right) {
            return left.lower_log_volume < right.lower_log_volume;
        }
    );
    result.scan_points.reserve(cache.size());
    for (const auto& [log_volume, point] : cache) {
        static_cast<void>(log_volume);
        result.scan_points.push_back(point.diagnostic);
    }
    return result;
}

}  // namespace epcsaft_equilibrium
