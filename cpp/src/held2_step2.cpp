#include "held2_step2.hpp"
#include "held2_tolerances.hpp"

#include <nlopt.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <functional>
#include <limits>
#include <stdexcept>

namespace epcsaft_equilibrium {
namespace {

struct Held2StageIJointEvaluation {
    std::vector<double> independent_modified_fractions;
    Held2StateEvaluation state;
    double normalized_volume = 0.0;
    double physical_total_ion_mole_fraction =
        std::numeric_limits<double>::quiet_NaN();
    double total_ion_mole_fraction_max =
        std::numeric_limits<double>::quiet_NaN();
    double tpd = std::numeric_limits<double>::infinity();
    std::string failure_reason;
};

using Held2StageIJointEvaluator = std::function<Held2StageIJointEvaluation(
    const std::vector<double>&
)>;

struct Held2StageIDirectResult {
    std::string failure_reason;
    int evaluation_count = 0;
    double minimum_tpd = std::numeric_limits<double>::infinity();
    std::optional<Held2StageIJointEvaluation> negative_candidate;
};

struct DirectContext {
    const Held2StageIJointEvaluator* evaluator = nullptr;
    Held2StageIDirectResult* result = nullptr;
    nlopt::opt* optimizer = nullptr;
    Held2ProgressObserver* observer = nullptr;
    bool stop_requested = false;
};

double pressure_residual(const Held2StageIJointEvaluation& evaluation) {
    return evaluation.state.pressure_stationarity_relative;
}

bool is_valid(const Held2StageIJointEvaluation& evaluation) {
    return std::isfinite(evaluation.tpd);
}

bool is_strict_negative(const Held2StageIJointEvaluation& evaluation) {
    return is_valid(evaluation)
        && audit_held2_tolerance(
            kHeld2TpdNegativeMargin, evaluation.tpd
        ).passed;
}

bool is_strict_stable_root(
    const Held2StageIJointEvaluation& evaluation
) {
    return is_valid(evaluation)
        && audit_held2_tolerance(
            kHeld2RootPressure,
            pressure_residual(evaluation)
        ).passed
        && audit_held2_tolerance(
            kHeld2MechanicalMargin,
            std::min(
                std::abs(
                    evaluation.state
                        .pressure_stationarity_derivative_log_volume
                ),
                std::abs(evaluation.state.hessian.back())
            )
        ).passed
        && evaluation.state.pressure_stationarity_derivative_log_volume < 0.0
        && evaluation.state.hessian.back() > 0.0;
}

double direct_objective(
    const std::vector<double>& chart_coordinates,
    std::vector<double>& gradient,
    void* opaque
) {
    auto& context = *static_cast<DirectContext*>(opaque);
    if (context.stop_requested) {
        // DIRECT may finish an already assembled sample batch after force_stop().
        // This return is ignored and is not a physical envelope evaluation.
        return 0.0;
    }
    if (!gradient.empty()) {
        throw std::invalid_argument("HELD2 DIRECT-L requested an unexpected gradient");
    }
    Held2StageIJointEvaluation evaluation =
        (*context.evaluator)(chart_coordinates);
    const int evaluation_index = context.result->evaluation_count++;
    Held2ProgressEvent progress;
    progress.kind = Held2ProgressKind::StageIEvaluation;
    progress.iteration = evaluation_index + 1;
    progress.objective = evaluation.tpd;
    progress.physical_total_ion_mole_fraction =
        evaluation.physical_total_ion_mole_fraction;
    progress.total_ion_mole_fraction_max =
        evaluation.total_ion_mole_fraction_max;
    progress.status = is_valid(evaluation) ? "admissible" : "failed";
    progress.reason = is_valid(evaluation) ? "" : "trial_evaluation_failed";
    if (evaluation.state.volume > 0.0) {
        progress.volume = evaluation.state.volume;
        progress.pressure_residual = pressure_residual(evaluation);
    }
    observe_held2(context.observer, progress);
    if (!is_valid(evaluation)) {
        Held2ProgressEvent failure;
        failure.kind = Held2ProgressKind::Failure;
        failure.stage = "STAGE I";
        failure.reason = "required_joint_evaluation_failed";
        observe_held2(context.observer, failure);
        context.result->failure_reason =
            "required_joint_evaluation_failed";
        context.stop_requested = true;
        context.optimizer->force_stop();
        return 0.0;
    }
    context.result->minimum_tpd = std::min(
        context.result->minimum_tpd,
        evaluation.tpd
    );
    if (is_strict_negative(evaluation)) {
        context.result->negative_candidate = evaluation;
        context.stop_requested = true;
        context.optimizer->force_stop();
        return evaluation.tpd;
    }
    return evaluation.tpd;
}

Held2StageIDirectResult run_direct_search(
    std::size_t joint_dimension,
    int provider_evaluation_budget,
    int provider_evaluations_per_trial,
    const Held2StageIJointEvaluator& evaluator,
    Held2ProgressObserver* observer
) {
    Held2StageIDirectResult result;
    nlopt::opt optimizer(nlopt::GN_DIRECT_L, joint_dimension);
    DirectContext context{&evaluator, &result, &optimizer, observer, false};
    optimizer.set_lower_bounds(std::vector<double>(joint_dimension, 0.0));
    optimizer.set_upper_bounds(std::vector<double>(joint_dimension, 1.0));
    optimizer.set_maxeval(
        provider_evaluation_budget / provider_evaluations_per_trial
    );
    optimizer.set_min_objective(direct_objective, &context);
    std::vector<double> initial(joint_dimension, 0.5);
    double minimum = std::numeric_limits<double>::infinity();
    try {
        const nlopt::result status = optimizer.optimize(initial, minimum);
        if (status != nlopt::MAXEVAL_REACHED) {
            result.failure_reason = "unexpected_solver_termination";
        }
    } catch (const nlopt::forced_stop&) {
    } catch (const std::exception& error) {
        result.failure_reason = std::string("solver_failure: ") + error.what();
    }
    return result;
}

}  // namespace

Held2TpdEvaluation evaluate_held2_tpd(
    const Held2StateEvaluation& reference,
    const std::vector<double>& feed,
    const Held2StateEvaluation& trial,
    const std::vector<double>& independent
) {
    const std::size_t dimension = feed.size();
    if (dimension == 0 || independent.size() != dimension
        || reference.gradient.size() != dimension + 1
        || trial.gradient.size() != dimension + 1
        || trial.hessian.size() != (dimension + 1) * (dimension + 1)
        || !std::isfinite(reference.objective)
        || !std::isfinite(trial.objective)) {
        throw std::invalid_argument("HELD2 Step-2 TPD evidence is incomplete");
    }
    Held2TpdEvaluation result;
    result.value = trial.objective - reference.objective;
    result.gradient = trial.gradient;
    for (std::size_t index = 0; index < dimension; ++index) {
        result.value -= reference.gradient[index]
            * (independent[index] - feed[index]);
        result.gradient[index] -= reference.gradient[index];
    }
    result.hessian = trial.hessian;
    return result;
}

Held2Step2Result run_held2_step2(
    const Held2Step1Result& step1,
    const Held2StateEvaluator& evaluator,
    int provider_evaluation_budget,
    Held2ProgressObserver* observer
) {
    constexpr int kPressureIntervals = 64;
    constexpr int kPressureSubdivisionDepth = 8;
    const auto wall_start = std::chrono::steady_clock::now();
    const std::clock_t cpu_start = std::clock();
    Held2Step2Result result;
    const auto finish = [&](Held2Step2Outcome outcome,
                            const std::string& reason) {
        result.outcome = outcome;
        result.reason = reason;
        result.timing.invocation_count = 1;
        result.timing.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start
        ).count();
        result.timing.cpu_seconds = static_cast<double>(
            std::clock() - cpu_start
        ) / static_cast<double>(CLOCKS_PER_SEC);
        result.timing.terminal_status =
            outcome == Held2Step2Outcome::Indeterminate
            ? "indeterminate"
            : "complete";
        result.timing.terminal_reason = reason;
        result.timing.next_step =
            outcome == Held2Step2Outcome::NegativeWitness ? 3 : 0;
        return result;
    };
    const auto record_minimum = [&](double tpd) {
        if (std::isfinite(tpd)
            && (!result.minimum_tpd || tpd < *result.minimum_tpd)) {
            result.minimum_tpd = tpd;
        }
    };
    if (step1.status != "complete" || !step1.coordinates
        || !step1.independent_feed || !step1.volume_bounds
        || !evaluator || provider_evaluation_budget < 2) {
        return finish(
            Held2Step2Outcome::Indeterminate, "invalid_step2_input"
        );
    }
    const std::vector<double>& feed = *step1.independent_feed;
    const Held2StateEvaluator counted = [&](const auto& composition,
                                             double log_volume) {
        ++result.timing.provider_evaluations;
        return evaluator(composition, log_volume);
    };
    const auto enumerate_pressure_roots = [&](
        const std::vector<double>& composition
    ) {
        ++result.timing.provider_evaluations;
        const std::array<double, 2> physical_bounds =
            (*step1.volume_bounds)(composition);
        return evaluate_held2_pressure_envelope(
            composition,
            physical_bounds,
            counted,
            kPressureIntervals,
            kPressureSubdivisionDepth
        );
    };

    Held2ProgressEvent progress;
    progress.kind = Held2ProgressKind::ReferenceStart;
    observe_held2(observer, progress);
    try {
        result.reference_envelope = enumerate_pressure_roots(feed);
    } catch (...) {
        return finish(
            Held2Step2Outcome::Indeterminate,
            "reference_evaluation_failed"
        );
    }
    for (std::size_t index = 0;
         index < result.reference_envelope->roots.size();
         ++index) {
        const Held2PressureRoot& root =
            result.reference_envelope->roots[index];
        progress = {};
        progress.kind = Held2ProgressKind::ReferenceRoot;
        progress.count = static_cast<int>(index);
        progress.volume = root.volume;
        progress.objective = root.objective;
        progress.pressure_residual = root.pressure_residual;
        progress.mechanical_class = root.mechanical_class;
        observe_held2(observer, progress);
    }
    const Held2PressureEnvelopeResult& reference_envelope =
        *result.reference_envelope;
    if (reference_envelope.outcome != "selected"
        || reference_envelope.selected_root_index < 0) {
        return finish(
            Held2Step2Outcome::Indeterminate,
            reference_envelope.failure_reason.empty()
            ? "reference_root_indeterminate"
            : reference_envelope.failure_reason
        );
    }
    const Held2PressureRoot& selected = reference_envelope.roots[
        static_cast<std::size_t>(reference_envelope.selected_root_index)
    ];
    for (const Held2PressureScanPoint& point :
         reference_envelope.scan_points) {
        const bool boundary = audit_held2_tolerance(
            kHeld2RootBoundary,
            point.log_volume - reference_envelope.lower_log_volume
        ).passed || audit_held2_tolerance(
            kHeld2RootBoundary,
            reference_envelope.upper_log_volume - point.log_volume
        ).passed;
        if (!boundary) {
            continue;
        }
        if (!point.valid) {
            return finish(
                Held2Step2Outcome::Indeterminate,
                "reference_boundary_evaluation_failed"
            );
        }
        const double scale = std::max(
            std::abs(point.objective), std::abs(selected.objective)
        );
        if (point.objective < selected.objective
            || audit_held2_tolerance(
                kHeld2StableObjectiveTie,
                point.objective - selected.objective,
                scale
            ).passed) {
            return finish(
                Held2Step2Outcome::Indeterminate,
                "reference_boundary_truncation"
            );
        }
    }
    result.reference = selected.state;
    if (!audit_held2_tolerance(
            kHeld2TpdReferenceZero,
            evaluate_held2_tpd(
                *result.reference, feed, *result.reference, feed
            ).value
        ).passed) {
        return finish(
            Held2Step2Outcome::Indeterminate, "reference_tpd_nonzero"
        );
    }
    progress = {};
    progress.kind = Held2ProgressKind::ReferenceSelected;
    progress.count = reference_envelope.selected_root_index;
    progress.status = "selected";
    observe_held2(observer, progress);
    const std::size_t dimension =
        step1.coordinates->independent_indices.size();
    const auto volume_bounds = [&](const auto& independent) {
        ++result.timing.provider_evaluations;
        return (*step1.volume_bounds)(independent);
    };
    const auto evaluate_at = [&](
        const auto& independent,
        double normalized_volume,
        const std::array<double, 2>& physical_bounds
    ) {
        Held2StageIJointEvaluation evaluation;
        evaluation.independent_modified_fractions = independent;
        evaluation.normalized_volume = normalized_volume;
        try {
            if (independent.size() != dimension
                || !std::isfinite(normalized_volume)
                || normalized_volume < 0.0
                || normalized_volume > 1.0
                || !std::isfinite(physical_bounds[0])
                || !std::isfinite(physical_bounds[1])
                || physical_bounds[0] <= 0.0
                || physical_bounds[1] <= physical_bounds[0]) {
                throw std::invalid_argument("invalid Step-2 search point");
            }
            const double lower_log_volume = std::log(physical_bounds[0]);
            const double log_volume = lower_log_volume
                + normalized_volume
                    * (std::log(physical_bounds[1]) - lower_log_volume);
            evaluation.state = counted(independent, log_volume);
            if (!std::isfinite(evaluation.state.volume)
                || evaluation.state.volume <= 0.0
                || !audit_held2_tolerance(
                    kHeld2JointVolumeConsistency,
                    std::log(evaluation.state.volume) - log_volume
                ).passed
                || !std::isfinite(
                    evaluation.state.pressure_stationarity_relative
                )) {
                throw std::invalid_argument(
                    "incomplete Step-2 joint state"
                );
            }
            const std::vector<double> physical =
                held2_lift_independent_fractions(
                    *step1.coordinates, independent
                );
            evaluation.physical_total_ion_mole_fraction = 0.0;
            for (std::size_t index = 0; index < physical.size(); ++index) {
                if (step1.coordinates->charges[index] != 0.0) {
                    evaluation.physical_total_ion_mole_fraction +=
                        physical[index];
                }
            }
            evaluation.total_ion_mole_fraction_max =
                step1.total_ion_mole_fraction_max;
            if (std::isfinite(step1.total_ion_mole_fraction_max)
                && !audit_held2_tolerance(
                    kHeld2PolytopeFeasibility,
                    std::max(
                        0.0,
                        evaluation.physical_total_ion_mole_fraction
                            - step1.total_ion_mole_fraction_max
                    )
                ).passed) {
                throw std::invalid_argument(
                    "Step-2 joint state exceeds Provider ion domain"
                );
            }
            evaluation.tpd = evaluate_held2_tpd(
                *result.reference,
                feed,
                evaluation.state,
                evaluation.independent_modified_fractions
            ).value;
        } catch (...) {}
        return evaluation;
    };
    const auto evaluate_independent = [&](
        const auto& independent,
        double normalized_volume
    ) {
        try {
            return evaluate_at(
                independent,
                normalized_volume,
                volume_bounds(independent)
            );
        } catch (...) {
            Held2StageIJointEvaluation evaluation;
            evaluation.independent_modified_fractions = independent;
            evaluation.normalized_volume = normalized_volume;
            return evaluation;
        }
    };
    const auto recertify = [&](const Held2StageIJointEvaluation& candidate) {
        return evaluate_independent(
            candidate.independent_modified_fractions,
            candidate.normalized_volume
        );
    };
    const auto refine_stationary_witness = [&](
        const Held2StageIJointEvaluation& certified
    ) {
        if (audit_held2_tolerance(
                kHeld2RootPressure,
                pressure_residual(certified)
            ).passed) {
            return certified;
        }
        try {
            const Held2PressureEnvelopeResult envelope =
                enumerate_pressure_roots(
                    certified.independent_modified_fractions
                );
            if (envelope.outcome != "selected"
                || envelope.selected_root_index < 0) {
                return certified;
            }
            const Held2PressureRoot& root = envelope.roots[
                static_cast<std::size_t>(
                    envelope.selected_root_index
                )
            ];
            const std::array<double, 2> bounds = volume_bounds(
                certified.independent_modified_fractions
            );
            const double lower_log_volume = std::log(bounds[0]);
            Held2StageIJointEvaluation stationary = evaluate_at(
                certified.independent_modified_fractions,
                (std::log(root.volume) - lower_log_volume)
                    / (std::log(bounds[1]) - lower_log_volume),
                bounds
            );
            if (is_strict_negative(stationary)) {
                return stationary;
            }
        } catch (...) {}
        return certified;
    };
    const auto make_witness = [&](const Held2StageIJointEvaluation& witness) {
        return Held2StageICandidate{
            held2_transform_physical_fractions(
                *step1.coordinates,
                held2_lift_independent_fractions(
                    *step1.coordinates,
                    witness.independent_modified_fractions
                )
            ),
            witness.state.volume,
            witness.tpd,
            witness.state.objective,
            witness.state.gradient,
        };
    };
    const auto accept_witness = [&](
        const Held2StageIJointEvaluation& witness,
        const char* reason
    ) {
        record_minimum(witness.tpd);
        result.negative_witness = make_witness(witness);
        return finish(Held2Step2Outcome::NegativeWitness, reason);
    };
    std::vector<double> trace = feed;
    bool has_independent_charge = false;
    for (std::size_t index = 0; index < dimension; ++index) {
        const std::size_t component =
            step1.coordinates->independent_indices[index];
        if (step1.coordinates->charges[component] == 0.0) {
            continue;
        }
        has_independent_charge = true;
        trace[index] = std::nextafter(
            step1.coordinates->independent_lower_bounds[index],
            step1.coordinates->independent_upper_bounds[index]
        );
    }
    if (has_independent_charge) {
        constexpr int kTraceVolumeIntervals = 128;
        std::array<double, 2> trace_bounds;
        try {
            trace_bounds = volume_bounds(trace);
        } catch (...) {
            return finish(
                Held2Step2Outcome::Indeterminate,
                "trace_boundary_evaluation_failed"
            );
        }
        for (int sample = 0; sample <= kTraceVolumeIntervals; ++sample) {
            const Held2StageIJointEvaluation boundary = evaluate_at(
                trace,
                static_cast<double>(sample)
                    / static_cast<double>(kTraceVolumeIntervals),
                trace_bounds
            );
            if (!is_valid(boundary)) {
                return finish(
                    Held2Step2Outcome::Indeterminate,
                    "trace_boundary_evaluation_failed"
                );
            }
            record_minimum(boundary.tpd);
            progress = {};
            progress.kind = Held2ProgressKind::StageIEvaluation;
            progress.iteration = sample + 1;
            progress.objective = boundary.tpd;
            progress.status = "admissible_trace_boundary";
            progress.volume = boundary.state.volume;
            progress.pressure_residual = pressure_residual(boundary);
            observe_held2(observer, progress);
            if (!is_strict_negative(boundary)) {
                continue;
            }
            const Held2StageIJointEvaluation certified =
                recertify(boundary);
            if (!is_strict_negative(certified)) {
                return finish(
                    Held2Step2Outcome::Indeterminate,
                    "negative_witness_recertification_failed"
                );
            }
            const Held2StageIJointEvaluation witness =
                refine_stationary_witness(certified);
            return accept_witness(
                witness,
                "trace_boundary_negative_tpd_witness"
            );
        }
    }
    const auto assess_search = [&](const Held2StageIDirectResult& search) {
        result.timing.optimizer_iterations +=
            static_cast<std::uint64_t>(search.evaluation_count);
        record_minimum(search.minimum_tpd);
        Held2StageIJointEvaluation assessment;
        if (search.negative_candidate) {
            const Held2StageIJointEvaluation certified = recertify(
                *search.negative_candidate
            );
            if (!is_strict_negative(certified)) {
                assessment.failure_reason =
                    "negative_witness_recertification_failed";
                return assessment;
            }
            return certified;
        }
        assessment.failure_reason = search.failure_reason;
        return assessment;
    };
    const Held2StageIJointEvaluator complete_joint = [
        &evaluate_independent,
        &step1
    ](const auto& cube) {
        const std::vector<double> composition_cube(
            cube.begin() + 1,
            cube.end()
        );
        return evaluate_independent(
            held2_map_unit_cube_to_independent_fractions(
                *step1.coordinates,
                composition_cube,
                step1.total_ion_mole_fraction_max
            ),
            cube.front()
        );
    };
    ++result.timing.optimizer_solves;
    const Held2StageIJointEvaluation complete_witness = assess_search(
        run_direct_search(
            dimension + 1,
            provider_evaluation_budget,
            2,
            complete_joint,
            observer
        )
    );
    if (is_valid(complete_witness)) {
        constexpr int kPreparationVolumeSamples = 65;
        constexpr int kPreparationRootRefinementLimit = 64;
        constexpr int kPreparationProviderEvaluations =
            1 + kPreparationVolumeSamples
            + kPreparationRootRefinementLimit;
        const int joint_provider_evaluations =
            static_cast<int>(result.timing.optimizer_iterations) * 2;
        const int preparation_budget =
            provider_evaluation_budget - joint_provider_evaluations;
        Held2StageIJointEvaluation prepared_witness;
        if (preparation_budget >= kPreparationProviderEvaluations) {
            constexpr double kPromisingPreparationTpd = 0.2;
            const Held2StageIJointEvaluator prepare_composition = [
                &evaluate_at,
                &volume_bounds,
                &step1
            ](const auto& cube) {
                Held2StageIJointEvaluation best;
                try {
                    const std::vector<double> independent =
                        held2_map_unit_cube_to_independent_fractions(
                            *step1.coordinates,
                            cube,
                            step1.total_ion_mole_fraction_max
                        );
                    const std::array<double, 2> bounds =
                        volume_bounds(independent);
                    std::optional<std::pair<
                        Held2StageIJointEvaluation,
                        Held2StageIJointEvaluation
                    >> promising_bracket;
                    Held2StageIJointEvaluation previous;
                    const auto brackets_zero = [](double left,
                                                  double right) {
                        return (left <= 0.0 && right >= 0.0)
                            || (left >= 0.0 && right <= 0.0);
                    };
                    for (int sample = 0;
                         sample < kPreparationVolumeSamples;
                         ++sample) {
                        const Held2StageIJointEvaluation evaluation =
                            evaluate_at(
                                independent,
                                static_cast<double>(sample)
                                    / static_cast<double>(
                                        kPreparationVolumeSamples - 1
                                    ),
                                bounds
                            );
                        if (!is_valid(evaluation)) {
                            return evaluation;
                        }
                        if (evaluation.tpd < best.tpd) {
                            best = evaluation;
                        }
                        if (sample > 0
                            && brackets_zero(
                                pressure_residual(previous),
                                pressure_residual(evaluation)
                            )
                            && (!promising_bracket
                                || std::min(previous.tpd, evaluation.tpd)
                                    < std::min(
                                        promising_bracket->first.tpd,
                                        promising_bracket->second.tpd
                                    ))) {
                            promising_bracket = {previous, evaluation};
                        }
                        previous = evaluation;
                    }
                    if (best.tpd <= kPromisingPreparationTpd
                        && promising_bracket) {
                        auto left = promising_bracket->first;
                        auto right = promising_bracket->second;
                        auto root = std::abs(pressure_residual(left))
                                <= std::abs(pressure_residual(right))
                            ? left : right;
                        for (int iteration = 0;
                             iteration < kPreparationRootRefinementLimit;
                             ++iteration) {
                            if (audit_held2_tolerance(
                                    kHeld2RootPressure,
                                    pressure_residual(root)
                                ).passed) {
                                break;
                            }
                            const double midpoint = 0.5 * (
                                left.normalized_volume
                                + right.normalized_volume
                            );
                            const auto middle = evaluate_at(
                                independent, midpoint, bounds
                            );
                            if (!is_valid(middle)) {
                                return middle;
                            }
                            if (std::abs(pressure_residual(middle))
                                < std::abs(pressure_residual(root))) {
                                root = middle;
                            }
                            if (brackets_zero(
                                    pressure_residual(left),
                                    pressure_residual(middle)
                                )) {
                                right = middle;
                            } else {
                                left = middle;
                            }
                        }
                        if (is_strict_stable_root(root)
                            && is_strict_negative(root)) {
                            return root;
                        }
                    }
                } catch (...) {}
                return best;
            };
            ++result.timing.optimizer_solves;
            const Held2StageIDirectResult preparation_search =
                run_direct_search(
                    dimension,
                    preparation_budget,
                    kPreparationProviderEvaluations,
                    prepare_composition,
                    observer
                );
            prepared_witness = assess_search(
                preparation_search
            );
            if (!prepared_witness.failure_reason.empty()) {
                return finish(
                    Held2Step2Outcome::Indeterminate,
                    prepared_witness.failure_reason
                );
            }
        }
        const Held2StageIJointEvaluation witness =
            refine_stationary_witness(
                is_valid(prepared_witness)
                ? prepared_witness
                : complete_witness
            );
        return accept_witness(witness, "negative_tpd_witness");
    }
    if (!complete_witness.failure_reason.empty()) {
        return finish(
            Held2Step2Outcome::Indeterminate,
            complete_witness.failure_reason
        );
    }
    return finish(
        Held2Step2Outcome::NoNegativeWitnessDetected,
        "declared_search_complete"
    );
}

}  // namespace epcsaft_equilibrium
