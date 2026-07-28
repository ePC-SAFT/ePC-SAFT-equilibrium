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

struct Held2StageIReducedEvaluation {
    std::vector<double> chart_coordinates;
    std::vector<double> independent_modified_fractions;
    Held2PressureEnvelopeResult pressure_envelope;
    double trial_volume = 0.0;
    double trial_pressure_residual = 0.0;
    double physical_total_ion_mole_fraction =
        std::numeric_limits<double>::quiet_NaN();
    double total_ion_mole_fraction_max =
        std::numeric_limits<double>::quiet_NaN();
    double tpd = std::numeric_limits<double>::infinity();
    bool certified = false;
    std::string failure_reason;
};

using Held2StageIReducedEvaluator = std::function<Held2StageIReducedEvaluation(
    const std::vector<double>&
)>;

struct Held2StageIDirectResult {
    std::string outcome = "indeterminate";
    std::string termination_reason;
    int failed_evaluation_count = 0;
    int negative_witness_index = -1;
    double minimum_tpd = std::numeric_limits<double>::infinity();
    std::vector<Held2StageIReducedEvaluation> evaluations;
};

struct DirectContext {
    const Held2StageIReducedEvaluator* evaluator = nullptr;
    Held2StageIDirectResult* result = nullptr;
    nlopt::opt* optimizer = nullptr;
    Held2ProgressObserver* observer = nullptr;
    bool stop_requested = false;
};

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
    Held2StageIReducedEvaluation evaluation =
        (*context.evaluator)(chart_coordinates);
    evaluation.chart_coordinates = chart_coordinates;
    const int evaluation_index = static_cast<int>(context.result->evaluations.size());
    context.result->evaluations.push_back(std::move(evaluation));
    const Held2StageIReducedEvaluation& retained =
        context.result->evaluations.back();
    Held2ProgressEvent progress;
    progress.kind = Held2ProgressKind::StageIEvaluation;
    progress.iteration = evaluation_index + 1;
    progress.objective = retained.tpd;
    progress.physical_total_ion_mole_fraction =
        retained.physical_total_ion_mole_fraction;
    progress.total_ion_mole_fraction_max =
        retained.total_ion_mole_fraction_max;
    progress.status = retained.certified ? "certified" : "failed";
    progress.reason = retained.failure_reason;
    if (retained.pressure_envelope.selected_root_index >= 0) {
        const Held2PressureRoot& selected = retained.pressure_envelope.roots[
            static_cast<std::size_t>(
                retained.pressure_envelope.selected_root_index
            )
        ];
        progress.volume = selected.volume;
        progress.pressure_residual = selected.pressure_residual;
    } else if (retained.trial_volume > 0.0) {
        progress.volume = retained.trial_volume;
        progress.pressure_residual = retained.trial_pressure_residual;
    }
    observe_held2(context.observer, progress);
    if (!retained.certified || !std::isfinite(retained.tpd)) {
        Held2ProgressEvent failure;
        failure.kind = Held2ProgressKind::Failure;
        failure.stage = "STAGE I";
        failure.reason = retained.failure_reason.empty()
            ? "required_envelope_evaluation_failed"
            : retained.failure_reason;
        observe_held2(context.observer, failure);
        ++context.result->failed_evaluation_count;
        context.result->termination_reason =
            "required_envelope_evaluation_failed";
        context.stop_requested = true;
        context.optimizer->force_stop();
        return 0.0;
    }
    context.result->minimum_tpd = std::min(
        context.result->minimum_tpd,
        retained.tpd
    );
    if (audit_held2_tolerance(kHeld2TpdNegativeMargin, retained.tpd).passed) {
        context.result->negative_witness_index = evaluation_index;
        context.result->termination_reason = "certified_negative_tpd";
        context.stop_requested = true;
        context.optimizer->force_stop();
        return retained.tpd;
    }
    return retained.tpd;
}

Held2StageIDirectResult run_direct_search(
    std::size_t composition_dimension,
    int evaluation_budget,
    const Held2StageIReducedEvaluator& evaluator,
    Held2ProgressObserver* observer
) {
    Held2StageIDirectResult result;
    nlopt::opt optimizer(nlopt::GN_DIRECT_L, composition_dimension);
    DirectContext context{&evaluator, &result, &optimizer, observer, false};
    optimizer.set_lower_bounds(std::vector<double>(composition_dimension, 0.0));
    optimizer.set_upper_bounds(std::vector<double>(composition_dimension, 1.0));
    optimizer.set_maxeval(evaluation_budget);
    optimizer.set_min_objective(direct_objective, &context);
    std::vector<double> initial(composition_dimension, 0.5);
    double minimum = std::numeric_limits<double>::infinity();
    try {
        const nlopt::result status = optimizer.optimize(initial, minimum);
        if (status == nlopt::MAXEVAL_REACHED
            && result.failed_evaluation_count == 0) {
            result.outcome = "no_negative_witness_detected";
            result.termination_reason = "declared_budget_exhausted";
        } else {
            result.termination_reason = "unexpected_solver_termination";
        }
    } catch (const nlopt::forced_stop&) {
        if (result.negative_witness_index >= 0) {
            result.outcome = "negative_witness_found";
        }
    } catch (const std::exception& error) {
        result.termination_reason = std::string("solver_failure: ") + error.what();
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
    int search_budget,
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
    if (step1.status != "complete" || !step1.coordinates
        || !step1.independent_feed || !step1.volume_bounds
        || !evaluator || search_budget < 1) {
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
    const auto envelope = [&](const std::vector<double>& composition) {
        ++result.timing.provider_evaluations;
        return evaluate_held2_pressure_envelope(
            composition,
            (*step1.volume_bounds)(composition),
            counted,
            kPressureIntervals,
            kPressureSubdivisionDepth
        );
    };

    Held2ProgressEvent progress;
    progress.kind = Held2ProgressKind::ReferenceStart;
    observe_held2(observer, progress);
    try {
        result.reference_envelope = envelope(feed);
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
    const auto evaluate_independent = [&](const auto& independent) {
        Held2StageIReducedEvaluation evaluation;
        evaluation.independent_modified_fractions = independent;
        try {
            if (independent.size() != dimension) {
                throw std::invalid_argument("invalid Step-2 search point");
            }
            evaluation.pressure_envelope = envelope(
                independent
            );
            if (evaluation.pressure_envelope.outcome != "selected"
                || evaluation.pressure_envelope.selected_root_index < 0) {
                throw std::runtime_error(
                    evaluation.pressure_envelope.failure_reason
                );
            }
            const Held2PressureRoot& trial =
                evaluation.pressure_envelope.roots[
                    static_cast<std::size_t>(
                        evaluation.pressure_envelope.selected_root_index
                    )
                ];
            evaluation.trial_volume = trial.volume;
            evaluation.trial_pressure_residual = trial.pressure_residual;
            evaluation.tpd = evaluate_held2_tpd(
                *result.reference,
                feed,
                trial.state,
                evaluation.independent_modified_fractions
            ).value;
            evaluation.certified = std::isfinite(evaluation.tpd);
        } catch (...) {
            evaluation.failure_reason = "trial_evaluation_failed";
        }
        return evaluation;
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
        const Held2StageIReducedEvaluation boundary =
            evaluate_independent(trace);
        if (!boundary.certified || !std::isfinite(boundary.tpd)) {
            return finish(
                Held2Step2Outcome::Indeterminate,
                "trace_boundary_evaluation_failed"
            );
        }
        result.minimum_tpd = boundary.tpd;
        progress = {};
        progress.kind = Held2ProgressKind::StageIEvaluation;
        progress.iteration = 1;
        progress.objective = boundary.tpd;
        progress.status = "certified_trace_boundary";
        progress.volume = boundary.trial_volume;
        progress.pressure_residual = boundary.trial_pressure_residual;
        observe_held2(observer, progress);
        if (audit_held2_tolerance(
                kHeld2TpdNegativeMargin, boundary.tpd
            ).passed) {
            result.negative_witness = Held2StageICandidate{
                held2_transform_physical_fractions(
                    *step1.coordinates,
                    held2_lift_independent_fractions(
                        *step1.coordinates, trace
                    )
                ),
                boundary.trial_volume,
                boundary.tpd,
            };
            return finish(
                Held2Step2Outcome::NegativeWitness,
                "trace_boundary_negative_tpd_witness"
            );
        }
    }
    const Held2StageIReducedEvaluator reduced = [
        &evaluate_independent,
        &step1
    ](const auto& cube) {
        return evaluate_independent(
            held2_map_unit_cube_to_independent_fractions(
                *step1.coordinates,
                cube,
                step1.total_ion_mole_fraction_max
            )
        );
    };
    ++result.timing.optimizer_solves;
    const Held2StageIDirectResult search = run_direct_search(
        dimension,
        search_budget,
        reduced,
        observer
    );
    result.timing.optimizer_iterations =
        static_cast<std::uint64_t>(search.evaluations.size());
    if (std::isfinite(search.minimum_tpd)
        && (!result.minimum_tpd
            || search.minimum_tpd < *result.minimum_tpd)) {
        result.minimum_tpd = search.minimum_tpd;
    }
    if (search.outcome == "negative_witness_found"
        && search.negative_witness_index >= 0) {
        const Held2StageIReducedEvaluation& witness = search.evaluations[
            static_cast<std::size_t>(search.negative_witness_index)
        ];
        result.negative_witness = Held2StageICandidate{
            held2_transform_physical_fractions(
                *step1.coordinates,
                held2_lift_independent_fractions(
                    *step1.coordinates,
                    witness.independent_modified_fractions
                )
            ),
            witness.trial_volume,
            witness.tpd,
        };
        return finish(
            Held2Step2Outcome::NegativeWitness,
            "negative_tpd_witness"
        );
    }
    if (search.outcome == "no_negative_witness_detected") {
        return finish(
            Held2Step2Outcome::NoNegativeWitnessDetected,
            "declared_search_complete"
        );
    }
    return finish(
        Held2Step2Outcome::Indeterminate,
        search.termination_reason.empty()
        ? "step2_search_failed"
        : search.termination_reason
    );
}

}  // namespace epcsaft_equilibrium
