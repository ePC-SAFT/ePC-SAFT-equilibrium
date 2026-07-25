#include "held2_step2.hpp"
#include "held2_tolerances.hpp"

#include <nlopt.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <stdexcept>

namespace epcsaft_equilibrium {
namespace {

std::string nlopt_version_string() {
    int major = 0;
    int minor = 0;
    int bugfix = 0;
    nlopt_version(&major, &minor, &bugfix);
    return std::to_string(major) + "." + std::to_string(minor) + "."
        + std::to_string(bugfix);
}

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
    ++context.result->completed_evaluation_count;
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

Held2StateEvaluation evaluate_simple_manufactured_state(
    double composition,
    double log_volume,
    const std::string& topology
) {
    const double volume = std::exp(log_volume);
    const double volume_delta = volume - 1.0;
    double composition_objective = 0.0;
    double composition_gradient = 0.0;
    if (topology == "no_negative") {
        const double delta = composition - 0.5;
        composition_objective = delta * delta;
        composition_gradient = 2.0 * delta;
    } else if (topology == "narrow_negative") {
        constexpr double center = 0.73;
        constexpr double width = 0.025;
        const double scaled = (composition - center) / width;
        composition_objective = -0.01 * std::exp(-scaled * scaled);
        composition_gradient = -2.0 * scaled / width * composition_objective;
    } else {
        throw std::invalid_argument("unknown simple manufactured Stage-I topology");
    }

    Held2StateEvaluation state;
    state.modified_fractions = {1.0 - composition, composition};
    state.physical_amounts = state.modified_fractions;
    state.volume = volume;
    state.objective = composition_objective + 2.5 * volume_delta * volume_delta;
    const double volume_gradient = 5.0 * volume_delta;
    state.gradient = {composition_gradient, volume * volume_gradient};
    state.hessian = {
        0.0,
        0.0,
        0.0,
        5.0 * volume * volume + volume * volume_gradient,
    };
    state.modified_potentials = {0.0, composition_gradient};
    state.pressure_stationarity_relative = -volume_gradient;
    state.pressure_stationarity_derivative_log_volume = -5.0 * volume;
    return state;
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
    const Held2StageIReducedEvaluator reduced = [&](const auto& cube) {
        Held2StageIReducedEvaluation evaluation;
        try {
            if (cube.size() != dimension) {
                throw std::invalid_argument("invalid Step-2 search point");
            }
            evaluation.independent_modified_fractions =
                held2_map_unit_cube_to_independent_fractions(
                    *step1.coordinates,
                    cube,
                    step1.total_ion_mole_fraction_max
                );
            evaluation.pressure_envelope = envelope(
                evaluation.independent_modified_fractions
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
    ++result.timing.optimizer_solves;
    const Held2StageIDirectResult search = solve_held2_stage_i_direct(
        dimension,
        search_budget,
        -kHeld2TpdNegativeMargin.atol,
        reduced,
        observer
    );
    result.timing.optimizer_iterations =
        static_cast<std::uint64_t>(search.evaluations.size());
    if (std::isfinite(search.minimum_tpd)) {
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

Held2StageIDirectResult solve_held2_stage_i_direct(
    std::size_t composition_dimension,
    int evaluation_budget,
    double negative_tpd_threshold,
    const Held2StageIReducedEvaluator& evaluator,
    Held2ProgressObserver* observer
) {
    if (composition_dimension == 0 || evaluation_budget < 1
        || !std::isfinite(negative_tpd_threshold)
        || negative_tpd_threshold != -kHeld2TpdNegativeMargin.atol) {
        throw std::invalid_argument("HELD2 DIRECT-L search policy is invalid");
    }
    Held2StageIDirectResult result;
    result.declared_evaluation_budget = evaluation_budget;
    result.solver_version = nlopt_version_string();
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

Held2StageIDirectResult solve_held2_manufactured_stage_i_direct(
    const std::string& topology,
    int evaluation_budget
) {
    const Held2Coordinates coordinates = make_held2_coordinates({0.0, 1.0, -1.0});
    const std::vector<double> feed = {0.5};
    const auto manufactured_envelope = [&coordinates, &topology](
        const std::vector<double>& independent
    ) {
        if (topology == "branch_switch") {
            return evaluate_held2_manufactured_pressure_envelope(
                "branch_switch",
                independent.front(),
                64
            );
        }
        if (topology == "boundary" || topology == "provider_failure") {
            return evaluate_held2_manufactured_pressure_envelope(
                topology == "boundary" ? "boundary" : "invalid",
                independent.front(),
                64
            );
        }
        const Held2StateEvaluator phase_evaluator = [&coordinates, &topology](
            const std::vector<double>& composition,
            double log_volume
        ) {
            if (topology == "negative") {
                Held2StateEvaluation state = evaluate_held2_manufactured_state(
                    coordinates,
                    composition,
                    log_volume
                );
                state.pressure_stationarity_relative *= -1.0;
                state.pressure_stationarity_derivative_log_volume *= -1.0;
                return state;
            }
            return evaluate_simple_manufactured_state(
                composition.front(),
                log_volume,
                topology
            );
        };
        return evaluate_held2_pressure_envelope(
            independent,
            {0.5, 1.5},
            phase_evaluator,
            64,
            8
        );
    };

    Held2PressureEnvelopeResult reference_envelope = manufactured_envelope(feed);
    Held2StateEvaluation reference;
    if (reference_envelope.outcome == "selected") {
        reference = reference_envelope.roots[static_cast<std::size_t>(
            reference_envelope.selected_root_index
        )].state;
    }
    const Held2StageIReducedEvaluator evaluator = [
        &coordinates,
        &topology,
        &manufactured_envelope,
        &reference,
        feed
    ](const std::vector<double>& chart_coordinates) {
        Held2StageIReducedEvaluation evaluation;
        evaluation.independent_modified_fractions =
            held2_map_unit_cube_to_independent_fractions(
                coordinates,
                chart_coordinates,
                std::numeric_limits<double>::quiet_NaN()
            );
        evaluation.pressure_envelope = manufactured_envelope(
            evaluation.independent_modified_fractions
        );
        if (evaluation.pressure_envelope.outcome != "selected") {
            evaluation.failure_reason =
                evaluation.pressure_envelope.failure_reason;
            return evaluation;
        }
        const Held2PressureRoot& selected = evaluation.pressure_envelope.roots[
            static_cast<std::size_t>(
                evaluation.pressure_envelope.selected_root_index
            )
        ];
        if (topology == "branch_switch") {
            evaluation.tpd = 10.0 + selected.objective;
        } else {
            evaluation.tpd = selected.objective - reference.objective;
            for (std::size_t index = 0; index < feed.size(); ++index) {
                evaluation.tpd -= reference.gradient[index]
                    * (evaluation.independent_modified_fractions[index] - feed[index]);
            }
        }
        evaluation.certified = true;
        return evaluation;
    };
    return solve_held2_stage_i_direct(
        coordinates.independent_indices.size(),
        evaluation_budget,
        -kHeld2TpdNegativeMargin.atol,
        evaluator
    );
}

}  // namespace epcsaft_equilibrium
