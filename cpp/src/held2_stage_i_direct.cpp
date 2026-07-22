#include "held2_stage_i_direct.hpp"
#include "held2_tolerances.hpp"

#include <nlopt.hpp>

#include <algorithm>
#include <cmath>
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
    if (!retained.certified || !std::isfinite(retained.tpd)) {
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

Held2StageIDirectResult solve_held2_stage_i_direct(
    std::size_t composition_dimension,
    int evaluation_budget,
    double negative_tpd_threshold,
    const Held2StageIReducedEvaluator& evaluator
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
    DirectContext context{&evaluator, &result, &optimizer, false};
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
                chart_coordinates
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
