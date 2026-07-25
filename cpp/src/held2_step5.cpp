#include "held2_step5.hpp"
#include "held2_tolerances.hpp"

#include <nlopt.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
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

double maximum_abs_difference(
    const std::vector<double>& left,
    const std::vector<double>& right
) {
    if (left.size() != right.size()) {
        return std::numeric_limits<double>::infinity();
    }
    double maximum = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        maximum = std::max(maximum, std::abs(left[index] - right[index]));
    }
    return maximum;
}

bool same_composition(
    const std::vector<double>& left,
    const std::vector<double>& right
) {
    return maximum_abs_difference(left, right)
        <= kHeld2BasinDuplicateComposition.atol;
}

bool same_physical_start(
    const Held2StageIIPhysicalStart& left,
    const Held2StageIIPhysicalStart& right
) {
    return same_composition(
               left.independent_modified_fractions,
               right.independent_modified_fractions
           )
        && std::abs(left.log_volume - right.log_volume)
            <= kHeld2BasinDuplicateLogVolume.atol;
}

Held2StageIIBasinEvaluation evaluate_fail_closed(
    const Held2StageIIBasinEvaluator& evaluator,
    const std::vector<double>& independent
) {
    try {
        return evaluator(independent);
    } catch (const std::exception& error) {
        Held2StageIIBasinEvaluation failed;
        failed.independent_modified_fractions = independent;
        failed.failure_reason = std::string("envelope_evaluator_exception: ")
            + error.what();
        return failed;
    }
}

void retain_evaluation(
    Held2StageIIBasinExplorationResult& result,
    Held2StageIIBasinEvaluation evaluation,
    const std::string& source
) {
    if (!evaluation.certified
        || !std::isfinite(evaluation.reduced_lower_value)) {
        ++result.failed_evaluation_count;
        result.evaluations.push_back(std::move(evaluation));
        return;
    }
    ++result.completed_evaluation_count;
    double minimum_stable_root_objective =
        std::numeric_limits<double>::infinity();
    for (const Held2PressureRoot& root : evaluation.pressure_envelope.roots) {
        if (root.mechanical_class == "strict_stable" && !root.boundary) {
            minimum_stable_root_objective = std::min(
                minimum_stable_root_objective,
                root.objective
            );
        }
    }
    const double composition_offset = evaluation.reduced_lower_value
        - minimum_stable_root_objective;
    int stable_branch_index = 0;
    for (const Held2PressureRoot& root : evaluation.pressure_envelope.roots) {
        if (root.mechanical_class != "strict_stable" || root.boundary) {
            continue;
        }
        Held2StageIIPhysicalStart start;
        start.independent_modified_fractions =
            evaluation.independent_modified_fractions;
        start.stable_branch_index = stable_branch_index++;
        start.log_volume = root.log_volume;
        start.volume = root.volume;
        start.reduced_lower_value = root.objective + composition_offset;
        start.source = source;
        start.root_origin = root.origin;
        start.root_completeness =
            evaluation.pressure_envelope.root_completeness;
        const bool duplicate = std::any_of(
            result.representatives.begin(),
            result.representatives.end(),
            [&start](const Held2StageIIPhysicalStart& known) {
                return same_physical_start(start, known);
            }
        );
        if (duplicate) {
            ++result.duplicate_start_count;
        } else {
            result.representatives.push_back(std::move(start));
        }
    }
    result.evaluations.push_back(std::move(evaluation));
}

struct DirectContext {
    const Held2Coordinates* coordinates = nullptr;
    const Held2StageIIBasinEvaluator* evaluator = nullptr;
    Held2StageIIBasinExplorationResult* result = nullptr;
    nlopt::opt* optimizer = nullptr;
    double total_ion_mole_fraction_max =
        std::numeric_limits<double>::quiet_NaN();
    bool stop_requested = false;
};

double direct_objective(
    const std::vector<double>& cube,
    std::vector<double>& gradient,
    void* opaque
) {
    auto& context = *static_cast<DirectContext*>(opaque);
    if (context.stop_requested) {
        return 0.0;
    }
    if (!gradient.empty()) {
        throw std::invalid_argument(
            "HELD2 Stage-II DIRECT-L requested an unexpected gradient"
        );
    }
    const std::vector<double> independent =
        held2_map_unit_cube_to_independent_fractions(
            *context.coordinates,
            cube,
            context.total_ion_mole_fraction_max
        );
    Held2StageIIBasinEvaluation evaluation = evaluate_fail_closed(
        *context.evaluator,
        independent
    );
    const double objective = evaluation.reduced_lower_value;
    retain_evaluation(*context.result, std::move(evaluation), "direct_l");
    if (context.result->failed_evaluation_count != 0) {
        context.stop_requested = true;
        context.result->termination_reason =
            "required_envelope_evaluation_failed";
        context.optimizer->force_stop();
        return 0.0;
    }
    return objective;
}

struct SobolParameters {
    int degree;
    std::uint32_t coefficient;
    std::array<std::uint32_t, 5> initial;
};

constexpr std::array<SobolParameters, 9> kSobolParameters = {{
    {1, 0, {1, 0, 0, 0, 0}},
    {2, 1, {1, 3, 0, 0, 0}},
    {3, 1, {1, 3, 1, 0, 0}},
    {3, 2, {1, 1, 1, 0, 0}},
    {4, 1, {1, 3, 5, 13, 0}},
    {4, 4, {1, 1, 5, 5, 0}},
    {5, 2, {1, 3, 3, 9, 7}},
    {5, 4, {1, 1, 5, 11, 27}},
    {5, 7, {1, 1, 7, 13, 3}},
}};

std::array<std::uint32_t, 32> sobol_directions(std::size_t dimension) {
    std::array<std::uint32_t, 32> directions{};
    if (dimension == 0) {
        for (std::size_t bit = 0; bit < directions.size(); ++bit) {
            directions[bit] = std::uint32_t{1} << (31U - bit);
        }
        return directions;
    }
    if (dimension > kSobolParameters.size()) {
        throw std::invalid_argument(
            "HELD2 Sobol explorer supports at most ten composition coordinates"
        );
    }
    const SobolParameters& parameters = kSobolParameters[dimension - 1];
    for (int bit = 0; bit < parameters.degree; ++bit) {
        directions[static_cast<std::size_t>(bit)] =
            parameters.initial[static_cast<std::size_t>(bit)]
            << (31 - bit);
    }
    for (int bit = parameters.degree; bit < 32; ++bit) {
        std::uint32_t value = directions[static_cast<std::size_t>(
            bit - parameters.degree
        )];
        value ^= value >> parameters.degree;
        for (int offset = 1; offset < parameters.degree; ++offset) {
            const int coefficient_bit = parameters.degree - 1 - offset;
            if ((parameters.coefficient >> coefficient_bit) & 1U) {
                value ^= directions[static_cast<std::size_t>(bit - offset)];
            }
        }
        directions[static_cast<std::size_t>(bit)] = value;
    }
    return directions;
}

}  // namespace

Held2Step5Result run_held2_step5(
    const Held2Step1Result& step1,
    const Held2Step4Result& step4,
    Held2PersistentState& state,
    const Held2StateEvaluator& evaluator,
    const Held2PackingFractionEvaluator& packing_fraction,
    const Held2ResourceProfile& resources,
    Held2ProgressObserver* observer
) {
    Held2Step5Result result;
    result.timing.invocation_count = 1;
    if (step4.status != "complete" || !step1.coordinates
        || !step1.volume_bounds || !evaluator || !packing_fraction
        || resources.step5_total_start_cap < 1) {
        result.reason = "invalid_step5_input";
        return result;
    }
    const Held2Coordinates& coordinates = *step1.coordinates;
    const std::size_t dimension = state.feed.size();
    std::vector<double> lower = coordinates.independent_lower_bounds;
    std::vector<double> upper = coordinates.independent_upper_bounds;
    const auto dependent = std::find(
        coordinates.retained_indices.begin(),
        coordinates.retained_indices.end(),
        coordinates.dependent_index
    );
    const double sum_upper = 1.0 - kHeld2ModifiedLowerScale
        * coordinates.modified_factors[static_cast<std::size_t>(
            dependent - coordinates.retained_indices.begin()
        )];
    const std::vector<std::vector<double>> sobol = held2_sobol_points(
        dimension, resources.step5_total_start_cap
    );
    const std::size_t begin = std::min(
        static_cast<std::size_t>(state.next_start_ordinal), sobol.size()
    );
    const std::size_t cap = std::min(
        sobol.size(),
        begin + static_cast<std::size_t>(
            std::max(1, resources.step5_start_epoch_size)
        )
    );
    state.starts_consumed_in_epoch = 0;
    state.start_epoch_added_member = false;
    ++state.start_epoch_index;
    double best = std::numeric_limits<double>::infinity();
    Held2MPoint best_point;
    for (std::size_t index = begin; index < cap; ++index) {
        const std::uint64_t ordinal = state.next_start_ordinal++;
        ++result.starts_consumed;
        const std::vector<double> start =
            held2_map_unit_cube_to_independent_fractions(
                coordinates,
                sobol[index],
                step1.total_ion_mole_fraction_max
            );
        Held2LocalCertificate certificate;
        certificate.start_ordinal = ordinal;
        Held2PressureEnvelopeResult envelope;
        try {
            envelope = evaluate_held2_pressure_envelope(
                start,
                (*step1.volume_bounds)(start),
                evaluator,
                64
            );
        } catch (...) {
            certificate.solver_status = "pressure_root_failed";
            result.attempts.push_back(certificate);
            continue;
        }
        ++result.timing.provider_evaluations;
        if (envelope.outcome != "selected") {
            certificate.solver_status = envelope.failure_reason;
            result.attempts.push_back(certificate);
            continue;
        }
        const int selected = envelope.selected_root_index;
        int branch = 0;
        for (int root = 0; root < selected; ++root) {
            branch += envelope.roots[static_cast<std::size_t>(root)]
                .mechanical_class == "strict_stable" ? 1 : 0;
        }
        const Held2PressureRoot& initial_root =
            envelope.roots[static_cast<std::size_t>(selected)];
        const Held2LocalSearchRun run =
            run_held2_local_pressure_root_search(
                evaluator,
                *step1.volume_bounds,
                state.feed,
                state.multipliers,
                start,
                branch,
                initial_root.log_volume,
                lower,
                upper,
                sum_upper,
                observer,
                state.major_iteration,
                static_cast<int>(ordinal)
            );
        ++result.timing.optimizer_solves;
        result.timing.optimizer_iterations +=
            static_cast<std::uint64_t>(run.optimizer_iterations);
        certificate.solver_status = run.solver_status;
        if (!run.solver_converged || !run.callback_error.empty()
            || run.variables.size() != dimension + 1) {
            result.attempts.push_back(certificate);
            continue;
        }
        const std::vector<double> independent(
            run.variables.begin(), run.variables.end() - 1
        );
        const Held2StateEvaluation terminal =
            evaluator(independent, run.variables.back());
        ++result.timing.provider_evaluations;
        const std::array<double, 2> volume_bounds =
            (*step1.volume_bounds)(independent);
        std::vector<double> full_lower = lower;
        std::vector<double> full_upper = upper;
        full_lower.push_back(std::log(volume_bounds[0]));
        full_upper.push_back(std::log(volume_bounds[1]));
        const Held2PhysicalKkt kkt = certify_held2_local_physical_kkt(
            terminal.gradient,
            state.multipliers,
            run.variables,
            full_lower,
            full_upper,
            sum_upper,
            run.coordinate_jacobian,
            run.lower_bound_multipliers,
            run.upper_bound_multipliers
        );
        double value = terminal.objective;
        for (std::size_t coordinate = 0;
             coordinate < dimension;
             ++coordinate) {
            value += state.multipliers[coordinate]
                * (state.feed[coordinate] - independent[coordinate]);
        }
        certificate.finite_and_in_domain =
            std::isfinite(value) && terminal.volume > 0.0;
        certificate.pressure_residual =
            terminal.pressure_stationarity_relative;
        certificate.primal_residual_inf = kkt.primal_inf_norm;
        certificate.stationarity_residual_inf = kkt.stationarity_inf_norm;
        certificate.dual_sign_violation_inf =
            kkt.dual_sign_violation_inf_norm;
        certificate.complementarity_inf = kkt.complementarity;
        certificate.dual_pullback_residual_inf =
            kkt.reconstruction_inf_norm;
        certificate.accepted = certificate.finite_and_in_domain
            && audit_held2_tolerance(
                kHeld2RootPressure, certificate.pressure_residual
            ).passed
            && audit_held2_tolerance(
                kHeld2Stage2Primal, kkt.primal_inf_norm
            ).passed
            && kkt.dual_signs_valid
            && audit_held2_tolerance(
                kHeld2Stage2Stationarity, kkt.stationarity_inf_norm
            ).passed
            && audit_held2_tolerance(
                kHeld2Stage2Complementarity, kkt.complementarity
            ).passed
            && audit_held2_tolerance(
                kHeld2Stage2DualPullback,
                kkt.reconstruction_inf_norm,
                kkt.reconstruction_scale
            ).passed;
        result.attempts.push_back(certificate);
        if (!assess_held2_stage_ii_step5(
                state.upper_bound, value, certificate.accepted
            ).qualified || value >= best) {
            continue;
        }
        best = value;
        best_point = {
            static_cast<std::uint64_t>(state.M.size()),
            independent,
            terminal.volume,
            packing_fraction(independent, terminal.volume),
            terminal.objective,
            "step5_local",
        };
        break;
    }
    if (!std::isfinite(best)) {
        result.reason = "step5_start_budget_exhausted";
        return result;
    }
    state.lower_value = best;
    state.starts_consumed_in_epoch += static_cast<int>(
        result.starts_consumed
    );
    result.lower_value = best;
    result.terminal = best_point;
    const bool duplicate = std::any_of(
        state.M.begin(), state.M.end(), [&](const Held2MPoint& point) {
            return maximum_abs_difference(
                point.independent_modified_fractions,
                best_point.independent_modified_fractions
            ) <= kHeld2BasinDuplicateComposition.atol
                && std::abs(
                    std::log(point.volume) - std::log(best_point.volume)
                ) <= kHeld2BasinDuplicateLogVolume.atol;
        }
    );
    if (!duplicate) {
        state.M.push_back(best_point);
        state.start_epoch_added_member = true;
        result.mathematical_set_changed = true;
    }
    result.status = "complete";
    result.reason = "step5_complete";
    result.timing.terminal_status = result.status;
    result.timing.terminal_reason = result.reason;
    result.timing.next_step = 6;
    return result;
}

std::vector<std::vector<double>> held2_sobol_points(
    std::size_t dimension,
    int count
) {
    if (dimension == 0 || count < 0 || dimension > 10) {
        throw std::invalid_argument("HELD2 Sobol policy is invalid");
    }
    std::vector<std::array<std::uint32_t, 32>> directions;
    directions.reserve(dimension);
    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
        directions.push_back(sobol_directions(coordinate));
    }
    std::vector<std::vector<double>> points;
    points.reserve(static_cast<std::size_t>(count));
    constexpr double scale = 1.0 / 4294967296.0;
    for (int index = 1; index <= count; ++index) {
        const std::uint32_t gray = static_cast<std::uint32_t>(index)
            ^ (static_cast<std::uint32_t>(index) >> 1U);
        std::vector<double> point(dimension, 0.0);
        for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
            std::uint32_t value = 0;
            for (std::size_t bit = 0; bit < 32; ++bit) {
                if ((gray >> bit) & 1U) {
                    value ^= directions[coordinate][bit];
                }
            }
            point[coordinate] = static_cast<double>(value) * scale;
        }
        points.push_back(std::move(point));
    }
    return points;
}

Held2StageIIBasinExplorationResult explore_held2_stage_ii_basins(
    const Held2Coordinates& coordinates,
    const std::vector<Held2StageIIBasinSeed>& seeds,
    int sobol_count,
    bool use_direct_escalation,
    int direct_evaluation_budget,
    double total_ion_mole_fraction_max,
    const Held2StageIIBasinEvaluator& evaluator
) {
    if (coordinates.independent_indices.empty() || sobol_count < 0
        || direct_evaluation_budget < 0
        || (use_direct_escalation && direct_evaluation_budget == 0)
        || (!std::isnan(total_ion_mole_fraction_max)
            && (!std::isfinite(total_ion_mole_fraction_max)
                || total_ion_mole_fraction_max <= 0.0
                || total_ion_mole_fraction_max > 1.0))) {
        throw std::invalid_argument("HELD2 Stage-II basin policy is invalid");
    }
    Held2StageIIBasinExplorationResult result;
    result.direct_solver_version = nlopt_version_string();
    result.declared_sobol_count = sobol_count;
    result.declared_direct_budget = use_direct_escalation
        ? direct_evaluation_budget
        : 0;

    std::vector<Held2StageIIBasinSeed> unique_seeds;
    for (const Held2StageIIBasinSeed& seed : seeds) {
        const bool duplicate = std::any_of(
            unique_seeds.begin(),
            unique_seeds.end(),
            [&seed](const Held2StageIIBasinSeed& known) {
                return same_composition(
                    seed.independent_modified_fractions,
                    known.independent_modified_fractions
                );
            }
        );
        if (duplicate) {
            ++result.duplicate_start_count;
        } else {
            unique_seeds.push_back(seed);
        }
    }
    for (const Held2StageIIBasinSeed& seed : unique_seeds) {
        retain_evaluation(
            result,
            evaluate_fail_closed(
                evaluator,
                seed.independent_modified_fractions
            ),
            seed.source
        );
        if (result.failed_evaluation_count != 0) {
            result.termination_reason = "required_envelope_evaluation_failed";
            return result;
        }
    }
    for (const std::vector<double>& cube : held2_sobol_points(
             coordinates.independent_indices.size(),
             sobol_count
         )) {
        const std::vector<double> independent =
            held2_map_unit_cube_to_independent_fractions(
                coordinates,
                cube,
                total_ion_mole_fraction_max
            );
        retain_evaluation(
            result,
            evaluate_fail_closed(evaluator, independent),
            "sobol"
        );
        if (result.failed_evaluation_count != 0) {
            result.termination_reason = "required_envelope_evaluation_failed";
            return result;
        }
    }

    if (use_direct_escalation) {
        result.direct_escalation_used = true;
        nlopt::opt optimizer(
            nlopt::GN_DIRECT_L,
            coordinates.independent_indices.size()
        );
        DirectContext context{
            &coordinates,
            &evaluator,
            &result,
            &optimizer,
            total_ion_mole_fraction_max,
            false,
        };
        optimizer.set_lower_bounds(std::vector<double>(
            coordinates.independent_indices.size(),
            0.0
        ));
        optimizer.set_upper_bounds(std::vector<double>(
            coordinates.independent_indices.size(),
            1.0
        ));
        optimizer.set_maxeval(direct_evaluation_budget);
        optimizer.set_min_objective(direct_objective, &context);
        std::vector<double> initial(
            coordinates.independent_indices.size(),
            0.5
        );
        double minimum = std::numeric_limits<double>::infinity();
        try {
            const nlopt::result status = optimizer.optimize(initial, minimum);
            if (status != nlopt::MAXEVAL_REACHED) {
                result.termination_reason = "unexpected_direct_termination";
                return result;
            }
        } catch (const nlopt::forced_stop&) {
            if (result.failed_evaluation_count != 0) {
                return result;
            }
            result.termination_reason = "unexpected_direct_forced_stop";
            return result;
        } catch (const std::exception& error) {
            result.termination_reason = std::string("direct_solver_failure: ")
                + error.what();
            return result;
        }
    }

    if (result.representatives.empty()) {
        result.termination_reason = "no_physical_basin_representatives";
        return result;
    }
    std::sort(
        result.representatives.begin(),
        result.representatives.end(),
        [](const Held2StageIIPhysicalStart& left,
           const Held2StageIIPhysicalStart& right) {
            if (left.reduced_lower_value != right.reduced_lower_value) {
                return left.reduced_lower_value < right.reduced_lower_value;
            }
            if (left.independent_modified_fractions
                != right.independent_modified_fractions) {
                return left.independent_modified_fractions
                    < right.independent_modified_fractions;
            }
            return left.log_volume < right.log_volume;
        }
    );
    result.outcome = "representatives_found";
    result.termination_reason = use_direct_escalation
        ? "declared_direct_budget_exhausted"
        : "deterministic_exploration_completed";
    return result;
}

}  // namespace epcsaft_equilibrium
