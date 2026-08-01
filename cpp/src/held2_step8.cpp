#include "held2_step8.hpp"
#include "held2_tolerances.hpp"

#include <algorithm>
#include <cmath>

namespace epcsaft_equilibrium {
namespace {

std::vector<double> independent(
    const Held2Coordinates& coordinates,
    const std::vector<double>& modified
) {
    std::vector<double> result;
    result.reserve(coordinates.independent_indices.size());
    for (std::size_t provider : coordinates.independent_indices) {
        const auto position = std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            provider
        );
        result.push_back(modified[static_cast<std::size_t>(
            position - coordinates.retained_indices.begin()
        )]);
    }
    return result;
}

}  // namespace

Held2Step8Result run_held2_step8(
    const Held2Step1Result& step1,
    const Held2Step6Result& step6,
    const Held2StateEvaluator& evaluator,
    const Held2PackingFractionEvaluator& packing_fraction,
    const Held2Step8Result* previous,
    const Held2StateValueEvaluator& value_evaluator
) {
    Held2Step8Result result;
    result.timing.invocation_count = 1;
    if (!step1.coordinates || !step1.independent_feed || !step1.volume_bounds
        || step6.status != "complete" || step6.candidates.size() < 2
        || !evaluator || !packing_fraction) {
        result.reason = "invalid_step8_input";
        return result;
    }
    std::vector<Held2MPoint> selected_points;
    const bool continue_certified_active_set = previous
        && previous->outcome == Held2Step8Outcome::CertifiedFeasible
        && previous->nlp && previous->nlp->accepted
        && !previous->candidate_ids.empty();
    const bool continue_rejected_candidate_problem = previous
        && previous->outcome == Held2Step8Outcome::Indeterminate
        && !previous->candidate_ids.empty();
    const bool continue_candidate_problem =
        continue_certified_active_set
        || continue_rejected_candidate_problem;
    const std::uint64_t newest_id = std::max_element(
        step6.candidates.begin(),
        step6.candidates.end(),
        [](const Held2MPoint& left, const Held2MPoint& right) {
            return left.insertion_id < right.insertion_id;
        }
    )->insertion_id;
    const std::vector<std::uint64_t>& previously_attempted =
        previous && !previous->attempted_candidate_ids.empty()
        ? previous->attempted_candidate_ids
        : previous && !previous->problem_candidate_ids.empty()
            ? previous->problem_candidate_ids
            : previous ? previous->candidate_ids
                       : result.attempted_candidate_ids;
    const bool newest_is_unattempted = continue_candidate_problem
        && std::find(
            previously_attempted.begin(),
            previously_attempted.end(),
            newest_id
        ) == previously_attempted.end();
    for (const Held2MPoint& point : step6.candidates) {
        const bool retained = continue_candidate_problem
            && std::find(
                previous->candidate_ids.begin(),
                previous->candidate_ids.end(),
                point.insertion_id
            ) != previous->candidate_ids.end();
        if (!continue_candidate_problem || retained
            || (continue_certified_active_set
                && point.insertion_id == newest_id)
            || (continue_rejected_candidate_problem
                && newest_is_unattempted
                && point.insertion_id == newest_id)) {
            selected_points.push_back(point);
        }
    }
    if (selected_points.size() < 2) {
        selected_points = step6.candidates;
    }
    for (const Held2MPoint& point : selected_points) {
        result.problem_candidate_ids.push_back(point.insertion_id);
    }
    result.attempted_candidate_ids = previously_attempted;
    for (std::uint64_t id : result.problem_candidate_ids) {
        if (std::find(
                result.attempted_candidate_ids.begin(),
                result.attempted_candidate_ids.end(),
                id
            ) == result.attempted_candidate_ids.end()) {
            result.attempted_candidate_ids.push_back(id);
        }
    }
    std::vector<double> previous_effective_problem_variables = previous
        ? previous->problem_candidate_variables : std::vector<double>{};
    if (continue_certified_active_set) {
        for (Held2MPoint& point : selected_points) {
            const auto refined = std::find_if(
                previous->active_phases.begin(),
                previous->active_phases.end(),
                [&](const Held2Phase& phase) {
                    return phase.stable_id == point.insertion_id;
                }
            );
            if (refined != previous->active_phases.end()) {
                bool neighborhood_boundary_active = false;
                for (std::size_t coordinate = 0;
                     coordinate
                        < point.independent_modified_fractions.size();
                     ++coordinate) {
                    neighborhood_boundary_active =
                        neighborhood_boundary_active
                        || std::abs(
                            refined->independent_modified_fractions[
                                coordinate
                            ]
                            - point.independent_modified_fractions[
                                coordinate
                            ]
                        ) >= kHeld2Problem67Radius
                            - kHeld2BoundActivity.atol;
                }
                if (neighborhood_boundary_active) {
                    point.independent_modified_fractions =
                        refined->independent_modified_fractions;
                    point.volume = refined->volume;
                    const auto terminal_id = std::find(
                        previous->candidate_ids.begin(),
                        previous->candidate_ids.end(),
                        point.insertion_id
                    );
                    const auto problem_id = std::find(
                        previous->problem_candidate_ids.begin(),
                        previous->problem_candidate_ids.end(),
                        point.insertion_id
                    );
                    const std::size_t block_size =
                        point.independent_modified_fractions.size() + 1;
                    if (terminal_id != previous->candidate_ids.end()
                        && problem_id
                            != previous->problem_candidate_ids.end()
                        && previous->candidate_variables.size()
                            == previous->candidate_ids.size()
                                * block_size
                        && previous_effective_problem_variables.size()
                            == previous->problem_candidate_ids.size()
                                * block_size) {
                        const std::size_t source = static_cast<std::size_t>(
                            terminal_id - previous->candidate_ids.begin()
                        ) * block_size;
                        const std::size_t destination =
                            static_cast<std::size_t>(
                                problem_id
                                - previous->problem_candidate_ids.begin()
                            ) * block_size;
                        std::copy_n(
                            previous->candidate_variables.begin() + source,
                            block_size,
                            previous_effective_problem_variables.begin()
                                + destination
                        );
                    }
                }
            }
        }
    }
    for (const Held2MPoint& point : selected_points) {
        result.problem_candidate_variables.insert(
            result.problem_candidate_variables.end(),
            point.independent_modified_fractions.begin(),
            point.independent_modified_fractions.end()
        );
        result.problem_candidate_variables.push_back(point.volume);
    }
    const bool same_as_previous_problem = previous
        && result.problem_candidate_ids == previous->problem_candidate_ids
        && result.problem_candidate_variables
            == previous->problem_candidate_variables;
    const bool same_as_rejected_terminal = previous
        && previous->outcome == Held2Step8Outcome::Indeterminate
        && result.problem_candidate_ids == previous->candidate_ids
        && result.problem_candidate_variables
            == previous->candidate_variables;
    const bool same_as_certified_effective_problem = previous
        && previous->outcome == Held2Step8Outcome::CertifiedFeasible
        && result.problem_candidate_ids
            == previous->problem_candidate_ids
        && result.problem_candidate_variables
            == previous_effective_problem_variables;
    if (same_as_previous_problem || same_as_rejected_terminal
        || same_as_certified_effective_problem) {
        Held2Step8Result unchanged = *previous;
        unchanged.problem_candidate_ids = result.problem_candidate_ids;
        unchanged.problem_candidate_variables =
            result.problem_candidate_variables;
        unchanged.attempted_candidate_ids = result.attempted_candidate_ids;
        unchanged.warm_start_used = false;
        unchanged.cold_fallback_used = false;
        unchanged.provider_state_evaluations = 0;
        unchanged.provider_value_evaluations = 0;
        unchanged.provider_volume_bound_evaluations = 0;
        unchanged.provider_packing_evaluations = 0;
        unchanged.timing = {};
        unchanged.timing.invocation_count = 1;
        unchanged.timing.terminal_status =
            unchanged.outcome == Held2Step8Outcome::Indeterminate
            ? "indeterminate" : "complete";
        unchanged.timing.terminal_reason = "unchanged_problem_67";
        unchanged.timing.next_step =
            unchanged.outcome == Held2Step8Outcome::CertifiedFeasible
            ? 9 : 7;
        return unchanged;
    }
    for (const Held2MPoint& point : selected_points) {
        result.candidate_ids.push_back(point.insertion_id);
    }
    result.candidate_variables = result.problem_candidate_variables;

    std::vector<Held2StageIICandidate> candidates;
    std::vector<std::array<double, 2>> bounds;
    candidates.reserve(selected_points.size());
    bounds.reserve(selected_points.size());
    std::uint64_t provider_evaluations = 0;
    std::uint64_t provider_state_evaluations = 0;
    std::uint64_t provider_value_evaluations = 0;
    std::uint64_t provider_volume_bound_evaluations = 0;
    std::uint64_t provider_packing_evaluations = 0;
    const auto record_provider_work = [&] {
        result.timing.provider_evaluations = provider_evaluations;
        result.provider_state_evaluations =
            provider_state_evaluations;
        result.provider_value_evaluations =
            provider_value_evaluations;
        result.provider_volume_bound_evaluations =
            provider_volume_bound_evaluations;
        result.provider_packing_evaluations =
            provider_packing_evaluations;
    };
    for (const Held2MPoint& point : selected_points) {
        const std::array<double, 2> physical_bounds =
            (*step1.volume_bounds)(point.independent_modified_fractions);
        ++provider_evaluations;
        ++provider_volume_bound_evaluations;
        candidates.push_back({
            point.independent_modified_fractions,
            point.volume,
            std::log(point.volume),
        });
        bounds.push_back({
            std::log(physical_bounds[0]), std::log(physical_bounds[1]),
        });
    }
    const Held2StateEvaluator counted_evaluator =
        [&evaluator, &provider_evaluations,
         &provider_state_evaluations](
            const std::vector<double>& composition,
            double log_volume
        ) {
            ++provider_evaluations;
            ++provider_state_evaluations;
            return evaluator(composition, log_volume);
        };
    const Held2VolumeBoundsEvaluator counted_volume_bounds =
        [&step1, &provider_evaluations,
         &provider_volume_bound_evaluations](
            const std::vector<double>& composition
        ) {
            ++provider_evaluations;
            ++provider_volume_bound_evaluations;
            return (*step1.volume_bounds)(composition);
        };
    const Held2StateValueEvaluator counted_value =
        value_evaluator
        ? Held2StateValueEvaluator(
            [&value_evaluator, &provider_evaluations,
             &provider_value_evaluations](
                const std::vector<double>& composition,
                double log_volume
            ) {
                ++provider_evaluations;
                ++provider_value_evaluations;
                return value_evaluator(composition, log_volume);
            }
        )
        : Held2StateValueEvaluator{};
    std::vector<double> initial;
    const std::size_t dimension =
        step1.coordinates->independent_indices.size();
    const std::size_t block_size = dimension + 2;
    if (continue_certified_active_set
        && previous->continuation_variables.size()
            == previous->candidate_ids.size() * block_size) {
        initial.resize(candidates.size() * block_size);
        for (std::size_t phase = 0; phase < candidates.size(); ++phase) {
            const std::size_t offset = phase * block_size;
            const auto known = std::find(
                previous->candidate_ids.begin(),
                previous->candidate_ids.end(),
                selected_points[phase].insertion_id
            );
            const bool found = known != previous->candidate_ids.end();
            const std::size_t previous_offset = found
                ? static_cast<std::size_t>(
                    known - previous->candidate_ids.begin()
                ) * block_size
                : 0;
            initial[offset] = found
                ? previous->continuation_variables[previous_offset]
                : 0.0;
            for (std::size_t coordinate = 0;
                 coordinate < dimension;
                 ++coordinate) {
                initial[offset + 1 + coordinate] = std::clamp(
                    !found
                    ? candidates[phase]
                        .independent_modified_fractions[coordinate]
                    : previous->continuation_variables[
                        previous_offset + 1 + coordinate
                    ],
                    std::max(
                        step1.coordinates
                            ->independent_lower_bounds[coordinate],
                        candidates[phase]
                            .independent_modified_fractions[coordinate]
                            - kHeld2Problem67Radius
                    ),
                    std::min(
                        step1.coordinates
                            ->independent_upper_bounds[coordinate],
                        candidates[phase]
                            .independent_modified_fractions[coordinate]
                            + kHeld2Problem67Radius
                    )
                );
            }
            initial[offset + 1 + dimension] = std::clamp(
                !found
                ? candidates[phase].phase_coordinate
                : previous->continuation_variables[
                    previous_offset + 1 + dimension
                ],
                bounds[phase][0],
                bounds[phase][1]
            );
        }
    }
    const bool warm_started = !initial.empty();
    result.warm_start_used = warm_started;
    const std::vector<double> physical_feed =
        held2_lift_independent_fractions(
            *step1.coordinates, *step1.independent_feed
        );
    const auto solve = [&](std::vector<double> start) {
        return solve_held2_problem67(
            *step1.coordinates, physical_feed,
            candidates, counted_evaluator, bounds,
            std::move(start), counted_value, 32, true,
            counted_volume_bounds
        );
    };
    const auto accepted_nlp_evidence = [](const Held2Problem67Result& value) {
        return value.numerical_status == "converged"
            && audit_held2_tolerance(
                kHeld2Stage3ModifiedBalance,
                value.modified_balance_inf_norm
            ).passed
            && audit_held2_tolerance(
                kHeld2Stage3Stationarity,
                value.kkt_stationarity_inf_norm
            ).passed
            && audit_held2_tolerance(
                kHeld2Stage3DualSign,
                value.dual_sign_violation_inf_norm
            ).passed
            && audit_held2_tolerance(
                kHeld2Stage3Complementarity,
                value.bound_complementarity_inf_norm
            ).passed;
    };
    Held2Problem67Result solved = solve(std::move(initial));
    if (warm_started
        && solved.numerical_status != "not_adjudicated"
        && (!accepted_nlp_evidence(solved) || solved.phases.size() < 2)) {
        result.cold_fallback_used = true;
        Held2Problem67Result cold = solve({});
        cold.stage_iii_solve_count += solved.stage_iii_solve_count;
        cold.optimizer_iteration_count +=
            solved.optimizer_iteration_count;
        solved = std::move(cold);
    }
    record_provider_work();
    result.timing.optimizer_solves =
        static_cast<std::uint64_t>(solved.stage_iii_solve_count);
    result.timing.optimizer_iterations =
        static_cast<std::uint64_t>(solved.optimizer_iteration_count);
    if (!solved.candidate_indices.empty()) {
        result.candidate_ids.clear();
        result.candidate_variables.clear();
        for (std::size_t index : solved.candidate_indices) {
            result.candidate_ids.push_back(
                selected_points.at(index).insertion_id
            );
            result.candidate_variables.insert(
                result.candidate_variables.end(),
                selected_points.at(index)
                    .independent_modified_fractions.begin(),
                selected_points.at(index)
                    .independent_modified_fractions.end()
            );
            result.candidate_variables.push_back(
                selected_points.at(index).volume
            );
        }
    }
    result.continuation_variables = solved.solution_variables;
    result.ordinary_balance_inf = solved.ordinary_balance_inf_norm;
    result.electroneutrality_inf = solved.phase_charge_inf_norm;
    result.electroneutrality_scale = solved.phase_charge_scale;
    result.pressure_residual_inf = solved.pressure_stationarity_inf_norm;
    result.phase_coalescences = solved.phase_coalescences;
    result.farkas = solved.feasibility_certificate;
    const bool nlp_attempted =
        solved.numerical_status != "not_adjudicated";
    const bool nlp_accepted = nlp_attempted
        && accepted_nlp_evidence(solved);
    if (nlp_attempted) {
        result.nlp = Held2NlpCertificate{
            solved.solver_status,
            solved.modified_balance_inf_norm,
            solved.kkt_stationarity_inf_norm,
            solved.dual_sign_violation_inf_norm,
            solved.bound_complementarity_inf_norm,
            nlp_accepted,
        };
    }

    if (solved.solver_status == "infeasible_problem_detected"
        && result.farkas && result.farkas->accepted) {
        result.outcome = Held2Step8Outcome::CertifiedInfeasible;
        result.reason = "problem_67_infeasible";
        result.timing.terminal_status = "complete";
        result.timing.terminal_reason = result.reason;
        result.timing.next_step = 7;
        return result;
    }
    if (solved.failure_reason == "collapsed_phase_set") {
        result.outcome = Held2Step8Outcome::InsufficientCandidates;
        result.reason = solved.failure_reason;
        result.timing.terminal_status = "complete";
        result.timing.terminal_reason = result.reason;
        result.timing.next_step = 7;
        return result;
    }
    if (solved.numerical_status != "converged"
        || !nlp_accepted
        || solved.phases.size() < 2) {
        result.reason = solved.failure_reason.empty()
            ? "problem_67_not_converged" : solved.failure_reason;
        result.timing.terminal_status = "indeterminate";
        result.timing.terminal_reason = result.reason;
        return result;
    }

    for (std::size_t phase_index = 0;
         phase_index < solved.phases.size(); ++phase_index) {
        const Held2Problem67Phase& phase = solved.phases[phase_index];
        const std::vector<double> composition =
            independent(*step1.coordinates, phase.modified_fractions);
        ++provider_evaluations;
        ++provider_state_evaluations;
        const Held2StateEvaluation state =
            evaluator(composition, std::log(phase.volume));
        ++provider_evaluations;
        ++provider_packing_evaluations;
        const double phase_packing_fraction =
            packing_fraction(composition, phase.volume);
        result.active_phases.push_back({
            result.candidate_ids.at(phase_index),
            phase.phase_fraction,
            composition,
            phase.physical_fractions,
            phase.volume,
            phase_packing_fraction,
            state.helmholtz_over_rt_reference_amount,
            state.pressure_pa,
            state.chemical_potentials_over_rt,
            state.objective,
            state.gradient,
        });
    }
    record_provider_work();
    std::sort(
        result.active_phases.begin(),
        result.active_phases.end(),
        [](const Held2Phase& left, const Held2Phase& right) {
            return left.stable_id < right.stable_id;
        }
    );
    result.candidate_variables.clear();
    for (const Held2Phase& phase : result.active_phases) {
        result.candidate_variables.insert(
            result.candidate_variables.end(),
            phase.independent_modified_fractions.begin(),
            phase.independent_modified_fractions.end()
        );
        result.candidate_variables.push_back(phase.volume);
    }
    result.outcome = Held2Step8Outcome::CertifiedFeasible;
    result.reason = "step8_complete";
    result.total_reduced_gibbs = solved.objective;
    result.timing.terminal_status = "complete";
    result.timing.terminal_reason = result.reason;
    result.timing.next_step = 9;
    return result;
}

}  // namespace epcsaft_equilibrium
