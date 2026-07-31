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
    const bool continue_active_set = previous
        && previous->outcome == Held2Step8Outcome::CertifiedFeasible
        && previous->nlp && previous->nlp->accepted
        && !previous->candidate_ids.empty();
    const std::uint64_t newest_id = std::max_element(
        step6.candidates.begin(),
        step6.candidates.end(),
        [](const Held2MPoint& left, const Held2MPoint& right) {
            return left.insertion_id < right.insertion_id;
        }
    )->insertion_id;
    for (const Held2MPoint& point : step6.candidates) {
        const bool active = continue_active_set
            && std::any_of(
                previous->active_phases.begin(),
                previous->active_phases.end(),
                [&point](const Held2Phase& phase) {
                    return phase.stable_id == point.insertion_id;
                }
            );
        if (!continue_active_set || active
            || point.insertion_id == newest_id) {
            selected_points.push_back(point);
        }
    }
    if (selected_points.size() < 2) {
        selected_points = step6.candidates;
    }
    if (continue_active_set) {
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
                }
            }
        }
    }
    for (const Held2MPoint& point : selected_points) {
        result.candidate_ids.push_back(point.insertion_id);
    }

    std::vector<Held2StageIICandidate> candidates;
    std::vector<std::array<double, 2>> bounds;
    candidates.reserve(selected_points.size());
    bounds.reserve(selected_points.size());
    std::uint64_t provider_evaluations = 0;
    for (const Held2MPoint& point : selected_points) {
        const std::array<double, 2> physical_bounds =
            (*step1.volume_bounds)(point.independent_modified_fractions);
        ++provider_evaluations;
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
        [&evaluator, &provider_evaluations](
            const std::vector<double>& composition,
            double log_volume
        ) {
            ++provider_evaluations;
            return evaluator(composition, log_volume);
        };
    const Held2StateValueEvaluator counted_value =
        value_evaluator
        ? Held2StateValueEvaluator(
            [&value_evaluator, &provider_evaluations](
                const std::vector<double>& composition,
                double log_volume
            ) {
                ++provider_evaluations;
                return value_evaluator(composition, log_volume);
            }
        )
        : Held2StateValueEvaluator{};
    std::vector<double> initial;
    const std::size_t dimension =
        step1.coordinates->independent_indices.size();
    const std::size_t block_size = dimension + 2;
    if (previous
        && previous->outcome == Held2Step8Outcome::CertifiedFeasible
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
    const std::vector<double> physical_feed =
        held2_lift_independent_fractions(
            *step1.coordinates, *step1.independent_feed
        );
    const auto solve = [&](std::vector<double> start) {
        return solve_held2_problem67(
            *step1.coordinates, physical_feed,
            candidates, counted_evaluator, bounds,
            std::move(start), counted_value
        );
    };
    Held2Problem67Result solved = solve(std::move(initial));
    if (warm_started
        && solved.failure_reason == "stage_iii_solver_not_converged") {
        Held2Problem67Result cold = solve({});
        cold.stage_iii_solve_count += solved.stage_iii_solve_count;
        cold.optimizer_iteration_count +=
            solved.optimizer_iteration_count;
        solved = std::move(cold);
    }
    result.timing.provider_evaluations = provider_evaluations;
    result.timing.optimizer_solves =
        static_cast<std::uint64_t>(solved.stage_iii_solve_count);
    result.timing.optimizer_iterations =
        static_cast<std::uint64_t>(solved.optimizer_iteration_count);
    if (!solved.candidate_indices.empty()) {
        result.candidate_ids.clear();
        for (std::size_t index : solved.candidate_indices) {
            result.candidate_ids.push_back(
                selected_points.at(index).insertion_id
            );
        }
    }
    result.continuation_variables = solved.solution_variables;
    result.ordinary_balance_inf = solved.ordinary_balance_inf_norm;
    result.electroneutrality_inf = solved.phase_charge_inf_norm;
    result.electroneutrality_scale = solved.phase_charge_scale;
    result.pressure_residual_inf = solved.pressure_stationarity_inf_norm;
    result.farkas = solved.feasibility_certificate;
    const bool nlp_attempted =
        solved.numerical_status != "not_adjudicated";
    const bool nlp_accepted = nlp_attempted
        && solved.numerical_status == "converged"
        && audit_held2_tolerance(
        kHeld2Stage3ModifiedBalance, solved.modified_balance_inf_norm
    ).passed && audit_held2_tolerance(
        kHeld2Stage3Stationarity, solved.kkt_stationarity_inf_norm
    ).passed && audit_held2_tolerance(
        kHeld2Stage3DualSign, solved.dual_sign_violation_inf_norm
    ).passed && audit_held2_tolerance(
        kHeld2Stage3Complementarity,
        solved.bound_complementarity_inf_norm
    ).passed;
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
        const Held2StateEvaluation state =
            evaluator(composition, std::log(phase.volume));
        ++provider_evaluations;
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
    result.timing.provider_evaluations = provider_evaluations;
    std::sort(
        result.active_phases.begin(),
        result.active_phases.end(),
        [](const Held2Phase& left, const Held2Phase& right) {
            return left.stable_id < right.stable_id;
        }
    );
    result.outcome = Held2Step8Outcome::CertifiedFeasible;
    result.reason = "step8_complete";
    result.total_reduced_gibbs = solved.objective;
    result.timing.terminal_status = "complete";
    result.timing.terminal_reason = result.reason;
    result.timing.next_step = 9;
    return result;
}

}  // namespace epcsaft_equilibrium
