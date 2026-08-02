#include "held2_step8.hpp"
#include "held2_tolerances.hpp"

#include <algorithm>
#include <cmath>
#include <functional>

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

bool all_finite(const std::vector<double>& values) {
    return std::all_of(
        values.begin(), values.end(),
        [](double value) { return std::isfinite(value); }
    );
}

bool all_finite(const std::vector<std::array<double, 2>>& values) {
    return std::all_of(
        values.begin(), values.end(),
        [](const auto& value) {
            return std::isfinite(value[0]) && std::isfinite(value[1]);
        }
    );
}

bool all_finite(const Held2Coordinates& coordinates) {
    if (!all_finite(coordinates.charges)
        || !all_finite(coordinates.modified_factors)
        || !all_finite(coordinates.independent_lower_bounds)
        || !all_finite(coordinates.independent_upper_bounds)) {
        return false;
    }
    return std::all_of(
        coordinates.polytope_constraints.begin(),
        coordinates.polytope_constraints.end(),
        [](const Held2PolytopeConstraint& constraint) {
            return all_finite(constraint.coefficients)
                && std::isfinite(constraint.upper_bound);
        }
    );
}

bool same_coordinates(
    const Held2Coordinates& left,
    const Held2Coordinates& right
) {
    if (left.charges != right.charges
        || left.eliminated_index != right.eliminated_index
        || left.dependent_index != right.dependent_index
        || left.paper_to_provider_indices
            != right.paper_to_provider_indices
        || left.provider_to_paper_indices
            != right.provider_to_paper_indices
        || left.compact_to_paper_indices
            != right.compact_to_paper_indices
        || left.retained_indices != right.retained_indices
        || left.independent_indices != right.independent_indices
        || left.modified_factors != right.modified_factors
        || left.independent_lower_bounds
            != right.independent_lower_bounds
        || left.independent_upper_bounds
            != right.independent_upper_bounds
        || left.polytope_constraints.size()
            != right.polytope_constraints.size()) {
        return false;
    }
    for (std::size_t index = 0;
         index < left.polytope_constraints.size(); ++index) {
        const Held2PolytopeConstraint& lhs =
            left.polytope_constraints[index];
        const Held2PolytopeConstraint& rhs =
            right.polytope_constraints[index];
        if (lhs.coefficients != rhs.coefficients
            || lhs.upper_bound != rhs.upper_bound) {
            return false;
        }
    }
    return true;
}

template <typename Value>
void mix_hash(std::size_t& seed, const Value& value) {
    const std::size_t hashed = std::hash<Value>{}(value);
    seed ^= hashed + 0x9e3779b9U + (seed << 6U) + (seed >> 2U);
}

template <typename Value>
void mix_values(std::size_t& seed, const std::vector<Value>& values) {
    mix_hash(seed, values.size());
    for (const Value& value : values) {
        mix_hash(seed, value);
    }
}

void mix_coordinates(std::size_t& seed, const Held2Coordinates& coordinates) {
    mix_values(seed, coordinates.charges);
    mix_hash(seed, coordinates.eliminated_index);
    mix_hash(seed, coordinates.dependent_index);
    mix_values(seed, coordinates.paper_to_provider_indices);
    mix_values(seed, coordinates.provider_to_paper_indices);
    mix_values(seed, coordinates.compact_to_paper_indices);
    mix_values(seed, coordinates.retained_indices);
    mix_values(seed, coordinates.independent_indices);
    mix_values(seed, coordinates.modified_factors);
    mix_values(seed, coordinates.independent_lower_bounds);
    mix_values(seed, coordinates.independent_upper_bounds);
    mix_hash(seed, coordinates.polytope_constraints.size());
    for (const Held2PolytopeConstraint& constraint
         : coordinates.polytope_constraints) {
        mix_values(seed, constraint.coefficients);
        mix_hash(seed, constraint.upper_bound);
    }
}

}  // namespace

bool Held2AlgorithmCache::ProblemKey::operator==(
    const ProblemKey& other
) const {
    return candidate_ids == other.candidate_ids
        && candidate_variables == other.candidate_variables
        && same_coordinates(coordinates, other.coordinates)
        && physical_feed == other.physical_feed
        && phase_coordinate_bounds == other.phase_coordinate_bounds
        && neighborhood_radius == other.neighborhood_radius;
}

std::size_t Held2AlgorithmCache::ProblemHash::operator()(
    const ProblemKey& problem
) const noexcept {
    std::size_t seed = problem.candidate_ids.size();
    for (std::uint64_t id : problem.candidate_ids) {
        mix_hash(seed, id);
    }
    mix_values(seed, problem.candidate_variables);
    mix_coordinates(seed, problem.coordinates);
    mix_values(seed, problem.physical_feed);
    mix_hash(seed, problem.phase_coordinate_bounds.size());
    for (const auto& bounds : problem.phase_coordinate_bounds) {
        mix_hash(seed, bounds[0]);
        mix_hash(seed, bounds[1]);
    }
    mix_hash(seed, problem.neighborhood_radius);
    return seed;
}

std::size_t Held2AlgorithmCache::StateHash::operator()(
    const StateKey& state
) const noexcept {
    std::size_t seed = 0;
    mix_values(seed, state.variables);
    return seed;
}

std::size_t Held2AlgorithmCache::VolumeBoundsHash::operator()(
    const VolumeBoundsKey& bounds
) const noexcept {
    std::size_t seed = 0;
    mix_values(seed, bounds.composition);
    return seed;
}

void Held2AlgorithmCache::begin_context(const void* context_token) {
    if (context_token && context_token == context_token_) {
        return;
    }
    states_.clear();
    volume_bounds_.clear();
    problems_.clear();
    provider_state_evaluations_ = 0;
    state_cache_hits_ = 0;
    context_token_ = context_token;
}

Held2StateEvaluation Held2AlgorithmCache::evaluate_state(
    const Held2StateEvaluator& evaluator,
    const std::vector<double>& composition,
    double log_volume,
    bool* provider_evaluated
) {
    if (provider_evaluated) {
        *provider_evaluated = false;
    }
    const bool cacheable = all_finite(composition)
        && std::isfinite(log_volume);
    StateKey key{composition};
    key.variables.push_back(log_volume);
    if (cacheable) {
        const auto retained = states_.find(key);
        if (retained != states_.end()) {
            ++state_cache_hits_;
            return retained->second;
        }
    }
    Held2StateEvaluation state = evaluator(composition, log_volume);
    ++provider_state_evaluations_;
    if (provider_evaluated) {
        *provider_evaluated = true;
    }
    if (cacheable) {
        states_.insert_or_assign(std::move(key), state);
    }
    return state;
}

const std::array<double, 2>*
Held2AlgorithmCache::find_volume_bounds(
    const std::vector<double>& composition
) const {
    if (!all_finite(composition)) {
        return nullptr;
    }
    const auto retained = volume_bounds_.find({composition});
    return retained == volume_bounds_.end() ? nullptr : &retained->second;
}

void Held2AlgorithmCache::retain_volume_bounds(
    const std::vector<double>& composition,
    std::array<double, 2> bounds
) {
    if (!all_finite(composition)) {
        return;
    }
    volume_bounds_.insert_or_assign({composition}, bounds);
}

const Held2Step8Result* Held2AlgorithmCache::find_problem(
    const std::vector<std::uint64_t>& candidate_ids,
    const std::vector<double>& candidate_variables,
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const std::vector<std::array<double, 2>>& phase_coordinate_bounds,
    double neighborhood_radius
) const {
    if (!all_finite(candidate_variables)
        || !all_finite(coordinates)
        || !all_finite(physical_feed)
        || !all_finite(phase_coordinate_bounds)
        || !std::isfinite(neighborhood_radius)) {
        return nullptr;
    }
    const auto retained = problems_.find({
        candidate_ids,
        candidate_variables,
        coordinates,
        physical_feed,
        phase_coordinate_bounds,
        neighborhood_radius,
    });
    return retained == problems_.end() ? nullptr : &retained->second;
}

void Held2AlgorithmCache::retain_problem(
    const std::vector<std::uint64_t>& candidate_ids,
    const std::vector<double>& candidate_variables,
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const std::vector<std::array<double, 2>>& phase_coordinate_bounds,
    double neighborhood_radius,
    Held2Step8Result result
) {
    if ((result.outcome != Held2Step8Outcome::CertifiedFeasible
         && result.outcome != Held2Step8Outcome::CertifiedInfeasible)
        || !all_finite(candidate_variables)
        || !all_finite(coordinates)
        || !all_finite(physical_feed)
        || !all_finite(phase_coordinate_bounds)
        || !std::isfinite(neighborhood_radius)) {
        return;
    }
    problems_.insert_or_assign(
        {
            candidate_ids,
            candidate_variables,
            coordinates,
            physical_feed,
            phase_coordinate_bounds,
            neighborhood_radius,
        },
        std::move(result)
    );
}

Held2Step8Result run_held2_step8(
    const Held2Step1Result& step1,
    const Held2Step6Result& step6,
    const Held2StateEvaluator& evaluator,
    const Held2PackingFractionEvaluator& packing_fraction,
    const Held2Step8Result* previous,
    double neighborhood_radius,
    Held2AlgorithmCache* shared_cache,
    const void* cache_context
) {
    Held2Step8Result result;
    const int local_context = 0;
    Held2AlgorithmCache local_cache;
    Held2AlgorithmCache& cache = shared_cache ? *shared_cache : local_cache;
    cache.begin_context(
        shared_cache ? cache_context : static_cast<const void*>(&local_context)
    );
    result.neighborhood_radius = neighborhood_radius;
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
                        ) >= neighborhood_radius
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
        result.problem_candidate_variables.insert(
            result.problem_candidate_variables.end(),
            point.independent_modified_fractions.begin(),
            point.independent_modified_fractions.end()
        );
        result.problem_candidate_variables.push_back(point.volume);
    }
    const auto reuse_problem = [&](
        const Held2Step8Result& retained,
        const char* reason
    ) {
        Held2Step8Result unchanged = retained;
        unchanged.problem_candidate_ids = result.problem_candidate_ids;
        unchanged.problem_candidate_variables =
            result.problem_candidate_variables;
        unchanged.attempted_candidate_ids = result.attempted_candidate_ids;
        unchanged.warm_start_used = false;
        unchanged.cold_fallback_used = false;
        unchanged.provider_state_evaluations = 0;
        unchanged.provider_volume_bound_evaluations = 0;
        unchanged.provider_packing_evaluations = 0;
        unchanged.timing = {};
        unchanged.timing.invocation_count = 1;
        unchanged.timing.terminal_status =
            unchanged.outcome == Held2Step8Outcome::Indeterminate
            ? "indeterminate" : "complete";
        unchanged.timing.terminal_reason = reason;
        unchanged.timing.next_step =
            unchanged.outcome == Held2Step8Outcome::CertifiedFeasible
            ? 9 : 7;
        return unchanged;
    };
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
    std::uint64_t provider_volume_bound_evaluations = 0;
    std::uint64_t provider_packing_evaluations = 0;
    const auto record_provider_work = [&] {
        result.timing.provider_evaluations = provider_evaluations;
        result.provider_state_evaluations =
            provider_state_evaluations;
        result.provider_volume_bound_evaluations =
            provider_volume_bound_evaluations;
        result.provider_packing_evaluations =
            provider_packing_evaluations;
    };
    const auto cached_volume_bounds = [&](
        const std::vector<double>& composition
    ) {
        if (const auto* retained = cache.find_volume_bounds(composition)) {
            return *retained;
        }
        ++provider_evaluations;
        ++provider_volume_bound_evaluations;
        const std::array<double, 2> bounds =
            (*step1.volume_bounds)(composition);
        cache.retain_volume_bounds(composition, bounds);
        return bounds;
    };
    for (const Held2MPoint& point : selected_points) {
        const std::array<double, 2> physical_bounds =
            cached_volume_bounds(point.independent_modified_fractions);
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
         &provider_state_evaluations, &cache](
            const std::vector<double>& composition,
            double log_volume
        ) {
            bool provider_evaluated = false;
            Held2StateEvaluation state = cache.evaluate_state(
                evaluator, composition, log_volume, &provider_evaluated
            );
            if (provider_evaluated) {
                ++provider_evaluations;
                ++provider_state_evaluations;
            }
            return state;
        };
    const Held2VolumeBoundsEvaluator counted_volume_bounds =
        [&cached_volume_bounds](
            const std::vector<double>& composition
        ) {
            return cached_volume_bounds(composition);
        };
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
                            - neighborhood_radius
                    ),
                    std::min(
                        step1.coordinates
                            ->independent_upper_bounds[coordinate],
                        candidates[phase]
                            .independent_modified_fractions[coordinate]
                            + neighborhood_radius
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
    if (const auto* retained = cache.find_problem(
            result.problem_candidate_ids,
            result.problem_candidate_variables,
            *step1.coordinates,
            physical_feed,
            bounds,
            neighborhood_radius
        )) {
        return reuse_problem(*retained, "cached_problem_67");
    }
    const auto retain_problem = [&] {
        cache.retain_problem(
            result.problem_candidate_ids,
            result.problem_candidate_variables,
            *step1.coordinates,
            physical_feed,
            bounds,
            result.neighborhood_radius,
            result
        );
    };
    const auto solve = [&](
        std::vector<double> start,
        double radius
    ) {
        return solve_held2_problem67(
            *step1.coordinates, physical_feed,
            candidates, counted_evaluator, bounds,
            std::move(start), 32, true,
            counted_volume_bounds, radius
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
    Held2Problem67Result solved = solve(
        std::move(initial), neighborhood_radius
    );
    if (solved.failure_reason == "problem_67_infeasible"
        && neighborhood_radius == kHeld2Problem67InitialRadius) {
        Held2Problem67Result expanded = solve(
            {}, kHeld2Problem67ExpandedRadius
        );
        if (accepted_nlp_evidence(expanded)
            && expanded.phases.size() >= 2) {
            expanded.stage_iii_solve_count +=
                solved.stage_iii_solve_count;
            expanded.optimizer_iteration_count +=
                solved.optimizer_iteration_count;
            solved = std::move(expanded);
            result.neighborhood_radius =
                kHeld2Problem67ExpandedRadius;
        } else {
            solved.stage_iii_solve_count +=
                expanded.stage_iii_solve_count;
            solved.optimizer_iteration_count +=
                expanded.optimizer_iteration_count;
        }
    }
    if (warm_started
        && solved.numerical_status != "not_adjudicated"
        && (!accepted_nlp_evidence(solved) || solved.phases.size() < 2)) {
        result.cold_fallback_used = true;
        Held2Problem67Result cold = solve(
            {}, result.neighborhood_radius
        );
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
        retain_problem();
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
        const Held2StateEvaluation state =
            counted_evaluator(composition, std::log(phase.volume));
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
    retain_problem();
    return result;
}

}  // namespace epcsaft_equilibrium
