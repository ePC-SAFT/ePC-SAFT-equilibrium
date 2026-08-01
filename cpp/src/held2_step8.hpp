#pragma once

#include "held2_step6.hpp"

#include <unordered_map>

namespace epcsaft_equilibrium {

enum class Held2Step8Outcome {
    CertifiedFeasible,
    CertifiedInfeasible,
    InsufficientCandidates,
    Indeterminate,
};

struct Held2Phase {
    std::uint64_t stable_id = 0;
    double phase_fraction = 0.0;
    std::vector<double> independent_modified_fractions;
    std::vector<double> physical_fractions_provider_order;
    double volume = 0.0;
    double packing_fraction = 0.0;
    double helmholtz_over_rt_reference_amount = 0.0;
    double pressure_pa = 0.0;
    std::vector<double> chemical_potentials_over_rt;
    double reduced_gibbs = 0.0;
    std::vector<double> reduced_gibbs_gradient;
};

struct Held2NlpCertificate {
    std::string solver_status = "not_run";
    double primal_residual_inf = std::numeric_limits<double>::infinity();
    double stationarity_residual_inf =
        std::numeric_limits<double>::infinity();
    double dual_sign_violation_inf =
        std::numeric_limits<double>::infinity();
    double complementarity_inf = std::numeric_limits<double>::infinity();
    bool accepted = false;
};

struct Held2Step8Result {
    Held2Step8Outcome outcome = Held2Step8Outcome::Indeterminate;
    std::string reason = "not_run";
    std::vector<std::uint64_t> problem_candidate_ids;
    std::vector<double> problem_candidate_variables;
    std::vector<std::uint64_t> attempted_candidate_ids;
    std::vector<std::uint64_t> candidate_ids;
    std::vector<double> candidate_variables;
    std::vector<double> continuation_variables;
    double neighborhood_radius = kHeld2Problem67InitialRadius;
    bool warm_start_used = false;
    bool cold_fallback_used = false;
    std::uint64_t provider_state_evaluations = 0;
    std::uint64_t provider_volume_bound_evaluations = 0;
    std::uint64_t provider_packing_evaluations = 0;
    std::optional<double> total_reduced_gibbs;
    double ordinary_balance_inf = std::numeric_limits<double>::infinity();
    double electroneutrality_inf = std::numeric_limits<double>::infinity();
    double electroneutrality_scale = 0.0;
    double pressure_residual_inf = std::numeric_limits<double>::infinity();
    std::vector<Held2Phase> active_phases;
    std::vector<Held2PhaseCoalescence> phase_coalescences;
    std::optional<Held2FarkasCertificate> farkas;
    std::optional<Held2NlpCertificate> nlp;
    Held2StepTiming timing;
};

class Held2AlgorithmCache {
public:
    void begin_context(const void* context_token);
    [[nodiscard]] Held2StateEvaluation evaluate_state(
        const Held2StateEvaluator& evaluator,
        const std::vector<double>& composition,
        double log_volume,
        bool* provider_evaluated = nullptr
    );
    [[nodiscard]] const std::array<double, 2>* find_volume_bounds(
        const std::vector<double>& composition
    ) const;
    void retain_volume_bounds(
        const std::vector<double>& composition,
        std::array<double, 2> bounds
    );
    [[nodiscard]] std::uint64_t provider_state_evaluations() const {
        return provider_state_evaluations_;
    }
    [[nodiscard]] std::uint64_t state_cache_hits() const {
        return state_cache_hits_;
    }
    [[nodiscard]] const Held2Step8Result* find_problem(
        const std::vector<std::uint64_t>& candidate_ids,
        const std::vector<double>& candidate_variables,
        const Held2Coordinates& coordinates,
        const std::vector<double>& physical_feed,
        const std::vector<std::array<double, 2>>& phase_coordinate_bounds,
        double neighborhood_radius
    ) const;
    void retain_problem(
        const std::vector<std::uint64_t>& candidate_ids,
        const std::vector<double>& candidate_variables,
        const Held2Coordinates& coordinates,
        const std::vector<double>& physical_feed,
        const std::vector<std::array<double, 2>>& phase_coordinate_bounds,
        double neighborhood_radius,
        Held2Step8Result result
    );

private:
    struct ProblemKey {
        std::vector<std::uint64_t> candidate_ids;
        std::vector<double> candidate_variables;
        Held2Coordinates coordinates;
        std::vector<double> physical_feed;
        std::vector<std::array<double, 2>> phase_coordinate_bounds;
        double neighborhood_radius = 0.0;

        [[nodiscard]] bool operator==(const ProblemKey& other) const;
    };
    struct ProblemHash {
        [[nodiscard]] std::size_t operator()(
            const ProblemKey& problem
        ) const noexcept;
    };

    struct StateKey {
        std::vector<double> variables;

        [[nodiscard]] bool operator==(const StateKey& other) const {
            return variables == other.variables;
        }
    };
    struct StateHash {
        [[nodiscard]] std::size_t operator()(
            const StateKey& state
        ) const noexcept;
    };
    struct VolumeBoundsKey {
        std::vector<double> composition;

        [[nodiscard]] bool operator==(
            const VolumeBoundsKey& other
        ) const {
            return composition == other.composition;
        }
    };
    struct VolumeBoundsHash {
        [[nodiscard]] std::size_t operator()(
            const VolumeBoundsKey& bounds
        ) const noexcept;
    };

    std::unordered_map<
        StateKey, Held2StateEvaluation, StateHash
    > states_;
    std::unordered_map<
        VolumeBoundsKey, std::array<double, 2>, VolumeBoundsHash
    > volume_bounds_;
    std::unordered_map<ProblemKey, Held2Step8Result, ProblemHash> problems_;
    std::uint64_t provider_state_evaluations_ = 0;
    std::uint64_t state_cache_hits_ = 0;
    const void* context_token_ = nullptr;
};

[[nodiscard]] Held2Step8Result run_held2_step8(
    const Held2Step1Result& step1,
    const Held2Step6Result& step6,
    const Held2StateEvaluator& evaluator,
    const Held2PackingFractionEvaluator& packing_fraction,
    const Held2Step8Result* previous = nullptr,
    double neighborhood_radius = kHeld2Problem67InitialRadius,
    Held2AlgorithmCache* cache = nullptr,
    const void* cache_context = nullptr
);

}  // namespace epcsaft_equilibrium
