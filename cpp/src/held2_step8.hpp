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
    Held2AlgorithmCache(
        std::string provider_fingerprint,
        Held2StateEvaluator state_evaluator,
        Held2VolumeBoundsEvaluator volume_bounds_evaluator,
        Held2PackingFractionEvaluator packing_fraction_evaluator
    );
    [[nodiscard]] const std::string& provider_fingerprint() const noexcept {
        return provider_fingerprint_;
    }
    [[nodiscard]] Held2StateEvaluation evaluate_state(
        const std::vector<double>& composition,
        double log_volume
    );
    [[nodiscard]] std::array<double, 2> evaluate_volume_bounds(
        const std::vector<double>& composition
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
    [[nodiscard]] std::uint64_t provider_volume_bound_evaluations() const {
        return provider_volume_bound_evaluations_;
    }
    [[nodiscard]] double evaluate_packing_fraction(
        const std::vector<double>& composition,
        double volume
    );
    [[nodiscard]] std::uint64_t provider_packing_evaluations() const {
        return provider_packing_evaluations_;
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
    struct PackingKey {
        std::vector<double> variables;

        [[nodiscard]] bool operator==(const PackingKey& other) const {
            return variables == other.variables;
        }
    };
    struct PackingHash {
        [[nodiscard]] std::size_t operator()(
            const PackingKey& packing
        ) const noexcept;
    };

    std::unordered_map<
        StateKey, Held2StateEvaluation, StateHash
    > states_;
    std::unordered_map<
        VolumeBoundsKey, std::array<double, 2>, VolumeBoundsHash
    > volume_bounds_;
    std::unordered_map<PackingKey, double, PackingHash> packing_fractions_;
    std::unordered_map<ProblemKey, Held2Step8Result, ProblemHash> problems_;
    std::uint64_t provider_state_evaluations_ = 0;
    std::uint64_t provider_volume_bound_evaluations_ = 0;
    std::uint64_t provider_packing_evaluations_ = 0;
    std::string provider_fingerprint_;
    Held2StateEvaluator state_evaluator_;
    Held2VolumeBoundsEvaluator volume_bounds_evaluator_;
    Held2PackingFractionEvaluator packing_fraction_evaluator_;
};

[[nodiscard]] Held2Step8Result run_held2_step8(
    const Held2Step1Result& step1,
    const Held2Step6Result& step6,
    Held2AlgorithmCache& cache,
    const Held2Step8Result* previous = nullptr,
    double neighborhood_radius = kHeld2Problem67InitialRadius
);

}  // namespace epcsaft_equilibrium
