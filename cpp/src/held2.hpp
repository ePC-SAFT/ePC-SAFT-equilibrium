#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace epcsaft_equilibrium {

inline constexpr const char* kHeld2ManufacturedFormulationId =
    "perdomo-held2.modified-mole.manufactured.v1";

struct Held2Coordinates {
    std::vector<double> charges;
    std::size_t eliminated_index = 0;
    std::size_t dependent_index = 0;
    std::vector<std::size_t> retained_indices;
    std::vector<std::size_t> independent_indices;
    std::vector<double> modified_factors;
    std::vector<double> independent_lower_bounds;
    std::vector<double> independent_upper_bounds;
};

struct Held2StateEvaluation;

[[nodiscard]] Held2Coordinates make_held2_coordinates(
    const std::vector<double>& charges
);

[[nodiscard]] std::vector<double> held2_transform_physical_fractions(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_fractions
);

[[nodiscard]] std::vector<double> held2_lift_modified_fractions(
    const Held2Coordinates& coordinates,
    const std::vector<double>& modified_fractions
);

[[nodiscard]] std::vector<double> held2_lift_independent_fractions(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent_modified_fractions
);

[[nodiscard]] std::vector<double> held2_transform_modified_potentials(
    const Held2Coordinates& coordinates,
    const std::vector<double>& chemical_potentials
);

[[nodiscard]] Held2StateEvaluation evaluate_held2_manufactured_state(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent_modified_fractions,
    double log_volume
);

[[nodiscard]] Held2StateEvaluation evaluate_held2_manufactured_search_objective(
    const Held2Coordinates& coordinates,
    const std::vector<double>& variables,
    const std::vector<double>& reference_variables,
    bool use_tpd
);

[[nodiscard]] double held2_manufactured_enumerated_objective(
    double feed_composition
);

struct Held2PhysicalPhaseBlock {
    double helmholtz_over_rt = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    double pressure_pa = 0.0;
};

struct Held2StateEvaluation {
    std::vector<double> modified_fractions;
    std::vector<double> physical_amounts;
    double volume = 0.0;
    double objective = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    std::vector<double> modified_potentials;
    double pressure_stationarity_relative = 0.0;
};

struct Held2StageICandidate {
    std::vector<double> modified_fractions;
    double volume = 0.0;
    double tpd = 0.0;
};

struct Held2StageIResult {
    std::string outcome;
    int declared_start_count = 0;
    int completed_start_count = 0;
    int failed_start_count = 0;
    std::vector<double> reference_modified_fractions;
    double reference_volume = 0.0;
    double minimum_tpd = 0.0;
    std::vector<Held2StageICandidate> candidates;
};

struct Held2StageIIBound {
    double lower_bound = 0.0;
    double upper_bound = 0.0;
    double multiplier = 0.0;
    int cut_count = 0;
};

struct Held2StageIICandidate {
    std::vector<double> modified_fractions;
    double volume = 0.0;
    double lower_gap = 0.0;
};

struct Held2StageIIAttempt {
    int attempt_id = 0;
    int major_iteration = 0;
    int start_index = 0;
    std::string start_source;
    std::vector<double> internal_start;
    std::vector<double> physical_start_modified_fractions;
    double physical_start_volume = 0.0;
    std::string solver_status = "not_run";
    bool solver_converged = false;
    std::string provider_status = "manufactured_oracle";
    std::string callback_error;
    std::vector<double> internal_terminal;
    std::vector<double> terminal_modified_fractions;
    double terminal_volume = 0.0;
    double objective = 0.0;
    double lower_value = 0.0;
    double pressure_residual = 0.0;
    std::vector<double> lower_bound_multipliers;
    std::vector<double> upper_bound_multipliers;
    double chart_jacobian_condition = 1.0;
    double dual_pullback_inf_norm = 0.0;
    double chart_kkt_inf_norm = 0.0;
    double physical_kkt_inf_norm = 0.0;
    double complementarity_inf_norm = 0.0;
    bool pressure_passed = false;
    bool physical_kkt_passed = false;
    bool cut_eligible = false;
    bool step6_eligible = false;
    int basin_id = -1;
    double same_major_upper_bound = 0.0;
    double same_major_multiplier = 0.0;
};

struct Held2StageIIResult {
    std::string outcome;
    std::string historical_dual_pullback_fixture_status = "not_assigned";
    int major_iterations = 0;
    int lower_starts_per_iteration = 0;
    int cut_count = 0;
    std::vector<Held2StageIIBound> bound_history;
    std::vector<Held2StageIIAttempt> attempt_trace;
    std::vector<Held2StageIICandidate> candidates;
};

struct Held2StageIIINlpEvaluation {
    double objective = 0.0;
    std::vector<double> objective_gradient;
    std::vector<double> constraints;
    std::vector<double> constraint_jacobian;
    std::vector<double> lagrangian_gradient;
    std::vector<double> lagrangian_hessian;
};

struct Held2StageIIIPhase {
    double phase_fraction = 0.0;
    std::vector<double> modified_fractions;
    std::vector<double> physical_fractions;
    double volume = 0.0;
};

struct Held2StageIIIResult {
    std::string solver_status = "not_run";
    std::string numerical_status = "not_adjudicated";
    std::string physical_status = "not_adjudicated";
    std::string feedback = "return_to_stage_ii";
    std::string failure_reason;
    std::string trace_refinement_status = "not_adjudicated";
    int input_candidate_count = 0;
    int retired_duplicate_count = 0;
    int trace_component_count = 0;
    int certified_modified_potential_count = 0;
    double objective = 0.0;
    double modified_balance_inf_norm = 0.0;
    double ordinary_balance_inf_norm = 0.0;
    double phase_charge_inf_norm = 0.0;
    double pressure_stationarity_inf_norm = 0.0;
    double modified_potential_mixed_gap = 0.0;
    double minimum_phase_distance = 0.0;
    double kkt_stationarity_inf_norm = 0.0;
    double enumeration_objective_gap = 0.0;
    std::vector<Held2StageIIIPhase> phases;
};

[[nodiscard]] Held2StateEvaluation evaluate_held2_phase_block(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent_modified_fractions,
    double log_volume,
    double pressure_over_rt,
    double target_pressure_pa,
    const Held2PhysicalPhaseBlock& block
);

[[nodiscard]] Held2StageIResult solve_held2_manufactured_stage_i(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed
);

[[nodiscard]] Held2StageIIResult solve_held2_manufactured_stage_ii(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed
);

[[nodiscard]] Held2StageIIINlpEvaluation evaluate_held2_manufactured_stage_iii_nlp(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::vector<std::array<double, 2>>& candidates,
    const std::vector<double>& variables,
    const std::vector<double>& equality_multipliers
);

[[nodiscard]] Held2StageIIIResult solve_held2_manufactured_stage_iii(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::vector<std::array<double, 2>>& candidates
);

struct Held2Certificate {
    bool accepted = false;
    bool independent_evidence = false;
    double modified_balance_abs = 0.0;
    double ordinary_balance_inf_norm = 0.0;
    double phase_charge_inf_norm = 0.0;
    double modified_potential_gap = 0.0;
    double pressure_stationarity_inf_norm = 0.0;
    double reduced_kkt_inf_norm = 0.0;
    double enumeration_objective_gap = 0.0;
    double independent_modified_composition_count = 0.0;
};

struct Held2ManufacturedEvaluation {
    Held2Coordinates coordinates;
    std::vector<double> modified_feed;
    double phase_fraction = 0.0;
    std::array<std::vector<double>, 2> modified_phases;
    std::array<std::vector<double>, 2> physical_phases;
    std::array<double, 2> phase_charge_residuals{};
    std::vector<double> modified_balance;
    std::vector<double> ordinary_balance;
    std::vector<double> transformed_modified_potentials;
    std::array<std::array<double, 2>, 2> phase_gibbs_gradients{};
    std::array<std::array<double, 2>, 2> phase_modified_potentials{};
    double modified_potential_gap = 0.0;
    double pressure_stationarity_inf_norm = 0.0;
    double objective = 0.0;
    std::array<double, 4> gradient{};
    Held2Certificate certificate;
};

[[nodiscard]] Held2ManufacturedEvaluation evaluate_held2_manufactured(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::array<double, 4>& variables,
    const std::vector<double>& chemical_potentials
);

}  // namespace epcsaft_equilibrium
