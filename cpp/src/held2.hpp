#pragma once

#include <array>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace epcsaft_equilibrium {

class Held2ProgressObserver;

inline constexpr const char* kHeld2ManufacturedFormulationId =
    "perdomo-held2.modified-mole.manufactured.v1";
inline constexpr double kHeld2PackingFractionMinimum = 1.0e-6;
inline constexpr double kHeld2PackingFractionMaximum = 0.74;

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

[[nodiscard]] std::vector<double> held2_map_unit_cube_to_independent_fractions(
    const Held2Coordinates& coordinates,
    const std::vector<double>& unit_cube_coordinates,
    double total_ion_mole_fraction_max
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
    double pressure_stationarity_derivative_log_volume = 0.0;
};

using Held2StateEvaluator = std::function<Held2StateEvaluation(
    const std::vector<double>&,
    double
)>;

using Held2VolumeBoundsEvaluator = std::function<std::array<double, 2>(
    const std::vector<double>&
)>;

struct Held2StageIIPressureRootReduction {
    double objective = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    double pressure_coordinate_gradient = 0.0;
    double pressure_coordinate_curvature = 0.0;
};

[[nodiscard]] Held2StageIIPressureRootReduction
reduce_held2_stage_ii_pressure_root(
    const std::vector<double>& independent,
    const std::vector<double>& feed,
    const std::vector<double>& multipliers,
    const Held2StateEvaluation& state
);

struct Held2StageIIStep5Assessment {
    bool qualified = false;
    double gap = 0.0;
    std::string reason;
};

[[nodiscard]] Held2StageIIStep5Assessment assess_held2_stage_ii_step5(
    double upper_bound,
    double local_value,
    bool local_state_certified
);

struct Held2PressureScanPoint {
    double log_volume = 0.0;
    double volume = 0.0;
    double pressure_residual = 0.0;
    double pressure_derivative_log_volume = 0.0;
    double objective = 0.0;
    bool valid = false;
    std::string failure;
};

struct Held2PressureScanInterval {
    double lower_log_volume = 0.0;
    double upper_log_volume = 0.0;
    int depth = 0;
    std::string status;
};

struct Held2PressureRoot {
    double log_volume = 0.0;
    double volume = 0.0;
    double objective = 0.0;
    double pressure_residual = 0.0;
    double pressure_derivative_log_volume = 0.0;
    double objective_curvature_log_volume = 0.0;
    std::string mechanical_class = "unclassified";
    std::string origin;
    bool boundary = false;
    Held2StateEvaluation state;
};

struct Held2PressureEnvelopeResult {
    std::string outcome = "indeterminate";
    std::string failure_reason;
    std::string root_completeness = "not_proven";
    std::string selection_scope =
        "lowest_among_detected_strict_stable_roots";
    int selected_root_index = -1;
    int evaluation_failure_count = 0;
    int refinement_failure_count = 0;
    int stationary_point_count = 0;
    int tangential_root_count = 0;
    int marginal_root_count = 0;
    int boundary_root_count = 0;
    int objective_tie_count = 0;
    int deduplicated_root_count = 0;
    double lower_log_volume = 0.0;
    double upper_log_volume = 0.0;
    std::vector<Held2PressureScanPoint> scan_points;
    std::vector<Held2PressureScanInterval> intervals;
    std::vector<Held2PressureRoot> roots;
};

[[nodiscard]] Held2PressureEnvelopeResult evaluate_held2_pressure_envelope(
    const std::vector<double>& independent_modified_fractions,
    const std::array<double, 2>& molar_volume_bounds,
    const Held2StateEvaluator& evaluator,
    int initial_interval_count,
    int maximum_subdivision_depth = 8
);

[[nodiscard]] Held2PressureEnvelopeResult
evaluate_held2_manufactured_pressure_envelope(
    const std::string& topology,
    double composition,
    int initial_interval_count
);

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
    bool lower_bound_available = false;
    double upper_bound = 0.0;
    std::vector<double> multipliers;
    int cut_count = 0;
    std::string upper_solver;
    std::string upper_solver_version;
    std::string upper_solver_status;
    bool upper_primal_feasible = false;
    bool upper_dual_feasible = false;
    double upper_primal_residual_inf = 0.0;
    double upper_primal_scale = 0.0;
    double upper_dual_residual_inf = 0.0;
    double upper_dual_scale = 0.0;
    double upper_complementarity_inf = 0.0;
    std::vector<double> cut_slacks;
    std::vector<double> cut_duals;
    std::vector<int> active_cut_ids;
};

struct Held2StageIILocalEvaluation {
    int sequence = 0;
    std::string callback;
    bool new_x = false;
    bool callback_succeeded = false;
    bool accepted_iterate = false;
    std::vector<double> raw_variables;
    std::vector<double> physical_variables;
    double maximum_bound_violation = 0.0;
    std::string mapping_status = "not_evaluated";
    std::string error;
};

struct Held2StageIIChartCoordinate {
    double raw = 0.0;
    double normalized = 0.0;
    bool normalized_boundary_contact = false;
};

[[nodiscard]] Held2StageIIChartCoordinate normalize_held2_stage_ii_chart_coordinate(
    double coordinate
);

struct Held2StageIICandidate {
    std::vector<double> modified_fractions;
    std::vector<double> independent_modified_fractions;
    double volume = 0.0;
    double phase_coordinate = 0.0;
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
    std::string provider_status = "not_evaluated";
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
    double dual_pullback_scale = 0.0;
    double chart_kkt_inf_norm = 0.0;
    double primal_inf_norm = 0.0;
    double dual_sign_violation_inf_norm = 0.0;
    double physical_kkt_inf_norm = 0.0;
    double complementarity_inf_norm = 0.0;
    bool pressure_passed = false;
    bool dual_signs_valid = false;
    bool physical_kkt_passed = false;
    bool cut_eligible = false;
    bool step5_qualified = false;
    std::string step5_reason = "not_evaluated";
    bool step6_eligible = false;
    double step6_gap = 0.0;
    bool step6_gap_passed = false;
    bool step6_gradient_passed = false;
    std::string step6_rejection_reason = "not_evaluated";
    double fixed_volume_gradient_inf_norm = 0.0;
    double fixed_volume_gradient_scale = 0.0;
    int basin_id = -1;
    double same_major_upper_bound = 0.0;
    std::vector<double> same_major_multipliers;
    std::vector<Held2StageIILocalEvaluation> evaluation_trace;
    std::vector<double> last_valid_physical_variables;
};

struct Held2StageIIResult {
    std::string outcome;
    std::string search_strategy = "continuation_sobol_direct_l_ipopt_v1";
    std::string global_explorer = "continuation_sobol_direct_l";
    std::string local_solver = "ipopt_exact_hessian";
    std::string globality_certificate = "not_guaranteed";
    int major_iterations = 0;
    int lower_starts_per_iteration = 0;
    int cut_count = 0;
    int exploration_evaluation_count = 0;
    int exploration_failure_count = 0;
    int exploration_representative_count = 0;
    int duplicate_representative_count = 0;
    int duplicate_terminal_count = 0;
    int distinct_basin_count = 0;
    int unresolved_candidate_identity_count = 0;
    int local_attempt_cap_per_major = 0;
    bool local_attempts_truncated = false;
    bool direct_escalation_used = false;
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

struct Held2StageIIIRetirementDecision {
    bool retire = false;
    std::string reason = "not_adjudicated";
    double complementarity_inf_norm = 0.0;
    double stationarity_residual = 0.0;
};

struct Held2StageIIILifecycleStep {
    int solve_index = 0;
    int active_candidate_count = 0;
    int removed_candidate_index = -1;
    std::string action;
    double phase_fraction = 0.0;
    double lower_bound_multiplier = 0.0;
    double reduced_derivative = 0.0;
    double complementarity_inf_norm = 0.0;
    std::vector<double> candidate_independent_modified_fractions;
    double candidate_volume = 0.0;
    std::string solver_status;
    std::string decision_reason;
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
    int retired_inactive_count = 0;
    int stage_iii_solve_count = 0;
    int active_set_resolve_count = 0;
    int trace_component_count = 0;
    int certified_modified_potential_count = 0;
    double objective = 0.0;
    double modified_balance_inf_norm = 0.0;
    double ordinary_balance_inf_norm = 0.0;
    double phase_charge_inf_norm = 0.0;
    double phase_charge_scale = 0.0;
    double pressure_stationarity_inf_norm = 0.0;
    double modified_potential_mixed_gap = 0.0;
    double modified_potential_scale = 0.0;
    double minimum_phase_distance = 0.0;
    std::string phase_identity_status = "not_adjudicated";
    double kkt_stationarity_inf_norm = 0.0;
    double dual_sign_violation_inf_norm = 0.0;
    double bound_complementarity_inf_norm = 0.0;
    double minimum_phase_fraction = 0.0;
    double enumeration_objective_gap = 0.0;
    bool kkt_evidence_available = false;
    bool physical_evidence_available = false;
    bool phase_identity_evidence_available = false;
    bool free_energy_gap_available = false;
    std::vector<Held2StageIIIPhase> phases;
    std::vector<Held2StageIIILifecycleStep> lifecycle;
};

[[nodiscard]] Held2StageIIIRetirementDecision held2_stage_iii_retirement_decision(
    double phase_fraction,
    double lower_bound_multiplier,
    double upper_bound_multiplier,
    double reduced_derivative,
    bool remaining_balance_feasible
);

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

[[nodiscard]] Held2StageIIResult solve_held2_stage_ii(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const Held2StateEvaluator& evaluator,
    const Held2VolumeBoundsEvaluator& volume_bounds_evaluator,
    const Held2StateEvaluation& reference,
    const std::vector<Held2StageICandidate>& stage_i_candidates,
    double total_ion_mole_fraction_max,
    int major_iteration_cap,
    int local_attempt_cap_per_major,
    Held2ProgressObserver* observer = nullptr
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

[[nodiscard]] Held2StageIIINlpEvaluation evaluate_held2_stage_iii_nlp(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const Held2StateEvaluator& evaluator,
    std::size_t phase_count,
    const std::vector<double>& variables,
    const std::vector<double>& equality_multipliers
);

[[nodiscard]] Held2StageIIIResult solve_held2_stage_iii(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const std::vector<Held2StageIICandidate>& candidates,
    const Held2StateEvaluator& evaluator,
    const std::vector<std::array<double, 2>>& phase_coordinate_bounds,
    Held2ProgressObserver* observer = nullptr
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
