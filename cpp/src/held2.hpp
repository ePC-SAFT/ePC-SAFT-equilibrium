#pragma once

#include <array>
#include <cstddef>
#include <functional>
#include <string>
#include <tuple>
#include <vector>

namespace epcsaft_equilibrium {

inline constexpr const char* kHeld2ManufacturedFormulationId =
    "perdomo-held2.modified-mole.manufactured.v1";
inline constexpr double kHeld2PackingFractionMinimum = 1.0e-6;
inline constexpr double kHeld2PackingFractionMaximum = 0.74;
inline constexpr double kHeld2ModifiedLowerScale = 1.0e-10;

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

[[nodiscard]] double held2_manufactured_enumerated_objective(
    double feed_composition
);

struct Held2PhysicalPhaseBlock {
    double helmholtz_over_rt = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    double pressure_pa = 0.0;
};

struct Held2PackingEvaluation {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
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
    double log_volume_gradient = 0.0;
    bool has_packing_evaluation = false;
    Held2PackingEvaluation packing;
};

using Held2StateEvaluator = std::function<Held2StateEvaluation(
    const std::vector<double>&,
    double
)>;

using Held2VolumeBoundsEvaluator = std::function<std::array<double, 2>(
    const std::vector<double>&
)>;

struct Held2StageICandidate {
    std::vector<double> modified_fractions;
    double volume = 0.0;
    double tpd = 0.0;
    std::array<double, 2> molar_volume_bounds{};
    double pressure_stationarity_relative = 0.0;
    double volume_gradient = 0.0;
    double packing_fraction = 0.0;
    bool lower_volume_bound_active = false;
    bool upper_volume_bound_active = false;
};

struct Held2ReferenceRoot {
    double log_volume = 0.0;
    double volume = 0.0;
    double objective = 0.0;
    double pressure_residual = 0.0;
    double curvature = 0.0;
    bool mechanically_stable = false;
};

struct Held2StageIResult {
    std::string outcome;
    int reference_scan_interval_count = 0;
    int reference_scan_point_count = 0;
    int reference_root_count = 0;
    int reference_stable_root_count = 0;
    int reference_evaluation_failure_count = 0;
    int reference_refinement_failure_count = 0;
    int declared_start_count = 0;
    int attempted_start_count = 0;
    int completed_start_count = 0;
    int failed_start_count = 0;
    int candidate_domain_evaluation_failure_count = 0;
    int candidate_domain_rejection_count = 0;
    int failed_start_index = -1;
    int failed_start_solver_status = 999;
    bool failed_start_solver_converged = false;
    std::string failed_start_reason;
    std::vector<double> failed_start_initial;
    bool volume_domain_search_complete = true;
    std::vector<double> reference_modified_fractions;
    double reference_volume = 0.0;
    double minimum_tpd = 0.0;
    std::vector<Held2ReferenceRoot> reference_roots;
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
    std::vector<double> independent_modified_fractions;
    double volume = 0.0;
    double phase_coordinate = 0.0;
    double lower_gap = 0.0;
};

struct Held2StageIIAttempt {
    int start_index = -1;
    int solver_status = 999;
    int iterations = -1;
    bool solver_converged = false;
    bool numerical_certified = false;
    bool provider_terminal_valid = false;
    bool improving = false;
    std::string callback_error;
    double lower_value = 0.0;
    double projected_kkt_inf_norm = 0.0;
    double constraint_violation = 0.0;
    double complementarity = 0.0;
    double final_step_norm = 0.0;
};

struct Held2StageIIResult {
    std::string outcome;
    int major_iterations = 0;
    int lower_starts_per_iteration = 0;
    int lower_attempted_start_count = 0;
    int lower_completed_start_count = 0;
    int lower_failed_start_count = 0;
    int first_failed_lower_start_index = -1;
    int first_failed_lower_solver_status = 999;
    std::string first_failed_lower_reason;
    std::vector<double> first_failed_lower_initial;
    int certified_improving_cut_count = 0;
    std::string search_completeness = "not_run";
    bool final_lower_search_complete = false;
    int cut_count = 0;
    std::vector<Held2StageIIBound> bound_history;
    std::vector<Held2StageIIAttempt> attempt_log;
    std::vector<Held2StageIICandidate> candidates;
};

struct Held2StageIILowerDecision {
    std::string decision;
    std::string search_completeness;
    bool lower_bound_certified = false;
};

[[nodiscard]] Held2StageIILowerDecision decide_held2_stage_ii_lower(
    double upper_bound,
    double best_certified_value,
    int certified_start_count,
    int declared_start_count
);

[[nodiscard]] std::vector<Held2StageIICandidate>
evaluate_held2_stage_ii_step6_test_adapter(
    const Held2Coordinates& coordinates,
    const std::vector<double>& feed_independent_modified_fractions,
    double upper_bound,
    const std::vector<double>& multipliers,
    const std::vector<Held2StateEvaluation>& states,
    const std::vector<std::vector<double>>& independent_modified_fractions,
    const std::vector<std::vector<double>>& fixed_volume_composition_gradients,
    const std::vector<double>& phase_coordinates
);

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

[[nodiscard]] Held2PackingEvaluation evaluate_held2_packing_block(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent_modified_fractions,
    double log_volume,
    double packing_fraction,
    const std::vector<double>& gradient,
    const std::vector<double>& hessian
);

[[nodiscard]] Held2StateEvaluation evaluate_held2_log_packing_state(
    const Held2StateEvaluator& log_volume_evaluator,
    const std::vector<double>& independent_modified_fractions,
    double log_packing_fraction,
    const std::array<double, 2>& molar_volume_bounds
);

[[nodiscard]] Held2StageIResult solve_held2_manufactured_stage_i(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed
);

[[nodiscard]] Held2StageIResult solve_held2_stage_i(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const Held2StateEvaluator& log_volume_evaluator,
    const Held2StateEvaluator& search_evaluator,
    const std::array<double, 2>& molar_volume_bounds,
    const Held2VolumeBoundsEvaluator& volume_bounds_evaluator,
    bool search_uses_log_packing,
    bool volume_domain_search_complete,
    int solve_start_limit
);

[[nodiscard]] Held2StageIIResult solve_held2_manufactured_stage_ii(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed
);

[[nodiscard]] Held2StageIIResult solve_held2_stage_ii(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const Held2StateEvaluator& search_evaluator,
    const Held2StateEvaluator& fixed_volume_evaluator,
    const Held2StateEvaluation& reference,
    const std::vector<Held2StageICandidate>& stage_i_candidates
);

[[nodiscard]] Held2StageIIINlpEvaluation evaluate_held2_manufactured_stage_iii_nlp(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::vector<std::array<double, 2>>& candidates,
    const std::vector<double>& variables,
    const std::vector<double>& equality_multipliers
);

[[nodiscard]] Held2StageIIINlpEvaluation evaluate_held2_stage_iii_nlp(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const Held2StateEvaluator& evaluator,
    std::size_t phase_count,
    const std::vector<double>& variables,
    const std::vector<double>& equality_multipliers
);

[[nodiscard]] std::tuple<
    int,
    int,
    int,
    std::vector<double>,
    std::vector<double>,
    std::vector<int>,
    std::vector<int>,
    std::vector<double>> inspect_held2_stage_iii_tnlp(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const std::vector<Held2StageIICandidate>& candidates,
    const std::array<double, 2>& phase_coordinate_bounds,
    const std::vector<double>& variables
);

[[nodiscard]] std::tuple<bool, int, std::string, std::string>
probe_held2_stage_iii_objective_trial(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const std::vector<Held2StageIICandidate>& candidates,
    const Held2StateEvaluator& evaluator,
    const std::array<double, 2>& phase_coordinate_bounds,
    const std::vector<double>& variables
);

[[nodiscard]] Held2StageIIIResult solve_held2_stage_iii(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const std::vector<Held2StageIICandidate>& candidates,
    const Held2StateEvaluator& evaluator,
    const std::array<double, 2>& phase_coordinate_bounds
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
