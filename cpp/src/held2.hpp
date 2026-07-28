#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace epcsaft_equilibrium {

class Held2ProgressObserver;

inline constexpr double kHeld2PackingFractionMinimum = 1.0e-6;
inline constexpr double kHeld2PackingFractionMaximum = 0.74;
inline constexpr double kHeld2ModifiedLowerScale = 1.0e-10;

enum class Held2CompositionDomain {
    FiniteSearch,
    TraceRefinement,
};

struct Held2PolytopeConstraint {
    std::string name;
    std::vector<double> coefficients;
    double upper_bound = 0.0;
};

struct Held2Coordinates {
    std::vector<double> charges;
    std::size_t eliminated_index = 0;
    std::size_t dependent_index = 0;
    std::vector<std::size_t> paper_to_provider_indices;
    std::vector<std::size_t> provider_to_paper_indices;
    std::vector<std::size_t> compact_to_paper_indices;
    std::vector<std::size_t> retained_indices;
    std::vector<std::size_t> independent_indices;
    std::vector<double> modified_factors;
    std::vector<double> independent_lower_bounds;
    std::vector<double> independent_upper_bounds;
    std::vector<Held2PolytopeConstraint> polytope_constraints;
};

struct Held2StepTiming {
    int step = 0;
    int invocation_count = 0;
    double wall_seconds = 0.0;
    double cpu_seconds = 0.0;
    std::uint64_t provider_evaluations = 0;
    std::uint64_t optimizer_solves = 0;
    std::uint64_t optimizer_iterations = 0;
    std::string terminal_status = "not_run";
    std::string terminal_reason = "not_run";
    int next_step = 0;
};

struct Held2StateEvaluation;

[[nodiscard]] Held2Coordinates make_held2_coordinates(
    const std::vector<double>& charges,
    const std::vector<std::string>& component_ids = {}
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
    const std::vector<double>& independent_modified_fractions,
    Held2CompositionDomain domain = Held2CompositionDomain::FiniteSearch
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

struct Held2PhysicalPhaseBlock {
    double helmholtz_over_rt = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    double pressure_pa = 0.0;
};

using Held2PhysicalVolumeBoundsEvaluator =
    std::function<std::array<double, 2>(const std::vector<double>&)>;

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
    double helmholtz_over_rt_reference_amount = 0.0;
    double pressure_pa = 0.0;
    std::vector<double> chemical_potentials_over_rt;
};

using Held2StateEvaluator = std::function<Held2StateEvaluation(
    const std::vector<double>&,
    double
)>;
using Held2StateValueEvaluator = std::function<double(
    const std::vector<double>&,
    double
)>;

using Held2VolumeBoundsEvaluator = std::function<std::array<double, 2>(
    const std::vector<double>&
)>;

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

struct Held2StageICandidate {
    std::vector<double> modified_fractions;
    double volume = 0.0;
    double tpd = 0.0;
};

struct Held2StageIICandidate {
    std::vector<double> modified_fractions;
    std::vector<double> independent_modified_fractions;
    double volume = 0.0;
    double phase_coordinate = 0.0;
    double lower_gap = 0.0;
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
    int retired_inactive_count = 0;
    int stage_iii_solve_count = 0;
    int optimizer_iteration_count = 0;
    int active_set_resolve_count = 0;
    int pressure_polish_iteration_count = 0;
    std::string pressure_polish_status = "not_run";
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
    double free_energy_upper_bound =
        std::numeric_limits<double>::quiet_NaN();
    double free_energy_gap =
        std::numeric_limits<double>::quiet_NaN();
    std::string free_energy_gap_provenance = "unavailable";
    bool kkt_evidence_available = false;
    bool physical_evidence_available = false;
    bool phase_identity_evidence_available = false;
    bool free_energy_gap_available = false;
    std::vector<Held2StageIIIPhase> phases;
    std::vector<double> solution_variables;
};

[[nodiscard]] Held2StateEvaluation evaluate_held2_phase_block(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent_modified_fractions,
    double log_volume,
    double pressure_over_rt,
    double target_pressure_pa,
    const Held2PhysicalPhaseBlock& block,
    Held2CompositionDomain domain = Held2CompositionDomain::FiniteSearch
);

[[nodiscard]] Held2StageIIIResult solve_held2_stage_iii(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const std::vector<Held2StageIICandidate>& candidates,
    const Held2StateEvaluator& evaluator,
    const std::vector<std::array<double, 2>>& phase_coordinate_bounds,
    double free_energy_upper_bound,
    const std::string& free_energy_gap_provenance,
    std::vector<double> initial = {},
    const Held2StateValueEvaluator& value_evaluator = {}
);

}  // namespace epcsaft_equilibrium
