#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace epcsaft_equilibrium {

struct DenseMatrix {
    std::size_t rows = 0;
    std::size_t columns = 0;
    std::vector<double> values;

    [[nodiscard]] double operator()(std::size_t row, std::size_t column) const;
    double& operator()(std::size_t row, std::size_t column);
};

struct EquilibriumConstantRecord {
    std::string source_id;
    std::string reference_id;
    std::string reaction_orientation;
    std::string conversion_id;
    bool dimensionless = false;
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
};

struct ReactionSystemInput {
    std::vector<std::string> species_ids;
    std::vector<int> charges;
    std::string provider_fingerprint;
    std::vector<double> molar_masses_kg_per_mol;
    DenseMatrix balance_matrix;
    std::vector<double> conserved_totals;
    DenseMatrix reaction_matrix;
    std::vector<double> feed_amounts;
    std::vector<double> ln_k;
    std::vector<EquilibriumConstantRecord> equilibrium_constant_records;
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
};

struct CompiledReactionSystem {
    std::string provider_fingerprint;
    std::size_t original_species_count = 0;
    std::vector<int> original_charges;
    std::vector<double> original_molar_masses_kg_per_mol;
    std::vector<double> original_feed_amounts;
    std::size_t species_count = 0;
    std::vector<int> charges;
    std::vector<double> molar_masses_kg_per_mol;
    DenseMatrix supplied_balance_matrix;
    std::vector<std::size_t> retained_species_indices;
    std::vector<std::size_t> removed_species_indices;
    DenseMatrix balance_matrix;
    DenseMatrix reaction_matrix;
    DenseMatrix supplied_reaction_transform;
    std::vector<double> balance_totals;
    std::vector<double> feed_amounts;
    std::vector<double> g_ref;
};

struct AmountChart {
    std::vector<int> charges;
    std::vector<std::size_t> cation_indices;
    std::vector<std::size_t> anion_indices;
    std::vector<std::size_t> neutral_indices;

    [[nodiscard]] std::size_t coordinate_count() const;
    [[nodiscard]] bool ionic() const;
};

struct AmountChartEvaluation {
    std::vector<double> amounts;
    std::vector<double> jacobian;
    std::vector<double> amount_hessians;
    double minimum_amount = 0.0;
    double charge_residual = 0.0;
};

struct MaxMinInitializationResult {
    std::string solver_status;
    std::vector<double> amounts;
    std::vector<double> amount_upper_bounds;
    bool strict_positive_feasible = false;
};

struct ChemicalSensitivityResult {
    std::string status = "unavailable";
    std::string failure_reason = "primal_solution_unavailable";
    std::vector<std::string> parameter_order;
    std::vector<double> amount_derivatives;
    std::vector<double> volume_derivatives;
    std::size_t kkt_dimension = 0;
    std::size_t kkt_rank = 0;
    double condition_number_inf = 0.0;
    std::vector<std::size_t> active_lower_bounds;
    std::vector<std::size_t> active_upper_bounds;
    std::vector<std::size_t> active_constraint_bounds;
    std::vector<std::size_t> active_trace_species;
    std::string chart_topology;
    std::string parameter_fingerprint;
    std::string provider_parameter_status = "not_applicable";
    std::string provider_parameter_failure_reason;
    std::string reference_parameter_status = "not_applicable";
    std::string reference_parameter_failure_reason;
};

struct ChemicalSearchAttempt {
    std::size_t ordinal = 0;
    std::size_t primary_ordinal = 0;
    std::string kind = "primary";
    long parent_ordinal = -1;
    std::string start_identity;
    std::string start_construction_status = "not_evaluated";
    std::string retraction_status = "not_needed";
    std::string continuation_status = "not_used";
    std::string provider_domain_status = "not_adjudicated";
    std::string solver_status;
    std::string callback_error;
    std::string terminal_status = "not_evaluated";
    std::vector<double> amounts;
    double volume_m3 = std::numeric_limits<double>::quiet_NaN();
    double objective = std::numeric_limits<double>::quiet_NaN();
    double balance_inf_norm = std::numeric_limits<double>::quiet_NaN();
    double charge_inf_norm = std::numeric_limits<double>::quiet_NaN();
    double pressure_relative_residual = std::numeric_limits<double>::quiet_NaN();
    double reaction_affinity_inf_norm = std::numeric_limits<double>::quiet_NaN();
    double kkt_stationarity_inf_norm = std::numeric_limits<double>::quiet_NaN();
    double complementarity_inf_norm = std::numeric_limits<double>::quiet_NaN();
    std::size_t kkt_dimension = 0;
    std::size_t kkt_rank = 0;
    double condition_number_inf = std::numeric_limits<double>::infinity();
    std::string local_minimum_status = "not_adjudicated";
    std::string trace_status = "not_adjudicated";
    long basin_ordinal = -1;
    std::size_t recovery_seed_count = 0;
    std::size_t recovery_solve_count = 0;
};

struct ChemicalSearchBasin {
    std::size_t ordinal = 0;
    std::size_t representative_attempt_ordinal = 0;
    std::vector<double> amounts;
    double volume_m3 = std::numeric_limits<double>::quiet_NaN();
    double objective = std::numeric_limits<double>::quiet_NaN();
};

struct ChemicalSearchBudgetPrefix {
    std::size_t primary_budget = 0;
    std::vector<std::size_t> attempted_primary_ordinals;
    std::vector<std::size_t> basin_ordinals;
    long selected_basin_ordinal = -1;
    bool selection_changed = false;
};

struct ChemicalSearchEvidence {
    std::string status = "not_evaluated";
    std::string continuation_status = "not_used";
    std::string continuation_blocker;
    std::string continuation_initial_model_fingerprint;
    double continuation_accepted_lambda =
        std::numeric_limits<double>::quiet_NaN();
    std::size_t continuation_attempt_count = 0;
    std::size_t primary_budget = 25;
    std::size_t primary_attempt_count = 0;
    std::size_t generated_start_count = 0;
    std::size_t budget_truncated_start_count = 0;
    std::size_t duplicate_start_count = 0;
    std::size_t infeasible_start_count = 0;
    std::size_t evaluated_start_count = 0;
    std::size_t domain_rejected_start_count = 0;
    std::size_t construction_rejected_start_count = 0;
    std::vector<ChemicalSearchAttempt> attempts;
    std::vector<ChemicalSearchBasin> basins;
    std::vector<ChemicalSearchBudgetPrefix> budget_prefixes;
    long selected_basin_ordinal = -1;
    double selected_objective = std::numeric_limits<double>::quiet_NaN();
    std::string selection_label = "lowest_observed_certified_local_value";
};

struct ChemicalSolveResult {
    bool accepted = false;
    std::string solver_status;
    std::string callback_error;
    std::string chemical_certification_level = "FEASIBLE_ONLY";
    std::string boundary_status = "not_adjudicated";
    std::vector<std::size_t> structural_zero_species_indices;
    std::string numerical_status = "not_adjudicated";
    std::string physical_status = "not_adjudicated";
    std::string provider_domain_status = "not_adjudicated";
    std::string local_minimum_status = "not_adjudicated";
    std::string negative_curvature_recovery_status = "not_needed";
    std::size_t negative_curvature_recovery_seed_count = 0;
    std::size_t negative_curvature_recovery_attempts = 0;
    int negative_curvature_recovery_selected_sign = 0;
    std::string trace_status = "not_adjudicated";
    std::vector<double> amounts;
    double volume_m3 = 0.0;
    double balance_inf_norm = 0.0;
    double charge_inf_norm = 0.0;
    double pressure_relative_residual = 0.0;
    double reaction_affinity_inf_norm = 0.0;
    std::vector<double> reaction_affinity_residuals;
    double packing_fraction = 0.0;
    double packing_fraction_min = std::numeric_limits<double>::quiet_NaN();
    double packing_fraction_max = std::numeric_limits<double>::quiet_NaN();
    double total_ion_fraction = std::numeric_limits<double>::quiet_NaN();
    double total_ion_fraction_max = std::numeric_limits<double>::quiet_NaN();
    double minimum_amount_mol = std::numeric_limits<double>::quiet_NaN();
    double trace_floor_mol = std::numeric_limits<double>::quiet_NaN();
    double kkt_stationarity_inf_norm = 0.0;
    std::vector<double> physical_stationarity_residuals;
    std::vector<double> physical_equality_multipliers;
    std::vector<double> physical_constraint_jacobian;
    std::size_t physical_constraint_rows = 0;
    std::size_t physical_constraint_columns = 0;
    std::vector<double> physical_lagrangian_gradient;
    std::vector<double> physical_to_chart_jacobian;
    std::size_t physical_to_chart_rows = 0;
    std::size_t physical_to_chart_columns = 0;
    double chart_physical_pullback_residual_inf_norm =
        std::numeric_limits<double>::quiet_NaN();
    double complementarity_inf_norm = 0.0;
    std::size_t kkt_dimension = 0;
    std::size_t kkt_rank = 0;
    double condition_number_inf = std::numeric_limits<double>::infinity();
    std::string failure_kind = "not_adjudicated";
    std::string failure_reason;
    std::vector<std::size_t> active_lower_bounds;
    std::vector<std::size_t> active_upper_bounds;
    std::vector<std::size_t> active_constraint_bounds;
    std::string reduced_hessian_status = "not_adjudicated";
    std::vector<double> reduced_hessian;
    std::vector<double> reduced_hessian_nullspace_basis;
    std::size_t reduced_hessian_nullspace_rows = 0;
    std::size_t reduced_hessian_nullspace_columns = 0;
    std::vector<double> reduced_hessian_eigenvalues;
    std::string reduced_hessian_spectrum_status = "not_adjudicated";
    std::size_t reduced_hessian_raw_positive_eigenvalues = 0;
    std::size_t reduced_hessian_raw_zero_eigenvalues = 0;
    std::size_t reduced_hessian_raw_negative_eigenvalues = 0;
    std::size_t reduced_hessian_positive_eigenvalues = 0;
    std::size_t reduced_hessian_zero_eigenvalues = 0;
    std::size_t reduced_hessian_negative_eigenvalues = 0;
    double reduced_hessian_scale = std::numeric_limits<double>::quiet_NaN();
    double reduced_hessian_eigenvalue_tolerance =
        std::numeric_limits<double>::quiet_NaN();
    std::vector<double> objective_gradient;
    std::vector<double> constraint_values;
    std::vector<double> constraint_jacobian;
    std::vector<double> lagrangian_gradient;
    std::vector<double> equality_multipliers;
    double chart_stationarity_inf_norm = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> lagrangian_hessian;
    std::vector<double> covariant_lagrangian_hessian;
    double derivative_check_step = std::numeric_limits<double>::quiet_NaN();
    double objective_gradient_check_relative_error =
        std::numeric_limits<double>::quiet_NaN();
    double constraint_jacobian_check_relative_error =
        std::numeric_limits<double>::quiet_NaN();
    double lagrangian_hessian_check_relative_error =
        std::numeric_limits<double>::quiet_NaN();
    std::string derivative_check_worst_entry;
    double derivative_check_worst_relative_error =
        std::numeric_limits<double>::quiet_NaN();
    double derivative_check_worst_analytic_value =
        std::numeric_limits<double>::quiet_NaN();
    double derivative_check_worst_finite_difference_value =
        std::numeric_limits<double>::quiet_NaN();
    double derivative_check_worst_step =
        std::numeric_limits<double>::quiet_NaN();
    std::vector<std::string> derivative_coordinate_order;
    std::string derivative_objective_basis =
        "dimensionless_fixed_TP_A_plus_PV_plus_reference_over_RT";
    std::string derivative_constraint_basis =
        "ordered_per_row_in_derivative_constraint_order";
    std::vector<std::string> derivative_constraint_order;
    std::vector<double> kkt_root_jacobian;
    std::size_t kkt_root_rows = 0;
    std::size_t kkt_root_columns = 0;
    std::string kkt_root_status = "not_evaluated";
    double kkt_root_jacobian_check_relative_error =
        std::numeric_limits<double>::quiet_NaN();
    ChemicalSearchEvidence search;
    ChemicalSensitivityResult sensitivities;
};

struct ManufacturedNlpEvaluation {
    double objective = 0.0;
    std::vector<double> objective_gradient;
    std::vector<double> constraints;
    std::vector<double> constraint_jacobian;
    std::vector<double> lagrangian_gradient;
    std::vector<double> lagrangian_hessian;
    std::vector<double> amounts;
    double volume_m3 = 0.0;
    // Populated only by the private inverse-chart evidence seam.
    std::vector<double> kkt_backtransform_rhs;
    std::vector<double> kkt_backtransform_solution;
};

struct ManufacturedReducedHessianEvidence {
    bool positive = false;
    std::string status = "second_order_inconclusive";
    double curvature = 0.0;
    std::vector<double> negative_direction;
    std::vector<double> reduced_hessian;
    std::vector<double> nullspace_basis;
    std::size_t nullspace_rows = 0;
    std::size_t nullspace_columns = 0;
    std::vector<double> eigenvalues;
    std::array<std::size_t, 3> inertia = {0, 0, 0};
    std::array<std::size_t, 3> raw_inertia = {0, 0, 0};
    std::string spectrum_status = "not_evaluated";
    double hessian_scale = std::numeric_limits<double>::quiet_NaN();
    double eigenvalue_tolerance = std::numeric_limits<double>::quiet_NaN();
};

struct ProviderPhaseBlockEvidence {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    double pressure_pa = 0.0;
    double packing_fraction = 0.0;
    std::vector<double> packing_gradient;
    std::vector<double> packing_hessian;
};

class ProviderContext;
struct ProviderActiveParameterSet;

struct SourceReferenceTransferEvidence {
    std::string status;
    std::string domain_status;
    std::string convergence_status;
    std::string reference_state_id;
    std::string activity_convention_id;
    std::vector<std::string> component_ids;
    std::size_t neutral_basis_row_count = 0;
    std::vector<double> neutral_basis;
    std::vector<double> log_fugacity_contractions;
    std::vector<double> activity_scale_log_contractions;
    std::vector<double> transfer_log_contractions;
    std::vector<double> source_composition;
    std::vector<std::string> derivative_availability;
    std::vector<double> pressure_derivatives_per_pa;
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    double source_reference_pressure_pa = 0.0;
    double standard_molality_mol_per_kg = 0.0;
    double reference_convergence_error = 0.0;
    std::string parameter_fingerprint;
    std::string topology_fingerprint;
    std::string component_order_fingerprint;
    std::string reference_state_fingerprint;
    std::string domain_fingerprint;
    std::string artifact_fingerprint;
    std::string helmholtz_basis_id;
};

struct SourceStandardStateResult {
    std::vector<double> standard_offsets;
    std::vector<double> ln_k_provider_basis;
    std::vector<double> pressure_derivatives_per_pa;
    std::vector<double> parameter_derivatives;
    std::size_t active_parameter_count = 0;
    double representation_residual_inf_norm = 0.0;
};

[[nodiscard]] CompiledReactionSystem compile_reaction_system(
    const ReactionSystemInput& input
);

[[nodiscard]] AmountChart make_amount_chart(const std::vector<int>& charges);

[[nodiscard]] AmountChartEvaluation evaluate_amount_chart(
    const AmountChart& chart,
    const std::vector<double>& coordinates
);

[[nodiscard]] std::vector<double> invert_amount_chart(
    const AmountChart& chart,
    const std::vector<double>& amounts
);

[[nodiscard]] MaxMinInitializationResult max_min_initialization(
    const DenseMatrix& balance_matrix,
    const std::vector<double>& feed_amounts,
    const std::vector<int>& charges,
    double trace_floor,
    double total_ion_fraction_max
);

[[nodiscard]] std::vector<std::size_t> homogeneous_structural_zeros(
    const DenseMatrix& balance_matrix,
    const std::vector<double>& balance_totals,
    const std::vector<int>& charges
);

[[nodiscard]] ChemicalSolveResult solve_manufactured_ideal_reaction(
    const CompiledReactionSystem& system,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& gauge_coefficients,
    double trace_floor,
    int max_iterations = 500
);

// Private manufactured nonconvex case used only to verify the generic
// negative-curvature recovery path. It is not a public model or capability.
[[nodiscard]] ChemicalSolveResult solve_manufactured_nonconvex_reaction(
    const CompiledReactionSystem& system,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& gauge_coefficients,
    double trace_floor,
    int max_iterations = 500,
    double quadratic_strength = 2.3
);

[[nodiscard]] ChemicalSolveResult solve_manufactured_inconsistent_derivative_reaction(
    const CompiledReactionSystem& system,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& gauge_coefficients,
    double trace_floor,
    int max_iterations = 200
);

[[nodiscard]] ManufacturedNlpEvaluation evaluate_manufactured_reaction_nlp(
    const CompiledReactionSystem& system,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& gauge_coefficients,
    const std::vector<double>& variables,
    const std::vector<double>& constraint_multipliers
);

// Private manufactured seams used to test the generic recovery mathematics.
[[nodiscard]] ManufacturedReducedHessianEvidence analyze_manufactured_reduced_hessian(
    const std::vector<double>& hessian,
    const std::vector<double>& constraint_jacobian = {},
    std::size_t constraint_count = 0
);

[[nodiscard]] std::vector<double> retract_manufactured_balance(
    const CompiledReactionSystem& system,
    const std::vector<double>& seed,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    double trace_floor
);

[[nodiscard]] std::vector<double> manufactured_recovery_displacement(
    const std::vector<double>& variables,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    const std::vector<double>& direction,
    int sign,
    std::size_t backtrack_index = 0
);

// Private derivative evidence seam. This is intentionally not part of the
// public equilibrium API; tests use it to exercise the inverse log-packing
// chart with the same exact chain rules and KKT pullback as Provider phases.
[[nodiscard]] ManufacturedNlpEvaluation
evaluate_manufactured_inverse_log_packing_nlp(
    const CompiledReactionSystem& system,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& gauge_coefficients,
    const std::vector<double>& variables,
    const std::vector<double>& constraint_multipliers,
    bool zero_kkt_rhs = false
);

[[nodiscard]] ProviderPhaseBlockEvidence evaluate_provider_phase_block(
    const ProviderContext& provider,
    double temperature_k,
    const std::vector<double>& amounts,
    double volume_m3
);

[[nodiscard]] SourceStandardStateResult transform_source_standard_state(
    const DenseMatrix& reaction_matrix,
    const std::vector<double>& source_ln_k,
    const std::vector<double>& log_activity_scale_factors,
    const std::vector<std::string>& species_ids,
    const std::vector<int>& charges,
    const std::string& provider_fingerprint,
    double temperature_k,
    double pressure_pa,
    const SourceReferenceTransferEvidence& reference
);

[[nodiscard]] ChemicalSolveResult solve_provider_reaction(
    const CompiledReactionSystem& system,
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    double packing_fraction_min,
    double packing_fraction_max,
    double total_ion_fraction_max,
    double trace_floor,
    const std::vector<double>& ln_k_pressure_derivatives_per_pa = {},
    const std::vector<double>& ln_k_parameter_derivatives = {},
    const ProviderActiveParameterSet* active_parameters = nullptr,
    int max_iterations = 500
);

[[nodiscard]] ChemicalSolveResult solve_provider_reaction_continuation(
    const CompiledReactionSystem& target_system,
    const ProviderContext& target_provider,
    double target_packing_fraction_min,
    double target_packing_fraction_max,
    double target_total_ion_fraction_max,
    const CompiledReactionSystem& initial_system,
    const ProviderContext& initial_provider,
    double initial_packing_fraction_min,
    double initial_packing_fraction_max,
    double initial_total_ion_fraction_max,
    double temperature_k,
    double pressure_pa,
    double trace_floor,
    int max_iterations = 500,
    bool continue_after_certified_target = false
);

}  // namespace epcsaft_equilibrium
