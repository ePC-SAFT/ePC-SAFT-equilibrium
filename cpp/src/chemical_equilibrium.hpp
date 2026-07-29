#pragma once

#include <cstddef>
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
    std::size_t negative_curvature_recovery_attempts = 0;
    int negative_curvature_recovery_selected_sign = 0;
    std::string trace_status = "not_adjudicated";
    std::vector<double> amounts;
    double volume_m3 = 0.0;
    double balance_inf_norm = 0.0;
    double charge_inf_norm = 0.0;
    double pressure_relative_residual = 0.0;
    double reaction_affinity_inf_norm = 0.0;
    double packing_fraction = 0.0;
    double kkt_stationarity_inf_norm = 0.0;
    double complementarity_inf_norm = 0.0;
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
    double curvature = 0.0;
    std::vector<double> negative_direction;
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
struct NeutralReferenceEvaluation;

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
    int max_iterations = 500
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
[[nodiscard]] ManufacturedReducedHessianEvidence
analyze_manufactured_reduced_hessian(const std::vector<double>& hessian);

[[nodiscard]] std::vector<double> manufactured_recovery_displacement(
    const std::vector<double>& variables,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    const std::vector<double>& direction,
    int sign
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
    const std::vector<int>& charges,
    const std::string& provider_fingerprint,
    double temperature_k,
    double pressure_pa,
    const NeutralReferenceEvaluation& reference
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

}  // namespace epcsaft_equilibrium
