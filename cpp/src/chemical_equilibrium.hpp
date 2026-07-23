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
    bool dimensionless = false;
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
};

struct ReactionSystemInput {
    std::vector<std::string> species_ids;
    std::vector<int> charges;
    std::string provider_fingerprint;
    DenseMatrix balance_matrix;
    DenseMatrix reaction_matrix;
    std::vector<double> feed_amounts;
    std::vector<double> ln_k;
    std::vector<EquilibriumConstantRecord> equilibrium_constant_records;
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
};

struct CompiledReactionSystem {
    std::vector<std::string> species_ids;
    std::vector<int> charges;
    DenseMatrix balance_matrix;
    DenseMatrix reaction_matrix;
    std::vector<double> balance_totals;
    std::vector<double> feed_amounts;
    std::vector<double> ln_k;
    std::vector<double> g_ref;
    std::string provider_fingerprint;
    std::size_t balance_rank = 0;
    std::size_t reaction_rank = 0;
    double reaction_qr_diagonal_ratio = 0.0;
    double reference_reconstruction_inf_norm = 0.0;
    double conservation_reaction_inf_norm = 0.0;
    double charge_reaction_inf_norm = 0.0;
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
    std::string reason;
    std::vector<double> amounts;
    std::vector<double> amount_upper_bounds;
    double max_min_amount = 0.0;
    double equality_inf_norm = 0.0;
    bool strict_positive_feasible = false;
};

struct ChemicalSolveResult {
    bool accepted = false;
    std::string solver_status;
    std::string callback_error;
    std::string numerical_status = "not_adjudicated";
    std::string physical_status = "not_adjudicated";
    std::string provider_domain_status = "not_adjudicated";
    std::string local_minimum_status = "not_adjudicated";
    std::string trace_status = "not_adjudicated";
    std::string predictive_status = "not_adjudicated";
    std::string finite_search_status = "not_applicable_single_phase_local_nlp";
    std::vector<double> amounts;
    double volume_m3 = 0.0;
    double objective = 0.0;
    double balance_inf_norm = 0.0;
    double charge_inf_norm = 0.0;
    double pressure_relative_residual = 0.0;
    double reaction_affinity_inf_norm = 0.0;
    double packing_fraction = 0.0;
    double kkt_stationarity_inf_norm = 0.0;
    double complementarity_inf_norm = 0.0;
    double final_lambda = 0.0;
    bool has_final_lambda = false;
    bool continuation_used = false;
    std::string kkt_scope = "not_adjudicated";
    std::vector<double> kkt_residual;
    std::vector<double> kkt_jacobian;
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
};

struct ProviderPhaseBlockEvidence {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    double pressure_pa = 0.0;
    double packing_fraction = 0.0;
    std::vector<double> packing_gradient;
    std::vector<double> packing_hessian;
    std::string parameter_fingerprint;
};

class ProviderContext;

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

[[nodiscard]] ChemicalSolveResult solve_manufactured_ideal_reaction(
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

[[nodiscard]] ProviderPhaseBlockEvidence evaluate_provider_phase_block(
    const ProviderContext& provider,
    double temperature_k,
    const std::vector<double>& amounts,
    double volume_m3
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
    int max_iterations = 500
);

}  // namespace epcsaft_equilibrium
