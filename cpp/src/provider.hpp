#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include <epcsaft/native_sdk_v1.h>

namespace epcsaft_equilibrium {

struct PhaseEvaluation {
    double amount_mol = 0.0;
    double volume_m3 = 0.0;
    epcsaft_phase_block_result_v1 provider{};
    std::string parameter_fingerprint;
};

struct MixturePhaseEvaluation {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    double pressure_pa = 0.0;
    std::string parameter_fingerprint;
};

struct PackingFractionEvaluation {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
};

struct NeutralReferenceEvaluation {
    std::size_t component_count = 0;
    std::size_t neutral_basis_row_count = 0;
    std::vector<double> neutral_basis;
    std::vector<double> log_fugacity_contractions;
    std::vector<double> pressure_derivatives_per_pa;
    std::vector<double> parameter_derivatives;
    std::size_t active_parameter_count = 0;
    std::uint32_t derivative_availability = 0;
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    double source_pressure_min_pa = 0.0;
    double source_pressure_max_pa = 0.0;
    double maximum_root_residual_pa = 0.0;
    double minimum_pressure_density_derivative_pa_m3_per_mol = 0.0;
    double maximum_density_condition_number = 0.0;
    double reference_derivative_convergence_error = 0.0;
    double maximum_relative_root_bracket_width = 0.0;
    double maximum_relative_root_density_step = 0.0;
    std::size_t stable_root_count = 0;
    std::size_t selected_stable_root_index = 0;
    double reference_convergence_error = 0.0;
    std::string reference_branch;
    std::string basis_id;
    std::string parameter_fingerprint;
    std::string topology_fingerprint;
};

struct ProviderActiveParameterInput {
    std::string family;
    std::string identity;
    std::vector<std::string> component_ids;
    double value = 0.0;
    std::string unit;
};

struct ProviderActiveParameter {
    std::string family;
    std::string identity;
    std::vector<std::string> component_ids;
    double value = 0.0;
    std::string unit;
    std::uint32_t family_code = 0;
    std::uint32_t identity_code = 0;
    std::int32_t component_index = -1;
    std::int32_t pair_component_index_a = -1;
    std::int32_t pair_component_index_b = -1;
};

struct ProviderActiveParameterSet {
    std::vector<ProviderActiveParameter> parameters;
    std::string topology_fingerprint;
};

struct ParameterizedPhaseEvaluation {
    MixturePhaseEvaluation phase;
    PackingFractionEvaluation packing;
    std::array<double, 2> molar_volume_bounds_m3_per_mol{};
    std::vector<double> state_parameter_derivatives;
    std::vector<double> pressure_parameter_derivatives_pa;
    std::vector<double> chemical_potential_parameter_derivatives_over_rt;
    std::string topology_fingerprint;
};

struct InversePackingGeometryEvaluation {
    double volume_m3 = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    std::string topology_fingerprint;
};

class ProviderContext {
public:
    ProviderContext(const epcsaft_native_sdk_v1& sdk, std::string fingerprint);

    [[nodiscard]] PhaseEvaluation evaluate(
        double temperature_k,
        double amount_mol,
        double volume_m3
    ) const;

    [[nodiscard]] MixturePhaseEvaluation evaluate_mixture(
        double temperature_k,
        const std::vector<double>& amounts_mol,
        double volume_m3
    ) const;

    [[nodiscard]] MixturePhaseEvaluation evaluate_electrolyte(
        double temperature_k,
        const std::vector<double>& amounts_mol,
        double volume_m3
    ) const;

    [[nodiscard]] double evaluate_electrolyte_value(
        double temperature_k,
        const std::vector<double>& amounts_mol,
        double volume_m3
    ) const;

    [[nodiscard]] std::array<double, 2> evaluate_molar_volume_bounds(
        double temperature_k,
        const std::vector<double>& mole_fractions,
        double packing_fraction_min,
        double packing_fraction_max
    ) const;

    [[nodiscard]] PackingFractionEvaluation evaluate_packing_fraction(
        double temperature_k,
        const std::vector<double>& amounts_mol,
        double volume_m3
    ) const;

    [[nodiscard]] NeutralReferenceEvaluation evaluate_neutral_reference(
        double temperature_k,
        double pressure_pa
    ) const;

    [[nodiscard]] NeutralReferenceEvaluation evaluate_neutral_reference_derivatives(
        double temperature_k,
        double pressure_pa,
        const ProviderActiveParameterSet& active_parameters = {}
    ) const;

    [[nodiscard]] ProviderActiveParameterSet resolve_active_parameters(
        double temperature_k,
        const std::vector<ProviderActiveParameterInput>& requests
    ) const;

    [[nodiscard]] ParameterizedPhaseEvaluation evaluate_reacting_phase_parameters(
        double temperature_k,
        const std::vector<double>& amounts_mol,
        double volume_m3,
        double packing_fraction_min,
        double packing_fraction_max,
        const ProviderActiveParameterSet& active_parameters
    ) const;

    [[nodiscard]] InversePackingGeometryEvaluation
    evaluate_inverse_packing_geometry(
        double temperature_k,
        const std::vector<double>& amounts_mol,
        double log_packing_fraction,
        const ProviderActiveParameterSet& active_parameters
    ) const;

    [[nodiscard]] const std::string& fingerprint() const;
    [[nodiscard]] const epcsaft_native_sdk_v1& sdk() const noexcept;

private:
    const epcsaft_native_sdk_v1& sdk_;
    std::string fingerprint_;
};

}  // namespace epcsaft_equilibrium
