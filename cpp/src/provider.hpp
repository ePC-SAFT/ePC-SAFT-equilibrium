#pragma once

#include <array>
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
    std::vector<double> reference_composition;
    std::uint32_t derivative_availability = EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_NONE_V1;
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    double solvent_molar_mass_kg_per_mol = 0.0;
    double reference_amount_mol = 0.0;
    double reference_number_density_mol_per_m3 = 0.0;
    double reference_molality_mol_per_kg = 0.0;
    double reference_convergence_error = 0.0;
    std::string basis_id;
    std::string parameter_fingerprint;
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

    [[nodiscard]] const std::string& fingerprint() const;

private:
    const epcsaft_native_sdk_v1& sdk_;
    std::string fingerprint_;
};

}  // namespace epcsaft_equilibrium
