#pragma once

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

    [[nodiscard]] const std::string& fingerprint() const;

private:
    const epcsaft_native_sdk_v1& sdk_;
    std::string fingerprint_;
};

}  // namespace epcsaft_equilibrium
