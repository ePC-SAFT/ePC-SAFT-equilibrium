#pragma once

#include <array>

#include "provider.hpp"

namespace epcsaft_equilibrium {

struct FlashPhaseEvaluation {
    std::array<double, 2> amounts_mol{};
    double volume_m3 = 0.0;
    MixturePhaseEvaluation provider{};
};

struct FlashNlpEvaluation {
    double objective = 0.0;
    std::array<double, 6> gradient{};
    std::array<double, 2> constraints{};
    std::array<double, 12> jacobian{};
    std::array<double, 21> hessian_lower{};
    FlashPhaseEvaluation liquid{};
    FlashPhaseEvaluation vapor{};
};

[[nodiscard]] FlashNlpEvaluation evaluate_two_phase_flash_nlp(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const std::array<double, 2>& overall_mole_fractions,
    const std::array<double, 6>& variables
);

}  // namespace epcsaft_equilibrium
