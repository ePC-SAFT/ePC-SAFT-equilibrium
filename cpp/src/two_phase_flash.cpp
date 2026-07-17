#include "two_phase_flash.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace epcsaft_equilibrium {
namespace {

constexpr double kGasConstantJPerMolK = 8.31446261815324;

void require_finite(double value, const char* name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
}

std::size_t lower_index(std::size_t row, std::size_t column) {
    return row * (row + 1) / 2 + column;
}

FlashPhaseEvaluation evaluate_phase(
    const ProviderContext& provider,
    double temperature_k,
    const std::array<double, 2>& amounts_mol,
    double log_volume
) {
    const double volume_m3 = std::exp(log_volume);
    require_finite(volume_m3, "phase volume");
    if (amounts_mol[0] <= 0.0 || amounts_mol[1] <= 0.0 || volume_m3 <= 0.0) {
        throw std::invalid_argument("phase amounts and volume must be positive");
    }
    MixturePhaseEvaluation evaluation = provider.evaluate_mixture(
        temperature_k,
        {amounts_mol[0], amounts_mol[1]},
        volume_m3
    );
    return {amounts_mol, volume_m3, std::move(evaluation)};
}

void add_phase_derivatives(
    const FlashPhaseEvaluation& phase,
    double pressure_over_rt,
    std::size_t offset,
    std::array<double, 6>& gradient,
    std::array<double, 21>& hessian_lower
) {
    const auto& provider_gradient = phase.provider.gradient;
    const auto& provider_hessian = phase.provider.hessian;
    const double volume = phase.volume_m3;
    gradient[offset] = provider_gradient[0];
    gradient[offset + 1] = provider_gradient[1];
    gradient[offset + 2] = (provider_gradient[2] + pressure_over_rt) * volume;

    for (std::size_t row = 0; row < 2; ++row) {
        for (std::size_t column = 0; column <= row; ++column) {
            hessian_lower[lower_index(offset + row, offset + column)] =
                provider_hessian[row * 3 + column];
        }
        hessian_lower[lower_index(offset + 2, offset + row)] =
            provider_hessian[2 * 3 + row] * volume;
    }
    hessian_lower[lower_index(offset + 2, offset + 2)] =
        provider_hessian[8] * volume * volume
        + (provider_gradient[2] + pressure_over_rt) * volume;
}

}  // namespace

FlashNlpEvaluation evaluate_two_phase_flash_nlp(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const std::array<double, 2>& overall_mole_fractions,
    const std::array<double, 6>& variables
) {
    require_finite(temperature_k, "temperature");
    require_finite(pressure_pa, "pressure");
    if (temperature_k <= 0.0 || pressure_pa <= 0.0) {
        throw std::invalid_argument("temperature and pressure must be positive");
    }
    for (double value : overall_mole_fractions) {
        require_finite(value, "overall mole fraction");
    }
    for (double value : variables) {
        require_finite(value, "flash NLP variable");
    }

    FlashPhaseEvaluation liquid = evaluate_phase(
        provider,
        temperature_k,
        {variables[0], variables[1]},
        variables[2]
    );
    FlashPhaseEvaluation vapor = evaluate_phase(
        provider,
        temperature_k,
        {variables[3], variables[4]},
        variables[5]
    );
    const double pressure_over_rt = pressure_pa / (kGasConstantJPerMolK * temperature_k);

    FlashNlpEvaluation result;
    result.objective = liquid.provider.value + vapor.provider.value
        + pressure_over_rt * (liquid.volume_m3 + vapor.volume_m3);
    result.constraints = {
        variables[0] + variables[3] - overall_mole_fractions[0],
        variables[1] + variables[4] - overall_mole_fractions[1],
    };
    result.jacobian = {
        1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
    };
    add_phase_derivatives(liquid, pressure_over_rt, 0, result.gradient, result.hessian_lower);
    add_phase_derivatives(vapor, pressure_over_rt, 3, result.gradient, result.hessian_lower);
    result.liquid = std::move(liquid);
    result.vapor = std::move(vapor);
    return result;
}

}  // namespace epcsaft_equilibrium
