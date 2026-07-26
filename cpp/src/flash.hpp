#pragma once

#include <array>
#include <string>
#include <vector>

#include "held.hpp"
#include "held2.hpp"
#include "held2_progress.hpp"
#include "held2_step2.hpp"
#include "provider.hpp"

namespace epcsaft_equilibrium {

struct FlashInput {
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    std::vector<double> overall_mole_fractions;
};

struct FlashResult {
    FlashInput input;
    std::string parameter_fingerprint;
    HeldResult held;
    std::string globality_certificate = "not_guaranteed";
};

struct Held2InstalledPressureEnvelopeDiagnostic {
    Held2PressureEnvelopeResult envelope;
    std::vector<std::string> component_ids;
    std::vector<double> charges;
    std::array<double, 2> molar_volume_bounds{};
    std::string parameter_fingerprint;
};

[[nodiscard]] Held2InstalledPressureEnvelopeDiagnostic
evaluate_held2_installed_pressure_envelope(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& independent_modified_fractions,
    int initial_interval_count
);

[[nodiscard]] Held2StageIIINlpEvaluation
evaluate_held2_installed_stage_iii_derivatives(
    const ProviderContext& provider,
    const FlashInput& input,
    std::size_t phase_count,
    const std::vector<double>& variables,
    const std::vector<double>& equality_multipliers
);

[[nodiscard]] FlashResult solve_tp_flash(
    const ProviderContext& provider,
    const FlashInput& input
);

}  // namespace epcsaft_equilibrium
