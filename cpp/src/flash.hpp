#pragma once

#include <string>
#include <variant>
#include <vector>

#include "held.hpp"
#include "held2_algorithm.hpp"
#include "held2_progress.hpp"
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
    std::variant<HeldResult, Held2AlgorithmResult> solve;
    std::string globality_certificate = "not_guaranteed";
};

[[nodiscard]] FlashResult solve_tp_flash(
    const ProviderContext& provider,
    const FlashInput& input,
    Held2ProgressObserver* observer = nullptr
);

}  // namespace epcsaft_equilibrium
