#pragma once

#include <string>
#include <variant>
#include <vector>

#include "held.hpp"
#include "held2_algorithm.hpp"
#include "held2_progress.hpp"
#include "provider.hpp"

namespace epcsaft_equilibrium {

struct FlashResult {
    FlashInput input;
    std::string parameter_fingerprint;
    std::variant<HeldResult, Held2AlgorithmResult> solve;
};

[[nodiscard]] FlashResult solve_tp_flash(
    const ProviderContext& provider,
    const FlashInput& input,
    Held2ProgressObserver* observer = nullptr
);

}  // namespace epcsaft_equilibrium
