#pragma once

#include <string>

#include "flash.hpp"
#include "held2_algorithm.hpp"

namespace epcsaft_equilibrium {

[[nodiscard]] std::string flash_result_to_json(const FlashResult& result);
[[nodiscard]] std::string held2_algorithm_result_to_json(
    const Held2AlgorithmResult& result,
    const Held2Input& input,
    const std::string& parameter_fingerprint
);

}  // namespace epcsaft_equilibrium
