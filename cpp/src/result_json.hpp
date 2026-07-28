#pragma once

#include <string>

#include "flash.hpp"

namespace epcsaft_equilibrium {

[[nodiscard]] std::string flash_result_to_json(const FlashResult& result);

}  // namespace epcsaft_equilibrium
