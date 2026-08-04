#pragma once

#include "held2.hpp"

#include <optional>

namespace epcsaft_equilibrium {

struct Held2Step1Result {
    std::string status = "indeterminate";
    std::string reason = "not_run";
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    double total_ion_mole_fraction_max =
        std::numeric_limits<double>::quiet_NaN();
    std::string provider_fingerprint;
    std::optional<Held2Coordinates> coordinates;
    std::optional<std::vector<double>> independent_feed;
    std::optional<std::array<double, 2>> feed_volume_bounds;
    std::optional<Held2VolumeBoundsEvaluator> volume_bounds;
    Held2StepTiming timing;
};

[[nodiscard]] Held2Step1Result run_held2_step1(
    const std::vector<std::string>& component_ids,
    const std::vector<double>& charges,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& physical_feed,
    const Held2PhysicalVolumeBoundsEvaluator& physical_volume_bounds,
    std::string provider_fingerprint,
    double total_ion_mole_fraction_max =
        std::numeric_limits<double>::quiet_NaN()
);

}  // namespace epcsaft_equilibrium
