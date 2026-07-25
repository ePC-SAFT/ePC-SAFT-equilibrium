#pragma once

#include "held2_step2.hpp"

#include <array>
#include <cstdint>
#include <functional>
#include <optional>

namespace epcsaft_equilibrium {

struct Held2MPoint {
    std::uint64_t insertion_id = 0;
    std::vector<double> independent_modified_fractions;
    double volume = 0.0;
    double packing_fraction = std::numeric_limits<double>::quiet_NaN();
    double reduced_gibbs = 0.0;
    std::vector<double> reduced_gibbs_gradient;
    std::string origin;
};

struct Held2PersistentState {
    int major_iteration = 0;
    int upper_solve_count = 0;
    double feed_reduced_gibbs = 0.0;
    double upper_bound = 0.0;
    double lower_value = -std::numeric_limits<double>::infinity();
    std::uint64_t next_start_ordinal = 0;
    std::optional<std::array<double, 2>> step5_volume_bounds;
    std::vector<double> feed;
    std::vector<double> multipliers;
    std::vector<Held2MPoint> M;
};

using Held2PressureRootEvaluator = std::function<Held2PressureEnvelopeResult(
    const std::vector<double>&
)>;

struct Held2Step3Result {
    std::string status = "indeterminate";
    std::string reason = "not_run";
    std::optional<Held2PersistentState> state;
    Held2StepTiming timing;
};

[[nodiscard]] Held2Step3Result run_held2_step3(
    const Held2Step1Result& step1,
    const Held2Step2Result& step2,
    const Held2PressureRootEvaluator& pressure_roots,
    Held2ProgressObserver* observer = nullptr
);

}  // namespace epcsaft_equilibrium
