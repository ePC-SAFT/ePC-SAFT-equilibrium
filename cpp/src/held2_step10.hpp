#pragma once

#include "held2_step9.hpp"

namespace epcsaft_equilibrium {

struct Held2TraceRefinement {
    std::uint64_t phase_id = 0;
    std::size_t component_index = 0;
    std::uint64_t reference_phase_id = 0;
    double initial_mole_fraction = 0.0;
    double refined_mole_fraction = 0.0;
    double final_potential_residual =
        std::numeric_limits<double>::infinity();
    std::string status = "not_run";
};

struct Held2Step10Result {
    std::string status = "indeterminate";
    std::string reason = "not_run";
    std::vector<Held2Phase> phases;
    std::vector<Held2TraceRefinement> refinements;
    std::optional<Held2PhysicalCertificate> final_certificate;
    Held2StepTiming timing;
};

[[nodiscard]] std::vector<double> held2_lift_trace_fractions(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent_modified_fractions
);

[[nodiscard]] Held2Step10Result run_held2_step10(
    const Held2Step1Result& step1,
    const Held2Step8Result& step8,
    const Held2Step9Result& step9,
    const Held2StateEvaluator& evaluator,
    Held2ProgressObserver* observer = nullptr
);

}  // namespace epcsaft_equilibrium
