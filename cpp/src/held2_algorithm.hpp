#pragma once

#include "held2_step7.hpp"
#include "held2_step10.hpp"

#include <cstddef>
#include <string_view>

namespace epcsaft_equilibrium {

struct FlashInput {
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    std::vector<double> overall_mole_fractions;
};

struct Held2ThermodynamicAccess {
    std::vector<std::string> component_ids;
    std::vector<double> charges;
    Held2StateEvaluator evaluate;
    Held2PhysicalVolumeBoundsEvaluator volume_bounds_physical;
    Held2PackingFractionEvaluator packing_fraction;
    double total_ion_mole_fraction_max =
        std::numeric_limits<double>::quiet_NaN();
    Held2StateEvaluator evaluate_trace;
};

struct Held2AlgorithmResult {
    std::string outcome = "indeterminate";
    std::string failure_stage;
    std::string failure_reason;
    std::string globality_certificate = "not_guaranteed";
    std::string phase_enumeration_certificate =
        "completeness_not_guaranteed";
    Held2Step1Result step1;
    std::optional<Held2Step2Result> step2;
    std::optional<Held2Step3Result> step3;
    std::vector<Held2Step4Result> step4_history;
    std::vector<Held2Step5Result> step5_history;
    std::vector<Held2Step6Result> step6_history;
    std::vector<Held2Step7Result> step7_history;
    std::vector<Held2Step8Result> step8_history;
    std::vector<Held2Step9Result> step9_history;
    std::vector<Held2Step10Result> step10_history;
    std::optional<Held2Step10Result> step10;
    std::optional<Held2PersistentState> final_state;
    std::vector<Held2Phase> phases;
    std::optional<double> total_free_energy_over_rt;
    std::vector<Held2StepTiming> step_timings;
    int upper_solve_count = 0;
};

[[nodiscard]] bool held2_step8_failure_returns_to_stage_ii(
    std::string_view reason
);

[[nodiscard]] bool held2_insufficient_candidate_recovery_required(
    std::string_view step5_reason,
    std::size_t candidate_count
);

[[nodiscard]] bool retain_held2_step8_feedback(
    Held2PersistentState& state,
    const Held2Step8Result& step8
);

[[nodiscard]] Held2AlgorithmResult run_held2_algorithm(
    const Held2ThermodynamicAccess& thermodynamics,
    const FlashInput& input,
    const Held2ResourceProfile& resources,
    Held2ProgressObserver* observer = nullptr
);

}  // namespace epcsaft_equilibrium
