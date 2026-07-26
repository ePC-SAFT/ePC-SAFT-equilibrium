#pragma once

#include "held2_step7.hpp"
#include "held2_step10.hpp"

namespace epcsaft_equilibrium {

class ProviderContext;

struct Held2Input {
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    std::vector<double> overall_mole_fractions_provider_order;
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
    Held2StateValueEvaluator evaluate_value;
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
    std::optional<Held2Step10Result> step10;
    std::optional<Held2PersistentState> final_state;
    std::vector<Held2Phase> phases;
    std::optional<double> total_free_energy_over_rt;
    std::vector<Held2StepTiming> step_timings;
    int upper_solve_count = 0;
};

[[nodiscard]] Held2AlgorithmResult run_held2_algorithm(
    const Held2ThermodynamicAccess& thermodynamics,
    const Held2Input& input,
    const Held2ResourceProfile& resources,
    Held2ProgressObserver* observer = nullptr
);

[[nodiscard]] Held2ThermodynamicAccess make_installed_held2_access(
    const ProviderContext& provider,
    const Held2Input& input
);

}  // namespace epcsaft_equilibrium
