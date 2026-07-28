#pragma once

#include "held2_step8.hpp"

namespace epcsaft_equilibrium {

enum class Held2Step9Outcome {
    Converged,
    PaperConvergenceFailed,
    Indeterminate,
};

enum class Held2Step9Action {
    TerminateIndeterminate,
    ReturnStageII,
    RunStep10,
};

struct Held2PotentialComparison {
    std::size_t component_index = 0;
    std::uint64_t left_phase_id = 0;
    std::uint64_t right_phase_id = 0;
    double numerator = 0.0;
    double denominator = 0.0;
    double ratio = 0.0;
    bool passed = false;
};

struct Held2PhysicalCertificate {
    double modified_balance_inf = std::numeric_limits<double>::infinity();
    double ordinary_balance_inf = std::numeric_limits<double>::infinity();
    double electroneutrality_inf = std::numeric_limits<double>::infinity();
    double pressure_residual_inf = std::numeric_limits<double>::infinity();
    double kkt_residual_inf = std::numeric_limits<double>::infinity();
    bool accepted = false;
};

struct Held2Step9Result {
    Held2Step9Outcome outcome = Held2Step9Outcome::Indeterminate;
    Held2Step9Action next_action =
        Held2Step9Action::TerminateIndeterminate;
    std::string reason = "not_run";
    std::optional<double> free_energy_gap;
    std::vector<Held2PotentialComparison> potential_comparisons;
    std::optional<Held2PhysicalCertificate> physical;
    Held2StepTiming timing;
};

[[nodiscard]] Held2Step9Result run_held2_step9(
    const Held2Step4Result& step4,
    const Held2Step8Result& step8,
    const Held2StateEvaluator& evaluator,
    Held2ProgressObserver* observer = nullptr
);

}  // namespace epcsaft_equilibrium
