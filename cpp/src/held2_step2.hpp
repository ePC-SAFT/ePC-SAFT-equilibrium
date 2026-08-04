#pragma once

#include <optional>
#include <string>
#include <vector>

#include "held2_progress.hpp"
#include "held2_step1.hpp"

namespace epcsaft_equilibrium {

enum class Held2Step2Outcome {
    NegativeWitness,
    NoNegativeWitnessDetected,
    Indeterminate,
};

struct Held2TpdEvaluation {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
};

struct Held2Step2Result {
    Held2Step2Outcome outcome = Held2Step2Outcome::Indeterminate;
    std::string reason = "not_run";
    std::string globality_certificate = "not_guaranteed";
    std::optional<Held2PressureEnvelopeResult> reference_envelope;
    std::optional<Held2StateEvaluation> reference;
    std::optional<Held2StageICandidate> negative_witness;
    std::optional<double> minimum_tpd;
    Held2StepTiming timing;
};

[[nodiscard]] Held2TpdEvaluation evaluate_held2_tpd(
    const Held2StateEvaluation& reference,
    const std::vector<double>& feed,
    const Held2StateEvaluation& trial,
    const std::vector<double>& independent
);

[[nodiscard]] Held2Step2Result run_held2_step2(
    const Held2Step1Result& step1,
    const Held2StateEvaluator& evaluator,
    int search_work_budget,
    Held2ProgressObserver* observer = nullptr
);

}  // namespace epcsaft_equilibrium
