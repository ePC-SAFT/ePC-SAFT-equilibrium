#pragma once

#include <functional>
#include <limits>
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
    int search_budget,
    Held2ProgressObserver* observer = nullptr
);

struct Held2StageIReducedEvaluation {
    std::vector<double> chart_coordinates;
    std::vector<double> independent_modified_fractions;
    Held2PressureEnvelopeResult pressure_envelope;
    double trial_volume = 0.0;
    double trial_pressure_residual = 0.0;
    double physical_total_ion_mole_fraction =
        std::numeric_limits<double>::quiet_NaN();
    double total_ion_mole_fraction_max =
        std::numeric_limits<double>::quiet_NaN();
    double tpd = std::numeric_limits<double>::infinity();
    bool certified = false;
    std::string failure_reason;
};

using Held2StageIReducedEvaluator = std::function<Held2StageIReducedEvaluation(
    const std::vector<double>&
)>;

struct Held2StageIDirectResult {
    std::string outcome = "indeterminate";
    std::string termination_reason;
    std::string search_strategy = "nlopt_direct_l_pressure_envelope_v1";
    std::string search_solver = "nlopt_gn_direct_l";
    std::string solver_version;
    std::string globality_certificate = "not_guaranteed";
    int declared_evaluation_budget = 0;
    int completed_evaluation_count = 0;
    int failed_evaluation_count = 0;
    int negative_witness_index = -1;
    double total_ion_mole_fraction_max =
        std::numeric_limits<double>::quiet_NaN();
    double minimum_tpd = std::numeric_limits<double>::infinity();
    std::vector<Held2StageIReducedEvaluation> evaluations;
};

[[nodiscard]] Held2StageIDirectResult solve_held2_stage_i_direct(
    std::size_t composition_dimension,
    int evaluation_budget,
    double negative_tpd_threshold,
    const Held2StageIReducedEvaluator& evaluator,
    Held2ProgressObserver* observer = nullptr
);

}  // namespace epcsaft_equilibrium
