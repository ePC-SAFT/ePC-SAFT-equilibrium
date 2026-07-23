#pragma once

#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "held2.hpp"

namespace epcsaft_equilibrium {

struct Held2StageIIBasinSeed {
    std::vector<double> independent_modified_fractions;
    std::string source;
};

struct Held2StageIIBasinEvaluation {
    std::vector<double> independent_modified_fractions;
    Held2PressureEnvelopeResult pressure_envelope;
    double reduced_lower_value = std::numeric_limits<double>::infinity();
    bool certified = false;
    std::string failure_reason;
};

using Held2StageIIBasinEvaluator = std::function<Held2StageIIBasinEvaluation(
    const std::vector<double>&
)>;

struct Held2StageIIPhysicalStart {
    std::vector<double> independent_modified_fractions;
    int stable_branch_index = -1;
    double log_volume = 0.0;
    double volume = 0.0;
    double reduced_lower_value = 0.0;
    std::string source;
    std::string root_origin;
    std::string root_completeness = "not_proven";
};

struct Held2StageIIBasinExplorationResult {
    std::string outcome = "indeterminate";
    std::string termination_reason;
    std::string strategy = "continuation_sobol_direct_l_v1";
    std::string direct_solver = "nlopt_gn_direct_l";
    std::string direct_solver_version;
    std::string globality_certificate = "not_guaranteed";
    int declared_sobol_count = 0;
    int declared_direct_budget = 0;
    int completed_evaluation_count = 0;
    int failed_evaluation_count = 0;
    int duplicate_start_count = 0;
    bool direct_escalation_used = false;
    std::vector<Held2StageIIBasinEvaluation> evaluations;
    std::vector<Held2StageIIPhysicalStart> representatives;
};

[[nodiscard]] std::vector<std::vector<double>> held2_sobol_points(
    std::size_t dimension,
    int count
);

[[nodiscard]] Held2StageIIBasinExplorationResult
explore_held2_stage_ii_basins(
    const Held2Coordinates& coordinates,
    const std::vector<Held2StageIIBasinSeed>& seeds,
    int sobol_count,
    bool use_direct_escalation,
    int direct_evaluation_budget,
    double total_ion_mole_fraction_max,
    const Held2StageIIBasinEvaluator& evaluator
);

}  // namespace epcsaft_equilibrium
