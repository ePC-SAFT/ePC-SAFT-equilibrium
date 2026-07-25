#pragma once

#include <functional>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "held2_step4.hpp"

namespace epcsaft_equilibrium {

using Held2PackingFractionEvaluator = std::function<double(
    const std::vector<double>&,
    double
)>;

struct Held2ResourceProfile {
    int step2_search_budget = 0;
    int step5_start_epoch_size = 0;
    int step5_total_start_cap = 0;
    int step7_major_iteration_cap = 0;
};

struct Held2LocalCertificate {
    std::uint64_t start_ordinal = 0;
    std::string solver_status;
    bool finite_and_in_domain = false;
    double pressure_residual = std::numeric_limits<double>::infinity();
    double primal_residual_inf = std::numeric_limits<double>::infinity();
    double stationarity_residual_inf = std::numeric_limits<double>::infinity();
    double dual_sign_violation_inf = std::numeric_limits<double>::infinity();
    double complementarity_inf = std::numeric_limits<double>::infinity();
    double dual_pullback_residual_inf =
        std::numeric_limits<double>::infinity();
    bool accepted = false;
};

struct Held2Step5Result {
    std::string status = "indeterminate";
    std::string reason = "not_run";
    std::optional<double> lower_value;
    std::optional<Held2MPoint> terminal;
    bool mathematical_set_changed = false;
    std::uint64_t starts_consumed = 0;
    std::vector<Held2LocalCertificate> attempts;
    Held2StepTiming timing;
};

[[nodiscard]] Held2Step5Result run_held2_step5(
    const Held2Step1Result& step1,
    const Held2Step4Result& step4,
    Held2PersistentState& state,
    const Held2StateEvaluator& evaluator,
    const Held2PackingFractionEvaluator& packing_fraction,
    const Held2ResourceProfile& resources,
    Held2ProgressObserver* observer = nullptr
);

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
