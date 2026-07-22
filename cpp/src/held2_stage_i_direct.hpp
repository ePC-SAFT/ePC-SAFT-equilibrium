#pragma once

#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "held2.hpp"
#include "held2_progress.hpp"

namespace epcsaft_equilibrium {

struct Held2StageIReducedEvaluation {
    std::vector<double> chart_coordinates;
    std::vector<double> independent_modified_fractions;
    Held2PressureEnvelopeResult pressure_envelope;
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

[[nodiscard]] Held2StageIDirectResult solve_held2_manufactured_stage_i_direct(
    const std::string& topology,
    int evaluation_budget
);

}  // namespace epcsaft_equilibrium
