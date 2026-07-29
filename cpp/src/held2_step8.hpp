#pragma once

#include "held2_step6.hpp"

namespace epcsaft_equilibrium {

enum class Held2Step8Outcome {
    CertifiedFeasible,
    CertifiedInfeasible,
    InsufficientCandidates,
    Indeterminate,
};

struct Held2Phase {
    std::uint64_t stable_id = 0;
    double phase_fraction = 0.0;
    std::vector<double> independent_modified_fractions;
    std::vector<double> physical_fractions_provider_order;
    double volume = 0.0;
    double packing_fraction = 0.0;
    double helmholtz_over_rt_reference_amount = 0.0;
    double pressure_pa = 0.0;
    std::vector<double> chemical_potentials_over_rt;
    double reduced_gibbs = 0.0;
    std::vector<double> reduced_gibbs_gradient;
};

struct Held2NlpCertificate {
    std::string solver_status = "not_run";
    double primal_residual_inf = std::numeric_limits<double>::infinity();
    double stationarity_residual_inf =
        std::numeric_limits<double>::infinity();
    double dual_sign_violation_inf =
        std::numeric_limits<double>::infinity();
    double complementarity_inf = std::numeric_limits<double>::infinity();
    bool accepted = false;
};

struct Held2Step8Result {
    Held2Step8Outcome outcome = Held2Step8Outcome::Indeterminate;
    std::string reason = "not_run";
    std::vector<std::uint64_t> candidate_ids;
    std::vector<double> continuation_variables;
    std::optional<double> total_reduced_gibbs;
    double ordinary_balance_inf = std::numeric_limits<double>::infinity();
    double electroneutrality_inf = std::numeric_limits<double>::infinity();
    double electroneutrality_scale = 0.0;
    double pressure_residual_inf = std::numeric_limits<double>::infinity();
    std::vector<Held2Phase> active_phases;
    std::optional<Held2NlpCertificate> nlp;
    Held2StepTiming timing;
};

[[nodiscard]] Held2Step8Result run_held2_step8(
    const Held2Step1Result& step1,
    const Held2Step6Result& step6,
    const Held2StateEvaluator& evaluator,
    const Held2PackingFractionEvaluator& packing_fraction,
    const Held2Step8Result* previous = nullptr,
    const Held2StateValueEvaluator& value_evaluator = {}
);

}  // namespace epcsaft_equilibrium
