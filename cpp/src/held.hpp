#pragma once

#include <array>
#include <string>
#include <vector>

#include "provider.hpp"

namespace epcsaft_equilibrium {

inline constexpr double kHeldTunnelPoleStrength = 1.0e-3;

[[nodiscard]] double held_tunnel_minimum_finite_exponential_distance();

struct HeldStateEvaluation {
    double x_methane = 0.0;
    double log_volume = 0.0;
    double volume_m3 = 0.0;
    std::array<double, 2> amounts_mol{};
    double g_bar = 0.0;
    std::array<double, 2> gradient{};
    std::array<double, 4> hessian{};
    double pressure_stationarity_relative = 0.0;
    MixturePhaseEvaluation provider{};
};

struct HeldTpdEvaluation {
    double d_bar = 0.0;
    std::array<double, 2> gradient{};
    std::array<double, 4> hessian{};
    HeldStateEvaluation state{};
};

struct HeldTunnelingEvaluation {
    double objective = 0.0;
    std::array<double, 2> gradient{};
    std::array<double, 4> hessian{};
    HeldTpdEvaluation tpd{};
};

[[nodiscard]] HeldStateEvaluation evaluate_held_state(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    double x_methane,
    double log_volume
);

[[nodiscard]] HeldTpdEvaluation evaluate_held_tpd(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const HeldStateEvaluation& reference,
    double x_methane,
    double log_volume
);

[[nodiscard]] HeldTunnelingEvaluation evaluate_held_tunneling(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const HeldStateEvaluation& reference,
    double minimum_x_methane,
    double minimum_tpd,
    double x_methane,
    double log_volume
);

struct HeldStageIStart {
    std::string role;
    double x_methane = 0.0;
    double log_volume = 0.0;
};

struct HeldReferenceAttempt {
    std::string role;
    double initial_log_volume = 0.0;
    bool solver_converged = false;
    std::string solver_status;
    int iterations = 0;
    bool accepted = false;
    double g_bar = 0.0;
    double log_volume = 0.0;
    double pressure_stationarity_relative = 0.0;
    std::string callback_error;
};

struct HeldStageIAttempt {
    std::string kind;
    std::string role;
    std::array<double, 2> initial_guess{};
    bool solver_converged = false;
    std::string solver_status;
    int iterations = 0;
    bool accepted = false;
    bool materially_perturbed = false;
    double objective = 0.0;
    double tpd = 0.0;
    double pressure_stationarity_relative = 0.0;
    std::string callback_error;
};

struct HeldStageIResult {
    std::string outcome = "indeterminate";
    std::string search_status = "not_started";
    std::string search_profile = "held-stage-i-binary-v1";
    std::string failure_reason;
    std::vector<HeldReferenceAttempt> reference_attempts;
    std::vector<HeldStageIStart> planned_starts;
    std::vector<HeldStageIAttempt> attempt_log;
    int starts_completed = 0;
    int negative_confirmations = 0;
    double confirmation_max_difference = 0.0;
    double best_tpd = 0.0;
    bool has_reference = false;
    HeldStateEvaluation reference{};
    HeldStateEvaluation best_state{};
};

struct HeldOuterCut {
    std::string identity;
    double intercept = 0.0;
    double slope = 0.0;
};

struct HeldOuterResult {
    std::string status = "infeasible";
    std::string failure_reason;
    double value = 0.0;
    double multiplier = 0.0;
    std::vector<std::string> active_cut_ids;
    std::vector<double> tied_multipliers;
};

struct HeldLowerEvaluation {
    double objective = 0.0;
    std::array<double, 2> gradient{};
    std::array<double, 4> hessian{};
    HeldStateEvaluation state{};
};

[[nodiscard]] HeldLowerEvaluation evaluate_held_lower(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane,
    double multiplier,
    double x_methane,
    double log_volume
);

struct HeldStageIIAttempt {
    std::string role;
    std::array<double, 2> initial_guess{};
    bool solver_converged = false;
    std::string solver_status;
    int iterations = 0;
    bool accepted = false;
    double objective = 0.0;
    double pressure_stationarity_relative = 0.0;
    std::string callback_error;
};

struct HeldStageIIStart {
    std::string start_class;
    std::string role;
    double x_methane = 0.0;
    double log_volume = 0.0;
};

struct HeldStageIICut {
    std::string identity;
    std::string endpoint;
    HeldStateEvaluation state{};
    HeldOuterCut outer{};
};

struct HeldStageIIInitialization {
    std::string status = "indeterminate";
    std::string failure_reason;
    std::vector<HeldStageIIAttempt> attempts;
    std::vector<HeldStageIICut> cuts;
    HeldOuterResult outer{};
};

struct HeldStageIILowerSearch {
    std::string status = "search_exhausted";
    double upper_bound = 0.0;
    double multiplier = 0.0;
    std::vector<HeldStageIIStart> planned_starts;
    std::vector<HeldStageIIAttempt> attempts;
    std::vector<HeldStageIICut> accepted_cuts;
    int starts_completed = 0;
};

struct HeldCandidateInput {
    std::string identity;
    double g_bar = 0.0;
    double x_methane = 0.0;
    double volume_m3 = 0.0;
    double composition_gradient = 0.0;
};

struct HeldStageIICandidate {
    std::string identity;
    double objective = 0.0;
    double x_methane = 0.0;
    double volume_m3 = 0.0;
};

struct HeldStageIIRejection {
    std::string identity;
    std::string reason;
};

struct HeldStageIICandidateResult {
    std::string status = "insufficient_candidates";
    std::vector<HeldStageIICandidate> candidates;
    std::vector<HeldStageIIRejection> rejections;
};

struct HeldStageIITrace {
    int major_iteration = 0;
    double outer_value = 0.0;
    double upper_bound = 0.0;
    double multiplier = 0.0;
    std::vector<std::string> active_cut_ids;
    std::vector<std::string> accepted_cut_ids;
    int lower_starts_completed = 0;
    std::vector<std::string> candidate_ids;
    std::vector<HeldStageIIRejection> rejections;
};

struct HeldStageIIResult {
    std::string outcome = "indeterminate";
    std::string search_status = "not_started";
    std::string search_profile = "held-stage-ii-binary-v1";
    std::string failure_reason;
    std::string stage_i_outcome;
    std::string stage_i_search_status;
    double best_tpd = 0.0;
    int major_iterations = 0;
    double upper_bound = 0.0;
    std::vector<HeldStageIIAttempt> endpoint_attempts;
    std::vector<HeldStageIICut> cuts;
    std::vector<HeldStageIICut> candidates;
    std::vector<HeldStageIITrace> trace;
};

[[nodiscard]] HeldOuterResult solve_held_outer_envelope(
    const std::vector<HeldOuterCut>& cuts
);

[[nodiscard]] HeldStageIIInitialization initialize_held_stage_ii_cuts(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane
);

[[nodiscard]] HeldStageIILowerSearch search_held_stage_ii_lower(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane,
    double multiplier,
    double upper_bound,
    const std::vector<HeldStateEvaluation>& previous_states,
    const std::string& identity_prefix = "lower"
);

[[nodiscard]] HeldStageIICandidateResult select_held_stage_ii_candidates(
    double feed_x_methane,
    double upper_bound,
    double multiplier,
    const std::vector<HeldCandidateInput>& points
);

[[nodiscard]] std::string held_stage_ii_budget_status(
    int major_iterations,
    int lower_starts,
    bool lower_satisfied
);

[[nodiscard]] HeldStageIIResult solve_held_stage_ii(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane
);

[[nodiscard]] HeldStageIResult solve_held_stage_i(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane
);

}  // namespace epcsaft_equilibrium
