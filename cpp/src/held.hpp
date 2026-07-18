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

[[nodiscard]] HeldOuterResult solve_held_outer_envelope(
    const std::vector<HeldOuterCut>& cuts
);

[[nodiscard]] HeldStageIIInitialization initialize_held_stage_ii_cuts(
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
