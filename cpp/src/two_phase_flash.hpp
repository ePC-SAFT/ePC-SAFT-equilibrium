#pragma once

#include <array>
#include <string>
#include <vector>

#include "provider.hpp"

namespace epcsaft_equilibrium {

struct FlashPhaseEvaluation {
    std::array<double, 2> amounts_mol{};
    double volume_m3 = 0.0;
    MixturePhaseEvaluation provider{};
};

struct FlashNlpEvaluation {
    double objective = 0.0;
    std::array<double, 6> gradient{};
    std::array<double, 2> constraints{};
    std::array<double, 12> jacobian{};
    std::array<double, 21> hessian_lower{};
    FlashPhaseEvaluation liquid{};
    FlashPhaseEvaluation vapor{};
};

[[nodiscard]] FlashNlpEvaluation evaluate_two_phase_flash_nlp(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const std::array<double, 2>& overall_mole_fractions,
    const std::array<double, 6>& variables
);

struct FlashAttemptRecord {
    std::string role;
    std::array<double, 6> initial_guess{};
    bool solver_converged = false;
    std::string solver_status;
    int iterations = 0;
    double constraint_violation = 0.0;
    std::string callback_error;
};

struct FlashSolveResult {
    bool accepted = false;
    bool solver_converged = false;
    bool numerical_converged = false;
    bool physical_accepted = false;
    std::string failure_reason;
    std::string solver_status;
    int iterations = 0;
    int attempts = 0;
    int confirmation_solves = 0;
    std::array<double, 6> solver_lower_bounds{};
    std::array<double, 6> solver_upper_bounds{};
    std::vector<FlashAttemptRecord> attempt_log;
    double solver_constraint_violation = 0.0;
    double confirmation_max_difference = 0.0;
    double material_balance_max_abs = 0.0;
    double pressure_stationarity_max_relative = 0.0;
    double chemical_potential_max_abs = 0.0;
    double kkt_stationarity_max_abs = 0.0;
    double phase_density_distance = 0.0;
    std::array<double, 2> equality_multipliers{};
    std::array<double, 6> lower_bound_multipliers{};
    std::array<double, 6> upper_bound_multipliers{};
    FlashNlpEvaluation evaluation{};
};

[[nodiscard]] FlashSolveResult solve_two_phase_flash(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const std::array<double, 2>& overall_mole_fractions
);

}  // namespace epcsaft_equilibrium
