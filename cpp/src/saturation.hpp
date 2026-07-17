#pragma once

#include <array>
#include <vector>

#include "provider.hpp"

namespace epcsaft_equilibrium {

struct NlpEvaluation {
    double objective = 0.0;
    std::array<double, 3> objective_gradient{};
    std::array<double, 3> constraints{};
    std::array<double, 9> jacobian{};
    std::array<double, 6> lagrangian_hessian_lower{};
    PhaseEvaluation vapor{};
    PhaseEvaluation liquid{};
    double saturation_pressure_pa = 0.0;
};

[[nodiscard]] NlpEvaluation evaluate_saturation_nlp(
    const ProviderContext& provider,
    double temperature_k,
    const std::array<double, 3>& variables,
    const std::array<double, 3>& multipliers
);

struct SaturationSolveResult {
    bool accepted = false;
    bool solver_converged = false;
    bool numerical_converged = false;
    bool physical_accepted = false;
    std::string failure_reason;
    std::string solver_status;
    int iterations = 0;
    int attempts = 0;
    int confirmation_solves = 0;
    std::array<double, 3> solver_lower_bounds{};
    std::array<double, 3> solver_upper_bounds{};
    struct AttemptRecord {
        std::string role;
        std::array<double, 3> initial_guess{};
        bool solver_converged = false;
        std::string solver_status;
        int iterations = 0;
        double constraint_violation = 0.0;
        std::string callback_error;
    };
    std::vector<AttemptRecord> attempt_log;
    double solver_constraint_violation = 0.0;
    double confirmation_max_relative_difference = 0.0;
    double pressure_relative_residual = 0.0;
    double chemical_potential_absolute_residual = 0.0;
    double phase_density_distance = 0.0;
    NlpEvaluation evaluation{};
};

[[nodiscard]] SaturationSolveResult solve_local_saturation(
    const ProviderContext& provider,
    double temperature_k,
    double liquid_density_upper_mol_m3
);

}  // namespace epcsaft_equilibrium
