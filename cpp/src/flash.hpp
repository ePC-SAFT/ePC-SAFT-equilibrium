#pragma once

#include <array>
#include <optional>
#include <string>
#include <vector>

#include "held.hpp"
#include "held2.hpp"
#include "held2_controller.hpp"
#include "held2_progress.hpp"
#include "held2_stage_i_direct.hpp"
#include "provider.hpp"

namespace epcsaft_equilibrium {

struct FlashInput {
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    std::vector<double> overall_mole_fractions;
};

struct SolverResourceProfile {
    int stage_i_evaluation_budget = 50;
    int stage_ii_major_iteration_cap = 24;
    int stage_ii_local_attempt_cap_per_major = 50;
};

struct Held2FlashResult {
    Held2PressureEnvelopeResult reference_pressure_envelope;
    Held2StageIDirectResult stage_i;
    std::optional<Held2StageIIResult> stage_ii;
    std::optional<Held2StageIIIResult> stage_iii;
    Held2WorkflowState workflow;
    std::string stage_ii_skip_reason;
    std::string stage_iii_skip_reason;
    std::string predictive_comparison_status =
        "not_allowed_before_physical_acceptance";
};

struct FlashResult {
    enum class Route { Held, Held2 };

    Route route = Route::Held;
    FlashInput input;
    std::string parameter_fingerprint;
    HeldResult held;
    Held2FlashResult held2;
    std::string globality_certificate = "not_guaranteed";
};

struct Held2InstalledPressureEnvelopeDiagnostic {
    Held2PressureEnvelopeResult envelope;
    std::vector<std::string> component_ids;
    std::vector<double> charges;
    std::array<double, 2> molar_volume_bounds{};
    std::string parameter_fingerprint;
};

struct Held2ManufacturedWorkflowResult {
    Held2StageIDirectResult stage_i;
    Held2StageIIResult stage_ii;
    Held2StageIIIResult stage_iii;
    Held2WorkflowState workflow;
    std::string globality_certificate = "not_guaranteed";
};

[[nodiscard]] Held2ManufacturedWorkflowResult
solve_held2_manufactured_workflow(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed
);

[[nodiscard]] Held2InstalledPressureEnvelopeDiagnostic
evaluate_held2_installed_pressure_envelope(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& independent_modified_fractions,
    int initial_interval_count
);

[[nodiscard]] Held2StageIIINlpEvaluation
evaluate_held2_installed_stage_iii_derivatives(
    const ProviderContext& provider,
    const FlashInput& input,
    std::size_t phase_count,
    const std::vector<double>& variables,
    const std::vector<double>& equality_multipliers
);

[[nodiscard]] FlashResult solve_tp_flash(
    const ProviderContext& provider,
    const FlashInput& input,
    const SolverResourceProfile& resources = {},
    Held2ProgressObserver* observer = nullptr
);

}  // namespace epcsaft_equilibrium
