#include "flash.hpp"

#include "held2_algorithm.hpp"
#include "held2_tolerances.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>

namespace epcsaft_equilibrium {
namespace {

constexpr double kGasConstantJPerMolK = 8.31446261815324;
constexpr const char* kNeutralFlashFingerprint =
    "sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286";
constexpr std::size_t kMixtureSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, evaluate_mixture_phase)
    + sizeof(epcsaft_evaluate_mixture_phase_v1);
constexpr std::size_t kElectrolyteSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, evaluate_electrolyte_phase)
    + sizeof(epcsaft_evaluate_mixture_phase_v1);
constexpr std::size_t kMolarVolumeSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, evaluate_molar_volume_bounds)
    + sizeof(epcsaft_evaluate_molar_volume_bounds_v1);
constexpr std::size_t kSourceDomainSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, total_ion_mole_fraction_max)
    + sizeof(double);

void require_mixture_sdk(const epcsaft_native_sdk_v1& sdk) {
    if (sdk.table_size < kMixtureSdkTableSize) {
        throw std::invalid_argument("provider capsule is missing the mixture SDK tail");
    }
    if (sdk.mixture_result_size != sizeof(epcsaft_mixture_phase_block_result_v1)) {
        throw std::invalid_argument(
            "provider capsule mixture result size does not match the v1 contract"
        );
    }
    if (sdk.evaluate_mixture_phase == nullptr) {
        throw std::invalid_argument("provider capsule is missing its mixture phase evaluator");
    }
}

void require_held2_sdk(const epcsaft_native_sdk_v1& sdk) {
    if (sdk.table_size < kElectrolyteSdkTableSize
        || sdk.component_count < 3
        || sdk.component_ids == nullptr
        || sdk.component_charges == nullptr
        || sdk.evaluate_electrolyte_phase == nullptr) {
        throw std::invalid_argument("provider capsule is missing the electrolyte phase contract");
    }
    if (sdk.mixture_result_size != sizeof(epcsaft_mixture_phase_block_result_v1)) {
        throw std::invalid_argument(
            "provider capsule mixture result size does not match the v1 contract"
        );
    }
    if (sdk.table_size < kMolarVolumeSdkTableSize
        || sdk.evaluate_molar_volume_bounds == nullptr) {
        throw std::invalid_argument(
            "provider capsule is missing the molar-volume domain contract"
        );
    }
    if (sdk.table_size < kSourceDomainSdkTableSize) {
        throw std::invalid_argument(
            "provider capsule is missing the electrolyte source-domain contract"
        );
    }
}

bool charged(const epcsaft_native_sdk_v1& sdk) {
    return sdk.table_size >= kElectrolyteSdkTableSize
        && sdk.component_charges != nullptr
        && std::any_of(
            sdk.component_charges,
            sdk.component_charges + sdk.component_count,
            [](int32_t charge) { return charge != 0; }
        );
}

void validate_input(const ProviderContext& provider, const FlashInput& input) {
    if (!std::isfinite(input.temperature_k) || input.temperature_k <= 0.0
        || !std::isfinite(input.pressure_pa) || input.pressure_pa <= 0.0
        || input.overall_mole_fractions.size() != provider.sdk().component_count
        || !std::all_of(
            input.overall_mole_fractions.begin(),
            input.overall_mole_fractions.end(),
            [](double value) { return std::isfinite(value) && value > 0.0; }
        )) {
        throw std::invalid_argument(
            "tp_flash inputs must be positive, finite, and match the Provider model"
        );
    }
    double total = 0.0;
    for (double value : input.overall_mole_fractions) {
        total += value;
    }
    if (std::abs(total - 1.0) > 1.0e-12) {
        throw std::invalid_argument("tp_flash composition must be normalized");
    }
}

double total_ion_mole_fraction(
    const std::vector<double>& charges,
    const std::vector<double>& physical_fractions
) {
    if (charges.size() != physical_fractions.size()) {
        throw std::invalid_argument(
            "HELD2 charge and physical-composition dimensions changed"
        );
    }
    double total = 0.0;
    for (std::size_t index = 0; index < charges.size(); ++index) {
        if (charges[index] != 0.0) {
            total += physical_fractions[index];
        }
    }
    return total;
}

class InstalledHeld2Problem {
public:
    InstalledHeld2Problem(
        const ProviderContext& provider,
        const FlashInput& input
    )
        : provider_(provider),
          input_(input),
          coordinates_(make_held2_coordinates(charges_from_sdk(provider.sdk()))),
          total_ion_mole_fraction_max_(
              provider.sdk().total_ion_mole_fraction_max
          ),
          pressure_over_rt_(
              input.pressure_pa / (kGasConstantJPerMolK * input.temperature_k)
          ) {
        const std::vector<double> modified_feed =
            held2_transform_physical_fractions(
                coordinates_, input_.overall_mole_fractions
            );
        independent_feed_.reserve(coordinates_.independent_indices.size());
        for (std::size_t component : coordinates_.independent_indices) {
            const auto retained = std::find(
                coordinates_.retained_indices.begin(),
                coordinates_.retained_indices.end(),
                component
            );
            if (retained == coordinates_.retained_indices.end()) {
                throw std::invalid_argument(
                    "HELD2 independent component is not retained"
                );
            }
            independent_feed_.push_back(modified_feed[static_cast<std::size_t>(
                retained - coordinates_.retained_indices.begin()
            )]);
        }
        static_cast<void>(held2_map_unit_cube_to_independent_fractions(
            coordinates_,
            std::vector<double>(coordinates_.independent_indices.size(), 0.0),
            total_ion_mole_fraction_max_
        ));
    }

    [[nodiscard]] Held2StateEvaluation evaluate(
        const std::vector<double>& independent,
        double log_volume
    ) const {
        return evaluate_physical(
            independent,
            held2_lift_independent_fractions(coordinates_, independent),
            log_volume
        );
    }

    [[nodiscard]] Held2StateEvaluation evaluate_trace(
        const std::vector<double>& independent,
        double log_volume
    ) const {
        return evaluate_physical(
            independent,
            held2_lift_trace_fractions(coordinates_, independent),
            log_volume
        );
    }

private:
    [[nodiscard]] Held2StateEvaluation evaluate_physical(
        const std::vector<double>& independent,
        const std::vector<double>& amounts,
        double log_volume
    ) const {
        MixturePhaseEvaluation provider_phase =
            provider_.evaluate_electrolyte(
                input_.temperature_k, amounts, std::exp(log_volume)
            );
        Held2PhysicalPhaseBlock block;
        block.helmholtz_over_rt = provider_phase.value;
        block.gradient = std::move(provider_phase.gradient);
        block.hessian = std::move(provider_phase.hessian);
        block.pressure_pa = provider_phase.pressure_pa;
        return evaluate_held2_phase_block(
            coordinates_,
            independent,
            log_volume,
            pressure_over_rt_,
            input_.pressure_pa,
            block
        );
    }

public:

    [[nodiscard]] std::array<double, 2> volume_bounds(
        const std::vector<double>& independent
    ) const {
        return provider_.evaluate_molar_volume_bounds(
            input_.temperature_k,
            held2_lift_independent_fractions(coordinates_, independent),
            kHeld2PackingFractionMinimum,
            kHeld2PackingFractionMaximum
        );
    }

    [[nodiscard]] Held2PressureEnvelopeResult envelope(
        const std::vector<double>& independent
    ) const {
        const auto evaluator = [this](
            const std::vector<double>& values,
            double log_volume
        ) {
            return evaluate(values, log_volume);
        };
        return evaluate_held2_pressure_envelope(
            independent, volume_bounds(independent), evaluator, 64, 8
        );
    }

    [[nodiscard]] const Held2Coordinates& coordinates() const {
        return coordinates_;
    }

    [[nodiscard]] const std::vector<double>& physical_feed() const {
        return input_.overall_mole_fractions;
    }

    [[nodiscard]] const std::vector<double>& independent_feed() const {
        return independent_feed_;
    }

    [[nodiscard]] std::vector<double> independent(
        const std::vector<double>& physical
    ) const {
        const std::vector<double> modified =
            held2_transform_physical_fractions(coordinates_, physical);
        std::vector<double> result;
        for (std::size_t provider : coordinates_.independent_indices) {
            const auto retained = std::find(
                coordinates_.retained_indices.begin(),
                coordinates_.retained_indices.end(),
                provider
            );
            result.push_back(modified[static_cast<std::size_t>(
                retained - coordinates_.retained_indices.begin()
            )]);
        }
        return result;
    }

    [[nodiscard]] double total_ion_mole_fraction_max() const {
        return total_ion_mole_fraction_max_;
    }

private:
    static std::vector<double> charges_from_sdk(
        const epcsaft_native_sdk_v1& sdk
    ) {
        if (sdk.component_charges == nullptr) {
            throw std::invalid_argument(
                "Provider SDK is missing electrolyte component charges"
            );
        }
        std::vector<double> charges;
        charges.reserve(sdk.component_count);
        for (std::size_t component = 0; component < sdk.component_count; ++component) {
            charges.push_back(static_cast<double>(sdk.component_charges[component]));
        }
        return charges;
    }

    const ProviderContext& provider_;
    const FlashInput& input_;
    Held2Coordinates coordinates_;
    double total_ion_mole_fraction_max_;
    double pressure_over_rt_;
    std::vector<double> independent_feed_;
};

struct InstalledStageI {
    Held2PressureEnvelopeResult reference_envelope;
    Held2StateEvaluation reference;
    Held2StageIDirectResult search;
};

InstalledStageI run_stage_i(
    const InstalledHeld2Problem& problem,
    int evaluation_budget,
    Held2ProgressObserver* observer
) {
    InstalledStageI result;
    result.search.total_ion_mole_fraction_max =
        problem.total_ion_mole_fraction_max();
    Held2ProgressEvent progress;
    progress.kind = Held2ProgressKind::ReferenceStart;
    observe_held2(observer, progress);
    try {
        result.reference_envelope = problem.envelope(problem.independent_feed());
    } catch (const std::exception& error) {
        result.search.declared_evaluation_budget = evaluation_budget;
        result.search.termination_reason =
            std::string("reference_envelope_failed: ") + error.what();
        progress = {};
        progress.kind = Held2ProgressKind::Failure;
        progress.stage = "REFERENCE";
        progress.reason = result.search.termination_reason;
        observe_held2(observer, progress);
        return result;
    }
    for (std::size_t index = 0; index < result.reference_envelope.roots.size();
         ++index) {
        const auto& root = result.reference_envelope.roots[index];
        progress = {};
        progress.kind = Held2ProgressKind::ReferenceRoot;
        progress.count = static_cast<int>(index + 1);
        progress.volume = root.volume;
        progress.pressure_residual = root.pressure_residual;
        progress.objective = root.objective;
        progress.mechanical_class = root.mechanical_class;
        observe_held2(observer, progress);
    }
    if (result.reference_envelope.outcome != "selected") {
        result.search.declared_evaluation_budget = evaluation_budget;
        result.search.termination_reason =
            "reference_envelope_failed: "
            + result.reference_envelope.failure_reason;
        progress = {};
        progress.kind = Held2ProgressKind::Failure;
        progress.stage = "REFERENCE";
        progress.reason = result.search.termination_reason;
        observe_held2(observer, progress);
        return result;
    }
    result.reference = result.reference_envelope.roots[
        static_cast<std::size_t>(
            result.reference_envelope.selected_root_index
        )
    ].state;
    progress = {};
    progress.kind = Held2ProgressKind::ReferenceSelected;
    progress.count = result.reference_envelope.selected_root_index + 1;
    progress.status = "certified";
    observe_held2(observer, progress);
    progress = {};
    progress.kind = Held2ProgressKind::StageStart;
    progress.stage = "STAGE I - DIRECT-L TPD SEARCH";
    observe_held2(observer, progress);
    const Held2StageIReducedEvaluator evaluator = [
        &problem,
        reference = result.reference
    ](const std::vector<double>& chart_coordinates) {
        Held2StageIReducedEvaluation evaluation;
        evaluation.chart_coordinates = chart_coordinates;
        try {
            evaluation.independent_modified_fractions =
                held2_map_unit_cube_to_independent_fractions(
                    problem.coordinates(),
                    chart_coordinates,
                    problem.total_ion_mole_fraction_max()
                );
            const std::vector<double> physical =
                held2_lift_independent_fractions(
                    problem.coordinates(),
                    evaluation.independent_modified_fractions
                );
            evaluation.physical_total_ion_mole_fraction =
                total_ion_mole_fraction(
                    problem.coordinates().charges, physical
                );
            evaluation.total_ion_mole_fraction_max =
                problem.total_ion_mole_fraction_max();
            evaluation.pressure_envelope = problem.envelope(
                evaluation.independent_modified_fractions
            );
            if (evaluation.pressure_envelope.outcome != "selected") {
                evaluation.failure_reason =
                    evaluation.pressure_envelope.failure_reason;
                return evaluation;
            }
            const Held2StateEvaluation& selected =
                evaluation.pressure_envelope.roots[
                    static_cast<std::size_t>(
                        evaluation.pressure_envelope.selected_root_index
                    )
                ].state;
            evaluation.tpd = selected.objective - reference.objective;
            for (std::size_t index = 0;
                 index < problem.independent_feed().size();
                 ++index) {
                evaluation.tpd -= reference.gradient[index]
                    * (evaluation.independent_modified_fractions[index]
                       - problem.independent_feed()[index]);
            }
            evaluation.certified = true;
        } catch (const std::exception& error) {
            evaluation.failure_reason =
                std::string("provider_evaluation_failed: ") + error.what();
        }
        return evaluation;
    };
    result.search = solve_held2_stage_i_direct(
        problem.coordinates().independent_indices.size(),
        evaluation_budget,
        -1.0e-8,
        evaluator,
        observer
    );
    result.search.total_ion_mole_fraction_max =
        problem.total_ion_mole_fraction_max();
    progress = {};
    progress.kind = Held2ProgressKind::Certificate;
    progress.stage = "STAGE I";
    progress.status = result.search.outcome;
    progress.reason = result.search.termination_reason;
    observe_held2(observer, progress);
    return result;
}

void observe_skipped(
    Held2ProgressObserver* observer,
    const std::string& reason,
    bool stage_ii
) {
    Held2ProgressEvent progress;
    progress.kind = Held2ProgressKind::StageSkipped;
    if (stage_ii) {
        progress.stage = "STAGE II";
        progress.reason = reason;
        observe_held2(observer, progress);
    }
    progress.stage = "STAGE III";
    progress.reason = reason;
    observe_held2(observer, progress);
}

void observe_final(
    Held2ProgressObserver* observer,
    const Held2WorkflowState& workflow
) {
    Held2ProgressEvent progress;
    progress.kind = Held2ProgressKind::Final;
    progress.status = workflow.outcome;
    progress.reason = workflow.failure_reason;
    observe_held2(observer, progress);
}

Held2FlashResult solve_held2(
    const ProviderContext& provider,
    const FlashInput& input,
    const SolverResourceProfile& resources,
    Held2ProgressObserver* observer
) {
    const InstalledHeld2Problem problem(provider, input);
    Held2ProgressEvent progress;
    progress.kind = Held2ProgressKind::CaseStart;
    progress.case_id = "installed-held2-controller";
    progress.temperature_k = input.temperature_k;
    progress.pressure_pa = input.pressure_pa;
    observe_held2(observer, progress);

    const InstalledStageI stage_i = run_stage_i(
        problem, resources.stage_i_evaluation_budget, observer
    );
    Held2WorkflowController workflow;
    Held2FlashResult result;
    result.reference_pressure_envelope = stage_i.reference_envelope;
    result.stage_i = stage_i.search;

    Held2NextAction next_action =
        workflow.complete_reference(stage_i.reference_envelope);
    if (next_action == Held2NextAction::TerminateIndeterminate) {
        result.stage_ii_skip_reason = workflow.state().failure_reason;
        result.stage_iii_skip_reason = result.stage_ii_skip_reason;
        observe_skipped(observer, result.stage_ii_skip_reason, true);
        observe_final(observer, workflow.state());
        result.workflow = workflow.state();
        return result;
    }
    next_action = workflow.complete_stage_i(stage_i.search);
    if (next_action == Held2NextAction::TerminateIndeterminate) {
        const std::string reason =
            workflow.state().outcome
                == "stage_i_finite_search_without_negative_witness"
            ? "stage_i_negative_witness_not_found"
            : workflow.state().failure_reason;
        result.stage_ii_skip_reason = reason;
        result.stage_iii_skip_reason = reason;
        observe_skipped(observer, reason, true);
        observe_final(observer, workflow.state());
        result.workflow = workflow.state();
        return result;
    }

    const auto& witness = stage_i.search.evaluations[
        static_cast<std::size_t>(stage_i.search.negative_witness_index)
    ];
    const auto& witness_state = witness.pressure_envelope.roots[
        static_cast<std::size_t>(
            witness.pressure_envelope.selected_root_index
        )
    ].state;
    const std::vector<Held2StageICandidate> witnesses = {{
        witness_state.modified_fractions,
        witness_state.volume,
        witness.tpd,
    }};
    progress = {};
    progress.kind = Held2ProgressKind::StageStart;
    progress.stage = "STAGE II - HiGHS LP / IPOPT LOWER SEARCH";
    observe_held2(observer, progress);
    const Held2StateEvaluator evaluator = [&problem](
        const std::vector<double>& independent,
        double log_volume
    ) {
        return problem.evaluate(independent, log_volume);
    };
    const Held2VolumeBoundsEvaluator bounds_evaluator = [&problem](
        const std::vector<double>& independent
    ) {
        return problem.volume_bounds(independent);
    };
    result.stage_ii = solve_held2_stage_ii(
        problem.coordinates(),
        problem.physical_feed(),
        evaluator,
        bounds_evaluator,
        stage_i.reference,
        witnesses,
        problem.total_ion_mole_fraction_max(),
        resources.stage_ii_major_iteration_cap,
        resources.stage_ii_local_attempt_cap_per_major,
        observer
    );
    next_action = workflow.complete_stage_ii(*result.stage_ii);
    if (next_action != Held2NextAction::EnterStageIII) {
        result.stage_iii_skip_reason = result.stage_ii->outcome;
        progress = {};
        progress.kind = Held2ProgressKind::Failure;
        progress.stage = "STAGE II";
        progress.reason = result.stage_ii->outcome;
        observe_held2(observer, progress);
        observe_skipped(observer, result.stage_iii_skip_reason, false);
        observe_final(observer, workflow.state());
        result.workflow = workflow.state();
        return result;
    }

    progress = {};
    progress.kind = Held2ProgressKind::StageStart;
    progress.stage = "STAGE III - IPOPT TOTAL FREE ENERGY";
    observe_held2(observer, progress);
    std::vector<std::array<double, 2>> phase_coordinate_bounds;
    phase_coordinate_bounds.reserve(result.stage_ii->candidates.size());
    for (const auto& candidate : result.stage_ii->candidates) {
        const std::array<double, 2> bounds =
            problem.volume_bounds(candidate.independent_modified_fractions);
        phase_coordinate_bounds.push_back({
            std::log(bounds[0]), std::log(bounds[1])
        });
    }
    const double free_energy_upper_bound =
        result.stage_ii->bound_history.empty()
        ? std::numeric_limits<double>::quiet_NaN()
        : result.stage_ii->bound_history.back().upper_bound;
    const std::string free_energy_gap_provenance =
        result.stage_ii->bound_history.empty()
        ? "unavailable"
        : "stage_ii_problem_64_same_major_upper_bound";
    result.stage_iii = solve_held2_stage_iii(
        problem.coordinates(),
        problem.physical_feed(),
        result.stage_ii->candidates,
        evaluator,
        phase_coordinate_bounds,
        free_energy_upper_bound,
        free_energy_gap_provenance
    );
    next_action = workflow.complete_stage_iii(*result.stage_iii);
    if (next_action != Held2NextAction::AcceptMultiphase) {
        progress = {};
        progress.kind = Held2ProgressKind::Failure;
        progress.stage = "STAGE III";
        progress.reason = result.stage_iii->failure_reason;
        observe_held2(observer, progress);
        observe_final(observer, workflow.state());
        result.workflow = workflow.state();
        return result;
    }
    result.predictive_comparison_status =
        "eligible_but_not_executed_private_controller";
    progress = {};
    progress.kind = Held2ProgressKind::Final;
    progress.status = "physical_equilibrium_accepted";
    observe_held2(observer, progress);
    result.workflow = workflow.state();
    return result;
}

}  // namespace

Held2ThermodynamicAccess make_installed_held2_access(
    const ProviderContext& provider,
    const Held2Input& input
) {
    require_held2_sdk(provider.sdk());
    const FlashInput flash_input{
        input.temperature_k,
        input.pressure_pa,
        input.overall_mole_fractions_provider_order,
    };
    validate_input(provider, flash_input);
    const auto problem = std::make_shared<InstalledHeld2Problem>(
        provider, flash_input
    );
    Held2ThermodynamicAccess result;
    result.component_ids.reserve(provider.sdk().component_count);
    result.charges.reserve(provider.sdk().component_count);
    for (std::size_t component = 0;
         component < provider.sdk().component_count;
         ++component) {
        const char* id = provider.sdk().component_ids[component];
        if (id == nullptr || id[0] == '\0') {
            throw std::invalid_argument(
                "Provider electrolyte component ID must not be empty"
            );
        }
        result.component_ids.emplace_back(id);
        result.charges.push_back(static_cast<double>(
            provider.sdk().component_charges[component]
        ));
    }
    result.evaluate = [problem](const auto& composition, double log_volume) {
        return problem->evaluate(composition, log_volume);
    };
    result.evaluate_trace = [problem](
        const auto& composition, double log_volume
    ) {
        return problem->evaluate_trace(composition, log_volume);
    };
    result.volume_bounds_physical = [problem](const auto& physical) {
        return problem->volume_bounds(problem->independent(physical));
    };
    result.packing_fraction = [problem](
        const auto& composition, double volume
    ) {
        return kHeld2PackingFractionMinimum
            * problem->volume_bounds(composition)[1] / volume;
    };
    result.total_ion_mole_fraction_max =
        problem->total_ion_mole_fraction_max();
    return result;
}

Held2ManufacturedWorkflowResult solve_held2_manufactured_workflow(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed
) {
    Held2WorkflowController workflow;
    Held2PressureEnvelopeResult reference;
    reference.outcome = "selected";
    reference.selected_root_index = 0;
    reference.roots.resize(1);
    Held2ManufacturedWorkflowResult result;
    result.stage_i = solve_held2_manufactured_stage_i_direct(
        "negative", 80
    );
    result.stage_ii = solve_held2_manufactured_stage_ii(
        charges, physical_feed
    );
    if (
        workflow.complete_reference(reference)
            != Held2NextAction::ContinueStageII
        || workflow.complete_stage_i(result.stage_i)
            != Held2NextAction::ContinueStageII
        || workflow.complete_stage_ii(result.stage_ii)
            != Held2NextAction::EnterStageIII
    ) {
        throw std::logic_error(
            "manufactured HELD2 evidence did not reach Stage III"
        );
    }

    std::vector<std::array<double, 2>> candidates;
    candidates.reserve(result.stage_ii.candidates.size());
    for (const auto& candidate : result.stage_ii.candidates) {
        if (candidate.independent_modified_fractions.size() != 1) {
            throw std::logic_error(
                "manufactured HELD2 controller requires one independent "
                "modified composition"
            );
        }
        candidates.push_back({
            candidate.independent_modified_fractions.front(),
            candidate.volume,
        });
    }
    result.stage_iii = solve_held2_manufactured_stage_iii(
        charges, physical_feed, candidates
    );
    if (
        workflow.complete_stage_iii(result.stage_iii)
        != Held2NextAction::AcceptMultiphase
    ) {
        throw std::logic_error(
            "manufactured HELD2 evidence did not pass final certification"
        );
    }
    result.workflow = workflow.state();
    return result;
}

Held2InstalledPressureEnvelopeDiagnostic
evaluate_held2_installed_pressure_envelope(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& independent_modified_fractions,
    int initial_interval_count
) {
    require_held2_sdk(provider.sdk());
    const Held2Coordinates coordinates = make_held2_coordinates(
        [&provider]() {
            std::vector<double> charges;
            charges.reserve(provider.sdk().component_count);
            for (std::size_t component = 0;
                 component < provider.sdk().component_count;
                 ++component) {
                charges.push_back(static_cast<double>(
                    provider.sdk().component_charges[component]
                ));
            }
            return charges;
        }()
    );
    FlashInput input{
        temperature_k,
        pressure_pa,
        held2_lift_independent_fractions(
            coordinates, independent_modified_fractions
        ),
    };
    validate_input(provider, input);
    const InstalledHeld2Problem problem(provider, input);

    Held2InstalledPressureEnvelopeDiagnostic diagnostic;
    diagnostic.charges = problem.coordinates().charges;
    diagnostic.component_ids.reserve(provider.sdk().component_count);
    for (std::size_t component = 0;
         component < provider.sdk().component_count;
         ++component) {
        const char* component_id = provider.sdk().component_ids[component];
        if (component_id == nullptr || component_id[0] == '\0') {
            throw std::invalid_argument(
                "provider electrolyte component ID must not be empty"
            );
        }
        diagnostic.component_ids.emplace_back(component_id);
    }
    diagnostic.molar_volume_bounds = problem.volume_bounds(
        independent_modified_fractions
    );
    const Held2StateEvaluator evaluator = [&problem](
        const std::vector<double>& independent,
        double log_volume
    ) {
        return problem.evaluate(independent, log_volume);
    };
    diagnostic.envelope = evaluate_held2_pressure_envelope(
        independent_modified_fractions,
        diagnostic.molar_volume_bounds,
        evaluator,
        initial_interval_count,
        8
    );
    diagnostic.parameter_fingerprint = provider.fingerprint();
    return diagnostic;
}

Held2StageIIINlpEvaluation evaluate_held2_installed_stage_iii_derivatives(
    const ProviderContext& provider,
    const FlashInput& input,
    std::size_t phase_count,
    const std::vector<double>& variables,
    const std::vector<double>& equality_multipliers
) {
    require_held2_sdk(provider.sdk());
    validate_input(provider, input);
    const InstalledHeld2Problem problem(provider, input);
    const Held2StateEvaluator evaluator = [&problem](
        const std::vector<double>& independent,
        double log_volume
    ) {
        return problem.evaluate(independent, log_volume);
    };
    return evaluate_held2_stage_iii_nlp(
        problem.coordinates(),
        problem.physical_feed(),
        evaluator,
        phase_count,
        variables,
        equality_multipliers
    );
}

FlashResult solve_tp_flash(
    const ProviderContext& provider,
    const FlashInput& input,
    const SolverResourceProfile& resources,
    Held2ProgressObserver* observer
) {
    require_mixture_sdk(provider.sdk());
    validate_input(provider, input);
    FlashResult result;
    result.input = input;
    result.parameter_fingerprint = provider.fingerprint();
    if (charged(provider.sdk())) {
        require_held2_sdk(provider.sdk());
        result.route = FlashResult::Route::Held2;
        result.held2 = solve_held2(provider, input, resources, observer);
        return result;
    }
    if (input.overall_mole_fractions.size() != 2
        || provider.sdk().component_count != 2
        || provider.fingerprint() != kNeutralFlashFingerprint) {
        throw std::invalid_argument(
            "neutral tp_flash requires the approved two-component fingerprint"
        );
    }
    result.route = FlashResult::Route::Held;
    result.held = solve_held(
        provider,
        input.temperature_k,
        input.pressure_pa,
        input.overall_mole_fractions.front()
    );
    return result;
}

}  // namespace epcsaft_equilibrium
