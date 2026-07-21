#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <epcsaft/native_sdk_v1.h>

#include "held.hpp"
#include "held2.hpp"
#include "saturation.hpp"

namespace py = pybind11;

namespace {

constexpr std::string_view kFlashFingerprint =
    "sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286";
constexpr std::size_t kPureSdkTableSize = offsetof(epcsaft_native_sdk_v1, component_count);
constexpr std::size_t kMixtureSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, evaluate_mixture_phase)
    + sizeof(epcsaft_evaluate_mixture_phase_v1);
constexpr std::size_t kElectrolyteSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, evaluate_electrolyte_phase)
    + sizeof(epcsaft_evaluate_mixture_phase_v1);
constexpr std::size_t kMolarVolumeSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, evaluate_molar_volume_bounds)
    + sizeof(epcsaft_evaluate_molar_volume_bounds_v1);
constexpr std::size_t kPackingSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, evaluate_packing_fraction)
    + sizeof(epcsaft_evaluate_packing_fraction_v1);
constexpr double kGasConstantJPerMolK = 8.31446261815324;

struct Held2ProviderMetadata {
    std::vector<double> charges;
    std::vector<std::string> component_ids;
};

const epcsaft_native_sdk_v1& checked_sdk(const py::capsule& capsule) {
    const char* name = capsule.name();
    if (name == nullptr || std::string_view(name) != EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME) {
        throw py::value_error("expected capsule epcsaft.native_sdk.v1");
    }
    const auto* sdk = capsule.get_pointer<epcsaft_native_sdk_v1>();
    if (sdk == nullptr) {
        throw py::value_error("provider capsule contains a null SDK table");
    }
    if (sdk->abi_version != EPCSAFT_NATIVE_SDK_V1_ABI_VERSION) {
        throw py::value_error("provider capsule ABI version mismatch");
    }
    if (sdk->table_size < kPureSdkTableSize) {
        throw py::value_error("provider capsule table is smaller than the v1 contract");
    }
    if (sdk->result_size != sizeof(epcsaft_phase_block_result_v1)) {
        throw py::value_error("provider capsule result size does not match the v1 contract");
    }
    if (sdk->model_context == nullptr || sdk->evaluate_pure_phase == nullptr) {
        throw py::value_error("provider capsule is missing its model context or phase evaluator");
    }
    return *sdk;
}

const epcsaft_native_sdk_v1& checked_mixture_sdk(const py::capsule& capsule) {
    const epcsaft_native_sdk_v1& sdk = checked_sdk(capsule);
    if (sdk.table_size < kMixtureSdkTableSize) {
        throw py::value_error("provider capsule is missing the mixture SDK tail");
    }
    if (sdk.component_count != 2) {
        throw py::value_error("provider capsule mixture component count must be two");
    }
    if (sdk.mixture_result_size != sizeof(epcsaft_mixture_phase_block_result_v1)) {
        throw py::value_error(
            "provider capsule mixture result size does not match the v1 contract"
        );
    }
    if (sdk.evaluate_mixture_phase == nullptr) {
        throw py::value_error("provider capsule is missing its mixture phase evaluator");
    }
    return sdk;
}

const epcsaft_native_sdk_v1& checked_electrolyte_sdk(const py::capsule& capsule) {
    const epcsaft_native_sdk_v1& sdk = checked_sdk(capsule);
    if (sdk.table_size < kElectrolyteSdkTableSize) {
        throw py::value_error("provider capsule is missing the electrolyte SDK tail");
    }
    if (sdk.component_count < 3 || sdk.component_ids == nullptr
        || sdk.component_charges == nullptr || sdk.evaluate_electrolyte_phase == nullptr) {
        throw py::value_error("provider capsule is missing the electrolyte phase contract");
    }
    if (sdk.mixture_result_size != sizeof(epcsaft_mixture_phase_block_result_v1)) {
        throw py::value_error(
            "provider capsule mixture result size does not match the v1 contract"
        );
    }
    return sdk;
}

const epcsaft_native_sdk_v1& checked_molar_volume_sdk(const py::capsule& capsule) {
    const epcsaft_native_sdk_v1& sdk = checked_electrolyte_sdk(capsule);
    if (sdk.table_size < kMolarVolumeSdkTableSize
        || sdk.evaluate_molar_volume_bounds == nullptr) {
        throw py::value_error("provider capsule is missing the molar-volume domain contract");
    }
    return sdk;
}

const epcsaft_native_sdk_v1& checked_packing_sdk(const py::capsule& capsule) {
    const epcsaft_native_sdk_v1& sdk = checked_molar_volume_sdk(capsule);
    if (sdk.table_size < kPackingSdkTableSize
        || sdk.evaluate_packing_fraction == nullptr) {
        throw py::value_error("provider capsule is missing the packing-fraction contract");
    }
    return sdk;
}

Held2ProviderMetadata held2_provider_metadata(const epcsaft_native_sdk_v1& sdk) {
    Held2ProviderMetadata result;
    result.charges.reserve(sdk.component_count);
    result.component_ids.reserve(sdk.component_count);
    for (std::size_t component = 0; component < sdk.component_count; ++component) {
        if (sdk.component_ids[component] == nullptr
            || sdk.component_ids[component][0] == '\0') {
            throw py::value_error("provider electrolyte component ID must not be empty");
        }
        result.charges.push_back(static_cast<double>(sdk.component_charges[component]));
        result.component_ids.emplace_back(sdk.component_ids[component]);
    }
    return result;
}

py::dict sdk_info(const py::capsule& capsule) {
    const epcsaft_native_sdk_v1& sdk = checked_sdk(capsule);
    const bool has_mixture_tail = sdk.table_size >= kMixtureSdkTableSize;
    py::dict result;
    result["capsule_name"] = EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME;
    result["abi_version"] = sdk.abi_version;
    result["table_size"] = sdk.table_size;
    result["mixture_prefix_size"] = kMixtureSdkTableSize;
    result["result_size"] = sdk.result_size;
    result["has_model_context"] = sdk.model_context != nullptr;
    result["has_evaluate_pure_phase"] = sdk.evaluate_pure_phase != nullptr;
    result["component_count"] = has_mixture_tail ? sdk.component_count : 0;
    result["mixture_result_size"] = has_mixture_tail ? sdk.mixture_result_size : 0;
    result["has_evaluate_mixture_phase"] =
        has_mixture_tail && sdk.evaluate_mixture_phase != nullptr;
    return result;
}

py::dict evaluate_mixture_phase(
    const py::capsule& capsule,
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3,
    const std::string& expected_fingerprint
) {
    const epcsaft_native_sdk_v1& sdk = checked_mixture_sdk(capsule);
    const epcsaft_equilibrium::ProviderContext provider(sdk, expected_fingerprint);
    const epcsaft_equilibrium::MixturePhaseEvaluation phase =
        provider.evaluate_mixture(temperature_k, amounts_mol, volume_m3);
    py::dict result;
    result["helmholtz_over_rt_reference_amount"] = phase.value;
    result["gradient"] = phase.gradient;
    result["hessian"] = phase.hessian;
    result["pressure_pa"] = phase.pressure_pa;
    result["parameter_fingerprint"] = phase.parameter_fingerprint;
    return result;
}

py::dict held_state_to_dict(const epcsaft_equilibrium::HeldStateEvaluation& state) {
    py::dict result;
    result["x_methane"] = state.x_methane;
    result["log_volume"] = state.log_volume;
    result["volume_m3"] = state.volume_m3;
    result["amounts_mol"] = state.amounts_mol;
    result["g_bar"] = state.g_bar;
    result["gradient"] = state.gradient;
    result["hessian"] = state.hessian;
    result["pressure_pa"] = state.provider.pressure_pa;
    result["pressure_stationarity_relative"] = state.pressure_stationarity_relative;
    result["parameter_fingerprint"] = state.provider.parameter_fingerprint;
    return result;
}

py::dict held_evaluate_state(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    double x_methane,
    double log_volume,
    const std::string& expected_fingerprint
) {
    const epcsaft_equilibrium::ProviderContext provider(
        checked_mixture_sdk(capsule),
        expected_fingerprint
    );
    return held_state_to_dict(epcsaft_equilibrium::evaluate_held_state(
        provider,
        temperature_k,
        pressure_pa,
        x_methane,
        log_volume
    ));
}

py::dict held_evaluate_tpd(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    double reference_x_methane,
    double reference_log_volume,
    double x_methane,
    double log_volume,
    const std::string& expected_fingerprint
) {
    const epcsaft_equilibrium::ProviderContext provider(
        checked_mixture_sdk(capsule),
        expected_fingerprint
    );
    const epcsaft_equilibrium::HeldStateEvaluation reference =
        epcsaft_equilibrium::evaluate_held_state(
            provider,
            temperature_k,
            pressure_pa,
            reference_x_methane,
            reference_log_volume
        );
    const epcsaft_equilibrium::HeldTpdEvaluation evaluation =
        epcsaft_equilibrium::evaluate_held_tpd(
            provider,
            temperature_k,
            pressure_pa,
            reference,
            x_methane,
            log_volume
        );
    py::dict result;
    result["d_bar"] = evaluation.d_bar;
    result["gradient"] = evaluation.gradient;
    result["hessian"] = evaluation.hessian;
    result["state"] = held_state_to_dict(evaluation.state);
    return result;
}

py::dict held_evaluate_tunneling(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    double reference_x_methane,
    double reference_log_volume,
    double minimum_x_methane,
    double minimum_log_volume,
    double x_methane,
    double log_volume,
    const std::string& expected_fingerprint
) {
    const epcsaft_equilibrium::ProviderContext provider(
        checked_mixture_sdk(capsule),
        expected_fingerprint
    );
    const epcsaft_equilibrium::HeldStateEvaluation reference =
        epcsaft_equilibrium::evaluate_held_state(
            provider,
            temperature_k,
            pressure_pa,
            reference_x_methane,
            reference_log_volume
        );
    const epcsaft_equilibrium::HeldTpdEvaluation minimum =
        epcsaft_equilibrium::evaluate_held_tpd(
            provider,
            temperature_k,
            pressure_pa,
            reference,
            minimum_x_methane,
            minimum_log_volume
        );
    const epcsaft_equilibrium::HeldTunnelingEvaluation evaluation =
        epcsaft_equilibrium::evaluate_held_tunneling(
            provider,
            temperature_k,
            pressure_pa,
            reference,
            minimum_x_methane,
            minimum.d_bar,
            x_methane,
            log_volume
        );
    py::dict result;
    result["objective"] = evaluation.objective;
    result["gradient"] = evaluation.gradient;
    result["hessian"] = evaluation.hessian;
    result["d_bar"] = evaluation.tpd.d_bar;
    return result;
}

py::dict held_stage_i(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane,
    const std::string& expected_fingerprint
) {
    const epcsaft_equilibrium::ProviderContext provider(
        checked_mixture_sdk(capsule),
        expected_fingerprint
    );
    epcsaft_equilibrium::HeldStageIResult solve;
    {
        py::gil_scoped_release release;
        solve = epcsaft_equilibrium::solve_held_stage_i(
            provider,
            temperature_k,
            pressure_pa,
            feed_x_methane
        );
    }

    py::list reference_attempts;
    for (const auto& attempt : solve.reference_attempts) {
        py::dict item;
        item["role"] = attempt.role;
        item["initial_log_volume"] = attempt.initial_log_volume;
        item["solver_converged"] = attempt.solver_converged;
        item["solver_status"] = attempt.solver_status;
        item["iterations"] = attempt.iterations;
        item["accepted"] = attempt.accepted;
        item["g_bar"] = attempt.g_bar;
        item["log_volume"] = attempt.log_volume;
        item["pressure_stationarity_relative"] =
            attempt.pressure_stationarity_relative;
        item["callback_error"] = attempt.callback_error;
        reference_attempts.append(std::move(item));
    }
    py::list planned_starts;
    for (const auto& start : solve.planned_starts) {
        py::dict item;
        item["role"] = start.role;
        item["x_methane"] = start.x_methane;
        item["log_volume"] = start.log_volume;
        planned_starts.append(std::move(item));
    }
    py::list attempt_log;
    for (const auto& attempt : solve.attempt_log) {
        py::dict item;
        item["kind"] = attempt.kind;
        item["role"] = attempt.role;
        item["initial_guess"] = attempt.initial_guess;
        item["solver_converged"] = attempt.solver_converged;
        item["solver_status"] = attempt.solver_status;
        item["iterations"] = attempt.iterations;
        item["accepted"] = attempt.accepted;
        item["materially_perturbed"] = attempt.materially_perturbed;
        item["objective"] = attempt.objective;
        item["tpd"] = attempt.tpd;
        item["pressure_stationarity_relative"] =
            attempt.pressure_stationarity_relative;
        item["callback_error"] = attempt.callback_error;
        attempt_log.append(std::move(item));
    }

    py::dict result;
    result["outcome"] = solve.outcome;
    result["search_status"] = solve.search_status;
    result["search_profile"] = solve.search_profile;
    result["globality_certificate"] = "not_guaranteed";
    result["failure_reason"] = solve.failure_reason;
    result["reference_attempts"] = reference_attempts;
    result["planned_starts"] = planned_starts;
    result["attempt_log"] = attempt_log;
    result["starts_completed"] = solve.starts_completed;
    result["negative_confirmations"] = solve.negative_confirmations;
    result["confirmation_max_difference"] = solve.confirmation_max_difference;
    result["best_tpd"] = solve.best_tpd;
    if (solve.has_reference) {
        result["reference"] = held_state_to_dict(solve.reference);
        result["best_state"] = held_state_to_dict(solve.best_state);
    }
    return result;
}

py::dict held_outer_to_dict(const epcsaft_equilibrium::HeldOuterResult& solve) {
    py::dict result;
    result["status"] = solve.status;
    result["failure_reason"] = solve.failure_reason;
    result["active_cut_ids"] = solve.active_cut_ids;
    result["tied_multipliers"] = solve.tied_multipliers;
    if (solve.status == "finite") {
        result["value"] = solve.value;
        result["multiplier"] = solve.multiplier;
    }
    return result;
}

py::dict held_outer_envelope(
    const std::vector<std::tuple<std::string, double, double>>& raw_cuts
) {
    std::vector<epcsaft_equilibrium::HeldOuterCut> cuts;
    cuts.reserve(raw_cuts.size());
    for (const auto& [identity, intercept, slope] : raw_cuts) {
        cuts.push_back({identity, intercept, slope});
    }
    return held_outer_to_dict(epcsaft_equilibrium::solve_held_outer_envelope(cuts));
}

py::dict held_stage_ii_initial_cuts(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane,
    const std::string& expected_fingerprint
) {
    const epcsaft_equilibrium::ProviderContext provider(
        checked_mixture_sdk(capsule),
        expected_fingerprint
    );
    epcsaft_equilibrium::HeldStageIIInitialization solve;
    {
        py::gil_scoped_release release;
        solve = epcsaft_equilibrium::initialize_held_stage_ii_cuts(
            provider,
            temperature_k,
            pressure_pa,
            feed_x_methane
        );
    }
    py::list attempts;
    for (const auto& attempt : solve.attempts) {
        py::dict item;
        item["role"] = attempt.role;
        item["initial_guess"] = attempt.initial_guess;
        item["solver_converged"] = attempt.solver_converged;
        item["solver_status"] = attempt.solver_status;
        item["iterations"] = attempt.iterations;
        item["accepted"] = attempt.accepted;
        item["objective"] = attempt.objective;
        item["pressure_stationarity_relative"] =
            attempt.pressure_stationarity_relative;
        item["callback_error"] = attempt.callback_error;
        attempts.append(std::move(item));
    }
    py::list cuts;
    for (const auto& cut : solve.cuts) {
        py::dict item;
        item["identity"] = cut.identity;
        item["endpoint"] = cut.endpoint;
        item["intercept"] = cut.outer.intercept;
        item["slope"] = cut.outer.slope;
        item["x_methane"] = cut.state.x_methane;
        item["log_volume"] = cut.state.log_volume;
        item["volume_m3"] = cut.state.volume_m3;
        item["pressure_stationarity_relative"] =
            cut.state.pressure_stationarity_relative;
        cuts.append(std::move(item));
    }
    py::dict result;
    result["status"] = solve.status;
    result["failure_reason"] = solve.failure_reason;
    result["attempts"] = attempts;
    result["cuts"] = cuts;
    result["outer"] = held_outer_to_dict(solve.outer);
    return result;
}

py::dict held_evaluate_lower(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane,
    double multiplier,
    double x_methane,
    double log_volume,
    const std::string& expected_fingerprint
) {
    const epcsaft_equilibrium::ProviderContext provider(
        checked_mixture_sdk(capsule),
        expected_fingerprint
    );
    const epcsaft_equilibrium::HeldLowerEvaluation evaluation =
        epcsaft_equilibrium::evaluate_held_lower(
            provider,
            temperature_k,
            pressure_pa,
            feed_x_methane,
            multiplier,
            x_methane,
            log_volume
        );
    py::dict result;
    result["objective"] = evaluation.objective;
    result["gradient"] = evaluation.gradient;
    result["hessian"] = evaluation.hessian;
    result["state"] = held_state_to_dict(evaluation.state);
    return result;
}

py::dict held_stage_ii_lower_search(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane,
    double multiplier,
    double upper_bound,
    const std::vector<std::array<double, 2>>& raw_previous_states,
    const std::string& expected_fingerprint
) {
    const epcsaft_equilibrium::ProviderContext provider(
        checked_mixture_sdk(capsule),
        expected_fingerprint
    );
    std::vector<epcsaft_equilibrium::HeldStateEvaluation> previous_states;
    previous_states.reserve(raw_previous_states.size());
    for (const auto& values : raw_previous_states) {
        epcsaft_equilibrium::HeldStateEvaluation state;
        state.x_methane = values[0];
        state.log_volume = values[1];
        previous_states.push_back(std::move(state));
    }
    epcsaft_equilibrium::HeldStageIILowerSearch solve;
    {
        py::gil_scoped_release release;
        solve = epcsaft_equilibrium::search_held_stage_ii_lower(
            provider,
            temperature_k,
            pressure_pa,
            feed_x_methane,
            multiplier,
            upper_bound,
            previous_states
        );
    }
    py::list planned_starts;
    for (const auto& start : solve.planned_starts) {
        py::dict item;
        item["class"] = start.start_class;
        item["role"] = start.role;
        item["x_methane"] = start.x_methane;
        item["log_volume"] = start.log_volume;
        planned_starts.append(std::move(item));
    }
    py::list attempts;
    for (const auto& attempt : solve.attempts) {
        py::dict item;
        item["role"] = attempt.role;
        item["initial_guess"] = attempt.initial_guess;
        item["solver_converged"] = attempt.solver_converged;
        item["solver_status"] = attempt.solver_status;
        item["iterations"] = attempt.iterations;
        item["accepted"] = attempt.accepted;
        item["objective"] = attempt.objective;
        item["pressure_stationarity_relative"] =
            attempt.pressure_stationarity_relative;
        item["callback_error"] = attempt.callback_error;
        attempts.append(std::move(item));
    }
    py::list cuts;
    for (const auto& cut : solve.accepted_cuts) {
        py::dict item;
        item["identity"] = cut.identity;
        item["g_bar"] = cut.state.g_bar;
        item["x_methane"] = cut.state.x_methane;
        item["log_volume"] = cut.state.log_volume;
        item["volume_m3"] = cut.state.volume_m3;
        item["composition_gradient"] = cut.state.gradient[0];
        item["intercept"] = cut.outer.intercept;
        item["slope"] = cut.outer.slope;
        item["objective"] = cut.outer.intercept + multiplier * cut.outer.slope;
        item["pressure_stationarity_relative"] =
            cut.state.pressure_stationarity_relative;
        cuts.append(std::move(item));
    }
    py::dict result;
    result["status"] = solve.status;
    result["upper_bound"] = solve.upper_bound;
    result["multiplier"] = solve.multiplier;
    result["starts_completed"] = solve.starts_completed;
    result["planned_starts"] = planned_starts;
    result["attempts"] = attempts;
    result["accepted_cuts"] = cuts;
    return result;
}

py::dict held_stage_ii_candidates(
    double feed_x_methane,
    double upper_bound,
    double multiplier,
    const std::vector<std::tuple<std::string, double, double, double, double>>&
        raw_points
) {
    std::vector<epcsaft_equilibrium::HeldCandidateInput> points;
    points.reserve(raw_points.size());
    for (const auto& [identity, g_bar, x_methane, volume_m3, gradient] : raw_points) {
        points.push_back({identity, g_bar, x_methane, volume_m3, gradient});
    }
    const epcsaft_equilibrium::HeldStageIICandidateResult solve =
        epcsaft_equilibrium::select_held_stage_ii_candidates(
            feed_x_methane,
            upper_bound,
            multiplier,
            points
        );
    py::list candidate_ids;
    for (const auto& candidate : solve.candidates) {
        candidate_ids.append(candidate.identity);
    }
    py::list rejections;
    for (const auto& rejection : solve.rejections) {
        py::dict item;
        item["identity"] = rejection.identity;
        item["reason"] = rejection.reason;
        rejections.append(std::move(item));
    }
    py::dict result;
    result["status"] = solve.status;
    result["candidate_ids"] = candidate_ids;
    result["rejections"] = rejections;
    return result;
}

py::dict held_stage_ii_trace_to_dict(
    const epcsaft_equilibrium::HeldStageIITrace& entry
) {
    py::list rejections;
    for (const auto& rejection : entry.rejections) {
        py::dict item;
        item["identity"] = rejection.identity;
        item["reason"] = rejection.reason;
        rejections.append(std::move(item));
    }
    py::dict result;
    result["major_iteration"] = entry.major_iteration;
    result["outer_value"] = entry.outer_value;
    result["upper_bound"] = entry.upper_bound;
    result["multiplier"] = entry.multiplier;
    result["active_cut_ids"] = entry.active_cut_ids;
    result["accepted_cut_ids"] = entry.accepted_cut_ids;
    result["lower_starts_completed"] = entry.lower_starts_completed;
    result["candidate_ids"] = entry.candidate_ids;
    result["rejections"] = rejections;
    result["stage_iii_outcome"] = entry.stage_iii_outcome;
    result["stage_iii_failure_reason"] = entry.stage_iii_failure_reason;
    return result;
}

py::dict held_stage_ii(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane,
    const std::string& expected_fingerprint
) {
    const epcsaft_equilibrium::ProviderContext provider(
        checked_mixture_sdk(capsule),
        expected_fingerprint
    );
    epcsaft_equilibrium::HeldStageIIResult solve;
    {
        py::gil_scoped_release release;
        solve = epcsaft_equilibrium::solve_held_stage_ii(
            provider,
            temperature_k,
            pressure_pa,
            feed_x_methane
        );
    }
    const auto cut_to_dict = [](const epcsaft_equilibrium::HeldStageIICut& cut) {
        py::dict item;
        item["identity"] = cut.identity;
        item["endpoint"] = cut.endpoint;
        item["g_bar"] = cut.state.g_bar;
        item["x_methane"] = cut.state.x_methane;
        item["log_volume"] = cut.state.log_volume;
        item["volume_m3"] = cut.state.volume_m3;
        item["composition_gradient"] = cut.state.gradient[0];
        item["pressure_stationarity_relative"] =
            cut.state.pressure_stationarity_relative;
        item["intercept"] = cut.outer.intercept;
        item["slope"] = cut.outer.slope;
        return item;
    };
    py::list endpoint_attempts;
    for (const auto& attempt : solve.endpoint_attempts) {
        py::dict item;
        item["role"] = attempt.role;
        item["solver_status"] = attempt.solver_status;
        item["iterations"] = attempt.iterations;
        item["accepted"] = attempt.accepted;
        item["callback_error"] = attempt.callback_error;
        endpoint_attempts.append(std::move(item));
    }
    py::list cuts;
    for (const auto& cut : solve.cuts) {
        cuts.append(cut_to_dict(cut));
    }
    py::list candidates;
    for (const auto& candidate : solve.candidates) {
        candidates.append(cut_to_dict(candidate));
    }
    py::list trace;
    for (const auto& entry : solve.trace) {
        trace.append(held_stage_ii_trace_to_dict(entry));
    }
    py::dict result;
    result["outcome"] = solve.outcome;
    result["search_status"] = solve.search_status;
    result["search_profile"] = solve.search_profile;
    result["globality_certificate"] = "not_guaranteed";
    result["failure_reason"] = solve.failure_reason;
    result["stage_i_outcome"] = solve.stage_i_outcome;
    result["stage_i_search_status"] = solve.stage_i_search_status;
    result["best_tpd"] = solve.best_tpd;
    result["major_iterations"] = solve.major_iterations;
    result["upper_bound"] = solve.upper_bound;
    result["endpoint_attempts"] = endpoint_attempts;
    result["cuts"] = cuts;
    result["candidates"] = candidates;
    result["trace"] = trace;
    return result;
}

std::vector<epcsaft_equilibrium::HeldStageIIICandidate> held_stage_iii_candidates(
    const std::vector<std::tuple<std::string, double, double>>& raw_candidates
) {
    std::vector<epcsaft_equilibrium::HeldStageIIICandidate> candidates;
    candidates.reserve(raw_candidates.size());
    for (const auto& [identity, x_methane, volume_m3_per_mol] : raw_candidates) {
        candidates.push_back({identity, x_methane, volume_m3_per_mol});
    }
    return candidates;
}

py::dict held_stage_iii_initialization_to_dict(
    const epcsaft_equilibrium::HeldStageIIIInitialization& initialization
) {
    py::dict result;
    result["status"] = initialization.status;
    result["failure_reason"] = initialization.failure_reason;
    result["phase_fractions"] = initialization.phase_fractions;
    result["initial_variables"] = initialization.has_variables
        ? py::cast(initialization.initial_variables)
        : py::cast(std::vector<double>{});
    result["lower_bounds"] = initialization.has_variables
        ? py::cast(initialization.lower_bounds)
        : py::cast(std::vector<double>{});
    result["upper_bounds"] = initialization.has_variables
        ? py::cast(initialization.upper_bounds)
        : py::cast(std::vector<double>{});
    return result;
}

py::dict held_stage_iii_initialize(
    double feed_x_methane,
    const std::vector<std::tuple<std::string, double, double>>& raw_candidates
) {
    return held_stage_iii_initialization_to_dict(
        epcsaft_equilibrium::initialize_held_stage_iii(
            feed_x_methane,
            held_stage_iii_candidates(raw_candidates)
        )
    );
}

py::dict held_stage_iii_classify(const py::dict& raw) {
    epcsaft_equilibrium::HeldStageIIIAcceptanceEvidence evidence;
    evidence.solver_converged = raw["solver_converged"].cast<bool>();
    evidence.callback_error = raw["callback_error"].cast<std::string>();
    evidence.solver_constraint_violation =
        raw["solver_constraint_violation"].cast<double>();
    evidence.material_balance_max_abs = raw["material_balance_max_abs"].cast<double>();
    evidence.pressure_stationarity_max_relative =
        raw["pressure_stationarity_max_relative"].cast<double>();
    evidence.kkt_stationarity_max_abs = raw["kkt_stationarity_max_abs"].cast<double>();
    evidence.inactive_bounds = raw["inactive_bounds"].cast<bool>();
    evidence.composition_distance = raw["composition_distance"].cast<double>();
    evidence.phase_density_distance = raw["phase_density_distance"].cast<double>();
    evidence.held_gap = raw["held_gap"].cast<double>();
    evidence.chemical_potential_max_relative =
        raw["chemical_potential_max_relative"].cast<double>();
    evidence.confirmation_succeeded = raw["confirmation_succeeded"].cast<bool>();
    evidence.confirmation_max_difference =
        raw["confirmation_max_difference"].cast<double>();
    const auto decision = epcsaft_equilibrium::classify_held_stage_iii(evidence);
    py::dict result;
    result["status"] = decision.status;
    result["failure_reason"] = decision.failure_reason;
    return result;
}

double held_stage_iii_mu_difference(
    const std::array<double, 2>& first,
    const std::array<double, 2>& second
) {
    return epcsaft_equilibrium::held_stage_iii_mu_difference(first, second);
}

py::dict held_evaluate_stage_iii(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane,
    const std::array<double, 6>& variables,
    const std::string& expected_fingerprint
) {
    const epcsaft_equilibrium::ProviderContext provider(
        checked_mixture_sdk(capsule),
        expected_fingerprint
    );
    const auto evaluation = epcsaft_equilibrium::evaluate_held_stage_iii_nlp(
        provider,
        temperature_k,
        pressure_pa,
        {feed_x_methane, 1.0 - feed_x_methane},
        variables
    );
    py::dict result;
    result["objective"] = evaluation.objective;
    result["gradient"] = evaluation.gradient;
    result["constraints"] = evaluation.constraints;
    result["jacobian"] = evaluation.jacobian;
    result["hessian_lower"] = evaluation.hessian_lower;
    return result;
}

py::dict held_stage_iii_phase_to_dict(
    const epcsaft_equilibrium::HeldStageIIIPhaseEvaluation& phase
);

py::list held_stage_iii_attempt_log_to_list(
    const std::vector<epcsaft_equilibrium::HeldStageIIILocalAttempt>& attempts
);

py::dict held_stage_iii_to_dict(
    const epcsaft_equilibrium::HeldStageIIIResult& solve
) {
    py::dict result;
    result["outcome"] = solve.outcome;
    result["search_profile"] = solve.search_profile;
    result["failure_reason"] = solve.failure_reason;
    result["initialization"] = held_stage_iii_initialization_to_dict(solve.initialization);
    result["solver_converged"] = solve.local.solver_converged;
    result["solver_status"] = solve.local.solver_status;
    result["iterations"] = solve.local.iterations;
    result["attempts"] = solve.local.attempts;
    result["attempt_log"] = held_stage_iii_attempt_log_to_list(solve.local.attempt_log);
    result["solver_lower_bounds"] = solve.local.solver_lower_bounds;
    result["solver_upper_bounds"] = solve.local.solver_upper_bounds;
    result["solver_constraint_violation"] = solve.local.solver_constraint_violation;
    result["confirmation_solves"] = solve.local.confirmation_solves;
    result["confirmation_succeeded"] = solve.confirmation_succeeded;
    result["confirmation_max_difference"] = solve.local.confirmation_max_difference;
    result["material_balance_max_abs"] = solve.local.material_balance_max_abs;
    result["pressure_stationarity_max_relative"] =
        solve.local.pressure_stationarity_max_relative;
    result["kkt_stationarity_max_abs"] = solve.local.kkt_stationarity_max_abs;
    result["inactive_bounds"] = solve.inactive_bounds;
    result["composition_distance"] = solve.composition_distance;
    result["phase_density_distance"] = solve.local.phase_density_distance;
    result["chemical_potential_max_relative"] = solve.chemical_potential_max_relative;
    result["held_gap"] = solve.held_gap;
    result["upper_bound"] = solve.upper_bound;
    result["equality_multipliers"] = solve.local.equality_multipliers;
    result["lower_bound_multipliers"] = solve.local.lower_bound_multipliers;
    result["upper_bound_multipliers"] = solve.local.upper_bound_multipliers;
    result["total_g_bar"] = solve.local.evaluation.objective;
    if (solve.local.accepted) {
        result["first_phase"] = held_stage_iii_phase_to_dict(
            solve.local.evaluation.liquid
        );
        result["second_phase"] = held_stage_iii_phase_to_dict(
            solve.local.evaluation.vapor
        );
    }
    return result;
}

py::dict held_stage_iii(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane,
    double upper_bound,
    const std::vector<std::tuple<std::string, double, double>>& raw_candidates,
    const std::string& expected_fingerprint
) {
    const epcsaft_equilibrium::ProviderContext provider(
        checked_mixture_sdk(capsule),
        expected_fingerprint
    );
    epcsaft_equilibrium::HeldStageIIIResult solve;
    {
        py::gil_scoped_release release;
        solve = epcsaft_equilibrium::solve_held_stage_iii(
            provider,
            temperature_k,
            pressure_pa,
            feed_x_methane,
            upper_bound,
            held_stage_iii_candidates(raw_candidates)
        );
    }
    return held_stage_iii_to_dict(solve);
}

py::dict solve_tp_flash(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::array<double, 2>& overall_mole_fractions
) {
    if (!std::isfinite(temperature_k) || temperature_k <= 0.0
        || !std::isfinite(pressure_pa) || pressure_pa <= 0.0
        || !std::isfinite(overall_mole_fractions[0])
        || !std::isfinite(overall_mole_fractions[1])
        || overall_mole_fractions[0] <= 0.0
        || overall_mole_fractions[1] <= 0.0
        || std::abs(overall_mole_fractions[0] + overall_mole_fractions[1] - 1.0)
            > 1.0e-12) {
        throw py::value_error("tp_flash native inputs must be positive, finite, and normalized");
    }
    const epcsaft_equilibrium::ProviderContext provider(
        checked_mixture_sdk(capsule),
        std::string(kFlashFingerprint)
    );
    epcsaft_equilibrium::HeldResult solve;
    {
        py::gil_scoped_release release;
        solve = epcsaft_equilibrium::solve_held(
            provider,
            temperature_k,
            pressure_pa,
            overall_mole_fractions[0]
        );
    }
    int attempt_count = static_cast<int>(
        solve.stage_i.reference_attempts.size() + solve.stage_i.attempt_log.size()
        + solve.stage_ii.endpoint_attempts.size()
    );
    for (const auto& entry : solve.stage_ii.trace) {
        attempt_count += entry.lower_starts_completed;
    }
    py::list attempts;
    for (const auto& attempt : solve.stage_iii_attempts) {
        attempt_count += attempt.result.local.attempts;
        py::dict item = held_stage_iii_to_dict(attempt.result);
        item["major_iteration"] = attempt.major_iteration;
        attempts.append(std::move(item));
    }
    py::list trace;
    for (const auto& entry : solve.stage_ii.trace) {
        trace.append(held_stage_ii_trace_to_dict(entry));
    }
    py::list search_profiles;
    search_profiles.append(solve.stage_i.search_profile);
    if (solve.stage_i.outcome == "negative_tpd") {
        search_profiles.append(solve.stage_ii.search_profile);
    }
    if (!solve.stage_iii_attempts.empty()) {
        search_profiles.append(solve.stage_iii_attempts.front().result.search_profile);
    }
    py::dict result;
    result["outcome"] = solve.outcome;
    result["search_status"] = solve.search_status;
    result["failure_reason"] = solve.failure_reason;
    result["stage_i_outcome"] = solve.stage_i.outcome;
    result["stage_i_search_status"] = solve.stage_i.search_status;
    result["attempts"] = attempt_count;
    result["major_iterations"] = solve.stage_ii.major_iterations;
    result["best_tpd"] = solve.stage_i.best_tpd;
    result["lower_bound"] = py::none();
    result["upper_bound"] = py::none();
    result["held_gap"] = py::none();
    result["material_balance_max_abs"] = py::none();
    result["pressure_stationarity_max_relative"] = py::none();
    result["kkt_stationarity_max_abs"] = py::none();
    result["chemical_potential_max_relative"] = py::none();
    result["confirmation_succeeded"] = false;
    result["confirmation_max_difference"] = py::none();
    result["search_profiles"] = search_profiles;
    result["stage_iii_attempts"] = attempts;
    result["trace"] = trace;
    result["globality_certificate"] = "not_guaranteed";

    result["solver_status"] = "not_adjudicated";
    result["numerical_status"] = "not_adjudicated";
    result["physical_status"] = "not_adjudicated";
    const bool accepted_reference = solve.stage_i.has_reference
        && std::any_of(
            solve.stage_i.reference_attempts.begin(),
            solve.stage_i.reference_attempts.end(),
            [](const auto& attempt) { return attempt.accepted; }
        );
    const bool failed_reference_solver_attempt = std::any_of(
        solve.stage_i.reference_attempts.begin(),
        solve.stage_i.reference_attempts.end(),
        [](const auto& attempt) {
            return !attempt.solver_converged || !attempt.callback_error.empty();
        }
    );
    const bool failed_reference_attempt = std::any_of(
        solve.stage_i.reference_attempts.begin(),
        solve.stage_i.reference_attempts.end(),
        [](const auto& attempt) { return !attempt.accepted; }
    );
    const bool failed_stage_i_solver_attempt = std::any_of(
        solve.stage_i.attempt_log.begin(),
        solve.stage_i.attempt_log.end(),
        [](const auto& attempt) {
            return !attempt.solver_converged || !attempt.callback_error.empty();
        }
    );
    const bool failed_stage_i_search_attempt = std::any_of(
        solve.stage_i.attempt_log.begin(),
        solve.stage_i.attempt_log.end(),
        [](const auto& attempt) {
            return !attempt.solver_converged
                || !attempt.callback_error.empty()
                || !attempt.accepted;
        }
    );
    const bool all_stage_i_starts_completed = solve.stage_i.planned_starts.size() == 20
        && solve.stage_i.starts_completed == 20;
    const bool no_stage_i_confirmation = solve.stage_i.negative_confirmations == 0
        && std::none_of(
            solve.stage_i.attempt_log.begin(),
            solve.stage_i.attempt_log.end(),
            [](const auto& attempt) { return attempt.kind == "confirmation"; }
        );
    const bool no_negative_tpd = std::isfinite(solve.stage_i.best_tpd)
        && solve.stage_i.best_tpd >= -1.0e-8;
    const bool completed_one_phase_search = accepted_reference
        && all_stage_i_starts_completed
        && no_stage_i_confirmation
        && no_negative_tpd;
    if (completed_one_phase_search) {
        const bool solver_passed = !failed_reference_solver_attempt
            && !failed_stage_i_solver_attempt;
        const bool numerical_passed = solver_passed
            && !failed_reference_attempt
            && !failed_stage_i_search_attempt;
        result["solver_status"] = solver_passed ? "passed" : "failed";
        result["numerical_status"] = numerical_passed ? "passed" : "failed";
        result["physical_status"] = numerical_passed ? "passed" : "not_adjudicated";
    } else if (solve.outcome == "accepted" && !solve.stage_iii_attempts.empty()) {
        const auto& local = solve.stage_iii_attempts.back().result.local;
        result["solver_status"] = local.solver_converged ? "passed" : "failed";
        result["numerical_status"] = local.numerical_converged ? "passed" : "failed";
        result["physical_status"] = local.physical_accepted ? "passed" : "failed";
    }
    result["temperature_k"] = temperature_k;
    result["pressure_pa"] = pressure_pa;
    result["overall_mole_fractions"] = overall_mole_fractions;
    result["parameter_fingerprint"] = std::string(kFlashFingerprint);

    if (solve.outcome == "one_phase" && solve.stage_i.has_reference) {
        const auto& state = solve.stage_i.reference;
        py::dict phase;
        phase["amount_mol"] = 1.0;
        phase["mole_fractions"] = state.amounts_mol;
        phase["volume_m3"] = state.volume_m3;
        phase["molar_density_mol_m3"] = 1.0 / state.volume_m3;
        phase["pressure_pa"] = state.provider.pressure_pa;
        phase["chemical_potential_over_rt"] = std::array<double, 2>{
            state.provider.gradient[0],
            state.provider.gradient[1],
        };
        py::list phases;
        phases.append(std::move(phase));
        result["phases"] = phases;
        result["phase_fractions"] = std::array<double, 1>{1.0};
        result["total_free_energy_over_rt"] = state.g_bar;
        result["pressure_stationarity_max_relative"] =
            std::abs(state.pressure_stationarity_relative);
    } else if (!solve.stage_iii_attempts.empty()) {
        const auto& refinement = solve.stage_iii_attempts.back().result;
        const auto& local = refinement.local;
        result["upper_bound"] = refinement.upper_bound;
        if (local.solver_converged) {
            result["lower_bound"] = local.evaluation.objective;
            result["held_gap"] = refinement.held_gap;
            result["material_balance_max_abs"] = local.material_balance_max_abs;
            result["pressure_stationarity_max_relative"] =
                local.pressure_stationarity_max_relative;
            result["kkt_stationarity_max_abs"] = local.kkt_stationarity_max_abs;
            result["chemical_potential_max_relative"] =
                refinement.chemical_potential_max_relative;
        }
        if (local.confirmation_solves > 0) {
            result["confirmation_succeeded"] = refinement.confirmation_succeeded;
            result["confirmation_max_difference"] = local.confirmation_max_difference;
        }
        if (solve.outcome == "accepted") {
            py::list phases;
            phases.append(held_stage_iii_phase_to_dict(local.evaluation.liquid));
            phases.append(held_stage_iii_phase_to_dict(local.evaluation.vapor));
            result["phases"] = phases;
            result["phase_fractions"] = std::array<double, 2>{
                local.evaluation.liquid.amounts_mol[0]
                    + local.evaluation.liquid.amounts_mol[1],
                local.evaluation.vapor.amounts_mol[0]
                    + local.evaluation.vapor.amounts_mol[1],
            };
            result["total_free_energy_over_rt"] = local.evaluation.objective;
            result["accepted_stage_iii"] = held_stage_iii_to_dict(refinement);
        }
    } else if (!solve.stage_ii.trace.empty()) {
        result["lower_bound"] = solve.stage_ii.trace.back().outer_value;
        result["upper_bound"] = solve.stage_ii.upper_bound;
    }
    return result;
}

py::dict evaluate_phase(
    const py::capsule& capsule,
    double temperature_k,
    double amount_mol,
    double volume_m3,
    const std::string& expected_fingerprint
) {
    const epcsaft_native_sdk_v1& sdk = checked_sdk(capsule);
    const epcsaft_equilibrium::ProviderContext provider(sdk, expected_fingerprint);
    const epcsaft_equilibrium::PhaseEvaluation evaluation =
        provider.evaluate(temperature_k, amount_mol, volume_m3);
    const epcsaft_phase_block_result_v1& phase = evaluation.provider;

    py::dict result;
    result["status"] = phase.status;
    result["amount_mol"] = amount_mol;
    result["volume_m3"] = volume_m3;
    result["helmholtz_over_rt_reference_amount"] = phase.helmholtz_over_rt_reference_amount;
    result["gradient"] = std::array<double, 2>{phase.gradient[0], phase.gradient[1]};
    result["hessian"] = std::array<double, 4>{
        phase.hessian[0], phase.hessian[1], phase.hessian[2], phase.hessian[3]
    };
    result["third"] = std::array<double, 8>{
        phase.third[0], phase.third[1], phase.third[2], phase.third[3],
        phase.third[4], phase.third[5], phase.third[6], phase.third[7]
    };
    result["pressure_pa"] = phase.pressure_pa;
    result["chemical_potential_over_rt"] = phase.chemical_potential_over_rt;
    result["parameter_fingerprint"] = evaluation.parameter_fingerprint;
    return result;
}

py::dict evaluate_nlp(
    const py::capsule& capsule,
    double temperature_k,
    const std::string& expected_fingerprint,
    const std::array<double, 3>& variables,
    const std::array<double, 3>& multipliers
) {
    const epcsaft_equilibrium::ProviderContext provider(
        checked_sdk(capsule),
        expected_fingerprint
    );
    const epcsaft_equilibrium::NlpEvaluation evaluation =
        epcsaft_equilibrium::evaluate_saturation_nlp(
            provider,
            temperature_k,
            variables,
            multipliers
        );
    py::dict result;
    result["objective"] = evaluation.objective;
    result["objective_gradient"] = evaluation.objective_gradient;
    result["constraints"] = evaluation.constraints;
    result["jacobian"] = evaluation.jacobian;
    result["lagrangian_hessian_lower"] = evaluation.lagrangian_hessian_lower;
    return result;
}

py::dict phase_to_dict(const epcsaft_equilibrium::PhaseEvaluation& phase) {
    py::dict result;
    result["amount_mol"] = phase.amount_mol;
    result["volume_m3"] = phase.volume_m3;
    result["molar_density_mol_m3"] = phase.amount_mol / phase.volume_m3;
    result["pressure_pa"] = phase.provider.pressure_pa;
    result["chemical_potential_over_rt"] = phase.provider.chemical_potential_over_rt;
    return result;
}

py::dict held_stage_iii_phase_to_dict(
    const epcsaft_equilibrium::HeldStageIIIPhaseEvaluation& phase
) {
    const double amount_mol = phase.amounts_mol[0] + phase.amounts_mol[1];
    py::dict result;
    result["amount_mol"] = amount_mol;
    result["mole_fractions"] = std::array<double, 2>{
        phase.amounts_mol[0] / amount_mol,
        phase.amounts_mol[1] / amount_mol,
    };
    result["volume_m3"] = phase.volume_m3;
    result["molar_density_mol_m3"] = amount_mol / phase.volume_m3;
    result["pressure_pa"] = phase.provider.pressure_pa;
    result["chemical_potential_over_rt"] = std::array<double, 2>{
        phase.provider.gradient[0],
        phase.provider.gradient[1],
    };
    return result;
}

py::list attempt_log_to_list(
    const std::vector<epcsaft_equilibrium::SaturationSolveResult::AttemptRecord>& attempts
) {
    py::list result;
    for (const auto& attempt : attempts) {
        py::dict item;
        item["role"] = attempt.role;
        item["initial_guess"] = attempt.initial_guess;
        item["solver_converged"] = attempt.solver_converged;
        item["solver_status"] = attempt.solver_status;
        item["iterations"] = attempt.iterations;
        item["constraint_violation"] = attempt.constraint_violation;
        item["callback_error"] = attempt.callback_error;
        result.append(std::move(item));
    }
    return result;
}

py::list held_stage_iii_attempt_log_to_list(
    const std::vector<epcsaft_equilibrium::HeldStageIIILocalAttempt>& attempts
) {
    py::list result;
    for (const auto& attempt : attempts) {
        py::dict item;
        item["role"] = attempt.role;
        item["initial_guess"] = attempt.initial_guess;
        item["solver_converged"] = attempt.solver_converged;
        item["solver_status"] = attempt.solver_status;
        item["iterations"] = attempt.iterations;
        item["constraint_violation"] = attempt.constraint_violation;
        item["callback_error"] = attempt.callback_error;
        result.append(std::move(item));
    }
    return result;
}

py::dict solve_saturation(
    const py::capsule& capsule,
    double temperature_k,
    const std::string& expected_fingerprint,
    double liquid_density_upper_mol_m3
) {
    const epcsaft_equilibrium::ProviderContext provider(
        checked_sdk(capsule),
        expected_fingerprint
    );
    epcsaft_equilibrium::SaturationSolveResult solve;
    {
        py::gil_scoped_release release;
        solve = epcsaft_equilibrium::solve_local_saturation(
            provider,
            temperature_k,
            liquid_density_upper_mol_m3
        );
    }

    py::dict diagnostics;
    diagnostics["solver_converged"] = solve.solver_converged;
    diagnostics["solver_status"] = solve.solver_status;
    diagnostics["iterations"] = solve.iterations;
    diagnostics["attempts"] = solve.attempts;
    diagnostics["attempt_log"] = attempt_log_to_list(solve.attempt_log);
    diagnostics["solver_lower_bounds"] = solve.solver_lower_bounds;
    diagnostics["solver_upper_bounds"] = solve.solver_upper_bounds;
    diagnostics["solver_constraint_violation"] = solve.solver_constraint_violation;
    diagnostics["numerical_converged"] = solve.numerical_converged;
    diagnostics["confirmation_solves"] = solve.confirmation_solves;
    diagnostics["confirmation_max_relative_difference"] =
        solve.confirmation_max_relative_difference;
    diagnostics["physical_accepted"] = solve.physical_accepted;
    diagnostics["pressure_relative_residual"] = solve.pressure_relative_residual;
    diagnostics["chemical_potential_absolute_residual"] =
        solve.chemical_potential_absolute_residual;
    diagnostics["phase_density_distance"] = solve.phase_density_distance;
    diagnostics["exact_derivatives"] = true;
    diagnostics["globality_certificate"] = false;
    diagnostics["failure_reason"] = solve.failure_reason;

    py::dict result;
    result["accepted"] = solve.accepted;
    result["temperature_k"] = temperature_k;
    result["parameter_fingerprint"] = expected_fingerprint;
    result["diagnostics"] = diagnostics;
    if (solve.accepted) {
        result["saturation_pressure_pa"] = solve.evaluation.saturation_pressure_pa;
        result["vapor"] = phase_to_dict(solve.evaluation.vapor);
        result["liquid"] = phase_to_dict(solve.evaluation.liquid);
    }
    return result;
}

py::dict held2_manufactured(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::array<double, 4>& variables,
    const std::vector<double>& chemical_potentials
) {
    const epcsaft_equilibrium::Held2ManufacturedEvaluation evaluation =
        epcsaft_equilibrium::evaluate_held2_manufactured(
            charges,
            physical_feed,
            variables,
            chemical_potentials
        );
    py::dict metrics;
    metrics["modified_balance_abs"] = evaluation.certificate.modified_balance_abs;
    metrics["ordinary_balance_inf_norm"] =
        evaluation.certificate.ordinary_balance_inf_norm;
    metrics["phase_charge_inf_norm"] = evaluation.certificate.phase_charge_inf_norm;
    metrics["modified_potential_gap"] =
        evaluation.certificate.modified_potential_gap;
    metrics["pressure_stationarity_inf_norm"] =
        evaluation.certificate.pressure_stationarity_inf_norm;
    metrics["reduced_kkt_inf_norm"] = evaluation.certificate.reduced_kkt_inf_norm;
    metrics["enumeration_objective_gap"] =
        evaluation.certificate.enumeration_objective_gap;
    metrics["independent_modified_composition_count"] =
        evaluation.certificate.independent_modified_composition_count;
    py::dict certificate;
    certificate["accepted"] = evaluation.certificate.accepted;
    certificate["independent_evidence"] = evaluation.certificate.independent_evidence;
    certificate["metrics"] = std::move(metrics);

    py::dict result;
    result["formulation_id"] = epcsaft_equilibrium::kHeld2ManufacturedFormulationId;
    result["globality_certificate"] = "not_guaranteed";
    result["eliminated_index"] = evaluation.coordinates.eliminated_index;
    result["dependent_index"] = evaluation.coordinates.dependent_index;
    result["retained_indices"] = evaluation.coordinates.retained_indices;
    result["independent_indices"] = evaluation.coordinates.independent_indices;
    result["modified_factors"] = evaluation.coordinates.modified_factors;
    result["independent_lower_bounds"] =
        evaluation.coordinates.independent_lower_bounds;
    result["independent_upper_bounds"] =
        evaluation.coordinates.independent_upper_bounds;
    result["modified_feed"] = evaluation.modified_feed;
    result["phase_fraction"] = evaluation.phase_fraction;
    result["modified_phases"] = evaluation.modified_phases;
    result["physical_phases"] = evaluation.physical_phases;
    result["phase_charge_residuals"] = evaluation.phase_charge_residuals;
    result["modified_balance"] = evaluation.modified_balance;
    result["ordinary_balance"] = evaluation.ordinary_balance;
    result["transformed_modified_potentials"] =
        evaluation.transformed_modified_potentials;
    result["phase_gibbs_gradients"] = evaluation.phase_gibbs_gradients;
    result["phase_modified_potentials"] = evaluation.phase_modified_potentials;
    result["modified_potential_gap"] = evaluation.modified_potential_gap;
    result["pressure_stationarity_inf_norm"] =
        evaluation.pressure_stationarity_inf_norm;
    result["objective"] = evaluation.objective;
    result["gradient"] = evaluation.gradient;
    result["certificate"] = std::move(certificate);
    return result;
}

py::dict held2_coordinate_evidence(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::vector<double>& chemical_potentials
) {
    const epcsaft_equilibrium::Held2Coordinates coordinates =
        epcsaft_equilibrium::make_held2_coordinates(charges);
    const std::vector<double> modified_feed =
        epcsaft_equilibrium::held2_transform_physical_fractions(
            coordinates,
            physical_feed
        );
    const std::vector<double> lifted_feed =
        epcsaft_equilibrium::held2_lift_modified_fractions(
            coordinates,
            modified_feed
        );
    py::dict result;
    result["eliminated_index"] = coordinates.eliminated_index;
    result["dependent_index"] = coordinates.dependent_index;
    result["retained_indices"] = coordinates.retained_indices;
    result["independent_indices"] = coordinates.independent_indices;
    result["modified_factors"] = coordinates.modified_factors;
    result["independent_lower_bounds"] = coordinates.independent_lower_bounds;
    result["independent_upper_bounds"] = coordinates.independent_upper_bounds;
    result["modified_feed"] = modified_feed;
    result["lifted_feed"] = lifted_feed;
    result["transformed_modified_potentials"] =
        epcsaft_equilibrium::held2_transform_modified_potentials(
            coordinates,
            chemical_potentials
        );
    return result;
}

py::dict held2_stage_ii_simplex_test_adapter(
    const std::vector<double>& independent_lower_bounds,
    const std::vector<double>& independent_upper_bounds,
    double composition_sum_upper,
    const std::vector<double>& values,
    const std::vector<double>& physical_gradient,
    const std::vector<double>& physical_hessian,
    const std::vector<double>& master_multiplier,
    const std::array<double, 2>& phase_bounds,
    const std::string& stage
) {
    if (stage != "stage_ii_simplex_forward"
        && stage != "stage_ii_simplex_inverse"
        && stage != "stage_ii_simplex_chain"
        && stage != "stage_ii_physical_kkt") {
        throw py::value_error("unsupported HELD2 Stage II simplex chart request");
    }
    std::vector<double> chart_values = values;
    if (stage == "stage_ii_simplex_chain") {
        if (values.size() < 2) {
            throw py::value_error("HELD2 Stage II simplex chain is incomplete");
        }
        chart_values.pop_back();
    }
    const auto [physical, chart, jacobian, component_hessians, gradient, hessian,
                stationarity, complementarity, singular] =
        epcsaft_equilibrium::evaluate_held2_stage_ii_simplex_test_adapter(
            independent_lower_bounds,
            independent_upper_bounds,
            composition_sum_upper,
            chart_values,
            physical_gradient,
            physical_hessian,
            master_multiplier,
            phase_bounds,
            stage == "stage_ii_simplex_inverse",
            stage == "stage_ii_physical_kkt"
        );
    py::dict result;
    result["physical"] = physical;
    result["chart"] = chart;
    result["jacobian"] = jacobian;
    result["component_hessians"] = component_hessians;
    result["gradient"] = gradient;
    result["hessian"] = hessian;
    result["stationarity_inf_norm"] = stationarity;
    result["complementarity"] = complementarity;
    result["singular"] = singular;
    return result;
}

py::dict held2_phase_block_evidence(
    const std::vector<double>& charges,
    const std::vector<double>& independent_modified_fractions,
    double log_volume,
    double pressure_over_rt,
    double target_pressure_pa,
    double helmholtz_over_rt,
    const std::vector<double>& gradient,
    const std::vector<double>& hessian,
    double provider_pressure_pa
) {
    const epcsaft_equilibrium::Held2Coordinates coordinates =
        epcsaft_equilibrium::make_held2_coordinates(charges);
    epcsaft_equilibrium::Held2PhysicalPhaseBlock block;
    block.helmholtz_over_rt = helmholtz_over_rt;
    block.gradient = gradient;
    block.hessian = hessian;
    block.pressure_pa = provider_pressure_pa;
    const epcsaft_equilibrium::Held2StateEvaluation evaluation =
        epcsaft_equilibrium::evaluate_held2_phase_block(
            coordinates,
            independent_modified_fractions,
            log_volume,
            pressure_over_rt,
            target_pressure_pa,
            block
        );
    py::dict result;
    result["modified_fractions"] = evaluation.modified_fractions;
    result["physical_amounts"] = evaluation.physical_amounts;
    result["volume"] = evaluation.volume;
    result["objective"] = evaluation.objective;
    result["gradient"] = evaluation.gradient;
    result["hessian"] = evaluation.hessian;
    result["modified_potentials"] = evaluation.modified_potentials;
    result["pressure_stationarity_relative"] =
        evaluation.pressure_stationarity_relative;
    return result;
}

py::dict held2_installed_phase_block(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& independent_modified_fractions,
    double log_volume,
    const std::string& expected_fingerprint
) {
    if (!std::isfinite(temperature_k) || !std::isfinite(pressure_pa)
        || temperature_k <= 0.0 || pressure_pa <= 0.0) {
        throw py::value_error("HELD2 temperature and pressure must be finite and positive");
    }
    const epcsaft_native_sdk_v1& sdk = checked_packing_sdk(capsule);
    Held2ProviderMetadata metadata = held2_provider_metadata(sdk);
    const epcsaft_equilibrium::Held2Coordinates coordinates =
        epcsaft_equilibrium::make_held2_coordinates(metadata.charges);
    const std::vector<double> physical_amounts =
        epcsaft_equilibrium::held2_lift_independent_fractions(
            coordinates,
            independent_modified_fractions
        );
    const double volume_m3 = std::exp(log_volume);
    if (!std::isfinite(volume_m3) || volume_m3 <= 0.0) {
        throw py::value_error("HELD2 phase volume must be finite and positive");
    }
    const epcsaft_equilibrium::ProviderContext provider(sdk, expected_fingerprint);
    const epcsaft_equilibrium::MixturePhaseEvaluation provider_phase =
        provider.evaluate_electrolyte(temperature_k, physical_amounts, volume_m3);
    const epcsaft_equilibrium::PackingFractionEvaluation provider_packing =
        provider.evaluate_packing_fraction(temperature_k, physical_amounts, volume_m3);
    epcsaft_equilibrium::Held2PhysicalPhaseBlock block;
    block.helmholtz_over_rt = provider_phase.value;
    block.gradient = provider_phase.gradient;
    block.hessian = provider_phase.hessian;
    block.pressure_pa = provider_phase.pressure_pa;
    const double pressure_over_rt = pressure_pa / (kGasConstantJPerMolK * temperature_k);
    const epcsaft_equilibrium::Held2StateEvaluation evaluation =
        epcsaft_equilibrium::evaluate_held2_phase_block(
            coordinates,
            independent_modified_fractions,
            log_volume,
            pressure_over_rt,
            pressure_pa,
            block
        );
    const epcsaft_equilibrium::Held2PackingEvaluation packing =
        epcsaft_equilibrium::evaluate_held2_packing_block(
            coordinates,
            independent_modified_fractions,
            log_volume,
            provider_packing.value,
            provider_packing.gradient,
            provider_packing.hessian
        );
    py::dict result;
    result["component_ids"] = std::move(metadata.component_ids);
    result["charges"] = std::move(metadata.charges);
    result["modified_fractions"] = evaluation.modified_fractions;
    result["physical_amounts"] = evaluation.physical_amounts;
    result["volume"] = evaluation.volume;
    result["objective"] = evaluation.objective;
    result["gradient"] = evaluation.gradient;
    result["hessian"] = evaluation.hessian;
    result["packing_fraction"] = packing.value;
    result["packing_gradient"] = packing.gradient;
    result["packing_hessian"] = packing.hessian;
    result["modified_potentials"] = evaluation.modified_potentials;
    result["pressure_stationarity_relative"] =
        evaluation.pressure_stationarity_relative;
    result["provider_pressure_pa"] = provider_phase.pressure_pa;
    result["pressure_over_rt"] = pressure_over_rt;
    result["parameter_fingerprint"] = provider_phase.parameter_fingerprint;
    result["globality_certificate"] = "not_guaranteed";
    return result;
}

py::dict held2_installed_log_packing_phase(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& independent_modified_fractions,
    double log_packing_fraction,
    const std::string& expected_fingerprint,
    const std::string& stage
) {
    if (stage != "log_packing_phase") {
        throw py::value_error("unsupported installed HELD2 coordinate request");
    }
    const epcsaft_native_sdk_v1& sdk = checked_packing_sdk(capsule);
    const Held2ProviderMetadata metadata = held2_provider_metadata(sdk);
    const epcsaft_equilibrium::Held2Coordinates coordinates =
        epcsaft_equilibrium::make_held2_coordinates(metadata.charges);
    const epcsaft_equilibrium::ProviderContext provider(sdk, expected_fingerprint);
    const double pressure_over_rt = pressure_pa / (kGasConstantJPerMolK * temperature_k);
    const epcsaft_equilibrium::Held2StateEvaluator log_volume_evaluator =
        [&provider, coordinates, temperature_k, pressure_pa, pressure_over_rt](
            const std::vector<double>& independent,
            double log_volume
        ) {
            const std::vector<double> physical_amounts =
                epcsaft_equilibrium::held2_lift_independent_fractions(
                    coordinates, independent
                );
            const double volume = std::exp(log_volume);
            const epcsaft_equilibrium::MixturePhaseEvaluation phase =
                provider.evaluate_electrolyte(temperature_k, physical_amounts, volume);
            epcsaft_equilibrium::Held2PhysicalPhaseBlock block;
            block.helmholtz_over_rt = phase.value;
            block.gradient = phase.gradient;
            block.hessian = phase.hessian;
            block.pressure_pa = phase.pressure_pa;
            epcsaft_equilibrium::Held2StateEvaluation state =
                epcsaft_equilibrium::evaluate_held2_phase_block(
                    coordinates,
                    independent,
                    log_volume,
                    pressure_over_rt,
                    pressure_pa,
                    block
                );
            const epcsaft_equilibrium::PackingFractionEvaluation packing =
                provider.evaluate_packing_fraction(
                    temperature_k, physical_amounts, volume
                );
            state.packing = epcsaft_equilibrium::evaluate_held2_packing_block(
                coordinates,
                independent,
                log_volume,
                packing.value,
                packing.gradient,
                packing.hessian
            );
            state.has_packing_evaluation = true;
            return state;
        };
    const std::vector<double> physical_amounts =
        epcsaft_equilibrium::held2_lift_independent_fractions(
            coordinates, independent_modified_fractions
        );
    const std::array<double, 2> volume_bounds = provider.evaluate_molar_volume_bounds(
        temperature_k,
        physical_amounts,
        epcsaft_equilibrium::kHeld2PackingFractionMinimum,
        epcsaft_equilibrium::kHeld2PackingFractionMaximum
    );
    const epcsaft_equilibrium::Held2StateEvaluation evaluation =
        epcsaft_equilibrium::evaluate_held2_log_packing_state(
            log_volume_evaluator,
            independent_modified_fractions,
            log_packing_fraction,
            volume_bounds
        );
    py::dict result;
    result["objective"] = evaluation.objective;
    result["gradient"] = evaluation.gradient;
    result["hessian"] = evaluation.hessian;
    result["volume"] = evaluation.volume;
    result["log_volume"] = std::log(evaluation.volume);
    result["packing_fraction"] = evaluation.packing.value;
    result["log_packing_residual"] =
        std::log(evaluation.packing.value) - log_packing_fraction;
    result["pressure_stationarity_relative"] =
        evaluation.pressure_stationarity_relative;
    result["modified_fractions"] = evaluation.modified_fractions;
    result["physical_amounts"] = evaluation.physical_amounts;
    result["globality_certificate"] = "not_guaranteed";
    return result;
}

py::dict held2_stage_i_to_dict(
    const epcsaft_equilibrium::Held2StageIResult& evaluation,
    const char* profile
) {
    py::list candidates;
    for (const epcsaft_equilibrium::Held2StageICandidate& candidate :
         evaluation.candidates) {
        py::dict item;
        item["modified_fractions"] = candidate.modified_fractions;
        item["volume"] = candidate.volume;
        item["tpd"] = candidate.tpd;
        item["molar_volume_bounds"] = candidate.molar_volume_bounds;
        item["pressure_stationarity_relative"] =
            candidate.pressure_stationarity_relative;
        item["volume_gradient"] = candidate.volume_gradient;
        item["packing_fraction"] = candidate.packing_fraction;
        item["lower_volume_bound_active"] =
            candidate.lower_volume_bound_active;
        item["upper_volume_bound_active"] =
            candidate.upper_volume_bound_active;
        candidates.append(std::move(item));
    }
    py::dict result;
    result["profile"] = profile;
    result["outcome"] = evaluation.outcome;
    result["globality_certificate"] = "not_guaranteed";
    result["reference_scan_interval_count"] =
        evaluation.reference_scan_interval_count;
    result["reference_scan_point_count"] =
        evaluation.reference_scan_point_count;
    result["reference_root_count"] = evaluation.reference_root_count;
    result["reference_stable_root_count"] =
        evaluation.reference_stable_root_count;
    result["reference_evaluation_failure_count"] =
        evaluation.reference_evaluation_failure_count;
    result["reference_refinement_failure_count"] =
        evaluation.reference_refinement_failure_count;
    result["declared_start_count"] = evaluation.declared_start_count;
    result["attempted_start_count"] = evaluation.attempted_start_count;
    result["completed_start_count"] = evaluation.completed_start_count;
    result["failed_start_count"] = evaluation.failed_start_count;
    result["search_completeness"] = evaluation.failed_start_count == 0
        ? "complete"
        : "incomplete";
    result["candidate_domain_evaluation_failure_count"] =
        evaluation.candidate_domain_evaluation_failure_count;
    result["candidate_domain_rejection_count"] =
        evaluation.candidate_domain_rejection_count;
    result["volume_domain_search_complete"] =
        evaluation.volume_domain_search_complete;
    if (evaluation.failed_start_index < 0) {
        result["failed_start_index"] = py::none();
        result["failed_start_solver_status"] = py::none();
        result["failed_start_solver_converged"] = py::none();
        result["failed_start_reason"] = py::none();
        result["failed_start_initial"] = py::none();
    } else {
        result["failed_start_index"] = evaluation.failed_start_index;
        result["failed_start_solver_status"] =
            evaluation.failed_start_solver_status;
        result["failed_start_solver_converged"] =
            evaluation.failed_start_solver_converged;
        result["failed_start_reason"] = evaluation.failed_start_reason;
        result["failed_start_initial"] = evaluation.failed_start_initial;
    }
    result["reference_modified_fractions"] =
        evaluation.reference_modified_fractions;
    if (evaluation.reference_modified_fractions.empty()) {
        result["reference_volume"] = py::none();
    } else {
        result["reference_volume"] = evaluation.reference_volume;
    }
    if (std::isfinite(evaluation.minimum_tpd)) {
        result["minimum_tpd"] = evaluation.minimum_tpd;
    } else {
        result["minimum_tpd"] = py::none();
    }
    py::list reference_roots;
    for (const epcsaft_equilibrium::Held2ReferenceRoot& root :
         evaluation.reference_roots) {
        py::dict item;
        item["log_volume"] = root.log_volume;
        item["volume"] = root.volume;
        item["objective"] = root.objective;
        item["pressure_residual"] = root.pressure_residual;
        item["curvature"] = root.curvature;
        item["mechanically_stable"] = root.mechanically_stable;
        reference_roots.append(std::move(item));
    }
    result["reference_roots"] = std::move(reference_roots);
    result["candidates"] = std::move(candidates);
    return result;
}

py::dict held2_installed_stage_i(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& physical_feed,
    const std::string& expected_fingerprint,
    const std::string& stage
) {
    if (stage != "stage_i" && stage != "stage_i_start_0"
        && stage != "controller") {
        throw py::value_error("unsupported installed HELD2 stage request");
    }
    if (!std::isfinite(temperature_k) || !std::isfinite(pressure_pa)
        || temperature_k <= 0.0 || pressure_pa <= 0.0) {
        throw py::value_error("HELD2 temperature and pressure must be finite and positive");
    }
    const epcsaft_native_sdk_v1& sdk = checked_packing_sdk(capsule);
    Held2ProviderMetadata metadata = held2_provider_metadata(sdk);
    const epcsaft_equilibrium::Held2Coordinates coordinates =
        epcsaft_equilibrium::make_held2_coordinates(metadata.charges);
    static_cast<void>(epcsaft_equilibrium::held2_transform_physical_fractions(
        coordinates,
        physical_feed
    ));
    const epcsaft_equilibrium::ProviderContext provider(sdk, expected_fingerprint);
    const std::array<double, 2> molar_volume_bounds =
        provider.evaluate_molar_volume_bounds(
            temperature_k,
            physical_feed,
            epcsaft_equilibrium::kHeld2PackingFractionMinimum,
            epcsaft_equilibrium::kHeld2PackingFractionMaximum
        );
    const double pressure_over_rt =
        pressure_pa / (kGasConstantJPerMolK * temperature_k);
    const epcsaft_equilibrium::Held2StateEvaluator evaluator =
        [&provider, coordinates, temperature_k, pressure_pa, pressure_over_rt](
            const std::vector<double>& independent_modified_fractions,
            double log_volume
        ) {
            const std::vector<double> physical_amounts =
                epcsaft_equilibrium::held2_lift_independent_fractions(
                    coordinates,
                    independent_modified_fractions
                );
            const epcsaft_equilibrium::MixturePhaseEvaluation phase =
                provider.evaluate_electrolyte(
                    temperature_k,
                    physical_amounts,
                    std::exp(log_volume)
                );
            epcsaft_equilibrium::Held2PhysicalPhaseBlock block;
            block.helmholtz_over_rt = phase.value;
            block.gradient = phase.gradient;
            block.hessian = phase.hessian;
            block.pressure_pa = phase.pressure_pa;
            epcsaft_equilibrium::Held2StateEvaluation state =
                epcsaft_equilibrium::evaluate_held2_phase_block(
                coordinates,
                independent_modified_fractions,
                log_volume,
                pressure_over_rt,
                pressure_pa,
                block
            );
            const epcsaft_equilibrium::PackingFractionEvaluation provider_packing =
                provider.evaluate_packing_fraction(
                    temperature_k,
                    physical_amounts,
                    std::exp(log_volume)
                );
            state.packing = epcsaft_equilibrium::evaluate_held2_packing_block(
                coordinates,
                independent_modified_fractions,
                log_volume,
                provider_packing.value,
                provider_packing.gradient,
                provider_packing.hessian
            );
            state.has_packing_evaluation = true;
            return state;
        };
    const epcsaft_equilibrium::Held2VolumeBoundsEvaluator volume_bounds_evaluator =
        [&provider, temperature_k](const std::vector<double>& physical_amounts) {
            return provider.evaluate_molar_volume_bounds(
                temperature_k,
                physical_amounts,
                epcsaft_equilibrium::kHeld2PackingFractionMinimum,
                epcsaft_equilibrium::kHeld2PackingFractionMaximum
            );
        };
    const epcsaft_equilibrium::Held2StateEvaluator search_evaluator =
        [&evaluator, &volume_bounds_evaluator, coordinates](
            const std::vector<double>& independent_modified_fractions,
            double log_packing_fraction
        ) {
            const std::vector<double> physical_amounts =
                epcsaft_equilibrium::held2_lift_independent_fractions(
                    coordinates, independent_modified_fractions
                );
            return epcsaft_equilibrium::evaluate_held2_log_packing_state(
                evaluator,
                independent_modified_fractions,
                log_packing_fraction,
                volume_bounds_evaluator(physical_amounts)
            );
        };
    epcsaft_equilibrium::Held2StageIResult evaluation;
    {
        py::gil_scoped_release release;
        const bool q_discriminator = stage == "stage_i_start_0";
        evaluation = epcsaft_equilibrium::solve_held2_stage_i(
            coordinates,
            physical_feed,
            evaluator,
            q_discriminator ? search_evaluator : evaluator,
            molar_volume_bounds,
            volume_bounds_evaluator,
            q_discriminator,
            true,
            q_discriminator ? 1 : -1
        );
    }
    py::dict result = held2_stage_i_to_dict(
        evaluation,
        "perdomo-held2-stage-i-installed-v1"
    );
    result["component_ids"] = std::move(metadata.component_ids);
    result["charges"] = std::move(metadata.charges);
    result["parameter_fingerprint"] = expected_fingerprint;
    result["packing_fraction_bounds"] = std::array<double, 2>{
        epcsaft_equilibrium::kHeld2PackingFractionMinimum,
        epcsaft_equilibrium::kHeld2PackingFractionMaximum,
    };
    result["molar_volume_bounds"] = molar_volume_bounds;
    result["log_molar_volume_bounds"] = std::array<double, 2>{
        std::log(molar_volume_bounds[0]),
        std::log(molar_volume_bounds[1]),
    };
    if (stage == "controller" && evaluation.outcome == "negative_tpd") {
        std::vector<double> feed_independent;
        const std::vector<double> modified_feed =
            epcsaft_equilibrium::held2_transform_physical_fractions(
                coordinates, physical_feed
            );
        for (std::size_t component : coordinates.independent_indices) {
            const auto retained = std::find(
                coordinates.retained_indices.begin(),
                coordinates.retained_indices.end(),
                component
            );
            feed_independent.push_back(modified_feed[static_cast<std::size_t>(
                retained - coordinates.retained_indices.begin()
            )]);
        }
        epcsaft_equilibrium::Held2StateEvaluation reference = evaluator(
            feed_independent, std::log(evaluation.reference_volume)
        );
        reference = search_evaluator(
            feed_independent, std::log(reference.packing.value)
        );
        epcsaft_equilibrium::Held2StageIIResult stage_ii;
        {
            py::gil_scoped_release release;
            stage_ii = epcsaft_equilibrium::solve_held2_stage_ii(
                coordinates,
                physical_feed,
                search_evaluator,
                evaluator,
                reference,
                evaluation.candidates
            );
        }
        py::dict stage_ii_payload;
        stage_ii_payload["outcome"] = stage_ii.outcome;
        stage_ii_payload["major_iterations"] = stage_ii.major_iterations;
        stage_ii_payload["lower_starts_per_iteration"] =
            stage_ii.lower_starts_per_iteration;
        stage_ii_payload["lower_attempted_start_count"] =
            stage_ii.lower_attempted_start_count;
        stage_ii_payload["lower_completed_start_count"] =
            stage_ii.lower_completed_start_count;
        stage_ii_payload["lower_failed_start_count"] =
            stage_ii.lower_failed_start_count;
        stage_ii_payload["certified_improving_cut_count"] =
            stage_ii.certified_improving_cut_count;
        stage_ii_payload["search_completeness"] =
            stage_ii.search_completeness;
        stage_ii_payload["final_lower_search_complete"] =
            stage_ii.final_lower_search_complete;
        if (stage_ii.first_failed_lower_start_index < 0) {
            stage_ii_payload["first_failed_lower_start_index"] = py::none();
            stage_ii_payload["first_failed_lower_solver_status"] = py::none();
            stage_ii_payload["first_failed_lower_reason"] = py::none();
            stage_ii_payload["first_failed_lower_initial"] = py::none();
        } else {
            stage_ii_payload["first_failed_lower_start_index"] =
                stage_ii.first_failed_lower_start_index;
            stage_ii_payload["first_failed_lower_solver_status"] =
                stage_ii.first_failed_lower_solver_status;
            stage_ii_payload["first_failed_lower_reason"] =
                stage_ii.first_failed_lower_reason;
            stage_ii_payload["first_failed_lower_initial"] =
                stage_ii.first_failed_lower_initial;
        }
        stage_ii_payload["cut_count"] = stage_ii.cut_count;
        py::list bound_history;
        for (const auto& bound : stage_ii.bound_history) {
            py::dict item;
            item["lower_bound"] = bound.lower_bound;
            item["upper_bound"] = bound.upper_bound;
            item["multiplier"] = bound.multiplier;
            item["cut_count"] = bound.cut_count;
            bound_history.append(std::move(item));
        }
        stage_ii_payload["bound_history"] = std::move(bound_history);
        py::list attempt_log;
        for (const auto& attempt : stage_ii.attempt_log) {
            py::dict item;
            item["start_index"] = attempt.start_index;
            item["solver_status"] = attempt.solver_status;
            item["iterations"] = attempt.iterations;
            item["solver_converged"] = attempt.solver_converged;
            item["numerical_certified"] = attempt.numerical_certified;
            item["provider_terminal_valid"] = attempt.provider_terminal_valid;
            item["improving"] = attempt.improving;
            item["callback_error"] = attempt.callback_error;
            item["lower_value"] = attempt.lower_value;
            item["projected_kkt_inf_norm"] =
                attempt.projected_kkt_inf_norm;
            item["constraint_violation"] = attempt.constraint_violation;
            item["complementarity"] = attempt.complementarity;
            item["final_step_norm"] = attempt.final_step_norm;
            attempt_log.append(std::move(item));
        }
        stage_ii_payload["attempt_log"] = std::move(attempt_log);
        stage_ii_payload["candidate_count"] = stage_ii.candidates.size();
        result["stage_ii"] = std::move(stage_ii_payload);
        if (stage_ii.outcome == "candidate_set") {
            epcsaft_equilibrium::Held2StageIIIResult stage_iii;
            {
                py::gil_scoped_release release;
                stage_iii = epcsaft_equilibrium::solve_held2_stage_iii(
                    coordinates,
                    physical_feed,
                    stage_ii.candidates,
                    search_evaluator,
                    {
                        std::log(epcsaft_equilibrium::kHeld2PackingFractionMinimum),
                        std::log(epcsaft_equilibrium::kHeld2PackingFractionMaximum),
                    }
                );
            }
            py::dict stage_iii_payload;
            stage_iii_payload["solver_status"] = stage_iii.solver_status;
            stage_iii_payload["numerical_status"] = stage_iii.numerical_status;
            stage_iii_payload["physical_status"] = stage_iii.physical_status;
            stage_iii_payload["feedback"] = stage_iii.feedback;
            stage_iii_payload["failure_reason"] = stage_iii.failure_reason;
            stage_iii_payload["input_candidate_count"] =
                stage_iii.input_candidate_count;
            stage_iii_payload["modified_balance_inf_norm"] =
                stage_iii.modified_balance_inf_norm;
            stage_iii_payload["ordinary_balance_inf_norm"] =
                stage_iii.ordinary_balance_inf_norm;
            stage_iii_payload["phase_charge_inf_norm"] =
                stage_iii.phase_charge_inf_norm;
            stage_iii_payload["pressure_stationarity_inf_norm"] =
                stage_iii.pressure_stationarity_inf_norm;
            stage_iii_payload["modified_potential_mixed_gap"] =
                stage_iii.modified_potential_mixed_gap;
            stage_iii_payload["minimum_phase_distance"] =
                stage_iii.minimum_phase_distance;
            stage_iii_payload["kkt_stationarity_inf_norm"] =
                stage_iii.kkt_stationarity_inf_norm;
            py::list phases;
            for (const auto& phase : stage_iii.phases) {
                py::dict item;
                item["phase_fraction"] = phase.phase_fraction;
                item["modified_fractions"] = phase.modified_fractions;
                item["physical_fractions"] = phase.physical_fractions;
                item["volume"] = phase.volume;
                phases.append(std::move(item));
            }
            stage_iii_payload["phases"] = std::move(phases);
            result["stage_iii"] = std::move(stage_iii_payload);
        }
    }
    return result;
}

py::dict held2_manufactured_stage_i(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::string& stage
) {
    if (stage == "stage_ii") {
        const epcsaft_equilibrium::Held2StageIIResult evaluation =
            epcsaft_equilibrium::solve_held2_manufactured_stage_ii(
                charges,
                physical_feed
            );
        py::list history;
        for (const epcsaft_equilibrium::Held2StageIIBound& bound :
             evaluation.bound_history) {
            py::dict item;
            item["lower_bound"] = bound.lower_bound;
            item["upper_bound"] = bound.upper_bound;
            item["multiplier"] = bound.multiplier;
            item["cut_count"] = bound.cut_count;
            history.append(std::move(item));
        }
        py::list candidates;
        for (const epcsaft_equilibrium::Held2StageIICandidate& candidate :
             evaluation.candidates) {
            py::dict item;
            item["modified_fractions"] = candidate.modified_fractions;
            item["volume"] = candidate.volume;
            item["lower_gap"] = candidate.lower_gap;
            candidates.append(std::move(item));
        }
        py::dict result;
        result["profile"] = "perdomo-held2-stage-ii-manufactured-v1";
        result["outcome"] = evaluation.outcome;
        result["globality_certificate"] = "not_guaranteed";
        result["major_iterations"] = evaluation.major_iterations;
        result["lower_starts_per_iteration"] =
            evaluation.lower_starts_per_iteration;
        result["cut_count"] = evaluation.cut_count;
        result["bound_history"] = std::move(history);
        result["candidates"] = std::move(candidates);
        return result;
    }
    if (stage != "stage_i") {
        throw py::value_error("unsupported manufactured HELD2 stage request");
    }
    const epcsaft_equilibrium::Held2StageIResult evaluation =
        epcsaft_equilibrium::solve_held2_manufactured_stage_i(
            charges,
            physical_feed
        );
    return held2_stage_i_to_dict(
        evaluation,
        "perdomo-held2-stage-i-manufactured-v1"
    );
}

py::dict held2_stage_ii_precedence(
    double upper_bound,
    double best_certified_value,
    int certified_start_count,
    int declared_start_count,
    const std::string& stage
) {
    if (stage != "stage_ii_precedence") {
        throw py::value_error("unsupported HELD2 Stage II precedence request");
    }
    const auto decision = epcsaft_equilibrium::decide_held2_stage_ii_lower(
        upper_bound,
        best_certified_value,
        certified_start_count,
        declared_start_count
    );
    py::dict result;
    result["decision"] = decision.decision;
    result["search_completeness"] = decision.search_completeness;
    result["lower_bound_certified"] = decision.lower_bound_certified;
    return result;
}

py::dict held2_stage_ii_step6(
    const std::vector<double>& charges,
    const std::vector<double>& feed_independent_modified_fractions,
    double upper_bound,
    const std::vector<double>& multipliers,
    const py::sequence& raw_cuts,
    const std::string& stage
) {
    if (stage != "stage_ii_step6") {
        throw py::value_error("unsupported HELD2 Stage II Step 6 request");
    }
    const epcsaft_equilibrium::Held2Coordinates coordinates =
        epcsaft_equilibrium::make_held2_coordinates(charges);
    std::vector<epcsaft_equilibrium::Held2StateEvaluation> states;
    std::vector<std::vector<double>> independent_modified_fractions;
    std::vector<std::vector<double>> fixed_volume_composition_gradients;
    std::vector<double> phase_coordinates;
    states.reserve(raw_cuts.size());
    independent_modified_fractions.reserve(raw_cuts.size());
    fixed_volume_composition_gradients.reserve(raw_cuts.size());
    phase_coordinates.reserve(raw_cuts.size());
    for (const py::handle raw : raw_cuts) {
        const py::sequence cut = py::reinterpret_borrow<py::sequence>(raw);
        if (cut.size() != 6) {
            throw py::value_error("HELD2 Stage II Step 6 cut must contain six fields");
        }
        const std::vector<double> independent =
            cut[0].cast<std::vector<double>>();
        const double packing_fraction = cut[1].cast<double>();
        if (!(packing_fraction > 0.0) || !std::isfinite(packing_fraction)) {
            throw py::value_error("HELD2 Stage II Step 6 packing fraction is invalid");
        }
        epcsaft_equilibrium::Held2StateEvaluation state;
        state.physical_amounts =
            epcsaft_equilibrium::held2_lift_independent_fractions(
                coordinates,
                independent
            );
        state.modified_fractions =
            epcsaft_equilibrium::held2_transform_physical_fractions(
                coordinates,
                state.physical_amounts
            );
        state.volume = cut[2].cast<double>();
        state.objective = cut[3].cast<double>();
        state.gradient = cut[4].cast<std::vector<double>>();
        state.has_packing_evaluation = true;
        state.packing.value = packing_fraction;
        states.push_back(std::move(state));
        independent_modified_fractions.push_back(independent);
        fixed_volume_composition_gradients.push_back(
            cut[5].cast<std::vector<double>>()
        );
        phase_coordinates.push_back(std::log(packing_fraction));
    }
    const std::vector<epcsaft_equilibrium::Held2StageIICandidate> candidates =
        epcsaft_equilibrium::evaluate_held2_stage_ii_step6_test_adapter(
            coordinates,
            feed_independent_modified_fractions,
            upper_bound,
            multipliers,
            states,
            independent_modified_fractions,
            fixed_volume_composition_gradients,
            phase_coordinates
        );
    std::vector<std::vector<double>> independent_candidates;
    independent_candidates.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        independent_candidates.push_back(
            candidate.independent_modified_fractions
        );
    }
    py::dict result;
    result["candidate_count"] = candidates.size();
    result["independent_modified_fractions"] =
        std::move(independent_candidates);
    return result;
}

py::dict held2_manufactured_stage_iii_derivatives(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::vector<std::array<double, 2>>& candidates,
    const std::vector<double>& variables,
    const std::vector<double>& equality_multipliers,
    const std::string& stage
) {
    if (stage != "stage_iii_derivatives") {
        throw py::value_error("unsupported manufactured HELD2 derivative request");
    }
    const epcsaft_equilibrium::Held2StageIIINlpEvaluation evaluation =
        epcsaft_equilibrium::evaluate_held2_manufactured_stage_iii_nlp(
            charges,
            physical_feed,
            candidates,
            variables,
            equality_multipliers
        );
    py::dict result;
    result["objective"] = evaluation.objective;
    result["objective_gradient"] = evaluation.objective_gradient;
    result["constraints"] = evaluation.constraints;
    result["constraint_jacobian"] = evaluation.constraint_jacobian;
    result["lagrangian_gradient"] = evaluation.lagrangian_gradient;
    result["lagrangian_hessian"] = evaluation.lagrangian_hessian;
    return result;
}

py::dict held2_general_stage_iii_test_seam(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::vector<std::vector<double>>& raw_candidates,
    const std::vector<double>& variables,
    const std::vector<double>& multipliers,
    const std::string& stage
) {
    if (stage != "stage_iii_general_schema"
        && stage != "stage_iii_invalid_trial") {
        throw py::value_error("unsupported general HELD2 Stage III test request");
    }
    const epcsaft_equilibrium::Held2Coordinates coordinates =
        epcsaft_equilibrium::make_held2_coordinates(charges);
    const std::size_t dimension = coordinates.independent_indices.size();
    std::vector<epcsaft_equilibrium::Held2StageIICandidate> candidates;
    candidates.reserve(raw_candidates.size());
    for (const std::vector<double>& raw : raw_candidates) {
        if (raw.size() != dimension + 1) {
            throw py::value_error("general HELD2 Stage III candidate dimension changed");
        }
        epcsaft_equilibrium::Held2StageIICandidate candidate;
        candidate.independent_modified_fractions.assign(raw.begin(), raw.end() - 1);
        candidate.phase_coordinate = raw.back();
        candidates.push_back(std::move(candidate));
    }
    int scientific_evaluator_call_count = 0;
    const epcsaft_equilibrium::Held2StateEvaluator evaluator =
        [&scientific_evaluator_call_count, coordinates](
            const std::vector<double>& independent,
            double phase_coordinate
        ) {
            ++scientific_evaluator_call_count;
            const std::vector<double> physical =
                epcsaft_equilibrium::held2_lift_independent_fractions(
                    coordinates, independent
                );
            epcsaft_equilibrium::Held2StateEvaluation result;
            result.physical_amounts = physical;
            result.modified_fractions =
                epcsaft_equilibrium::held2_transform_physical_fractions(
                    coordinates, physical
                );
            result.volume = std::exp(phase_coordinate);
            result.gradient.reserve(independent.size() + 1);
            result.objective = 0.5 * phase_coordinate * phase_coordinate;
            for (double value : independent) {
                result.objective += 0.5 * value * value;
                result.gradient.push_back(value);
            }
            result.gradient.push_back(phase_coordinate);
            const std::size_t state_dimension = independent.size() + 1;
            result.hessian.assign(state_dimension * state_dimension, 0.0);
            for (std::size_t index = 0; index < state_dimension; ++index) {
                result.hessian[index * state_dimension + index] = 1.0;
            }
            result.modified_potentials.assign(
                coordinates.retained_indices.size(), 0.0
            );
            return result;
        };
    const std::array<double, 2> phase_coordinate_bounds{
        std::log(epcsaft_equilibrium::kHeld2PackingFractionMinimum),
        std::log(epcsaft_equilibrium::kHeld2PackingFractionMaximum),
    };
    py::dict result;
    if (stage == "stage_iii_invalid_trial") {
        const auto [accepted, domain_rejections, last_domain_rejection,
                    fatal_callback_error] =
            epcsaft_equilibrium::probe_held2_stage_iii_objective_trial(
                coordinates,
                physical_feed,
                candidates,
                evaluator,
                phase_coordinate_bounds,
                variables
            );
        result["objective_accepted"] = accepted;
        result["scientific_evaluator_call_count"] =
            scientific_evaluator_call_count;
        result["recoverable_domain_rejection_count"] = domain_rejections;
        result["last_domain_rejection"] = last_domain_rejection;
        result["fatal_callback_error"] = fatal_callback_error;
        return result;
    }

    const auto [variable_count, constraint_count, jacobian_nonzero_count,
                constraint_lower_bounds, constraint_upper_bounds,
                jacobian_rows, jacobian_columns, jacobian_values] =
        epcsaft_equilibrium::inspect_held2_stage_iii_tnlp(
            coordinates,
            physical_feed,
            candidates,
            phase_coordinate_bounds,
            variables
        );
    const epcsaft_equilibrium::Held2StageIIINlpEvaluation evaluation =
        epcsaft_equilibrium::evaluate_held2_stage_iii_nlp(
            coordinates,
            physical_feed,
            evaluator,
            candidates.size(),
            variables,
            multipliers
        );
    result["variable_count"] = variable_count;
    result["constraint_count"] = constraint_count;
    result["jacobian_nonzero_count"] = jacobian_nonzero_count;
    result["constraint_lower_bounds"] = constraint_lower_bounds;
    result["constraint_upper_bounds"] = constraint_upper_bounds;
    result["jacobian_rows"] = jacobian_rows;
    result["jacobian_columns"] = jacobian_columns;
    result["jacobian_values"] = jacobian_values;
    result["objective"] = evaluation.objective;
    result["objective_gradient"] = evaluation.objective_gradient;
    result["constraints"] = evaluation.constraints;
    result["constraint_jacobian"] = evaluation.constraint_jacobian;
    result["lagrangian_gradient"] = evaluation.lagrangian_gradient;
    result["lagrangian_hessian"] = evaluation.lagrangian_hessian;
    return result;
}

py::dict held2_manufactured_stage_iii(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::vector<std::array<double, 2>>& candidates,
    const std::string& stage
) {
    if (stage != "stage_iii") {
        throw py::value_error("unsupported manufactured HELD2 Stage III request");
    }
    const epcsaft_equilibrium::Held2StageIIIResult evaluation =
        epcsaft_equilibrium::solve_held2_manufactured_stage_iii(
            charges,
            physical_feed,
            candidates
        );
    py::list phases;
    for (const epcsaft_equilibrium::Held2StageIIIPhase& phase : evaluation.phases) {
        py::dict item;
        item["phase_fraction"] = phase.phase_fraction;
        item["modified_fractions"] = phase.modified_fractions;
        item["physical_fractions"] = phase.physical_fractions;
        item["volume"] = phase.volume;
        phases.append(std::move(item));
    }
    py::dict result;
    result["profile"] = "perdomo-held2-stage-iii-manufactured-v1";
    result["solver_status"] = evaluation.solver_status;
    result["numerical_status"] = evaluation.numerical_status;
    result["physical_status"] = evaluation.physical_status;
    result["globality_certificate"] = "not_guaranteed";
    result["feedback"] = evaluation.feedback;
    result["failure_reason"] = evaluation.failure_reason;
    result["trace_refinement_status"] = evaluation.trace_refinement_status;
    result["input_candidate_count"] = evaluation.input_candidate_count;
    result["retired_duplicate_count"] = evaluation.retired_duplicate_count;
    result["trace_component_count"] = evaluation.trace_component_count;
    result["certified_modified_potential_count"] =
        evaluation.certified_modified_potential_count;
    result["objective"] = evaluation.objective;
    result["modified_balance_inf_norm"] = evaluation.modified_balance_inf_norm;
    result["ordinary_balance_inf_norm"] = evaluation.ordinary_balance_inf_norm;
    result["phase_charge_inf_norm"] = evaluation.phase_charge_inf_norm;
    result["pressure_stationarity_inf_norm"] =
        evaluation.pressure_stationarity_inf_norm;
    result["modified_potential_mixed_gap"] =
        evaluation.modified_potential_mixed_gap;
    result["minimum_phase_distance"] = evaluation.minimum_phase_distance;
    result["kkt_stationarity_inf_norm"] = evaluation.kkt_stationarity_inf_norm;
    result["enumeration_objective_gap"] = evaluation.enumeration_objective_gap;
    result["phases"] = std::move(phases);
    return result;
}

}  // namespace

PYBIND11_MODULE(_equilibrium, module) {
    module.doc() = "Native local equilibrium formulation over epcsaft.native_sdk.v1";
    module.def("sdk_info", &sdk_info, py::arg("capsule"));
    module.def(
        "evaluate_mixture_phase",
        &evaluate_mixture_phase,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("amounts_mol"),
        py::arg("volume_m3"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "_held_evaluate_state",
        &held_evaluate_state,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("x_methane"),
        py::arg("log_volume"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "_held_evaluate_tpd",
        &held_evaluate_tpd,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("reference_x_methane"),
        py::arg("reference_log_volume"),
        py::arg("x_methane"),
        py::arg("log_volume"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "_held_evaluate_tunneling",
        &held_evaluate_tunneling,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("reference_x_methane"),
        py::arg("reference_log_volume"),
        py::arg("minimum_x_methane"),
        py::arg("minimum_log_volume"),
        py::arg("x_methane"),
        py::arg("log_volume"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "_held_stage_i",
        &held_stage_i,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("feed_x_methane"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "_held_outer_envelope",
        &held_outer_envelope,
        py::arg("cuts")
    );
    module.def(
        "_held_stage_ii_initial_cuts",
        &held_stage_ii_initial_cuts,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("feed_x_methane"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "_held_evaluate_lower",
        &held_evaluate_lower,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("feed_x_methane"),
        py::arg("multiplier"),
        py::arg("x_methane"),
        py::arg("log_volume"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "_held_stage_ii_lower_search",
        &held_stage_ii_lower_search,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("feed_x_methane"),
        py::arg("multiplier"),
        py::arg("upper_bound"),
        py::arg("previous_states"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "_held_stage_ii_candidates",
        &held_stage_ii_candidates,
        py::arg("feed_x_methane"),
        py::arg("upper_bound"),
        py::arg("multiplier"),
        py::arg("points")
    );
    module.def(
        "_held_stage_ii_budget_status",
        &epcsaft_equilibrium::held_stage_ii_budget_status,
        py::arg("major_iterations"),
        py::arg("lower_starts"),
        py::arg("lower_satisfied")
    );
    module.def(
        "_held_stage_ii",
        &held_stage_ii,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("feed_x_methane"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "_held_stage_iii_initialize",
        &held_stage_iii_initialize,
        py::arg("feed_x_methane"),
        py::arg("candidates")
    );
    module.def(
        "_held_stage_iii_classify",
        &held_stage_iii_classify,
        py::arg("evidence")
    );
    module.def(
        "_held_stage_iii_mu_difference",
        &held_stage_iii_mu_difference,
        py::arg("first"),
        py::arg("second")
    );
    module.def(
        "_held_evaluate_stage_iii",
        &held_evaluate_stage_iii,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("feed_x_methane"),
        py::arg("variables"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "_held_stage_iii",
        &held_stage_iii,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("feed_x_methane"),
        py::arg("upper_bound"),
        py::arg("candidates"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "_held2_adapter",
        &held2_coordinate_evidence,
        py::arg("charges"),
        py::arg("physical_feed"),
        py::arg("chemical_potentials")
    );
    module.def(
        "_held2_adapter",
        &held2_manufactured,
        py::arg("charges"),
        py::arg("physical_feed"),
        py::arg("variables"),
        py::arg("chemical_potentials")
    );
    module.def(
        "_held2_adapter",
        &held2_phase_block_evidence,
        py::arg("charges"),
        py::arg("independent_modified_fractions"),
        py::arg("log_volume"),
        py::arg("pressure_over_rt"),
        py::arg("target_pressure_pa"),
        py::arg("helmholtz_over_rt"),
        py::arg("gradient"),
        py::arg("hessian"),
        py::arg("provider_pressure_pa")
    );
    module.def(
        "_held2_adapter",
        &held2_installed_phase_block,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("independent_modified_fractions"),
        py::arg("log_volume"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "_held2_adapter",
        &held2_installed_log_packing_phase,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("independent_modified_fractions"),
        py::arg("log_packing_fraction"),
        py::arg("expected_fingerprint"),
        py::arg("stage")
    );
    module.def(
        "_held2_adapter",
        &held2_installed_stage_i,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("physical_feed"),
        py::arg("expected_fingerprint"),
        py::arg("stage")
    );
    module.def(
        "_held2_adapter",
        &held2_manufactured_stage_i,
        py::arg("charges"),
        py::arg("physical_feed"),
        py::arg("stage")
    );
    module.def(
        "_held2_adapter",
        &held2_stage_ii_simplex_test_adapter,
        py::arg("independent_lower_bounds"),
        py::arg("independent_upper_bounds"),
        py::arg("composition_sum_upper"),
        py::arg("values"),
        py::arg("physical_gradient"),
        py::arg("physical_hessian"),
        py::arg("master_multiplier"),
        py::arg("phase_bounds"),
        py::arg("stage")
    );
    module.def(
        "_held2_adapter",
        &held2_stage_ii_precedence,
        py::arg("upper_bound"),
        py::arg("best_certified_value"),
        py::arg("certified_start_count"),
        py::arg("declared_start_count"),
        py::arg("stage")
    );
    module.def(
        "_held2_adapter",
        &held2_stage_ii_step6,
        py::arg("charges"),
        py::arg("feed_independent_modified_fractions"),
        py::arg("upper_bound"),
        py::arg("multipliers"),
        py::arg("cuts"),
        py::arg("stage")
    );
    module.def(
        "_held2_adapter",
        &held2_manufactured_stage_iii,
        py::arg("charges"),
        py::arg("physical_feed"),
        py::arg("candidates"),
        py::arg("stage")
    );
    module.def(
        "_held2_adapter",
        &held2_manufactured_stage_iii_derivatives,
        py::arg("charges"),
        py::arg("physical_feed"),
        py::arg("candidates"),
        py::arg("variables"),
        py::arg("equality_multipliers"),
        py::arg("stage")
    );
    module.def(
        "_held2_adapter",
        &held2_general_stage_iii_test_seam,
        py::arg("charges"),
        py::arg("physical_feed"),
        py::arg("candidates"),
        py::arg("variables"),
        py::arg("multipliers"),
        py::arg("stage")
    );
    module.def(
        "_solve_tp_flash",
        &solve_tp_flash,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("overall_mole_fractions")
    );
    module.def(
        "evaluate_phase",
        &evaluate_phase,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("amount_mol"),
        py::arg("volume_m3"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "evaluate_nlp",
        &evaluate_nlp,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("expected_fingerprint"),
        py::arg("variables"),
        py::arg("multipliers")
    );
    module.def(
        "solve_saturation",
        &solve_saturation,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("expected_fingerprint"),
        py::arg("liquid_density_upper_mol_m3")
    );
}
