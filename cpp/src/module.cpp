#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <epcsaft/native_sdk_v1.h>

#include "held.hpp"
#include "held2.hpp"
#include "held2_stage_i_direct.hpp"
#include "held2_stage_ii_basin.hpp"
#include "held2_stage_ii_upper.hpp"
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
constexpr double kGasConstantJPerMolK = 8.31446261815324;

py::dict held2_stage_ii_to_dict(
    const epcsaft_equilibrium::Held2StageIIResult& evaluation
);
py::dict held2_stage_iii_to_dict(
    const epcsaft_equilibrium::Held2StageIIIResult& evaluation
);

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

const epcsaft_native_sdk_v1& checked_molar_volume_sdk(
    const py::capsule& capsule
) {
    const epcsaft_native_sdk_v1& sdk = checked_electrolyte_sdk(capsule);
    if (sdk.table_size < kMolarVolumeSdkTableSize
        || sdk.evaluate_molar_volume_bounds == nullptr) {
        throw py::value_error(
            "provider capsule is missing the molar-volume domain contract"
        );
    }
    return sdk;
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

py::dict held2_state_evaluation_to_dict(
    const epcsaft_equilibrium::Held2StateEvaluation& evaluation
) {
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
    result["pressure_stationarity_derivative_log_volume"] =
        evaluation.pressure_stationarity_derivative_log_volume;
    return result;
}

py::dict held2_manufactured_search_objective_evidence(
    const std::vector<double>& charges,
    const std::vector<double>& variables,
    const std::vector<double>& reference_variables,
    bool use_tpd,
    const std::string& stage
) {
    if (stage != "search_objective") {
        throw py::value_error("unsupported HELD2 search-objective adapter stage");
    }
    const epcsaft_equilibrium::Held2Coordinates coordinates =
        epcsaft_equilibrium::make_held2_coordinates(charges);
    return held2_state_evaluation_to_dict(
        epcsaft_equilibrium::evaluate_held2_manufactured_search_objective(
            coordinates,
            variables,
            reference_variables,
            use_tpd
        )
    );
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
    return held2_state_evaluation_to_dict(evaluation);
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
    const epcsaft_native_sdk_v1& sdk = checked_electrolyte_sdk(capsule);
    std::vector<double> charges;
    std::vector<std::string> component_ids;
    charges.reserve(sdk.component_count);
    component_ids.reserve(sdk.component_count);
    for (std::size_t component = 0; component < sdk.component_count; ++component) {
        if (sdk.component_ids[component] == nullptr
            || sdk.component_ids[component][0] == '\0') {
            throw py::value_error("provider electrolyte component ID must not be empty");
        }
        charges.push_back(static_cast<double>(sdk.component_charges[component]));
        component_ids.emplace_back(sdk.component_ids[component]);
    }
    const epcsaft_equilibrium::Held2Coordinates coordinates =
        epcsaft_equilibrium::make_held2_coordinates(charges);
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
    py::dict result;
    result["component_ids"] = std::move(component_ids);
    result["charges"] = std::move(charges);
    result["modified_fractions"] = evaluation.modified_fractions;
    result["physical_amounts"] = evaluation.physical_amounts;
    result["volume"] = evaluation.volume;
    result["objective"] = evaluation.objective;
    result["gradient"] = evaluation.gradient;
    result["hessian"] = evaluation.hessian;
    result["modified_potentials"] = evaluation.modified_potentials;
    result["pressure_stationarity_relative"] =
        evaluation.pressure_stationarity_relative;
    result["pressure_stationarity_derivative_log_volume"] =
        evaluation.pressure_stationarity_derivative_log_volume;
    result["provider_pressure_pa"] = provider_phase.pressure_pa;
    result["pressure_over_rt"] = pressure_over_rt;
    result["parameter_fingerprint"] = provider_phase.parameter_fingerprint;
    result["globality_certificate"] = "not_guaranteed";
    return result;
}

py::dict held2_manufactured_stage_i(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::string& stage
) {
    if (stage == "stage_ii") {
        py::dict result = held2_stage_ii_to_dict(
            epcsaft_equilibrium::solve_held2_manufactured_stage_ii(
                charges, physical_feed
            )
        );
        result["profile"] = "perdomo-held2-stage-ii-manufactured-v1";
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
    py::list candidates;
    for (const epcsaft_equilibrium::Held2StageICandidate& candidate :
         evaluation.candidates) {
        py::dict item;
        item["modified_fractions"] = candidate.modified_fractions;
        item["volume"] = candidate.volume;
        item["tpd"] = candidate.tpd;
        candidates.append(std::move(item));
    }
    py::dict result;
    result["profile"] = "perdomo-held2-stage-i-manufactured-v1";
    result["outcome"] = evaluation.outcome;
    result["globality_certificate"] = "not_guaranteed";
    result["declared_start_count"] = evaluation.declared_start_count;
    result["completed_start_count"] = evaluation.completed_start_count;
    result["failed_start_count"] = evaluation.failed_start_count;
    result["reference_modified_fractions"] =
        evaluation.reference_modified_fractions;
    result["reference_volume"] = evaluation.reference_volume;
    result["minimum_tpd"] = evaluation.minimum_tpd;
    result["candidates"] = std::move(candidates);
    return result;
}

py::dict held2_pressure_envelope_to_dict(
    const epcsaft_equilibrium::Held2PressureEnvelopeResult& evaluation
) {
    py::list points;
    for (const epcsaft_equilibrium::Held2PressureScanPoint& point :
         evaluation.scan_points) {
        py::dict item;
        item["log_volume"] = point.log_volume;
        item["volume"] = point.volume;
        item["pressure_residual"] = point.pressure_residual;
        item["pressure_derivative_log_volume"] =
            point.pressure_derivative_log_volume;
        item["objective"] = point.objective;
        item["valid"] = point.valid;
        item["failure"] = point.failure;
        points.append(std::move(item));
    }
    py::list intervals;
    for (const epcsaft_equilibrium::Held2PressureScanInterval& interval :
         evaluation.intervals) {
        py::dict item;
        item["lower_log_volume"] = interval.lower_log_volume;
        item["upper_log_volume"] = interval.upper_log_volume;
        item["depth"] = interval.depth;
        item["status"] = interval.status;
        intervals.append(std::move(item));
    }
    py::list roots;
    for (const epcsaft_equilibrium::Held2PressureRoot& root : evaluation.roots) {
        py::dict item;
        item["log_volume"] = root.log_volume;
        item["volume"] = root.volume;
        item["objective"] = root.objective;
        item["pressure_residual"] = root.pressure_residual;
        item["pressure_derivative_log_volume"] =
            root.pressure_derivative_log_volume;
        item["objective_curvature_log_volume"] =
            root.objective_curvature_log_volume;
        item["mechanical_class"] = root.mechanical_class;
        item["origin"] = root.origin;
        item["boundary"] = root.boundary;
        roots.append(std::move(item));
    }
    py::dict result;
    result["outcome"] = evaluation.outcome;
    result["failure_reason"] = evaluation.failure_reason;
    result["root_completeness"] = evaluation.root_completeness;
    result["selection_scope"] = evaluation.selection_scope;
    result["selected_root_index"] = evaluation.selected_root_index;
    result["evaluation_failure_count"] = evaluation.evaluation_failure_count;
    result["refinement_failure_count"] = evaluation.refinement_failure_count;
    result["stationary_point_count"] = evaluation.stationary_point_count;
    result["tangential_root_count"] = evaluation.tangential_root_count;
    result["marginal_root_count"] = evaluation.marginal_root_count;
    result["boundary_root_count"] = evaluation.boundary_root_count;
    result["objective_tie_count"] = evaluation.objective_tie_count;
    result["deduplicated_root_count"] = evaluation.deduplicated_root_count;
    result["scan_points"] = std::move(points);
    result["intervals"] = std::move(intervals);
    result["roots"] = std::move(roots);
    result["globality_certificate"] = "not_guaranteed";
    return result;
}

py::dict held2_stage_i_reduced_evaluation_to_dict(
    const epcsaft_equilibrium::Held2StageIReducedEvaluation& evaluation
) {
    py::dict result;
    result["chart_coordinates"] = evaluation.chart_coordinates;
    result["independent_modified_fractions"] =
        evaluation.independent_modified_fractions;
    result["tpd"] = std::isfinite(evaluation.tpd)
        ? py::cast(evaluation.tpd)
        : py::none();
    result["certified"] = evaluation.certified;
    result["failure_reason"] = evaluation.failure_reason;
    result["pressure_envelope"] = held2_pressure_envelope_to_dict(
        evaluation.pressure_envelope
    );
    result["root_completeness"] =
        evaluation.pressure_envelope.root_completeness;
    result["pressure_certified"] = false;
    result["mechanical_class"] = py::none();
    result["root_origin"] = py::none();
    result["selected_root_log_volume"] = py::none();
    const int selected = evaluation.pressure_envelope.selected_root_index;
    if (selected >= 0
        && static_cast<std::size_t>(selected)
            < evaluation.pressure_envelope.roots.size()) {
        const epcsaft_equilibrium::Held2PressureRoot& root =
            evaluation.pressure_envelope.roots[static_cast<std::size_t>(selected)];
        result["pressure_certified"] = std::abs(root.pressure_residual) <= 1.0e-8;
        result["mechanical_class"] = root.mechanical_class;
        result["root_origin"] = root.origin;
        result["selected_root_log_volume"] = root.log_volume;
    }
    return result;
}

py::dict held2_stage_i_direct_to_dict(
    const epcsaft_equilibrium::Held2StageIDirectResult& evaluation
) {
    py::list evaluations;
    for (const epcsaft_equilibrium::Held2StageIReducedEvaluation& item :
         evaluation.evaluations) {
        evaluations.append(held2_stage_i_reduced_evaluation_to_dict(item));
    }
    py::dict result;
    result["outcome"] = evaluation.outcome;
    result["termination_reason"] = evaluation.termination_reason;
    result["search_strategy"] = evaluation.search_strategy;
    result["search_solver"] = evaluation.search_solver;
    result["solver_version"] = evaluation.solver_version;
    result["declared_evaluation_budget"] = evaluation.declared_evaluation_budget;
    result["completed_evaluation_count"] = evaluation.completed_evaluation_count;
    result["failed_evaluation_count"] = evaluation.failed_evaluation_count;
    result["minimum_tpd"] = std::isfinite(evaluation.minimum_tpd)
        ? py::cast(evaluation.minimum_tpd)
        : py::none();
    result["evaluations"] = std::move(evaluations);
    result["negative_witness"] = py::none();
    if (evaluation.negative_witness_index >= 0) {
        result["negative_witness"] = held2_stage_i_reduced_evaluation_to_dict(
            evaluation.evaluations[static_cast<std::size_t>(
                evaluation.negative_witness_index
            )]
        );
    }
    result["globality_certificate"] = evaluation.globality_certificate;
    return result;
}

py::dict held2_manufactured_stage_i_direct(
    const std::string& topology,
    int evaluation_budget
) {
    return held2_stage_i_direct_to_dict(
        epcsaft_equilibrium::solve_held2_manufactured_stage_i_direct(
            topology,
            evaluation_budget
        )
    );
}

py::dict held2_manufactured_pressure_envelope(
    const std::string& topology,
    double composition,
    int initial_interval_count
) {
    return held2_pressure_envelope_to_dict(
        epcsaft_equilibrium::evaluate_held2_manufactured_pressure_envelope(
            topology,
            composition,
            initial_interval_count
        )
    );
}

py::dict held2_stage_ii_basin_exploration_to_dict(
    const epcsaft_equilibrium::Held2StageIIBasinExplorationResult& evaluation
) {
    py::list evaluations;
    for (const epcsaft_equilibrium::Held2StageIIBasinEvaluation& item :
         evaluation.evaluations) {
        py::dict serialized;
        serialized["independent_modified_fractions"] =
            item.independent_modified_fractions;
        serialized["reduced_lower_value"] = std::isfinite(
            item.reduced_lower_value
        ) ? py::cast(item.reduced_lower_value) : py::none();
        serialized["certified"] = item.certified;
        serialized["failure_reason"] = item.failure_reason;
        serialized["pressure_envelope"] =
            held2_pressure_envelope_to_dict(item.pressure_envelope);
        evaluations.append(std::move(serialized));
    }
    py::list representatives;
    for (const epcsaft_equilibrium::Held2StageIIPhysicalStart& start :
         evaluation.representatives) {
        py::dict serialized;
        serialized["independent_modified_fractions"] =
            start.independent_modified_fractions;
        serialized["log_volume"] = start.log_volume;
        serialized["volume"] = start.volume;
        serialized["reduced_lower_value"] = start.reduced_lower_value;
        serialized["source"] = start.source;
        serialized["root_origin"] = start.root_origin;
        serialized["root_completeness"] = start.root_completeness;
        representatives.append(std::move(serialized));
    }
    py::dict result;
    result["outcome"] = evaluation.outcome;
    result["termination_reason"] = evaluation.termination_reason;
    result["strategy"] = evaluation.strategy;
    result["direct_solver"] = evaluation.direct_solver;
    result["direct_solver_version"] = evaluation.direct_solver_version;
    result["declared_sobol_count"] = evaluation.declared_sobol_count;
    result["declared_direct_budget"] = evaluation.declared_direct_budget;
    result["completed_evaluation_count"] =
        evaluation.completed_evaluation_count;
    result["failed_evaluation_count"] = evaluation.failed_evaluation_count;
    result["duplicate_start_count"] = evaluation.duplicate_start_count;
    result["direct_escalation_used"] = evaluation.direct_escalation_used;
    result["evaluations"] = std::move(evaluations);
    result["representatives"] = std::move(representatives);
    result["globality_certificate"] = evaluation.globality_certificate;
    return result;
}

py::dict held2_manufactured_stage_ii_basin_explorer(
    const std::string& topology,
    bool use_direct_escalation,
    int direct_evaluation_budget
) {
    const epcsaft_equilibrium::Held2Coordinates coordinates =
        epcsaft_equilibrium::make_held2_coordinates({0.0, 1.0, -1.0});
    std::vector<epcsaft_equilibrium::Held2StageIIBasinSeed> seeds;
    std::string pressure_topology;
    if (topology == "same_composition_different_density") {
        seeds = {{{0.5}, "reference"}, {{0.5}, "duplicate_reference"}};
        pressure_topology = "three_root";
    } else if (topology == "tied_density_branches") {
        seeds = {{{0.5}, "reference"}};
        pressure_topology = "tied";
    } else if (topology == "different_composition_same_density") {
        seeds = {{{0.25}, "left"}, {{0.75}, "right"}};
        pressure_topology = "one_root";
    } else if (topology == "duplicates" || topology == "direct") {
        seeds = {
            {{0.5}, "reference"},
            {{0.5}, "duplicate_reference"},
            {{0.5 + 0.25e-7}, "near_duplicate_reference"},
        };
        pressure_topology = "one_root";
    } else if (topology == "provider_failure") {
        seeds = {{{0.5}, "reference"}};
        pressure_topology = "invalid";
    } else {
        throw py::value_error("unknown manufactured Stage-II basin topology");
    }
    const epcsaft_equilibrium::Held2StageIIBasinEvaluator evaluator = [
        pressure_topology
    ](const std::vector<double>& independent) {
        epcsaft_equilibrium::Held2StageIIBasinEvaluation evaluation;
        evaluation.independent_modified_fractions = independent;
        evaluation.pressure_envelope =
            epcsaft_equilibrium::evaluate_held2_manufactured_pressure_envelope(
                pressure_topology,
                independent.front(),
                64
            );
        if (evaluation.pressure_envelope.outcome != "selected"
            && evaluation.pressure_envelope.failure_reason
                != "stable_objective_tie") {
            evaluation.failure_reason = evaluation.pressure_envelope.failure_reason;
            return evaluation;
        }
        for (const epcsaft_equilibrium::Held2PressureRoot& root :
             evaluation.pressure_envelope.roots) {
            if (root.mechanical_class == "strict_stable" && !root.boundary) {
                evaluation.reduced_lower_value = std::min(
                    evaluation.reduced_lower_value,
                    root.objective + std::pow(independent.front() - 0.5, 2)
                );
            }
        }
        if (!std::isfinite(evaluation.reduced_lower_value)) {
            evaluation.failure_reason = "no_strict_stable_root";
            return evaluation;
        }
        evaluation.certified = true;
        return evaluation;
    };
    const epcsaft_equilibrium::Held2StageIIBasinExplorationResult evaluation =
        epcsaft_equilibrium::explore_held2_stage_ii_basins(
            coordinates,
            seeds,
            0,
            use_direct_escalation,
            direct_evaluation_budget,
            evaluator
        );
    return held2_stage_ii_basin_exploration_to_dict(evaluation);
}

py::dict held2_installed_pressure_envelope(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& independent_modified_fractions,
    const std::string& expected_fingerprint,
    int initial_interval_count
) {
    if (!std::isfinite(temperature_k) || !std::isfinite(pressure_pa)
        || temperature_k <= 0.0 || pressure_pa <= 0.0) {
        throw py::value_error(
            "HELD2 temperature and pressure must be finite and positive"
        );
    }
    const epcsaft_native_sdk_v1& sdk = checked_molar_volume_sdk(capsule);
    std::vector<double> charges;
    std::vector<std::string> component_ids;
    charges.reserve(sdk.component_count);
    component_ids.reserve(sdk.component_count);
    for (std::size_t component = 0; component < sdk.component_count; ++component) {
        if (sdk.component_ids[component] == nullptr
            || sdk.component_ids[component][0] == '\0') {
            throw py::value_error(
                "provider electrolyte component ID must not be empty"
            );
        }
        charges.push_back(static_cast<double>(sdk.component_charges[component]));
        component_ids.emplace_back(sdk.component_ids[component]);
    }
    const epcsaft_equilibrium::Held2Coordinates coordinates =
        epcsaft_equilibrium::make_held2_coordinates(charges);
    const std::vector<double> physical_amounts =
        epcsaft_equilibrium::held2_lift_independent_fractions(
            coordinates,
            independent_modified_fractions
        );
    const epcsaft_equilibrium::ProviderContext provider(
        sdk,
        expected_fingerprint
    );
    const std::array<double, 2> molar_volume_bounds =
        provider.evaluate_molar_volume_bounds(
            temperature_k,
            physical_amounts,
            epcsaft_equilibrium::kHeld2PackingFractionMinimum,
            epcsaft_equilibrium::kHeld2PackingFractionMaximum
        );
    const double pressure_over_rt = pressure_pa
        / (kGasConstantJPerMolK * temperature_k);
    const epcsaft_equilibrium::Held2StateEvaluator evaluator = [
        &provider,
        coordinates,
        temperature_k,
        pressure_pa,
        pressure_over_rt
    ](
        const std::vector<double>& independent,
        double log_volume
    ) {
        const std::vector<double> amounts =
            epcsaft_equilibrium::held2_lift_independent_fractions(
                coordinates,
                independent
            );
        const double volume = std::exp(log_volume);
        const epcsaft_equilibrium::MixturePhaseEvaluation provider_phase =
            provider.evaluate_electrolyte(
                temperature_k,
                amounts,
                volume
            );
        epcsaft_equilibrium::Held2PhysicalPhaseBlock block;
        block.helmholtz_over_rt = provider_phase.value;
        block.gradient = provider_phase.gradient;
        block.hessian = provider_phase.hessian;
        block.pressure_pa = provider_phase.pressure_pa;
        return epcsaft_equilibrium::evaluate_held2_phase_block(
            coordinates,
            independent,
            log_volume,
            pressure_over_rt,
            pressure_pa,
            block
        );
    };
    py::dict result = held2_pressure_envelope_to_dict(
        epcsaft_equilibrium::evaluate_held2_pressure_envelope(
            independent_modified_fractions,
            molar_volume_bounds,
            evaluator,
            initial_interval_count,
            8
        )
    );
    result["component_ids"] = std::move(component_ids);
    result["charges"] = std::move(charges);
    result["molar_volume_bounds"] = molar_volume_bounds;
    result["parameter_fingerprint"] = expected_fingerprint;
    return result;
}

py::dict held2_installed_stage_ii_basin_explorer(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::vector<double>>& independent_seeds,
    const std::string& expected_fingerprint,
    int sobol_count
) {
    if (!std::isfinite(temperature_k) || !std::isfinite(pressure_pa)
        || temperature_k <= 0.0 || pressure_pa <= 0.0) {
        throw py::value_error(
            "HELD2 temperature and pressure must be finite and positive"
        );
    }
    const epcsaft_native_sdk_v1& sdk = checked_molar_volume_sdk(capsule);
    std::vector<double> charges;
    charges.reserve(sdk.component_count);
    for (std::size_t component = 0; component < sdk.component_count; ++component) {
        charges.push_back(static_cast<double>(sdk.component_charges[component]));
    }
    const epcsaft_equilibrium::Held2Coordinates coordinates =
        epcsaft_equilibrium::make_held2_coordinates(charges);
    const epcsaft_equilibrium::ProviderContext provider(
        sdk,
        expected_fingerprint
    );
    const double pressure_over_rt = pressure_pa
        / (kGasConstantJPerMolK * temperature_k);
    const epcsaft_equilibrium::Held2StateEvaluator phase_evaluator = [
        &provider,
        coordinates,
        temperature_k,
        pressure_pa,
        pressure_over_rt
    ](
        const std::vector<double>& independent,
        double log_volume
    ) {
        const std::vector<double> amounts =
            epcsaft_equilibrium::held2_lift_independent_fractions(
                coordinates,
                independent
            );
        const epcsaft_equilibrium::MixturePhaseEvaluation provider_phase =
            provider.evaluate_electrolyte(
                temperature_k,
                amounts,
                std::exp(log_volume)
            );
        epcsaft_equilibrium::Held2PhysicalPhaseBlock block;
        block.helmholtz_over_rt = provider_phase.value;
        block.gradient = provider_phase.gradient;
        block.hessian = provider_phase.hessian;
        block.pressure_pa = provider_phase.pressure_pa;
        return epcsaft_equilibrium::evaluate_held2_phase_block(
            coordinates,
            independent,
            log_volume,
            pressure_over_rt,
            pressure_pa,
            block
        );
    };
    const epcsaft_equilibrium::Held2StageIIBasinEvaluator evaluator = [
        &provider,
        &phase_evaluator,
        coordinates,
        temperature_k
    ](const std::vector<double>& independent) {
        epcsaft_equilibrium::Held2StageIIBasinEvaluation evaluation;
        evaluation.independent_modified_fractions = independent;
        const std::vector<double> amounts =
            epcsaft_equilibrium::held2_lift_independent_fractions(
                coordinates,
                independent
            );
        const std::array<double, 2> bounds =
            provider.evaluate_molar_volume_bounds(
                temperature_k,
                amounts,
                epcsaft_equilibrium::kHeld2PackingFractionMinimum,
                epcsaft_equilibrium::kHeld2PackingFractionMaximum
            );
        evaluation.pressure_envelope =
            epcsaft_equilibrium::evaluate_held2_pressure_envelope(
                independent,
                bounds,
                phase_evaluator,
                64,
                8
            );
        if (evaluation.pressure_envelope.outcome != "selected"
            && evaluation.pressure_envelope.failure_reason
                != "stable_objective_tie") {
            evaluation.failure_reason = evaluation.pressure_envelope.failure_reason;
            return evaluation;
        }
        for (const epcsaft_equilibrium::Held2PressureRoot& root :
             evaluation.pressure_envelope.roots) {
            if (root.mechanical_class == "strict_stable" && !root.boundary) {
                evaluation.reduced_lower_value = std::min(
                    evaluation.reduced_lower_value,
                    root.objective
                );
            }
        }
        if (!std::isfinite(evaluation.reduced_lower_value)) {
            evaluation.failure_reason = "no_strict_stable_root";
            return evaluation;
        }
        evaluation.certified = true;
        return evaluation;
    };
    std::vector<epcsaft_equilibrium::Held2StageIIBasinSeed> seeds;
    seeds.reserve(independent_seeds.size());
    for (const std::vector<double>& independent : independent_seeds) {
        seeds.push_back({independent, "external_seed"});
    }
    const epcsaft_equilibrium::Held2StageIIBasinExplorationResult evaluation =
        epcsaft_equilibrium::explore_held2_stage_ii_basins(
            coordinates,
            seeds,
            sobol_count,
            false,
            0,
            evaluator
        );
    py::dict result = held2_stage_ii_basin_exploration_to_dict(evaluation);
    result["parameter_fingerprint"] = expected_fingerprint;
    return result;
}

class InstalledHeld2Problem {
public:
    InstalledHeld2Problem(
        const epcsaft_native_sdk_v1& sdk,
        double temperature_k,
        double pressure_pa,
        std::vector<double> physical_feed,
        std::string expected_fingerprint
    )
        : temperature_k_(temperature_k),
          pressure_pa_(pressure_pa),
          physical_feed_(std::move(physical_feed)),
          coordinates_(epcsaft_equilibrium::make_held2_coordinates(
              charges_from_sdk(sdk)
          )),
          provider_(sdk, std::move(expected_fingerprint)),
          pressure_over_rt_(pressure_pa_
              / (kGasConstantJPerMolK * temperature_k_)) {
        if (!std::isfinite(temperature_k_) || !std::isfinite(pressure_pa_)
            || temperature_k_ <= 0.0 || pressure_pa_ <= 0.0) {
            throw py::value_error(
                "HELD2 temperature and pressure must be finite and positive"
            );
        }
        const std::vector<double> modified_feed =
            epcsaft_equilibrium::held2_transform_physical_fractions(
                coordinates_, physical_feed_
            );
        independent_feed_.reserve(coordinates_.independent_indices.size());
        for (std::size_t component : coordinates_.independent_indices) {
            const auto retained = std::find(
                coordinates_.retained_indices.begin(),
                coordinates_.retained_indices.end(),
                component
            );
            if (retained == coordinates_.retained_indices.end()) {
                throw py::value_error(
                    "HELD2 independent component is not retained"
                );
            }
            independent_feed_.push_back(modified_feed[static_cast<std::size_t>(
                retained - coordinates_.retained_indices.begin()
            )]);
        }
    }

    [[nodiscard]] epcsaft_equilibrium::Held2StateEvaluation evaluate(
        const std::vector<double>& independent,
        double log_volume
    ) const {
        const std::vector<double> amounts =
            epcsaft_equilibrium::held2_lift_independent_fractions(
                coordinates_, independent
            );
        const epcsaft_equilibrium::MixturePhaseEvaluation provider_phase =
            provider_.evaluate_electrolyte(
                temperature_k_, amounts, std::exp(log_volume)
            );
        epcsaft_equilibrium::Held2PhysicalPhaseBlock block;
        block.helmholtz_over_rt = provider_phase.value;
        block.gradient = provider_phase.gradient;
        block.hessian = provider_phase.hessian;
        block.pressure_pa = provider_phase.pressure_pa;
        return epcsaft_equilibrium::evaluate_held2_phase_block(
            coordinates_,
            independent,
            log_volume,
            pressure_over_rt_,
            pressure_pa_,
            block
        );
    }

    [[nodiscard]] std::array<double, 2> volume_bounds(
        const std::vector<double>& independent
    ) const {
        const std::vector<double> amounts =
            epcsaft_equilibrium::held2_lift_independent_fractions(
                coordinates_, independent
            );
        return provider_.evaluate_molar_volume_bounds(
            temperature_k_,
            amounts,
            epcsaft_equilibrium::kHeld2PackingFractionMinimum,
            epcsaft_equilibrium::kHeld2PackingFractionMaximum
        );
    }

    [[nodiscard]] epcsaft_equilibrium::Held2PressureEnvelopeResult envelope(
        const std::vector<double>& independent
    ) const {
        const auto evaluator = [this](const std::vector<double>& values, double q) {
            return evaluate(values, q);
        };
        return epcsaft_equilibrium::evaluate_held2_pressure_envelope(
            independent, volume_bounds(independent), evaluator, 64, 8
        );
    }

    [[nodiscard]] const epcsaft_equilibrium::Held2Coordinates& coordinates() const {
        return coordinates_;
    }

    [[nodiscard]] const std::vector<double>& physical_feed() const {
        return physical_feed_;
    }

    [[nodiscard]] const std::vector<double>& independent_feed() const {
        return independent_feed_;
    }

    [[nodiscard]] const std::string& fingerprint() const {
        return provider_.fingerprint();
    }

private:
    static std::vector<double> charges_from_sdk(
        const epcsaft_native_sdk_v1& sdk
    ) {
        std::vector<double> charges;
        charges.reserve(sdk.component_count);
        for (std::size_t component = 0; component < sdk.component_count; ++component) {
            charges.push_back(static_cast<double>(sdk.component_charges[component]));
        }
        return charges;
    }

    double temperature_k_;
    double pressure_pa_;
    std::vector<double> physical_feed_;
    epcsaft_equilibrium::Held2Coordinates coordinates_;
    epcsaft_equilibrium::ProviderContext provider_;
    double pressure_over_rt_;
    std::vector<double> independent_feed_;
};

struct InstalledHeld2StageI {
    epcsaft_equilibrium::Held2PressureEnvelopeResult reference_envelope;
    epcsaft_equilibrium::Held2StateEvaluation reference;
    epcsaft_equilibrium::Held2StageIDirectResult search;
};

InstalledHeld2StageI run_installed_held2_stage_i(
    const InstalledHeld2Problem& problem,
    int evaluation_budget
) {
    InstalledHeld2StageI result;
    try {
        result.reference_envelope = problem.envelope(problem.independent_feed());
    } catch (const std::exception& error) {
        result.search.declared_evaluation_budget = evaluation_budget;
        result.search.termination_reason = std::string(
            "reference_envelope_failed: "
        ) + error.what();
        return result;
    }
    if (result.reference_envelope.outcome != "selected") {
        result.search.declared_evaluation_budget = evaluation_budget;
        result.search.termination_reason = "reference_envelope_failed: "
            + result.reference_envelope.failure_reason;
        return result;
    }
    result.reference = result.reference_envelope.roots[static_cast<std::size_t>(
        result.reference_envelope.selected_root_index
    )].state;
    const epcsaft_equilibrium::Held2StageIReducedEvaluator evaluator = [
        &problem,
        reference = result.reference
    ](const std::vector<double>& chart_coordinates) {
        epcsaft_equilibrium::Held2StageIReducedEvaluation evaluation;
        evaluation.chart_coordinates = chart_coordinates;
        try {
            evaluation.independent_modified_fractions =
                epcsaft_equilibrium::held2_map_unit_cube_to_independent_fractions(
                    problem.coordinates(), chart_coordinates
                );
            evaluation.pressure_envelope = problem.envelope(
                evaluation.independent_modified_fractions
            );
            if (evaluation.pressure_envelope.outcome != "selected") {
                evaluation.failure_reason =
                    evaluation.pressure_envelope.failure_reason;
                return evaluation;
            }
            const epcsaft_equilibrium::Held2StateEvaluation& selected =
                evaluation.pressure_envelope.roots[static_cast<std::size_t>(
                    evaluation.pressure_envelope.selected_root_index
                )].state;
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
            evaluation.failure_reason = std::string(
                "provider_evaluation_failed: "
            ) + error.what();
        }
        return evaluation;
    };
    result.search = epcsaft_equilibrium::solve_held2_stage_i_direct(
        problem.coordinates().independent_indices.size(),
        evaluation_budget,
        -1.0e-8,
        evaluator
    );
    return result;
}

py::dict held2_installed_stage_i_direct(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& physical_feed,
    const std::string& expected_fingerprint,
    int evaluation_budget
) {
    const epcsaft_native_sdk_v1& sdk = checked_molar_volume_sdk(capsule);
    const InstalledHeld2Problem problem(
        sdk,
        temperature_k,
        pressure_pa,
        physical_feed,
        expected_fingerprint
    );
    const InstalledHeld2StageI stage_i = run_installed_held2_stage_i(
        problem, evaluation_budget
    );
    py::dict result = held2_stage_i_direct_to_dict(stage_i.search);
    result["reference_pressure_envelope"] =
        held2_pressure_envelope_to_dict(stage_i.reference_envelope);
    result["reference_independent_modified_fractions"] =
        problem.independent_feed();
    result["parameter_fingerprint"] = problem.fingerprint();
    return result;
}

py::dict held2_stage_ii_to_dict(
    const epcsaft_equilibrium::Held2StageIIResult& evaluation
) {
    py::list bounds;
    for (const auto& bound : evaluation.bound_history) {
        py::dict item;
        item["lower_bound"] = bound.lower_bound;
        item["upper_bound"] = bound.upper_bound;
        item["multipliers"] = bound.multipliers;
        item["cut_count"] = bound.cut_count;
        item["upper_solver"] = bound.upper_solver;
        item["upper_solver_version"] = bound.upper_solver_version;
        item["upper_solver_status"] = bound.upper_solver_status;
        item["upper_primal_feasible"] = bound.upper_primal_feasible;
        item["upper_dual_feasible"] = bound.upper_dual_feasible;
        item["upper_primal_residual_inf"] = bound.upper_primal_residual_inf;
        item["upper_dual_residual_inf"] = bound.upper_dual_residual_inf;
        item["cut_slacks"] = bound.cut_slacks;
        item["cut_duals"] = bound.cut_duals;
        item["active_cut_ids"] = bound.active_cut_ids;
        bounds.append(std::move(item));
    }
    py::list attempts;
    int solver_converged_count = 0;
    int physical_kkt_passed_count = 0;
    int step6_eligible_count = 0;
    for (const auto& attempt : evaluation.attempt_trace) {
        py::dict item;
        item["attempt_id"] = attempt.attempt_id;
        item["major_iteration"] = attempt.major_iteration;
        item["start_index"] = attempt.start_index;
        item["start_source"] = attempt.start_source;
        item["internal_start"] = attempt.internal_start;
        item["physical_start_modified_fractions"] =
            attempt.physical_start_modified_fractions;
        item["physical_start_volume"] = attempt.physical_start_volume;
        item["solver_status"] = attempt.solver_status;
        item["solver_converged"] = attempt.solver_converged;
        item["provider_status"] = attempt.provider_status;
        item["callback_error"] = attempt.callback_error;
        item["internal_terminal"] = attempt.internal_terminal;
        item["terminal_modified_fractions"] =
            attempt.terminal_modified_fractions;
        item["terminal_volume"] = attempt.terminal_volume;
        item["objective"] = attempt.objective;
        item["lower_value"] = attempt.lower_value;
        item["pressure_residual"] = attempt.pressure_residual;
        item["lower_bound_multipliers"] = attempt.lower_bound_multipliers;
        item["upper_bound_multipliers"] = attempt.upper_bound_multipliers;
        item["chart_jacobian_condition"] = attempt.chart_jacobian_condition;
        item["dual_pullback_inf_norm"] = attempt.dual_pullback_inf_norm;
        item["chart_kkt_inf_norm"] = attempt.chart_kkt_inf_norm;
        item["physical_kkt_inf_norm"] = attempt.physical_kkt_inf_norm;
        item["complementarity_inf_norm"] = attempt.complementarity_inf_norm;
        item["pressure_passed"] = attempt.pressure_passed;
        item["dual_signs_valid"] = attempt.dual_signs_valid;
        item["physical_kkt_passed"] = attempt.physical_kkt_passed;
        item["cut_eligible"] = attempt.cut_eligible;
        item["step6_eligible"] = attempt.step6_eligible;
        item["basin_id"] = attempt.basin_id;
        item["same_major_upper_bound"] = attempt.same_major_upper_bound;
        item["same_major_multipliers"] = attempt.same_major_multipliers;
        attempts.append(std::move(item));
        solver_converged_count += attempt.solver_converged ? 1 : 0;
        physical_kkt_passed_count += attempt.physical_kkt_passed ? 1 : 0;
        step6_eligible_count += attempt.step6_eligible ? 1 : 0;
    }
    py::list candidates;
    for (const auto& candidate : evaluation.candidates) {
        py::dict item;
        item["modified_fractions"] = candidate.modified_fractions;
        item["independent_modified_fractions"] =
            candidate.independent_modified_fractions;
        item["volume"] = candidate.volume;
        item["phase_coordinate"] = candidate.phase_coordinate;
        item["lower_gap"] = candidate.lower_gap;
        candidates.append(std::move(item));
    }
    py::dict result;
    result["outcome"] = evaluation.outcome;
    result["search_strategy"] = evaluation.search_strategy;
    result["global_explorer"] = evaluation.global_explorer;
    result["local_solver"] = evaluation.local_solver;
    result["globality_certificate"] = evaluation.globality_certificate;
    result["major_iterations"] = evaluation.major_iterations;
    result["lower_starts_per_iteration"] =
        evaluation.lower_starts_per_iteration;
    result["cut_count"] = evaluation.cut_count;
    result["exploration_evaluation_count"] =
        evaluation.exploration_evaluation_count;
    result["exploration_failure_count"] =
        evaluation.exploration_failure_count;
    result["exploration_representative_count"] =
        evaluation.exploration_representative_count;
    result["duplicate_representative_count"] =
        evaluation.duplicate_representative_count;
    result["duplicate_terminal_count"] = evaluation.duplicate_terminal_count;
    result["distinct_basin_count"] = evaluation.distinct_basin_count;
    result["local_attempt_cap_per_major"] =
        evaluation.local_attempt_cap_per_major;
    result["local_attempts_truncated"] = evaluation.local_attempts_truncated;
    result["direct_escalation_used"] = evaluation.direct_escalation_used;
    result["bound_history"] = std::move(bounds);
    result["attempt_trace"] = std::move(attempts);
    py::dict attempt_classification;
    attempt_classification["declared"] = evaluation.attempt_trace.size();
    attempt_classification["solver_converged"] = solver_converged_count;
    attempt_classification["solver_failed"] =
        static_cast<int>(evaluation.attempt_trace.size())
        - solver_converged_count;
    attempt_classification["physical_kkt_passed"] =
        physical_kkt_passed_count;
    attempt_classification["step6_eligible"] = step6_eligible_count;
    result["attempt_classification"] = std::move(attempt_classification);
    result["candidates"] = std::move(candidates);
    return result;
}

py::dict held2_stage_iii_to_dict(
    const epcsaft_equilibrium::Held2StageIIIResult& evaluation
) {
    py::list phases;
    for (const auto& phase : evaluation.phases) {
        py::dict item;
        item["phase_fraction"] = phase.phase_fraction;
        item["modified_fractions"] = phase.modified_fractions;
        item["physical_fractions"] = phase.physical_fractions;
        item["volume"] = phase.volume;
        phases.append(std::move(item));
    }
    py::list lifecycle;
    for (const auto& step : evaluation.lifecycle) {
        py::dict item;
        item["solve_index"] = step.solve_index;
        item["active_candidate_count"] = step.active_candidate_count;
        item["removed_candidate_index"] = step.removed_candidate_index;
        item["action"] = step.action;
        item["phase_fraction"] = step.phase_fraction;
        item["lower_bound_multiplier"] = step.lower_bound_multiplier;
        item["reduced_derivative"] = step.reduced_derivative;
        item["complementarity_inf_norm"] = step.complementarity_inf_norm;
        item["candidate_independent_modified_fractions"] =
            step.candidate_independent_modified_fractions;
        item["candidate_volume"] = step.candidate_volume;
        item["solver_status"] = step.solver_status;
        item["decision_reason"] = step.decision_reason;
        lifecycle.append(std::move(item));
    }
    py::dict result;
    result["solver_status"] = evaluation.solver_status;
    result["numerical_status"] = evaluation.numerical_status;
    result["physical_status"] = evaluation.physical_status;
    result["feedback"] = evaluation.feedback;
    result["failure_reason"] = evaluation.failure_reason;
    result["trace_refinement_status"] = evaluation.trace_refinement_status;
    result["input_candidate_count"] = evaluation.input_candidate_count;
    result["retired_duplicate_count"] = evaluation.retired_duplicate_count;
    result["retired_inactive_count"] = evaluation.retired_inactive_count;
    result["stage_iii_solve_count"] = evaluation.stage_iii_solve_count;
    result["active_set_resolve_count"] = evaluation.active_set_resolve_count;
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
    result["kkt_stationarity_inf_norm"] =
        evaluation.kkt_stationarity_inf_norm;
    result["dual_sign_violation_inf_norm"] =
        evaluation.dual_sign_violation_inf_norm;
    result["bound_complementarity_inf_norm"] =
        evaluation.bound_complementarity_inf_norm;
    result["minimum_phase_fraction"] = evaluation.minimum_phase_fraction;
    result["enumeration_objective_gap"] = evaluation.enumeration_objective_gap;
    result["phases"] = std::move(phases);
    result["lifecycle"] = std::move(lifecycle);
    result["globality_certificate"] = "not_guaranteed";
    return result;
}

py::dict held2_installed_controller(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& physical_feed,
    const std::string& expected_fingerprint,
    int stage_i_evaluation_budget,
    int stage_ii_major_iteration_cap,
    int stage_ii_local_attempt_cap_per_major
) {
    const epcsaft_native_sdk_v1& sdk = checked_molar_volume_sdk(capsule);
    const InstalledHeld2Problem problem(
        sdk,
        temperature_k,
        pressure_pa,
        physical_feed,
        expected_fingerprint
    );
    const InstalledHeld2StageI stage_i = run_installed_held2_stage_i(
        problem, stage_i_evaluation_budget
    );
    py::dict result;
    result["controller"] = "perdomo_held2_integrated_private_v1";
    result["stage_order"] = py::make_tuple(
        "homogeneous_reference", "stage_i", "stage_ii", "stage_iii"
    );
    result["parameter_fingerprint"] = problem.fingerprint();
    result["globality_certificate"] = "not_guaranteed";
    result["reference_pressure_envelope"] =
        held2_pressure_envelope_to_dict(stage_i.reference_envelope);
    result["stage_i"] = held2_stage_i_direct_to_dict(stage_i.search);
    result["stage_ii"] = py::none();
    result["stage_iii"] = py::none();
    result["predictive_comparison_status"] =
        "not_allowed_before_physical_acceptance";

    if (stage_i.search.outcome != "negative_witness_found"
        || stage_i.search.negative_witness_index < 0) {
        result["outcome"] = "stage_i_finite_search_without_negative_witness";
        result["failure_stage"] = "stage_i";
        return result;
    }
    const auto& witness = stage_i.search.evaluations[static_cast<std::size_t>(
        stage_i.search.negative_witness_index
    )];
    if (witness.pressure_envelope.selected_root_index < 0) {
        result["outcome"] = "indeterminate_stage_i_witness";
        result["failure_stage"] = "stage_i";
        return result;
    }
    const auto& witness_state = witness.pressure_envelope.roots[
        static_cast<std::size_t>(witness.pressure_envelope.selected_root_index)
    ].state;
    const std::vector<epcsaft_equilibrium::Held2StageICandidate> witnesses = {{
        witness_state.modified_fractions,
        witness_state.volume,
        witness.tpd,
    }};
    const epcsaft_equilibrium::Held2StateEvaluator evaluator = [
        &problem
    ](const std::vector<double>& independent, double q) {
        return problem.evaluate(independent, q);
    };
    const epcsaft_equilibrium::Held2VolumeBoundsEvaluator bounds_evaluator = [
        &problem
    ](const std::vector<double>& independent) {
        return problem.volume_bounds(independent);
    };
    const epcsaft_equilibrium::Held2StageIIResult stage_ii =
        epcsaft_equilibrium::solve_held2_stage_ii(
            problem.coordinates(),
            problem.physical_feed(),
            evaluator,
            bounds_evaluator,
            stage_i.reference,
            witnesses,
            stage_ii_major_iteration_cap,
            stage_ii_local_attempt_cap_per_major
        );
    result["stage_ii"] = held2_stage_ii_to_dict(stage_ii);
    if (stage_ii.outcome != "candidate_set" || stage_ii.candidates.size() < 2) {
        result["outcome"] = "indeterminate_stage_ii";
        result["failure_stage"] = "stage_ii";
        return result;
    }

    std::vector<std::array<double, 2>> phase_coordinate_bounds;
    phase_coordinate_bounds.reserve(stage_ii.candidates.size());
    for (const auto& candidate : stage_ii.candidates) {
        const std::array<double, 2> volume_bounds = problem.volume_bounds(
            candidate.independent_modified_fractions
        );
        phase_coordinate_bounds.push_back({
            std::log(volume_bounds[0]), std::log(volume_bounds[1])
        });
    }
    const epcsaft_equilibrium::Held2StageIIIResult stage_iii =
        epcsaft_equilibrium::solve_held2_stage_iii(
            problem.coordinates(),
            problem.physical_feed(),
            stage_ii.candidates,
            evaluator,
            phase_coordinate_bounds
        );
    result["stage_iii"] = held2_stage_iii_to_dict(stage_iii);
    if (stage_iii.physical_status != "accepted") {
        result["outcome"] = "indeterminate_stage_iii";
        result["failure_stage"] = "stage_iii";
        return result;
    }
    result["outcome"] = "physical_equilibrium_accepted";
    result["failure_stage"] = py::none();
    result["predictive_comparison_status"] =
        "eligible_but_not_executed_private_controller";
    return result;
}

py::dict held2_stage_ii_upper_lp(
    const std::vector<double>& intercepts,
    const std::vector<std::vector<double>>& slopes,
    double value_lower_bound,
    const std::string& solver
) {
    if (intercepts.size() != slopes.size()) {
        throw py::value_error("HELD2 upper LP cut intercept and slope counts differ");
    }
    const std::size_t dimension = slopes.empty() ? 1 : slopes.front().size();
    epcsaft_equilibrium::Held2StageIIUpperProblem problem;
    problem.multiplier_lower_bounds.assign(
        dimension,
        -std::numeric_limits<double>::infinity()
    );
    problem.multiplier_upper_bounds.assign(
        dimension,
        std::numeric_limits<double>::infinity()
    );
    problem.value_lower_bound = value_lower_bound;
    for (std::size_t index = 0; index < intercepts.size(); ++index) {
        problem.cuts.push_back({
            static_cast<int>(index),
            intercepts[index],
            slopes[index],
        });
    }
    epcsaft_equilibrium::Held2StageIIUpperResult evaluation;
    if (solver == "highs") {
        evaluation = epcsaft_equilibrium::solve_held2_stage_ii_upper_highs(problem);
    } else if (solver == "analytic_1d_test_oracle") {
        evaluation =
            epcsaft_equilibrium::solve_held2_stage_ii_upper_analytic_1d(problem);
    } else {
        throw py::value_error("unsupported HELD2 upper LP solver request");
    }
    py::dict result;
    result["outcome"] = evaluation.outcome;
    result["solver"] = evaluation.solver;
    result["solver_status"] = evaluation.solver_status;
    result["solver_version"] = evaluation.solver_version;
    result["solver_finished"] = evaluation.solver_finished;
    result["primal_feasible"] = evaluation.primal_feasible;
    result["dual_feasible"] = evaluation.dual_feasible;
    result["upper_bound"] = evaluation.outcome == "optimal"
        ? py::cast(evaluation.upper_bound)
        : py::none();
    result["multipliers"] = evaluation.multipliers;
    result["cut_slacks"] = evaluation.cut_slacks;
    result["cut_duals"] = evaluation.cut_duals;
    result["active_cut_ids"] = evaluation.active_cut_ids;
    result["primal_residual_inf"] = evaluation.primal_residual_inf;
    result["dual_residual_inf"] = evaluation.dual_residual_inf;
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

py::dict held2_manufactured_stage_iii(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::vector<std::array<double, 2>>& candidates,
    const std::string& stage
) {
    if (stage != "stage_iii") {
        throw py::value_error("unsupported manufactured HELD2 Stage III request");
    }
    py::dict result = held2_stage_iii_to_dict(
        epcsaft_equilibrium::solve_held2_manufactured_stage_iii(
            charges, physical_feed, candidates
        )
    );
    result["profile"] = "perdomo-held2-stage-iii-manufactured-v1";
    return result;
}

py::dict held2_stage_iii_retirement_decision(
    double phase_fraction,
    double lower_bound_multiplier,
    double upper_bound_multiplier,
    double reduced_derivative,
    bool remaining_balance_feasible
) {
    const epcsaft_equilibrium::Held2StageIIIRetirementDecision decision =
        epcsaft_equilibrium::held2_stage_iii_retirement_decision(
            phase_fraction,
            lower_bound_multiplier,
            upper_bound_multiplier,
            reduced_derivative,
            remaining_balance_feasible
        );
    py::dict result;
    result["retire"] = decision.retire;
    result["reason"] = decision.reason;
    result["complementarity_inf_norm"] = decision.complementarity_inf_norm;
    result["stationarity_residual"] = decision.stationarity_residual;
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
        &held2_manufactured_search_objective_evidence,
        py::arg("charges"),
        py::arg("variables"),
        py::arg("reference_variables"),
        py::arg("use_tpd"),
        py::arg("stage")
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
        &held2_manufactured_stage_i,
        py::arg("charges"),
        py::arg("physical_feed"),
        py::arg("stage")
    );
    module.def(
        "_held2_pressure_envelope",
        &held2_manufactured_pressure_envelope,
        py::arg("topology"),
        py::arg("composition"),
        py::arg("initial_interval_count")
    );
    module.def(
        "_held2_stage_i_direct",
        &held2_manufactured_stage_i_direct,
        py::arg("topology"),
        py::arg("evaluation_budget")
    );
    module.def(
        "_held2_stage_i_direct",
        &held2_installed_stage_i_direct,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("physical_feed"),
        py::arg("expected_fingerprint"),
        py::arg("evaluation_budget")
    );
    module.def(
        "_held2_controller",
        &held2_installed_controller,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("physical_feed"),
        py::arg("expected_fingerprint"),
        py::arg("stage_i_evaluation_budget") = 50,
        py::arg("stage_ii_major_iteration_cap") = 8,
        py::arg("stage_ii_local_attempt_cap_per_major") = 50
    );
    module.def(
        "_held2_stage_ii_upper_lp",
        &held2_stage_ii_upper_lp,
        py::arg("intercepts"),
        py::arg("slopes"),
        py::arg("value_lower_bound") =
            -std::numeric_limits<double>::infinity(),
        py::arg("solver") = "highs"
    );
    module.def(
        "_held2_stage_ii_basin_explorer",
        &held2_manufactured_stage_ii_basin_explorer,
        py::arg("topology"),
        py::arg("use_direct_escalation") = false,
        py::arg("direct_evaluation_budget") = 0
    );
    module.def(
        "_held2_stage_ii_basin_explorer",
        &held2_installed_stage_ii_basin_explorer,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("independent_seeds"),
        py::arg("expected_fingerprint"),
        py::arg("sobol_count") = 0
    );
    module.def(
        "_held2_pressure_envelope",
        &held2_installed_pressure_envelope,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("independent_modified_fractions"),
        py::arg("expected_fingerprint"),
        py::arg("initial_interval_count")
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
        "_held2_stage_iii_retirement_decision",
        &held2_stage_iii_retirement_decision,
        py::arg("phase_fraction"),
        py::arg("lower_bound_multiplier"),
        py::arg("upper_bound_multiplier"),
        py::arg("reduced_derivative"),
        py::arg("remaining_balance_feasible")
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
