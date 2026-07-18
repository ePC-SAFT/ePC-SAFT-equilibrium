#include <array>
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
#include "saturation.hpp"
#include "two_phase_flash.hpp"

namespace py = pybind11;

namespace {

constexpr std::string_view kFlashFingerprint =
    "sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286";
constexpr std::size_t kPureSdkTableSize = offsetof(epcsaft_native_sdk_v1, component_count);
constexpr std::size_t kMixtureSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, evaluate_mixture_phase)
    + sizeof(epcsaft_evaluate_mixture_phase_v1);

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
        py::list rejections;
        for (const auto& rejection : entry.rejections) {
            py::dict item;
            item["identity"] = rejection.identity;
            item["reason"] = rejection.reason;
            rejections.append(std::move(item));
        }
        py::dict item;
        item["major_iteration"] = entry.major_iteration;
        item["outer_value"] = entry.outer_value;
        item["upper_bound"] = entry.upper_bound;
        item["multiplier"] = entry.multiplier;
        item["active_cut_ids"] = entry.active_cut_ids;
        item["accepted_cut_ids"] = entry.accepted_cut_ids;
        item["lower_starts_completed"] = entry.lower_starts_completed;
        item["candidate_ids"] = entry.candidate_ids;
        item["rejections"] = rejections;
        trace.append(std::move(item));
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

py::dict evaluate_two_phase_flash_nlp(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::array<double, 2>& overall_mole_fractions,
    const std::array<double, 6>& variables,
    const std::string& expected_fingerprint
) {
    const epcsaft_equilibrium::ProviderContext provider(
        checked_mixture_sdk(capsule),
        expected_fingerprint
    );
    const epcsaft_equilibrium::FlashNlpEvaluation evaluation =
        epcsaft_equilibrium::evaluate_two_phase_flash_nlp(
            provider,
            temperature_k,
            pressure_pa,
            overall_mole_fractions,
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

py::dict phase_to_dict(const epcsaft_equilibrium::PhaseEvaluation& phase) {
    py::dict result;
    result["amount_mol"] = phase.amount_mol;
    result["volume_m3"] = phase.volume_m3;
    result["molar_density_mol_m3"] = phase.amount_mol / phase.volume_m3;
    result["pressure_pa"] = phase.provider.pressure_pa;
    result["chemical_potential_over_rt"] = phase.provider.chemical_potential_over_rt;
    return result;
}

py::dict flash_phase_to_dict(const epcsaft_equilibrium::FlashPhaseEvaluation& phase) {
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

py::list flash_attempt_log_to_list(
    const std::vector<epcsaft_equilibrium::FlashAttemptRecord>& attempts
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

py::dict solve_two_phase_flash(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::array<double, 2>& overall_mole_fractions
) {
    const epcsaft_equilibrium::ProviderContext provider(
        checked_mixture_sdk(capsule),
        std::string(kFlashFingerprint)
    );
    epcsaft_equilibrium::FlashSolveResult solve;
    {
        py::gil_scoped_release release;
        solve = epcsaft_equilibrium::solve_two_phase_flash(
            provider,
            temperature_k,
            pressure_pa,
            overall_mole_fractions
        );
    }

    py::dict diagnostics;
    diagnostics["solver_converged"] = solve.solver_converged;
    diagnostics["solver_status"] = solve.solver_status;
    diagnostics["iterations"] = solve.iterations;
    diagnostics["attempts"] = solve.attempts;
    diagnostics["attempt_log"] = flash_attempt_log_to_list(solve.attempt_log);
    diagnostics["solver_lower_bounds"] = solve.solver_lower_bounds;
    diagnostics["solver_upper_bounds"] = solve.solver_upper_bounds;
    diagnostics["solver_constraint_violation"] = solve.solver_constraint_violation;
    diagnostics["numerical_converged"] = solve.numerical_converged;
    diagnostics["confirmation_solves"] = solve.confirmation_solves;
    diagnostics["confirmation_max_difference"] = solve.confirmation_max_difference;
    diagnostics["physical_accepted"] = solve.physical_accepted;
    diagnostics["material_balance_max_abs"] = solve.material_balance_max_abs;
    diagnostics["pressure_stationarity_max_relative"] =
        solve.pressure_stationarity_max_relative;
    diagnostics["chemical_potential_max_abs"] = solve.chemical_potential_max_abs;
    diagnostics["kkt_stationarity_max_abs"] = solve.kkt_stationarity_max_abs;
    diagnostics["phase_density_distance"] = solve.phase_density_distance;
    diagnostics["equality_multipliers"] = solve.equality_multipliers;
    diagnostics["lower_bound_multipliers"] = solve.lower_bound_multipliers;
    diagnostics["upper_bound_multipliers"] = solve.upper_bound_multipliers;
    diagnostics["exact_derivatives"] = true;
    diagnostics["globality_certificate"] = false;
    diagnostics["failure_reason"] = solve.failure_reason;

    py::dict result;
    result["accepted"] = solve.accepted;
    result["temperature_k"] = temperature_k;
    result["pressure_pa"] = pressure_pa;
    result["overall_mole_fractions"] = overall_mole_fractions;
    result["parameter_fingerprint"] = std::string(kFlashFingerprint);
    result["diagnostics"] = diagnostics;
    if (solve.accepted) {
        result["liquid"] = flash_phase_to_dict(solve.evaluation.liquid);
        result["vapor"] = flash_phase_to_dict(solve.evaluation.vapor);
        result["liquid_phase_fraction"] =
            solve.evaluation.liquid.amounts_mol[0] + solve.evaluation.liquid.amounts_mol[1];
        result["vapor_phase_fraction"] =
            solve.evaluation.vapor.amounts_mol[0] + solve.evaluation.vapor.amounts_mol[1];
        result["total_free_energy_over_rt"] = solve.evaluation.objective;
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
        "evaluate_phase",
        &evaluate_phase,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("amount_mol"),
        py::arg("volume_m3"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "evaluate_two_phase_flash_nlp",
        &evaluate_two_phase_flash_nlp,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("overall_mole_fractions"),
        py::arg("variables"),
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
        "_solve_two_phase_flash",
        &solve_two_phase_flash,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("overall_mole_fractions")
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
