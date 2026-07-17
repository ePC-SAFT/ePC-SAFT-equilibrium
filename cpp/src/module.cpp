#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <epcsaft/native_sdk_v1.h>

#include "saturation.hpp"
#include "two_phase_flash.hpp"

namespace py = pybind11;

namespace {

constexpr std::string_view kFlashFingerprint =
    "sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286";
constexpr std::size_t kPureSdkTableSize = offsetof(epcsaft_native_sdk_v1, component_count);

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
    if (sdk.table_size < sizeof(epcsaft_native_sdk_v1)) {
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
    const bool has_mixture_tail = sdk.table_size >= sizeof(epcsaft_native_sdk_v1);
    py::dict result;
    result["capsule_name"] = EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME;
    result["abi_version"] = sdk.abi_version;
    result["table_size"] = sdk.table_size;
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
