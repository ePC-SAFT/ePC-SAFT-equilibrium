#include "provider.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <utility>

namespace epcsaft_equilibrium {
namespace {

std::string decode_provider_char_array(
    const char* value,
    std::size_t capacity,
    const char* field_name
) {
    const void* terminator = std::memchr(value, '\0', capacity);
    if (terminator == nullptr) {
        throw std::invalid_argument(std::string(field_name) + " is missing a NUL terminator");
    }
    const auto* end = static_cast<const char*>(terminator);
    return std::string(value, static_cast<std::size_t>(end - value));
}

void require_finite(double value, const char* name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
}

MixturePhaseEvaluation evaluate_mixture_callback(
    const epcsaft_native_sdk_v1& sdk,
    epcsaft_evaluate_mixture_phase_v1 callback,
    const std::string& fingerprint,
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3
) {
    if (sdk.component_count < 2 || amounts_mol.size() != sdk.component_count) {
        throw std::invalid_argument("provider mixture component count mismatch");
    }
    if (sdk.mixture_result_size != sizeof(epcsaft_mixture_phase_block_result_v1)
        || callback == nullptr) {
        throw std::invalid_argument("provider capsule is missing the mixture phase contract");
    }
    const std::size_t coordinate_count = sdk.component_count + 1;
    std::vector<double> gradient(coordinate_count, 0.0);
    std::vector<double> hessian(coordinate_count * coordinate_count, 0.0);
    epcsaft_mixture_phase_block_result_v1 phase{};
    phase.struct_size = sizeof(phase);
    phase.coordinate_count = coordinate_count;
    phase.gradient_capacity = gradient.size();
    phase.hessian_capacity = hessian.size();
    phase.gradient = gradient.data();
    phase.hessian = hessian.data();
    const int status = callback(
        sdk.model_context,
        temperature_k,
        amounts_mol.data(),
        amounts_mol.size(),
        volume_m3,
        &phase
    );
    if (phase.struct_size != sizeof(epcsaft_mixture_phase_block_result_v1)) {
        throw std::invalid_argument("provider mixture result struct size mismatch");
    }
    if (status != phase.status) {
        throw std::runtime_error("provider mixture evaluation returned inconsistent status values");
    }
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1) {
        const std::string error = decode_provider_char_array(
            phase.error,
            sizeof(phase.error),
            "provider mixture error"
        );
        throw std::domain_error("provider mixture phase evaluation failed: " + error);
    }
    if (phase.coordinate_count != coordinate_count
        || phase.gradient_capacity != gradient.size()
        || phase.hessian_capacity != hessian.size()
        || phase.gradient != gradient.data()
        || phase.hessian != hessian.data()) {
        throw std::invalid_argument("provider mixture result buffer contract changed");
    }
    require_finite(phase.helmholtz_over_rt_reference_amount, "provider mixture value");
    require_finite(phase.pressure_pa, "provider mixture pressure");
    if (!std::all_of(gradient.begin(), gradient.end(), [](double value) {
            return std::isfinite(value);
        })
        || !std::all_of(hessian.begin(), hessian.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::invalid_argument("provider mixture tensors must be finite");
    }
    const std::string parameter_fingerprint = decode_provider_char_array(
        phase.parameter_fingerprint,
        sizeof(phase.parameter_fingerprint),
        "provider mixture parameter fingerprint"
    );
    if (parameter_fingerprint != fingerprint) {
        throw std::invalid_argument(
            "provider mixture result fingerprint does not match the requested model"
        );
    }
    return {
        phase.helmholtz_over_rt_reference_amount,
        std::move(gradient),
        std::move(hessian),
        phase.pressure_pa,
        std::move(parameter_fingerprint),
    };
}

}  // namespace

ProviderContext::ProviderContext(const epcsaft_native_sdk_v1& sdk, std::string fingerprint)
    : sdk_(sdk), fingerprint_(std::move(fingerprint)) {
    if (fingerprint_.empty()) {
        throw std::invalid_argument("expected provider fingerprint must not be empty");
    }
}

PhaseEvaluation ProviderContext::evaluate(
    double temperature_k,
    double amount_mol,
    double volume_m3
) const {
    epcsaft_phase_block_result_v1 phase{};
    phase.struct_size = sizeof(phase);
    const int status = sdk_.evaluate_pure_phase(
        sdk_.model_context,
        temperature_k,
        amount_mol,
        volume_m3,
        &phase
    );
    if (phase.struct_size != sizeof(epcsaft_phase_block_result_v1)) {
        throw std::invalid_argument("provider phase result struct size mismatch");
    }
    if (status != phase.status) {
        throw std::runtime_error("provider phase evaluation returned inconsistent status values");
    }
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1) {
        const std::string error = decode_provider_char_array(
            phase.error,
            sizeof(phase.error),
            "provider error"
        );
        throw std::domain_error("provider phase evaluation failed: " + error);
    }
    const std::string parameter_fingerprint = decode_provider_char_array(
        phase.parameter_fingerprint,
        sizeof(phase.parameter_fingerprint),
        "provider parameter fingerprint"
    );
    if (parameter_fingerprint != fingerprint_) {
        throw std::invalid_argument(
            "provider phase result fingerprint does not match the requested model"
        );
    }
    return {amount_mol, volume_m3, phase, std::move(parameter_fingerprint)};
}

MixturePhaseEvaluation ProviderContext::evaluate_mixture(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3
) const {
    return evaluate_mixture_callback(
        sdk_, sdk_.evaluate_mixture_phase, fingerprint_, temperature_k, amounts_mol, volume_m3
    );
}

MixturePhaseEvaluation ProviderContext::evaluate_electrolyte(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3
) const {
    return evaluate_mixture_callback(
        sdk_, sdk_.evaluate_electrolyte_phase, fingerprint_, temperature_k, amounts_mol, volume_m3
    );
}

std::array<double, 2> ProviderContext::evaluate_molar_volume_bounds(
    double temperature_k,
    const std::vector<double>& mole_fractions,
    double packing_fraction_min,
    double packing_fraction_max
) const {
    if (sdk_.evaluate_molar_volume_bounds == nullptr
        || mole_fractions.size() != sdk_.component_count) {
        throw std::invalid_argument("provider capsule is missing the molar-volume domain contract");
    }
    require_finite(temperature_k, "temperature");
    require_finite(packing_fraction_min, "minimum packing fraction");
    require_finite(packing_fraction_max, "maximum packing fraction");
    if (temperature_k <= 0.0 || packing_fraction_min <= 0.0
        || packing_fraction_max <= packing_fraction_min) {
        throw std::invalid_argument(
            "temperature and packing-fraction bounds must be positive and ordered"
        );
    }
    if (!std::all_of(mole_fractions.begin(), mole_fractions.end(), [](double value) {
            return std::isfinite(value) && value >= 0.0;
        })) {
        throw std::invalid_argument("mole fractions must be finite and nonnegative");
    }
    std::array<double, 2> bounds{};
    const int status = sdk_.evaluate_molar_volume_bounds(
        sdk_.model_context,
        fingerprint_.c_str(),
        temperature_k,
        mole_fractions.data(),
        mole_fractions.size(),
        packing_fraction_min,
        packing_fraction_max,
        bounds.data(),
        bounds.size()
    );
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1) {
        throw std::domain_error(
            "provider molar-volume domain/fingerprint evaluation failed with status "
            + std::to_string(status)
        );
    }
    if (!std::isfinite(bounds[0]) || !std::isfinite(bounds[1])
        || bounds[0] <= 0.0 || bounds[1] <= bounds[0]) {
        throw std::invalid_argument(
            "provider molar-volume bounds must be finite, positive, and ordered"
        );
    }
    return bounds;
}

PackingFractionEvaluation ProviderContext::evaluate_packing_fraction(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3
) const {
    if (sdk_.evaluate_packing_fraction == nullptr
        || amounts_mol.size() != sdk_.component_count) {
        throw std::invalid_argument("provider capsule is missing the packing-fraction contract");
    }
    require_finite(temperature_k, "temperature");
    require_finite(volume_m3, "volume");
    if (temperature_k <= 0.0 || volume_m3 <= 0.0
        || !std::all_of(amounts_mol.begin(), amounts_mol.end(), [](double value) {
            return std::isfinite(value) && value >= 0.0;
        })) {
        throw std::invalid_argument(
            "packing-fraction state must be finite, nonnegative, and positive-volume"
        );
    }
    const std::size_t coordinate_count = amounts_mol.size() + 1;
    PackingFractionEvaluation result;
    result.gradient.assign(coordinate_count, 0.0);
    result.hessian.assign(coordinate_count * coordinate_count, 0.0);
    const int status = sdk_.evaluate_packing_fraction(
        sdk_.model_context,
        fingerprint_.c_str(),
        temperature_k,
        amounts_mol.data(),
        amounts_mol.size(),
        volume_m3,
        &result.value,
        result.gradient.data(),
        result.gradient.size(),
        result.hessian.data(),
        result.hessian.size()
    );
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1) {
        throw std::domain_error(
            "provider packing-fraction evaluation failed with status "
            + std::to_string(status)
        );
    }
    require_finite(result.value, "provider packing fraction");
    if (!std::all_of(result.gradient.begin(), result.gradient.end(), [](double value) {
            return std::isfinite(value);
        })
        || !std::all_of(result.hessian.begin(), result.hessian.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::invalid_argument("provider packing-fraction tensors must be finite");
    }
    return result;
}

const std::string& ProviderContext::fingerprint() const {
    return fingerprint_;
}

}  // namespace epcsaft_equilibrium
