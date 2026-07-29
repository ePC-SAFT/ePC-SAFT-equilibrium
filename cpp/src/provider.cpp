#include "provider.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <set>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace epcsaft_equilibrium {
namespace {

constexpr std::size_t kNeutralReferenceSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, evaluate_neutral_reference)
    + sizeof(epcsaft_evaluate_neutral_reference_v1);
constexpr std::size_t kCapabilitySdkTableSize =
    offsetof(epcsaft_native_sdk_v1, capabilities)
    + sizeof(const epcsaft_native_capability_descriptor_v1*);
constexpr std::size_t kNeutralReferenceDerivativeSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, evaluate_neutral_reference_derivatives)
    + sizeof(epcsaft_evaluate_neutral_reference_derivatives_v1);
constexpr std::size_t kReactingPhaseParameterSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, evaluate_reacting_phase_parameters)
    + sizeof(epcsaft_evaluate_reacting_phase_parameters_v1);
constexpr std::size_t kInversePackingGeometrySdkTableSize =
    offsetof(epcsaft_native_sdk_v1, evaluate_inverse_packing_geometry)
    + sizeof(epcsaft_evaluate_inverse_packing_geometry_v1);
constexpr double kNeutralReferenceConvergenceErrorMax = 5.0e-5;
constexpr double kNeutralReferenceRootResidualRelativeMax = 1.0e-9;
constexpr double kNeutralReferenceRootBracketRelativeMax = 1.0e-8;
constexpr double kNeutralReferenceRootDensityStepRelativeMax = 1.0e-3;

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

std::uint32_t parameter_family_code(std::string_view family) {
    if (family == "binary_interaction_kij") {
        return EPCSAFT_NATIVE_PARAMETER_FAMILY_BINARY_INTERACTION_KIJ_V1;
    }
    if (family == "binary_interaction_lij") {
        return EPCSAFT_NATIVE_PARAMETER_FAMILY_BINARY_INTERACTION_LIJ_V1;
    }
    if (family == "segment_count") {
        return EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_COUNT_V1;
    }
    if (family == "segment_diameter") {
        return EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_DIAMETER_V1;
    }
    if (family == "dispersion_energy_over_k") {
        return EPCSAFT_NATIVE_PARAMETER_FAMILY_DISPERSION_ENERGY_OVER_K_V1;
    }
    if (family == "born_diameter") {
        return EPCSAFT_NATIVE_PARAMETER_FAMILY_BORN_DIAMETER_V1;
    }
    if (family == "solvation_factor") {
        return EPCSAFT_NATIVE_PARAMETER_FAMILY_SOLVATION_FACTOR_V1;
    }
    if (family == "ion_fraction_suppression") {
        return EPCSAFT_NATIVE_PARAMETER_FAMILY_ION_FRACTION_SUPPRESSION_V1;
    }
    if (family == "association_energy_over_k") {
        return EPCSAFT_NATIVE_PARAMETER_FAMILY_ASSOCIATION_ENERGY_OVER_K_V1;
    }
    if (family == "association_volume") {
        return EPCSAFT_NATIVE_PARAMETER_FAMILY_ASSOCIATION_VOLUME_V1;
    }
    if (family == "ionic_region_relative_permittivity") {
        return EPCSAFT_NATIVE_PARAMETER_FAMILY_IONIC_REGION_RELATIVE_PERMITTIVITY_V1;
    }
    if (family == "relative_permittivity") {
        return EPCSAFT_NATIVE_PARAMETER_FAMILY_RELATIVE_PERMITTIVITY_V1;
    }
    throw std::invalid_argument("Provider active-parameter family is unsupported");
}

std::uint32_t parameter_identity_code(std::string_view identity) {
    if (identity == "unordered_component_pair") {
        return EPCSAFT_NATIVE_PARAMETER_IDENTITY_UNORDERED_COMPONENT_PAIR_V1;
    }
    if (identity == "component") {
        return EPCSAFT_NATIVE_PARAMETER_IDENTITY_COMPONENT_V1;
    }
    if (identity == "model") {
        return EPCSAFT_NATIVE_PARAMETER_IDENTITY_MODEL_V1;
    }
    throw std::invalid_argument("Provider active-parameter identity is unsupported");
}

std::uint32_t parameter_coordinate_kind(std::uint32_t family) {
    return family + (
        EPCSAFT_NATIVE_CAPABILITY_COORDINATE_BINARY_INTERACTION_KIJ_V1
        - EPCSAFT_NATIVE_PARAMETER_FAMILY_BINARY_INTERACTION_KIJ_V1
    );
}

std::int32_t component_index(
    const epcsaft_native_sdk_v1& sdk,
    const std::string& component_id
) {
    std::int32_t found = -1;
    for (std::size_t index = 0; index < sdk.component_count; ++index) {
        if (sdk.component_ids[index] != nullptr
            && component_id == sdk.component_ids[index]) {
            if (found >= 0) {
                throw std::invalid_argument(
                    "Provider component identity is ambiguous"
                );
            }
            found = static_cast<std::int32_t>(index);
        }
    }
    if (found < 0) {
        throw std::invalid_argument(
            "Provider active-parameter component identity is missing"
        );
    }
    return found;
}

std::vector<epcsaft_active_parameter_request_v1> native_active_requests(
    const ProviderActiveParameterSet& active_parameters
) {
    std::vector<epcsaft_active_parameter_request_v1> result(
        active_parameters.parameters.size()
    );
    for (std::size_t index = 0; index < result.size(); ++index) {
        const ProviderActiveParameter& source = active_parameters.parameters[index];
        result[index] = {
            sizeof(epcsaft_active_parameter_request_v1),
            source.family_code,
            source.identity_code,
            source.component_index,
            source.pair_component_index_a,
            source.pair_component_index_b,
            source.value,
            source.unit.c_str(),
        };
    }
    return result;
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

NeutralReferenceEvaluation checked_neutral_reference_value(
    const epcsaft_native_sdk_v1& sdk,
    const std::string& fingerprint,
    double temperature_k,
    double pressure_pa,
    epcsaft_neutral_reference_result_v1& native,
    std::vector<double> neutral_basis,
    std::vector<double> contractions,
    const std::vector<double>& reference_composition
) {
    NeutralReferenceEvaluation result;
    result.component_count = sdk.component_count;
    result.neutral_basis_row_count = sdk.neutral_reference_basis_row_count;
    if (native.struct_size != sizeof(native)
        || native.component_count != result.component_count
        || native.neutral_basis_row_count != result.neutral_basis_row_count
        || native.neutral_basis_capacity != neutral_basis.size()
        || native.contraction_capacity != contractions.size()
        || native.reference_composition_capacity != reference_composition.size()
        || native.neutral_basis != neutral_basis.data()
        || native.log_fugacity_contractions != contractions.data()
        || native.reference_composition != reference_composition.data()) {
        throw std::invalid_argument("Provider neutral-reference result contract changed");
    }
    if (native.derivative_availability
        != EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_NONE_V1) {
        throw std::invalid_argument(
            "Provider neutral-reference value derivative availability is invalid"
        );
    }
    if (!std::all_of(neutral_basis.begin(), neutral_basis.end(),
            [](double value) { return std::isfinite(value); })
        || !std::all_of(
            contractions.begin(),
            contractions.end(),
            [](double value) { return std::isfinite(value); }
        )) {
        throw std::invalid_argument("Provider neutral-reference values must be finite");
    }
    double composition_sum = 0.0;
    for (std::size_t component = 0; component < result.component_count; ++component) {
        const double value = reference_composition[component];
        if (!std::isfinite(value) || value < 0.0
            || (sdk.component_charges[component] != 0 && value != 0.0)) {
            throw std::invalid_argument(
                "Provider neutral-reference composition must be salt-free"
            );
        }
        composition_sum += value;
    }
    for (std::size_t row = 0; row < result.neutral_basis_row_count; ++row) {
        double charge = 0.0;
        for (std::size_t component = 0; component < result.component_count; ++component) {
            charge += neutral_basis[row * result.component_count + component]
                * static_cast<double>(sdk.component_charges[component]);
        }
        if (std::abs(charge) > 1.0e-12) {
            throw std::invalid_argument(
                "Provider neutral-reference basis row is not charge neutral"
            );
        }
    }
    if (std::abs(composition_sum - 1.0) > 1.0e-12
        || native.temperature_k != temperature_k
        || native.pressure_pa != pressure_pa
        || native.reference_amount_mol != 1.0
        || native.reference_number_density_mol_per_m3 != 1.0
        || !std::isfinite(native.solvent_molar_mass_kg_per_mol)
        || native.solvent_molar_mass_kg_per_mol <= 0.0
        || !std::isfinite(native.reference_molality_mol_per_kg)
        || native.reference_molality_mol_per_kg <= 0.0
        || !std::isfinite(native.reference_convergence_error)
        || native.reference_convergence_error < 0.0
        || native.reference_convergence_error > kNeutralReferenceConvergenceErrorMax) {
        throw std::invalid_argument("Provider neutral-reference scalar identity is invalid");
    }
    result.basis_id = decode_provider_char_array(
        native.helmholtz_basis_id,
        sizeof(native.helmholtz_basis_id),
        "Provider neutral-reference basis"
    );
    result.parameter_fingerprint = decode_provider_char_array(
        native.parameter_fingerprint,
        sizeof(native.parameter_fingerprint),
        "Provider neutral-reference fingerprint"
    );
    if (result.basis_id != EPCSAFT_NATIVE_HELMHOLTZ_BASIS_ID_V1
        || result.parameter_fingerprint != fingerprint) {
        throw std::invalid_argument("Provider neutral-reference identity mismatch");
    }
    result.neutral_basis = std::move(neutral_basis);
    result.log_fugacity_contractions = std::move(contractions);
    result.derivative_availability = native.derivative_availability;
    result.temperature_k = native.temperature_k;
    result.pressure_pa = native.pressure_pa;
    result.reference_convergence_error = native.reference_convergence_error;
    return result;
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

double ProviderContext::evaluate_electrolyte_value(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3
) const {
    if (sdk_.electrolyte_phase_value_result_size
            != sizeof(epcsaft_electrolyte_phase_value_result_v1)
        || sdk_.evaluate_electrolyte_phase_value == nullptr
        || amounts_mol.size() != sdk_.component_count) {
        throw std::invalid_argument(
            "provider capsule is missing the electrolyte value contract"
        );
    }
    epcsaft_electrolyte_phase_value_result_v1 result{};
    result.struct_size = sizeof(result);
    const int status = sdk_.evaluate_electrolyte_phase_value(
        sdk_.model_context,
        temperature_k,
        amounts_mol.data(),
        amounts_mol.size(),
        volume_m3,
        &result
    );
    if (status != result.status) {
        throw std::runtime_error(
            "provider electrolyte value returned inconsistent status values"
        );
    }
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1) {
        throw std::domain_error(
            "provider electrolyte value evaluation failed: "
            + decode_provider_char_array(
                result.error,
                sizeof(result.error),
                "provider electrolyte value error"
            )
        );
    }
    require_finite(
        result.helmholtz_over_rt_reference_amount,
        "provider electrolyte value"
    );
    if (decode_provider_char_array(
            result.parameter_fingerprint,
            sizeof(result.parameter_fingerprint),
            "provider electrolyte value fingerprint"
        ) != fingerprint_) {
        throw std::invalid_argument(
            "provider electrolyte value fingerprint does not match the requested model"
        );
    }
    return result.helmholtz_over_rt_reference_amount;
}

std::array<double, 2> ProviderContext::evaluate_molar_volume_bounds(
    double temperature_k,
    const std::vector<double>& mole_fractions,
    double packing_fraction_min,
    double packing_fraction_max
) const {
    if (sdk_.evaluate_molar_volume_bounds == nullptr
        || mole_fractions.size() != sdk_.component_count) {
        throw std::invalid_argument(
            "provider capsule is missing the molar-volume domain contract"
        );
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
    if (!std::all_of(
            mole_fractions.begin(),
            mole_fractions.end(),
            [](double value) {
                return std::isfinite(value) && value >= 0.0;
            }
        )) {
        throw std::invalid_argument(
            "mole fractions must be finite and nonnegative"
        );
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
    result.gradient.resize(coordinate_count);
    result.hessian.resize(coordinate_count * coordinate_count);
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

NeutralReferenceEvaluation ProviderContext::evaluate_neutral_reference(
    double temperature_k,
    double pressure_pa
) const {
    if (sdk_.table_size < kNeutralReferenceSdkTableSize
        || sdk_.neutral_reference_result_size
            != sizeof(epcsaft_neutral_reference_result_v1)
        || sdk_.evaluate_neutral_reference == nullptr
        || sdk_.neutral_reference_basis_row_count == 0
        || sdk_.neutral_reference_basis_row_count + 1 != sdk_.component_count) {
        throw std::invalid_argument("Provider neutral-reference ABI contract is incomplete");
    }
    require_finite(temperature_k, "neutral-reference temperature");
    require_finite(pressure_pa, "neutral-reference pressure");
    if (temperature_k <= 0.0 || pressure_pa <= 0.0
        || temperature_k < sdk_.source_temperature_min_k
        || temperature_k > sdk_.source_temperature_max_k) {
        throw std::invalid_argument("Provider neutral-reference source domain is incompatible");
    }

    NeutralReferenceEvaluation result;
    result.component_count = sdk_.component_count;
    result.neutral_basis_row_count = sdk_.neutral_reference_basis_row_count;
    result.neutral_basis.resize(
        result.component_count * result.neutral_basis_row_count
    );
    result.log_fugacity_contractions.resize(result.neutral_basis_row_count);
    std::vector<double> reference_composition(result.component_count, 0.0);

    epcsaft_neutral_reference_result_v1 native{};
    native.struct_size = sizeof(native);
    native.component_count = result.component_count;
    native.neutral_basis_row_count = result.neutral_basis_row_count;
    native.neutral_basis_capacity = result.neutral_basis.size();
    native.contraction_capacity = result.log_fugacity_contractions.size();
    native.reference_composition_capacity = reference_composition.size();
    native.neutral_basis = result.neutral_basis.data();
    native.log_fugacity_contractions = result.log_fugacity_contractions.data();
    native.reference_composition = reference_composition.data();
    const int status = sdk_.evaluate_neutral_reference(
        sdk_.model_context,
        fingerprint_.c_str(),
        temperature_k,
        pressure_pa,
        &native
    );
    if (status != native.status) {
        throw std::runtime_error(
            "Provider neutral-reference evaluation returned inconsistent status values"
        );
    }
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1) {
        throw std::domain_error(
            "Provider neutral-reference evaluation failed: "
            + decode_provider_char_array(
                native.error,
                sizeof(native.error),
                "Provider neutral-reference error"
            )
        );
    }
    return checked_neutral_reference_value(
        sdk_,
        fingerprint_,
        temperature_k,
        pressure_pa,
        native,
        std::move(result.neutral_basis),
        std::move(result.log_fugacity_contractions),
        reference_composition
    );
}

ProviderActiveParameterSet ProviderContext::resolve_active_parameters(
    double temperature_k,
    const std::vector<ProviderActiveParameterInput>& requests
) const {
    ProviderActiveParameterSet result;
    if (requests.empty()) {
        return result;
    }
    if (sdk_.table_size < kCapabilitySdkTableSize
        || sdk_.table_size < kReactingPhaseParameterSdkTableSize
        || sdk_.capability_count == 0
        || sdk_.capabilities == nullptr
        || sdk_.component_ids == nullptr
        || sdk_.reacting_phase_parameter_result_size
            != sizeof(epcsaft_reacting_phase_parameter_result_v1)
        || sdk_.evaluate_reacting_phase_parameters == nullptr) {
        throw std::invalid_argument(
            "Provider active-parameter ABI contract is incomplete"
        );
    }
    require_finite(temperature_k, "active-parameter temperature");
    std::set<std::string> selected_coordinates;
    result.parameters.reserve(requests.size());
    for (const ProviderActiveParameterInput& request : requests) {
        require_finite(request.value, "active-parameter value");
        if (request.family.empty() || request.identity.empty() || request.unit.empty()) {
            throw std::invalid_argument("Provider active-parameter request is incomplete");
        }
        ProviderActiveParameter parameter;
        parameter.family = request.family;
        parameter.identity = request.identity;
        parameter.component_ids = request.component_ids;
        parameter.value = request.value;
        parameter.unit = request.unit;
        parameter.family_code = parameter_family_code(parameter.family);
        parameter.identity_code = parameter_identity_code(parameter.identity);
        if (parameter.identity_code
            == EPCSAFT_NATIVE_PARAMETER_IDENTITY_COMPONENT_V1) {
            if (parameter.component_ids.size() != 1) {
                throw std::invalid_argument(
                    "Provider component parameter requires exactly one component id"
                );
            }
            parameter.component_index = component_index(
                sdk_, parameter.component_ids.front()
            );
        } else if (parameter.identity_code
            == EPCSAFT_NATIVE_PARAMETER_IDENTITY_UNORDERED_COMPONENT_PAIR_V1) {
            if (parameter.component_ids.size() != 2
                || parameter.component_ids[0] == parameter.component_ids[1]) {
                throw std::invalid_argument(
                    "Provider pair parameter requires two distinct component ids"
                );
            }
            parameter.pair_component_index_a = component_index(
                sdk_, parameter.component_ids[0]
            );
            parameter.pair_component_index_b = component_index(
                sdk_, parameter.component_ids[1]
            );
            if (parameter.pair_component_index_b
                < parameter.pair_component_index_a) {
                std::swap(
                    parameter.pair_component_index_a,
                    parameter.pair_component_index_b
                );
            }
        } else if (!parameter.component_ids.empty()) {
            throw std::invalid_argument(
                "Provider model parameter must not contain component ids"
            );
        }

        const epcsaft_native_capability_descriptor_v1* match = nullptr;
        for (std::size_t descriptor_index = 0;
             descriptor_index < sdk_.capability_count;
             ++descriptor_index) {
            const auto& descriptor = sdk_.capabilities[descriptor_index];
            if (descriptor.struct_size != sizeof(descriptor)
                || descriptor.schema_version
                    != EPCSAFT_NATIVE_CAPABILITY_SCHEMA_VERSION_V1
                || descriptor.model_domain
                    != EPCSAFT_NATIVE_MODEL_DOMAIN_REACTING_ELECTROLYTE_PHASE_V1
                || descriptor.tensor_layout
                    != EPCSAFT_NATIVE_TENSOR_LAYOUT_ROW_MAJOR_V1
                || descriptor.derivative_order < 2
                || descriptor.maturity
                    != EPCSAFT_NATIVE_CAPABILITY_DERIVATIVE_READY_V1
                || descriptor.authority_effect
                    != EPCSAFT_NATIVE_AUTHORITY_EFFECT_NONE_V1
                || descriptor.parameter_family != parameter.family_code
                || descriptor.parameter_identity != parameter.identity_code
                || descriptor.state_coordinate_count != sdk_.component_count + 1
                || descriptor.active_parameter_count != 1
                || descriptor.coordinate_count != 1
                || descriptor.coordinates == nullptr
                || descriptor.component_count != sdk_.component_count
                || descriptor.component_ids == nullptr
                || !std::isfinite(descriptor.temperature_min_k)
                || !std::isfinite(descriptor.temperature_max_k)
                || descriptor.temperature_min_k > descriptor.temperature_max_k
                || temperature_k < descriptor.temperature_min_k
                || temperature_k > descriptor.temperature_max_k) {
                continue;
            }
            bool identity_matches = true;
            for (std::size_t component = 0;
                 component < sdk_.component_count;
                 ++component) {
                identity_matches = identity_matches
                    && descriptor.component_ids[component] != nullptr
                    && sdk_.component_ids[component] != nullptr
                    && std::string_view(descriptor.component_ids[component])
                        == sdk_.component_ids[component];
            }
            const auto& active_coordinate = descriptor.coordinates[0];
            identity_matches = identity_matches
                && active_coordinate.struct_size == sizeof(active_coordinate)
                && active_coordinate.kind
                    == parameter_coordinate_kind(parameter.family_code)
                && active_coordinate.component_index == parameter.component_index
                && active_coordinate.pair_component_index_a
                    == parameter.pair_component_index_a
                && active_coordinate.pair_component_index_b
                    == parameter.pair_component_index_b
                && active_coordinate.unit != nullptr
                && parameter.unit == active_coordinate.unit;
            const std::string descriptor_fingerprint =
                decode_provider_char_array(
                    descriptor.parameter_fingerprint,
                    sizeof(descriptor.parameter_fingerprint),
                    "Provider active-parameter descriptor fingerprint"
                );
            const std::string descriptor_topology = decode_provider_char_array(
                descriptor.topology_fingerprint,
                sizeof(descriptor.topology_fingerprint),
                "Provider active-parameter descriptor topology"
            );
            const std::string descriptor_basis = decode_provider_char_array(
                descriptor.helmholtz_basis_id,
                sizeof(descriptor.helmholtz_basis_id),
                "Provider active-parameter descriptor basis"
            );
            identity_matches = identity_matches
                && descriptor_fingerprint == fingerprint_
                && !descriptor_topology.empty()
                && descriptor_basis == EPCSAFT_NATIVE_HELMHOLTZ_BASIS_ID_V1;
            if (!identity_matches) {
                continue;
            }
            if (match != nullptr) {
                throw std::invalid_argument(
                    "Provider active-parameter descriptor is ambiguous"
                );
            }
            match = &descriptor;
        }
        if (match == nullptr) {
            throw std::invalid_argument(
                "Provider active-parameter descriptor is missing or incompatible"
            );
        }
        const std::string topology = decode_provider_char_array(
            match->topology_fingerprint,
            sizeof(match->topology_fingerprint),
            "Provider active-parameter topology"
        );
        if (!result.topology_fingerprint.empty()
            && result.topology_fingerprint != topology) {
            throw std::invalid_argument(
                "Provider active-parameter requests have incompatible topologies"
            );
        }
        result.topology_fingerprint = topology;
        const std::string coordinate_key =
            std::to_string(parameter.family_code) + ":"
            + std::to_string(parameter.identity_code) + ":"
            + std::to_string(parameter.component_index) + ":"
            + std::to_string(parameter.pair_component_index_a) + ":"
            + std::to_string(parameter.pair_component_index_b);
        if (!selected_coordinates.insert(coordinate_key).second) {
            throw std::invalid_argument(
                "Provider active-parameter coordinate is duplicated"
            );
        }
        result.parameters.push_back(std::move(parameter));
    }
    return result;
}

ParameterizedPhaseEvaluation ProviderContext::evaluate_reacting_phase_parameters(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3,
    double packing_fraction_min,
    double packing_fraction_max,
    const ProviderActiveParameterSet& active_parameters
) const {
    if (active_parameters.parameters.empty()
        || active_parameters.topology_fingerprint.empty()
        || sdk_.table_size < kReactingPhaseParameterSdkTableSize
        || sdk_.reacting_phase_parameter_result_size
            != sizeof(epcsaft_reacting_phase_parameter_result_v1)
        || sdk_.evaluate_reacting_phase_parameters == nullptr
        || amounts_mol.size() != sdk_.component_count) {
        throw std::invalid_argument(
            "Provider atomic active reacting-phase contract is incomplete"
        );
    }
    const std::size_t state_count = amounts_mol.size() + 1;
    const std::size_t parameter_count = active_parameters.parameters.size();
    ParameterizedPhaseEvaluation result;
    result.phase.gradient.resize(state_count);
    result.phase.hessian.resize(state_count * state_count);
    result.packing.gradient.resize(state_count);
    result.packing.hessian.resize(state_count * state_count);
    result.state_parameter_derivatives.resize(state_count * parameter_count);
    std::vector<double> parameter_gradient(parameter_count);
    std::vector<double> packing_parameter_derivatives(parameter_count);
    std::vector<double> bound_derivatives(2 * parameter_count);
    result.pressure_parameter_derivatives_pa.resize(parameter_count);
    result.chemical_potential_parameter_derivatives_over_rt.resize(
        amounts_mol.size() * parameter_count
    );
    const std::vector<epcsaft_active_parameter_request_v1> native_requests =
        native_active_requests(active_parameters);
    epcsaft_reacting_phase_parameter_result_v1 native{};
    native.struct_size = sizeof(native);
    native.state_coordinate_count = state_count;
    native.active_parameter_count = parameter_count;
    native.state_gradient_capacity = result.phase.gradient.size();
    native.state_hessian_capacity = result.phase.hessian.size();
    native.parameter_gradient_capacity = parameter_gradient.size();
    native.state_parameter_capacity = result.state_parameter_derivatives.size();
    native.packing_parameter_capacity = packing_parameter_derivatives.size();
    native.molar_volume_bound_parameter_capacity = bound_derivatives.size();
    native.pressure_parameter_capacity =
        result.pressure_parameter_derivatives_pa.size();
    native.chemical_potential_parameter_capacity =
        result.chemical_potential_parameter_derivatives_over_rt.size();
    native.state_gradient = result.phase.gradient.data();
    native.state_hessian = result.phase.hessian.data();
    native.parameter_gradient = parameter_gradient.data();
    native.state_parameter_derivatives =
        result.state_parameter_derivatives.data();
    native.packing_parameter_derivatives =
        packing_parameter_derivatives.data();
    native.molar_volume_bound_parameter_derivatives_m3_per_mol =
        bound_derivatives.data();
    native.pressure_parameter_derivatives_pa =
        result.pressure_parameter_derivatives_pa.data();
    native.chemical_potential_parameter_derivatives_over_rt =
        result.chemical_potential_parameter_derivatives_over_rt.data();
    native.packing_state_gradient_capacity = result.packing.gradient.size();
    native.packing_state_hessian_capacity = result.packing.hessian.size();
    native.packing_state_gradient = result.packing.gradient.data();
    native.packing_state_hessian = result.packing.hessian.data();
    const int status = sdk_.evaluate_reacting_phase_parameters(
        sdk_.model_context,
        fingerprint_.c_str(),
        active_parameters.topology_fingerprint.c_str(),
        temperature_k,
        amounts_mol.data(),
        amounts_mol.size(),
        volume_m3,
        packing_fraction_min,
        packing_fraction_max,
        native_requests.data(),
        native_requests.size(),
        &native
    );
    if (status != native.status) {
        throw std::runtime_error(
            "Provider active reacting-phase evaluation returned inconsistent status values"
        );
    }
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1) {
        throw std::domain_error(
            "Provider active reacting-phase evaluation failed: "
            + decode_provider_char_array(
                native.error,
                sizeof(native.error),
                "Provider active reacting-phase error"
            )
        );
    }
    if (native.struct_size != sizeof(native)
        || native.state_coordinate_count != state_count
        || native.active_parameter_count != parameter_count
        || native.state_gradient_capacity != result.phase.gradient.size()
        || native.state_hessian_capacity != result.phase.hessian.size()
        || native.parameter_gradient_capacity != parameter_gradient.size()
        || native.state_parameter_capacity
            != result.state_parameter_derivatives.size()
        || native.packing_parameter_capacity
            != packing_parameter_derivatives.size()
        || native.molar_volume_bound_parameter_capacity
            != bound_derivatives.size()
        || native.pressure_parameter_capacity
            != result.pressure_parameter_derivatives_pa.size()
        || native.chemical_potential_parameter_capacity
            != result.chemical_potential_parameter_derivatives_over_rt.size()
        || native.packing_state_gradient_capacity
            != result.packing.gradient.size()
        || native.packing_state_hessian_capacity
            != result.packing.hessian.size()
        || native.state_gradient != result.phase.gradient.data()
        || native.state_hessian != result.phase.hessian.data()
        || native.parameter_gradient != parameter_gradient.data()
        || native.state_parameter_derivatives
            != result.state_parameter_derivatives.data()
        || native.packing_parameter_derivatives
            != packing_parameter_derivatives.data()
        || native.molar_volume_bound_parameter_derivatives_m3_per_mol
            != bound_derivatives.data()
        || native.pressure_parameter_derivatives_pa
            != result.pressure_parameter_derivatives_pa.data()
        || native.chemical_potential_parameter_derivatives_over_rt
            != result.chemical_potential_parameter_derivatives_over_rt.data()
        || native.packing_state_gradient != result.packing.gradient.data()
        || native.packing_state_hessian != result.packing.hessian.data()) {
        throw std::invalid_argument(
            "Provider active reacting-phase result buffer contract changed"
        );
    }
    const auto finite = [](const std::vector<double>& values) {
        return std::all_of(values.begin(), values.end(), [](double value) {
            return std::isfinite(value);
        });
    };
    if (!finite(result.phase.gradient) || !finite(result.phase.hessian)
        || !finite(result.packing.gradient) || !finite(result.packing.hessian)
        || !finite(parameter_gradient)
        || !finite(result.state_parameter_derivatives)
        || !finite(packing_parameter_derivatives) || !finite(bound_derivatives)
        || !finite(result.pressure_parameter_derivatives_pa)
        || !finite(result.chemical_potential_parameter_derivatives_over_rt)
        || !std::isfinite(native.helmholtz_over_rt_reference_amount)
        || !std::isfinite(native.pressure_pa)
        || !std::isfinite(native.packing_fraction)
        || !std::isfinite(native.molar_volume_bounds_m3_per_mol[0])
        || !std::isfinite(native.molar_volume_bounds_m3_per_mol[1])
        || native.molar_volume_bounds_m3_per_mol[0] <= 0.0
        || native.molar_volume_bounds_m3_per_mol[1]
            <= native.molar_volume_bounds_m3_per_mol[0]) {
        throw std::invalid_argument(
            "Provider active reacting-phase tensors or bounds are invalid"
        );
    }
    result.phase.value = native.helmholtz_over_rt_reference_amount;
    result.phase.pressure_pa = native.pressure_pa;
    result.packing.value = native.packing_fraction;
    result.phase.parameter_fingerprint = decode_provider_char_array(
        native.parameter_fingerprint,
        sizeof(native.parameter_fingerprint),
        "Provider active reacting-phase fingerprint"
    );
    result.topology_fingerprint = decode_provider_char_array(
        native.topology_fingerprint,
        sizeof(native.topology_fingerprint),
        "Provider active reacting-phase topology"
    );
    const std::string basis = decode_provider_char_array(
        native.helmholtz_basis_id,
        sizeof(native.helmholtz_basis_id),
        "Provider active reacting-phase basis"
    );
    if (result.phase.parameter_fingerprint != fingerprint_
        || result.topology_fingerprint
            != active_parameters.topology_fingerprint
        || basis != EPCSAFT_NATIVE_HELMHOLTZ_BASIS_ID_V1) {
        throw std::invalid_argument(
            "Provider active reacting-phase identity changed"
        );
    }
    result.molar_volume_bounds_m3_per_mol = {
        native.molar_volume_bounds_m3_per_mol[0],
        native.molar_volume_bounds_m3_per_mol[1],
    };
    constexpr double kProjectionRelativeTolerance = 2.0e-11;
    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        for (std::size_t species = 0; species < amounts_mol.size(); ++species) {
            const double mixed =
                result.state_parameter_derivatives[
                    species * parameter_count + parameter
                ];
            const double projected =
                result.chemical_potential_parameter_derivatives_over_rt[
                    species * parameter_count + parameter
                ];
            const double scale = std::max({1.0, std::abs(mixed), std::abs(projected)});
            if (std::abs(mixed - projected) > kProjectionRelativeTolerance * scale) {
                throw std::invalid_argument(
                    "Provider active chemical-potential projection is inconsistent"
                );
            }
        }
        const double projected_pressure =
            -8.31446261815324 * temperature_k
            * result.state_parameter_derivatives[
                amounts_mol.size() * parameter_count + parameter
            ];
        const double returned_pressure =
            result.pressure_parameter_derivatives_pa[parameter];
        const double scale = std::max(
            {1.0, std::abs(projected_pressure), std::abs(returned_pressure)}
        );
        if (std::abs(projected_pressure - returned_pressure)
            > kProjectionRelativeTolerance * scale) {
            throw std::invalid_argument(
                "Provider active pressure projection is inconsistent"
            );
        }
    }
    return result;
}

InversePackingGeometryEvaluation ProviderContext::evaluate_inverse_packing_geometry(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double log_packing_fraction,
    const ProviderActiveParameterSet& active_parameters
) const {
    if (sdk_.table_size < kInversePackingGeometrySdkTableSize
        || sdk_.inverse_packing_geometry_result_size
            != sizeof(epcsaft_inverse_packing_geometry_result_v1)
        || sdk_.evaluate_inverse_packing_geometry == nullptr
        || amounts_mol.size() != sdk_.component_count) {
        throw std::invalid_argument(
            "Provider inverse-packing geometry ABI contract is incomplete"
        );
    }
    require_finite(temperature_k, "inverse-packing temperature");
    require_finite(log_packing_fraction, "log packing fraction");
    if (temperature_k <= 0.0
        || !std::all_of(amounts_mol.begin(), amounts_mol.end(), [](double value) {
            return std::isfinite(value) && value > 0.0;
        })) {
        throw std::invalid_argument(
            "inverse-packing geometry state must be finite and positive"
        );
    }
    std::string topology = active_parameters.topology_fingerprint;
    if (topology.empty()) {
        if (sdk_.capability_count == 0 || sdk_.capabilities == nullptr
            || sdk_.table_size < kCapabilitySdkTableSize) {
            throw std::invalid_argument(
                "Provider inverse-packing topology capability is incomplete"
            );
        }
        std::size_t matching_descriptor_count = 0;
        for (std::size_t index = 0; index < sdk_.capability_count; ++index) {
            const auto& descriptor = sdk_.capabilities[index];
            if (descriptor.struct_size != sizeof(descriptor)
                || descriptor.schema_version
                    != EPCSAFT_NATIVE_CAPABILITY_SCHEMA_VERSION_V1
                || descriptor.model_domain
                    != EPCSAFT_NATIVE_MODEL_DOMAIN_REACTING_ELECTROLYTE_PHASE_V1
                || descriptor.tensor_layout
                    != EPCSAFT_NATIVE_TENSOR_LAYOUT_ROW_MAJOR_V1
                || descriptor.derivative_order < 2
                || descriptor.maturity
                    != EPCSAFT_NATIVE_CAPABILITY_DERIVATIVE_READY_V1
                || descriptor.authority_effect
                    != EPCSAFT_NATIVE_AUTHORITY_EFFECT_NONE_V1
                || descriptor.state_coordinate_count
                    != sdk_.component_count + 1
                || descriptor.coordinates == nullptr
                || descriptor.component_count != sdk_.component_count) {
                continue;
            }
            const std::string descriptor_fingerprint = decode_provider_char_array(
                descriptor.parameter_fingerprint,
                sizeof(descriptor.parameter_fingerprint),
                "Provider inverse-packing capability fingerprint"
            );
            if (descriptor_fingerprint != fingerprint_) {
                continue;
            }
            const std::string descriptor_topology = decode_provider_char_array(
                descriptor.topology_fingerprint,
                sizeof(descriptor.topology_fingerprint),
                "Provider inverse-packing capability topology"
            );
            if (descriptor_topology.empty()) {
                continue;
            }
            const std::string descriptor_basis = decode_provider_char_array(
                descriptor.helmholtz_basis_id,
                sizeof(descriptor.helmholtz_basis_id),
                "Provider inverse-packing capability basis"
            );
            if (descriptor_basis != EPCSAFT_NATIVE_HELMHOLTZ_BASIS_ID_V1
                || descriptor.component_ids == nullptr) {
                continue;
            }
            bool component_order_matches = true;
            for (std::size_t component = 0;
                 component < sdk_.component_count;
                 ++component) {
                component_order_matches = component_order_matches
                    && descriptor.component_ids[component] != nullptr
                    && sdk_.component_ids[component] != nullptr
                    && std::string_view(descriptor.component_ids[component])
                        == sdk_.component_ids[component];
            }
            if (!component_order_matches) {
                continue;
            }
            ++matching_descriptor_count;
            if (topology.empty()) {
                topology = descriptor_topology;
            } else if (topology != descriptor_topology) {
                throw std::invalid_argument(
                    "Provider inverse-packing capability topologies disagree"
                );
            }
        }
        if (topology.empty()) {
            throw std::invalid_argument(
                "Provider inverse-packing topology capability is missing"
            );
        }
        if (matching_descriptor_count == 0) {
            throw std::invalid_argument(
                "Provider inverse-packing topology capability is not derivative-ready"
            );
        }
    }
    const std::size_t coordinate_count = amounts_mol.size() + 1
        + active_parameters.parameters.size();
    InversePackingGeometryEvaluation result;
    result.gradient.resize(coordinate_count);
    result.hessian.resize(coordinate_count * coordinate_count);
    const std::vector<epcsaft_active_parameter_request_v1> native_requests =
        native_active_requests(active_parameters);
    epcsaft_inverse_packing_geometry_result_v1 native{};
    native.struct_size = sizeof(native);
    native.coordinate_count = coordinate_count;
    native.active_parameter_count = active_parameters.parameters.size();
    native.gradient_capacity = result.gradient.size();
    native.hessian_capacity = result.hessian.size();
    native.gradient = result.gradient.data();
    native.hessian = result.hessian.data();
    const int status = sdk_.evaluate_inverse_packing_geometry(
        sdk_.model_context,
        fingerprint_.c_str(),
        topology.c_str(),
        temperature_k,
        amounts_mol.data(),
        amounts_mol.size(),
        log_packing_fraction,
        native_requests.empty() ? nullptr : native_requests.data(),
        native_requests.size(),
        &native
    );
    if (status != native.status) {
        throw std::runtime_error(
            "Provider inverse-packing geometry returned inconsistent status values"
        );
    }
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1) {
        throw std::domain_error(
            "Provider inverse-packing geometry failed: "
            + decode_provider_char_array(
                native.error,
                sizeof(native.error),
                "Provider inverse-packing geometry error"
            )
        );
    }
    if (native.struct_size != sizeof(native)
        || native.coordinate_count != coordinate_count
        || native.active_parameter_count != active_parameters.parameters.size()
        || native.gradient_capacity != result.gradient.size()
        || native.hessian_capacity != result.hessian.size()
        || native.gradient != result.gradient.data()
        || native.hessian != result.hessian.data()) {
        throw std::invalid_argument(
            "Provider inverse-packing geometry result buffer contract changed"
        );
    }
    const std::string returned_parameter_fingerprint = decode_provider_char_array(
        native.parameter_fingerprint,
        sizeof(native.parameter_fingerprint),
        "Provider inverse-packing geometry fingerprint"
    );
    const std::string returned_topology = decode_provider_char_array(
        native.topology_fingerprint,
        sizeof(native.topology_fingerprint),
        "Provider inverse-packing geometry topology"
    );
    if (returned_parameter_fingerprint != fingerprint_ || returned_topology != topology) {
        throw std::invalid_argument(
            "Provider inverse-packing geometry identity changed"
        );
    }
    if (!std::isfinite(native.volume_m3) || native.volume_m3 <= 0.0
        || !std::all_of(result.gradient.begin(), result.gradient.end(),
            [](double value) { return std::isfinite(value); })
        || !std::all_of(result.hessian.begin(), result.hessian.end(),
            [](double value) { return std::isfinite(value); })) {
        throw std::invalid_argument(
            "Provider inverse-packing geometry tensors are invalid"
        );
    }
    result.volume_m3 = native.volume_m3;
    result.topology_fingerprint = returned_topology;
    return result;
}

NeutralReferenceEvaluation ProviderContext::evaluate_neutral_reference_derivatives(
    double temperature_k,
    double pressure_pa,
    const ProviderActiveParameterSet& active_parameters
) const {
    if (sdk_.table_size < kCapabilitySdkTableSize
        || sdk_.table_size < kNeutralReferenceDerivativeSdkTableSize
        || sdk_.capability_count == 0
        || sdk_.capabilities == nullptr
        || sdk_.neutral_reference_derivative_result_size
            != sizeof(epcsaft_neutral_reference_derivative_result_v1)
        || sdk_.evaluate_neutral_reference_derivatives == nullptr) {
        throw std::invalid_argument(
            "Provider neutral-reference derivative ABI contract is incomplete"
        );
    }
    const epcsaft_native_capability_descriptor_v1* topology_descriptor = nullptr;
    for (std::size_t index = 0; index < sdk_.capability_count; ++index) {
        const auto& descriptor = sdk_.capabilities[index];
        if (descriptor.struct_size == sizeof(descriptor)
            && descriptor.schema_version == EPCSAFT_NATIVE_CAPABILITY_SCHEMA_VERSION_V1
            && descriptor.model_domain
                == EPCSAFT_NATIVE_MODEL_DOMAIN_REACTING_ELECTROLYTE_PHASE_V1
            && descriptor.maturity
                == EPCSAFT_NATIVE_CAPABILITY_DERIVATIVE_READY_V1
            && descriptor.authority_effect
                == EPCSAFT_NATIVE_AUTHORITY_EFFECT_NONE_V1) {
            const std::string candidate_topology = decode_provider_char_array(
                descriptor.topology_fingerprint,
                sizeof(descriptor.topology_fingerprint),
                "Provider derivative topology fingerprint"
            );
            if (active_parameters.topology_fingerprint.empty()
                || candidate_topology
                    == active_parameters.topology_fingerprint) {
                topology_descriptor = &descriptor;
                break;
            }
        }
    }
    if (topology_descriptor == nullptr) {
        throw std::invalid_argument(
            "Provider neutral-reference derivative topology is not advertised"
        );
    }
    const std::string topology_fingerprint = decode_provider_char_array(
        topology_descriptor->topology_fingerprint,
        sizeof(topology_descriptor->topology_fingerprint),
        "Provider derivative topology fingerprint"
    );
    const std::string descriptor_fingerprint = decode_provider_char_array(
        topology_descriptor->parameter_fingerprint,
        sizeof(topology_descriptor->parameter_fingerprint),
        "Provider derivative parameter fingerprint"
    );
    const std::string descriptor_basis = decode_provider_char_array(
        topology_descriptor->helmholtz_basis_id,
        sizeof(topology_descriptor->helmholtz_basis_id),
        "Provider derivative basis"
    );
    if (topology_fingerprint.empty()
        || descriptor_fingerprint != fingerprint_
        || descriptor_basis != EPCSAFT_NATIVE_HELMHOLTZ_BASIS_ID_V1) {
        throw std::invalid_argument(
            "Provider neutral-reference derivative descriptor identity mismatch"
        );
    }

    const std::size_t component_count = sdk_.component_count;
    const std::size_t row_count = sdk_.neutral_reference_basis_row_count;
    std::vector<double> neutral_basis(component_count * row_count, 0.0);
    std::vector<double> contractions(row_count, 0.0);
    std::vector<double> reference_composition(component_count, 0.0);
    std::vector<double> pressure_derivatives(row_count, 0.0);
    const std::size_t parameter_count = active_parameters.parameters.size();
    std::vector<double> parameter_derivatives(row_count * parameter_count, 0.0);
    const std::vector<epcsaft_active_parameter_request_v1> native_requests =
        native_active_requests(active_parameters);
    epcsaft_neutral_reference_derivative_result_v1 native{};
    native.struct_size = sizeof(native);
    native.value.struct_size = sizeof(native.value);
    native.value.component_count = component_count;
    native.value.neutral_basis_row_count = row_count;
    native.value.neutral_basis_capacity = neutral_basis.size();
    native.value.contraction_capacity = contractions.size();
    native.value.reference_composition_capacity = reference_composition.size();
    native.value.neutral_basis = neutral_basis.data();
    native.value.log_fugacity_contractions = contractions.data();
    native.value.reference_composition = reference_composition.data();
    native.active_parameter_count = parameter_count;
    native.pressure_derivative_capacity = pressure_derivatives.size();
    native.parameter_derivative_capacity = parameter_derivatives.size();
    native.pressure_derivatives_per_pa = pressure_derivatives.data();
    native.parameter_derivatives = parameter_derivatives.empty()
        ? nullptr
        : parameter_derivatives.data();
    const int status = sdk_.evaluate_neutral_reference_derivatives(
        sdk_.model_context,
        fingerprint_.c_str(),
        topology_fingerprint.c_str(),
        temperature_k,
        pressure_pa,
        native_requests.empty() ? nullptr : native_requests.data(),
        native_requests.size(),
        &native
    );
    if (status != native.status || status != native.value.status) {
        throw std::runtime_error(
            "Provider neutral-reference derivative returned inconsistent status values"
        );
    }
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1) {
        throw std::domain_error(
            "Provider neutral-reference derivative evaluation failed: "
            + decode_provider_char_array(
                native.error,
                sizeof(native.error),
                "Provider neutral-reference derivative error"
            )
        );
    }
    if (native.struct_size != sizeof(native)
        || native.value.struct_size != sizeof(native.value)
        || native.active_parameter_count != parameter_count
        || native.pressure_derivative_capacity != pressure_derivatives.size()
        || native.parameter_derivative_capacity != parameter_derivatives.size()
        || native.pressure_derivatives_per_pa != pressure_derivatives.data()
        || native.parameter_derivatives
            != (parameter_derivatives.empty()
                ? nullptr
                : parameter_derivatives.data())
        || native.derivative_availability
            != (
                EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_PRESSURE_V1
                | (parameter_count == 0
                    ? EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_NONE_V1
                    : EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_PARAMETERS_V1)
            )) {
        throw std::invalid_argument(
            "Provider neutral-reference derivative result contract changed"
        );
    }

    NeutralReferenceEvaluation result = checked_neutral_reference_value(
        sdk_,
        fingerprint_,
        temperature_k,
        pressure_pa,
        native.value,
        std::move(neutral_basis),
        std::move(contractions),
        reference_composition
    );
    const std::string returned_fingerprint = decode_provider_char_array(
        native.parameter_fingerprint,
        sizeof(native.parameter_fingerprint),
        "Provider neutral-reference derivative fingerprint"
    );
    const std::string returned_topology = decode_provider_char_array(
        native.topology_fingerprint,
        sizeof(native.topology_fingerprint),
        "Provider neutral-reference derivative topology"
    );
    const std::string returned_basis = decode_provider_char_array(
        native.helmholtz_basis_id,
        sizeof(native.helmholtz_basis_id),
        "Provider neutral-reference derivative basis"
    );
    const std::string reference_branch = decode_provider_char_array(
        native.reference_branch,
        sizeof(native.reference_branch),
        "Provider neutral-reference derivative branch"
    );
    if (returned_fingerprint != fingerprint_
        || returned_topology != topology_fingerprint
        || returned_basis != EPCSAFT_NATIVE_HELMHOLTZ_BASIS_ID_V1
        || reference_branch.empty()
        || !std::isfinite(native.source_pressure_min_pa)
        || !std::isfinite(native.source_pressure_max_pa)
        || native.source_pressure_min_pa <= 0.0
        || native.source_pressure_min_pa > pressure_pa
        || native.source_pressure_max_pa < pressure_pa
        || native.source_pressure_max_pa < native.source_pressure_min_pa
        || !std::isfinite(native.maximum_root_residual_pa)
        || native.maximum_root_residual_pa < 0.0
        || native.maximum_root_residual_pa
            > kNeutralReferenceRootResidualRelativeMax * pressure_pa
        || !std::isfinite(native.minimum_pressure_density_derivative_pa_m3_per_mol)
        || native.minimum_pressure_density_derivative_pa_m3_per_mol <= 0.0
        || !std::isfinite(native.maximum_density_condition_number)
        || native.maximum_density_condition_number < 0.0
        || !std::isfinite(native.reference_derivative_convergence_error)
        || native.reference_derivative_convergence_error < 0.0
        || native.reference_derivative_convergence_error
            > kNeutralReferenceConvergenceErrorMax
        || !std::isfinite(native.maximum_relative_root_bracket_width)
        || native.maximum_relative_root_bracket_width < 0.0
        || native.maximum_relative_root_bracket_width
            > kNeutralReferenceRootBracketRelativeMax
        || !std::isfinite(native.maximum_relative_root_density_step)
        || native.maximum_relative_root_density_step < 0.0
        || native.maximum_relative_root_density_step
            > kNeutralReferenceRootDensityStepRelativeMax
        || native.stable_root_count == 0
        || native.selected_stable_root_index >= native.stable_root_count
        || !std::all_of(
            pressure_derivatives.begin(),
            pressure_derivatives.end(),
            [](double value) { return std::isfinite(value); }
        )
        || !std::all_of(
            parameter_derivatives.begin(),
            parameter_derivatives.end(),
            [](double value) { return std::isfinite(value); }
        )) {
        throw std::invalid_argument(
            "Provider neutral-reference derivative certificate is invalid"
        );
    }
    result.pressure_derivatives_per_pa = std::move(pressure_derivatives);
    result.parameter_derivatives = std::move(parameter_derivatives);
    result.active_parameter_count = parameter_count;
    result.derivative_availability = native.derivative_availability;
    result.source_pressure_min_pa = native.source_pressure_min_pa;
    result.source_pressure_max_pa = native.source_pressure_max_pa;
    result.maximum_root_residual_pa = native.maximum_root_residual_pa;
    result.minimum_pressure_density_derivative_pa_m3_per_mol =
        native.minimum_pressure_density_derivative_pa_m3_per_mol;
    result.maximum_density_condition_number =
        native.maximum_density_condition_number;
    result.reference_derivative_convergence_error =
        native.reference_derivative_convergence_error;
    result.maximum_relative_root_bracket_width =
        native.maximum_relative_root_bracket_width;
    result.maximum_relative_root_density_step =
        native.maximum_relative_root_density_step;
    result.stable_root_count = native.stable_root_count;
    result.selected_stable_root_index = native.selected_stable_root_index;
    result.reference_branch = std::move(reference_branch);
    result.topology_fingerprint = std::move(returned_topology);
    return result;
}

const std::string& ProviderContext::fingerprint() const {
    return fingerprint_;
}

const epcsaft_native_sdk_v1& ProviderContext::sdk() const noexcept {
    return sdk_;
}

}  // namespace epcsaft_equilibrium
