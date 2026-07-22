#include "provider.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>

namespace epcsaft_equilibrium {
namespace {

constexpr std::size_t kNeutralReferenceSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, evaluate_neutral_reference)
    + sizeof(epcsaft_evaluate_neutral_reference_v1);
constexpr double kNeutralReferenceConvergenceErrorMax = 5.0e-5;
constexpr double kRankToleranceMultiplier = 4096.0;

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

std::size_t rank_of_rows(
    const std::vector<double>& values,
    std::size_t rows,
    std::size_t columns
) {
    if (rows == 0 || columns == 0 || values.size() != rows * columns) {
        return 0;
    }
    std::vector<double> matrix = values;
    double scale = 0.0;
    for (double value : matrix) {
        scale = std::max(scale, std::abs(value));
    }
    const double tolerance = kRankToleranceMultiplier * std::numeric_limits<double>::epsilon()
        * std::max(1.0, scale) * static_cast<double>(std::max(rows, columns));
    std::size_t pivot_row = 0;
    for (std::size_t column = 0; column < columns && pivot_row < rows; ++column) {
        std::size_t pivot = pivot_row;
        for (std::size_t row = pivot_row + 1; row < rows; ++row) {
            if (std::abs(matrix[row * columns + column])
                > std::abs(matrix[pivot * columns + column])) {
                pivot = row;
            }
        }
        if (std::abs(matrix[pivot * columns + column]) <= tolerance) {
            continue;
        }
        if (pivot != pivot_row) {
            for (std::size_t index = column; index < columns; ++index) {
                std::swap(
                    matrix[pivot * columns + index],
                    matrix[pivot_row * columns + index]
                );
            }
        }
        const double pivot_value = matrix[pivot_row * columns + column];
        for (std::size_t row = pivot_row + 1; row < rows; ++row) {
            const double factor = matrix[row * columns + column] / pivot_value;
            for (std::size_t index = column; index < columns; ++index) {
                matrix[row * columns + index] -= factor
                    * matrix[pivot_row * columns + index];
            }
        }
        ++pivot_row;
    }
    return pivot_row;
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

NeutralReferenceEvaluation ProviderContext::evaluate_neutral_reference(
    double temperature_k,
    double pressure_pa
) const {
    if (sdk_.table_size < kNeutralReferenceSdkTableSize
        || sdk_.neutral_reference_result_size
            != sizeof(epcsaft_neutral_reference_result_v1)
        || sdk_.evaluate_neutral_reference == nullptr
        || sdk_.component_count == 0
        || sdk_.neutral_reference_basis_row_count == 0
        || sdk_.neutral_reference_basis_row_count > sdk_.component_count
        || sdk_.component_ids == nullptr
        || sdk_.component_charges == nullptr) {
        throw std::invalid_argument("Provider neutral-reference ABI contract is incomplete");
    }
    if (!std::isfinite(temperature_k) || temperature_k <= 0.0
        || !std::isfinite(pressure_pa) || pressure_pa <= 0.0
        || !std::isfinite(sdk_.source_temperature_min_k)
        || !std::isfinite(sdk_.source_temperature_max_k)
        || sdk_.source_temperature_min_k <= 0.0
        || sdk_.source_temperature_max_k < sdk_.source_temperature_min_k
        || temperature_k < sdk_.source_temperature_min_k
        || temperature_k > sdk_.source_temperature_max_k) {
        throw std::invalid_argument("Provider neutral-reference source domain is incompatible");
    }

    const std::size_t component_count = sdk_.component_count;
    const std::size_t row_count = sdk_.neutral_reference_basis_row_count;
    NeutralReferenceEvaluation result;
    result.component_count = component_count;
    result.neutral_basis_row_count = row_count;
    result.neutral_basis.assign(row_count * component_count, 0.0);
    result.log_fugacity_contractions.assign(row_count, 0.0);
    result.reference_composition.assign(component_count, 0.0);

    epcsaft_neutral_reference_result_v1 native{};
    native.struct_size = sizeof(native);
    native.component_count = component_count;
    native.neutral_basis_row_count = row_count;
    native.neutral_basis_capacity = result.neutral_basis.size();
    native.contraction_capacity = result.log_fugacity_contractions.size();
    native.reference_composition_capacity = result.reference_composition.size();
    native.neutral_basis = result.neutral_basis.data();
    native.log_fugacity_contractions = result.log_fugacity_contractions.data();
    native.reference_composition = result.reference_composition.data();
    const int status = sdk_.evaluate_neutral_reference(
        sdk_.model_context,
        fingerprint_.c_str(),
        temperature_k,
        pressure_pa,
        &native
    );
    if (native.struct_size != sizeof(native)) {
        throw std::invalid_argument("Provider neutral-reference result size changed");
    }
    if (status != native.status) {
        throw std::runtime_error(
            "Provider neutral-reference evaluation returned inconsistent status values"
        );
    }
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1) {
        const std::string error = decode_provider_char_array(
            native.error,
            sizeof(native.error),
            "Provider neutral-reference error"
        );
        throw std::domain_error("Provider neutral-reference evaluation failed: " + error);
    }
    if (native.component_count != component_count
        || native.neutral_basis_row_count != row_count
        || native.neutral_basis_capacity != result.neutral_basis.size()
        || native.contraction_capacity != result.log_fugacity_contractions.size()
        || native.reference_composition_capacity != result.reference_composition.size()
        || native.neutral_basis != result.neutral_basis.data()
        || native.log_fugacity_contractions != result.log_fugacity_contractions.data()
        || native.reference_composition != result.reference_composition.data()) {
        throw std::invalid_argument("Provider neutral-reference result buffers changed");
    }
    if (native.derivative_availability
        != EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_NONE_V1) {
        throw std::invalid_argument(
            "Provider neutral-reference derivative availability is unsupported"
        );
    }
    for (double value : result.neutral_basis) {
        require_finite(value, "Provider neutral-reference basis");
    }
    for (double value : result.log_fugacity_contractions) {
        require_finite(value, "Provider neutral-reference contraction");
    }
    for (double value : result.reference_composition) {
        if (!std::isfinite(value) || value < 0.0) {
            throw std::invalid_argument(
                "Provider neutral-reference composition must be finite and nonnegative"
            );
        }
    }
    if (rank_of_rows(result.neutral_basis, row_count, component_count) != row_count) {
        throw std::invalid_argument("Provider neutral-reference basis is rank deficient");
    }
    double composition_sum = 0.0;
    double neutral_composition_sum = 0.0;
    for (std::size_t component = 0; component < component_count; ++component) {
        composition_sum += result.reference_composition[component];
        if (sdk_.component_charges[component] == 0) {
            neutral_composition_sum += result.reference_composition[component];
        } else if (result.reference_composition[component] != 0.0) {
            throw std::invalid_argument(
                "Provider neutral-reference composition must be salt-free"
            );
        }
    }
    const double composition_tolerance = kRankToleranceMultiplier
        * std::numeric_limits<double>::epsilon() * static_cast<double>(component_count);
    if (std::abs(composition_sum - 1.0) > composition_tolerance
        || neutral_composition_sum <= 0.0) {
        throw std::invalid_argument(
            "Provider neutral-reference composition is not a salt-free simplex"
        );
    }
    for (std::size_t row = 0; row < row_count; ++row) {
        double charge = 0.0;
        for (std::size_t component = 0; component < component_count; ++component) {
            charge += result.neutral_basis[row * component_count + component]
                * static_cast<double>(sdk_.component_charges[component]);
        }
        if (std::abs(charge) > composition_tolerance) {
            throw std::invalid_argument("Provider neutral-reference basis row is not charge neutral");
        }
    }
    if (native.temperature_k != temperature_k || native.pressure_pa != pressure_pa
        || !std::isfinite(native.solvent_molar_mass_kg_per_mol)
        || native.solvent_molar_mass_kg_per_mol <= 0.0
        || !std::isfinite(native.reference_amount_mol)
        || native.reference_amount_mol != 1.0
        || !std::isfinite(native.reference_number_density_mol_per_m3)
        || native.reference_number_density_mol_per_m3 != 1.0
        || !std::isfinite(native.reference_molality_mol_per_kg)
        || native.reference_molality_mol_per_kg <= 0.0
        || !std::isfinite(native.reference_convergence_error)
        || native.reference_convergence_error < 0.0
        || native.reference_convergence_error > kNeutralReferenceConvergenceErrorMax) {
        throw std::invalid_argument(
            "Provider neutral-reference evaluation returned invalid scalars"
        );
    }
    std::string parameter_fingerprint = decode_provider_char_array(
        native.parameter_fingerprint,
        sizeof(native.parameter_fingerprint),
        "Provider neutral-reference fingerprint"
    );
    if (parameter_fingerprint != fingerprint_) {
        throw std::invalid_argument("Provider neutral-reference fingerprint mismatch");
    }
    std::string basis_id = decode_provider_char_array(
        native.helmholtz_basis_id,
        sizeof(native.helmholtz_basis_id),
        "Provider neutral-reference basis"
    );
    if (basis_id != EPCSAFT_NATIVE_HELMHOLTZ_BASIS_ID_V1) {
        throw std::invalid_argument("Provider neutral-reference basis identity mismatch");
    }
    result.derivative_availability = native.derivative_availability;
    result.temperature_k = native.temperature_k;
    result.pressure_pa = native.pressure_pa;
    result.solvent_molar_mass_kg_per_mol = native.solvent_molar_mass_kg_per_mol;
    result.reference_amount_mol = native.reference_amount_mol;
    result.reference_number_density_mol_per_m3 = native.reference_number_density_mol_per_m3;
    result.reference_molality_mol_per_kg = native.reference_molality_mol_per_kg;
    result.reference_convergence_error = native.reference_convergence_error;
    result.basis_id = std::move(basis_id);
    result.parameter_fingerprint = std::move(parameter_fingerprint);
    return result;
}

const std::string& ProviderContext::fingerprint() const {
    return fingerprint_;
}

}  // namespace epcsaft_equilibrium
