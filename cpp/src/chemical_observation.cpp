#include "chemical_observation.hpp"

#include "chemical_equilibrium.hpp"
#include "provider.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <pybind11/stl.h>
#include <epcsaft/native_sdk_v1.h>
#include <epcsaft/regression/evaluator_v1.h>

namespace py = pybind11;

namespace epcsaft_equilibrium {
namespace {

constexpr std::size_t kChemicalSourceDomainSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, total_ion_mole_fraction_max) + sizeof(double);
struct ChemicalProviderMetadata {
    std::vector<std::string> component_ids;
    std::vector<int> charges;
};

DenseMatrix observation_dense_matrix(const py::handle& value, const char* field) {
    const auto rows = py::cast<std::vector<std::vector<double>>>(value);
    if (rows.empty()) return {};
    const auto columns = rows.front().size();
    if (columns == 0) throw py::value_error(std::string(field) + " rows must not be empty");
    DenseMatrix result{rows.size(), columns, {}};
    result.values.reserve(rows.size() * columns);
    for (const auto& row : rows) {
        if (row.size() != columns) throw py::value_error(std::string(field) + " must be rectangular");
        result.values.insert(result.values.end(), row.begin(), row.end());
    }
    return result;
}

std::vector<EquilibriumConstantRecord> observation_records(const py::handle& value) {
    std::vector<EquilibriumConstantRecord> result;
    for (const auto item : py::cast<py::tuple>(value)) {
        const auto record = py::cast<py::dict>(item);
        result.push_back({
            py::cast<std::string>(record["source_id"]),
            py::cast<std::string>(record["reference_id"]),
            py::cast<std::string>(record["reaction_orientation"]),
            py::cast<std::string>(record["conversion_id"]),
            py::cast<bool>(record["dimensionless"]),
            py::cast<double>(record["temperature_k"]),
            py::cast<double>(record["pressure_pa"]),
        });
    }
    return result;
}

ReactionSystemInput observation_input(const py::dict& spec) {
    ReactionSystemInput input;
    input.species_ids = py::cast<std::vector<std::string>>(spec["species_ids"]);
    input.charges = py::cast<std::vector<int>>(spec["charges"]);
    input.provider_fingerprint = py::cast<std::string>(spec["provider_fingerprint"]);
    input.molar_masses_kg_per_mol = py::cast<std::vector<double>>(spec["molar_masses_kg_per_mol"]);
    input.balance_matrix = observation_dense_matrix(spec["balance_matrix"], "balance matrix");
    input.conserved_totals = py::cast<std::vector<double>>(spec["conserved_totals"]);
    input.reaction_matrix = observation_dense_matrix(spec["reaction_matrix"], "reaction matrix");
    input.feed_amounts = py::cast<std::vector<double>>(spec["feed_amounts"]);
    input.ln_k = py::cast<std::vector<double>>(spec["ln_k"]);
    input.equilibrium_constant_records = observation_records(spec["equilibrium_constant_records"]);
    input.temperature_k = py::cast<double>(spec["temperature_k"]);
    input.pressure_pa = py::cast<double>(spec["pressure_pa"]);
    return input;
}

void validate_observation_identity(
    const ChemicalProviderMetadata& metadata, const ReactionSystemInput& input
) {
    if (metadata.component_ids != input.species_ids || metadata.charges != input.charges) {
        throw py::value_error("problem component identity does not match installed Provider");
    }
    if (input.provider_fingerprint.empty()) {
        throw py::value_error("problem Provider fingerprint is empty");
    }
}

const epcsaft_native_sdk_v1& checked_chemical_sdk(const py::capsule& capsule) {
    const char* name = capsule.name();
    if (name == nullptr || std::string_view(name) != EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME) {
        throw py::value_error("expected capsule epcsaft.native_sdk.v1");
    }
    const auto* sdk = capsule.get_pointer<epcsaft_native_sdk_v1>();
    if (sdk == nullptr || sdk->abi_version != EPCSAFT_NATIVE_SDK_V1_ABI_VERSION
        || sdk->table_size < kChemicalSourceDomainSdkTableSize
        || sdk->model_context == nullptr || sdk->component_count < 3
        || sdk->component_ids == nullptr || sdk->component_charges == nullptr
        || sdk->mixture_result_size != sizeof(epcsaft_mixture_phase_block_result_v1)
        || sdk->evaluate_electrolyte_phase == nullptr
        || sdk->evaluate_molar_volume_bounds == nullptr
        || sdk->evaluate_packing_fraction == nullptr) {
        throw py::value_error("Provider capsule is missing the reacting-phase SDK contract");
    }
    return *sdk;
}

ChemicalProviderMetadata chemical_provider_metadata(const epcsaft_native_sdk_v1& sdk) {
    ChemicalProviderMetadata result;
    result.component_ids.reserve(sdk.component_count);
    result.charges.reserve(sdk.component_count);
    for (std::size_t component = 0; component < sdk.component_count; ++component) {
        if (sdk.component_ids[component] == nullptr || sdk.component_ids[component][0] == '\0') {
            throw py::value_error("Provider component identity is incomplete");
        }
        result.component_ids.emplace_back(sdk.component_ids[component]);
        result.charges.push_back(static_cast<int>(sdk.component_charges[component]));
    }
    return result;
}

void validate_provider_identity(
    const ChemicalProviderMetadata& metadata, const ReactionSystemInput& input
) {
    validate_observation_identity(metadata, input);
}

constexpr char kChemicalObservationCapsuleName[] =
    "epcsaft.regression.evaluator.v1";
constexpr double kUniversalGasConstant = 8.31446261815324;

static_assert(
    std::string_view(EPCSAFT_REGRESSION_EVALUATOR_V1_SOLVER_STATUS_SOLVE_SUCCEEDED)
        == "solve_succeeded"
);
static_assert(
    std::string_view(EPCSAFT_REGRESSION_EVALUATOR_V1_NUMERICAL_STATUS_PASSED)
        == "passed"
);
static_assert(
    std::string_view(EPCSAFT_REGRESSION_EVALUATOR_V1_PHYSICAL_STATUS_PASSED)
        == "passed"
);
static_assert(
    std::string_view(EPCSAFT_REGRESSION_EVALUATOR_V1_DERIVATIVE_STATUS_AVAILABLE)
        == "available"
);

struct ChemicalObservationPrimitive {
    std::string kind;
    std::string component_id;
    std::size_t component_index = 0;
};

struct ChemicalObservationRow {
    std::string row_id;
    std::string state_id;
    std::string state_schema_id;
    std::string source_id;
    std::string transform_id;
    std::string reference_id;
    std::string reference_fingerprint;
    ChemicalObservationPrimitive primitive;
    ReactionSystemInput input;
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
};

struct ChemicalObservationContext {
    epcsaft_regression_evaluator_sdk_v1 sdk{};
    const epcsaft_native_sdk_v1* provider_sdk = nullptr;
    py::capsule provider_capsule;
    std::vector<ChemicalObservationRow> rows;
    std::vector<double> packing_bounds;
    double trace_floor = 0.0;
    std::vector<ProviderActiveParameterInput> parameter_templates;
    std::string provider_fingerprint;
    std::string evaluator_identity;
    std::string capability_id;
    std::string capability_fingerprint;
    std::string provider_artifact_identity;
    std::string owner_artifact_identity;
    std::string contract_fingerprint;
    std::string expected_provider_topology_fingerprint;
    std::vector<std::string> primitive_ids;
    std::vector<const char*> primitive_id_ptrs;
    std::vector<std::string> primitive_units;
    std::vector<const char*> primitive_unit_ptrs;
    std::vector<std::string> transform_ids;
    std::vector<const char*> transform_id_ptrs;
    std::vector<std::string> reference_ids;
    std::vector<const char*> reference_id_ptrs;
    std::vector<std::string> reference_fingerprints;
    std::vector<const char*> reference_fingerprint_ptrs;
    std::vector<std::string> row_ids;
    std::vector<const char*> row_id_ptrs;
    std::vector<std::string> state_ids;
    std::vector<const char*> state_id_ptrs;
    std::vector<std::string> state_schema_ids;
    std::vector<const char*> state_schema_id_ptrs;
    std::vector<std::string> source_ids;
    std::vector<const char*> source_id_ptrs;
    std::vector<std::string> parameter_ids;
    std::vector<const char*> parameter_id_ptrs;
    std::vector<std::string> parameter_units;
    std::vector<const char*> parameter_unit_ptrs;
    std::string artifact_identity;
};

void destroy_chemical_observation_context(PyObject* capsule) {
    void* pointer = PyCapsule_GetPointer(
        capsule, kChemicalObservationCapsuleName
    );
    if (pointer != nullptr) {
        auto* sdk = static_cast<epcsaft_regression_evaluator_sdk_v1*>(
            pointer
        );
        delete static_cast<ChemicalObservationContext*>(sdk->model_context);
    }
}

epcsaft_regression_evaluator_sdk_v1& checked_chemical_observation_sdk(
    const py::capsule& capsule
) {
    const char* name = capsule.name();
    if (name == nullptr || std::string_view(name) != kChemicalObservationCapsuleName) {
        throw py::value_error(
            "expected capsule epcsaft.regression.evaluator.v1"
        );
    }
    auto* sdk = capsule.get_pointer<epcsaft_regression_evaluator_sdk_v1>();
    if (sdk == nullptr
        || sdk->abi_version != EPCSAFT_REGRESSION_EVALUATOR_V1_ABI_VERSION
        || sdk->table_size < sizeof(epcsaft_regression_evaluator_sdk_v1)
        || sdk->model_context == nullptr
        || sdk->evaluate == nullptr
        || sdk->row_ids == nullptr
        || sdk->state_ids == nullptr
        || sdk->state_schema_ids == nullptr
        || sdk->observation_source_ids == nullptr
        || sdk->primitive_ids == nullptr
        || sdk->primitive_units == nullptr
        || sdk->transform_ids == nullptr
        || sdk->reference_ids == nullptr
        || sdk->reference_fingerprints == nullptr
        || sdk->parameter_ids == nullptr
        || sdk->parameter_units == nullptr
        || sdk->evaluator_identity == nullptr
        || sdk->capability_id == nullptr
        || sdk->capability_fingerprint == nullptr
        || sdk->provider_artifact_identity == nullptr
        || sdk->owner_artifact_identity == nullptr
        || sdk->contract_fingerprint == nullptr
        || sdk->model_fingerprint == nullptr
        || sdk->provider_parameter_fingerprint == nullptr
        || sdk->expected_provider_topology_fingerprint == nullptr
        || sdk->provider_sdk_capsule_name == nullptr
        || sdk->provider_sdk_abi_version != EPCSAFT_NATIVE_SDK_V1_ABI_VERSION
        || sdk->provider_sdk_table_size == 0
        || sdk->provider_sdk_result_size == 0
        || sdk->provider_sdk_mixture_result_size == 0
        || sdk->provider_sdk_reacting_phase_parameter_result_size == 0
        || sdk->artifact_identity == nullptr
        || sdk->result_size != sizeof(epcsaft_regression_evaluator_result_v1)
        || sdk->row_result_size != sizeof(epcsaft_regression_evaluator_row_result_v1)) {
        throw py::value_error("chemical-observation capsule has an invalid SDK table");
    }
    return *sdk;
}

std::string provider_parameter_name(const ProviderActiveParameter& parameter) {
    std::string name = "provider_parameter[" + parameter.family + ";"
        + parameter.identity + ";";
    for (std::size_t index = 0; index < parameter.component_ids.size(); ++index) {
        if (index != 0) name += ",";
        name += parameter.component_ids[index];
    }
    return name + "]";
}

std::vector<ProviderActiveParameterInput> observation_parameter_templates(
    const py::tuple& values
) {
    std::vector<ProviderActiveParameterInput> result;
    result.reserve(values.size());
    for (const py::handle item : values) {
        const py::dict request = py::cast<py::dict>(item);
        ProviderActiveParameterInput parameter{
            py::cast<std::string>(request["family"]),
            py::cast<std::string>(request["identity"]),
            py::cast<std::vector<std::string>>(request["component_ids"]),
            py::cast<double>(request["value"]),
            py::cast<std::string>(request["unit"]),
        };
        if (parameter.family.empty() || parameter.identity.empty()
            || parameter.unit.empty()) {
            throw py::value_error("chemical-observation parameter template is incomplete");
        }
        result.push_back(std::move(parameter));
    }
    return result;
}

struct ChemicalObservationSolve {
    ChemicalSolveResult result;
    ProviderActiveParameterSet active_parameters;
};

ChemicalObservationSolve solve_chemical_observation(
    ChemicalObservationContext& context,
    ChemicalObservationRow& row,
    const std::vector<double>& parameter_values
) {
    if (parameter_values.size() != context.parameter_templates.size()) {
        throw std::invalid_argument("chemical-observation parameter count is incorrect");
    }
    if (context.provider_sdk == nullptr) {
        throw std::runtime_error("Provider SDK lifetime is unavailable");
    }
    const epcsaft_native_sdk_v1& sdk = *context.provider_sdk;
    const auto metadata = chemical_provider_metadata(sdk);
    validate_provider_identity(metadata, row.input);
    const ProviderContext provider(sdk, context.provider_fingerprint);
    std::vector<ProviderActiveParameterInput> requests = context.parameter_templates;
    for (std::size_t index = 0; index < requests.size(); ++index) {
        requests[index].value = parameter_values[index];
    }
    ChemicalObservationSolve solved;
    solved.active_parameters = provider.resolve_active_parameters(
        row.temperature_k, requests
    );
    if (solved.active_parameters.parameters.size() != requests.size()) {
        throw std::invalid_argument("Provider active-parameter resolution is incomplete");
    }
    ReactionSystemInput input = row.input;
    std::vector<double> pressure_derivatives;
    std::vector<double> parameter_derivatives;
    parameter_derivatives.assign(
        input.reaction_matrix.rows * solved.active_parameters.parameters.size(),
        0.0
    );
    solved.result = solve_provider_reaction(
        compile_reaction_system(input),
        provider,
        row.temperature_k,
        row.pressure_pa,
        context.packing_bounds[0],
        context.packing_bounds[1],
        sdk.total_ion_mole_fraction_max,
        context.trace_floor,
        pressure_derivatives,
        parameter_derivatives,
        &solved.active_parameters
    );
    if (solved.result.sensitivities.status == "available") {
        solved.result.sensitivities.provider_parameter_status = "available";
        solved.result.sensitivities.provider_parameter_failure_reason.clear();
    } else {
        solved.result.sensitivities.provider_parameter_status = "unavailable";
        solved.result.sensitivities.provider_parameter_failure_reason =
            solved.result.sensitivities.failure_reason;
    }
    return solved;
}

void observation_text(char* target, const std::string& value) {
    constexpr std::size_t capacity = EPCSAFT_REGRESSION_EVALUATOR_V1_TEXT_CAPACITY;
    if (value.size() >= capacity) throw std::invalid_argument("diagnostic text is too long");
    std::fill(target, target + capacity, '\0');
    std::copy(value.begin(), value.end(), target);
}

bool evaluate_observation_row(
    ChemicalObservationContext& context,
    ChemicalObservationRow& row,
    const std::vector<double>& parameter_values,
    bool jacobian_requested,
    double& value,
    double* jacobian,
    epcsaft_regression_evaluator_row_result_v1& row_result
) {
    row_result = {};
    row_result.struct_size = sizeof(row_result);
    const auto solved = solve_chemical_observation(context, row, parameter_values);
    observation_text(row_result.solver_status, solved.result.solver_status);
    observation_text(row_result.numerical_status, solved.result.numerical_status);
    observation_text(row_result.physical_status, solved.result.physical_status);
    observation_text(row_result.derivative_status, solved.result.sensitivities.status);
    row_result.kkt_dimension = solved.result.sensitivities.kkt_dimension;
    row_result.kkt_rank = solved.result.sensitivities.kkt_rank;
    row_result.kkt_condition_number_inf = solved.result.sensitivities.condition_number_inf;
    if (!solved.result.accepted
        || (jacobian_requested && solved.result.sensitivities.status != "available")) {
        row_result.status = EPCSAFT_REGRESSION_EVALUATOR_STATUS_UNAVAILABLE_V1;
        observation_text(
            row_result.reason,
            solved.result.callback_error.empty()
                ? (solved.result.sensitivities.failure_reason.empty()
                    ? "homogeneous chemical state was not certified"
                    : solved.result.sensitivities.failure_reason)
                : solved.result.callback_error
        );
        return false;
    }
    const ProviderContext provider(*context.provider_sdk, context.provider_fingerprint);
    const auto phase = provider.evaluate_reacting_phase_parameters(
        row.temperature_k, solved.result.amounts, solved.result.volume_m3,
        context.packing_bounds[0], context.packing_bounds[1], solved.active_parameters
    );
    if (phase.topology_fingerprint
        != context.expected_provider_topology_fingerprint) {
        throw std::runtime_error("Provider active-parameter topology changed");
    }
    const std::size_t species_count = solved.result.amounts.size();
    const std::size_t state_count = species_count + 1;
    std::vector<std::size_t> columns;
    if (jacobian_requested) {
        for (const auto& parameter : solved.active_parameters.parameters) {
            const auto expected = provider_parameter_name(parameter);
            const auto it = std::find(
                solved.result.sensitivities.parameter_order.begin(),
                solved.result.sensitivities.parameter_order.end(), expected
            );
            if (it == solved.result.sensitivities.parameter_order.end()) {
                row_result.status = EPCSAFT_REGRESSION_EVALUATOR_STATUS_UNAVAILABLE_V1;
                observation_text(
                    row_result.reason,
                    "chemical sensitivity parameter order is incomplete"
                );
                return false;
            }
            columns.push_back(static_cast<std::size_t>(
                std::distance(solved.result.sensitivities.parameter_order.begin(), it)
            ));
        }
    }
    const std::size_t component = row.primitive.component_index;
    if (row.primitive.kind == "neutral_component_fugacity_pa") {
        const double mu = phase.phase.gradient[component];
        if (!std::isfinite(mu)
            || mu > std::log(std::numeric_limits<double>::max())
            || mu < std::log(std::numeric_limits<double>::min())) {
            throw std::domain_error("neutral component fugacity is outside binary64 range");
        }
        value = kUniversalGasConstant * row.temperature_k * std::exp(mu);
        for (std::size_t parameter = 0; parameter < columns.size(); ++parameter) {
            const auto source = columns[parameter];
            double dmu = phase.chemical_potential_parameter_derivatives_over_rt[
                component * columns.size() + parameter
            ];
            for (std::size_t species = 0; species < species_count; ++species) {
                dmu += phase.phase.hessian[component * state_count + species]
                    * solved.result.sensitivities.amount_derivatives[
                        source * species_count + species
                    ];
            }
            dmu += phase.phase.hessian[component * state_count + species_count]
                * solved.result.sensitivities.volume_derivatives[source];
            jacobian[parameter] = value * dmu;
        }
    } else {
        const double total = std::accumulate(
            solved.result.amounts.begin(), solved.result.amounts.end(), 0.0
        );
        if (!(total > 0.0) || !std::isfinite(total)) throw std::domain_error("chemical amount total is invalid");
        value = solved.result.amounts[component] / total;
        for (std::size_t parameter = 0; parameter < columns.size(); ++parameter) {
            const auto source = columns[parameter];
            double dtotal = 0.0;
            for (std::size_t species = 0; species < species_count; ++species) {
                dtotal += solved.result.sensitivities.amount_derivatives[
                    source * species_count + species
                ];
            }
            jacobian[parameter] = (
                solved.result.sensitivities.amount_derivatives[
                    source * species_count + component
                ] * total - solved.result.amounts[component] * dtotal
            ) / (total * total);
        }
    }
    if (!std::isfinite(value)
        || (jacobian_requested
            && !std::all_of(
                jacobian,
                jacobian + context.parameter_templates.size(),
                [](double derivative) { return std::isfinite(derivative); }
            ))) {
        throw std::domain_error("chemical-observation value or derivative is nonfinite");
    }
    observation_text(
        row_result.chart_topology,
        solved.result.sensitivities.chart_topology
    );
    observation_text(
        row_result.provider_topology_fingerprint,
        phase.topology_fingerprint
    );
    row_result.status = EPCSAFT_REGRESSION_EVALUATOR_STATUS_OK_V1;
    row_result.reason[0] = '\0';
    return true;
}

int evaluate_chemical_observation_v1(
    void* opaque,
    const double* parameter_values,
    std::size_t parameter_count,
    int32_t request_mode,
    epcsaft_regression_evaluator_result_v1* output
) {
    if (opaque == nullptr || output == nullptr) return EPCSAFT_REGRESSION_EVALUATOR_STATUS_INVALID_INPUT_V1;
    if (output->struct_size != sizeof(*output)) {
        return EPCSAFT_REGRESSION_EVALUATOR_STATUS_INVALID_INPUT_V1;
    }
    auto* context = static_cast<ChemicalObservationContext*>(opaque);
    const bool jacobian_requested = request_mode == EPCSAFT_REGRESSION_EVALUATOR_REQUEST_VALUES_AND_JACOBIAN_V1;
    try {
        if ((request_mode != EPCSAFT_REGRESSION_EVALUATOR_REQUEST_VALUES_ONLY_V1
                && request_mode
                    != EPCSAFT_REGRESSION_EVALUATOR_REQUEST_VALUES_AND_JACOBIAN_V1)
            || output->row_count != context->rows.size()
            || output->parameter_count != parameter_count
            || parameter_count != context->parameter_templates.size()
            || output->value_capacity < context->rows.size() || output->values == nullptr
            || output->row_result_capacity < context->rows.size() || output->row_results == nullptr
            || (parameter_count != 0 && parameter_values == nullptr)
            || (jacobian_requested && (output->jacobian == nullptr
                || output->jacobian_capacity < context->rows.size() * parameter_count))) {
            output->status = EPCSAFT_REGRESSION_EVALUATOR_STATUS_INVALID_INPUT_V1;
            observation_text(output->error, "invalid evaluator request or result buffers");
            return output->status;
        }
        for (std::size_t index = 0; index < context->rows.size(); ++index) {
            if (output->row_results[index].struct_size
                != sizeof(epcsaft_regression_evaluator_row_result_v1)) {
                output->status = EPCSAFT_REGRESSION_EVALUATOR_STATUS_INVALID_INPUT_V1;
                observation_text(output->error, "invalid evaluator row-result buffer");
                return output->status;
            }
        }
        const std::vector<double> values(parameter_values, parameter_values + parameter_count);
        std::fill(
            output->values,
            output->values + context->rows.size(),
            std::numeric_limits<double>::quiet_NaN()
        );
        if (jacobian_requested) {
            std::fill(
                output->jacobian,
                output->jacobian + context->rows.size() * parameter_count,
                std::numeric_limits<double>::quiet_NaN()
            );
        }
        bool all_ok = true;
        for (std::size_t index = 0; index < context->rows.size(); ++index) {
            try {
                all_ok = evaluate_observation_row(
                    *context, context->rows[index], values, jacobian_requested,
                    output->values[index],
                    jacobian_requested ? output->jacobian + index * parameter_count : nullptr,
                    output->row_results[index]
                ) && all_ok;
            } catch (const std::exception& error) {
                auto& result = output->row_results[index];
                result.status = EPCSAFT_REGRESSION_EVALUATOR_STATUS_UNAVAILABLE_V1;
                observation_text(result.reason, error.what());
                all_ok = false;
            } catch (...) {
                auto& result = output->row_results[index];
                result.status = EPCSAFT_REGRESSION_EVALUATOR_STATUS_UNAVAILABLE_V1;
                observation_text(result.reason, "unknown evaluator row failure");
                all_ok = false;
            }
        }
        observation_text(output->provider_parameter_fingerprint, context->provider_fingerprint);
        observation_text(output->artifact_identity, context->artifact_identity);
        output->request_mode = request_mode;
        output->status = all_ok
            ? EPCSAFT_REGRESSION_EVALUATOR_STATUS_OK_V1
            : EPCSAFT_REGRESSION_EVALUATOR_STATUS_UNAVAILABLE_V1;
        if (all_ok) {
            output->error[0] = '\0';
        } else {
            observation_text(output->error, "one or more evaluator rows failed");
        }
        return output->status;
    } catch (const std::exception& error) {
        output->status = EPCSAFT_REGRESSION_EVALUATOR_STATUS_UNAVAILABLE_V1;
        try { observation_text(output->error, error.what()); } catch (...) { output->error[0] = '\0'; }
        return output->status;
    } catch (...) {
        output->status = EPCSAFT_REGRESSION_EVALUATOR_STATUS_UNAVAILABLE_V1;
        try {
            observation_text(output->error, "unknown evaluator failure");
        } catch (...) {
            output->error[0] = '\0';
        }
        return output->status;
    }
}

ChemicalObservationRow observation_row(
    const py::dict& record,
    const ChemicalProviderMetadata& metadata
) {
    ChemicalObservationRow row;
    row.row_id = py::cast<std::string>(record["row_id"]);
    row.state_id = py::cast<std::string>(record["state_id"]);
    row.state_schema_id = py::cast<std::string>(record["state_schema_id"]);
    row.source_id = py::cast<std::string>(record["source_id"]);
    row.transform_id = py::cast<std::string>(record["transform_id"]);
    row.reference_id = py::cast<std::string>(record["reference_id"]);
    row.reference_fingerprint = py::cast<std::string>(
        record["reference_fingerprint"]
    );
    row.input = observation_input(py::cast<py::dict>(record["problem"]));
    row.temperature_k = row.input.temperature_k;
    row.pressure_pa = row.input.pressure_pa;
    const auto primitive = py::cast<py::dict>(record["primitive"]);
    row.primitive.kind = py::cast<std::string>(primitive["kind"]);
    row.primitive.component_id = py::cast<std::string>(primitive["component_id"]);
    if (row.row_id.empty() || row.state_id.empty() || row.state_schema_id.empty()
        || row.source_id.empty() || row.reference_id.empty()
        || row.reference_fingerprint.empty()) {
        throw py::value_error("chemical-observation row identity is incomplete");
    }
    if (row.transform_id != "identity" && row.transform_id != "natural_log") {
        throw py::value_error("unsupported chemical-observation transform");
    }
    if (row.primitive.kind != "neutral_component_fugacity_pa"
        && row.primitive.kind != "species_mole_fraction") {
        throw py::value_error("unsupported chemical-observation primitive");
    }
    const auto component = std::find(
        metadata.component_ids.begin(),
        metadata.component_ids.end(),
        row.primitive.component_id
    );
    if (component == metadata.component_ids.end()) {
        throw py::value_error("chemical-observation component is not in Provider order");
    }
    row.primitive.component_index = static_cast<std::size_t>(
        std::distance(metadata.component_ids.begin(), component)
    );
    if (row.primitive.kind == "neutral_component_fugacity_pa"
        && metadata.charges[row.primitive.component_index] != 0) {
        throw py::value_error(
            "neutral_component_fugacity_pa requires a neutral Provider component"
        );
    }
    validate_provider_identity(metadata, row.input);
    return row;
}

void populate_pointer_table(
    const std::vector<std::string>& values,
    std::vector<const char*>& pointers
) {
    pointers.reserve(values.size());
    for (const auto& value : values) {
        pointers.push_back(value.c_str());
    }
}

py::capsule make_chemical_observation_context(
    const py::capsule& provider_capsule,
    const py::tuple& row_records,
    const std::vector<double>& packing_bounds,
    double trace_floor,
    const py::tuple& parameter_templates,
    const std::string& evaluator_identity,
    const std::string& capability_id,
    const std::string& capability_fingerprint,
    const std::string& provider_artifact_identity,
    const std::string& owner_artifact_identity,
    const std::string& contract_fingerprint,
    const std::string& artifact_identity
) {
    const auto& provider_sdk = checked_chemical_sdk(provider_capsule);
    const auto metadata = chemical_provider_metadata(provider_sdk);
    if (row_records.empty()) {
        throw py::value_error("chemical-observation context requires rows");
    }
    if (packing_bounds.size() != 2
        || !std::isfinite(packing_bounds[0])
        || !std::isfinite(packing_bounds[1])
        || !(0.0 < packing_bounds[0] && packing_bounds[0] < packing_bounds[1])) {
        throw py::value_error("packing-fraction bounds are invalid");
    }
    if (!std::isfinite(trace_floor) || trace_floor <= 0.0) {
        throw py::value_error("trace floor must be finite and positive");
    }
    if (evaluator_identity.empty() || capability_id.empty()
        || capability_fingerprint.empty() || provider_artifact_identity.empty()
        || owner_artifact_identity.empty() || contract_fingerprint.empty()
        || artifact_identity.empty()) {
        throw py::value_error("installed evaluator artifact identity is incomplete");
    }
    auto context = std::make_unique<ChemicalObservationContext>();
    context->provider_sdk = &provider_sdk;
    context->provider_capsule = provider_capsule;
    context->packing_bounds = packing_bounds;
    context->trace_floor = trace_floor;
    context->parameter_templates = observation_parameter_templates(parameter_templates);
    if (context->parameter_templates.empty()) {
        throw py::value_error("chemical-observation context requires active parameters");
    }
    context->provider_fingerprint = py::cast<std::string>(
        py::cast<py::dict>(
            py::cast<py::dict>(row_records[0])["problem"]
        )["provider_fingerprint"]
    );
    context->evaluator_identity = evaluator_identity;
    context->capability_id = capability_id;
    context->capability_fingerprint = capability_fingerprint;
    context->provider_artifact_identity = provider_artifact_identity;
    context->owner_artifact_identity = owner_artifact_identity;
    context->contract_fingerprint = contract_fingerprint;
    context->artifact_identity = artifact_identity;

    context->rows.reserve(row_records.size());
    for (const auto item : row_records) {
        auto row = observation_row(py::cast<py::dict>(item), metadata);
        if (row.input.provider_fingerprint != context->provider_fingerprint) {
            throw py::value_error("chemical-observation rows have inconsistent fingerprints");
        }
        if (std::find(
                context->row_ids.begin(), context->row_ids.end(), row.row_id
            ) != context->row_ids.end()) {
            throw py::value_error("chemical-observation row identities must be unique");
        }
        context->row_ids.push_back(row.row_id);
        context->state_ids.push_back(row.state_id);
        context->state_schema_ids.push_back(row.state_schema_id);
        context->source_ids.push_back(row.source_id);
        context->primitive_ids.push_back(
            row.primitive.kind + ";" + row.primitive.component_id
        );
        context->primitive_units.push_back(
            row.primitive.kind == "neutral_component_fugacity_pa"
                ? "Pa"
                : "dimensionless"
        );
        context->transform_ids.push_back(row.transform_id);
        context->reference_ids.push_back(row.reference_id);
        context->reference_fingerprints.push_back(row.reference_fingerprint);
        context->rows.push_back(std::move(row));
    }
    for (const auto& parameter : context->parameter_templates) {
        std::string components;
        for (std::size_t index = 0; index < parameter.component_ids.size(); ++index) {
            if (index != 0) {
                components += ",";
            }
            components += parameter.component_ids[index];
        }
        context->parameter_ids.push_back(
            parameter.family + ";" + parameter.identity + ";" + components
        );
        context->parameter_units.push_back(parameter.unit);
    }
    const ProviderContext provider(provider_sdk, context->provider_fingerprint);
    for (const auto& row : context->rows) {
        const auto active = provider.resolve_active_parameters(
            row.temperature_k,
            context->parameter_templates
        );
        if (active.parameters.size() != context->parameter_templates.size()
            || active.topology_fingerprint.empty()) {
            throw py::value_error(
                "Provider active-parameter topology is incomplete"
            );
        }
        if (context->expected_provider_topology_fingerprint.empty()) {
            context->expected_provider_topology_fingerprint =
                active.topology_fingerprint;
        } else if (active.topology_fingerprint
            != context->expected_provider_topology_fingerprint) {
            throw py::value_error(
                "chemical-observation rows require one Provider topology"
            );
        }
    }
    populate_pointer_table(context->row_ids, context->row_id_ptrs);
    populate_pointer_table(context->state_ids, context->state_id_ptrs);
    populate_pointer_table(context->state_schema_ids, context->state_schema_id_ptrs);
    populate_pointer_table(context->source_ids, context->source_id_ptrs);
    populate_pointer_table(context->primitive_ids, context->primitive_id_ptrs);
    populate_pointer_table(context->primitive_units, context->primitive_unit_ptrs);
    populate_pointer_table(context->transform_ids, context->transform_id_ptrs);
    populate_pointer_table(context->reference_ids, context->reference_id_ptrs);
    populate_pointer_table(
        context->reference_fingerprints,
        context->reference_fingerprint_ptrs
    );
    populate_pointer_table(context->parameter_ids, context->parameter_id_ptrs);
    populate_pointer_table(context->parameter_units, context->parameter_unit_ptrs);

    auto& sdk = context->sdk;
    sdk.abi_version = EPCSAFT_REGRESSION_EVALUATOR_V1_ABI_VERSION;
    sdk.table_size = sizeof(sdk);
    sdk.model_context = context.get();
    sdk.row_count = context->rows.size();
    sdk.parameter_count = context->parameter_templates.size();
    sdk.row_ids = context->row_id_ptrs.data();
    sdk.state_ids = context->state_id_ptrs.data();
    sdk.state_schema_ids = context->state_schema_id_ptrs.data();
    sdk.observation_source_ids = context->source_id_ptrs.data();
    sdk.primitive_ids = context->primitive_id_ptrs.data();
    sdk.primitive_units = context->primitive_unit_ptrs.data();
    sdk.transform_ids = context->transform_id_ptrs.data();
    sdk.reference_ids = context->reference_id_ptrs.data();
    sdk.reference_fingerprints = context->reference_fingerprint_ptrs.data();
    sdk.parameter_ids = context->parameter_id_ptrs.data();
    sdk.parameter_units = context->parameter_unit_ptrs.data();
    sdk.evaluator_identity = context->evaluator_identity.c_str();
    sdk.capability_id = context->capability_id.c_str();
    sdk.capability_fingerprint = context->capability_fingerprint.c_str();
    sdk.provider_artifact_identity = context->provider_artifact_identity.c_str();
    sdk.owner_artifact_identity = context->owner_artifact_identity.c_str();
    sdk.contract_fingerprint = context->contract_fingerprint.c_str();
    sdk.model_fingerprint = context->provider_fingerprint.c_str();
    sdk.provider_parameter_fingerprint = context->provider_fingerprint.c_str();
    sdk.expected_provider_topology_fingerprint =
        context->expected_provider_topology_fingerprint.c_str();
    sdk.provider_sdk_capsule_name = EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME;
    sdk.provider_sdk_abi_version = provider_sdk.abi_version;
    sdk.provider_sdk_table_size = provider_sdk.table_size;
    sdk.provider_sdk_result_size = provider_sdk.result_size;
    sdk.provider_sdk_mixture_result_size = provider_sdk.mixture_result_size;
    sdk.provider_sdk_neutral_reference_result_size =
        provider_sdk.neutral_reference_result_size;
    sdk.provider_sdk_neutral_reference_derivative_result_size =
        provider_sdk.neutral_reference_derivative_result_size;
    sdk.provider_sdk_reacting_phase_parameter_result_size =
        provider_sdk.reacting_phase_parameter_result_size;
    sdk.artifact_identity = context->artifact_identity.c_str();
    sdk.single_thread_non_reentrant = 1;
    sdk.value_only_avoids_derivative_work = 0;
    sdk.result_size = sizeof(epcsaft_regression_evaluator_result_v1);
    sdk.row_result_size = sizeof(epcsaft_regression_evaluator_row_result_v1);
    sdk.evaluate = &evaluate_chemical_observation_v1;
    auto* sdk_pointer = &sdk;
    context.release();
    return py::capsule(
        sdk_pointer,
        kChemicalObservationCapsuleName,
        &destroy_chemical_observation_context
    );
}

std::vector<std::string> text_table(
    const char* const* values,
    std::size_t count
) {
    std::vector<std::string> result;
    result.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        result.emplace_back(values[index]);
    }
    return result;
}

py::dict evaluate_chemical_observation(
    const py::capsule& capsule,
    const std::vector<double>& parameter_values,
    bool with_jacobian
) {
    auto& sdk = checked_chemical_observation_sdk(capsule);
    std::vector<double> values(sdk.row_count);
    std::vector<double> jacobian(
        with_jacobian ? sdk.row_count * sdk.parameter_count : 0
    );
    std::vector<epcsaft_regression_evaluator_row_result_v1> rows(sdk.row_count);
    for (auto& row : rows) {
        row.struct_size = sizeof(row);
    }
    epcsaft_regression_evaluator_result_v1 result{};
    result.struct_size = sizeof(result);
    result.row_count = sdk.row_count;
    result.parameter_count = sdk.parameter_count;
    result.value_capacity = values.size();
    result.jacobian_capacity = jacobian.size();
    result.row_result_capacity = rows.size();
    result.values = values.data();
    result.jacobian = with_jacobian ? jacobian.data() : nullptr;
    result.row_results = rows.data();
    int status = EPCSAFT_REGRESSION_EVALUATOR_STATUS_UNAVAILABLE_V1;
    {
        py::gil_scoped_release release;
        status = sdk.evaluate(
            sdk.model_context,
            parameter_values.data(),
            parameter_values.size(),
            with_jacobian
                ? EPCSAFT_REGRESSION_EVALUATOR_REQUEST_VALUES_AND_JACOBIAN_V1
                : EPCSAFT_REGRESSION_EVALUATOR_REQUEST_VALUES_ONLY_V1,
            &result
        );
    }
    if (status != result.status) {
        throw std::runtime_error("evaluator callback returned inconsistent status");
    }
    py::list row_results;
    for (std::size_t index = 0; index < sdk.row_count; ++index) {
        const auto& row = rows[index];
        py::dict record;
        record["row_id"] = sdk.row_ids[index];
        record["state_id"] = sdk.state_ids[index];
        record["status"] = row.status;
        record["reason"] = row.reason;
        record["solver_status"] = row.solver_status;
        record["numerical_status"] = row.numerical_status;
        record["physical_status"] = row.physical_status;
        record["derivative_status"] = row.derivative_status;
        record["chart_topology"] = row.chart_topology;
        record["provider_topology_fingerprint"] =
            row.provider_topology_fingerprint;
        record["kkt_dimension"] = row.kkt_dimension;
        record["kkt_rank"] = row.kkt_rank;
        record["kkt_condition_number_inf"] = row.kkt_condition_number_inf;
        row_results.append(std::move(record));
    }
    py::dict output;
    output["status"] = result.status;
    output["error"] = result.error;
    output["request_mode"] = result.request_mode;
    output["values"] = std::move(values);
    output["jacobian"] = std::move(jacobian);
    output["row_results"] = std::move(row_results);
    output["row_ids"] = text_table(sdk.row_ids, sdk.row_count);
    output["state_ids"] = text_table(sdk.state_ids, sdk.row_count);
    output["state_schema_ids"] = text_table(sdk.state_schema_ids, sdk.row_count);
    output["source_ids"] = text_table(
        sdk.observation_source_ids, sdk.row_count
    );
    output["primitive_ids"] = text_table(sdk.primitive_ids, sdk.row_count);
    output["primitive_units"] = text_table(sdk.primitive_units, sdk.row_count);
    output["transform_ids"] = text_table(sdk.transform_ids, sdk.row_count);
    output["reference_ids"] = text_table(sdk.reference_ids, sdk.row_count);
    output["reference_fingerprints"] = text_table(
        sdk.reference_fingerprints, sdk.row_count
    );
    output["parameter_ids"] = text_table(sdk.parameter_ids, sdk.parameter_count);
    output["parameter_units"] = text_table(sdk.parameter_units, sdk.parameter_count);
    output["evaluator_identity"] = sdk.evaluator_identity;
    output["capability_id"] = sdk.capability_id;
    output["capability_fingerprint"] = sdk.capability_fingerprint;
    output["provider_artifact_identity"] = sdk.provider_artifact_identity;
    output["owner_artifact_identity"] = sdk.owner_artifact_identity;
    output["contract_fingerprint"] = sdk.contract_fingerprint;
    output["model_fingerprint"] = sdk.model_fingerprint;
    output["provider_parameter_fingerprint"] = result.provider_parameter_fingerprint;
    output["expected_provider_topology_fingerprint"] =
        sdk.expected_provider_topology_fingerprint;
    output["provider_sdk_capsule_name"] = sdk.provider_sdk_capsule_name;
    output["provider_sdk_abi_version"] = sdk.provider_sdk_abi_version;
    output["provider_sdk_table_size"] = sdk.provider_sdk_table_size;
    output["provider_sdk_result_size"] = sdk.provider_sdk_result_size;
    output["provider_sdk_mixture_result_size"] =
        sdk.provider_sdk_mixture_result_size;
    output["provider_sdk_neutral_reference_result_size"] =
        sdk.provider_sdk_neutral_reference_result_size;
    output["provider_sdk_neutral_reference_derivative_result_size"] =
        sdk.provider_sdk_neutral_reference_derivative_result_size;
    output["provider_sdk_reacting_phase_parameter_result_size"] =
        sdk.provider_sdk_reacting_phase_parameter_result_size;
    output["artifact_identity"] = result.artifact_identity;
    output["single_thread_non_reentrant"] = sdk.single_thread_non_reentrant != 0;
    output["value_only_avoids_derivative_work"] =
        sdk.value_only_avoids_derivative_work != 0;
    return output;
}

}  // namespace

void bind_chemical_observation(py::module_& module) {
    module.def(
        "_chemical_observation_context",
        &make_chemical_observation_context,
        py::arg("provider_capsule"),
        py::arg("rows"),
        py::arg("packing_fraction_bounds"),
        py::arg("trace_floor"),
        py::arg("parameter_templates"),
        py::arg("evaluator_identity"),
        py::arg("capability_id"),
        py::arg("capability_fingerprint"),
        py::arg("provider_artifact_identity"),
        py::arg("owner_artifact_identity"),
        py::arg("contract_fingerprint"),
        py::arg("artifact_identity")
    );
    module.def(
        "_chemical_observation_evaluate",
        &evaluate_chemical_observation,
        py::arg("context"),
        py::arg("parameter_values"),
        py::arg("with_jacobian") = true
    );
}

}  // namespace epcsaft_equilibrium
