#include "chemical_equilibrium.hpp"
#include "provider.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <epcsaft/native_sdk_v1.h>

namespace py = pybind11;

namespace epcsaft_equilibrium {
namespace {

constexpr std::size_t kChemicalSourceDomainSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, total_ion_mole_fraction_max) + sizeof(double);

struct ChemicalProviderMetadata {
    std::vector<std::string> component_ids;
    std::vector<int> charges;
};

const epcsaft_native_sdk_v1& checked_chemical_sdk(const py::capsule& capsule) {
    const char* name = capsule.name();
    if (name == nullptr || std::string_view(name) != EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME) {
        throw py::value_error("expected capsule epcsaft.native_sdk.v1");
    }
    const auto* sdk = capsule.get_pointer<epcsaft_native_sdk_v1>();
    if (sdk == nullptr || sdk->abi_version != EPCSAFT_NATIVE_SDK_V1_ABI_VERSION
        || sdk->table_size < kChemicalSourceDomainSdkTableSize) {
        throw py::value_error("Provider capsule is missing the reacting-phase SDK contract");
    }
    if (sdk->model_context == nullptr || sdk->component_count < 3
        || sdk->component_ids == nullptr || sdk->component_charges == nullptr
        || sdk->mixture_result_size != sizeof(epcsaft_mixture_phase_block_result_v1)
        || sdk->evaluate_electrolyte_phase == nullptr
        || sdk->evaluate_molar_volume_bounds == nullptr
        || sdk->evaluate_packing_fraction == nullptr) {
        throw py::value_error("Provider capsule is missing the reacting-phase callbacks");
    }
    if (!std::isfinite(sdk->source_temperature_min_k)
        || !std::isfinite(sdk->source_temperature_max_k)
        || sdk->source_temperature_min_k <= 0.0
        || sdk->source_temperature_max_k < sdk->source_temperature_min_k
        || (!std::isnan(sdk->total_ion_mole_fraction_max)
            && (!std::isfinite(sdk->total_ion_mole_fraction_max)
                || sdk->total_ion_mole_fraction_max < 0.0
                || sdk->total_ion_mole_fraction_max > 1.0))) {
        throw py::value_error("Provider source-domain metadata is invalid");
    }
    return *sdk;
}

ChemicalProviderMetadata chemical_provider_metadata(const epcsaft_native_sdk_v1& sdk) {
    ChemicalProviderMetadata result;
    result.component_ids.reserve(sdk.component_count);
    result.charges.reserve(sdk.component_count);
    for (std::size_t component = 0; component < sdk.component_count; ++component) {
        if (sdk.component_ids[component] == nullptr
            || sdk.component_ids[component][0] == '\0') {
            throw py::value_error("Provider component identity is incomplete");
        }
        result.component_ids.emplace_back(sdk.component_ids[component]);
        result.charges.push_back(static_cast<int>(sdk.component_charges[component]));
    }
    return result;
}

DenseMatrix dense_matrix(const py::handle& value, const char* field) {
    const std::vector<std::vector<double>> rows = py::cast<std::vector<std::vector<double>>>(
        value
    );
    if (rows.empty()) {
        return {};
    }
    const std::size_t columns = rows.front().size();
    if (columns == 0) {
        throw py::value_error(std::string(field) + " rows must not be empty");
    }
    DenseMatrix result{rows.size(), columns, {}};
    result.values.reserve(rows.size() * columns);
    for (const std::vector<double>& row : rows) {
        if (row.size() != columns) {
            throw py::value_error(std::string(field) + " must be rectangular");
        }
        result.values.insert(result.values.end(), row.begin(), row.end());
    }
    return result;
}

std::vector<EquilibriumConstantRecord> equilibrium_constant_records(
    const py::handle& value
) {
    std::vector<EquilibriumConstantRecord> result;
    for (const py::handle item : py::cast<py::tuple>(value)) {
        const py::dict record = py::cast<py::dict>(item);
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

ReactionSystemInput reaction_system_input(const py::dict& spec) {
    ReactionSystemInput input;
    input.species_ids = py::cast<std::vector<std::string>>(spec["species_ids"]);
    input.charges = py::cast<std::vector<int>>(spec["charges"]);
    input.provider_fingerprint = py::cast<std::string>(spec["provider_fingerprint"]);
    input.molar_masses_kg_per_mol = py::cast<std::vector<double>>(
        spec["molar_masses_kg_per_mol"]
    );
    input.balance_matrix = dense_matrix(spec["balance_matrix"], "balance matrix");
    input.reaction_matrix = dense_matrix(spec["reaction_matrix"], "reaction matrix");
    input.feed_amounts = py::cast<std::vector<double>>(spec["feed_amounts"]);
    input.ln_k = py::cast<std::vector<double>>(spec["ln_k"]);
    input.equilibrium_constant_records = equilibrium_constant_records(
        spec["equilibrium_constant_records"]
    );
    input.temperature_k = py::cast<double>(spec["temperature_k"]);
    input.pressure_pa = py::cast<double>(spec["pressure_pa"]);
    return input;
}

py::dict amount_chart_evidence(
    const std::vector<int>& charges,
    const std::vector<double>& coordinates,
    double trace_floor
) {
    if (!std::isfinite(trace_floor) || trace_floor <= 0.0) {
        throw py::value_error("trace floor must be finite and positive");
    }
    const AmountChart chart = make_amount_chart(charges);
    const AmountChartEvaluation evaluation = evaluate_amount_chart(chart, coordinates);
    py::dict result;
    result["amounts"] = evaluation.amounts;
    result["coordinate_count"] = chart.coordinate_count();
    result["jacobian"] = evaluation.jacobian;
    result["amount_hessians"] = evaluation.amount_hessians;
    result["minimum_amount"] = evaluation.minimum_amount;
    result["charge_residual"] = evaluation.charge_residual;
    result["trace_status"] = evaluation.minimum_amount <= trace_floor
        ? "at_or_below_floor"
        : "interior";
    return result;
}

py::dict chemical_result(
    const char* profile,
    const ChemicalSolveResult& evaluation
) {
    py::dict result;
    result["profile"] = profile;
    result["accepted"] = evaluation.accepted;
    result["solver_status"] = evaluation.solver_status;
    result["callback_error"] = evaluation.callback_error;
    result["chemical_certification_level"] =
        evaluation.chemical_certification_level;
    result["boundary_status"] = evaluation.boundary_status;
    result["structural_zero_species_indices"] =
        evaluation.structural_zero_species_indices;
    result["numerical_status"] = evaluation.numerical_status;
    result["physical_status"] = evaluation.physical_status;
    result["provider_domain_status"] = evaluation.provider_domain_status;
    result["local_minimum_status"] = evaluation.local_minimum_status;
    result["trace_status"] = evaluation.trace_status;
    result["globality_certificate"] = "not_guaranteed";
    result["amounts"] = evaluation.amounts;
    result["volume_m3"] = evaluation.volume_m3;
    result["balance_inf_norm"] = evaluation.balance_inf_norm;
    result["charge_inf_norm"] = evaluation.charge_inf_norm;
    result["pressure_relative_residual"] = evaluation.pressure_relative_residual;
    result["reaction_affinity_inf_norm"] = evaluation.reaction_affinity_inf_norm;
    result["packing_fraction"] = evaluation.packing_fraction;
    result["kkt_stationarity_inf_norm"] = evaluation.kkt_stationarity_inf_norm;
    result["complementarity_inf_norm"] = evaluation.complementarity_inf_norm;
    return result;
}

py::dict solve_manufactured(const py::dict& spec, const py::dict& options) {
    const ReactionSystemInput input = reaction_system_input(spec);
    const CompiledReactionSystem compiled = compile_reaction_system(input);
    const std::vector<double> gauge_coefficients = options.contains("gauge_coefficients")
        ? py::cast<std::vector<double>>(options["gauge_coefficients"])
        : std::vector<double>{};
    const double trace_floor = options.contains("trace_floor")
        ? py::cast<double>(options["trace_floor"])
        : 1.0e-12;
    const int max_iterations = options.contains("test_max_iterations")
        ? py::cast<int>(options["test_max_iterations"])
        : 500;
    return chemical_result(
        "manufactured_ideal_nonpredictive",
        solve_manufactured_ideal_reaction(
            compiled,
            input.temperature_k,
            input.pressure_pa,
            gauge_coefficients,
            trace_floor,
            max_iterations
        )
    );
}

py::dict provider_block_evidence(
    const py::capsule& capsule,
    double temperature_k,
    const std::vector<double>& amounts,
    double volume_m3,
    const std::string& expected_fingerprint
) {
    const epcsaft_native_sdk_v1& sdk = checked_chemical_sdk(capsule);
    const ChemicalProviderMetadata metadata = chemical_provider_metadata(sdk);
    if (amounts.size() != metadata.component_ids.size()) {
        throw py::value_error("Provider component order has the wrong length");
    }
    if (temperature_k < sdk.source_temperature_min_k
        || temperature_k > sdk.source_temperature_max_k) {
        throw py::value_error("temperature is outside the Provider source domain");
    }
    const ProviderContext provider(sdk, expected_fingerprint);
    const ProviderPhaseBlockEvidence evaluation = evaluate_provider_phase_block(
        provider, temperature_k, amounts, volume_m3
    );
    py::dict result;
    result["gradient"] = evaluation.gradient;
    result["pressure_pa"] = evaluation.pressure_pa;
    return result;
}

void validate_provider_identity(
    const ChemicalProviderMetadata& metadata,
    const ReactionSystemInput& input
) {
    if (input.species_ids != metadata.component_ids) {
        throw py::value_error("Provider capsule component order does not match the reaction system");
    }
    if (input.charges != metadata.charges) {
        throw py::value_error("Provider capsule charges do not match the reaction system");
    }
}

py::dict solve_provider_input(
    const epcsaft_native_sdk_v1& sdk,
    const ChemicalProviderMetadata& metadata,
    const ReactionSystemInput& input,
    const ProviderContext& provider,
    const char* profile,
    const py::dict& options
) {
    validate_provider_identity(metadata, input);
    if (input.temperature_k < sdk.source_temperature_min_k
        || input.temperature_k > sdk.source_temperature_max_k) {
        throw py::value_error("temperature is outside the Provider source domain");
    }
    const CompiledReactionSystem compiled = compile_reaction_system(input);
    double ionic_feed = 0.0;
    double total_feed = 0.0;
    for (std::size_t species = 0; species < input.feed_amounts.size(); ++species) {
        total_feed += input.feed_amounts[species];
        if (input.charges[species] != 0) {
            ionic_feed += input.feed_amounts[species];
        }
    }
    if (std::isfinite(sdk.total_ion_mole_fraction_max)
        && (total_feed <= 0.0
            || ionic_feed / total_feed > sdk.total_ion_mole_fraction_max)) {
        throw py::value_error("feed composition exceeds the Provider source domain");
    }
    if (!options.contains("packing_fraction_bounds")) {
        throw py::value_error("packing_fraction_bounds are required from the calling formulation");
    }
    const std::vector<double> packing_bounds = py::cast<std::vector<double>>(
        options["packing_fraction_bounds"]
    );
    if (packing_bounds.size() != 2) {
        throw py::value_error("packing_fraction_bounds must contain two values");
    }
    const double trace_floor = options.contains("trace_floor")
        ? py::cast<double>(options["trace_floor"])
        : 1.0e-12;
    const int max_iterations = options.contains("test_max_iterations")
        ? py::cast<int>(options["test_max_iterations"])
        : 500;
    py::dict result = chemical_result(
        profile,
        solve_provider_reaction(
            compiled,
            provider,
            input.temperature_k,
            input.pressure_pa,
            packing_bounds[0],
            packing_bounds[1],
            sdk.total_ion_mole_fraction_max,
            trace_floor,
            max_iterations
        )
    );
    result["parameter_fingerprint"] = input.provider_fingerprint;
    result["packing_fraction_bounds"] = packing_bounds;
    return result;
}

py::dict solve_provider_manufactured(
    const py::capsule& capsule,
    const py::dict& spec,
    const py::dict& options
) {
    const epcsaft_native_sdk_v1& sdk = checked_chemical_sdk(capsule);
    const ChemicalProviderMetadata metadata = chemical_provider_metadata(sdk);
    const ReactionSystemInput input = reaction_system_input(spec);
    const ProviderContext provider(sdk, input.provider_fingerprint);
    return solve_provider_input(
        sdk,
        metadata,
        input,
        provider,
        "installed_provider_manufactured_nonpredictive",
        options
    );
}

py::dict solve_provider_source(
    const py::capsule& capsule,
    const py::dict& spec,
    const py::dict& source_standard_state,
    const py::dict& options
) {
    const epcsaft_native_sdk_v1& sdk = checked_chemical_sdk(capsule);
    const ChemicalProviderMetadata metadata = chemical_provider_metadata(sdk);
    ReactionSystemInput input = reaction_system_input(spec);
    validate_provider_identity(metadata, input);
    const std::string source_standard_state_id = py::cast<std::string>(
        source_standard_state["id"]
    );
    const std::string source_activity_scale_id = py::cast<std::string>(
        source_standard_state["activity_scale_id"]
    );
    const std::vector<double> log_activity_scale_factors =
        py::cast<std::vector<double>>(
            source_standard_state["log_activity_scale_factors"]
        );
    if (source_standard_state_id.empty() || source_activity_scale_id.empty()) {
        throw py::value_error("source standard-state identity is incomplete");
    }
    for (const EquilibriumConstantRecord& record : input.equilibrium_constant_records) {
        if (record.conversion_id
                != "source-standard-state-to-provider-neutral-reference"
            || record.reaction_orientation != "products_positive"
            || !record.dimensionless
            || record.reference_id.empty()
            || record.temperature_k != input.temperature_k
            || record.pressure_pa != input.pressure_pa) {
            throw py::value_error(
                "source equilibrium-constant provenance is incompatible"
            );
        }
        if (source_standard_state_id != record.reference_id) {
            throw py::value_error(
                "source standard-state identity is inconsistent"
            );
        }
    }
    const ProviderContext provider(sdk, input.provider_fingerprint);
    const NeutralReferenceEvaluation reference = provider.evaluate_neutral_reference(
        input.temperature_k, input.pressure_pa
    );
    const SourceStandardStateResult transformed = transform_source_standard_state(
        input.reaction_matrix,
        input.ln_k,
        log_activity_scale_factors,
        input.charges,
        input.provider_fingerprint,
        input.temperature_k,
        input.pressure_pa,
        reference
    );
    input.ln_k = transformed.ln_k_provider_basis;
    for (EquilibriumConstantRecord& record : input.equilibrium_constant_records) {
        record.reference_id = "provider-helmholtz-coordinate-basis";
        record.conversion_id = "already-provider-basis";
    }
    py::dict result = solve_provider_input(
        sdk,
        metadata,
        input,
        provider,
        "installed_provider_source_bound",
        options
    );
    result["standard_offsets"] = transformed.standard_offsets;
    result["ln_k_provider_basis"] = transformed.ln_k_provider_basis;
    result["reference_representation_residual_inf_norm"] =
        transformed.representation_residual_inf_norm;
    result["source_standard_state_id"] = source_standard_state_id;
    result["source_activity_scale_id"] = source_activity_scale_id;
    result["provider_reference_id"] = reference.basis_id;
    result["reference_derivative_availability"] =
        reference.derivative_availability;
    result["reference_convergence_error"] =
        reference.reference_convergence_error;
    return result;
}

py::dict manufactured_nlp_evidence(
    const py::dict& spec,
    const std::vector<double>& variables,
    const std::vector<double>& constraint_multipliers,
    const std::vector<double>& gauge_coefficients
) {
    const ReactionSystemInput input = reaction_system_input(spec);
    const ManufacturedNlpEvaluation evaluation = evaluate_manufactured_reaction_nlp(
        compile_reaction_system(input),
        input.temperature_k,
        input.pressure_pa,
        gauge_coefficients,
        variables,
        constraint_multipliers
    );
    py::dict result;
    result["objective"] = evaluation.objective;
    result["objective_gradient"] = evaluation.objective_gradient;
    result["constraints"] = evaluation.constraints;
    result["constraint_jacobian"] = evaluation.constraint_jacobian;
    result["lagrangian_gradient"] = evaluation.lagrangian_gradient;
    result["lagrangian_hessian"] = evaluation.lagrangian_hessian;
    result["amounts"] = evaluation.amounts;
    result["volume_m3"] = evaluation.volume_m3;
    return result;
}

}  // namespace

void bind_chemical_equilibrium(py::module_& module) {
    module.def(
        "_chemical_amount_chart",
        &amount_chart_evidence,
        py::arg("charges"),
        py::arg("coordinates"),
        py::arg("trace_floor")
    );
    module.def(
        "_chemical_solve_manufactured",
        &solve_manufactured,
        py::arg("spec"),
        py::arg("options")
    );
    module.def(
        "_chemical_evaluate_provider_block",
        &provider_block_evidence,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("amounts"),
        py::arg("volume_m3"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "_chemical_solve_provider_manufactured",
        &solve_provider_manufactured,
        py::arg("capsule"),
        py::arg("spec"),
        py::arg("options")
    );
    module.def(
        "_chemical_solve_provider_source",
        &solve_provider_source,
        py::arg("capsule"),
        py::arg("spec"),
        py::arg("source_standard_state"),
        py::arg("options")
    );
    module.def(
        "_chemical_evaluate_manufactured_nlp",
        &manufactured_nlp_evidence,
        py::arg("spec"),
        py::arg("variables"),
        py::arg("constraint_multipliers"),
        py::arg("gauge_coefficients") = std::vector<double>{}
    );
}

}  // namespace epcsaft_equilibrium
