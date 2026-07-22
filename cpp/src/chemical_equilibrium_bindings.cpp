#include "chemical_equilibrium.hpp"
#include "provider.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
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
            py::cast<bool>(record["dimensionless"]),
            py::cast<double>(record["temperature_k"]),
            py::cast<double>(record["pressure_pa"]),
            py::cast<std::string>(record["pressure_binding"]),
        });
    }
    return result;
}

ReactionSystemInput reaction_system_input(const py::dict& spec) {
    ReactionSystemInput input;
    input.species_ids = py::cast<std::vector<std::string>>(spec["species_ids"]);
    input.charges = py::cast<std::vector<int>>(spec["charges"]);
    input.provider_component_ids = py::cast<std::vector<std::string>>(
        spec["provider_component_ids"]
    );
    input.provider_charges = py::cast<std::vector<int>>(spec["provider_charges"]);
    input.provider_fingerprint = py::cast<std::string>(spec["provider_fingerprint"]);
    input.expected_provider_component_ids = py::cast<std::vector<std::string>>(
        spec["expected_provider_component_ids"]
    );
    input.expected_provider_charges = py::cast<std::vector<int>>(
        spec["expected_provider_charges"]
    );
    input.expected_provider_fingerprint = py::cast<std::string>(
        spec["expected_provider_fingerprint"]
    );
    input.balance_matrix = dense_matrix(spec["balance_matrix"], "balance matrix");
    input.declared_balance_rank = py::cast<std::size_t>(spec["declared_balance_rank"]);
    input.reaction_matrix = dense_matrix(spec["reaction_matrix"], "reaction matrix");
    input.feed_amounts = py::cast<std::vector<double>>(spec["feed_amounts"]);
    input.ln_k = py::cast<std::vector<double>>(spec["ln_k"]);
    input.equilibrium_constant_records = equilibrium_constant_records(
        spec["equilibrium_constant_records"]
    );
    input.temperature_k = py::cast<double>(spec["temperature_k"]);
    input.pressure_pa = py::cast<double>(spec["pressure_pa"]);
    input.complete_closed_system = py::cast<bool>(spec["complete_closed_system"]);
    return input;
}

py::dict compile_system(const py::dict& spec) {
    const CompiledReactionSystem compiled = compile_reaction_system(
        reaction_system_input(spec)
    );
    py::dict result;
    result["species_ids"] = compiled.species_ids;
    result["charges"] = compiled.charges;
    result["balance_rank"] = compiled.balance_rank;
    result["reaction_rank"] = compiled.reaction_rank;
    result["reaction_qr_diagonal_ratio"] = compiled.reaction_qr_diagonal_ratio;
    result["balance_totals"] = compiled.balance_totals;
    result["g_ref"] = compiled.g_ref;
    result["reference_reconstruction_inf_norm"] =
        compiled.reference_reconstruction_inf_norm;
    result["conservation_reaction_inf_norm"] =
        compiled.conservation_reaction_inf_norm;
    result["charge_reaction_inf_norm"] = compiled.charge_reaction_inf_norm;
    result["provider_fingerprint"] = compiled.provider_fingerprint;
    return result;
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

std::vector<double> amount_chart_inverse(
    const std::vector<int>& charges,
    const std::vector<double>& amounts
) {
    return invert_amount_chart(make_amount_chart(charges), amounts);
}

py::dict max_min_evidence(
    const py::handle& balance_matrix,
    const std::vector<double>& feed_amounts,
    const std::vector<int>& charges,
    double trace_floor
) {
    const MaxMinInitializationResult initialization = max_min_initialization(
        dense_matrix(balance_matrix, "balance matrix"),
        feed_amounts,
        charges,
        trace_floor,
        std::numeric_limits<double>::quiet_NaN()
    );
    py::dict result;
    result["solver_status"] = initialization.solver_status;
    result["reason"] = initialization.reason;
    result["amounts"] = initialization.amounts;
    result["max_min_amount"] = initialization.max_min_amount;
    result["equality_inf_norm"] = initialization.equality_inf_norm;
    result["strict_positive_feasible"] = initialization.strict_positive_feasible;
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
    result["artifact_input_status"] = evaluation.artifact_input_status;
    result["numerical_status"] = evaluation.numerical_status;
    result["physical_status"] = evaluation.physical_status;
    result["provider_domain_status"] = evaluation.provider_domain_status;
    result["local_minimum_status"] = evaluation.local_minimum_status;
    result["trace_status"] = evaluation.trace_status;
    result["predictive_status"] = evaluation.predictive_status;
    result["finite_search_status"] = evaluation.finite_search_status;
    result["globality_certificate"] = "not_guaranteed";
    result["amounts"] = evaluation.amounts;
    result["volume_m3"] = evaluation.volume_m3;
    result["objective"] = evaluation.objective;
    result["balance_inf_norm"] = evaluation.balance_inf_norm;
    result["charge_inf_norm"] = evaluation.charge_inf_norm;
    result["pressure_relative_residual"] = evaluation.pressure_relative_residual;
    result["reaction_affinity_inf_norm"] = evaluation.reaction_affinity_inf_norm;
    result["packing_fraction"] = evaluation.packing_fraction;
    result["kkt_stationarity_inf_norm"] = evaluation.kkt_stationarity_inf_norm;
    result["complementarity_inf_norm"] = evaluation.complementarity_inf_norm;
    result["kkt_scope"] = evaluation.kkt_scope;
    result["final_lambda"] = evaluation.has_final_lambda
        ? py::cast(evaluation.final_lambda)
        : py::none();
    result["continuation_used"] = evaluation.continuation_used;
    result["max_min_solve_count"] = evaluation.max_min_solve_count;
    result["active_balance_constraint_count"] =
        evaluation.active_balance_constraint_count;
    result["kkt_residual"] = evaluation.kkt_residual;
    result["kkt_jacobian"] = evaluation.kkt_jacobian;
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
    result["value"] = evaluation.value;
    result["gradient"] = evaluation.gradient;
    result["hessian"] = evaluation.hessian;
    result["pressure_pa"] = evaluation.pressure_pa;
    result["packing_fraction"] = evaluation.packing_fraction;
    result["packing_gradient"] = evaluation.packing_gradient;
    result["packing_hessian"] = evaluation.packing_hessian;
    result["parameter_fingerprint"] = evaluation.parameter_fingerprint;
    result["component_ids"] = metadata.component_ids;
    result["charges"] = metadata.charges;
    return result;
}

py::dict neutral_reference_record(
    const NeutralReferenceEvaluation& evaluation,
    const ChemicalProviderMetadata& metadata,
    double temperature_k,
    double pressure_pa
) {
    py::dict result;
    result["component_count"] = evaluation.component_count;
    result["neutral_basis_row_count"] = evaluation.neutral_basis_row_count;
    py::list basis;
    for (std::size_t row = 0; row < evaluation.neutral_basis_row_count; ++row) {
        std::vector<double> values(
            evaluation.neutral_basis.begin()
                + static_cast<std::ptrdiff_t>(row * evaluation.component_count),
            evaluation.neutral_basis.begin()
                + static_cast<std::ptrdiff_t>((row + 1) * evaluation.component_count)
        );
        basis.append(values);
    }
    result["neutral_basis"] = basis;
    result["log_fugacity_contractions"] = evaluation.log_fugacity_contractions;
    result["reference_composition"] = evaluation.reference_composition;
    result["derivative_availability"] = evaluation.derivative_availability;
    result["temperature_k"] = evaluation.temperature_k;
    result["pressure_pa"] = evaluation.pressure_pa;
    result["solvent_molar_mass_kg_per_mol"] =
        evaluation.solvent_molar_mass_kg_per_mol;
    result["reference_amount_mol"] = evaluation.reference_amount_mol;
    result["reference_number_density_mol_per_m3"] =
        evaluation.reference_number_density_mol_per_m3;
    result["reference_molality_mol_per_kg"] =
        evaluation.reference_molality_mol_per_kg;
    result["reference_convergence_error"] = evaluation.reference_convergence_error;
    result["basis_id"] = evaluation.basis_id;
    result["parameter_fingerprint"] = evaluation.parameter_fingerprint;
    result["component_ids"] = metadata.component_ids;
    result["charges"] = metadata.charges;
    if (evaluation.temperature_k != temperature_k || evaluation.pressure_pa != pressure_pa) {
        throw py::value_error("Provider neutral-reference T/P binding changed");
    }
    return result;
}

py::dict provider_neutral_reference_evidence(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::string& expected_fingerprint
) {
    const epcsaft_native_sdk_v1& sdk = checked_chemical_sdk(capsule);
    const ChemicalProviderMetadata metadata = chemical_provider_metadata(sdk);
    const ProviderContext provider(sdk, expected_fingerprint);
    return neutral_reference_record(
        provider.evaluate_neutral_reference(temperature_k, pressure_pa),
        metadata,
        temperature_k,
        pressure_pa
    );
}

py::dict source_standard_state_evidence(const py::dict& source) {
    const std::vector<std::string> component_ids = py::cast<std::vector<std::string>>(
        source["component_ids"]
    );
    const std::vector<int> charges = py::cast<std::vector<int>>(source["charges"]);
    const std::string fingerprint = py::cast<std::string>(source["provider_fingerprint"]);
    const double temperature_k = py::cast<double>(source["temperature_k"]);
    const double pressure_pa = py::cast<double>(source["pressure_pa"]);
    const DenseMatrix reaction_matrix = dense_matrix(source["reaction_matrix"], "reaction matrix");
    const std::vector<double> source_ln_k = py::cast<std::vector<double>>(source["ln_k"]);
    const std::vector<double> scales = py::cast<std::vector<double>>(
        source["log_activity_scale_factors"]
    );
    const py::dict reference_record = py::cast<py::dict>(source["reference"]);
    const std::vector<std::string> reference_component_ids =
        py::cast<std::vector<std::string>>(reference_record["component_ids"]);
    const std::vector<int> reference_charges = py::cast<std::vector<int>>(
        reference_record["charges"]
    );
    if (reference_component_ids != component_ids || reference_charges != charges) {
        throw py::value_error("Provider neutral-reference component identity mismatch");
    }
    NeutralReferenceEvaluation reference;
    reference.component_count = py::cast<std::size_t>(reference_record["component_count"]);
    reference.neutral_basis_row_count = py::cast<std::size_t>(
        reference_record["neutral_basis_row_count"]
    );
    const std::vector<std::vector<double>> basis = py::cast<std::vector<std::vector<double>>>(
        reference_record["neutral_basis"]
    );
    for (const std::vector<double>& row : basis) {
        reference.neutral_basis.insert(
            reference.neutral_basis.end(), row.begin(), row.end()
        );
    }
    reference.log_fugacity_contractions = py::cast<std::vector<double>>(
        reference_record["log_fugacity_contractions"]
    );
    reference.reference_composition = py::cast<std::vector<double>>(
        reference_record["reference_composition"]
    );
    reference.derivative_availability = py::cast<std::uint32_t>(
        reference_record["derivative_availability"]
    );
    reference.temperature_k = py::cast<double>(reference_record["temperature_k"]);
    reference.pressure_pa = py::cast<double>(reference_record["pressure_pa"]);
    reference.solvent_molar_mass_kg_per_mol = py::cast<double>(
        reference_record["solvent_molar_mass_kg_per_mol"]
    );
    reference.reference_amount_mol = py::cast<double>(reference_record["reference_amount_mol"]);
    reference.reference_number_density_mol_per_m3 = py::cast<double>(
        reference_record["reference_number_density_mol_per_m3"]
    );
    reference.reference_molality_mol_per_kg = py::cast<double>(
        reference_record["reference_molality_mol_per_kg"]
    );
    reference.reference_convergence_error = py::cast<double>(
        reference_record["reference_convergence_error"]
    );
    reference.basis_id = py::cast<std::string>(reference_record["basis_id"]);
    reference.parameter_fingerprint = py::cast<std::string>(
        reference_record["parameter_fingerprint"]
    );
    if (!source.contains("equilibrium_constant_records")) {
        throw py::value_error("source equilibrium-constant records are incomplete");
    }
    const py::tuple records = py::cast<py::tuple>(source["equilibrium_constant_records"]);
    if (records.size() != reaction_matrix.rows) {
        throw py::value_error("source equilibrium-constant records do not match reactions");
    }
    std::string source_standard_state_id;
    for (const py::handle item : records) {
        const py::dict record = py::cast<py::dict>(item);
        const std::string source_id = py::cast<std::string>(record["source_id"]);
        const std::string standard_state_id = py::cast<std::string>(
            record["source_standard_state_id"]
        );
        if (source_id.empty() || standard_state_id.empty()) {
            throw py::value_error("source identity is incomplete");
        }
        if (source_standard_state_id.empty()) {
            source_standard_state_id = standard_state_id;
        } else if (source_standard_state_id != standard_state_id) {
            throw py::value_error("source standard-state identity is inconsistent");
        }
        if (!py::cast<bool>(record["dimensionless"])
            || py::cast<std::string>(record["reference_id"]) != reference.basis_id
            || py::cast<double>(record["temperature_k"]) != temperature_k
            || py::cast<double>(record["pressure_pa"]) != pressure_pa
            || py::cast<std::string>(record["pressure_binding"]) != "fixed") {
            throw py::value_error("source/reference standard-state identity is incompatible");
        }
    }
    const SourceStandardStateResult transformed = transform_source_standard_state(
        reaction_matrix,
        source_ln_k,
        scales,
        charges,
        component_ids,
        fingerprint,
        temperature_k,
        pressure_pa,
        reference
    );
    py::dict result;
    result["standard_offsets"] = transformed.standard_offsets;
    result["ln_k_provider_basis"] = transformed.ln_k_provider_basis;
    py::list coordinates;
    for (std::size_t row = 0; row < transformed.reaction_to_neutral_basis.rows; ++row) {
        std::vector<double> values(
            transformed.reaction_to_neutral_basis.values.begin()
                + static_cast<std::ptrdiff_t>(
                    row * transformed.reaction_to_neutral_basis.columns
                ),
            transformed.reaction_to_neutral_basis.values.begin()
                + static_cast<std::ptrdiff_t>(
                    (row + 1) * transformed.reaction_to_neutral_basis.columns
                )
        );
        coordinates.append(values);
    }
    result["reaction_to_neutral_basis"] = coordinates;
    result["representation_residual_inf_norm"] =
        transformed.representation_residual_inf_norm;
    result["basis_condition_ratio"] = transformed.basis_condition_ratio;
    result["derivative_availability"] = transformed.derivative_availability;
    result["basis_id"] = transformed.basis_id;
    result["parameter_fingerprint"] = transformed.parameter_fingerprint;
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
    if (input.provider_component_ids != metadata.component_ids) {
        throw py::value_error("Provider capsule component order does not match the reaction system");
    }
    if (input.provider_charges != metadata.charges) {
        throw py::value_error("Provider capsule charges do not match the reaction system");
    }
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
    const ProviderContext provider(sdk, input.provider_fingerprint);
    py::dict result = chemical_result(
        "installed_provider_manufactured_nonpredictive",
        solve_provider_reaction(
            compiled,
            provider,
            input.temperature_k,
            input.pressure_pa,
            packing_bounds[0],
            packing_bounds[1],
            sdk.total_ion_mole_fraction_max,
            trace_floor,
            std::numeric_limits<double>::quiet_NaN(),
            {},
            max_iterations
        )
    );
    result["parameter_fingerprint"] = input.provider_fingerprint;
    result["packing_fraction_bounds"] = packing_bounds;
    result["predictive_status"] = "manufactured_nonpredictive";
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
    module.def("_chemical_compile_system", &compile_system, py::arg("spec"));
    module.def(
        "_chemical_amount_chart",
        &amount_chart_evidence,
        py::arg("charges"),
        py::arg("coordinates"),
        py::arg("trace_floor")
    );
    module.def(
        "_chemical_amount_chart_inverse",
        &amount_chart_inverse,
        py::arg("charges"),
        py::arg("amounts")
    );
    module.def(
        "_chemical_max_min_initialization",
        &max_min_evidence,
        py::arg("balance_matrix"),
        py::arg("feed_amounts"),
        py::arg("charges"),
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
        "_chemical_evaluate_provider_neutral_reference",
        &provider_neutral_reference_evidence,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "_chemical_transform_source_standard_state",
        &source_standard_state_evidence,
        py::arg("source")
    );
    module.def(
        "_chemical_solve_provider_manufactured",
        &solve_provider_manufactured,
        py::arg("capsule"),
        py::arg("spec"),
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
