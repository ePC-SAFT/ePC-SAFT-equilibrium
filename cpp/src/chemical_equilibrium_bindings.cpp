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

constexpr std::string_view kWaterReferenceBasis =
    "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3";
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
    input.has_standard_state_transformation =
        spec.contains("has_standard_state_transformation")
        && py::cast<bool>(spec["has_standard_state_transformation"]);
    if (input.has_standard_state_transformation) {
        input.provider_basis_id = py::cast<std::string>(spec["provider_basis_id"]);
        input.standard_state_transformation_residual = py::cast<double>(
            spec["standard_state_transformation_residual"]
        );
    }
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
    result["balance_totals"] = compiled.balance_totals;
    result["g_ref"] = compiled.g_ref;
    result["reference_reconstruction_inf_norm"] =
        compiled.reference_reconstruction_inf_norm;
    result["conservation_reaction_inf_norm"] =
        compiled.conservation_reaction_inf_norm;
    result["charge_reaction_inf_norm"] = compiled.charge_reaction_inf_norm;
    result["provider_fingerprint"] = compiled.provider_fingerprint;
    result["has_standard_state_transformation"] =
        compiled.has_standard_state_transformation;
    result["provider_basis_id"] = compiled.provider_basis_id;
    result["standard_state_transformation_residual"] =
        compiled.standard_state_transformation_residual;
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

py::dict provider_standard_reference_evidence(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::string& expected_fingerprint
) {
    const epcsaft_native_sdk_v1& sdk = checked_chemical_sdk(capsule);
    const ChemicalProviderMetadata metadata = chemical_provider_metadata(sdk);
    if (metadata.component_ids
            != std::vector<std::string>{
                "water", "hydronium-cation", "hydroxide-anion"
            }
        || metadata.charges != std::vector<int>{0, 1, -1}) {
        throw py::value_error(
            "Provider standard-reference identity component order or charges mismatch"
        );
    }
    const ProviderContext provider(sdk, expected_fingerprint);
    const StandardReferenceEvaluation evaluation = provider.evaluate_standard_reference(
        temperature_k, pressure_pa
    );
    py::dict result;
    result["formula_unit_log_fugacity"] = evaluation.formula_unit_log_fugacity;
    result["pure_solvent_log_fugacity"] = evaluation.pure_solvent_log_fugacity;
    result["solvent_molar_mass_kg_per_mol"] =
        evaluation.solvent_molar_mass_kg_per_mol;
    result["reference_molality_mol_per_kg"] =
        evaluation.reference_molality_mol_per_kg;
    result["convergence_error"] = evaluation.convergence_error;
    result["basis_id"] = evaluation.basis_id;
    result["parameter_fingerprint"] = evaluation.parameter_fingerprint;
    result["component_ids"] = metadata.component_ids;
    result["charges"] = metadata.charges;
    result["temperature_k"] = temperature_k;
    result["pressure_pa"] = pressure_pa;
    return result;
}

py::dict transform_water_self_ionization_standard_state_evidence(
    const py::dict& record,
    const py::dict& reference_record
) {
    if (py::cast<int>(record["schema_version"]) != 1
        || py::cast<std::string>(record["capability"])
            != "private-water-self-ionization-reference-v1") {
        throw py::value_error("IAPWS source identity is incomplete");
    }
    const py::dict source = py::cast<py::dict>(record["source"]);
    if (py::cast<std::string>(source["id"]) != "iapws-r11-24"
        || py::cast<std::string>(source["url"])
            != "https://iapws.org/documents/release/Ionization.download"
        || py::cast<std::string>(source["equation_locator"])
            != "equations-1-through-4-and-table-1"
        || py::cast<std::string>(source["reaction"])
            != "2 H2O = H3O+ + OH-") {
        throw py::value_error("IAPWS source identity is incomplete");
    }

    const py::dict state = py::cast<py::dict>(record["state"]);
    const double temperature_k = py::cast<double>(state["temperature_k"]);
    const double pressure_pa = py::cast<double>(state["pressure_pa"]);
    if (temperature_k != 298.15) {
        throw py::value_error("IAPWS temperature binding does not match");
    }
    if (pressure_pa != 100000.0
        || py::cast<std::string>(state["pressure_binding"]) != "fixed") {
        throw py::value_error("IAPWS pressure binding does not match");
    }
    if (py::cast<std::string>(state["phase"]) != "ordinary-liquid-water") {
        throw py::value_error("IAPWS state identity does not match");
    }
    if (py::cast<double>(state["density_kg_per_m3"]) != 997.047039) {
        throw py::value_error("independent density record does not match");
    }
    const py::dict density_source = py::cast<py::dict>(state["density_source"]);
    if (py::cast<std::string>(density_source["id"]) != "llnl-tr-805304"
        || py::cast<std::string>(density_source["url"])
            != "https://www.osti.gov/servlets/purl/1630404"
        || py::cast<std::string>(density_source["locator"])
            != "section-4.2.3-page-35-high-density-solution") {
        throw py::value_error("independent density source identity is incomplete");
    }

    const py::dict standard_state = py::cast<py::dict>(record["standard_state"]);
    if (!py::cast<bool>(standard_state["equilibrium_constant_dimensionless"])) {
        throw py::value_error("IAPWS lnKw must be explicitly dimensionless");
    }
    if (py::cast<std::string>(standard_state["ionic_basis"]) != "molality"
        || py::cast<std::string>(standard_state["solvent_basis"]) != "mole-fraction") {
        throw py::value_error("IAPWS mixed standard-state identity does not match");
    }
    if (py::cast<double>(standard_state["standard_molality_mol_per_kg"]) != 1.0) {
        throw py::value_error("IAPWS thermodynamic standard molality must be 1 mol/kg");
    }
    if (py::cast<std::string>(standard_state["solvent_molar_mass_unit"])
        != "kg/mol") {
        throw py::value_error("standard-state molar-mass unit must be kg/mol");
    }
    if (py::cast<std::string>(standard_state["transformation_sign"])
        != "lnK_provider=lnKw_mixed+delta_standard_offset") {
        throw py::value_error("standard-state transformation sign does not match");
    }

    const py::dict provider_binding = py::cast<py::dict>(record["provider_binding"]);
    const std::string provider_basis_id = py::cast<std::string>(
        provider_binding["basis_id"]
    );
    if (provider_basis_id != kWaterReferenceBasis) {
        throw py::value_error("Provider basis identity does not match");
    }
    const std::vector<std::string> component_ids = py::cast<std::vector<std::string>>(
        provider_binding["component_ids"]
    );
    const std::vector<int> charges = py::cast<std::vector<int>>(
        provider_binding["charges"]
    );
    if (component_ids
        != std::vector<std::string>{"water", "hydronium-cation", "hydroxide-anion"}) {
        throw py::value_error("Provider component order does not match");
    }
    if (charges != std::vector<int>{0, 1, -1}) {
        throw py::value_error("Provider charges do not match");
    }
    if (py::cast<std::vector<double>>(provider_binding["reaction_stoichiometry"])
        != std::vector<double>{-2.0, 1.0, 1.0}) {
        throw py::value_error("water self-ionization requires the exact 1:1 ionic reaction");
    }
    const double maximum_convergence_error = py::cast<double>(
        provider_binding["maximum_reference_convergence_error"]
    );
    if (!std::isfinite(maximum_convergence_error)
        || maximum_convergence_error < 0.0) {
        throw py::value_error("reference convergence bound is invalid");
    }

    const py::dict values = py::cast<py::dict>(record["values"]);
    const double p_kw = py::cast<double>(values["p_kw"]);
    const double ln_kw = py::cast<double>(values["ln_kw"]);
    if (!std::isfinite(p_kw) || !std::isfinite(ln_kw)
        || std::abs(ln_kw + std::log(10.0) * p_kw) > 1.0e-12) {
        throw py::value_error("IAPWS pKw and lnKw are inconsistent");
    }

    if (py::cast<std::string>(reference_record["basis_id"]) != provider_basis_id) {
        throw py::value_error("Provider standard-reference basis identity mismatch");
    }
    if (py::cast<std::vector<std::string>>(reference_record["component_ids"])
        != component_ids) {
        throw py::value_error("Provider standard-reference component order mismatch");
    }
    if (py::cast<std::vector<int>>(reference_record["charges"]) != charges) {
        throw py::value_error("Provider standard-reference charges mismatch");
    }
    if (py::cast<double>(reference_record["temperature_k"]) != temperature_k) {
        throw py::value_error("Provider standard-reference temperature binding mismatch");
    }
    if (py::cast<double>(reference_record["pressure_pa"]) != pressure_pa) {
        throw py::value_error("Provider standard-reference pressure binding mismatch");
    }

    StandardReferenceEvaluation reference;
    reference.formula_unit_log_fugacity = py::cast<double>(
        reference_record["formula_unit_log_fugacity"]
    );
    reference.pure_solvent_log_fugacity = py::cast<double>(
        reference_record["pure_solvent_log_fugacity"]
    );
    reference.solvent_molar_mass_kg_per_mol = py::cast<double>(
        reference_record["solvent_molar_mass_kg_per_mol"]
    );
    reference.reference_molality_mol_per_kg = py::cast<double>(
        reference_record["reference_molality_mol_per_kg"]
    );
    reference.convergence_error = py::cast<double>(reference_record["convergence_error"]);
    reference.basis_id = py::cast<std::string>(reference_record["basis_id"]);
    reference.parameter_fingerprint = py::cast<std::string>(
        reference_record["parameter_fingerprint"]
    );
    if (std::abs(reference.solvent_molar_mass_kg_per_mol - 0.01801528) > 1.0e-12) {
        throw py::value_error("Provider water molar mass does not match kg/mol identity");
    }
    if (!std::isfinite(reference.reference_molality_mol_per_kg)
        || reference.reference_molality_mol_per_kg <= 0.0) {
        throw py::value_error("Provider terminal molality must be finite and positive");
    }
    if (!std::isfinite(reference.convergence_error)
        || reference.convergence_error < 0.0
        || reference.convergence_error > maximum_convergence_error) {
        throw py::value_error("Provider standard-reference convergence is unacceptable");
    }

    const MixedStandardStateResult transformed =
        transform_water_self_ionization_standard_state(ln_kw, reference);
    const double independently_reconstructed = ln_kw
        + 2.0 * std::log(reference.solvent_molar_mass_kg_per_mol)
        + reference.formula_unit_log_fugacity
        - 2.0 * reference.pure_solvent_log_fugacity;
    const double transformation_residual = transformed.ln_k_provider_basis
        - independently_reconstructed;
    py::dict result;
    result["delta_standard_offset"] = transformed.delta_standard_offset;
    result["ln_k_provider_basis"] = transformed.ln_k_provider_basis;
    result["transformation_residual"] = transformation_residual;
    result["provider_basis_id"] = provider_basis_id;
    result["source_id"] = py::cast<std::string>(source["id"]);
    result["parameter_fingerprint"] = reference.parameter_fingerprint;
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
        solve_provider_manufactured_reaction(
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
        "_chemical_evaluate_provider_standard_reference",
        &provider_standard_reference_evidence,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("expected_fingerprint")
    );
    module.def(
        "_chemical_transform_water_self_ionization_standard_state",
        &transform_water_self_ionization_standard_state_evidence,
        py::arg("record"),
        py::arg("provider_reference")
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
