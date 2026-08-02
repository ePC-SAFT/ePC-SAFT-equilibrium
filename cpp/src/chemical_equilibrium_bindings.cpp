#include "chemical_equilibrium.hpp"
#include "chemical_observation.hpp"
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
constexpr std::size_t kChemicalNeutralReferenceSizeSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, neutral_reference_result_size)
    + sizeof(std::size_t);
constexpr std::size_t kChemicalNeutralReferenceDerivativeSizeSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, neutral_reference_derivative_result_size)
    + sizeof(std::size_t);
constexpr std::size_t kChemicalReactingPhaseParameterSizeSdkTableSize =
    offsetof(epcsaft_native_sdk_v1, reacting_phase_parameter_result_size)
    + sizeof(std::size_t);

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

std::vector<ProviderActiveParameterInput> active_parameter_inputs(
    const py::object& value
) {
    std::vector<ProviderActiveParameterInput> result;
    if (value.is_none()) {
        return result;
    }
    for (const py::handle item : py::cast<py::tuple>(value)) {
        const py::dict request = py::cast<py::dict>(item);
        result.push_back({
            py::cast<std::string>(request["family"]),
            py::cast<std::string>(request["identity"]),
            py::cast<std::vector<std::string>>(request["component_ids"]),
            py::cast<double>(request["value"]),
            py::cast<std::string>(request["unit"]),
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
    input.conserved_totals = py::cast<std::vector<double>>(spec["conserved_totals"]);
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

py::dict chemical_result(const ChemicalSolveResult& evaluation) {
    py::dict result;
    const auto optional_float = [](double value) -> py::object {
        return std::isfinite(value) ? py::cast(value) : py::none();
    };
    const auto optional_index = [](long value) -> py::object {
        return value >= 0 ? py::cast(value) : py::none();
    };
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
    result["negative_curvature_recovery_status"] =
        evaluation.negative_curvature_recovery_status;
    result["negative_curvature_recovery_attempts"] =
        evaluation.negative_curvature_recovery_attempts;
    result["negative_curvature_recovery_selected_sign"] =
        evaluation.negative_curvature_recovery_selected_sign;
    result["trace_status"] = evaluation.trace_status;
    result["globality_certificate"] = "not_guaranteed";
    result["amounts"] = evaluation.amounts;
    result["volume_m3"] = evaluation.volume_m3;
    result["balance_inf_norm"] = evaluation.balance_inf_norm;
    result["charge_inf_norm"] = evaluation.charge_inf_norm;
    result["pressure_relative_residual"] = evaluation.pressure_relative_residual;
    result["reaction_affinity_inf_norm"] = evaluation.reaction_affinity_inf_norm;
    result["reaction_affinity_residuals"] =
        evaluation.reaction_affinity_residuals;
    result["packing_fraction"] = evaluation.packing_fraction;
    result["packing_fraction_min"] = optional_float(
        evaluation.packing_fraction_min
    );
    result["packing_fraction_max"] = optional_float(
        evaluation.packing_fraction_max
    );
    result["total_ion_fraction"] = optional_float(
        evaluation.total_ion_fraction
    );
    result["total_ion_fraction_max"] = optional_float(
        evaluation.total_ion_fraction_max
    );
    result["minimum_amount_mol"] = optional_float(
        evaluation.minimum_amount_mol
    );
    result["trace_floor_mol"] = optional_float(evaluation.trace_floor_mol);
    result["kkt_stationarity_inf_norm"] = evaluation.kkt_stationarity_inf_norm;
    result["physical_stationarity_residuals"] =
        evaluation.physical_stationarity_residuals;
    result["complementarity_inf_norm"] = evaluation.complementarity_inf_norm;
    result["kkt_dimension"] = evaluation.kkt_dimension;
    result["kkt_rank"] = evaluation.kkt_rank;
    result["condition_number_inf"] =
        optional_float(evaluation.condition_number_inf);
    result["failure_kind"] = evaluation.failure_kind;
    result["failure_reason"] = evaluation.failure_reason;
    result["active_lower_bounds"] = evaluation.active_lower_bounds;
    result["active_upper_bounds"] = evaluation.active_upper_bounds;
    result["active_constraint_bounds"] =
        evaluation.active_constraint_bounds;
    result["reduced_hessian_status"] = evaluation.reduced_hessian_status;
    result["reduced_hessian"] = evaluation.reduced_hessian;
    result["reduced_hessian_nullspace_basis"] =
        evaluation.reduced_hessian_nullspace_basis;
    result["reduced_hessian_nullspace_shape"] = py::make_tuple(
        evaluation.reduced_hessian_nullspace_rows,
        evaluation.reduced_hessian_nullspace_columns
    );
    result["reduced_hessian_eigenvalues"] =
        evaluation.reduced_hessian_eigenvalues;
    result["reduced_hessian_spectrum_status"] =
        evaluation.reduced_hessian_spectrum_status;
    result["reduced_hessian_raw_inertia"] = py::make_tuple(
        evaluation.reduced_hessian_raw_positive_eigenvalues,
        evaluation.reduced_hessian_raw_zero_eigenvalues,
        evaluation.reduced_hessian_raw_negative_eigenvalues
    );
    result["reduced_hessian_inertia"] = py::make_tuple(
        evaluation.reduced_hessian_positive_eigenvalues,
        evaluation.reduced_hessian_zero_eigenvalues,
        evaluation.reduced_hessian_negative_eigenvalues
    );
    result["reduced_hessian_scale"] =
        optional_float(evaluation.reduced_hessian_scale);
    result["reduced_hessian_eigenvalue_tolerance"] =
        optional_float(evaluation.reduced_hessian_eigenvalue_tolerance);
    result["objective_gradient"] = evaluation.objective_gradient;
    result["constraint_values"] = evaluation.constraint_values;
    result["constraint_jacobian"] = evaluation.constraint_jacobian;
    result["lagrangian_gradient"] = evaluation.lagrangian_gradient;
    result["equality_multipliers"] = evaluation.equality_multipliers;
    result["chart_stationarity_inf_norm"] =
        optional_float(evaluation.chart_stationarity_inf_norm);
    result["lagrangian_hessian"] = evaluation.lagrangian_hessian;
    result["covariant_lagrangian_hessian"] =
        evaluation.covariant_lagrangian_hessian;
    const std::size_t derivative_dimension =
        evaluation.objective_gradient.size();
    result["derivative_coordinate_order"] =
        evaluation.derivative_coordinate_order.size() == derivative_dimension
        ? py::cast(evaluation.derivative_coordinate_order)
        : py::cast(std::vector<std::string>{});
    result["derivative_objective_basis"] =
        evaluation.derivative_objective_basis;
    result["derivative_constraint_basis"] =
        evaluation.derivative_constraint_basis;
    result["derivative_constraint_order"] =
        evaluation.derivative_constraint_order.size()
            == evaluation.constraint_values.size()
        ? py::cast(evaluation.derivative_constraint_order)
        : py::cast(std::vector<std::string>{});
    result["objective_gradient_shape"] = py::make_tuple(
        derivative_dimension
    );
    result["constraint_jacobian_shape"] = py::make_tuple(
        evaluation.constraint_values.size(), derivative_dimension
    );
    result["lagrangian_hessian_shape"] = py::make_tuple(
        derivative_dimension, derivative_dimension
    );
    result["kkt_root_jacobian"] = evaluation.kkt_root_jacobian;
    result["kkt_root_shape"] = py::make_tuple(
        evaluation.kkt_root_rows,
        evaluation.kkt_root_columns
    );
    result["kkt_root_status"] = evaluation.kkt_root_status;
    const auto criterion = [](
        const char* name,
        double value,
        double limit,
        const char* comparison,
        bool passed
    ) {
        py::dict record;
        record["name"] = name;
        record["status"] = std::isfinite(value) && passed ? "passed" : "failed";
        record["value"] = value;
        record["limit"] = limit;
        record["comparison"] = comparison;
        return record;
    };
    py::list numerical_criteria;
    numerical_criteria.append(criterion(
        "balance_inf_norm", evaluation.balance_inf_norm, 1.0e-9, "<=",
        evaluation.balance_inf_norm <= 1.0e-9
    ));
    numerical_criteria.append(criterion(
        "complementarity_inf_norm", evaluation.complementarity_inf_norm, 1.0e-7, "<=",
        evaluation.complementarity_inf_norm <= 1.0e-7
    ));
    numerical_criteria.append(criterion(
        "kkt_stationarity_inf_norm", evaluation.kkt_stationarity_inf_norm, 1.0e-7, "<=",
        evaluation.kkt_stationarity_inf_norm <= 1.0e-7
    ));
    numerical_criteria.append(criterion(
        "chart_stationarity_inf_norm", evaluation.chart_stationarity_inf_norm, 1.0e-7, "<=",
        evaluation.chart_stationarity_inf_norm <= 1.0e-7
    ));
    numerical_criteria.append(criterion(
        "derivative_evidence_finite", evaluation.callback_error.empty() ? 1.0 : 0.0,
        1.0, ">=", evaluation.callback_error.empty()
    ));
    result["numerical_criteria"] = std::move(numerical_criteria);
    py::list physical_criteria;
    physical_criteria.append(criterion(
        "balance_inf_norm", evaluation.balance_inf_norm, 1.0e-9, "<=",
        evaluation.balance_inf_norm <= 1.0e-9
    ));
    physical_criteria.append(criterion(
        "charge_inf_norm", evaluation.charge_inf_norm, 1.0e-9, "<=",
        evaluation.charge_inf_norm <= 1.0e-9
    ));
    physical_criteria.append(criterion(
        "pressure_relative_residual",
        evaluation.pressure_relative_residual,
        1.0e-8,
        "<=",
        evaluation.pressure_relative_residual <= 1.0e-8
    ));
    physical_criteria.append(criterion(
        "reaction_affinity_inf_norm",
        evaluation.reaction_affinity_inf_norm,
        1.0e-7,
        "<=",
        evaluation.reaction_affinity_inf_norm <= 1.0e-7
    ));
    if (std::isfinite(evaluation.minimum_amount_mol)
        && std::isfinite(evaluation.trace_floor_mol)) {
        physical_criteria.append(criterion(
            "minimum_amount_mol",
            evaluation.minimum_amount_mol,
            evaluation.trace_floor_mol,
            ">",
            evaluation.minimum_amount_mol > evaluation.trace_floor_mol
        ));
    }
    if (std::isfinite(evaluation.packing_fraction_min)
        && std::isfinite(evaluation.packing_fraction_max)) {
        physical_criteria.append(criterion(
            "packing_fraction_lower_bound",
            evaluation.packing_fraction,
            evaluation.packing_fraction_min,
            ">=",
            evaluation.packing_fraction >= evaluation.packing_fraction_min
        ));
        physical_criteria.append(criterion(
            "packing_fraction_upper_bound",
            evaluation.packing_fraction,
            evaluation.packing_fraction_max,
            "<=",
            evaluation.packing_fraction <= evaluation.packing_fraction_max
        ));
    }
    if (std::isfinite(evaluation.total_ion_fraction_max)) {
        physical_criteria.append(criterion(
            "total_ion_fraction",
            evaluation.total_ion_fraction,
            evaluation.total_ion_fraction_max + 1.0e-12,
            "<=",
            evaluation.total_ion_fraction
                <= evaluation.total_ion_fraction_max + 1.0e-12
        ));
    }
    result["physical_criteria"] = std::move(physical_criteria);
    py::dict search;
    search["status"] = evaluation.search.status;
    search["continuation_status"] = evaluation.search.continuation_status;
    search["primary_budget"] = evaluation.search.primary_budget;
    search["primary_attempt_count"] =
        evaluation.search.primary_attempt_count;
    search["selected_basin_ordinal"] =
        optional_index(evaluation.search.selected_basin_ordinal);
    search["selected_objective"] =
        optional_float(evaluation.search.selected_objective);
    search["selection_label"] = evaluation.search.selection_label;
    py::list attempts;
    for (const ChemicalSearchAttempt& attempt : evaluation.search.attempts) {
        py::dict record;
        record["ordinal"] = attempt.ordinal;
        record["primary_ordinal"] = attempt.primary_ordinal;
        record["kind"] = attempt.kind;
        record["parent_ordinal"] = optional_index(attempt.parent_ordinal);
        record["start_identity"] = attempt.start_identity;
        record["start_construction_status"] =
            attempt.start_construction_status;
        record["retraction_status"] = attempt.retraction_status;
        record["continuation_status"] = attempt.continuation_status;
        record["provider_domain_status"] =
            attempt.provider_domain_status;
        record["solver_status"] = attempt.solver_status;
        record["callback_error"] = attempt.callback_error;
        record["terminal_status"] = attempt.terminal_status;
        record["amounts"] = attempt.amounts;
        record["volume_m3"] = optional_float(attempt.volume_m3);
        record["objective"] = optional_float(attempt.objective);
        record["balance_inf_norm"] =
            optional_float(attempt.balance_inf_norm);
        record["charge_inf_norm"] =
            optional_float(attempt.charge_inf_norm);
        record["pressure_relative_residual"] =
            optional_float(attempt.pressure_relative_residual);
        record["reaction_affinity_inf_norm"] =
            optional_float(attempt.reaction_affinity_inf_norm);
        record["kkt_stationarity_inf_norm"] =
            optional_float(attempt.kkt_stationarity_inf_norm);
        record["complementarity_inf_norm"] =
            optional_float(attempt.complementarity_inf_norm);
        record["kkt_dimension"] = attempt.kkt_dimension;
        record["kkt_rank"] = attempt.kkt_rank;
        record["condition_number_inf"] =
            optional_float(attempt.condition_number_inf);
        record["local_minimum_status"] =
            attempt.local_minimum_status;
        record["trace_status"] = attempt.trace_status;
        record["basin_ordinal"] = optional_index(attempt.basin_ordinal);
        record["recovery_seed_count"] = attempt.recovery_seed_count;
        record["recovery_solve_count"] = attempt.recovery_solve_count;
        attempts.append(std::move(record));
    }
    search["attempts"] = std::move(attempts);
    py::list basins;
    for (const ChemicalSearchBasin& basin : evaluation.search.basins) {
        py::dict record;
        record["ordinal"] = basin.ordinal;
        record["representative_attempt_ordinal"] =
            basin.representative_attempt_ordinal;
        record["amounts"] = basin.amounts;
        record["volume_m3"] = optional_float(basin.volume_m3);
        record["objective"] = optional_float(basin.objective);
        basins.append(std::move(record));
    }
    search["basins"] = std::move(basins);
    py::list prefixes;
    for (const ChemicalSearchBudgetPrefix& prefix
         : evaluation.search.budget_prefixes) {
        py::dict record;
        record["primary_budget"] = prefix.primary_budget;
        record["attempted_primary_ordinals"] =
            prefix.attempted_primary_ordinals;
        record["basin_ordinals"] = prefix.basin_ordinals;
        record["selected_basin_ordinal"] =
            optional_index(prefix.selected_basin_ordinal);
        record["selection_changed"] = prefix.selection_changed;
        prefixes.append(std::move(record));
    }
    search["budget_prefixes"] = std::move(prefixes);
    result["search"] = std::move(search);
    py::dict sensitivities;
    sensitivities["status"] = evaluation.sensitivities.status;
    sensitivities["failure_reason"] = evaluation.sensitivities.failure_reason;
    const std::size_t species_count = evaluation.amounts.size();
    if (evaluation.sensitivities.amount_derivatives.size()
            != evaluation.sensitivities.parameter_order.size() * species_count
        || evaluation.sensitivities.volume_derivatives.size()
            != evaluation.sensitivities.parameter_order.size()) {
        throw std::runtime_error(
            "internal chemical-sensitivity result dimensions are inconsistent"
        );
    }
    py::tuple parameter_order(evaluation.sensitivities.parameter_order.size());
    for (std::size_t parameter = 0;
         parameter < evaluation.sensitivities.parameter_order.size();
         ++parameter) {
        parameter_order[parameter] =
            evaluation.sensitivities.parameter_order[parameter];
    }
    sensitivities["parameter_order"] = std::move(parameter_order);
    py::tuple amount_derivatives(evaluation.sensitivities.parameter_order.size());
    for (std::size_t parameter = 0;
         parameter < evaluation.sensitivities.parameter_order.size();
         ++parameter) {
        py::tuple row(species_count);
        for (std::size_t species = 0; species < species_count; ++species) {
            row[species] = evaluation.sensitivities.amount_derivatives[
                parameter * species_count + species
            ];
        }
        amount_derivatives[parameter] = std::move(row);
    }
    sensitivities["amount_derivatives"] = std::move(amount_derivatives);
    sensitivities["volume_derivatives"] = py::tuple(py::cast(
        evaluation.sensitivities.volume_derivatives
    ));
    sensitivities["kkt_dimension"] = evaluation.sensitivities.kkt_dimension;
    sensitivities["kkt_rank"] = evaluation.sensitivities.kkt_rank;
    sensitivities["condition_number_inf"] =
        evaluation.sensitivities.condition_number_inf;
    sensitivities["active_lower_bounds"] = py::tuple(py::cast(
        evaluation.sensitivities.active_lower_bounds
    ));
    sensitivities["active_upper_bounds"] = py::tuple(py::cast(
        evaluation.sensitivities.active_upper_bounds
    ));
    sensitivities["active_constraint_bounds"] = py::tuple(py::cast(
        evaluation.sensitivities.active_constraint_bounds
    ));
    sensitivities["active_trace_species"] = py::tuple(py::cast(
        evaluation.sensitivities.active_trace_species
    ));
    sensitivities["chart_topology"] =
        evaluation.sensitivities.chart_topology;
    sensitivities["parameter_fingerprint"] =
        evaluation.sensitivities.parameter_fingerprint;
    sensitivities["provider_parameter_status"] =
        evaluation.sensitivities.provider_parameter_status;
    sensitivities["provider_parameter_failure_reason"] =
        evaluation.sensitivities.provider_parameter_failure_reason;
    sensitivities["reference_parameter_status"] =
        evaluation.sensitivities.reference_parameter_status;
    sensitivities["reference_parameter_failure_reason"] =
        evaluation.sensitivities.reference_parameter_failure_reason;
    result["sensitivities"] = std::move(sensitivities);
    return result;
}

py::dict solve_manufactured(const py::dict& spec, double trace_floor) {
    const ReactionSystemInput input = reaction_system_input(spec);
    const CompiledReactionSystem compiled = compile_reaction_system(input);
    return chemical_result(
        solve_manufactured_ideal_reaction(
            compiled,
            input.temperature_k,
            input.pressure_pa,
            {},
            trace_floor
        )
    );
}

py::dict solve_manufactured_nonconvex(
    const py::dict& spec,
    double trace_floor,
    int max_iterations,
    double quadratic_strength
) {
    const ReactionSystemInput input = reaction_system_input(spec);
    const CompiledReactionSystem compiled = compile_reaction_system(input);
    return chemical_result(
        solve_manufactured_nonconvex_reaction(
            compiled,
            input.temperature_k,
            input.pressure_pa,
            {},
            trace_floor,
            max_iterations,
            quadratic_strength
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
    result["component_ids"] = metadata.component_ids;
    result["coordinate_order"] = py::make_tuple(
        "component_amounts_mol_in_provider_order", "volume_m3"
    );
    result["density_transformation"] = "rho_i_mol_per_m3=amount_i_mol/volume_m3";
    result["helmholtz_basis"] = "dimensionless_A_over_RT";
    result["reference_transformation"] =
        "provider_phase_block_only; source-standard offsets are linear external terms";
    result["value"] = evaluation.value;
    result["gradient"] = evaluation.gradient;
    result["hessian"] = evaluation.hessian;
    result["pressure_pa"] = evaluation.pressure_pa;
    result["packing_fraction"] = evaluation.packing_fraction;
    result["packing_gradient"] = evaluation.packing_gradient;
    result["packing_hessian"] = evaluation.packing_hessian;
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
    const std::vector<double>& packing_bounds,
    double trace_floor,
    const std::vector<double>& ln_k_pressure_derivatives_per_pa = {},
    const std::vector<double>& ln_k_parameter_derivatives = {},
    const ProviderActiveParameterSet* active_parameters = nullptr
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
    if (packing_bounds.size() != 2) {
        throw py::value_error("packing_fraction_bounds must contain two values");
    }
    std::vector<double> compiled_ln_k_pressure_derivatives_per_pa;
    if (!ln_k_pressure_derivatives_per_pa.empty()) {
        if (ln_k_pressure_derivatives_per_pa.size()
                != input.reaction_matrix.rows
            || compiled.supplied_reaction_transform.rows
                != compiled.reaction_matrix.rows
            || compiled.supplied_reaction_transform.columns
                != input.reaction_matrix.rows) {
            throw py::value_error(
                "source pressure derivatives do not match the compiled reaction basis"
            );
        }
        compiled_ln_k_pressure_derivatives_per_pa.assign(
            compiled.reaction_matrix.rows, 0.0
        );
        for (std::size_t reaction = 0;
             reaction < compiled.reaction_matrix.rows;
             ++reaction) {
            for (std::size_t supplied = 0;
                 supplied < input.reaction_matrix.rows;
                 ++supplied) {
                compiled_ln_k_pressure_derivatives_per_pa[reaction] +=
                    compiled.supplied_reaction_transform(reaction, supplied)
                    * ln_k_pressure_derivatives_per_pa[supplied];
            }
        }
    }
    std::vector<double> compiled_ln_k_parameter_derivatives;
    const std::size_t active_parameter_count = active_parameters == nullptr
        ? 0
        : active_parameters->parameters.size();
    if (!ln_k_parameter_derivatives.empty()) {
        if (active_parameter_count == 0
            || ln_k_parameter_derivatives.size()
                != input.reaction_matrix.rows * active_parameter_count
            || compiled.supplied_reaction_transform.rows
                != compiled.reaction_matrix.rows
            || compiled.supplied_reaction_transform.columns
                != input.reaction_matrix.rows) {
            throw py::value_error(
                "source parameter derivatives do not match the compiled reaction basis"
            );
        }
        compiled_ln_k_parameter_derivatives.assign(
            compiled.reaction_matrix.rows * active_parameter_count, 0.0
        );
        for (std::size_t reaction = 0;
             reaction < compiled.reaction_matrix.rows;
             ++reaction) {
            for (std::size_t supplied = 0;
                 supplied < input.reaction_matrix.rows;
                 ++supplied) {
                for (std::size_t parameter = 0;
                     parameter < active_parameter_count;
                     ++parameter) {
                    compiled_ln_k_parameter_derivatives[
                        reaction * active_parameter_count + parameter
                    ] += compiled.supplied_reaction_transform(reaction, supplied)
                        * ln_k_parameter_derivatives[
                            supplied * active_parameter_count + parameter
                        ];
                }
            }
        }
    } else if (active_parameter_count != 0) {
        compiled_ln_k_parameter_derivatives.assign(
            compiled.reaction_matrix.rows * active_parameter_count, 0.0
        );
    }
    ChemicalSolveResult evaluation = solve_provider_reaction(
        compiled,
        provider,
        input.temperature_k,
        input.pressure_pa,
        packing_bounds[0],
        packing_bounds[1],
        sdk.total_ion_mole_fraction_max,
        trace_floor,
        compiled_ln_k_pressure_derivatives_per_pa,
        compiled_ln_k_parameter_derivatives,
        active_parameters
    );
    if (active_parameter_count == 0) {
        evaluation.sensitivities.provider_parameter_status = "not_applicable";
        evaluation.sensitivities.provider_parameter_failure_reason.clear();
    } else if (evaluation.sensitivities.status == "available") {
        evaluation.sensitivities.provider_parameter_status = "available";
        evaluation.sensitivities.provider_parameter_failure_reason.clear();
    } else {
        evaluation.sensitivities.provider_parameter_status = "unavailable";
        evaluation.sensitivities.provider_parameter_failure_reason =
            evaluation.sensitivities.failure_reason;
    }
    py::dict result = chemical_result(evaluation);
    result["parameter_fingerprint"] = input.provider_fingerprint;
    result["packing_fraction_bounds"] = packing_bounds;
    result["provider_sdk_capsule_name"] = EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME;
    result["provider_sdk_abi_version"] = sdk.abi_version;
    result["provider_sdk_table_size"] = sdk.table_size;
    result["provider_sdk_result_size"] = sdk.result_size;
    result["provider_sdk_mixture_result_size"] = sdk.mixture_result_size;
    result["provider_sdk_neutral_reference_result_size"] =
        sdk.table_size >= kChemicalNeutralReferenceSizeSdkTableSize
        ? sdk.neutral_reference_result_size
        : 0;
    result["provider_sdk_neutral_reference_derivative_result_size"] =
        sdk.table_size >= kChemicalNeutralReferenceDerivativeSizeSdkTableSize
        ? sdk.neutral_reference_derivative_result_size
        : 0;
    result["provider_sdk_reacting_phase_parameter_result_size"] =
        sdk.table_size >= kChemicalReactingPhaseParameterSizeSdkTableSize
        ? sdk.reacting_phase_parameter_result_size
        : 0;
    return result;
}

py::dict solve_provider_source(
    const py::capsule& capsule,
    const py::dict& spec,
    const py::dict& source_standard_state,
    const std::vector<double>& packing_bounds,
    double trace_floor,
    bool sensitivities_requested,
    const ProviderActiveParameterSet* active_parameters
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
    const double source_reference_pressure_pa = py::cast<double>(
        source_standard_state["reference_pressure_pa"]
    );
    if (source_standard_state_id.empty() || source_activity_scale_id.empty()) {
        throw py::value_error("source standard-state identity is incomplete");
    }
    if (!std::isfinite(source_reference_pressure_pa)
        || source_reference_pressure_pa <= 0.0) {
        throw py::value_error(
            "source standard-state reference pressure is invalid"
        );
    }
    for (const EquilibriumConstantRecord& record : input.equilibrium_constant_records) {
        if (record.conversion_id
                != "source-standard-state-to-provider-neutral-reference"
            || record.reaction_orientation != "products_positive"
            || !record.dimensionless
            || record.reference_id.empty()
            || record.temperature_k != input.temperature_k
            || record.pressure_pa != source_reference_pressure_pa) {
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
    const NeutralReferenceEvaluation reference = sensitivities_requested
        ? provider.evaluate_neutral_reference_derivatives(
            input.temperature_k,
            input.pressure_pa,
            active_parameters == nullptr
                ? ProviderActiveParameterSet{}
                : *active_parameters
        )
        : provider.evaluate_neutral_reference(
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
        record.pressure_pa = input.pressure_pa;
    }
    py::dict result = solve_provider_input(
        sdk,
        metadata,
        input,
        provider,
        packing_bounds,
        trace_floor,
        sensitivities_requested
            ? transformed.pressure_derivatives_per_pa
            : std::vector<double>{},
        sensitivities_requested
            ? transformed.parameter_derivatives
            : std::vector<double>{},
        active_parameters
    );
    py::dict sensitivities = py::cast<py::dict>(result["sensitivities"]);
    sensitivities["reference_parameter_status"] =
        sensitivities_requested ? "available" : "not_applicable";
    sensitivities["reference_parameter_failure_reason"] = "";
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

py::dict solve_chemical_equilibrium(
    const py::object& capsule,
    const py::dict& spec,
    const py::object& source_standard_state,
    const py::object& packing_fraction_bounds,
    double trace_floor,
    const py::object& active_parameters
) {
    const std::vector<ProviderActiveParameterInput> active_requests =
        active_parameter_inputs(active_parameters);
    if (capsule.is_none()) {
        if (!source_standard_state.is_none() || !packing_fraction_bounds.is_none()) {
            throw py::value_error(
                "ideal-gas chemical equilibrium cannot consume Provider metadata"
            );
        }
        return solve_manufactured(spec, trace_floor);
    }
    if (packing_fraction_bounds.is_none()) {
        throw py::value_error(
            "packing_fraction_bounds are required from the calling formulation"
        );
    }
    const std::vector<double> packing_bounds =
        py::cast<std::vector<double>>(packing_fraction_bounds);
    const py::capsule provider_capsule = py::cast<py::capsule>(capsule);
    if (source_standard_state.is_none()) {
        const epcsaft_native_sdk_v1& sdk = checked_chemical_sdk(provider_capsule);
        const ChemicalProviderMetadata metadata = chemical_provider_metadata(sdk);
        const ReactionSystemInput input = reaction_system_input(spec);
        const ProviderContext provider(sdk, input.provider_fingerprint);
        const ProviderActiveParameterSet resolved =
            provider.resolve_active_parameters(input.temperature_k, active_requests);
        return solve_provider_input(
            sdk,
            metadata,
            input,
            provider,
            packing_bounds,
            trace_floor,
            {},
            {},
            resolved.parameters.empty() ? nullptr : &resolved
        );
    }
    const epcsaft_native_sdk_v1& sdk = checked_chemical_sdk(provider_capsule);
    const ReactionSystemInput input = reaction_system_input(spec);
    const ProviderContext provider(sdk, input.provider_fingerprint);
    const ProviderActiveParameterSet resolved =
        provider.resolve_active_parameters(input.temperature_k, active_requests);
    return solve_provider_source(
        provider_capsule,
        spec,
        py::cast<py::dict>(source_standard_state),
        packing_bounds,
        trace_floor,
        !active_parameters.is_none(),
        resolved.parameters.empty() ? nullptr : &resolved
    );
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

py::dict manufactured_inverse_log_packing_nlp_evidence(
    const py::dict& spec,
    const std::vector<double>& variables,
    const std::vector<double>& constraint_multipliers,
    const std::vector<double>& gauge_coefficients,
    bool zero_kkt_rhs
) {
    const ReactionSystemInput input = reaction_system_input(spec);
    const ManufacturedNlpEvaluation evaluation =
        evaluate_manufactured_inverse_log_packing_nlp(
            compile_reaction_system(input),
            input.temperature_k,
            input.pressure_pa,
            gauge_coefficients,
            variables,
            constraint_multipliers,
            zero_kkt_rhs
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
    result["kkt_backtransform_rhs"] = evaluation.kkt_backtransform_rhs;
    result["kkt_backtransform_solution"] = evaluation.kkt_backtransform_solution;
    return result;
}

py::dict manufactured_reduced_hessian_evidence(
    const std::vector<double>& hessian,
    const std::vector<double>& constraint_jacobian,
    std::size_t constraint_count
) {
    const ManufacturedReducedHessianEvidence evidence =
        analyze_manufactured_reduced_hessian(
            hessian, constraint_jacobian, constraint_count
        );
    py::dict result;
    result["positive"] = evidence.positive;
    result["status"] = evidence.status;
    result["curvature"] = evidence.curvature;
    result["negative_direction"] = evidence.negative_direction;
    result["reduced_hessian"] = evidence.reduced_hessian;
    result["nullspace_basis"] = evidence.nullspace_basis;
    result["nullspace_shape"] = py::make_tuple(
        evidence.nullspace_rows, evidence.nullspace_columns
    );
    result["eigenvalues"] = evidence.eigenvalues;
    result["inertia"] = evidence.inertia;
    result["raw_inertia"] = evidence.raw_inertia;
    result["spectrum_status"] = evidence.spectrum_status;
    result["hessian_scale"] = evidence.hessian_scale;
    result["eigenvalue_tolerance"] = evidence.eigenvalue_tolerance;
    return result;
}

std::vector<double> manufactured_balance_retraction_evidence(
    const py::dict& spec,
    const std::vector<double>& seed,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    double trace_floor
) {
    const ReactionSystemInput input = reaction_system_input(spec);
    return retract_manufactured_balance(
        compile_reaction_system(input),
        seed,
        lower,
        upper,
        trace_floor
    );
}

std::vector<double> manufactured_recovery_displacement_evidence(
    const std::vector<double>& variables,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    const std::vector<double>& direction,
    int sign,
    std::size_t backtrack_index
) {
    return manufactured_recovery_displacement(
        variables, lower, upper, direction, sign, backtrack_index
    );
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
        "_chemical_equilibrium",
        &solve_chemical_equilibrium,
        py::arg("capsule"),
        py::arg("spec"),
        py::arg("source_standard_state"),
        py::arg("packing_fraction_bounds"),
        py::arg("trace_floor"),
        py::arg("active_parameters") = py::none()
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
        "_chemical_evaluate_manufactured_nlp",
        &manufactured_nlp_evidence,
        py::arg("spec"),
        py::arg("variables"),
        py::arg("constraint_multipliers"),
        py::arg("gauge_coefficients") = std::vector<double>{}
    );
    module.def(
        "_chemical_evaluate_manufactured_inverse_log_packing_nlp",
        &manufactured_inverse_log_packing_nlp_evidence,
        py::arg("spec"),
        py::arg("variables"),
        py::arg("constraint_multipliers"),
        py::arg("gauge_coefficients") = std::vector<double>{},
        py::arg("zero_kkt_rhs") = false
    );
    module.def(
        "_chemical_analyze_manufactured_reduced_hessian",
        &manufactured_reduced_hessian_evidence,
        py::arg("hessian"),
        py::arg("constraint_jacobian") = std::vector<double>{},
        py::arg("constraint_count") = 0
    );
    module.def(
        "_chemical_retract_manufactured_balance",
        &manufactured_balance_retraction_evidence,
        py::arg("spec"),
        py::arg("seed"),
        py::arg("lower"),
        py::arg("upper"),
        py::arg("trace_floor")
    );
    module.def(
        "_chemical_manufactured_recovery_displacement",
        &manufactured_recovery_displacement_evidence,
        py::arg("variables"),
        py::arg("lower"),
        py::arg("upper"),
        py::arg("direction"),
        py::arg("sign"),
        py::arg("backtrack_index") = 0
    );
    module.def(
        "_chemical_solve_manufactured_nonconvex",
        &solve_manufactured_nonconvex,
        py::arg("spec"),
        py::arg("trace_floor") = 1.0e-12,
        py::arg("max_iterations") = 500,
        py::arg("quadratic_strength") = 2.3
    );
    bind_chemical_observation(module);
}

}  // namespace epcsaft_equilibrium
