#include "chemical_equilibrium.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace epcsaft_equilibrium {
namespace {

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

py::dict compile_system(const py::dict& spec) {
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

    const CompiledReactionSystem compiled = compile_reaction_system(input);
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
}

}  // namespace epcsaft_equilibrium
