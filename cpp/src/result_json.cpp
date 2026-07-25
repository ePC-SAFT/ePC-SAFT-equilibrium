#include "result_json.hpp"

#include "held2_tolerances.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace epcsaft_equilibrium {
namespace {

class JsonValue {
public:
    using Array = std::vector<JsonValue>;
    using Object = std::vector<std::pair<std::string, JsonValue>>;

    JsonValue() : value_(nullptr) {}
    JsonValue(std::nullptr_t) : value_(nullptr) {}
    JsonValue(bool value) : value_(value) {}
    JsonValue(double value) : value_(value) {}
    JsonValue(const char* value) : value_(std::string(value)) {}
    JsonValue(std::string value) : value_(std::move(value)) {}
    JsonValue(std::string_view value) : value_(std::string(value)) {}

    template <
        typename Integer,
        std::enable_if_t<
            std::is_integral_v<Integer>
                && !std::is_same_v<std::remove_cv_t<Integer>, bool>,
            int
        > = 0
    >
    JsonValue(Integer value) : value_(static_cast<std::int64_t>(value)) {}

    template <typename Value>
    JsonValue(const std::vector<Value>& values) : value_(Array{}) {
        Array& array = std::get<Array>(value_);
        array.reserve(values.size());
        for (const Value& value : values) {
            array.emplace_back(value);
        }
    }

    template <typename Value, std::size_t Size>
    JsonValue(const std::array<Value, Size>& values) : value_(Array{}) {
        Array& array = std::get<Array>(value_);
        array.reserve(Size);
        for (const Value& value : values) {
            array.emplace_back(value);
        }
    }

    static JsonValue array() {
        return JsonValue(Array{});
    }

    static JsonValue object() {
        return JsonValue(Object{});
    }

    JsonValue& operator[](std::string key) {
        Object& object = std::get<Object>(value_);
        for (auto& [existing_key, existing_value] : object) {
            if (existing_key == key) {
                return existing_value;
            }
        }
        object.emplace_back(std::move(key), JsonValue{});
        return object.back().second;
    }

    void append(JsonValue value) {
        std::get<Array>(value_).push_back(std::move(value));
    }

    [[nodiscard]] const auto& value() const {
        return value_;
    }

private:
    explicit JsonValue(Array value) : value_(std::move(value)) {}
    explicit JsonValue(Object value) : value_(std::move(value)) {}

    std::variant<
        std::nullptr_t,
        bool,
        std::int64_t,
        double,
        std::string,
        Array,
        Object
    > value_;
};

void write_json_string(std::ostream& output, std::string_view value) {
    output << '"';
    for (unsigned char character : value) {
        switch (character) {
        case '"':
            output << "\\\"";
            break;
        case '\\':
            output << "\\\\";
            break;
        case '\b':
            output << "\\b";
            break;
        case '\f':
            output << "\\f";
            break;
        case '\n':
            output << "\\n";
            break;
        case '\r':
            output << "\\r";
            break;
        case '\t':
            output << "\\t";
            break;
        default:
            if (character < 0x20) {
                output << "\\u"
                       << std::hex << std::setw(4) << std::setfill('0')
                       << static_cast<int>(character)
                       << std::dec << std::setw(0) << std::setfill(' ');
            } else {
                output << static_cast<char>(character);
            }
        }
    }
    output << '"';
}

void write_json(std::ostream& output, const JsonValue& value) {
    std::visit(
        [&output](const auto& item) {
            using Value = std::decay_t<decltype(item)>;
            if constexpr (std::is_same_v<Value, std::nullptr_t>) {
                output << "null";
            } else if constexpr (std::is_same_v<Value, bool>) {
                output << (item ? "true" : "false");
            } else if constexpr (std::is_same_v<Value, std::int64_t>) {
                output << item;
            } else if constexpr (std::is_same_v<Value, double>) {
                if (std::isfinite(item)) {
                    output << std::setprecision(
                        std::numeric_limits<double>::max_digits10
                    ) << item;
                } else {
                    output << "null";
                }
            } else if constexpr (std::is_same_v<Value, std::string>) {
                write_json_string(output, item);
            } else if constexpr (std::is_same_v<Value, JsonValue::Array>) {
                output << '[';
                for (std::size_t index = 0; index < item.size(); ++index) {
                    if (index != 0) {
                        output << ',';
                    }
                    write_json(output, item[index]);
                }
                output << ']';
            } else {
                output << '{';
                for (std::size_t index = 0; index < item.size(); ++index) {
                    if (index != 0) {
                        output << ',';
                    }
                    write_json_string(output, item[index].first);
                    output << ':';
                    write_json(output, item[index].second);
                }
                output << '}';
            }
        },
        value.value()
    );
}

JsonValue held_stage_ii_trace_to_json(const HeldStageIITrace& entry) {
    JsonValue rejections = JsonValue::array();
    for (const auto& rejection : entry.rejections) {
        JsonValue item = JsonValue::object();
        item["identity"] = rejection.identity;
        item["reason"] = rejection.reason;
        rejections.append(std::move(item));
    }
    JsonValue result = JsonValue::object();
    result["major_iteration"] = entry.major_iteration;
    result["outer_value"] = entry.outer_value;
    result["upper_bound"] = entry.upper_bound;
    result["multiplier"] = entry.multiplier;
    result["active_cut_ids"] = entry.active_cut_ids;
    result["accepted_cut_ids"] = entry.accepted_cut_ids;
    result["lower_starts_completed"] = entry.lower_starts_completed;
    result["candidate_ids"] = entry.candidate_ids;
    result["rejections"] = std::move(rejections);
    result["stage_iii_outcome"] = entry.stage_iii_outcome;
    result["stage_iii_failure_reason"] = entry.stage_iii_failure_reason;
    return result;
}

JsonValue held_stage_iii_initialization_to_json(
    const HeldStageIIIInitialization& initialization
) {
    JsonValue result = JsonValue::object();
    result["status"] = initialization.status;
    result["failure_reason"] = initialization.failure_reason;
    result["phase_fractions"] = initialization.phase_fractions;
    result["initial_variables"] = initialization.has_variables
        ? JsonValue(initialization.initial_variables)
        : JsonValue(std::vector<double>{});
    result["lower_bounds"] = initialization.has_variables
        ? JsonValue(initialization.lower_bounds)
        : JsonValue(std::vector<double>{});
    result["upper_bounds"] = initialization.has_variables
        ? JsonValue(initialization.upper_bounds)
        : JsonValue(std::vector<double>{});
    return result;
}

JsonValue held_stage_iii_phase_to_json(
    const HeldStageIIIPhaseEvaluation& phase
) {
    const double amount_mol = phase.amounts_mol[0] + phase.amounts_mol[1];
    JsonValue result = JsonValue::object();
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

JsonValue held_stage_iii_attempt_log_to_json(
    const std::vector<HeldStageIIILocalAttempt>& attempts
) {
    JsonValue result = JsonValue::array();
    for (const auto& attempt : attempts) {
        JsonValue item = JsonValue::object();
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

JsonValue held_stage_iii_to_json(const HeldStageIIIResult& solve) {
    JsonValue result = JsonValue::object();
    result["outcome"] = solve.outcome;
    result["search_profile"] = solve.search_profile;
    result["failure_reason"] = solve.failure_reason;
    result["initialization"] =
        held_stage_iii_initialization_to_json(solve.initialization);
    result["solver_converged"] = solve.local.solver_converged;
    result["solver_status"] = solve.local.solver_status;
    result["iterations"] = solve.local.iterations;
    result["attempts"] = solve.local.attempts;
    result["attempt_log"] =
        held_stage_iii_attempt_log_to_json(solve.local.attempt_log);
    result["solver_lower_bounds"] = solve.local.solver_lower_bounds;
    result["solver_upper_bounds"] = solve.local.solver_upper_bounds;
    result["solver_constraint_violation"] =
        solve.local.solver_constraint_violation;
    result["confirmation_solves"] = solve.local.confirmation_solves;
    result["confirmation_succeeded"] = solve.confirmation_succeeded;
    result["confirmation_max_difference"] =
        solve.local.confirmation_max_difference;
    result["material_balance_max_abs"] = solve.local.material_balance_max_abs;
    result["pressure_stationarity_max_relative"] =
        solve.local.pressure_stationarity_max_relative;
    result["kkt_stationarity_max_abs"] = solve.local.kkt_stationarity_max_abs;
    result["inactive_bounds"] = solve.inactive_bounds;
    result["composition_distance"] = solve.composition_distance;
    result["phase_density_distance"] = solve.local.phase_density_distance;
    result["chemical_potential_max_relative"] =
        solve.chemical_potential_max_relative;
    result["held_gap"] = solve.held_gap;
    result["upper_bound"] = solve.upper_bound;
    result["equality_multipliers"] = solve.local.equality_multipliers;
    result["lower_bound_multipliers"] = solve.local.lower_bound_multipliers;
    result["upper_bound_multipliers"] = solve.local.upper_bound_multipliers;
    result["total_g_bar"] = solve.local.evaluation.objective;
    if (solve.local.accepted) {
        result["first_phase"] =
            held_stage_iii_phase_to_json(solve.local.evaluation.liquid);
        result["second_phase"] =
            held_stage_iii_phase_to_json(solve.local.evaluation.vapor);
    }
    return result;
}

JsonValue neutral_flash_to_json(const FlashResult& flash) {
    const HeldResult& solve = flash.held;
    int attempt_count = static_cast<int>(
        solve.stage_i.reference_attempts.size()
        + solve.stage_i.attempt_log.size()
        + solve.stage_ii.endpoint_attempts.size()
    );
    for (const auto& entry : solve.stage_ii.trace) {
        attempt_count += entry.lower_starts_completed;
    }
    JsonValue attempts = JsonValue::array();
    for (const auto& attempt : solve.stage_iii_attempts) {
        attempt_count += attempt.result.local.attempts;
        JsonValue item = held_stage_iii_to_json(attempt.result);
        item["major_iteration"] = attempt.major_iteration;
        attempts.append(std::move(item));
    }
    JsonValue trace = JsonValue::array();
    for (const auto& entry : solve.stage_ii.trace) {
        trace.append(held_stage_ii_trace_to_json(entry));
    }
    JsonValue search_profiles = JsonValue::array();
    search_profiles.append(solve.stage_i.search_profile);
    if (solve.stage_i.outcome == "negative_tpd") {
        search_profiles.append(solve.stage_ii.search_profile);
    }
    if (!solve.stage_iii_attempts.empty()) {
        search_profiles.append(
            solve.stage_iii_attempts.front().result.search_profile
        );
    }

    JsonValue result = JsonValue::object();
    result["outcome"] = solve.outcome;
    result["search_status"] = solve.search_status;
    result["failure_reason"] = solve.failure_reason;
    result["stage_i_outcome"] = solve.stage_i.outcome;
    result["stage_i_search_status"] = solve.stage_i.search_status;
    result["attempts"] = attempt_count;
    result["major_iterations"] = solve.stage_ii.major_iterations;
    result["best_tpd"] = solve.stage_i.best_tpd;
    result["lower_bound"] = nullptr;
    result["upper_bound"] = nullptr;
    result["held_gap"] = nullptr;
    result["material_balance_max_abs"] = nullptr;
    result["pressure_stationarity_max_relative"] = nullptr;
    result["kkt_stationarity_max_abs"] = nullptr;
    result["chemical_potential_max_relative"] = nullptr;
    result["confirmation_succeeded"] = false;
    result["confirmation_max_difference"] = nullptr;
    result["search_profiles"] = std::move(search_profiles);
    result["stage_iii_attempts"] = std::move(attempts);
    result["trace"] = std::move(trace);
    result["globality_certificate"] = flash.globality_certificate;
    result["solver_status"] = "not_adjudicated";
    result["numerical_status"] = "not_adjudicated";
    result["physical_status"] = "not_adjudicated";

    const bool accepted_reference = solve.stage_i.has_reference
        && std::any_of(
            solve.stage_i.reference_attempts.begin(),
            solve.stage_i.reference_attempts.end(),
            [](const auto& attempt) { return attempt.accepted; }
        );
    const bool failed_reference_solver_attempt = std::any_of(
        solve.stage_i.reference_attempts.begin(),
        solve.stage_i.reference_attempts.end(),
        [](const auto& attempt) {
            return !attempt.solver_converged
                || !attempt.callback_error.empty();
        }
    );
    const bool failed_reference_attempt = std::any_of(
        solve.stage_i.reference_attempts.begin(),
        solve.stage_i.reference_attempts.end(),
        [](const auto& attempt) { return !attempt.accepted; }
    );
    const bool failed_stage_i_solver_attempt = std::any_of(
        solve.stage_i.attempt_log.begin(),
        solve.stage_i.attempt_log.end(),
        [](const auto& attempt) {
            return !attempt.solver_converged
                || !attempt.callback_error.empty();
        }
    );
    const bool failed_stage_i_search_attempt = std::any_of(
        solve.stage_i.attempt_log.begin(),
        solve.stage_i.attempt_log.end(),
        [](const auto& attempt) {
            return !attempt.solver_converged
                || !attempt.callback_error.empty()
                || !attempt.accepted;
        }
    );
    const bool all_stage_i_starts_completed =
        solve.stage_i.planned_starts.size() == 20
        && solve.stage_i.starts_completed == 20;
    const bool no_stage_i_confirmation =
        solve.stage_i.negative_confirmations == 0
        && std::none_of(
            solve.stage_i.attempt_log.begin(),
            solve.stage_i.attempt_log.end(),
            [](const auto& attempt) {
                return attempt.kind == "confirmation";
            }
        );
    const bool no_negative_tpd = std::isfinite(solve.stage_i.best_tpd)
        && solve.stage_i.best_tpd >= -1.0e-8;
    const bool completed_one_phase_search = accepted_reference
        && all_stage_i_starts_completed
        && no_stage_i_confirmation
        && no_negative_tpd;
    if (completed_one_phase_search) {
        const bool solver_passed = !failed_reference_solver_attempt
            && !failed_stage_i_solver_attempt;
        const bool numerical_passed = solver_passed
            && !failed_reference_attempt
            && !failed_stage_i_search_attempt;
        result["solver_status"] = solver_passed ? "passed" : "failed";
        result["numerical_status"] = numerical_passed ? "passed" : "failed";
        result["physical_status"] =
            numerical_passed ? "passed" : "not_adjudicated";
    } else if (
        solve.outcome == "accepted" && !solve.stage_iii_attempts.empty()
    ) {
        const auto& local = solve.stage_iii_attempts.back().result.local;
        result["solver_status"] =
            local.solver_converged ? "passed" : "failed";
        result["numerical_status"] =
            local.numerical_converged ? "passed" : "failed";
        result["physical_status"] =
            local.physical_accepted ? "passed" : "failed";
    }

    result["temperature_k"] = flash.input.temperature_k;
    result["pressure_pa"] = flash.input.pressure_pa;
    result["overall_mole_fractions"] =
        flash.input.overall_mole_fractions;
    result["parameter_fingerprint"] = flash.parameter_fingerprint;

    if (solve.outcome == "one_phase" && solve.stage_i.has_reference) {
        const auto& state = solve.stage_i.reference;
        JsonValue phase = JsonValue::object();
        phase["amount_mol"] = 1.0;
        phase["mole_fractions"] = state.amounts_mol;
        phase["volume_m3"] = state.volume_m3;
        phase["molar_density_mol_m3"] = 1.0 / state.volume_m3;
        phase["pressure_pa"] = state.provider.pressure_pa;
        phase["chemical_potential_over_rt"] = std::array<double, 2>{
            state.provider.gradient[0],
            state.provider.gradient[1],
        };
        JsonValue phases = JsonValue::array();
        phases.append(std::move(phase));
        result["phases"] = std::move(phases);
        result["phase_fractions"] = std::array<double, 1>{1.0};
        result["total_free_energy_over_rt"] = state.g_bar;
        result["pressure_stationarity_max_relative"] =
            std::abs(state.pressure_stationarity_relative);
    } else if (!solve.stage_iii_attempts.empty()) {
        const auto& refinement = solve.stage_iii_attempts.back().result;
        const auto& local = refinement.local;
        result["upper_bound"] = refinement.upper_bound;
        if (local.solver_converged) {
            result["lower_bound"] = local.evaluation.objective;
            result["held_gap"] = refinement.held_gap;
            result["material_balance_max_abs"] =
                local.material_balance_max_abs;
            result["pressure_stationarity_max_relative"] =
                local.pressure_stationarity_max_relative;
            result["kkt_stationarity_max_abs"] =
                local.kkt_stationarity_max_abs;
            result["chemical_potential_max_relative"] =
                refinement.chemical_potential_max_relative;
        }
        if (local.confirmation_solves > 0) {
            result["confirmation_succeeded"] =
                refinement.confirmation_succeeded;
            result["confirmation_max_difference"] =
                local.confirmation_max_difference;
        }
        if (solve.outcome == "accepted") {
            JsonValue phases = JsonValue::array();
            phases.append(
                held_stage_iii_phase_to_json(local.evaluation.liquid)
            );
            phases.append(
                held_stage_iii_phase_to_json(local.evaluation.vapor)
            );
            result["phases"] = std::move(phases);
            result["phase_fractions"] = std::array<double, 2>{
                local.evaluation.liquid.amounts_mol[0]
                    + local.evaluation.liquid.amounts_mol[1],
                local.evaluation.vapor.amounts_mol[0]
                    + local.evaluation.vapor.amounts_mol[1],
            };
            result["total_free_energy_over_rt"] =
                local.evaluation.objective;
            result["accepted_stage_iii"] =
                held_stage_iii_to_json(refinement);
        }
    } else if (!solve.stage_ii.trace.empty()) {
        result["lower_bound"] = solve.stage_ii.trace.back().outer_value;
        result["upper_bound"] = solve.stage_ii.upper_bound;
    }
    return result;
}

JsonValue held2_tolerance_audit_to_json(const Held2ToleranceAudit& audit) {
    JsonValue result = JsonValue::object();
    result["name"] = audit.tolerance->name;
    result["category"] = audit.tolerance->category;
    result["failure_meaning"] = audit.tolerance->failure_meaning;
    result["relation"] = held2_tolerance_relation_name(
        audit.tolerance->relation
    );
    result["residual"] = audit.residual;
    result["scale"] = audit.scale;
    result["atol"] = audit.tolerance->atol;
    result["rtol"] = audit.tolerance->rtol;
    result["threshold"] = audit.threshold;
    result["passed"] = audit.passed;
    return result;
}

JsonValue held2_pressure_envelope_to_json(
    const Held2PressureEnvelopeResult& evaluation
) {
    JsonValue points = JsonValue::array();
    for (const Held2PressureScanPoint& point : evaluation.scan_points) {
        JsonValue item = JsonValue::object();
        item["log_volume"] = point.log_volume;
        item["volume"] = point.volume;
        item["pressure_residual"] = point.pressure_residual;
        item["pressure_derivative_log_volume"] =
            point.pressure_derivative_log_volume;
        item["objective"] = point.objective;
        item["valid"] = point.valid;
        item["failure"] = point.failure;
        points.append(std::move(item));
    }
    JsonValue intervals = JsonValue::array();
    for (const Held2PressureScanInterval& interval : evaluation.intervals) {
        JsonValue item = JsonValue::object();
        item["lower_log_volume"] = interval.lower_log_volume;
        item["upper_log_volume"] = interval.upper_log_volume;
        item["depth"] = interval.depth;
        item["status"] = interval.status;
        intervals.append(std::move(item));
    }
    JsonValue roots = JsonValue::array();
    for (const Held2PressureRoot& root : evaluation.roots) {
        JsonValue item = JsonValue::object();
        item["log_volume"] = root.log_volume;
        item["volume"] = root.volume;
        item["objective"] = root.objective;
        item["pressure_residual"] = root.pressure_residual;
        item["pressure_derivative_log_volume"] =
            root.pressure_derivative_log_volume;
        item["objective_curvature_log_volume"] =
            root.objective_curvature_log_volume;
        item["mechanical_class"] = root.mechanical_class;
        item["origin"] = root.origin;
        item["boundary"] = root.boundary;
        JsonValue tolerance_audits = JsonValue::array();
        tolerance_audits.append(held2_tolerance_audit_to_json(
            audit_held2_tolerance(kHeld2RootPressure, root.pressure_residual)
        ));
        const double mechanical_margin = std::min(
            std::abs(root.pressure_derivative_log_volume),
            std::abs(root.objective_curvature_log_volume)
        );
        tolerance_audits.append(held2_tolerance_audit_to_json(
            audit_held2_tolerance(kHeld2MechanicalMargin, mechanical_margin)
        ));
        const double boundary_distance = std::min(
            std::abs(root.log_volume - evaluation.lower_log_volume),
            std::abs(evaluation.upper_log_volume - root.log_volume)
        );
        tolerance_audits.append(held2_tolerance_audit_to_json(
            audit_held2_tolerance(kHeld2RootBoundary, boundary_distance)
        ));
        item["tolerance_audits"] = std::move(tolerance_audits);
        roots.append(std::move(item));
    }
    JsonValue result = JsonValue::object();
    result["outcome"] = evaluation.outcome;
    result["failure_reason"] = evaluation.failure_reason;
    result["root_completeness"] = evaluation.root_completeness;
    result["selection_scope"] = evaluation.selection_scope;
    result["selected_root_index"] = evaluation.selected_root_index;
    result["evaluation_failure_count"] = evaluation.evaluation_failure_count;
    result["refinement_failure_count"] = evaluation.refinement_failure_count;
    result["stationary_point_count"] = evaluation.stationary_point_count;
    result["tangential_root_count"] = evaluation.tangential_root_count;
    result["marginal_root_count"] = evaluation.marginal_root_count;
    result["boundary_root_count"] = evaluation.boundary_root_count;
    result["objective_tie_count"] = evaluation.objective_tie_count;
    result["deduplicated_root_count"] = evaluation.deduplicated_root_count;
    result["lower_log_volume"] = evaluation.lower_log_volume;
    result["upper_log_volume"] = evaluation.upper_log_volume;
    result["scan_points"] = std::move(points);
    result["intervals"] = std::move(intervals);
    result["roots"] = std::move(roots);
    result["globality_certificate"] = "not_guaranteed";
    return result;
}

JsonValue held2_stage_i_reduced_evaluation_to_json(
    const Held2StageIReducedEvaluation& evaluation
) {
    JsonValue result = JsonValue::object();
    result["chart_coordinates"] = evaluation.chart_coordinates;
    result["independent_modified_fractions"] =
        evaluation.independent_modified_fractions;
    result["physical_total_ion_mole_fraction"] =
        std::isfinite(evaluation.physical_total_ion_mole_fraction)
        ? JsonValue(evaluation.physical_total_ion_mole_fraction)
        : JsonValue(nullptr);
    result["total_ion_mole_fraction_max"] =
        std::isfinite(evaluation.total_ion_mole_fraction_max)
        ? JsonValue(evaluation.total_ion_mole_fraction_max)
        : JsonValue(nullptr);
    result["tpd"] = std::isfinite(evaluation.tpd)
        ? JsonValue(evaluation.tpd)
        : JsonValue(nullptr);
    result["certified"] = evaluation.certified;
    result["failure_reason"] = evaluation.failure_reason;
    result["pressure_envelope"] =
        held2_pressure_envelope_to_json(evaluation.pressure_envelope);
    result["root_completeness"] =
        evaluation.pressure_envelope.root_completeness;
    result["tpd_tolerance_audit"] = nullptr;
    if (std::isfinite(evaluation.tpd)) {
        result["tpd_tolerance_audit"] = held2_tolerance_audit_to_json(
            audit_held2_tolerance(kHeld2TpdNegativeMargin, evaluation.tpd)
        );
    }
    result["pressure_certified"] = false;
    result["mechanical_class"] = nullptr;
    result["root_origin"] = nullptr;
    result["selected_root_log_volume"] = nullptr;
    const int selected = evaluation.pressure_envelope.selected_root_index;
    if (
        selected >= 0
        && static_cast<std::size_t>(selected)
            < evaluation.pressure_envelope.roots.size()
    ) {
        const Held2PressureRoot& root =
            evaluation.pressure_envelope.roots[
                static_cast<std::size_t>(selected)
            ];
        result["pressure_certified"] = audit_held2_tolerance(
            kHeld2RootPressure, root.pressure_residual
        ).passed;
        result["mechanical_class"] = root.mechanical_class;
        result["root_origin"] = root.origin;
        result["selected_root_log_volume"] = root.log_volume;
    }
    return result;
}

JsonValue held2_stage_i_direct_to_json(
    const Held2StageIDirectResult& evaluation
) {
    JsonValue evaluations = JsonValue::array();
    for (const Held2StageIReducedEvaluation& item : evaluation.evaluations) {
        evaluations.append(held2_stage_i_reduced_evaluation_to_json(item));
    }
    JsonValue result = JsonValue::object();
    result["outcome"] = evaluation.outcome;
    result["termination_reason"] = evaluation.termination_reason;
    result["search_strategy"] = evaluation.search_strategy;
    result["search_solver"] = evaluation.search_solver;
    result["solver_version"] = evaluation.solver_version;
    result["declared_evaluation_budget"] =
        evaluation.declared_evaluation_budget;
    result["completed_evaluation_count"] =
        evaluation.completed_evaluation_count;
    result["failed_evaluation_count"] = evaluation.failed_evaluation_count;
    result["total_ion_mole_fraction_max"] =
        std::isfinite(evaluation.total_ion_mole_fraction_max)
        ? JsonValue(evaluation.total_ion_mole_fraction_max)
        : JsonValue(nullptr);
    result["minimum_tpd"] = std::isfinite(evaluation.minimum_tpd)
        ? JsonValue(evaluation.minimum_tpd)
        : JsonValue(nullptr);
    result["evaluations"] = std::move(evaluations);
    result["negative_witness"] = nullptr;
    if (evaluation.negative_witness_index >= 0) {
        result["negative_witness"] =
            held2_stage_i_reduced_evaluation_to_json(
                evaluation.evaluations[static_cast<std::size_t>(
                    evaluation.negative_witness_index
                )]
            );
    }
    result["globality_certificate"] = evaluation.globality_certificate;
    return result;
}

JsonValue held2_stage_ii_to_json(const Held2StageIIResult& evaluation) {
    JsonValue bounds = JsonValue::array();
    for (const auto& bound : evaluation.bound_history) {
        JsonValue item = JsonValue::object();
        item["lower_bound"] = bound.lower_bound;
        item["lower_bound_available"] = bound.lower_bound_available;
        item["upper_bound"] = bound.upper_bound;
        item["multipliers"] = bound.multipliers;
        item["cut_count"] = bound.cut_count;
        item["upper_solver"] = bound.upper_solver;
        item["upper_solver_version"] = bound.upper_solver_version;
        item["upper_solver_status"] = bound.upper_solver_status;
        item["upper_primal_feasible"] = bound.upper_primal_feasible;
        item["upper_dual_feasible"] = bound.upper_dual_feasible;
        item["upper_primal_residual_inf"] =
            bound.upper_primal_residual_inf;
        item["upper_primal_scale"] = bound.upper_primal_scale;
        item["upper_dual_residual_inf"] = bound.upper_dual_residual_inf;
        item["upper_dual_scale"] = bound.upper_dual_scale;
        item["upper_complementarity_inf"] =
            bound.upper_complementarity_inf;
        JsonValue tolerance_audits = JsonValue::array();
        if (
            std::isfinite(bound.upper_primal_residual_inf)
            && std::isfinite(bound.upper_primal_scale)
        ) {
            tolerance_audits.append(held2_tolerance_audit_to_json(
                audit_held2_tolerance(
                    kHeld2LpPrimal,
                    bound.upper_primal_residual_inf,
                    bound.upper_primal_scale
                )
            ));
        }
        if (
            std::isfinite(bound.upper_dual_residual_inf)
            && std::isfinite(bound.upper_dual_scale)
        ) {
            tolerance_audits.append(held2_tolerance_audit_to_json(
                audit_held2_tolerance(
                    kHeld2LpDual,
                    bound.upper_dual_residual_inf,
                    bound.upper_dual_scale
                )
            ));
        }
        if (std::isfinite(bound.upper_complementarity_inf)) {
            tolerance_audits.append(held2_tolerance_audit_to_json(
                audit_held2_tolerance(
                    kHeld2LpComplementarity,
                    bound.upper_complementarity_inf
                )
            ));
        }
        item["tolerance_audits"] = std::move(tolerance_audits);
        item["cut_slacks"] = bound.cut_slacks;
        item["cut_duals"] = bound.cut_duals;
        item["active_cut_ids"] = bound.active_cut_ids;
        bounds.append(std::move(item));
    }

    JsonValue attempts = JsonValue::array();
    int solver_converged_count = 0;
    int physical_kkt_passed_count = 0;
    int step6_eligible_count = 0;
    for (const auto& attempt : evaluation.attempt_trace) {
        JsonValue item = JsonValue::object();
        item["attempt_id"] = attempt.attempt_id;
        item["major_iteration"] = attempt.major_iteration;
        item["start_index"] = attempt.start_index;
        item["start_source"] = attempt.start_source;
        item["internal_start"] = attempt.internal_start;
        item["physical_start_modified_fractions"] =
            attempt.physical_start_modified_fractions;
        item["physical_start_volume"] = attempt.physical_start_volume;
        item["solver_status"] = attempt.solver_status;
        item["solver_converged"] = attempt.solver_converged;
        item["provider_status"] = attempt.provider_status;
        item["callback_error"] = attempt.callback_error;
        item["internal_terminal"] = attempt.internal_terminal;
        item["terminal_modified_fractions"] =
            attempt.terminal_modified_fractions;
        item["terminal_volume"] = attempt.terminal_volume;
        item["objective"] = attempt.objective;
        item["lower_value"] = attempt.lower_value;
        item["pressure_residual"] = attempt.pressure_residual;
        item["lower_bound_multipliers"] =
            attempt.lower_bound_multipliers;
        item["upper_bound_multipliers"] =
            attempt.upper_bound_multipliers;
        item["chart_jacobian_condition"] =
            attempt.chart_jacobian_condition;
        item["dual_pullback_inf_norm"] = attempt.dual_pullback_inf_norm;
        item["dual_pullback_scale"] = attempt.dual_pullback_scale;
        item["primal_inf_norm"] = attempt.primal_inf_norm;
        item["dual_sign_violation_inf_norm"] =
            attempt.dual_sign_violation_inf_norm;
        item["chart_kkt_inf_norm"] = attempt.chart_kkt_inf_norm;
        item["physical_kkt_inf_norm"] = attempt.physical_kkt_inf_norm;
        item["complementarity_inf_norm"] =
            attempt.complementarity_inf_norm;
        item["pressure_passed"] = attempt.pressure_passed;
        item["dual_signs_valid"] = attempt.dual_signs_valid;
        item["physical_kkt_passed"] = attempt.physical_kkt_passed;
        item["cut_eligible"] = attempt.cut_eligible;
        item["step5_qualified"] = attempt.step5_qualified;
        item["step5_reason"] = attempt.step5_reason;
        item["step6_eligible"] = attempt.step6_eligible;
        item["step6_gap"] = attempt.step6_gap;
        item["step6_gap_passed"] = attempt.step6_gap_passed;
        item["step6_gradient_passed"] =
            attempt.step6_gradient_passed;
        item["step6_rejection_reason"] =
            attempt.step6_rejection_reason;
        item["fixed_volume_gradient_inf_norm"] =
            attempt.fixed_volume_gradient_inf_norm;
        item["fixed_volume_gradient_scale"] =
            attempt.fixed_volume_gradient_scale;

        JsonValue tolerance_audits = JsonValue::array();
        const bool certificate_evaluated =
            attempt.provider_status == "provider_exact"
            || attempt.provider_status == "manufactured_oracle";
        if (
            certificate_evaluated
            && std::isfinite(attempt.pressure_residual)
        ) {
            tolerance_audits.append(held2_tolerance_audit_to_json(
                audit_held2_tolerance(
                    kHeld2RootPressure, attempt.pressure_residual
                )
            ));
        }
        if (
            certificate_evaluated
            && std::isfinite(attempt.primal_inf_norm)
        ) {
            tolerance_audits.append(held2_tolerance_audit_to_json(
                audit_held2_tolerance(
                    kHeld2Stage2Primal, attempt.primal_inf_norm
                )
            ));
        }
        if (
            certificate_evaluated
            && std::isfinite(attempt.dual_sign_violation_inf_norm)
        ) {
            tolerance_audits.append(held2_tolerance_audit_to_json(
                audit_held2_tolerance(
                    kHeld2Stage2DualSign,
                    attempt.dual_sign_violation_inf_norm
                )
            ));
        }
        if (
            certificate_evaluated
            && std::isfinite(attempt.dual_pullback_inf_norm)
            && std::isfinite(attempt.dual_pullback_scale)
        ) {
            tolerance_audits.append(held2_tolerance_audit_to_json(
                audit_held2_tolerance(
                    kHeld2Stage2DualPullback,
                    attempt.dual_pullback_inf_norm,
                    attempt.dual_pullback_scale
                )
            ));
        }
        if (
            certificate_evaluated
            && std::isfinite(attempt.physical_kkt_inf_norm)
        ) {
            tolerance_audits.append(held2_tolerance_audit_to_json(
                audit_held2_tolerance(
                    kHeld2Stage2Stationarity,
                    attempt.physical_kkt_inf_norm
                )
            ));
        }
        if (
            certificate_evaluated
            && std::isfinite(attempt.complementarity_inf_norm)
        ) {
            tolerance_audits.append(held2_tolerance_audit_to_json(
                audit_held2_tolerance(
                    kHeld2Stage2Complementarity,
                    attempt.complementarity_inf_norm
                )
            ));
        }
        if (certificate_evaluated && std::isfinite(attempt.step6_gap)) {
            tolerance_audits.append(held2_tolerance_audit_to_json(
                audit_held2_tolerance(kHeld2Step6Gap, attempt.step6_gap)
            ));
        }
        if (
            certificate_evaluated
            && std::isfinite(attempt.fixed_volume_gradient_inf_norm)
            && std::isfinite(attempt.fixed_volume_gradient_scale)
        ) {
            tolerance_audits.append(held2_tolerance_audit_to_json(
                audit_held2_tolerance(
                    kHeld2Step6Gradient,
                    attempt.fixed_volume_gradient_inf_norm,
                    attempt.fixed_volume_gradient_scale
                )
            ));
        }
        item["tolerance_audits"] = std::move(tolerance_audits);
        item["basin_id"] = attempt.basin_id;
        item["same_major_upper_bound"] =
            attempt.same_major_upper_bound;
        item["same_major_multipliers"] =
            attempt.same_major_multipliers;
        attempts.append(std::move(item));
        solver_converged_count += attempt.solver_converged ? 1 : 0;
        physical_kkt_passed_count +=
            attempt.physical_kkt_passed ? 1 : 0;
        step6_eligible_count += attempt.step6_eligible ? 1 : 0;
    }

    JsonValue candidates = JsonValue::array();
    for (const auto& candidate : evaluation.candidates) {
        JsonValue item = JsonValue::object();
        item["modified_fractions"] = candidate.modified_fractions;
        item["independent_modified_fractions"] =
            candidate.independent_modified_fractions;
        item["volume"] = candidate.volume;
        item["phase_coordinate"] = candidate.phase_coordinate;
        item["lower_gap"] = candidate.lower_gap;
        candidates.append(std::move(item));
    }

    JsonValue major_contexts = JsonValue::array();
    for (const auto& context : evaluation.major_contexts) {
        JsonValue item = JsonValue::object();
        item["major_id"] = context.major_id;
        item["upper_solve_id"] = context.upper_solve_id;
        item["upper_bound"] = context.upper_bound;
        item["multipliers"] = context.multipliers;
        item["active_cut_ids"] = context.active_cut_ids;
        item["lower_attempt_ids"] = context.lower_attempt_ids;
        item["current_basin_ids"] = context.current_basin_ids;
        item["pressure_branch_ids"] = context.pressure_branch_ids;
        item["step5_qualified_attempt_ids"] =
            context.step5_qualified_attempt_ids;
        item["step6_eligible_attempt_ids"] =
            context.step6_eligible_attempt_ids;
        item["certificate_failed_attempt_ids"] =
            context.certificate_failed_attempt_ids;
        major_contexts.append(std::move(item));
    }

    JsonValue result = JsonValue::object();
    result["outcome"] = evaluation.outcome;
    result["search_strategy"] = evaluation.search_strategy;
    result["global_explorer"] = evaluation.global_explorer;
    result["local_solver"] = evaluation.local_solver;
    result["globality_certificate"] = evaluation.globality_certificate;
    result["major_iterations"] = evaluation.major_iterations;
    result["lower_starts_per_iteration"] =
        evaluation.lower_starts_per_iteration;
    result["cut_count"] = evaluation.cut_count;
    result["exploration_evaluation_count"] =
        evaluation.exploration_evaluation_count;
    result["exploration_failure_count"] =
        evaluation.exploration_failure_count;
    result["exploration_representative_count"] =
        evaluation.exploration_representative_count;
    result["duplicate_representative_count"] =
        evaluation.duplicate_representative_count;
    result["duplicate_terminal_count"] =
        evaluation.duplicate_terminal_count;
    result["distinct_basin_count"] = evaluation.distinct_basin_count;
    result["unresolved_candidate_identity_count"] =
        evaluation.unresolved_candidate_identity_count;
    result["local_attempt_cap_per_major"] =
        evaluation.local_attempt_cap_per_major;
    result["local_attempts_truncated"] =
        evaluation.local_attempts_truncated;
    result["direct_escalation_used"] = evaluation.direct_escalation_used;
    result["next_action"] = held2_next_action_name(evaluation.next_action);
    result["bound_history"] = std::move(bounds);
    result["major_contexts"] = std::move(major_contexts);
    result["attempt_trace"] = std::move(attempts);
    JsonValue attempt_classification = JsonValue::object();
    attempt_classification["declared"] = evaluation.attempt_trace.size();
    attempt_classification["solver_converged"] =
        solver_converged_count;
    attempt_classification["solver_failed"] =
        static_cast<int>(evaluation.attempt_trace.size())
        - solver_converged_count;
    attempt_classification["physical_kkt_passed"] =
        physical_kkt_passed_count;
    attempt_classification["step6_eligible"] = step6_eligible_count;
    result["attempt_classification"] =
        std::move(attempt_classification);
    result["candidates"] = std::move(candidates);
    return result;
}

JsonValue held2_stage_iii_to_json(
    const Held2StageIIIResult& evaluation
) {
    JsonValue phases = JsonValue::array();
    for (const auto& phase : evaluation.phases) {
        JsonValue item = JsonValue::object();
        item["phase_fraction"] = phase.phase_fraction;
        item["modified_fractions"] = phase.modified_fractions;
        item["physical_fractions"] = phase.physical_fractions;
        item["volume"] = phase.volume;
        phases.append(std::move(item));
    }
    JsonValue lifecycle = JsonValue::array();
    for (const auto& step : evaluation.lifecycle) {
        JsonValue item = JsonValue::object();
        item["solve_index"] = step.solve_index;
        item["active_candidate_count"] = step.active_candidate_count;
        item["removed_candidate_index"] = step.removed_candidate_index;
        item["action"] = step.action;
        item["phase_fraction"] = step.phase_fraction;
        item["lower_bound_multiplier"] = step.lower_bound_multiplier;
        item["reduced_derivative"] = step.reduced_derivative;
        item["complementarity_inf_norm"] =
            step.complementarity_inf_norm;
        item["candidate_independent_modified_fractions"] =
            step.candidate_independent_modified_fractions;
        item["candidate_volume"] = step.candidate_volume;
        item["solver_status"] = step.solver_status;
        item["decision_reason"] = step.decision_reason;
        lifecycle.append(std::move(item));
    }

    JsonValue result = JsonValue::object();
    result["solver_status"] = evaluation.solver_status;
    result["numerical_status"] = evaluation.numerical_status;
    result["physical_status"] = evaluation.physical_status;
    result["feedback"] = evaluation.feedback;
    result["failure_reason"] = evaluation.failure_reason;
    result["trace_refinement_status"] =
        evaluation.trace_refinement_status;
    result["input_candidate_count"] = evaluation.input_candidate_count;
    result["retired_duplicate_count"] =
        evaluation.retired_duplicate_count;
    result["retired_inactive_count"] = evaluation.retired_inactive_count;
    result["stage_iii_solve_count"] = evaluation.stage_iii_solve_count;
    result["active_set_resolve_count"] =
        evaluation.active_set_resolve_count;
    result["pressure_polish_iteration_count"] =
        evaluation.pressure_polish_iteration_count;
    result["pressure_polish_status"] = evaluation.pressure_polish_status;
    result["trace_component_count"] = evaluation.trace_component_count;
    result["certified_modified_potential_count"] =
        evaluation.certified_modified_potential_count;
    result["objective"] = evaluation.objective;
    result["modified_balance_inf_norm"] =
        evaluation.modified_balance_inf_norm;
    result["ordinary_balance_inf_norm"] =
        evaluation.ordinary_balance_inf_norm;
    result["phase_charge_inf_norm"] = evaluation.phase_charge_inf_norm;
    result["phase_charge_scale"] = evaluation.phase_charge_scale;
    result["pressure_stationarity_inf_norm"] =
        evaluation.pressure_stationarity_inf_norm;
    result["modified_potential_mixed_gap"] =
        evaluation.modified_potential_mixed_gap;
    result["modified_potential_scale"] =
        evaluation.modified_potential_scale;
    result["minimum_phase_distance"] = evaluation.minimum_phase_distance;
    result["phase_identity_status"] = evaluation.phase_identity_status;
    result["kkt_stationarity_inf_norm"] =
        evaluation.kkt_stationarity_inf_norm;
    result["dual_sign_violation_inf_norm"] =
        evaluation.dual_sign_violation_inf_norm;
    result["bound_complementarity_inf_norm"] =
        evaluation.bound_complementarity_inf_norm;
    result["minimum_phase_fraction"] = evaluation.minimum_phase_fraction;
    result["free_energy_upper_bound"] =
        evaluation.free_energy_upper_bound;
    result["free_energy_gap"] = evaluation.free_energy_gap;
    result["free_energy_gap_provenance"] =
        evaluation.free_energy_gap_provenance;
    result["kkt_evidence_available"] =
        evaluation.kkt_evidence_available;
    result["physical_evidence_available"] =
        evaluation.physical_evidence_available;
    result["phase_identity_evidence_available"] =
        evaluation.phase_identity_evidence_available;
    result["free_energy_gap_available"] =
        evaluation.free_energy_gap_available;

    JsonValue tolerance_audits = JsonValue::array();
    if (evaluation.physical_evidence_available) {
        tolerance_audits.append(held2_tolerance_audit_to_json(
            audit_held2_tolerance(
                kHeld2Stage3ModifiedBalance,
                evaluation.modified_balance_inf_norm
            )
        ));
        tolerance_audits.append(held2_tolerance_audit_to_json(
            audit_held2_tolerance(
                kHeld2Stage3ExplicitBalance,
                evaluation.ordinary_balance_inf_norm
            )
        ));
        tolerance_audits.append(held2_tolerance_audit_to_json(
            audit_held2_tolerance(
                kHeld2Stage3Charge,
                evaluation.phase_charge_inf_norm,
                evaluation.phase_charge_scale
            )
        ));
        tolerance_audits.append(held2_tolerance_audit_to_json(
            audit_held2_tolerance(
                kHeld2Stage3Pressure,
                evaluation.pressure_stationarity_inf_norm
            )
        ));
        tolerance_audits.append(held2_tolerance_audit_to_json(
            audit_held2_tolerance(
                kHeld2Stage3Potential,
                evaluation.modified_potential_mixed_gap,
                evaluation.modified_potential_scale
            )
        ));
    }
    if (evaluation.kkt_evidence_available) {
        tolerance_audits.append(held2_tolerance_audit_to_json(
            audit_held2_tolerance(
                kHeld2Stage3Stationarity,
                evaluation.kkt_stationarity_inf_norm
            )
        ));
        tolerance_audits.append(held2_tolerance_audit_to_json(
            audit_held2_tolerance(
                kHeld2Stage3DualSign,
                evaluation.dual_sign_violation_inf_norm
            )
        ));
        tolerance_audits.append(held2_tolerance_audit_to_json(
            audit_held2_tolerance(
                kHeld2Stage3Complementarity,
                evaluation.bound_complementarity_inf_norm
            )
        ));
        tolerance_audits.append(held2_tolerance_audit_to_json(
            audit_held2_tolerance(
                kHeld2PhaseActivity,
                evaluation.minimum_phase_fraction
            )
        ));
    }
    if (evaluation.phase_identity_evidence_available) {
        tolerance_audits.append(held2_tolerance_audit_to_json(
            audit_held2_tolerance(
                kHeld2PhaseDistinct,
                evaluation.minimum_phase_distance
            )
        ));
    }
    if (evaluation.free_energy_gap_available) {
        tolerance_audits.append(held2_tolerance_audit_to_json(
            audit_held2_tolerance(
                kHeld2Stage3FreeEnergyGap,
                evaluation.free_energy_gap
            )
        ));
    }
    result["tolerance_audits"] = std::move(tolerance_audits);
    result["phases"] = std::move(phases);
    result["lifecycle"] = std::move(lifecycle);
    result["globality_certificate"] = "not_guaranteed";
    return result;
}

void update_held2_workflow_json(
    JsonValue& result,
    const Held2WorkflowState& workflow
) {
    JsonValue transitions = JsonValue::array();
    for (const auto& transition : workflow.transitions) {
        JsonValue item = JsonValue::object();
        item["completed_step"] = transition.completed_step;
        item["next_action"] =
            held2_next_action_name(transition.next_action);
        item["outcome"] = transition.outcome;
        item["reason"] = transition.reason;
        item["feedback_owner_step"] = transition.feedback_owner_step;
        transitions.append(std::move(item));
    }
    result["outcome"] = workflow.outcome;
    result["failure_stage"] = workflow.failure_stage.empty()
        ? JsonValue(nullptr)
        : JsonValue(workflow.failure_stage);
    result["failure_reason"] = workflow.failure_reason.empty()
        ? JsonValue(nullptr)
        : JsonValue(workflow.failure_reason);
    result["completed_step"] = workflow.completed_step;
    result["feedback_owner_step"] = workflow.feedback_owner_step;
    result["next_action"] = held2_next_action_name(workflow.next_action);
    result["transitions"] = std::move(transitions);
}

JsonValue held2_flash_to_json(const FlashResult& flash) {
    const Held2FlashResult& solve = flash.held2;
    JsonValue result = JsonValue::object();
    result["controller"] = "perdomo_held2_steps_1_to_10_v1";
    result["stage_order"] = std::array<const char*, 4>{
        "step_1_reference",
        "steps_2_3_stage_i",
        "steps_4_7_stage_ii",
        "steps_8_10_stage_iii",
    };
    result["parameter_fingerprint"] = flash.parameter_fingerprint;
    result["globality_certificate"] = flash.globality_certificate;
    result["reference_pressure_envelope"] =
        held2_pressure_envelope_to_json(
            solve.reference_pressure_envelope
        );
    result["stage_i"] = held2_stage_i_direct_to_json(solve.stage_i);
    result["stage_ii"] = solve.stage_ii.has_value()
        ? held2_stage_ii_to_json(*solve.stage_ii)
        : JsonValue(nullptr);
    result["stage_iii"] = solve.stage_iii.has_value()
        ? held2_stage_iii_to_json(*solve.stage_iii)
        : JsonValue(nullptr);
    result["stage_ii_skip_reason"] = solve.stage_ii_skip_reason.empty()
        ? JsonValue(nullptr)
        : JsonValue(solve.stage_ii_skip_reason);
    result["stage_iii_skip_reason"] =
        solve.stage_iii_skip_reason.empty()
        ? JsonValue(nullptr)
        : JsonValue(solve.stage_iii_skip_reason);
    result["predictive_comparison_status"] =
        solve.predictive_comparison_status;
    update_held2_workflow_json(result, solve.workflow);
    return result;
}

JsonValue paper_timing_to_json(const Held2StepTiming& timing) {
    JsonValue result = JsonValue::object();
    result["step"] = timing.step;
    result["invocation_count"] = timing.invocation_count;
    result["wall_seconds"] = timing.wall_seconds;
    result["cpu_seconds"] = timing.cpu_seconds;
    result["provider_evaluations"] = timing.provider_evaluations;
    result["optimizer_solves"] = timing.optimizer_solves;
    result["optimizer_iterations"] = timing.optimizer_iterations;
    result["terminal_status"] = timing.terminal_status;
    result["terminal_reason"] = timing.terminal_reason;
    result["next_step"] = timing.next_step;
    return result;
}

JsonValue paper_phase_to_json(const Held2Phase& phase) {
    JsonValue result = JsonValue::object();
    result["stable_id"] = phase.stable_id;
    result["phase_fraction"] = phase.phase_fraction;
    result["independent_modified_fractions"] =
        phase.independent_modified_fractions;
    result["mole_fractions"] =
        phase.physical_fractions_provider_order;
    result["molar_volume_m3_mol"] = phase.volume;
    result["packing_fraction"] = phase.packing_fraction;
    return result;
}

JsonValue paper_physical_certificate_to_json(
    const Held2PhysicalCertificate& value
) {
    JsonValue result = JsonValue::object();
    result["modified_balance_inf"] = value.modified_balance_inf;
    result["ordinary_balance_inf"] = value.ordinary_balance_inf;
    result["electroneutrality_inf"] = value.electroneutrality_inf;
    result["pressure_residual_inf"] = value.pressure_residual_inf;
    result["kkt_residual_inf"] = value.kkt_residual_inf;
    result["accepted"] = value.accepted;
    return result;
}

JsonValue paper_step8_to_json(const Held2Step8Result& step) {
    JsonValue result = JsonValue::object();
    result["outcome"] =
        step.outcome == Held2Step8Outcome::CertifiedFeasible
        ? "certified_feasible"
        : step.outcome == Held2Step8Outcome::CertifiedInfeasible
            ? "certified_infeasible"
            : step.outcome == Held2Step8Outcome::InsufficientCandidates
                ? "insufficient_candidates" : "indeterminate";
    result["reason"] = step.reason;
    result["total_reduced_gibbs"] = step.total_reduced_gibbs
        ? JsonValue(*step.total_reduced_gibbs) : JsonValue(nullptr);
    result["ordinary_balance_inf"] = step.ordinary_balance_inf;
    result["electroneutrality_inf"] = step.electroneutrality_inf;
    result["electroneutrality_scale"] = step.electroneutrality_scale;
    result["pressure_residual_inf"] = step.pressure_residual_inf;
    JsonValue phases = JsonValue::array();
    for (const Held2Phase& phase : step.active_phases) {
        phases.append(paper_phase_to_json(phase));
    }
    result["active_phases"] = std::move(phases);
    if (step.feasibility) {
        JsonValue certificate = JsonValue::object();
        certificate["solver_status"] = step.feasibility->solver_status;
        certificate["feasible"] = step.feasibility->feasible;
        certificate["infeasible"] = step.feasibility->infeasible;
        certificate["farkas_certificate_valid"] =
            step.feasibility->farkas_certificate_valid;
        certificate["primal_residual_inf"] =
            step.feasibility->primal_residual_inf;
        certificate["certificate_residual_inf"] =
            step.feasibility->certificate_residual_inf;
        result["feasibility"] = std::move(certificate);
    } else {
        result["feasibility"] = nullptr;
    }
    if (step.nlp) {
        JsonValue certificate = JsonValue::object();
        certificate["solver_status"] = step.nlp->solver_status;
        certificate["primal_residual_inf"] = step.nlp->primal_residual_inf;
        certificate["stationarity_residual_inf"] =
            step.nlp->stationarity_residual_inf;
        certificate["dual_sign_violation_inf"] =
            step.nlp->dual_sign_violation_inf;
        certificate["complementarity_inf"] =
            step.nlp->complementarity_inf;
        certificate["accepted"] = step.nlp->accepted;
        result["nlp"] = std::move(certificate);
    } else {
        result["nlp"] = nullptr;
    }
    JsonValue lifecycle = JsonValue::array();
    for (const Held2LifecycleDecision& value : step.lifecycle) {
        JsonValue decision = JsonValue::object();
        decision["stable_id"] = value.stable_id;
        decision["action"] = value.action;
        decision["reason"] = value.reason;
        decision["reduced_resolve_accepted"] =
            value.reduced_resolve_accepted;
        lifecycle.append(std::move(decision));
    }
    result["lifecycle"] = std::move(lifecycle);
    result["timing"] = paper_timing_to_json(step.timing);
    return result;
}

JsonValue paper_step9_to_json(const Held2Step9Result& step) {
    JsonValue result = JsonValue::object();
    result["outcome"] = step.outcome == Held2Step9Outcome::Converged
        ? "converged"
        : step.outcome == Held2Step9Outcome::PaperConvergenceFailed
            ? "paper_convergence_failed" : "indeterminate";
    result["reason"] = step.reason;
    result["free_energy_gap"] = step.free_energy_gap
        ? JsonValue(*step.free_energy_gap) : JsonValue(nullptr);
    result["physical"] = step.physical
        ? paper_physical_certificate_to_json(*step.physical)
        : JsonValue(nullptr);
    JsonValue comparisons = JsonValue::array();
    for (const Held2PotentialComparison& value : step.potential_comparisons) {
        JsonValue comparison = JsonValue::object();
        comparison["component_index"] = value.component_index;
        comparison["left_phase_id"] = value.left_phase_id;
        comparison["right_phase_id"] = value.right_phase_id;
        comparison["numerator"] = value.numerator;
        comparison["denominator"] = value.denominator;
        comparison["ratio"] = value.ratio;
        comparison["passed"] = value.passed;
        comparisons.append(std::move(comparison));
    }
    result["potential_comparisons"] = std::move(comparisons);
    result["timing"] = paper_timing_to_json(step.timing);
    return result;
}

JsonValue paper_step10_to_json(const Held2Step10Result& step) {
    JsonValue result = JsonValue::object();
    result["status"] = step.status;
    result["reason"] = step.reason;
    JsonValue phases = JsonValue::array();
    for (const Held2Phase& phase : step.phases) {
        phases.append(paper_phase_to_json(phase));
    }
    result["phases"] = std::move(phases);
    JsonValue refinements = JsonValue::array();
    for (const Held2TraceRefinement& value : step.refinements) {
        JsonValue refinement = JsonValue::object();
        refinement["phase_id"] = value.phase_id;
        refinement["component_index"] = value.component_index;
        refinement["reference_phase_id"] = value.reference_phase_id;
        refinement["initial_mole_fraction"] = value.initial_mole_fraction;
        refinement["refined_mole_fraction"] = value.refined_mole_fraction;
        refinement["final_potential_residual"] =
            value.final_potential_residual;
        refinement["status"] = value.status;
        refinements.append(std::move(refinement));
    }
    result["refinements"] = std::move(refinements);
    result["final_certificate"] = step.final_certificate
        ? paper_physical_certificate_to_json(*step.final_certificate)
        : JsonValue(nullptr);
    result["timing"] = paper_timing_to_json(step.timing);
    return result;
}

JsonValue paper_algorithm_to_json(
    const Held2AlgorithmResult& solve,
    const Held2Input& input,
    const std::string& fingerprint
) {
    JsonValue result = JsonValue::object();
    result["controller"] = "perdomo_held2_paper_steps_1_to_10_v1";
    result["route"] = "held2";
    result["parameter_fingerprint"] = fingerprint;
    result["temperature_k"] = input.temperature_k;
    result["pressure_pa"] = input.pressure_pa;
    result["overall_mole_fractions"] =
        input.overall_mole_fractions_provider_order;
    result["outcome"] = solve.outcome;
    result["failure_stage"] = solve.failure_stage.empty()
        ? JsonValue(nullptr) : JsonValue(solve.failure_stage);
    result["failure_reason"] = solve.failure_reason.empty()
        ? JsonValue(nullptr) : JsonValue(solve.failure_reason);
    result["globality_certificate"] = solve.globality_certificate;
    result["phase_enumeration_certificate"] =
        solve.phase_enumeration_certificate;
    result["upper_solve_count"] = solve.upper_solve_count;
    JsonValue tolerances = JsonValue::object();
    tolerances["epsilon_b"] = kHeld2PaperStep6Gap.atol;
    tolerances["epsilon_lambda"] = kHeld2PaperStep6Derivative.rtol;
    tolerances["epsilon_lambda_floor"] =
        kHeld2PaperStep6Derivative.atol;
    tolerances["epsilon_eta"] = kHeld2PaperStep6PackingDistinct.atol;
    tolerances["epsilon_x"] = kHeld2PaperStep6CompositionDistinct.atol;
    tolerances["epsilon_g"] = kHeld2PaperFreeEnergyGap.atol;
    tolerances["epsilon_mu"] = kHeld2PaperPotentialRatio.atol;
    result["paper_default_tolerances"] = std::move(tolerances);

    JsonValue phases = JsonValue::array();
    for (const Held2Phase& phase : solve.phases) {
        phases.append(paper_phase_to_json(phase));
    }
    result["phases"] = std::move(phases);
    JsonValue timings = JsonValue::array();
    for (const Held2StepTiming& timing : solve.step_timings) {
        timings.append(paper_timing_to_json(timing));
    }
    result["step_timings"] = std::move(timings);

    JsonValue step1 = JsonValue::object();
    step1["status"] = solve.step1.status;
    step1["reason"] = solve.step1.reason;
    step1["independent_feed"] = solve.step1.independent_feed
        ? JsonValue(*solve.step1.independent_feed) : JsonValue(nullptr);
    step1["timing"] = paper_timing_to_json(solve.step1.timing);
    result["step1"] = std::move(step1);

    if (solve.step2) {
        JsonValue step2 = JsonValue::object();
        step2["reason"] = solve.step2->reason;
        step2["globality_certificate"] =
            solve.step2->globality_certificate;
        step2["minimum_tpd"] = solve.step2->minimum_tpd
            ? JsonValue(*solve.step2->minimum_tpd) : JsonValue(nullptr);
        step2["reference_pressure_envelope"] =
            solve.step2->reference_envelope
            ? held2_pressure_envelope_to_json(
                *solve.step2->reference_envelope
            ) : JsonValue(nullptr);
        step2["timing"] = paper_timing_to_json(solve.step2->timing);
        result["step2"] = std::move(step2);
    } else {
        result["step2"] = nullptr;
    }
    const auto history = [](const auto& values) {
        JsonValue items = JsonValue::array();
        for (const auto& value : values) {
            JsonValue item = JsonValue::object();
            item["status"] = value.status;
            item["reason"] = value.reason;
            item["timing"] = paper_timing_to_json(value.timing);
            items.append(std::move(item));
        }
        return items;
    };
    result["step4_history"] = history(solve.step4_history);
    result["step5_history"] = history(solve.step5_history);
    result["step6_history"] = history(solve.step6_history);
    result["step7_history"] = history(solve.step7_history);

    JsonValue step8_history = JsonValue::array();
    for (const Held2Step8Result& value : solve.step8_history) {
        step8_history.append(paper_step8_to_json(value));
    }
    result["step8_history"] = std::move(step8_history);
    JsonValue step9_history = JsonValue::array();
    for (const Held2Step9Result& value : solve.step9_history) {
        step9_history.append(paper_step9_to_json(value));
    }
    result["step9_history"] = std::move(step9_history);
    if (solve.step10) {
        result["step10"] = paper_step10_to_json(*solve.step10);
    } else {
        result["step10"] = nullptr;
    }
    return result;
}

}  // namespace

std::string flash_result_to_json(const FlashResult& result) {
    std::ostringstream output;
    write_json(
        output,
        result.route == FlashResult::Route::Held2
            ? held2_flash_to_json(result)
            : neutral_flash_to_json(result)
    );
    output << '\n';
    return output.str();
}

std::string held2_algorithm_result_to_json(
    const Held2AlgorithmResult& result,
    const Held2Input& input,
    const std::string& parameter_fingerprint
) {
    std::ostringstream output;
    write_json(
        output,
        paper_algorithm_to_json(result, input, parameter_fingerprint)
    );
    output << '\n';
    return output.str();
}

}  // namespace epcsaft_equilibrium
