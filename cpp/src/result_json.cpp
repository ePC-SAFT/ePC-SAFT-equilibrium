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
    const HeldResult& solve = std::get<HeldResult>(flash.solve);
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

const char* tolerance_relation_name(Held2ToleranceRelation relation) {
    switch (relation) {
        case Held2ToleranceRelation::AbsAtMost: return "abs_at_most";
        case Held2ToleranceRelation::AtLeast: return "at_least";
        case Held2ToleranceRelation::GreaterThan: return "greater_than";
        case Held2ToleranceRelation::LessThanNegative:
            return "less_than_negative";
        case Held2ToleranceRelation::SolverTarget: return "solver_target";
    }
    return "unknown";
}

JsonValue held2_tolerance_audit_to_json(const Held2ToleranceAudit& audit) {
    JsonValue result = JsonValue::object();
    result["name"] = audit.tolerance->name;
    result["category"] = audit.tolerance->category;
    result["failure_meaning"] = audit.tolerance->failure_meaning;
    result["relation"] = tolerance_relation_name(
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
    result["helmholtz_over_rt_reference_amount"] =
        phase.helmholtz_over_rt_reference_amount;
    result["pressure_pa"] = phase.pressure_pa;
    result["chemical_potential_over_rt"] =
        phase.chemical_potentials_over_rt;
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
    result["candidate_ids"] = step.candidate_ids;
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
    result["total_free_energy_over_rt"] =
        solve.total_free_energy_over_rt
        ? JsonValue(*solve.total_free_energy_over_rt)
        : JsonValue(nullptr);
    JsonValue tolerances = JsonValue::object();
    tolerances["epsilon_b"] = kHeld2PaperStep6Gap.atol;
    tolerances["epsilon_lambda"] = kHeld2PaperStep6Derivative.rtol;
    tolerances["epsilon_lambda_floor"] =
        kHeld2PaperStep6Derivative.atol;
    tolerances["epsilon_eta"] = kHeld2PaperStep6PackingDistinct.atol;
    tolerances["epsilon_x"] = kHeld2PaperStep6CompositionDistinct.atol;
    tolerances["epsilon_g"] = kHeld2PaperFreeEnergyGap.atol;
    tolerances["epsilon_mu"] = kHeld2PaperPotentialRatio.atol;
    result["effective_tolerances"] = std::move(tolerances);

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
    if (const auto* held2 = std::get_if<Held2AlgorithmResult>(&result.solve)) {
        write_json(
            output,
            paper_algorithm_to_json(
                *held2,
                {
                    result.input.temperature_k,
                    result.input.pressure_pa,
                    result.input.overall_mole_fractions,
                },
                result.parameter_fingerprint
            )
        );
    } else {
        write_json(output, neutral_flash_to_json(result));
    }
    output << '\n';
    return output.str();
}

}  // namespace epcsaft_equilibrium
