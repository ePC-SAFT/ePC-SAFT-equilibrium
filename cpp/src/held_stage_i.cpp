#include "held.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>
#include <coin/IpTNLP.hpp>

namespace epcsaft_equilibrium {
namespace {

constexpr double kGasConstantJPerMolK = 8.31446261815324;
constexpr double kCompositionLower = 1.0e-8;
constexpr double kCompositionUpper = 1.0 - 1.0e-8;
constexpr double kVolumeLowerM3PerMol = 1.0e-5;
constexpr double kVolumeUpperM3PerMol = 1.0e-1;
constexpr double kLocalPhysicalTolerance = 1.0e-8;
constexpr double kConfirmationTolerance = 1.0e-7;
constexpr double kTpdNegativeThreshold = -1.0e-8;
constexpr double kBoundMargin = 1.0e-7;
constexpr int kIpoptIterations = 300;

std::string ipopt_status_name(Ipopt::ApplicationReturnStatus status) {
    switch (status) {
        case Ipopt::Solve_Succeeded:
            return "solve_succeeded";
        case Ipopt::Solved_To_Acceptable_Level:
            return "solved_to_acceptable_level";
        case Ipopt::Infeasible_Problem_Detected:
            return "infeasible_problem_detected";
        case Ipopt::Search_Direction_Becomes_Too_Small:
            return "search_direction_too_small";
        case Ipopt::Diverging_Iterates:
            return "diverging_iterates";
        case Ipopt::Feasible_Point_Found:
            return "feasible_point_found";
        case Ipopt::Maximum_Iterations_Exceeded:
            return "maximum_iterations_exceeded";
        case Ipopt::Restoration_Failed:
            return "restoration_failed";
        case Ipopt::Error_In_Step_Computation:
            return "error_in_step_computation";
        case Ipopt::Not_Enough_Degrees_Of_Freedom:
            return "not_enough_degrees_of_freedom";
        case Ipopt::Invalid_Problem_Definition:
            return "invalid_problem_definition";
        case Ipopt::Invalid_Option:
            return "invalid_option";
        case Ipopt::Invalid_Number_Detected:
            return "invalid_number_detected";
        case Ipopt::Unrecoverable_Exception:
            return "unrecoverable_exception";
        case Ipopt::NonIpopt_Exception_Thrown:
            return "non_ipopt_exception_thrown";
        case Ipopt::Insufficient_Memory:
            return "insufficient_memory";
        case Ipopt::Internal_Error:
            return "internal_error";
        default:
            return "ipopt_status_" + std::to_string(static_cast<int>(status));
    }
}

void configure_ipopt(const Ipopt::SmartPtr<Ipopt::IpoptApplication>& application) {
    application->Options()->SetStringValue("option_file_name", "");
    application->Options()->SetIntegerValue("print_level", 0);
    application->Options()->SetStringValue("sb", "yes");
    application->Options()->SetIntegerValue("max_iter", kIpoptIterations);
    application->Options()->SetNumericValue("tol", 1.0e-10);
    application->Options()->SetNumericValue("acceptable_tol", 1.0e-9);
    application->Options()->SetIntegerValue("acceptable_iter", 0);
    application->Options()->SetNumericValue("constr_viol_tol", 1.0e-10);
    application->Options()->SetStringValue("jacobian_approximation", "exact");
    application->Options()->SetStringValue("hessian_approximation", "exact");
    application->Options()->SetStringValue("nlp_scaling_method", "none");
    application->Options()->SetNumericValue("bound_relax_factor", 0.0);
    application->Options()->SetStringValue("honor_original_bounds", "yes");
    application->Options()->SetStringValue("check_derivatives_for_naninf", "yes");
}

enum class StageIObjective { homogeneous, tpd, tunneling };

struct SearchRun {
    bool solver_converged = false;
    std::string solver_status;
    std::string callback_error;
    int iterations = 0;
    std::array<double, 2> initial_guess{};
    std::array<double, 2> variables{};
};

class SearchTnlp final : public Ipopt::TNLP {
public:
    SearchTnlp(
        const ProviderContext& provider,
        double temperature_k,
        double pressure_pa,
        const HeldStateEvaluation& reference,
        double minimum_x_methane,
        double minimum_tpd,
        StageIObjective objective,
        std::array<double, 2> initial,
        std::array<double, 2> lower,
        std::array<double, 2> upper
    )
        : provider_(provider),
          temperature_k_(temperature_k),
          pressure_pa_(pressure_pa),
          reference_(reference),
          minimum_x_methane_(minimum_x_methane),
          minimum_tpd_(minimum_tpd),
          objective_(objective),
          initial_(initial),
          lower_(lower),
          upper_(upper) {}

    bool get_nlp_info(
        Ipopt::Index& n,
        Ipopt::Index& m,
        Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag,
        IndexStyleEnum& index_style
    ) override {
        n = 2;
        m = 0;
        nnz_jac_g = 0;
        nnz_h_lag = 3;
        index_style = TNLP::C_STYLE;
        return true;
    }

    bool get_bounds_info(
        Ipopt::Index n,
        Ipopt::Number* x_l,
        Ipopt::Number* x_u,
        Ipopt::Index m,
        Ipopt::Number*,
        Ipopt::Number*
    ) override {
        if (n != 2 || m != 0) {
            return false;
        }
        std::copy(lower_.begin(), lower_.end(), x_l);
        std::copy(upper_.begin(), upper_.end(), x_u);
        return true;
    }

    bool get_starting_point(
        Ipopt::Index n,
        bool init_x,
        Ipopt::Number* x,
        bool init_z,
        Ipopt::Number*,
        Ipopt::Number*,
        Ipopt::Index m,
        bool init_lambda,
        Ipopt::Number*
    ) override {
        if (n != 2 || m != 0 || !init_x || init_z || init_lambda) {
            return false;
        }
        std::copy(initial_.begin(), initial_.end(), x);
        return true;
    }

    bool eval_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number& objective
    ) override {
        if (n != 2) {
            return false;
        }
        try {
            objective = evaluate(x).objective;
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
    }

    bool eval_grad_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number* gradient
    ) override {
        if (n != 2) {
            return false;
        }
        try {
            const ObjectiveEvaluation evaluation = evaluate(x);
            std::copy(evaluation.gradient.begin(), evaluation.gradient.end(), gradient);
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
    }

    bool eval_g(
        Ipopt::Index,
        const Ipopt::Number*,
        bool,
        Ipopt::Index m,
        Ipopt::Number*
    ) override {
        return m == 0;
    }

    bool eval_jac_g(
        Ipopt::Index,
        const Ipopt::Number*,
        bool,
        Ipopt::Index m,
        Ipopt::Index nonzero_count,
        Ipopt::Index*,
        Ipopt::Index*,
        Ipopt::Number*
    ) override {
        return m == 0 && nonzero_count == 0;
    }

    bool eval_h(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number objective_factor,
        Ipopt::Index m,
        const Ipopt::Number*,
        bool,
        Ipopt::Index nonzero_count,
        Ipopt::Index* rows,
        Ipopt::Index* columns,
        Ipopt::Number* values
    ) override {
        if (n != 2 || m != 0 || nonzero_count != 3) {
            return false;
        }
        if (values == nullptr) {
            rows[0] = 0;
            columns[0] = 0;
            rows[1] = 1;
            columns[1] = 0;
            rows[2] = 1;
            columns[2] = 1;
            return true;
        }
        try {
            const ObjectiveEvaluation evaluation = evaluate(x);
            values[0] = objective_factor * evaluation.hessian[0];
            values[1] = objective_factor * evaluation.hessian[2];
            values[2] = objective_factor * evaluation.hessian[3];
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
    }

    void finalize_solution(
        Ipopt::SolverReturn,
        Ipopt::Index n,
        const Ipopt::Number* x,
        const Ipopt::Number*,
        const Ipopt::Number*,
        Ipopt::Index m,
        const Ipopt::Number*,
        const Ipopt::Number*,
        Ipopt::Number,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*
    ) override {
        if (n == 2 && m == 0 && x != nullptr) {
            std::copy(x, x + n, solution_.begin());
        }
    }

    [[nodiscard]] const std::array<double, 2>& solution() const { return solution_; }
    [[nodiscard]] const std::string& callback_error() const { return callback_error_; }

private:
    struct ObjectiveEvaluation {
        double objective = 0.0;
        std::array<double, 2> gradient{};
        std::array<double, 4> hessian{};
    };

    [[nodiscard]] ObjectiveEvaluation evaluate(const Ipopt::Number* x) const {
        if (objective_ == StageIObjective::homogeneous) {
            const HeldStateEvaluation evaluation = evaluate_held_state(
                provider_,
                temperature_k_,
                pressure_pa_,
                x[0],
                x[1]
            );
            return {evaluation.g_bar, evaluation.gradient, evaluation.hessian};
        }
        if (objective_ == StageIObjective::tunneling) {
            const HeldTunnelingEvaluation evaluation = evaluate_held_tunneling(
                provider_,
                temperature_k_,
                pressure_pa_,
                reference_,
                minimum_x_methane_,
                minimum_tpd_,
                x[0],
                x[1]
            );
            return {evaluation.objective, evaluation.gradient, evaluation.hessian};
        }
        const HeldTpdEvaluation evaluation = evaluate_held_tpd(
            provider_,
            temperature_k_,
            pressure_pa_,
            reference_,
            x[0],
            x[1]
        );
        return {evaluation.d_bar, evaluation.gradient, evaluation.hessian};
    }

    const ProviderContext& provider_;
    double temperature_k_;
    double pressure_pa_;
    HeldStateEvaluation reference_;
    double minimum_x_methane_;
    double minimum_tpd_;
    StageIObjective objective_;
    std::array<double, 2> initial_;
    std::array<double, 2> lower_;
    std::array<double, 2> upper_;
    std::array<double, 2> solution_{};
    std::string callback_error_;
};

SearchRun run_search_ipopt(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const HeldStateEvaluation& reference,
    double minimum_x_methane,
    double minimum_tpd,
    StageIObjective objective,
    const std::array<double, 2>& initial,
    const std::array<double, 2>& lower,
    const std::array<double, 2>& upper
) {
    auto* raw_problem = new SearchTnlp(
        provider,
        temperature_k,
        pressure_pa,
        reference,
        minimum_x_methane,
        minimum_tpd,
        objective,
        initial,
        lower,
        upper
    );
    Ipopt::SmartPtr<Ipopt::TNLP> problem = raw_problem;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application = IpoptApplicationFactory();
    configure_ipopt(application);

    SearchRun result;
    result.initial_guess = initial;
    const Ipopt::ApplicationReturnStatus initialize_status = application->Initialize();
    if (initialize_status != Ipopt::Solve_Succeeded) {
        result.solver_status = "initialization_" + ipopt_status_name(initialize_status);
        return result;
    }
    const Ipopt::ApplicationReturnStatus status = application->OptimizeTNLP(problem);
    result.solver_status = ipopt_status_name(status);
    result.solver_converged = status == Ipopt::Solve_Succeeded;
    result.variables = raw_problem->solution();
    result.callback_error = raw_problem->callback_error();
    const Ipopt::SmartPtr<Ipopt::SolveStatistics> statistics = application->Statistics();
    if (IsValid(statistics)) {
        result.iterations = statistics->IterationCount();
    }
    return result;
}

bool volume_is_interior(double log_volume) {
    return log_volume > std::log(kVolumeLowerM3PerMol) + kBoundMargin
        && log_volume < std::log(kVolumeUpperM3PerMol) - kBoundMargin;
}

bool locally_accepts(const SearchRun& run, const HeldTpdEvaluation& evaluation) {
    return run.solver_converged
        && run.callback_error.empty()
        && std::isfinite(evaluation.d_bar)
        && std::abs(evaluation.state.pressure_stationarity_relative)
            <= kLocalPhysicalTolerance
        && volume_is_interior(evaluation.state.log_volume);
}

std::vector<HeldStageIStart> make_planned_starts(
    double temperature_k,
    double pressure_pa,
    double feed_x_methane
) {
    const double liquid_log_volume = std::log(1.0 / 20'000.0);
    const double vapor_volume = std::clamp(
        2.0 * kGasConstantJPerMolK * temperature_k / pressure_pa,
        1.01 * kVolumeLowerM3PerMol,
        0.99 * kVolumeUpperM3PerMol
    );
    const double vapor_log_volume = std::log(vapor_volume);
    std::vector<HeldStageIStart> starts;
    starts.reserve(20);
    const auto append_pair = [&](const std::string& role, double x_methane) {
        const double bounded_x = std::clamp(x_methane, kCompositionLower, kCompositionUpper);
        starts.push_back({role + "_liquid", bounded_x, liquid_log_volume});
        starts.push_back({role + "_vapor", bounded_x, vapor_log_volume});
    };
    append_pair("feed_near_lower", feed_x_methane - 0.02);
    append_pair("feed_near_upper", feed_x_methane + 0.02);
    append_pair("component_1_rich", 1.0e-4);
    append_pair("component_2_rich", 1.0 - 1.0e-4);
    append_pair("stratified_low", 0.1);
    append_pair("stratified_high", 0.9);
    append_pair("stratified_quarter", 0.25);
    append_pair("stratified_three_quarter", 0.75);
    append_pair("stratified_inner_low", 0.4);
    append_pair("stratified_inner_high", 0.6);
    return starts;
}

HeldStageIAttempt make_attempt(
    const std::string& kind,
    const std::string& role,
    const SearchRun& run,
    bool accepted,
    bool materially_perturbed,
    double objective,
    double tpd,
    double pressure_stationarity_relative
) {
    return {
        kind,
        role,
        run.initial_guess,
        run.solver_converged,
        run.solver_status,
        run.iterations,
        accepted,
        materially_perturbed,
        objective,
        tpd,
        pressure_stationarity_relative,
        run.callback_error,
    };
}

double confirmation_difference(
    const HeldTpdEvaluation& candidate,
    const HeldTpdEvaluation& confirmation
) {
    return std::max({
        std::abs(candidate.d_bar - confirmation.d_bar),
        std::abs(candidate.state.x_methane - confirmation.state.x_methane),
        std::abs(candidate.state.volume_m3 - confirmation.state.volume_m3)
            / std::max({candidate.state.volume_m3, confirmation.state.volume_m3, 1.0e-12}),
    });
}

}  // namespace

HeldStageIResult solve_held_stage_i(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane
) {
    HeldStageIResult result;
    result.planned_starts = make_planned_starts(
        temperature_k,
        pressure_pa,
        feed_x_methane
    );
    if (!std::isfinite(temperature_k) || !std::isfinite(pressure_pa)
        || !std::isfinite(feed_x_methane) || temperature_k <= 0.0 || pressure_pa <= 0.0
        || feed_x_methane <= kCompositionLower
        || feed_x_methane >= kCompositionUpper) {
        result.search_status = "invalid_input";
        result.failure_reason = "Stage I requires finite positive T and P and an interior feed";
        return result;
    }

    const std::array<std::pair<std::string, double>, 2> reference_starts{{
        {"liquid", std::log(1.0 / 20'000.0)},
        {"vapor", std::log(std::clamp(
            2.0 * kGasConstantJPerMolK * temperature_k / pressure_pa,
            1.01 * kVolumeLowerM3PerMol,
            0.99 * kVolumeUpperM3PerMol
        ))},
    }};
    double lowest_reference = std::numeric_limits<double>::infinity();
    for (const auto& [role, initial] : reference_starts) {
        const SearchRun run = run_search_ipopt(
            provider,
            temperature_k,
            pressure_pa,
            {},
            feed_x_methane,
            0.0,
            StageIObjective::homogeneous,
            {feed_x_methane, initial},
            {feed_x_methane, std::log(kVolumeLowerM3PerMol)},
            {feed_x_methane, std::log(kVolumeUpperM3PerMol)}
        );
        HeldReferenceAttempt record;
        record.role = role;
        record.initial_log_volume = initial;
        record.solver_converged = run.solver_converged;
        record.solver_status = run.solver_status;
        record.iterations = run.iterations;
        record.g_bar = std::numeric_limits<double>::infinity();
        record.log_volume = run.variables[1];
        record.pressure_stationarity_relative = std::numeric_limits<double>::infinity();
        record.callback_error = run.callback_error;
        if (run.solver_converged && run.callback_error.empty()) {
            try {
                HeldStateEvaluation candidate = evaluate_held_state(
                    provider,
                    temperature_k,
                    pressure_pa,
                    feed_x_methane,
                    run.variables[1]
                );
                record.g_bar = candidate.g_bar;
                record.pressure_stationarity_relative =
                    candidate.pressure_stationarity_relative;
                record.accepted = std::abs(candidate.pressure_stationarity_relative)
                        <= kLocalPhysicalTolerance
                    && volume_is_interior(candidate.log_volume);
                if (record.accepted && candidate.g_bar < lowest_reference) {
                    lowest_reference = candidate.g_bar;
                    result.reference = std::move(candidate);
                    result.has_reference = true;
                }
            } catch (const std::exception& error) {
                record.callback_error = error.what();
            }
        }
        result.reference_attempts.push_back(std::move(record));
    }
    if (!result.has_reference) {
        result.search_status = "reference_failed";
        result.failure_reason = "homogeneous feed reference search failed";
        return result;
    }

    result.search_status = "searching";
    result.best_tpd = 0.0;
    result.best_state = result.reference;
    bool search_failed = false;
    for (const HeldStageIStart& start : result.planned_starts) {
        const bool lower_side = start.x_methane < result.best_state.x_methane;
        std::array<double, 2> tunnel_lower{
            lower_side
                ? kCompositionLower
                : result.best_state.x_methane + kHeldTunnelPoleExclusion,
            std::log(kVolumeLowerM3PerMol),
        };
        std::array<double, 2> tunnel_upper{
            lower_side
                ? result.best_state.x_methane - kHeldTunnelPoleExclusion
                : kCompositionUpper,
            std::log(kVolumeUpperM3PerMol),
        };
        if (tunnel_lower[0] >= tunnel_upper[0]
            || start.x_methane <= tunnel_lower[0]
            || start.x_methane >= tunnel_upper[0]) {
            SearchRun singular;
            singular.initial_guess = {start.x_methane, start.log_volume};
            singular.solver_status = "singular_restart";
            singular.callback_error =
                "tunneling restart lies in the singular composition neighborhood";
            result.attempt_log.push_back(make_attempt(
                "tunneling",
                start.role,
                singular,
                false,
                false,
                std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity()
            ));
            search_failed = true;
            continue;
        }

        const SearchRun tunnel = run_search_ipopt(
            provider,
            temperature_k,
            pressure_pa,
            result.reference,
            result.best_state.x_methane,
            result.best_tpd,
            StageIObjective::tunneling,
            {start.x_methane, start.log_volume},
            tunnel_lower,
            tunnel_upper
        );
        bool tunnel_accepted = false;
        HeldTunnelingEvaluation tunnel_evaluation;
        if (tunnel.solver_converged && tunnel.callback_error.empty()) {
            try {
                tunnel_evaluation = evaluate_held_tunneling(
                    provider,
                    temperature_k,
                    pressure_pa,
                    result.reference,
                    result.best_state.x_methane,
                    result.best_tpd,
                    tunnel.variables[0],
                    tunnel.variables[1]
                );
                tunnel_accepted = std::isfinite(tunnel_evaluation.objective)
                    && std::abs(
                           tunnel_evaluation.tpd.state.pressure_stationarity_relative
                       )
                        <= kLocalPhysicalTolerance
                    && volume_is_interior(tunnel_evaluation.tpd.state.log_volume);
            } catch (const std::exception&) {
                tunnel_accepted = false;
            }
        }
        result.attempt_log.push_back(make_attempt(
            "tunneling",
            start.role,
            tunnel,
            tunnel_accepted,
            false,
            tunnel_accepted ? tunnel_evaluation.objective
                            : std::numeric_limits<double>::infinity(),
            tunnel_accepted ? tunnel_evaluation.tpd.d_bar
                            : std::numeric_limits<double>::infinity(),
            tunnel_accepted ? tunnel_evaluation.tpd.state.pressure_stationarity_relative
                            : std::numeric_limits<double>::infinity()
        ));
        if (!tunnel_accepted) {
            search_failed = true;
            continue;
        }

        const std::array<double, 2> search_lower{
            kCompositionLower,
            std::log(kVolumeLowerM3PerMol),
        };
        const std::array<double, 2> search_upper{
            kCompositionUpper,
            std::log(kVolumeUpperM3PerMol),
        };
        const SearchRun polish = run_search_ipopt(
            provider,
            temperature_k,
            pressure_pa,
            result.reference,
            result.best_state.x_methane,
            result.best_tpd,
            StageIObjective::tpd,
            tunnel.variables,
            search_lower,
            search_upper
        );
        HeldTpdEvaluation candidate;
        bool polish_accepted = false;
        if (polish.solver_converged && polish.callback_error.empty()) {
            try {
                candidate = evaluate_held_tpd(
                    provider,
                    temperature_k,
                    pressure_pa,
                    result.reference,
                    polish.variables[0],
                    polish.variables[1]
                );
                polish_accepted = locally_accepts(polish, candidate);
            } catch (const std::exception&) {
                polish_accepted = false;
            }
        }
        result.attempt_log.push_back(make_attempt(
            "tpd_polish",
            start.role,
            polish,
            polish_accepted,
            false,
            polish_accepted ? candidate.d_bar : std::numeric_limits<double>::infinity(),
            polish_accepted ? candidate.d_bar : std::numeric_limits<double>::infinity(),
            polish_accepted ? candidate.state.pressure_stationarity_relative
                            : std::numeric_limits<double>::infinity()
        ));
        if (!polish_accepted) {
            search_failed = true;
            continue;
        }
        ++result.starts_completed;
        if (candidate.d_bar < result.best_tpd) {
            result.best_tpd = candidate.d_bar;
            result.best_state = candidate.state;
        }
        if (candidate.d_bar >= kTpdNegativeThreshold) {
            continue;
        }

        std::array<double, 2> confirmation_initial = polish.variables;
        const double composition_shift = confirmation_initial[0] + 0.02 < kCompositionUpper
            ? 0.02
            : -0.02;
        const double volume_shift = confirmation_initial[1] + 0.05 < search_upper[1]
            ? 0.05
            : -0.05;
        confirmation_initial[0] += composition_shift;
        confirmation_initial[1] += volume_shift;
        const bool materially_perturbed = std::abs(composition_shift) >= 0.01
            && std::abs(volume_shift) >= 0.04;
        const SearchRun confirmation = run_search_ipopt(
            provider,
            temperature_k,
            pressure_pa,
            result.reference,
            result.best_state.x_methane,
            result.best_tpd,
            StageIObjective::tpd,
            confirmation_initial,
            search_lower,
            search_upper
        );
        ++result.negative_confirmations;
        HeldTpdEvaluation confirmed;
        bool confirmation_accepted = false;
        if (confirmation.solver_converged && confirmation.callback_error.empty()) {
            try {
                confirmed = evaluate_held_tpd(
                    provider,
                    temperature_k,
                    pressure_pa,
                    result.reference,
                    confirmation.variables[0],
                    confirmation.variables[1]
                );
                confirmation_accepted = locally_accepts(confirmation, confirmed)
                    && confirmed.d_bar < kTpdNegativeThreshold;
                if (confirmation_accepted) {
                    result.confirmation_max_difference = confirmation_difference(
                        candidate,
                        confirmed
                    );
                    confirmation_accepted = result.confirmation_max_difference
                        <= kConfirmationTolerance;
                }
            } catch (const std::exception&) {
                confirmation_accepted = false;
            }
        }
        result.attempt_log.push_back(make_attempt(
            "confirmation",
            start.role,
            confirmation,
            confirmation_accepted,
            materially_perturbed,
            confirmation_accepted ? confirmed.d_bar
                                  : std::numeric_limits<double>::infinity(),
            confirmation_accepted ? confirmed.d_bar
                                  : std::numeric_limits<double>::infinity(),
            confirmation_accepted ? confirmed.state.pressure_stationarity_relative
                                  : std::numeric_limits<double>::infinity()
        ));
        if (confirmation_accepted) {
            result.outcome = "negative_tpd";
            result.search_status = "confirmed_negative";
            result.failure_reason.clear();
            return result;
        }
        search_failed = true;
    }

    if (search_failed) {
        result.outcome = "indeterminate";
        result.search_status = "search_failed";
        result.failure_reason = "one or more declared Stage I attempts failed";
        return result;
    }
    result.outcome = "no_negative_found";
    result.search_status = "source_heuristic_complete";
    result.failure_reason.clear();
    return result;
}

}  // namespace epcsaft_equilibrium
