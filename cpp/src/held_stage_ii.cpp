#include "held.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>
#include <coin/IpTNLP.hpp>

namespace epcsaft_equilibrium {
namespace {

constexpr double kFloatingPointAllowance = 256.0 * std::numeric_limits<double>::epsilon();
constexpr double kGasConstantJPerMolK = 8.31446261815324;
constexpr double kCompositionLower = 1.0e-8;
constexpr double kCompositionUpper = 1.0 - 1.0e-8;
constexpr double kVolumeLowerM3PerMol = 1.0e-5;
constexpr double kVolumeUpperM3PerMol = 1.0e-1;
constexpr double kLocalPhysicalTolerance = 1.0e-8;
constexpr double kBoundMargin = 1.0e-7;
constexpr double kRelativeDensitySeparation = 1.0e-3;
constexpr double kCandidateGapTolerance = 1.0e-2;
constexpr double kCandidateMultiplierTolerance = 0.5;
constexpr double kCandidateCompositionSeparation = 1.0e-3;
constexpr int kStageIIMajorIterations = 100;
constexpr int kStageIILowerStarts = 20;
constexpr int kIpoptIterations = 300;

double numerical_allowance(double left, double right) {
    return kFloatingPointAllowance * std::max({1.0, std::abs(left), std::abs(right)});
}

bool nearly_equal(double left, double right) {
    return std::abs(left - right) <= numerical_allowance(left, right);
}

struct OuterVertex {
    double value = 0.0;
    double multiplier = 0.0;
};

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

struct LowerRun {
    bool solver_converged = false;
    std::string solver_status;
    std::string callback_error;
    int iterations = 0;
    std::array<double, 2> initial_guess{};
    std::array<double, 2> variables{};
};

class LowerTnlp final : public Ipopt::TNLP {
public:
    LowerTnlp(
        const ProviderContext& provider,
        double temperature_k,
        double pressure_pa,
        double feed_x_methane,
        double multiplier,
        std::array<double, 2> initial,
        std::array<double, 2> lower,
        std::array<double, 2> upper
    )
        : provider_(provider),
          temperature_k_(temperature_k),
          pressure_pa_(pressure_pa),
          feed_x_methane_(feed_x_methane),
          multiplier_(multiplier),
          initial_(initial),
          lower_(lower),
          upper_(upper),
          solution_(initial) {}

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
            const HeldLowerEvaluation evaluation = evaluate(x);
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
            const HeldLowerEvaluation evaluation = evaluate(x);
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
    [[nodiscard]] HeldLowerEvaluation evaluate(const Ipopt::Number* x) const {
        return evaluate_held_lower(
            provider_,
            temperature_k_,
            pressure_pa_,
            feed_x_methane_,
            multiplier_,
            x[0],
            x[1]
        );
    }

    const ProviderContext& provider_;
    double temperature_k_;
    double pressure_pa_;
    double feed_x_methane_;
    double multiplier_;
    std::array<double, 2> initial_;
    std::array<double, 2> lower_;
    std::array<double, 2> upper_;
    std::array<double, 2> solution_{};
    std::string callback_error_;
};

LowerRun run_lower_ipopt(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane,
    double multiplier,
    const std::array<double, 2>& initial,
    const std::array<double, 2>& lower,
    const std::array<double, 2>& upper
) {
    auto* raw_problem = new LowerTnlp(
        provider,
        temperature_k,
        pressure_pa,
        feed_x_methane,
        multiplier,
        initial,
        lower,
        upper
    );
    Ipopt::SmartPtr<Ipopt::TNLP> problem = raw_problem;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application = IpoptApplicationFactory();
    configure_ipopt(application);

    LowerRun result;
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

bool locally_accepts(const LowerRun& run, const HeldLowerEvaluation& evaluation) {
    return run.solver_converged && run.callback_error.empty()
        && std::isfinite(evaluation.objective)
        && std::abs(evaluation.state.pressure_stationarity_relative)
            <= kLocalPhysicalTolerance
        && volume_is_interior(evaluation.state.log_volume);
}

bool density_is_distinct(
    const HeldStateEvaluation& candidate,
    const std::string& endpoint,
    const std::vector<HeldStageIICut>& retained
) {
    const double candidate_density = 1.0 / candidate.volume_m3;
    return std::all_of(
        retained.begin(),
        retained.end(),
        [candidate_density, &endpoint](const HeldStageIICut& cut) {
            if (cut.endpoint != endpoint) {
                return true;
            }
            const double retained_density = 1.0 / cut.state.volume_m3;
            return std::abs(candidate_density - retained_density)
                    / std::max(candidate_density, retained_density)
                >= kRelativeDensitySeparation;
        }
    );
}

std::string indexed_role(const std::string& start_class, int index) {
    return start_class + "_" + (index < 10 ? "0" : "") + std::to_string(index);
}

std::vector<HeldStageIIStart> make_lower_starts(
    double temperature_k,
    double pressure_pa,
    double feed_x_methane,
    const std::vector<HeldStateEvaluation>& previous_states
) {
    const double liquid_log_volume = std::log(1.0 / 20'000.0);
    const double vapor_log_volume = std::log(std::clamp(
        2.0 * kGasConstantJPerMolK * temperature_k / pressure_pa,
        1.01 * kVolumeLowerM3PerMol,
        0.99 * kVolumeUpperM3PerMol
    ));
    constexpr std::array<double, 7> stratified_compositions{
        0.1, 0.9, 0.25, 0.75, 0.4, 0.6, 0.5,
    };
    std::vector<HeldStageIIStart> starts;
    starts.reserve(kStageIILowerStarts);
    for (int index = 0; index < kStageIILowerStarts; ++index) {
        const int group = index / 3;
        if (index % 3 == 0) {
            const HeldStateEvaluation* previous = previous_states.empty()
                ? nullptr
                : &previous_states[
                      previous_states.size() - 1
                      - static_cast<std::size_t>(group) % previous_states.size()
                  ];
            const double base_x = previous == nullptr ? feed_x_methane : previous->x_methane;
            const double base_log_volume =
                previous == nullptr ? liquid_log_volume : previous->log_volume;
            const double shift = group % 2 == 0 ? 0.01 : -0.01;
            starts.push_back({
                "shifted_previous",
                indexed_role("shifted_previous", index),
                std::clamp(base_x + shift, kCompositionLower, kCompositionUpper),
                base_log_volume,
            });
        } else if (index % 3 == 1) {
            starts.push_back({
                "component_near",
                indexed_role("component_near", index),
                group % 2 == 0 ? 1.0 - 1.0e-4 : 1.0e-4,
                group % 2 == 0 ? liquid_log_volume : vapor_log_volume,
            });
        } else {
            starts.push_back({
                "stratified",
                indexed_role("stratified", index),
                stratified_compositions[static_cast<std::size_t>(group)
                                        % stratified_compositions.size()],
                group % 2 == 0 ? vapor_log_volume : liquid_log_volume,
            });
        }
    }
    return starts;
}

double relative_density_difference(double left_volume, double right_volume) {
    const double left_density = 1.0 / left_volume;
    const double right_density = 1.0 / right_volume;
    return std::abs(left_density - right_density)
        / std::max(left_density, right_density);
}

}  // namespace

HeldLowerEvaluation evaluate_held_lower(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane,
    double multiplier,
    double x_methane,
    double log_volume
) {
    if (!std::isfinite(feed_x_methane) || !std::isfinite(multiplier)) {
        throw std::invalid_argument("HELD lower feed and multiplier must be finite");
    }
    HeldStateEvaluation state = evaluate_held_state(
        provider,
        temperature_k,
        pressure_pa,
        x_methane,
        log_volume
    );
    HeldLowerEvaluation result;
    result.objective = state.g_bar + multiplier * (feed_x_methane - x_methane);
    result.gradient = {state.gradient[0] - multiplier, state.gradient[1]};
    result.hessian = state.hessian;
    result.state = std::move(state);
    if (!std::isfinite(result.objective)
        || !std::all_of(result.gradient.begin(), result.gradient.end(), [](double value) {
               return std::isfinite(value);
           })) {
        throw std::overflow_error("HELD lower evaluation must be finite");
    }
    return result;
}

HeldOuterResult solve_held_outer_envelope(const std::vector<HeldOuterCut>& cuts) {
    HeldOuterResult result;
    if (std::any_of(cuts.begin(), cuts.end(), [](const HeldOuterCut& cut) {
            return cut.intercept == -std::numeric_limits<double>::infinity();
        })) {
        result.failure_reason = "a cut has no finite feasible upper value";
        return result;
    }
    if (cuts.empty()
        || std::any_of(cuts.begin(), cuts.end(), [](const HeldOuterCut& cut) {
               return cut.identity.empty() || !std::isfinite(cut.intercept)
                   || !std::isfinite(cut.slope);
           })) {
        result.status = "invalid_input";
        result.failure_reason = "outer cuts require identities and finite coefficients";
        return result;
    }
    const bool has_positive_slope = std::any_of(
        cuts.begin(), cuts.end(), [](const HeldOuterCut& cut) { return cut.slope > 0.0; }
    );
    const bool has_negative_slope = std::any_of(
        cuts.begin(), cuts.end(), [](const HeldOuterCut& cut) { return cut.slope < 0.0; }
    );
    if (!has_positive_slope || !has_negative_slope) {
        result.status = "unbounded";
        result.failure_reason =
            "cuts do not bracket the feed with opposite endpoint slopes";
        return result;
    }

    std::vector<OuterVertex> feasible_vertices;
    for (std::size_t left_index = 0; left_index < cuts.size(); ++left_index) {
        const HeldOuterCut& left = cuts[left_index];
        for (std::size_t right_index = left_index + 1; right_index < cuts.size();
             ++right_index) {
            const HeldOuterCut& right = cuts[right_index];
            const double slope_delta = left.slope - right.slope;
            if (nearly_equal(slope_delta, 0.0)) {
                continue;
            }
            const double multiplier = (right.intercept - left.intercept) / slope_delta;
            const double value = left.intercept + left.slope * multiplier;
            if (!std::isfinite(multiplier) || !std::isfinite(value)) {
                continue;
            }
            const bool feasible = std::all_of(
                cuts.begin(), cuts.end(), [value, multiplier](const HeldOuterCut& cut) {
                    const double upper = cut.intercept + cut.slope * multiplier;
                    return std::isfinite(upper)
                        && value <= upper + numerical_allowance(value, upper);
                }
            );
            if (feasible) {
                feasible_vertices.push_back({value, multiplier});
            }
        }
    }
    if (feasible_vertices.empty()) {
        result.failure_reason = "no feasible finite cut intersection was found";
        return result;
    }

    const double best_value = std::max_element(
        feasible_vertices.begin(),
        feasible_vertices.end(),
        [](const OuterVertex& left, const OuterVertex& right) {
            return left.value < right.value;
        }
    )->value;
    for (const OuterVertex& vertex : feasible_vertices) {
        if (nearly_equal(vertex.value, best_value)) {
            result.tied_multipliers.push_back(vertex.multiplier);
        }
    }
    std::sort(result.tied_multipliers.begin(), result.tied_multipliers.end());
    result.tied_multipliers.erase(
        std::unique(
            result.tied_multipliers.begin(),
            result.tied_multipliers.end(),
            [](double left, double right) { return nearly_equal(left, right); }
        ),
        result.tied_multipliers.end()
    );

    result.status = "finite";
    result.value = best_value;
    result.multiplier = result.tied_multipliers.front();
    for (const HeldOuterCut& cut : cuts) {
        const double upper = cut.intercept + cut.slope * result.multiplier;
        if (nearly_equal(result.value, upper)) {
            result.active_cut_ids.push_back(cut.identity);
        }
    }
    return result;
}

HeldStageIIInitialization initialize_held_stage_ii_cuts(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane
) {
    HeldStageIIInitialization result;
    if (!std::isfinite(temperature_k) || !std::isfinite(pressure_pa)
        || !std::isfinite(feed_x_methane) || temperature_k <= 0.0 || pressure_pa <= 0.0
        || feed_x_methane <= kCompositionLower
        || feed_x_methane >= kCompositionUpper) {
        result.status = "invalid_input";
        result.failure_reason =
            "Stage II requires finite positive T and P and an interior feed";
        return result;
    }

    const double liquid_log_volume = std::log(1.0 / 20'000.0);
    const double vapor_log_volume = std::log(std::clamp(
        2.0 * kGasConstantJPerMolK * temperature_k / pressure_pa,
        1.01 * kVolumeLowerM3PerMol,
        0.99 * kVolumeUpperM3PerMol
    ));
    const std::array<std::tuple<std::string, std::string, double, double>, 4> starts{{
        {"lower_endpoint_liquid", "lower", kCompositionLower, liquid_log_volume},
        {"lower_endpoint_vapor", "lower", kCompositionLower, vapor_log_volume},
        {"upper_endpoint_liquid", "upper", kCompositionUpper, liquid_log_volume},
        {"upper_endpoint_vapor", "upper", kCompositionUpper, vapor_log_volume},
    }};

    for (const auto& [role, endpoint, composition, initial_log_volume] : starts) {
        LowerRun run = run_lower_ipopt(
            provider,
            temperature_k,
            pressure_pa,
            feed_x_methane,
            0.0,
            {composition, initial_log_volume},
            {composition, std::log(kVolumeLowerM3PerMol)},
            {composition, std::log(kVolumeUpperM3PerMol)}
        );
        HeldStageIIAttempt attempt;
        attempt.role = role;
        attempt.initial_guess = run.initial_guess;
        attempt.solver_converged = run.solver_converged;
        attempt.solver_status = run.solver_status;
        attempt.iterations = run.iterations;
        attempt.callback_error = run.callback_error;
        HeldLowerEvaluation evaluation;
        if (run.solver_converged && run.callback_error.empty()) {
            try {
                evaluation = evaluate_held_lower(
                    provider,
                    temperature_k,
                    pressure_pa,
                    feed_x_methane,
                    0.0,
                    composition,
                    run.variables[1]
                );
                attempt.accepted = locally_accepts(run, evaluation);
                if (attempt.accepted) {
                    attempt.objective = evaluation.objective;
                    attempt.pressure_stationarity_relative =
                        evaluation.state.pressure_stationarity_relative;
                    if (density_is_distinct(evaluation.state, endpoint, result.cuts)) {
                        HeldOuterCut outer{
                            role,
                            evaluation.state.g_bar,
                            feed_x_methane - evaluation.state.x_methane,
                        };
                        result.cuts.push_back({
                            role,
                            endpoint,
                            std::move(evaluation.state),
                            std::move(outer),
                        });
                    }
                }
            } catch (const std::exception& error) {
                attempt.callback_error = error.what();
            }
        }
        result.attempts.push_back(std::move(attempt));
    }

    const bool has_lower = std::any_of(
        result.cuts.begin(), result.cuts.end(), [](const HeldStageIICut& cut) {
            return cut.endpoint == "lower";
        }
    );
    const bool has_upper = std::any_of(
        result.cuts.begin(), result.cuts.end(), [](const HeldStageIICut& cut) {
            return cut.endpoint == "upper";
        }
    );
    if (!has_lower || !has_upper) {
        result.status = "endpoint_failed";
        result.failure_reason = "valid endpoint cuts were not found on both feed sides";
        return result;
    }
    std::vector<HeldOuterCut> outer_cuts;
    outer_cuts.reserve(result.cuts.size());
    for (const HeldStageIICut& cut : result.cuts) {
        outer_cuts.push_back(cut.outer);
    }
    result.outer = solve_held_outer_envelope(outer_cuts);
    if (result.outer.status != "finite") {
        result.status = "endpoint_failed";
        result.failure_reason = result.outer.failure_reason;
        return result;
    }
    result.status = "initialized";
    return result;
}

HeldStageIILowerSearch search_held_stage_ii_lower(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane,
    double multiplier,
    double upper_bound,
    const std::vector<HeldStateEvaluation>& previous_states,
    const std::string& identity_prefix
) {
    HeldStageIILowerSearch result;
    result.upper_bound = upper_bound;
    result.multiplier = multiplier;
    if (!std::isfinite(temperature_k) || !std::isfinite(pressure_pa)
        || !std::isfinite(feed_x_methane) || !std::isfinite(multiplier)
        || !std::isfinite(upper_bound) || temperature_k <= 0.0 || pressure_pa <= 0.0
        || feed_x_methane <= kCompositionLower
        || feed_x_methane >= kCompositionUpper || identity_prefix.empty()) {
        result.status = "invalid_input";
        return result;
    }
    result.planned_starts = make_lower_starts(
        temperature_k,
        pressure_pa,
        feed_x_methane,
        previous_states
    );
    const std::array<double, 2> lower{
        kCompositionLower,
        std::log(kVolumeLowerM3PerMol),
    };
    const std::array<double, 2> upper{
        kCompositionUpper,
        std::log(kVolumeUpperM3PerMol),
    };
    for (std::size_t index = 0; index < result.planned_starts.size(); ++index) {
        const HeldStageIIStart& start = result.planned_starts[index];
        LowerRun run = run_lower_ipopt(
            provider,
            temperature_k,
            pressure_pa,
            feed_x_methane,
            multiplier,
            {start.x_methane, start.log_volume},
            lower,
            upper
        );
        HeldStageIIAttempt attempt;
        attempt.role = start.role;
        attempt.initial_guess = run.initial_guess;
        attempt.solver_converged = run.solver_converged;
        attempt.solver_status = run.solver_status;
        attempt.iterations = run.iterations;
        attempt.callback_error = run.callback_error;
        HeldLowerEvaluation evaluation;
        if (run.solver_converged && run.callback_error.empty()) {
            try {
                evaluation = evaluate_held_lower(
                    provider,
                    temperature_k,
                    pressure_pa,
                    feed_x_methane,
                    multiplier,
                    run.variables[0],
                    run.variables[1]
                );
                attempt.accepted = locally_accepts(run, evaluation);
                if (attempt.accepted) {
                    attempt.objective = evaluation.objective;
                    attempt.pressure_stationarity_relative =
                        evaluation.state.pressure_stationarity_relative;
                    const std::string identity = identity_prefix + "_s"
                        + (index < 10 ? "0" : "") + std::to_string(index);
                    HeldOuterCut outer{
                        identity,
                        evaluation.state.g_bar,
                        feed_x_methane - evaluation.state.x_methane,
                    };
                    result.accepted_cuts.push_back({
                        identity,
                        "",
                        std::move(evaluation.state),
                        std::move(outer),
                    });
                } else {
                    attempt.callback_error = "local physical acceptance failed";
                }
            } catch (const std::exception& error) {
                attempt.callback_error = error.what();
            }
        }
        result.attempts.push_back(std::move(attempt));
        ++result.starts_completed;
        if (!result.accepted_cuts.empty()
            && result.accepted_cuts.back().outer.intercept
                    + multiplier * result.accepted_cuts.back().outer.slope
                <= upper_bound) {
            result.status = "below_upper_found";
            return result;
        }
    }
    result.status = "search_exhausted";
    return result;
}

HeldStageIICandidateResult select_held_stage_ii_candidates(
    double feed_x_methane,
    double upper_bound,
    double multiplier,
    const std::vector<HeldCandidateInput>& points
) {
    HeldStageIICandidateResult result;
    for (const HeldCandidateInput& point : points) {
        const double lower_objective =
            point.g_bar + multiplier * (feed_x_methane - point.x_methane);
        if (!std::isfinite(lower_objective) || !std::isfinite(point.volume_m3)
            || point.volume_m3 <= 0.0) {
            result.rejections.push_back({point.identity, "invalid_point"});
            continue;
        }
        const double upper_lower_gap = upper_bound - lower_objective;
        // A pre-stop local point remains a valid outer cut, but only a lower
        // point satisfying the declared L <= UBD inner-stop condition may
        // enter the candidate set.
        if (upper_lower_gap < -numerical_allowance(upper_bound, lower_objective)) {
            result.rejections.push_back({point.identity, "above_upper"});
            continue;
        }
        if (upper_lower_gap > kCandidateGapTolerance) {
            result.rejections.push_back({point.identity, "upper_lower_gap"});
            continue;
        }
        const double multiplier_residual =
            std::abs(point.composition_gradient - multiplier)
            / std::max(std::abs(multiplier), 1.0);
        if (multiplier_residual > kCandidateMultiplierTolerance) {
            result.rejections.push_back({point.identity, "multiplier_residual"});
            continue;
        }

        HeldStageIICandidate candidate{
            point.identity,
            lower_objective,
            point.x_methane,
            point.volume_m3,
        };
        auto duplicate = std::find_if(
            result.candidates.begin(),
            result.candidates.end(),
            [&candidate](const HeldStageIICandidate& retained) {
                return std::abs(candidate.x_methane - retained.x_methane)
                        < kCandidateCompositionSeparation
                    && relative_density_difference(
                           candidate.volume_m3,
                           retained.volume_m3
                       )
                        < kRelativeDensitySeparation;
            }
        );
        if (duplicate == result.candidates.end()) {
            result.candidates.push_back(std::move(candidate));
            continue;
        }
        if (candidate.objective < duplicate->objective
            || (candidate.objective == duplicate->objective
                && candidate.identity < duplicate->identity)) {
            result.rejections.push_back({
                duplicate->identity,
                "duplicate_replaced_by_lower_objective",
            });
            *duplicate = std::move(candidate);
        } else {
            result.rejections.push_back({
                candidate.identity,
                "duplicate_of_lower_objective",
            });
        }
    }
    if (result.candidates.size() > 2) {
        result.status = "scope_exceeded";
    } else if (result.candidates.size() == 2) {
        result.status = "stage_iii_ready";
    }
    return result;
}

std::string held_stage_ii_budget_status(
    int major_iterations,
    int lower_starts,
    bool lower_satisfied
) {
    if (major_iterations >= kStageIIMajorIterations
        || (!lower_satisfied && lower_starts >= kStageIILowerStarts)) {
        return "search_exhausted";
    }
    return "searching";
}

HeldStageIIResult solve_held_stage_ii(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane
) {
    HeldStageIIResult result;
    const HeldStageIResult stage_i = solve_held_stage_i(
        provider,
        temperature_k,
        pressure_pa,
        feed_x_methane
    );
    result.stage_i_outcome = stage_i.outcome;
    result.stage_i_search_status = stage_i.search_status;
    result.best_tpd = stage_i.best_tpd;
    if (stage_i.outcome == "no_negative_found") {
        result.outcome = "not_required";
        result.search_status = "stage_i_no_negative";
        return result;
    }
    if (stage_i.outcome != "negative_tpd") {
        result.search_status = "stage_i_indeterminate";
        result.failure_reason = stage_i.failure_reason;
        return result;
    }

    HeldStageIIInitialization initialization = initialize_held_stage_ii_cuts(
        provider,
        temperature_k,
        pressure_pa,
        feed_x_methane
    );
    result.endpoint_attempts = std::move(initialization.attempts);
    result.cuts = std::move(initialization.cuts);
    if (initialization.status != "initialized") {
        result.search_status = "endpoint_initialization_failed";
        result.failure_reason = initialization.failure_reason;
        return result;
    }
    result.upper_bound = stage_i.reference.g_bar;
    result.search_status = "searching";

    for (int major_iteration = 0; major_iteration < kStageIIMajorIterations;
         ++major_iteration) {
        std::vector<HeldOuterCut> outer_cuts;
        outer_cuts.reserve(result.cuts.size());
        for (const HeldStageIICut& cut : result.cuts) {
            outer_cuts.push_back(cut.outer);
        }
        const HeldOuterResult outer = solve_held_outer_envelope(outer_cuts);
        if (outer.status != "finite") {
            result.search_status = "outer_failed";
            result.failure_reason = outer.failure_reason;
            return result;
        }
        result.upper_bound = std::min(result.upper_bound, outer.value);

        std::vector<HeldStateEvaluation> previous_states;
        previous_states.reserve(result.cuts.size());
        for (const HeldStageIICut& cut : result.cuts) {
            previous_states.push_back(cut.state);
        }
        const std::string identity_prefix =
            "major_" + (major_iteration < 10 ? std::string("0") : std::string())
            + std::to_string(major_iteration);
        HeldStageIILowerSearch lower = search_held_stage_ii_lower(
            provider,
            temperature_k,
            pressure_pa,
            feed_x_methane,
            outer.multiplier,
            result.upper_bound,
            previous_states,
            identity_prefix
        );
        std::vector<std::string> accepted_cut_ids;
        accepted_cut_ids.reserve(lower.accepted_cuts.size());
        for (HeldStageIICut& cut : lower.accepted_cuts) {
            accepted_cut_ids.push_back(cut.identity);
            result.cuts.push_back(std::move(cut));
        }

        std::vector<HeldCandidateInput> candidate_points;
        candidate_points.reserve(result.cuts.size());
        for (const HeldStageIICut& cut : result.cuts) {
            candidate_points.push_back({
                cut.identity,
                cut.state.g_bar,
                cut.state.x_methane,
                cut.state.volume_m3,
                cut.state.gradient[0],
            });
        }
        const HeldStageIICandidateResult selection = select_held_stage_ii_candidates(
            feed_x_methane,
            result.upper_bound,
            outer.multiplier,
            candidate_points
        );
        HeldStageIITrace trace;
        trace.major_iteration = major_iteration;
        trace.outer_value = outer.value;
        trace.upper_bound = result.upper_bound;
        trace.multiplier = outer.multiplier;
        trace.active_cut_ids = outer.active_cut_ids;
        trace.accepted_cut_ids = std::move(accepted_cut_ids);
        trace.lower_starts_completed = lower.starts_completed;
        trace.rejections = selection.rejections;
        for (const HeldStageIICandidate& candidate : selection.candidates) {
            trace.candidate_ids.push_back(candidate.identity);
        }
        result.trace.push_back(std::move(trace));
        ++result.major_iterations;

        if (held_stage_ii_budget_status(
                result.major_iterations,
                lower.starts_completed,
                lower.status == "below_upper_found"
            )
            == "search_exhausted") {
            result.outcome = "search_exhausted";
            result.search_status = "search_exhausted";
            result.failure_reason = lower.status == "search_exhausted"
                ? "Stage II lower-search start budget exhausted"
                : "Stage II major-iteration budget exhausted";
            return result;
        }
        if (selection.status == "scope_exceeded") {
            result.outcome = "scope_exceeded";
            result.search_status = "third_candidate";
            result.failure_reason = "Stage II found more than two distinct candidates";
            return result;
        }
        if (selection.status == "stage_iii_ready") {
            for (const HeldStageIICandidate& candidate : selection.candidates) {
                const auto retained = std::find_if(
                    result.cuts.begin(),
                    result.cuts.end(),
                    [&candidate](const HeldStageIICut& cut) {
                        return cut.identity == candidate.identity;
                    }
                );
                if (retained != result.cuts.end()) {
                    result.candidates.push_back(*retained);
                }
            }
            if (result.candidates.size() != 2) {
                result.search_status = "candidate_lookup_failed";
                result.failure_reason = "Stage II candidate diagnostics lost a retained cut";
                return result;
            }
            result.outcome = "stage_iii_ready";
            result.search_status = "two_candidates";
            return result;
        }
    }

    result.outcome = "search_exhausted";
    result.search_status = "search_exhausted";
    result.failure_reason = "Stage II major-iteration budget exhausted";
    return result;
}

}  // namespace epcsaft_equilibrium
