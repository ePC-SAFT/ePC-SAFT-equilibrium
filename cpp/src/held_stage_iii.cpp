#include "held.hpp"
#include "two_phase_flash.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>
#include <coin/IpTNLP.hpp>

namespace epcsaft_equilibrium {
namespace {

constexpr double kGasConstantJPerMolK = 8.31446261815324;
constexpr double kAmountLowerMol = 1.0e-10;
constexpr double kVolumeLowerM3 = 1.0e-8;
constexpr double kVolumeUpperM3 = 1.0e-1;
constexpr double kSolverConstraintTolerance = 1.0e-10;
constexpr double kPhysicalTolerance = 1.0e-8;
constexpr double kConfirmationTolerance = 1.0e-7;
constexpr double kPhaseDistanceTolerance = 1.0e-3;
constexpr double kBoundMargin = 1.0e-7;
constexpr double kStageIIICandidateRadiusMol = 1.0e-3;
constexpr double kHeldChemicalPotentialTolerance = 1.0e-6;
constexpr double kHeldFreeEnergyTolerance = 1.0e-6;

void require_finite(double value, const char* name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
}

std::size_t lower_index(std::size_t row, std::size_t column) {
    return row * (row + 1) / 2 + column;
}

FlashPhaseEvaluation evaluate_phase(
    const ProviderContext& provider,
    double temperature_k,
    const std::array<double, 2>& amounts_mol,
    double log_volume
) {
    const double volume_m3 = std::exp(log_volume);
    require_finite(volume_m3, "phase volume");
    if (amounts_mol[0] <= 0.0 || amounts_mol[1] <= 0.0 || volume_m3 <= 0.0) {
        throw std::invalid_argument("phase amounts and volume must be positive");
    }
    MixturePhaseEvaluation evaluation = provider.evaluate_mixture(
        temperature_k,
        {amounts_mol[0], amounts_mol[1]},
        volume_m3
    );
    return {amounts_mol, volume_m3, std::move(evaluation)};
}

void add_phase_derivatives(
    const FlashPhaseEvaluation& phase,
    double pressure_over_rt,
    std::size_t offset,
    std::array<double, 6>& gradient,
    std::array<double, 21>& hessian_lower
) {
    const auto& provider_gradient = phase.provider.gradient;
    const auto& provider_hessian = phase.provider.hessian;
    const double volume = phase.volume_m3;
    gradient[offset] = provider_gradient[0];
    gradient[offset + 1] = provider_gradient[1];
    gradient[offset + 2] = (provider_gradient[2] + pressure_over_rt) * volume;

    for (std::size_t row = 0; row < 2; ++row) {
        for (std::size_t column = 0; column <= row; ++column) {
            hessian_lower[lower_index(offset + row, offset + column)] =
                provider_hessian[row * 3 + column];
        }
        hessian_lower[lower_index(offset + 2, offset + row)] =
            provider_hessian[2 * 3 + row] * volume;
    }
    hessian_lower[lower_index(offset + 2, offset + 2)] =
        provider_hessian[8] * volume * volume
        + (provider_gradient[2] + pressure_over_rt) * volume;
}

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

std::array<double, 6> solver_lower_bounds() {
    return {
        kAmountLowerMol,
        kAmountLowerMol,
        std::log(kVolumeLowerM3),
        kAmountLowerMol,
        kAmountLowerMol,
        std::log(kVolumeLowerM3),
    };
}

std::array<double, 6> solver_upper_bounds(
    const std::array<double, 2>& overall_mole_fractions
) {
    return {
        overall_mole_fractions[0],
        overall_mole_fractions[1],
        std::log(kVolumeUpperM3),
        overall_mole_fractions[0],
        overall_mole_fractions[1],
        std::log(kVolumeUpperM3),
    };
}

struct AttemptResult {
    bool solver_converged = false;
    std::string solver_status;
    std::string callback_error;
    int iterations = 0;
    double constraint_violation = std::numeric_limits<double>::infinity();
    std::array<double, 6> initial_guess{};
    std::array<double, 6> variables{};
    std::array<double, 2> equality_multipliers{};
    std::array<double, 6> lower_bound_multipliers{};
    std::array<double, 6> upper_bound_multipliers{};
};

class FlashTnlp final : public Ipopt::TNLP {
public:
    FlashTnlp(
        const ProviderContext& provider,
        double temperature_k,
        double pressure_pa,
        std::array<double, 2> overall_mole_fractions,
        std::array<double, 6> initial,
        std::array<double, 6> lower,
        std::array<double, 6> upper
    )
        : provider_(provider),
          temperature_k_(temperature_k),
          pressure_pa_(pressure_pa),
          overall_mole_fractions_(overall_mole_fractions),
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
        n = 6;
        m = 2;
        nnz_jac_g = 12;
        nnz_h_lag = 21;
        index_style = TNLP::C_STYLE;
        return true;
    }

    bool get_bounds_info(
        Ipopt::Index n,
        Ipopt::Number* x_l,
        Ipopt::Number* x_u,
        Ipopt::Index m,
        Ipopt::Number* g_l,
        Ipopt::Number* g_u
    ) override {
        if (n != 6 || m != 2) {
            return false;
        }
        std::copy(lower_.begin(), lower_.end(), x_l);
        std::copy(upper_.begin(), upper_.end(), x_u);
        std::fill(g_l, g_l + m, 0.0);
        std::fill(g_u, g_u + m, 0.0);
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
        if (n != 6 || m != 2 || !init_x || init_z || init_lambda) {
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
        if (n != 6) {
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
        if (n != 6) {
            return false;
        }
        try {
            const FlashNlpEvaluation evaluation = evaluate(x);
            std::copy(evaluation.gradient.begin(), evaluation.gradient.end(), gradient);
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
    }

    bool eval_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Index m,
        Ipopt::Number* constraints
    ) override {
        if (n != 6 || m != 2) {
            return false;
        }
        try {
            const FlashNlpEvaluation evaluation = evaluate(x);
            std::copy(evaluation.constraints.begin(), evaluation.constraints.end(), constraints);
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
    }

    bool eval_jac_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Index m,
        Ipopt::Index nonzero_count,
        Ipopt::Index* rows,
        Ipopt::Index* columns,
        Ipopt::Number* values
    ) override {
        if (n != 6 || m != 2 || nonzero_count != 12) {
            return false;
        }
        if (values == nullptr) {
            for (Ipopt::Index row = 0; row < 2; ++row) {
                for (Ipopt::Index column = 0; column < 6; ++column) {
                    rows[row * 6 + column] = row;
                    columns[row * 6 + column] = column;
                }
            }
            return true;
        }
        try {
            const FlashNlpEvaluation evaluation = evaluate(x);
            std::copy(evaluation.jacobian.begin(), evaluation.jacobian.end(), values);
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
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
        if (n != 6 || m != 2 || nonzero_count != 21) {
            return false;
        }
        if (values == nullptr) {
            std::size_t index = 0;
            for (Ipopt::Index row = 0; row < 6; ++row) {
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    rows[index] = row;
                    columns[index] = column;
                    ++index;
                }
            }
            return true;
        }
        try {
            const FlashNlpEvaluation evaluation = evaluate(x);
            for (std::size_t index = 0; index < evaluation.hessian_lower.size(); ++index) {
                values[index] = objective_factor * evaluation.hessian_lower[index];
            }
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
        const Ipopt::Number* z_l,
        const Ipopt::Number* z_u,
        Ipopt::Index m,
        const Ipopt::Number*,
        const Ipopt::Number* lambda,
        Ipopt::Number,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*
    ) override {
        if (n == 6 && m == 2 && x != nullptr && z_l != nullptr && z_u != nullptr
            && lambda != nullptr) {
            std::copy(x, x + n, solution_.begin());
            std::copy(z_l, z_l + n, lower_bound_multipliers_.begin());
            std::copy(z_u, z_u + n, upper_bound_multipliers_.begin());
            std::copy(lambda, lambda + m, equality_multipliers_.begin());
        }
    }

    [[nodiscard]] const std::array<double, 6>& solution() const { return solution_; }
    [[nodiscard]] const std::array<double, 2>& equality_multipliers() const {
        return equality_multipliers_;
    }
    [[nodiscard]] const std::array<double, 6>& lower_bound_multipliers() const {
        return lower_bound_multipliers_;
    }
    [[nodiscard]] const std::array<double, 6>& upper_bound_multipliers() const {
        return upper_bound_multipliers_;
    }
    [[nodiscard]] const std::string& callback_error() const { return callback_error_; }

private:
    [[nodiscard]] FlashNlpEvaluation evaluate(const Ipopt::Number* x) const {
        return evaluate_two_phase_flash_nlp(
            provider_,
            temperature_k_,
            pressure_pa_,
            overall_mole_fractions_,
            {x[0], x[1], x[2], x[3], x[4], x[5]}
        );
    }

    const ProviderContext& provider_;
    double temperature_k_;
    double pressure_pa_;
    std::array<double, 2> overall_mole_fractions_;
    std::array<double, 6> initial_;
    std::array<double, 6> lower_;
    std::array<double, 6> upper_;
    std::array<double, 6> solution_{};
    std::array<double, 2> equality_multipliers_{};
    std::array<double, 6> lower_bound_multipliers_{};
    std::array<double, 6> upper_bound_multipliers_{};
    std::string callback_error_;
};

AttemptResult run_ipopt(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const std::array<double, 2>& overall_mole_fractions,
    const std::array<double, 6>& initial,
    const std::array<double, 6>& lower,
    const std::array<double, 6>& upper
) {
    auto* raw_problem = new FlashTnlp(
        provider,
        temperature_k,
        pressure_pa,
        overall_mole_fractions,
        initial,
        lower,
        upper
    );
    Ipopt::SmartPtr<Ipopt::TNLP> problem = raw_problem;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application = IpoptApplicationFactory();
    application->Options()->SetIntegerValue("print_level", 0);
    application->Options()->SetStringValue("sb", "yes");
    application->Options()->SetIntegerValue("max_iter", 300);
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

    AttemptResult result;
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
    result.equality_multipliers = raw_problem->equality_multipliers();
    result.lower_bound_multipliers = raw_problem->lower_bound_multipliers();
    result.upper_bound_multipliers = raw_problem->upper_bound_multipliers();
    result.callback_error = raw_problem->callback_error();
    const Ipopt::SmartPtr<Ipopt::SolveStatistics> statistics = application->Statistics();
    if (IsValid(statistics)) {
        result.iterations = statistics->IterationCount();
        double dual_infeasibility = 0.0;
        double complementarity = 0.0;
        double kkt_error = 0.0;
        statistics->Infeasibilities(
            dual_infeasibility,
            result.constraint_violation,
            complementarity,
            kkt_error
        );
    }
    return result;
}

void canonicalize_phases(FlashNlpEvaluation& evaluation) {
    const double liquid_density =
        (evaluation.liquid.amounts_mol[0] + evaluation.liquid.amounts_mol[1])
        / evaluation.liquid.volume_m3;
    const double vapor_density =
        (evaluation.vapor.amounts_mol[0] + evaluation.vapor.amounts_mol[1])
        / evaluation.vapor.volume_m3;
    if (liquid_density < vapor_density) {
        std::swap(evaluation.liquid, evaluation.vapor);
    }
}

bool inside_bounds(
    const std::array<double, 6>& variables,
    const std::array<double, 6>& lower,
    const std::array<double, 6>& upper
) {
    for (std::size_t index = 0; index < variables.size(); ++index) {
        if (variables[index] <= lower[index] + kBoundMargin
            || variables[index] >= upper[index] - kBoundMargin) {
            return false;
        }
    }
    return true;
}

double relative_difference(double first, double second) {
    return std::abs(first - second) / std::max({std::abs(first), std::abs(second), 1.0e-12});
}

double phase_amount(const FlashPhaseEvaluation& phase) {
    return phase.amounts_mol[0] + phase.amounts_mol[1];
}

double phase_composition(const FlashPhaseEvaluation& phase, std::size_t component) {
    return phase.amounts_mol[component] / phase_amount(phase);
}

bool physical_acceptance(
    const AttemptResult& attempt,
    FlashNlpEvaluation& evaluation,
    double pressure_pa,
    const std::array<double, 6>& lower,
    const std::array<double, 6>& upper,
    bool require_legacy_chemical_potential_tolerance,
    bool require_density_separation,
    FlashSolveResult& result
) {
    canonicalize_phases(evaluation);
    result.material_balance_max_abs = std::max(
        std::abs(evaluation.constraints[0]),
        std::abs(evaluation.constraints[1])
    );
    result.pressure_stationarity_max_relative = std::max(
        std::abs(evaluation.liquid.provider.pressure_pa - pressure_pa) / pressure_pa,
        std::abs(evaluation.vapor.provider.pressure_pa - pressure_pa) / pressure_pa
    );
    result.chemical_potential_max_abs = 0.0;
    for (std::size_t component = 0; component < 2; ++component) {
        result.chemical_potential_max_abs = std::max(
            result.chemical_potential_max_abs,
            std::abs(
                evaluation.liquid.provider.gradient[component]
                - evaluation.vapor.provider.gradient[component]
            )
        );
    }
    result.kkt_stationarity_max_abs = 0.0;
    for (std::size_t column = 0; column < 6; ++column) {
        double stationarity = evaluation.gradient[column]
            - attempt.lower_bound_multipliers[column]
            + attempt.upper_bound_multipliers[column];
        for (std::size_t row = 0; row < 2; ++row) {
            stationarity += attempt.equality_multipliers[row] * evaluation.jacobian[row * 6 + column];
        }
        result.kkt_stationarity_max_abs = std::max(
            result.kkt_stationarity_max_abs,
            std::abs(stationarity)
        );
    }
    const double liquid_density = phase_amount(evaluation.liquid) / evaluation.liquid.volume_m3;
    const double vapor_density = phase_amount(evaluation.vapor) / evaluation.vapor.volume_m3;
    result.phase_density_distance = (liquid_density - vapor_density) / liquid_density;
    result.equality_multipliers = attempt.equality_multipliers;
    result.lower_bound_multipliers = attempt.lower_bound_multipliers;
    result.upper_bound_multipliers = attempt.upper_bound_multipliers;
    const bool finite = std::isfinite(evaluation.objective)
        && std::isfinite(result.material_balance_max_abs)
        && std::isfinite(result.pressure_stationarity_max_relative)
        && std::isfinite(result.chemical_potential_max_abs)
        && std::isfinite(result.kkt_stationarity_max_abs)
        && std::isfinite(result.phase_density_distance);
    return attempt.solver_converged
        && attempt.callback_error.empty()
        && std::isfinite(attempt.constraint_violation)
        && attempt.constraint_violation <= kSolverConstraintTolerance
        && finite
        && result.material_balance_max_abs <= kSolverConstraintTolerance
        && result.pressure_stationarity_max_relative <= kPhysicalTolerance
        && (!require_legacy_chemical_potential_tolerance
            || result.chemical_potential_max_abs <= kPhysicalTolerance)
        && result.kkt_stationarity_max_abs <= kPhysicalTolerance
        && (!require_density_separation
            || result.phase_density_distance >= kPhaseDistanceTolerance)
        && phase_amount(evaluation.liquid) > 0.0
        && phase_amount(evaluation.vapor) > 0.0
        && inside_bounds(attempt.variables, lower, upper);
}

double confirmation_difference(
    const FlashNlpEvaluation& first,
    const FlashNlpEvaluation& second
) {
    double difference = relative_difference(first.objective, second.objective);
    difference = std::max(
        difference,
        relative_difference(phase_amount(first.liquid), phase_amount(second.liquid))
    );
    difference = std::max(
        difference,
        relative_difference(phase_amount(first.vapor), phase_amount(second.vapor))
    );
    difference = std::max(
        difference,
        relative_difference(first.liquid.volume_m3, second.liquid.volume_m3)
    );
    difference = std::max(
        difference,
        relative_difference(first.vapor.volume_m3, second.vapor.volume_m3)
    );
    for (std::size_t component = 0; component < 2; ++component) {
        difference = std::max(
            difference,
            std::abs(
                phase_composition(first.liquid, component)
                - phase_composition(second.liquid, component)
            )
        );
        difference = std::max(
            difference,
            std::abs(
                phase_composition(first.vapor, component)
                - phase_composition(second.vapor, component)
            )
        );
    }
    return difference;
}

}  // namespace

FlashNlpEvaluation evaluate_two_phase_flash_nlp(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const std::array<double, 2>& overall_mole_fractions,
    const std::array<double, 6>& variables
) {
    require_finite(temperature_k, "temperature");
    require_finite(pressure_pa, "pressure");
    if (temperature_k <= 0.0 || pressure_pa <= 0.0) {
        throw std::invalid_argument("temperature and pressure must be positive");
    }
    for (double value : overall_mole_fractions) {
        require_finite(value, "overall mole fraction");
    }
    for (double value : variables) {
        require_finite(value, "flash NLP variable");
    }

    FlashPhaseEvaluation liquid = evaluate_phase(
        provider,
        temperature_k,
        {variables[0], variables[1]},
        variables[2]
    );
    FlashPhaseEvaluation vapor = evaluate_phase(
        provider,
        temperature_k,
        {variables[3], variables[4]},
        variables[5]
    );
    const double pressure_over_rt = pressure_pa / (kGasConstantJPerMolK * temperature_k);

    FlashNlpEvaluation result;
    result.objective = liquid.provider.value + vapor.provider.value
        + pressure_over_rt * (liquid.volume_m3 + vapor.volume_m3);
    result.constraints = {
        variables[0] + variables[3] - overall_mole_fractions[0],
        variables[1] + variables[4] - overall_mole_fractions[1],
    };
    result.jacobian = {
        1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
    };
    add_phase_derivatives(liquid, pressure_over_rt, 0, result.gradient, result.hessian_lower);
    add_phase_derivatives(vapor, pressure_over_rt, 3, result.gradient, result.hessian_lower);
    result.liquid = std::move(liquid);
    result.vapor = std::move(vapor);
    return result;
}

HeldStageIIIInitialization initialize_held_stage_iii(
    double feed_x_methane,
    const std::vector<HeldStageIIICandidate>& candidates
) {
    HeldStageIIIInitialization result;
    if (!std::isfinite(feed_x_methane) || feed_x_methane <= 0.0
        || feed_x_methane >= 1.0) {
        result.failure_reason = "feed composition must be finite and strictly interior";
        return result;
    }
    if (candidates.size() > 2) {
        result.status = "scope_exceeded";
        result.failure_reason = "more than two Stage III candidates exceed binary v1 scope";
        return result;
    }
    if (candidates.size() != 2) {
        result.failure_reason = "two Stage III candidates are required";
        return result;
    }
    for (const HeldStageIIICandidate& candidate : candidates) {
        if (!std::isfinite(candidate.x_methane)
            || !std::isfinite(candidate.volume_m3_per_mol)
            || candidate.volume_m3_per_mol <= 0.0) {
            result.failure_reason = "candidate composition and volume must be finite and positive";
            return result;
        }
        if (candidate.x_methane <= 1.0e-8 || candidate.x_methane >= 1.0 - 1.0e-8) {
            result.status = "scope_exceeded";
            result.failure_reason = "trace-bound candidate requires a later slice";
            return result;
        }
    }

    const double composition_difference =
        candidates[0].x_methane - candidates[1].x_methane;
    if (composition_difference == 0.0) {
        result.failure_reason = "candidate lever rule is singular";
        return result;
    }
    const double first_fraction =
        (feed_x_methane - candidates[1].x_methane) / composition_difference;
    const double second_fraction = 1.0 - first_fraction;
    if (!std::isfinite(first_fraction) || !std::isfinite(second_fraction)
        || first_fraction <= 0.0 || second_fraction <= 0.0) {
        result.failure_reason = "candidate lever rule does not contain the feed";
        return result;
    }

    result.phase_fractions = {first_fraction, second_fraction};
    result.initial_variables = {
        first_fraction * candidates[0].x_methane,
        first_fraction * (1.0 - candidates[0].x_methane),
        std::log(first_fraction * candidates[0].volume_m3_per_mol),
        second_fraction * candidates[1].x_methane,
        second_fraction * (1.0 - candidates[1].x_methane),
        std::log(second_fraction * candidates[1].volume_m3_per_mol),
    };
    result.lower_bounds = {
        std::max(kAmountLowerMol, result.initial_variables[0] - kStageIIICandidateRadiusMol),
        std::max(kAmountLowerMol, result.initial_variables[1] - kStageIIICandidateRadiusMol),
        std::log(kVolumeLowerM3),
        std::max(kAmountLowerMol, result.initial_variables[3] - kStageIIICandidateRadiusMol),
        std::max(kAmountLowerMol, result.initial_variables[4] - kStageIIICandidateRadiusMol),
        std::log(kVolumeLowerM3),
    };
    result.upper_bounds = {
        result.initial_variables[0] + kStageIIICandidateRadiusMol,
        result.initial_variables[1] + kStageIIICandidateRadiusMol,
        std::log(kVolumeUpperM3),
        result.initial_variables[3] + kStageIIICandidateRadiusMol,
        result.initial_variables[4] + kStageIIICandidateRadiusMol,
        std::log(kVolumeUpperM3),
    };
    const std::array<double, 2> feed{feed_x_methane, 1.0 - feed_x_methane};
    for (std::size_t component = 0; component < 2; ++component) {
        if (result.lower_bounds[component] + result.lower_bounds[component + 3]
                > feed[component]
            || result.upper_bounds[component] + result.upper_bounds[component + 3]
                < feed[component]) {
            result.failure_reason = "component balances are infeasible in candidate neighborhoods";
            return result;
        }
    }
    if (!inside_bounds(
            result.initial_variables,
            result.lower_bounds,
            result.upper_bounds
        )) {
        result.failure_reason = "lever-rule initial state is not interior to Stage III bounds";
        return result;
    }
    result.status = "initialized";
    result.failure_reason.clear();
    result.has_variables = true;
    return result;
}

double held_stage_iii_mu_difference(
    const std::array<double, 2>& first,
    const std::array<double, 2>& second
) {
    double maximum = 0.0;
    for (std::size_t component = 0; component < 2; ++component) {
        const double denominator =
            std::max({std::abs(first[component]), std::abs(second[component]), 1.0});
        maximum = std::max(
            maximum,
            std::abs(first[component] - second[component]) / denominator
        );
    }
    return maximum;
}

HeldStageIIIAcceptanceDecision classify_held_stage_iii(
    const HeldStageIIIAcceptanceEvidence& evidence
) {
    HeldStageIIIAcceptanceDecision result;
    if (!evidence.solver_converged || !evidence.callback_error.empty()
        || !std::isfinite(evidence.solver_constraint_violation)
        || evidence.solver_constraint_violation > kSolverConstraintTolerance) {
        result.failure_reason = "solver_failure";
        return result;
    }
    if (!std::isfinite(evidence.material_balance_max_abs)
        || evidence.material_balance_max_abs > kSolverConstraintTolerance) {
        result.failure_reason = "material_balance_failure";
        return result;
    }
    if (!std::isfinite(evidence.pressure_stationarity_max_relative)
        || evidence.pressure_stationarity_max_relative > kPhysicalTolerance) {
        result.failure_reason = "pressure_failure";
        return result;
    }
    if (!std::isfinite(evidence.kkt_stationarity_max_abs)
        || evidence.kkt_stationarity_max_abs > kPhysicalTolerance) {
        result.failure_reason = "kkt_failure";
        return result;
    }
    if (!evidence.inactive_bounds) {
        result.failure_reason = "active_bound_failure";
        return result;
    }
    if ((!std::isfinite(evidence.composition_distance)
         || evidence.composition_distance < kPhaseDistanceTolerance)
        && (!std::isfinite(evidence.phase_density_distance)
            || evidence.phase_density_distance < kPhaseDistanceTolerance)) {
        result.failure_reason = "phase_collapse";
        return result;
    }
    if (!std::isfinite(evidence.held_gap) || evidence.held_gap < 0.0
        || evidence.held_gap > kHeldFreeEnergyTolerance) {
        result.failure_reason = "held_gap_failure";
        return result;
    }
    if (!std::isfinite(evidence.chemical_potential_max_relative)
        || evidence.chemical_potential_max_relative > kHeldChemicalPotentialTolerance) {
        result.failure_reason = "chemical_potential_failure";
        return result;
    }
    if (!evidence.confirmation_succeeded
        || !std::isfinite(evidence.confirmation_max_difference)
        || evidence.confirmation_max_difference > kConfirmationTolerance) {
        result.failure_reason = "confirmation_failure";
        return result;
    }
    result.status = "accepted";
    return result;
}

HeldStageIIIResult solve_held_stage_iii(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    double feed_x_methane,
    double upper_bound,
    const std::vector<HeldStageIIICandidate>& candidates
) {
    HeldStageIIIResult result;
    result.upper_bound = upper_bound;
    result.initialization = initialize_held_stage_iii(feed_x_methane, candidates);
    if (result.initialization.status != "initialized") {
        result.outcome = result.initialization.status;
        result.failure_reason = result.initialization.failure_reason;
        return result;
    }
    if (!std::isfinite(temperature_k) || !std::isfinite(pressure_pa)
        || temperature_k <= 0.0 || pressure_pa <= 0.0 || !std::isfinite(upper_bound)) {
        result.failure_reason = "Stage III inputs must be finite and positive";
        return result;
    }

    const std::array<double, 2> feed{feed_x_methane, 1.0 - feed_x_methane};
    const auto& lower = result.initialization.lower_bounds;
    const auto& upper = result.initialization.upper_bounds;
    AttemptResult attempt = run_ipopt(
        provider,
        temperature_k,
        pressure_pa,
        feed,
        result.initialization.initial_variables,
        lower,
        upper
    );
    result.local.attempts = 1;
    result.local.solver_lower_bounds = lower;
    result.local.solver_upper_bounds = upper;
    result.local.solver_converged = attempt.solver_converged;
    result.local.solver_status = attempt.solver_status;
    result.local.iterations = attempt.iterations;
    result.local.solver_constraint_violation = attempt.constraint_violation;
    result.local.attempt_log.push_back({
        "refinement",
        attempt.initial_guess,
        attempt.solver_converged,
        attempt.solver_status,
        attempt.iterations,
        attempt.constraint_violation,
        attempt.callback_error,
    });
    if (!attempt.solver_converged || !attempt.callback_error.empty()) {
        result.failure_reason = attempt.callback_error.empty()
            ? "Stage III refinement solve failed"
            : attempt.callback_error;
        return result;
    }

    try {
        result.local.evaluation = evaluate_two_phase_flash_nlp(
            provider,
            temperature_k,
            pressure_pa,
            feed,
            attempt.variables
        );
    } catch (const std::exception& error) {
        result.failure_reason = error.what();
        result.local.attempt_log.back().callback_error = error.what();
        return result;
    }
    result.local.physical_accepted = physical_acceptance(
        attempt,
        result.local.evaluation,
        pressure_pa,
        lower,
        upper,
        false,
        false,
        result.local
    );
    result.local.equality_multipliers = attempt.equality_multipliers;
    result.local.lower_bound_multipliers = attempt.lower_bound_multipliers;
    result.local.upper_bound_multipliers = attempt.upper_bound_multipliers;
    result.inactive_bounds = inside_bounds(attempt.variables, lower, upper);
    result.held_gap = upper_bound - result.local.evaluation.objective;
    const double gap_roundoff = 64.0 * std::numeric_limits<double>::epsilon()
        * std::max({1.0, std::abs(upper_bound), std::abs(result.local.evaluation.objective)});
    if (result.held_gap < 0.0 && std::abs(result.held_gap) <= gap_roundoff) {
        result.held_gap = 0.0;
    }
    result.chemical_potential_max_relative = held_stage_iii_mu_difference(
        {
            result.local.evaluation.liquid.provider.gradient[0],
            result.local.evaluation.liquid.provider.gradient[1],
        },
        {
            result.local.evaluation.vapor.provider.gradient[0],
            result.local.evaluation.vapor.provider.gradient[1],
        }
    );
    result.composition_distance = std::abs(
        phase_composition(result.local.evaluation.liquid, 0)
        - phase_composition(result.local.evaluation.vapor, 0)
    );

    HeldStageIIIAcceptanceEvidence evidence{
        attempt.solver_converged,
        attempt.callback_error,
        attempt.constraint_violation,
        result.local.material_balance_max_abs,
        result.local.pressure_stationarity_max_relative,
        result.local.kkt_stationarity_max_abs,
        result.inactive_bounds,
        result.composition_distance,
        result.local.phase_density_distance,
        result.held_gap,
        result.chemical_potential_max_relative,
        true,
        0.0,
    };
    const HeldStageIIIAcceptanceDecision provisional = classify_held_stage_iii(evidence);
    if (provisional.status != "accepted") {
        result.failure_reason = provisional.failure_reason;
        return result;
    }

    std::array<double, 6> confirmation_seed = attempt.variables;
    for (std::size_t component = 0; component < 2; ++component) {
        const double forward = std::min(
            upper[component] - confirmation_seed[component],
            confirmation_seed[component + 3] - lower[component + 3]
        );
        const double reverse = std::min(
            confirmation_seed[component] - lower[component],
            upper[component + 3] - confirmation_seed[component + 3]
        );
        const double direction = forward >= reverse ? 1.0 : -1.0;
        const double available = std::max(forward, reverse);
        const double transfer = direction * std::min(2.0e-4, 0.5 * available);
        confirmation_seed[component] += transfer;
        confirmation_seed[component + 3] -= transfer;
    }
    confirmation_seed[2] += 0.05;
    confirmation_seed[5] -= 0.05;
    AttemptResult confirmation = run_ipopt(
        provider,
        temperature_k,
        pressure_pa,
        feed,
        confirmation_seed,
        lower,
        upper
    );
    ++result.local.attempts;
    ++result.local.confirmation_solves;
    result.local.attempt_log.push_back({
        "confirmation",
        confirmation.initial_guess,
        confirmation.solver_converged,
        confirmation.solver_status,
        confirmation.iterations,
        confirmation.constraint_violation,
        confirmation.callback_error,
    });
    FlashNlpEvaluation confirmed;
    FlashSolveResult confirmed_local;
    if (confirmation.solver_converged && confirmation.callback_error.empty()) {
        try {
            confirmed = evaluate_two_phase_flash_nlp(
                provider,
                temperature_k,
                pressure_pa,
                feed,
                confirmation.variables
            );
            result.confirmation_succeeded = physical_acceptance(
                confirmation,
                confirmed,
                pressure_pa,
                lower,
                upper,
                false,
                false,
                confirmed_local
            );
        } catch (const std::exception& error) {
            result.local.attempt_log.back().callback_error = error.what();
        }
    }
    result.local.confirmation_max_difference = result.confirmation_succeeded
        ? confirmation_difference(result.local.evaluation, confirmed)
        : std::numeric_limits<double>::infinity();
    evidence.confirmation_succeeded = result.confirmation_succeeded;
    evidence.confirmation_max_difference = result.local.confirmation_max_difference;
    const HeldStageIIIAcceptanceDecision decision = classify_held_stage_iii(evidence);
    if (decision.status != "accepted") {
        result.failure_reason = decision.failure_reason;
        return result;
    }
    result.outcome = "accepted";
    result.local.accepted = true;
    result.local.numerical_converged = true;
    result.failure_reason.clear();
    return result;
}

FlashSolveResult solve_two_phase_flash(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const std::array<double, 2>& overall_mole_fractions
) {
    FlashSolveResult last_result;
    last_result.failure_reason = "no deterministic seed produced an accepted local two-phase result";
    const std::array<double, 6> lower = solver_lower_bounds();
    const std::array<double, 6> upper = solver_upper_bounds(overall_mole_fractions);
    last_result.solver_lower_bounds = lower;
    last_result.solver_upper_bounds = upper;
    const double delta = std::min(
        0.20,
        0.45 * std::min(overall_mole_fractions[0], 1.0 - overall_mole_fractions[0])
    );
    for (double sign : {1.0, -1.0}) {
        const double liquid_x0 = overall_mole_fractions[0] - sign * delta;
        const double vapor_x0 = overall_mole_fractions[0] + sign * delta;
        const std::array<double, 6> initial{
            0.5 * liquid_x0,
            0.5 * (1.0 - liquid_x0),
            std::log(0.5 / 20'000.0),
            0.5 * vapor_x0,
            0.5 * (1.0 - vapor_x0),
            std::log(0.5 * kGasConstantJPerMolK * temperature_k / pressure_pa),
        };
        AttemptResult attempt = run_ipopt(
            provider,
            temperature_k,
            pressure_pa,
            overall_mole_fractions,
            initial,
            lower,
            upper
        );
        ++last_result.attempts;
        last_result.attempt_log.push_back({
            "search",
            attempt.initial_guess,
            attempt.solver_converged,
            attempt.solver_status,
            attempt.iterations,
            attempt.constraint_violation,
            attempt.callback_error,
        });
        last_result.solver_converged = attempt.solver_converged;
        last_result.solver_status = attempt.solver_status;
        last_result.iterations = attempt.iterations;
        last_result.solver_constraint_violation = attempt.constraint_violation;
        if (!attempt.callback_error.empty()) {
            last_result.failure_reason = attempt.callback_error;
        }
        if (!attempt.solver_converged) {
            continue;
        }
        FlashNlpEvaluation candidate;
        try {
            candidate = evaluate_two_phase_flash_nlp(
                provider,
                temperature_k,
                pressure_pa,
                overall_mole_fractions,
                attempt.variables
            );
        } catch (const std::exception& error) {
            last_result.failure_reason = error.what();
            continue;
        }
        FlashSolveResult candidate_result = last_result;
        candidate_result.evaluation = candidate;
        candidate_result.physical_accepted = physical_acceptance(
            attempt,
            candidate_result.evaluation,
            pressure_pa,
            lower,
            upper,
            true,
            true,
            candidate_result
        );
        if (!candidate_result.physical_accepted) {
            candidate_result.failure_reason = "Ipopt solution failed local physical acceptance";
            last_result = candidate_result;
            continue;
        }

        std::array<double, 6> confirmation_seed = attempt.variables;
        double transfer_sign = 1.0;
        if (confirmation_seed[3] - 0.02 * overall_mole_fractions[0] <= kAmountLowerMol
            || confirmation_seed[4] - 0.02 * overall_mole_fractions[1] <= kAmountLowerMol) {
            transfer_sign = -1.0;
        }
        for (std::size_t component = 0; component < 2; ++component) {
            const double transfer = transfer_sign * 0.02 * overall_mole_fractions[component];
            confirmation_seed[component] += transfer;
            confirmation_seed[component + 3] -= transfer;
        }
        confirmation_seed[2] += 0.05;
        confirmation_seed[5] -= 0.05;
        AttemptResult confirmation = run_ipopt(
            provider,
            temperature_k,
            pressure_pa,
            overall_mole_fractions,
            confirmation_seed,
            lower,
            upper
        );
        ++candidate_result.attempts;
        ++candidate_result.confirmation_solves;
        candidate_result.attempt_log.push_back({
            "confirmation",
            confirmation.initial_guess,
            confirmation.solver_converged,
            confirmation.solver_status,
            confirmation.iterations,
            confirmation.constraint_violation,
            confirmation.callback_error,
        });
        if (!confirmation.solver_converged || !confirmation.callback_error.empty()) {
            candidate_result.failure_reason = "perturbed confirmation solve failed";
            last_result = candidate_result;
            continue;
        }
        FlashNlpEvaluation confirmed;
        try {
            confirmed = evaluate_two_phase_flash_nlp(
                provider,
                temperature_k,
                pressure_pa,
                overall_mole_fractions,
                confirmation.variables
            );
        } catch (const std::exception& error) {
            candidate_result.failure_reason = error.what();
            last_result = candidate_result;
            continue;
        }
        FlashSolveResult confirmed_result;
        if (!physical_acceptance(
                confirmation,
                confirmed,
                pressure_pa,
                lower,
                upper,
                true,
                true,
                confirmed_result
            )) {
            candidate_result.failure_reason = "perturbed confirmation failed local acceptance";
            last_result = candidate_result;
            continue;
        }
        candidate_result.confirmation_max_difference = confirmation_difference(
            candidate_result.evaluation,
            confirmed
        );
        candidate_result.numerical_converged =
            candidate_result.confirmation_max_difference <= kConfirmationTolerance;
        if (!candidate_result.numerical_converged) {
            candidate_result.failure_reason = "perturbed confirmation disagreed with the local result";
            last_result = candidate_result;
            continue;
        }
        candidate_result.accepted = true;
        candidate_result.failure_reason.clear();
        return candidate_result;
    }
    return last_result;
}

}  // namespace epcsaft_equilibrium
