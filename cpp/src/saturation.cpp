#include "saturation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>
#include <coin/IpTNLP.hpp>

namespace epcsaft_equilibrium {

namespace {

constexpr double kGasConstantJPerMolK = 8.31446261815324;
constexpr double kVaporDensityLower = 1.0e-8;
constexpr double kVaporDensityUpper = 5.0e3;
constexpr double kLiquidDensityLower = 5.001e3;
constexpr double kPressureLowerPa = 1.0e-12;
constexpr double kPressureUpperPa = 1.0e8;
constexpr double kConstraintTolerance = 1.0e-8;
constexpr double kPhaseDistanceTolerance = 1.0e-3;

void require_finite(double value, const char* name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
}

struct PhaseDerivatives {
    PhaseEvaluation phase;
    double pressure_q;
    double pressure_qq;
    double chemical_potential_q;
    double chemical_potential_qq;
};

PhaseDerivatives evaluate_phase_derivatives(
    const ProviderContext& provider,
    double temperature_k,
    double log_density
) {
    const double volume_m3 = std::exp(-log_density);
    require_finite(volume_m3, "phase volume");
    if (volume_m3 <= 0.0) {
        throw std::invalid_argument("phase volume must be positive");
    }
    PhaseEvaluation phase = provider.evaluate(temperature_k, 1.0, volume_m3);
    const double rt = kGasConstantJPerMolK * temperature_k;
    const double h_nv = phase.provider.hessian[1];
    const double h_vv = phase.provider.hessian[3];
    const double t_nvv = phase.provider.third[3];
    const double t_vvv = phase.provider.third[7];
    return {
        std::move(phase),
        rt * h_vv * volume_m3,
        -rt * (t_vvv * volume_m3 * volume_m3 + h_vv * volume_m3),
        -h_nv * volume_m3,
        h_nv * volume_m3 + t_nvv * volume_m3 * volume_m3,
    };
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

struct AttemptResult {
    bool solver_converged = false;
    std::string solver_status;
    std::string callback_error;
    int iterations = 0;
    double solver_constraint_violation = std::numeric_limits<double>::infinity();
    std::array<double, 3> initial_guess{};
    std::array<double, 3> variables{};
};

std::array<double, 3> solver_lower_bounds() {
    return {
        std::log(kVaporDensityLower),
        std::log(kLiquidDensityLower),
        std::log(kPressureLowerPa),
    };
}

std::array<double, 3> solver_upper_bounds(double liquid_density_upper_mol_m3) {
    return {
        std::log(kVaporDensityUpper),
        std::log(liquid_density_upper_mol_m3),
        std::log(kPressureUpperPa),
    };
}

class SaturationTnlp final : public Ipopt::TNLP {
public:
    SaturationTnlp(
        const ProviderContext& provider,
        double temperature_k,
        std::array<double, 3> initial,
        double liquid_density_upper_mol_m3
    )
        : provider_(provider),
          temperature_k_(temperature_k),
          initial_(initial),
          lower_(solver_lower_bounds()),
          upper_(solver_upper_bounds(liquid_density_upper_mol_m3)) {}

    bool get_nlp_info(
        Ipopt::Index& n,
        Ipopt::Index& m,
        Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag,
        IndexStyleEnum& index_style
    ) override {
        n = 3;
        m = 3;
        nnz_jac_g = 9;
        nnz_h_lag = 6;
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
        if (n != 3 || m != 3) {
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
        if (n != 3 || m != 3 || !init_x || init_z || init_lambda) {
            return false;
        }
        std::copy(initial_.begin(), initial_.end(), x);
        return true;
    }

    bool eval_f(
        Ipopt::Index n,
        const Ipopt::Number*,
        bool,
        Ipopt::Number& objective
    ) override {
        if (n != 3) {
            return false;
        }
        objective = 0.0;
        return true;
    }

    bool eval_grad_f(
        Ipopt::Index n,
        const Ipopt::Number*,
        bool,
        Ipopt::Number* gradient
    ) override {
        if (n != 3) {
            return false;
        }
        std::fill(gradient, gradient + n, 0.0);
        return true;
    }

    bool eval_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Index m,
        Ipopt::Number* constraints
    ) override {
        if (n != 3 || m != 3) {
            return false;
        }
        try {
            const NlpEvaluation evaluation = evaluate(x, {0.0, 0.0, 0.0});
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
        if (n != 3 || m != 3 || nonzero_count != 9) {
            return false;
        }
        if (values == nullptr) {
            for (Ipopt::Index row = 0; row < 3; ++row) {
                for (Ipopt::Index column = 0; column < 3; ++column) {
                    rows[row * 3 + column] = row;
                    columns[row * 3 + column] = column;
                }
            }
            return true;
        }
        try {
            const NlpEvaluation evaluation = evaluate(x, {0.0, 0.0, 0.0});
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
        Ipopt::Number,
        Ipopt::Index m,
        const Ipopt::Number* lambda,
        bool,
        Ipopt::Index nonzero_count,
        Ipopt::Index* rows,
        Ipopt::Index* columns,
        Ipopt::Number* values
    ) override {
        if (n != 3 || m != 3 || nonzero_count != 6) {
            return false;
        }
        if (values == nullptr) {
            constexpr std::array<Ipopt::Index, 6> kRows{0, 1, 1, 2, 2, 2};
            constexpr std::array<Ipopt::Index, 6> kColumns{0, 0, 1, 0, 1, 2};
            std::copy(kRows.begin(), kRows.end(), rows);
            std::copy(kColumns.begin(), kColumns.end(), columns);
            return true;
        }
        try {
            const NlpEvaluation evaluation = evaluate(
                x,
                {lambda[0], lambda[1], lambda[2]}
            );
            std::copy(
                evaluation.lagrangian_hessian_lower.begin(),
                evaluation.lagrangian_hessian_lower.end(),
                values
            );
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
        Ipopt::Index,
        const Ipopt::Number*,
        const Ipopt::Number*,
        Ipopt::Number,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*
    ) override {
        if (n == 3 && x != nullptr) {
            std::copy(x, x + n, solution_.begin());
        }
    }

    [[nodiscard]] const std::array<double, 3>& solution() const {
        return solution_;
    }

    [[nodiscard]] const std::string& callback_error() const {
        return callback_error_;
    }

private:
    [[nodiscard]] NlpEvaluation evaluate(
        const Ipopt::Number* x,
        const std::array<double, 3>& multipliers
    ) const {
        return evaluate_saturation_nlp(
            provider_,
            temperature_k_,
            {x[0], x[1], x[2]},
            multipliers
        );
    }

    const ProviderContext& provider_;
    double temperature_k_;
    std::array<double, 3> initial_;
    std::array<double, 3> lower_;
    std::array<double, 3> upper_;
    std::array<double, 3> solution_{};
    std::string callback_error_;
};

AttemptResult run_ipopt(
    const ProviderContext& provider,
    double temperature_k,
    const std::array<double, 3>& initial,
    double liquid_density_upper_mol_m3
) {
    auto* raw_problem = new SaturationTnlp(
        provider,
        temperature_k,
        initial,
        liquid_density_upper_mol_m3
    );
    Ipopt::SmartPtr<Ipopt::TNLP> problem = raw_problem;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application = IpoptApplicationFactory();
    application->Options()->SetIntegerValue("print_level", 0);
    application->Options()->SetStringValue("sb", "yes");
    application->Options()->SetIntegerValue("max_iter", 200);
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

    const Ipopt::ApplicationReturnStatus initialize_status = application->Initialize();
    AttemptResult result;
    result.initial_guess = initial;
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
        double dual_infeasibility = 0.0;
        double complementarity = 0.0;
        double kkt_error = 0.0;
        statistics->Infeasibilities(
            dual_infeasibility,
            result.solver_constraint_violation,
            complementarity,
            kkt_error
        );
    }
    return result;
}

double relative_difference(double first, double second) {
    return std::abs(first - second) / std::max({std::abs(first), std::abs(second), 1.0e-300});
}

bool inside_bounds(
    const std::array<double, 3>& variables,
    double liquid_density_upper_mol_m3
) {
    constexpr double kMargin = 1.0e-7;
    const std::array<double, 3> lower = solver_lower_bounds();
    const std::array<double, 3> upper = solver_upper_bounds(liquid_density_upper_mol_m3);
    for (std::size_t index = 0; index < variables.size(); ++index) {
        if (variables[index] <= lower[index] + kMargin
            || variables[index] >= upper[index] - kMargin) {
            return false;
        }
    }
    return true;
}

bool physical_acceptance(
    const AttemptResult& attempt,
    const NlpEvaluation& evaluation,
    double liquid_density_upper_mol_m3,
    SaturationSolveResult& result
) {
    const double vapor_density = 1.0 / evaluation.vapor.volume_m3;
    const double liquid_density = 1.0 / evaluation.liquid.volume_m3;
    result.pressure_relative_residual = std::max(
        std::abs(evaluation.constraints[0]),
        std::abs(evaluation.constraints[1])
    );
    result.chemical_potential_absolute_residual = std::abs(evaluation.constraints[2]);
    result.phase_density_distance = (liquid_density - vapor_density) / liquid_density;
    const bool finite = std::isfinite(evaluation.saturation_pressure_pa)
        && std::isfinite(vapor_density)
        && std::isfinite(liquid_density)
        && std::isfinite(result.pressure_relative_residual)
        && std::isfinite(result.chemical_potential_absolute_residual)
        && std::isfinite(result.phase_density_distance);
    const bool locally_stable = evaluation.jacobian[0] > 0.0 && evaluation.jacobian[4] > 0.0;
    return attempt.solver_converged
        && finite
        && evaluation.saturation_pressure_pa > 0.0
        && vapor_density > 0.0
        && liquid_density > vapor_density
        && result.pressure_relative_residual <= kConstraintTolerance
        && result.chemical_potential_absolute_residual <= kConstraintTolerance
        && result.phase_density_distance > kPhaseDistanceTolerance
        && locally_stable
        && inside_bounds(attempt.variables, liquid_density_upper_mol_m3);
}

}  // namespace

ProviderContext::ProviderContext(const epcsaft_native_sdk_v1& sdk, std::string fingerprint)
    : sdk_(sdk), fingerprint_(std::move(fingerprint)) {
    if (fingerprint_.empty()) {
        throw std::invalid_argument("expected provider fingerprint must not be empty");
    }
}

PhaseEvaluation ProviderContext::evaluate(
    double temperature_k,
    double amount_mol,
    double volume_m3
) const {
    epcsaft_phase_block_result_v1 phase{};
    phase.struct_size = sizeof(phase);
    const int status = sdk_.evaluate_pure_phase(
        sdk_.model_context,
        temperature_k,
        amount_mol,
        volume_m3,
        &phase
    );
    if (status != phase.status) {
        throw std::runtime_error("provider phase evaluation returned inconsistent status values");
    }
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1) {
        throw std::domain_error(
            "provider phase evaluation failed: " + std::string(phase.error)
        );
    }
    if (phase.struct_size != sizeof(epcsaft_phase_block_result_v1)) {
        throw std::invalid_argument("provider phase result struct size mismatch");
    }
    if (std::string(phase.parameter_fingerprint) != fingerprint_) {
        throw std::invalid_argument(
            "provider phase result fingerprint does not match the requested model"
        );
    }
    return {amount_mol, volume_m3, phase};
}

const std::string& ProviderContext::fingerprint() const {
    return fingerprint_;
}

NlpEvaluation evaluate_saturation_nlp(
    const ProviderContext& provider,
    double temperature_k,
    const std::array<double, 3>& variables,
    const std::array<double, 3>& multipliers
) {
    require_finite(temperature_k, "temperature");
    if (temperature_k <= 0.0) {
        throw std::invalid_argument("temperature must be positive");
    }
    for (double value : variables) {
        require_finite(value, "NLP variable");
    }
    for (double value : multipliers) {
        require_finite(value, "constraint multiplier");
    }

    PhaseDerivatives vapor = evaluate_phase_derivatives(provider, temperature_k, variables[0]);
    PhaseDerivatives liquid = evaluate_phase_derivatives(provider, temperature_k, variables[1]);
    const double saturation_pressure_pa = std::exp(variables[2]);
    require_finite(saturation_pressure_pa, "saturation pressure");
    if (saturation_pressure_pa <= 0.0) {
        throw std::invalid_argument("saturation pressure must be positive");
    }

    const double vapor_pressure_ratio =
        vapor.phase.provider.pressure_pa / saturation_pressure_pa;
    const double liquid_pressure_ratio =
        liquid.phase.provider.pressure_pa / saturation_pressure_pa;
    const double vapor_pressure_q_ratio = vapor.pressure_q / saturation_pressure_pa;
    const double liquid_pressure_q_ratio = liquid.pressure_q / saturation_pressure_pa;

    NlpEvaluation result{
        0.0,
        {0.0, 0.0, 0.0},
        {
            vapor_pressure_ratio - 1.0,
            liquid_pressure_ratio - 1.0,
            vapor.phase.provider.chemical_potential_over_rt
                - liquid.phase.provider.chemical_potential_over_rt,
        },
        {
            vapor_pressure_q_ratio, 0.0, -vapor_pressure_ratio,
            0.0, liquid_pressure_q_ratio, -liquid_pressure_ratio,
            vapor.chemical_potential_q, -liquid.chemical_potential_q, 0.0,
        },
        {
            multipliers[0] * vapor.pressure_qq / saturation_pressure_pa
                + multipliers[2] * vapor.chemical_potential_qq,
            0.0,
            multipliers[1] * liquid.pressure_qq / saturation_pressure_pa
                - multipliers[2] * liquid.chemical_potential_qq,
            -multipliers[0] * vapor_pressure_q_ratio,
            -multipliers[1] * liquid_pressure_q_ratio,
            multipliers[0] * vapor_pressure_ratio
                + multipliers[1] * liquid_pressure_ratio,
        },
        std::move(vapor.phase),
        std::move(liquid.phase),
        saturation_pressure_pa,
    };
    return result;
}

SaturationSolveResult solve_local_saturation(
    const ProviderContext& provider,
    double temperature_k,
    double liquid_density_upper_mol_m3
) {
    require_finite(liquid_density_upper_mol_m3, "liquid density upper bound");
    if (liquid_density_upper_mol_m3 <= kLiquidDensityLower) {
        throw std::invalid_argument("liquid density upper bound must exceed 5001 mol m^-3");
    }

    const std::array<double, 7> vapor_seeds{1.0e-5, 1.0e-3, 0.1, 10.0, 100.0, 1000.0, 3000.0};
    const std::array<double, 3> liquid_seeds{18'000.0, 12'000.0, 24'000.0};
    SaturationSolveResult last_result;
    last_result.failure_reason = "no Ipopt seed produced a certified local boundary";
    last_result.solver_lower_bounds = solver_lower_bounds();
    last_result.solver_upper_bounds = solver_upper_bounds(liquid_density_upper_mol_m3);

    for (double vapor_density : vapor_seeds) {
        for (double liquid_density : liquid_seeds) {
            const double pressure_seed = kGasConstantJPerMolK * temperature_k * vapor_density;
            const std::array<double, 3> initial{
                std::log(vapor_density),
                std::log(liquid_density),
                std::log(pressure_seed),
            };
            AttemptResult attempt;
            attempt.initial_guess = initial;
            try {
                const NlpEvaluation seed_evaluation = evaluate_saturation_nlp(
                    provider,
                    temperature_k,
                    initial,
                    {0.0, 0.0, 0.0}
                );
                if (seed_evaluation.jacobian[0] <= 0.0
                    || seed_evaluation.jacobian[4] <= 0.0) {
                    attempt.solver_status = "seed_rejected_local_instability";
                } else {
                    attempt = run_ipopt(
                        provider,
                        temperature_k,
                        initial,
                        liquid_density_upper_mol_m3
                    );
                }
            } catch (const std::exception& error) {
                attempt.solver_status = "seed_evaluation_failed";
                attempt.callback_error = error.what();
            }
            ++last_result.attempts;
            last_result.attempt_log.push_back({
                "search",
                attempt.initial_guess,
                attempt.solver_converged,
                attempt.solver_status,
                attempt.iterations,
                attempt.solver_constraint_violation,
                attempt.callback_error,
            });
            last_result.solver_status = attempt.solver_status;
            last_result.solver_converged = attempt.solver_converged;
            last_result.iterations = attempt.iterations;
            last_result.solver_constraint_violation = attempt.solver_constraint_violation;
            if (!attempt.callback_error.empty()) {
                last_result.failure_reason = attempt.callback_error;
            }
            if (!attempt.solver_converged) {
                continue;
            }

            NlpEvaluation candidate;
            try {
                candidate = evaluate_saturation_nlp(
                    provider,
                    temperature_k,
                    attempt.variables,
                    {0.0, 0.0, 0.0}
                );
            } catch (const std::exception& error) {
                last_result.failure_reason = error.what();
                continue;
            }
            SaturationSolveResult candidate_result = last_result;
            candidate_result.evaluation = candidate;
            candidate_result.physical_accepted = physical_acceptance(
                    attempt,
                    candidate,
                    liquid_density_upper_mol_m3,
                    candidate_result
                );
            if (!candidate_result.physical_accepted) {
                candidate_result.failure_reason = "Ipopt solution failed local physical acceptance";
                last_result = candidate_result;
                continue;
            }

            const std::array<std::array<double, 3>, 2> perturbations{{
                {0.05, -0.05, 0.02},
                {-0.05, 0.05, -0.02},
            }};
            bool confirmations_agree = true;
            double maximum_difference = 0.0;
            for (const auto& perturbation : perturbations) {
                std::array<double, 3> confirmation_seed{};
                for (std::size_t index = 0; index < confirmation_seed.size(); ++index) {
                    confirmation_seed[index] = attempt.variables[index] + perturbation[index];
                }
                AttemptResult confirmation = run_ipopt(
                    provider,
                    temperature_k,
                    confirmation_seed,
                    liquid_density_upper_mol_m3
                );
                ++candidate_result.attempts;
                candidate_result.attempt_log.push_back({
                    "confirmation",
                    confirmation.initial_guess,
                    confirmation.solver_converged,
                    confirmation.solver_status,
                    confirmation.iterations,
                    confirmation.solver_constraint_violation,
                    confirmation.callback_error,
                });
                if (!confirmation.solver_converged) {
                    confirmations_agree = false;
                    break;
                }
                NlpEvaluation confirmed;
                try {
                    confirmed = evaluate_saturation_nlp(
                        provider,
                        temperature_k,
                        confirmation.variables,
                        {0.0, 0.0, 0.0}
                    );
                } catch (const std::exception&) {
                    confirmations_agree = false;
                    break;
                }
                SaturationSolveResult confirmed_result;
                if (!physical_acceptance(
                        confirmation,
                        confirmed,
                        liquid_density_upper_mol_m3,
                        confirmed_result
                    )) {
                    confirmations_agree = false;
                    break;
                }
                maximum_difference = std::max({
                    maximum_difference,
                    relative_difference(
                        candidate.saturation_pressure_pa,
                        confirmed.saturation_pressure_pa
                    ),
                    relative_difference(candidate.vapor.volume_m3, confirmed.vapor.volume_m3),
                    relative_difference(candidate.liquid.volume_m3, confirmed.liquid.volume_m3),
                });
                ++candidate_result.confirmation_solves;
            }
            candidate_result.confirmation_max_relative_difference = maximum_difference;
            candidate_result.numerical_converged =
                confirmations_agree && maximum_difference <= 1.0e-7;
            if (!confirmations_agree || maximum_difference > 1.0e-7) {
                candidate_result.failure_reason = "perturbed local solves did not confirm the boundary";
                last_result = candidate_result;
                continue;
            }
            candidate_result.accepted = true;
            candidate_result.failure_reason.clear();
            return candidate_result;
        }
    }
    return last_result;
}

}  // namespace epcsaft_equilibrium
