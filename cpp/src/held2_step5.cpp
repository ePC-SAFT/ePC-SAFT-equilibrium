#include "held2_step5.hpp"
#include "held2_tolerances.hpp"

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>

namespace epcsaft_equilibrium {
namespace {

struct Step5Assessment {
    bool qualified = false;
    double gap = 0.0;
    std::string reason;
};

Step5Assessment assess_step5(
    double upper_bound,
    double local_value,
    bool local_state_certified
) {
    Step5Assessment result;
    result.gap = upper_bound - local_value;
    if (!std::isfinite(upper_bound) || !std::isfinite(local_value)) {
        result.reason = "nonfinite_lower_upper_comparison";
    } else if (!local_state_certified) {
        result.reason = "local_state_not_certified";
    } else if (local_value <= upper_bound) {
        result.qualified = true;
        result.reason = "lower_not_above_upper";
    } else if (audit_held2_tolerance(kHeld2Step6Gap, result.gap).passed) {
        result.qualified = true;
        result.reason = "lower_equal_within_step6_gate";
    } else {
        result.reason = "certified_local_above_upper";
    }
    return result;
}

struct Step5LocalRun {
    bool converged = false;
    std::string status = "not_run";
    std::vector<double> variables;
    int iterations = 0;
};

class Step5Tnlp final : public Ipopt::TNLP {
public:
    Step5Tnlp(
        Held2StateEvaluator objective,
        std::vector<Held2PolytopeConstraint> constraints,
        std::vector<double> initial,
        std::vector<double> lower,
        std::vector<double> upper
    )
        : objective_(std::move(objective)),
          constraints_(std::move(constraints)),
          initial_(std::move(initial)),
          lower_(std::move(lower)),
          upper_(std::move(upper)) {}

    bool get_nlp_info(
        Ipopt::Index& n,
        Ipopt::Index& m,
        Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag,
        IndexStyleEnum& index_style
    ) override {
        n = static_cast<Ipopt::Index>(initial_.size());
        m = static_cast<Ipopt::Index>(constraints_.size());
        nnz_jac_g = m * (n - 1);
        nnz_h_lag = n * (n + 1) / 2;
        index_style = C_STYLE;
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
        if (n != static_cast<Ipopt::Index>(initial_.size())
            || m != static_cast<Ipopt::Index>(constraints_.size())) {
            return false;
        }
        std::copy(lower_.begin(), lower_.end(), x_l);
        std::copy(upper_.begin(), upper_.end(), x_u);
        for (Ipopt::Index row = 0; row < m; ++row) {
            g_l[row] = -1.0e20;
            g_u[row] = constraints_[static_cast<std::size_t>(row)]
                .upper_bound;
        }
        return true;
    }

    bool get_starting_point(
        Ipopt::Index n,
        bool init_x,
        Ipopt::Number* x,
        bool init_z,
        Ipopt::Number*,
        Ipopt::Number*,
        Ipopt::Index,
        bool init_lambda,
        Ipopt::Number*
    ) override {
        if (n != static_cast<Ipopt::Index>(initial_.size())
            || !init_x || init_z || init_lambda) {
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
        try {
            objective = evaluate(n, x).objective;
            return true;
        } catch (const std::exception& error) {
            static_cast<void>(error);
            return false;
        }
    }

    bool eval_grad_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number* gradient
    ) override {
        try {
            const Held2StateEvaluation& state = evaluate(n, x);
            std::copy(state.gradient.begin(), state.gradient.end(), gradient);
            return true;
        } catch (const std::exception& error) {
            static_cast<void>(error);
            return false;
        }
    }

    bool eval_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Index m,
        Ipopt::Number* values
    ) override {
        if (n != static_cast<Ipopt::Index>(initial_.size())
            || m != static_cast<Ipopt::Index>(constraints_.size())) {
            return false;
        }
        for (Ipopt::Index row = 0; row < m; ++row) {
            const auto& coefficients =
                constraints_[static_cast<std::size_t>(row)].coefficients;
            values[row] = std::inner_product(
                coefficients.begin(), coefficients.end(), x, 0.0
            );
        }
        return true;
    }

    bool eval_jac_g(
        Ipopt::Index n,
        const Ipopt::Number*,
        bool,
        Ipopt::Index m,
        Ipopt::Index nonzero_count,
        Ipopt::Index* rows,
        Ipopt::Index* columns,
        Ipopt::Number* values
    ) override {
        if (nonzero_count != m * (n - 1)) {
            return false;
        }
        Ipopt::Index position = 0;
        for (Ipopt::Index row = 0; row < m; ++row) {
            for (Ipopt::Index column = 0; column < n - 1; ++column) {
                if (values == nullptr) {
                    rows[position] = row;
                    columns[position] = column;
                } else {
                    values[position] =
                        constraints_[static_cast<std::size_t>(row)]
                            .coefficients[
                                static_cast<std::size_t>(column)
                            ];
                }
                ++position;
            }
        }
        return true;
    }

    bool eval_h(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number objective_factor,
        Ipopt::Index,
        const Ipopt::Number*,
        bool,
        Ipopt::Index nonzero_count,
        Ipopt::Index* rows,
        Ipopt::Index* columns,
        Ipopt::Number* values
    ) override {
        if (nonzero_count != n * (n + 1) / 2) {
            return false;
        }
        Ipopt::Index position = 0;
        if (values == nullptr) {
            for (Ipopt::Index row = 0; row < n; ++row) {
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    rows[position] = row;
                    columns[position++] = column;
                }
            }
            return true;
        }
        try {
            const Held2StateEvaluation& state = evaluate(n, x);
            for (Ipopt::Index row = 0; row < n; ++row) {
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    values[position++] = objective_factor * state.hessian[
                        static_cast<std::size_t>(row * n + column)
                    ];
                }
            }
            return true;
        } catch (const std::exception& error) {
            static_cast<void>(error);
            return false;
        }
    }

    bool intermediate_callback(
        Ipopt::AlgorithmMode,
        Ipopt::Index iteration,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Index,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*
    ) override {
        iterations_ = std::max(iterations_, static_cast<int>(iteration) + 1);
        return true;
    }

    void finalize_solution(
        Ipopt::SolverReturn status,
        Ipopt::Index n,
        const Ipopt::Number* x,
        const Ipopt::Number* z_lower,
        const Ipopt::Number* z_upper,
        Ipopt::Index m,
        const Ipopt::Number*,
        const Ipopt::Number* lambda,
        Ipopt::Number,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*
    ) override {
        static_cast<void>(z_lower);
        static_cast<void>(z_upper);
        static_cast<void>(m);
        static_cast<void>(lambda);
        variables_.assign(x, x + n);
        converged_ = status == Ipopt::SUCCESS
            || status == Ipopt::STOP_AT_ACCEPTABLE_POINT;
        status_ = status == Ipopt::SUCCESS
            ? "success"
            : status == Ipopt::STOP_AT_ACCEPTABLE_POINT
                ? "acceptable_point"
                : "failed_" + std::to_string(static_cast<int>(status));
    }

    [[nodiscard]] Step5LocalRun result() const {
        return {
            converged_,
            status_,
            variables_,
            iterations_,
        };
    }

private:
    const Held2StateEvaluation& evaluate(
        Ipopt::Index n,
        const Ipopt::Number* x
    ) const {
        if (n != static_cast<Ipopt::Index>(initial_.size())) {
            throw std::invalid_argument(
                "HELD2 Step-5 variable count changed"
            );
        }
        const std::vector<double> variables(x, x + n);
        if (!cached_state_ || variables != cached_variables_) {
            cached_variables_ = variables;
            cached_state_ = objective_(
                std::vector<double>(variables.begin(), variables.end() - 1),
                variables.back()
            );
        }
        return *cached_state_;
    }

    Held2StateEvaluator objective_;
    std::vector<Held2PolytopeConstraint> constraints_;
    std::vector<double> initial_;
    std::vector<double> lower_;
    std::vector<double> upper_;
    mutable std::vector<double> cached_variables_;
    mutable std::optional<Held2StateEvaluation> cached_state_;
    bool converged_ = false;
    std::string status_ = "not_run";
    std::vector<double> variables_;
    int iterations_ = 0;
};

Step5LocalRun solve_step5_local(
    Held2StateEvaluator objective,
    const std::vector<Held2PolytopeConstraint>& constraints,
    const std::vector<double>& initial,
    const std::vector<double>& lower,
    const std::vector<double>& upper
) {
    Ipopt::SmartPtr<Step5Tnlp> problem = new Step5Tnlp(
        std::move(objective),
        constraints,
        initial,
        lower,
        upper
    );
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application =
        IpoptApplicationFactory();
    application->Options()->SetStringValue("option_file_name", "");
    application->Options()->SetIntegerValue("print_level", 0);
    application->Options()->SetStringValue("sb", "yes");
    application->Options()->SetIntegerValue("max_iter", 300);
    application->Options()->SetNumericValue(
        "tol", kHeld2IpoptTarget.atol
    );
    application->Options()->SetStringValue(
        "jacobian_approximation", "exact"
    );
    application->Options()->SetStringValue(
        "hessian_approximation", "exact"
    );
    application->Options()->SetNumericValue(
        "constr_viol_tol", kHeld2IpoptConstraint.atol
    );
    application->Options()->SetNumericValue("bound_relax_factor", 0.0);
    if (application->Initialize() != Ipopt::Solve_Succeeded) {
        return {};
    }
    application->OptimizeTNLP(problem);
    return problem->result();
}

}  // namespace

Held2Step5Result run_held2_step5(
    const Held2Step1Result& step1,
    const Held2Step4Result& step4,
    Held2PersistentState& state,
    const Held2StateEvaluator& evaluator,
    const Held2PackingFractionEvaluator& packing_fraction,
    const Held2ResourceProfile& resources,
    Held2ProgressObserver* observer
) {
    Held2Step5Result result;
    result.timing.invocation_count = 1;
    if (step4.status != "complete" || !step1.coordinates
        || !step1.volume_bounds || !evaluator || !packing_fraction
        || resources.step5_start_cap < 1) {
        result.reason = "invalid_step5_input";
        return result;
    }
    const Held2Coordinates& coordinates = *step1.coordinates;
    const std::size_t dimension = state.feed.size();
    std::vector<double> lower = coordinates.independent_lower_bounds;
    std::vector<double> upper = coordinates.independent_upper_bounds;
    const std::size_t begin =
        static_cast<std::size_t>(state.next_start_ordinal);
    std::mt19937_64 random(0);
    random.discard(begin * (dimension + 1));
    std::uniform_real_distribution<double> unit;
    std::vector<std::vector<double>> starts(
        static_cast<std::size_t>(resources.step5_start_cap),
        std::vector<double>(dimension + 1)
    );
    for (std::vector<double>& start : starts) {
        for (double& value : start) {
            value = unit(random);
        }
    }
    const Held2StateEvaluator counted_evaluator =
        [&evaluator, &result](
            const std::vector<double>& composition,
            double log_volume
        ) {
            ++result.timing.provider_evaluations;
            return evaluator(composition, log_volume);
        };
    const Held2VolumeBoundsEvaluator counted_volume_bounds =
        [&step1, &result](const std::vector<double>& composition) {
            ++result.timing.provider_evaluations;
            return (*step1.volume_bounds)(composition);
        };
    if (!state.step5_volume_bounds) {
        std::array<double, 2> bounds{
            std::numeric_limits<double>::infinity(), 0.0
        };
        try {
            for (int index = 0; index < resources.step5_start_cap; ++index) {
                const std::vector<double>& point =
                    starts[static_cast<std::size_t>(index)];
                const std::array<double, 2> local = counted_volume_bounds(
                    held2_map_unit_cube_to_independent_fractions(
                        coordinates,
                        std::vector<double>(
                            point.begin(),
                            point.begin()
                                + static_cast<std::ptrdiff_t>(dimension)
                        ),
                        step1.total_ion_mole_fraction_max
                    )
                );
                bounds[0] = std::min(bounds[0], local[0]);
                bounds[1] = std::max(bounds[1], local[1]);
            }
        } catch (...) {
            result.reason = "volume_bounds_failed";
            return result;
        }
        state.step5_volume_bounds = bounds;
    }
    const std::array<double, 2>& solve_volume_bounds =
        *state.step5_volume_bounds;
    double best = std::numeric_limits<double>::infinity();
    Held2MPoint best_point;
    for (int attempt = 0;
         attempt < resources.step5_start_cap;
        ++attempt) {
        const std::vector<double>& point =
            starts[static_cast<std::size_t>(attempt)];
        const std::uint64_t ordinal = state.next_start_ordinal++;
        ++result.starts_consumed;
        const std::vector<double> start =
            held2_map_unit_cube_to_independent_fractions(
                coordinates,
                std::vector<double>(
                    point.begin(),
                    point.begin()
                        + static_cast<std::ptrdiff_t>(dimension)
                ),
                step1.total_ion_mole_fraction_max
            );
        Held2LocalCertificate certificate;
        certificate.start_ordinal = ordinal;
        const double initial_log_volume =
            std::log(solve_volume_bounds[0])
            + point.back()
                * std::log(
                    solve_volume_bounds[1] / solve_volume_bounds[0]
                );
        std::vector<double> initial = start;
        initial.push_back(initial_log_volume);
        std::vector<double> solver_lower = lower;
        std::vector<double> solver_upper = upper;
        solver_lower.push_back(std::log(solve_volume_bounds[0]));
        solver_upper.push_back(std::log(solve_volume_bounds[1]));
        const Held2StateEvaluator objective = [
            &counted_evaluator,
            feed = state.feed,
            multipliers = state.multipliers
        ](const std::vector<double>& composition, double log_volume) {
            Held2StateEvaluation evaluated =
                counted_evaluator(composition, log_volume);
            for (std::size_t coordinate = 0;
                 coordinate < composition.size();
                 ++coordinate) {
                evaluated.objective += multipliers[coordinate]
                    * (feed[coordinate] - composition[coordinate]);
                evaluated.gradient[coordinate] -= multipliers[coordinate];
            }
            return evaluated;
        };
        const Step5LocalRun run = solve_step5_local(
            objective,
            coordinates.polytope_constraints,
            initial,
            solver_lower,
            solver_upper
        );
        ++result.timing.optimizer_solves;
        result.timing.optimizer_iterations +=
            static_cast<std::uint64_t>(run.iterations);
        certificate.solver_status = run.status;
        if (!run.converged || run.variables.size() != dimension + 1) {
            result.attempts.push_back(certificate);
            continue;
        }
        const std::vector<double> independent(
            run.variables.begin(), run.variables.end() - 1
        );
        const Held2StateEvaluation terminal =
            counted_evaluator(independent, run.variables.back());
        double value = terminal.objective;
        for (std::size_t coordinate = 0;
             coordinate < dimension;
             ++coordinate) {
            value += state.multipliers[coordinate]
                * (state.feed[coordinate] - independent[coordinate]);
        }
        certificate.finite_and_in_domain =
            std::isfinite(value)
            && terminal.volume >= solve_volume_bounds[0]
            && terminal.volume <= solve_volume_bounds[1];
        certificate.accepted = certificate.finite_and_in_domain;
        result.attempts.push_back(certificate);
        Held2ProgressEvent progress;
        progress.kind = Held2ProgressKind::Certificate;
        progress.stage = "STEP 5 LOCAL";
        progress.major_iteration = state.major_iteration;
        progress.attempt = static_cast<int>(ordinal);
        progress.status = certificate.accepted ? "accepted" : "rejected";
        progress.reason = run.status;
        progress.pressure_residual =
            terminal.pressure_stationarity_relative;
        progress.objective = value;
        progress.upper_bound = state.upper_bound;
        progress.gap = state.upper_bound - value;
        observe_held2(observer, progress);
        if (!assess_step5(
                state.upper_bound, value, certificate.accepted
            ).qualified || value >= best) {
            continue;
        }
        best = value;
        ++result.timing.provider_evaluations;
        const double best_packing_fraction =
            packing_fraction(independent, terminal.volume);
        best_point = {
            static_cast<std::uint64_t>(state.M.size()),
            independent,
            terminal.volume,
            best_packing_fraction,
            terminal.objective,
            terminal.gradient,
            "step5_local",
        };
        break;
    }
    if (!std::isfinite(best)) {
        result.reason = "step5_start_budget_exhausted";
        return result;
    }
    state.lower_value = best;
    result.lower_value = best;
    result.terminal = best_point;
    state.M.push_back(best_point);
    result.status = "complete";
    result.reason = "step5_complete";
    result.timing.terminal_status = result.status;
    result.timing.terminal_reason = result.reason;
    result.timing.next_step = 6;
    return result;
}

}  // namespace epcsaft_equilibrium
