#include "held2_step5.hpp"
#include "held2_tolerances.hpp"

#include <nlopt.hpp>
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

std::string nlopt_version_string() {
    int major = 0;
    int minor = 0;
    int bugfix = 0;
    nlopt_version(&major, &minor, &bugfix);
    return std::to_string(major) + "." + std::to_string(minor) + "."
        + std::to_string(bugfix);
}

double maximum_abs_difference(
    const std::vector<double>& left,
    const std::vector<double>& right
) {
    if (left.size() != right.size()) {
        return std::numeric_limits<double>::infinity();
    }
    double maximum = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        maximum = std::max(maximum, std::abs(left[index] - right[index]));
    }
    return maximum;
}

bool same_composition(
    const std::vector<double>& left,
    const std::vector<double>& right
) {
    return maximum_abs_difference(left, right)
        <= kHeld2BasinDuplicateComposition.atol;
}

bool same_physical_start(
    const Held2StageIIPhysicalStart& left,
    const Held2StageIIPhysicalStart& right
) {
    return same_composition(
               left.independent_modified_fractions,
               right.independent_modified_fractions
           )
        && std::abs(left.log_volume - right.log_volume)
            <= kHeld2BasinDuplicateLogVolume.atol;
}

Held2StageIIBasinEvaluation evaluate_fail_closed(
    const Held2StageIIBasinEvaluator& evaluator,
    const std::vector<double>& independent
) {
    try {
        return evaluator(independent);
    } catch (const std::exception& error) {
        Held2StageIIBasinEvaluation failed;
        failed.independent_modified_fractions = independent;
        failed.failure_reason = std::string("envelope_evaluator_exception: ")
            + error.what();
        return failed;
    }
}

void retain_evaluation(
    Held2StageIIBasinExplorationResult& result,
    Held2StageIIBasinEvaluation evaluation,
    const std::string& source
) {
    if (!evaluation.certified
        || !std::isfinite(evaluation.reduced_lower_value)) {
        ++result.failed_evaluation_count;
        result.evaluations.push_back(std::move(evaluation));
        return;
    }
    ++result.completed_evaluation_count;
    double minimum_stable_root_objective =
        std::numeric_limits<double>::infinity();
    for (const Held2PressureRoot& root : evaluation.pressure_envelope.roots) {
        if (root.mechanical_class == "strict_stable" && !root.boundary) {
            minimum_stable_root_objective = std::min(
                minimum_stable_root_objective,
                root.objective
            );
        }
    }
    const double composition_offset = evaluation.reduced_lower_value
        - minimum_stable_root_objective;
    int stable_branch_index = 0;
    for (const Held2PressureRoot& root : evaluation.pressure_envelope.roots) {
        if (root.mechanical_class != "strict_stable" || root.boundary) {
            continue;
        }
        Held2StageIIPhysicalStart start;
        start.independent_modified_fractions =
            evaluation.independent_modified_fractions;
        start.stable_branch_index = stable_branch_index++;
        start.log_volume = root.log_volume;
        start.volume = root.volume;
        start.reduced_lower_value = root.objective + composition_offset;
        start.source = source;
        start.root_origin = root.origin;
        start.root_completeness =
            evaluation.pressure_envelope.root_completeness;
        const bool duplicate = std::any_of(
            result.representatives.begin(),
            result.representatives.end(),
            [&start](const Held2StageIIPhysicalStart& known) {
                return same_physical_start(start, known);
            }
        );
        if (duplicate) {
            ++result.duplicate_start_count;
        } else {
            result.representatives.push_back(std::move(start));
        }
    }
    result.evaluations.push_back(std::move(evaluation));
}

struct DirectContext {
    const Held2Coordinates* coordinates = nullptr;
    const Held2StageIIBasinEvaluator* evaluator = nullptr;
    Held2StageIIBasinExplorationResult* result = nullptr;
    nlopt::opt* optimizer = nullptr;
    double total_ion_mole_fraction_max =
        std::numeric_limits<double>::quiet_NaN();
    bool stop_requested = false;
};

double direct_objective(
    const std::vector<double>& cube,
    std::vector<double>& gradient,
    void* opaque
) {
    auto& context = *static_cast<DirectContext*>(opaque);
    if (context.stop_requested) {
        return 0.0;
    }
    if (!gradient.empty()) {
        throw std::invalid_argument(
            "HELD2 Stage-II DIRECT-L requested an unexpected gradient"
        );
    }
    const std::vector<double> independent =
        held2_map_unit_cube_to_independent_fractions(
            *context.coordinates,
            cube,
            context.total_ion_mole_fraction_max
        );
    Held2StageIIBasinEvaluation evaluation = evaluate_fail_closed(
        *context.evaluator,
        independent
    );
    const double objective = evaluation.reduced_lower_value;
    retain_evaluation(*context.result, std::move(evaluation), "direct_l");
    if (context.result->failed_evaluation_count != 0) {
        context.stop_requested = true;
        context.result->termination_reason =
            "required_envelope_evaluation_failed";
        context.optimizer->force_stop();
        return 0.0;
    }
    return objective;
}

struct SobolParameters {
    int degree;
    std::uint32_t coefficient;
    std::array<std::uint32_t, 5> initial;
};

constexpr std::array<SobolParameters, 9> kSobolParameters = {{
    {1, 0, {1, 0, 0, 0, 0}},
    {2, 1, {1, 3, 0, 0, 0}},
    {3, 1, {1, 3, 1, 0, 0}},
    {3, 2, {1, 1, 1, 0, 0}},
    {4, 1, {1, 3, 5, 13, 0}},
    {4, 4, {1, 1, 5, 5, 0}},
    {5, 2, {1, 3, 3, 9, 7}},
    {5, 4, {1, 1, 5, 11, 27}},
    {5, 7, {1, 1, 7, 13, 3}},
}};

std::array<std::uint32_t, 32> sobol_directions(std::size_t dimension) {
    std::array<std::uint32_t, 32> directions{};
    if (dimension == 0) {
        for (std::size_t bit = 0; bit < directions.size(); ++bit) {
            directions[bit] = std::uint32_t{1} << (31U - bit);
        }
        return directions;
    }
    if (dimension > kSobolParameters.size()) {
        throw std::invalid_argument(
            "HELD2 Sobol explorer supports at most ten composition coordinates"
        );
    }
    const SobolParameters& parameters = kSobolParameters[dimension - 1];
    for (int bit = 0; bit < parameters.degree; ++bit) {
        directions[static_cast<std::size_t>(bit)] =
            parameters.initial[static_cast<std::size_t>(bit)]
            << (31 - bit);
    }
    for (int bit = parameters.degree; bit < 32; ++bit) {
        std::uint32_t value = directions[static_cast<std::size_t>(
            bit - parameters.degree
        )];
        value ^= value >> parameters.degree;
        for (int offset = 1; offset < parameters.degree; ++offset) {
            const int coefficient_bit = parameters.degree - 1 - offset;
            if ((parameters.coefficient >> coefficient_bit) & 1U) {
                value ^= directions[static_cast<std::size_t>(bit - offset)];
            }
        }
        directions[static_cast<std::size_t>(bit)] = value;
    }
    return directions;
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
        if (!assess_held2_stage_ii_step5(
                state.upper_bound, value, certificate.accepted
            ).qualified || value >= best) {
            continue;
        }
        best = value;
        best_point = {
            static_cast<std::uint64_t>(state.M.size()),
            independent,
            terminal.volume,
            packing_fraction(independent, terminal.volume),
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

std::vector<std::vector<double>> held2_sobol_points(
    std::size_t dimension,
    int count
) {
    if (dimension == 0 || count < 0 || dimension > 10) {
        throw std::invalid_argument("HELD2 Sobol policy is invalid");
    }
    std::vector<std::array<std::uint32_t, 32>> directions;
    directions.reserve(dimension);
    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
        directions.push_back(sobol_directions(coordinate));
    }
    std::vector<std::vector<double>> points;
    points.reserve(static_cast<std::size_t>(count));
    constexpr double scale = 1.0 / 4294967296.0;
    for (int index = 1; index <= count; ++index) {
        const std::uint32_t gray = static_cast<std::uint32_t>(index)
            ^ (static_cast<std::uint32_t>(index) >> 1U);
        std::vector<double> point(dimension, 0.0);
        for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
            std::uint32_t value = 0;
            for (std::size_t bit = 0; bit < 32; ++bit) {
                if ((gray >> bit) & 1U) {
                    value ^= directions[coordinate][bit];
                }
            }
            point[coordinate] = static_cast<double>(value) * scale;
        }
        points.push_back(std::move(point));
    }
    return points;
}

Held2StageIIBasinExplorationResult explore_held2_stage_ii_basins(
    const Held2Coordinates& coordinates,
    const std::vector<Held2StageIIBasinSeed>& seeds,
    int sobol_count,
    bool use_direct_escalation,
    int direct_evaluation_budget,
    double total_ion_mole_fraction_max,
    const Held2StageIIBasinEvaluator& evaluator
) {
    if (coordinates.independent_indices.empty() || sobol_count < 0
        || direct_evaluation_budget < 0
        || (use_direct_escalation && direct_evaluation_budget == 0)
        || (!std::isnan(total_ion_mole_fraction_max)
            && (!std::isfinite(total_ion_mole_fraction_max)
                || total_ion_mole_fraction_max <= 0.0
                || total_ion_mole_fraction_max > 1.0))) {
        throw std::invalid_argument("HELD2 Stage-II basin policy is invalid");
    }
    Held2StageIIBasinExplorationResult result;
    result.direct_solver_version = nlopt_version_string();
    result.declared_sobol_count = sobol_count;
    result.declared_direct_budget = use_direct_escalation
        ? direct_evaluation_budget
        : 0;

    std::vector<Held2StageIIBasinSeed> unique_seeds;
    for (const Held2StageIIBasinSeed& seed : seeds) {
        const bool duplicate = std::any_of(
            unique_seeds.begin(),
            unique_seeds.end(),
            [&seed](const Held2StageIIBasinSeed& known) {
                return same_composition(
                    seed.independent_modified_fractions,
                    known.independent_modified_fractions
                );
            }
        );
        if (duplicate) {
            ++result.duplicate_start_count;
        } else {
            unique_seeds.push_back(seed);
        }
    }
    for (const Held2StageIIBasinSeed& seed : unique_seeds) {
        retain_evaluation(
            result,
            evaluate_fail_closed(
                evaluator,
                seed.independent_modified_fractions
            ),
            seed.source
        );
        if (result.failed_evaluation_count != 0) {
            result.termination_reason = "required_envelope_evaluation_failed";
            return result;
        }
    }
    for (const std::vector<double>& cube : held2_sobol_points(
             coordinates.independent_indices.size(),
             sobol_count
         )) {
        const std::vector<double> independent =
            held2_map_unit_cube_to_independent_fractions(
                coordinates,
                cube,
                total_ion_mole_fraction_max
            );
        retain_evaluation(
            result,
            evaluate_fail_closed(evaluator, independent),
            "sobol"
        );
        if (result.failed_evaluation_count != 0) {
            result.termination_reason = "required_envelope_evaluation_failed";
            return result;
        }
    }

    if (use_direct_escalation) {
        result.direct_escalation_used = true;
        nlopt::opt optimizer(
            nlopt::GN_DIRECT_L,
            coordinates.independent_indices.size()
        );
        DirectContext context{
            &coordinates,
            &evaluator,
            &result,
            &optimizer,
            total_ion_mole_fraction_max,
            false,
        };
        optimizer.set_lower_bounds(std::vector<double>(
            coordinates.independent_indices.size(),
            0.0
        ));
        optimizer.set_upper_bounds(std::vector<double>(
            coordinates.independent_indices.size(),
            1.0
        ));
        optimizer.set_maxeval(direct_evaluation_budget);
        optimizer.set_min_objective(direct_objective, &context);
        std::vector<double> initial(
            coordinates.independent_indices.size(),
            0.5
        );
        double minimum = std::numeric_limits<double>::infinity();
        try {
            const nlopt::result status = optimizer.optimize(initial, minimum);
            if (status != nlopt::MAXEVAL_REACHED) {
                result.termination_reason = "unexpected_direct_termination";
                return result;
            }
        } catch (const nlopt::forced_stop&) {
            if (result.failed_evaluation_count != 0) {
                return result;
            }
            result.termination_reason = "unexpected_direct_forced_stop";
            return result;
        } catch (const std::exception& error) {
            result.termination_reason = std::string("direct_solver_failure: ")
                + error.what();
            return result;
        }
    }

    if (result.representatives.empty()) {
        result.termination_reason = "no_physical_basin_representatives";
        return result;
    }
    std::sort(
        result.representatives.begin(),
        result.representatives.end(),
        [](const Held2StageIIPhysicalStart& left,
           const Held2StageIIPhysicalStart& right) {
            if (left.reduced_lower_value != right.reduced_lower_value) {
                return left.reduced_lower_value < right.reduced_lower_value;
            }
            if (left.independent_modified_fractions
                != right.independent_modified_fractions) {
                return left.independent_modified_fractions
                    < right.independent_modified_fractions;
            }
            return left.log_volume < right.log_volume;
        }
    );
    result.outcome = "representatives_found";
    result.termination_reason = use_direct_escalation
        ? "declared_direct_budget_exhausted"
        : "deterministic_exploration_completed";
    return result;
}

}  // namespace epcsaft_equilibrium
