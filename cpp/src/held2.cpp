#include "held2.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>

namespace epcsaft_equilibrium {
namespace {

constexpr double kCoordinateTolerance = 1.0e-12;
// The approved manufactured oracle uses a strict 1e-8 direct certificate.
constexpr double kCertificateTolerance = 1.0e-8;
// Shared with the existing Stage-III projected-KKT certificate.
constexpr double kStageIIKktTolerance = 1.0e-7;
constexpr double kStageITpdThreshold = -1.0e-8;
constexpr double kStageIIStep6EpsilonB = 1.0e-8;
constexpr double kStageIIStep6EpsilonLambda = 1.0e-8;
constexpr double kStageIIStep6EpsilonEta = 1.0e-3;
constexpr double kStageIIStep6EpsilonX = 1.0e-3;
constexpr double kStageIIStep6QStationarityTolerance = 1.0e-8;
constexpr int kStageIIpoptIterations = 300;
constexpr unsigned int kStageISeed = 2025;

void require_finite_vector(const std::vector<double>& values, const char* name) {
    if (!std::all_of(values.begin(), values.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
}

double maximum_abs_difference(
    const std::vector<double>& left,
    const std::vector<double>& right
) {
    if (left.size() != right.size()) {
        throw std::invalid_argument("vectors must have matching sizes");
    }
    double result = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        result = std::max(result, std::abs(left[index] - right[index]));
    }
    return result;
}

double charge_residual(
    const std::vector<double>& charges,
    const std::vector<double>& fractions
) {
    double result = 0.0;
    for (std::size_t index = 0; index < charges.size(); ++index) {
        result += charges[index] * fractions[index];
    }
    return result;
}

double manufactured_helmholtz(double composition, double molar_volume) {
    const double shifted = composition - 0.5;
    constexpr double inner_squared = 0.15 * 0.15;
    constexpr double outer_squared = 0.30 * 0.30;
    const double composition_energy = 1.0e4
        * (std::pow(shifted, 6) / 6.0
           - (inner_squared + outer_squared) * std::pow(shifted, 4) / 4.0
           + inner_squared * outer_squared * shifted * shifted / 2.0);
    const double volume_delta = molar_volume - 1.2;
    return composition_energy + 0.5 * 5.0 * volume_delta * volume_delta;
}

std::array<double, 2> manufactured_gibbs_gradient(
    double composition,
    double molar_volume
) {
    const double shifted = composition - 0.5;
    constexpr double inner_squared = 0.15 * 0.15;
    constexpr double outer_squared = 0.30 * 0.30;
    return {
        1.0e4
            * (std::pow(shifted, 5)
               - (inner_squared + outer_squared) * std::pow(shifted, 3)
               + inner_squared * outer_squared * shifted),
        5.0 * (molar_volume - 1.2) + 1.0,
    };
}

double manufactured_gibbs(double composition, double molar_volume) {
    return manufactured_helmholtz(composition, molar_volume) + molar_volume;
}

double manufactured_composition_hessian(double composition) {
    const double shifted = composition - 0.5;
    constexpr double inner_squared = 0.15 * 0.15;
    constexpr double outer_squared = 0.30 * 0.30;
    return 1.0e4
        * (5.0 * std::pow(shifted, 4)
           - 3.0 * (inner_squared + outer_squared) * shifted * shifted
           + inner_squared * outer_squared);
}

std::array<double, 2> manufactured_modified_potentials(
    double composition,
    double molar_volume
);

Held2StateEvaluation evaluate_manufactured_state_impl(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent,
    double log_volume
) {
    if (independent.size() != 1 || !std::isfinite(log_volume)) {
        throw std::invalid_argument("manufactured HELD2 phase state has the wrong size");
    }
    Held2StateEvaluation state;
    state.physical_amounts = held2_lift_independent_fractions(coordinates, independent);
    state.modified_fractions = held2_transform_physical_fractions(
        coordinates,
        state.physical_amounts
    );
    state.volume = std::exp(log_volume);
    if (!std::isfinite(state.volume) || state.volume <= 0.0) {
        throw std::invalid_argument("manufactured HELD2 phase volume must be positive");
    }
    const double composition = independent.front();
    const auto phase_gradient = manufactured_gibbs_gradient(composition, state.volume);
    state.objective = manufactured_gibbs(composition, state.volume);
    state.gradient = {phase_gradient[0], state.volume * phase_gradient[1]};
    state.hessian = {
        manufactured_composition_hessian(composition),
        0.0,
        0.0,
        5.0 * state.volume * state.volume + state.volume * phase_gradient[1],
    };
    const auto potentials = manufactured_modified_potentials(composition, state.volume);
    state.modified_potentials.assign(potentials.begin(), potentials.end());
    state.pressure_stationarity_relative = phase_gradient[1];
    state.pressure_stationarity_derivative_log_volume = 5.0 * state.volume;
    return state;
}

void configure_held2_ipopt(const Ipopt::SmartPtr<Ipopt::IpoptApplication>& application) {
    application->Options()->SetStringValue("option_file_name", "");
    application->Options()->SetIntegerValue("print_level", 0);
    application->Options()->SetStringValue("sb", "yes");
    application->Options()->SetIntegerValue("max_iter", kStageIIpoptIterations);
    application->Options()->SetNumericValue("tol", 1.0e-10);
    application->Options()->SetNumericValue("acceptable_tol", 1.0e-9);
    application->Options()->SetIntegerValue("acceptable_iter", 0);
    application->Options()->SetStringValue("jacobian_approximation", "exact");
    application->Options()->SetStringValue("hessian_approximation", "exact");
    application->Options()->SetStringValue("nlp_scaling_method", "none");
    application->Options()->SetNumericValue("bound_relax_factor", 0.0);
    application->Options()->SetStringValue("honor_original_bounds", "yes");
    application->Options()->SetStringValue("check_derivatives_for_naninf", "yes");
}

struct Held2SearchRun {
    bool solver_converged = false;
    std::string callback_error;
    std::vector<double> variables;
    int solver_status = 999;
    int iterations = -1;
    double final_step_norm = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> lower_bound_multipliers;
    std::vector<double> upper_bound_multipliers;
    std::vector<double> constraint_multipliers;
    std::vector<double> coordinate_jacobian;
};

using Held2SearchDomainPredicate = std::function<bool(
    const std::vector<double>&,
    double
)>;

std::array<double, 2> strict_interior_log_bounds(
    const std::array<double, 2>& raw_bounds
) {
    if (!std::isfinite(raw_bounds[0]) || !std::isfinite(raw_bounds[1])
        || raw_bounds[0] <= 0.0 || raw_bounds[1] <= raw_bounds[0]) {
        throw std::invalid_argument(
            "HELD2 physical bounds must be finite, positive, and ordered"
        );
    }

    double lower = std::log(raw_bounds[0]);
    double upper = std::log(raw_bounds[1]);
    while (!(std::exp(lower) > raw_bounds[0])) {
        const double next = std::nextafter(
            lower, std::numeric_limits<double>::infinity()
        );
        if (!(next > lower) || !(next < upper)) {
            throw std::invalid_argument(
                "HELD2 physical bounds have no representable logarithmic interior"
            );
        }
        lower = next;
    }
    while (!(std::exp(upper) < raw_bounds[1])) {
        const double next = std::nextafter(
            upper, -std::numeric_limits<double>::infinity()
        );
        if (!(next < upper) || !(next > lower)) {
            throw std::invalid_argument(
                "HELD2 physical bounds have no representable logarithmic interior"
            );
        }
        upper = next;
    }
    return {lower, upper};
}

struct Held2StageIISimplexChart {
    std::vector<double> physical_coordinates;
    std::vector<double> chart_coordinates;
    std::vector<double> jacobian;
    std::vector<double> component_hessians;
    bool singular = false;
};

struct Held2PhysicalKkt {
    double stationarity_inf_norm = std::numeric_limits<double>::infinity();
    double complementarity = std::numeric_limits<double>::infinity();
    double reconstruction_inf_norm = std::numeric_limits<double>::infinity();
    bool dual_signs_valid = false;
};

[[nodiscard]] Held2StageIISimplexChart evaluate_held2_stage_ii_simplex_chart(
    const std::vector<double>& independent_lower_bounds,
    const std::vector<double>& independent_upper_bounds,
    double composition_sum_upper,
    const std::vector<double>& coordinates,
    bool inverse
);

[[nodiscard]] Held2StateEvaluation transform_held2_stage_ii_simplex_state(
    const Held2StateEvaluation& physical_state,
    const Held2StageIISimplexChart& chart
);

class Held2SearchTnlp final : public Ipopt::TNLP {
public:
    Held2SearchTnlp(
        Held2StateEvaluator evaluator,
        std::vector<double> reference_variables,
        Held2StateEvaluation reference,
        bool use_tpd,
        std::vector<double> initial,
        std::vector<double> lower,
        std::vector<double> upper,
        double composition_sum_upper,
        std::vector<double> source_domain_coefficients,
        double source_domain_upper,
        Held2SearchDomainPredicate trial_domain_predicate
    )
        : evaluator_(std::move(evaluator)),
          reference_variables_(std::move(reference_variables)),
          reference_(std::move(reference)),
          use_tpd_(use_tpd),
          initial_(std::move(initial)),
          lower_(std::move(lower)),
          upper_(std::move(upper)),
          composition_sum_upper_(composition_sum_upper),
          source_domain_coefficients_(std::move(source_domain_coefficients)),
          source_domain_upper_(source_domain_upper),
          trial_domain_predicate_(std::move(trial_domain_predicate)) {}

    bool get_nlp_info(
        Ipopt::Index& n,
        Ipopt::Index& m,
        Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag,
        IndexStyleEnum& index_style
    ) override {
        n = static_cast<Ipopt::Index>(initial_.size());
        m = constraint_count();
        nnz_jac_g = has_composition_constraint() ? n - 1 : 0;
        if (has_source_domain_constraint()) {
            nnz_jac_g += static_cast<Ipopt::Index>(std::count_if(
                source_domain_coefficients_.begin(),
                source_domain_coefficients_.end(),
                [](double value) { return value != 0.0; }
            ));
        }
        nnz_h_lag = n * (n + 1) / 2;
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
        if (n != static_cast<Ipopt::Index>(initial_.size()) || m != constraint_count()) {
            return false;
        }
        std::copy(lower_.begin(), lower_.end(), x_l);
        std::copy(upper_.begin(), upper_.end(), x_u);
        if (has_composition_constraint()) {
            g_l[composition_constraint_row()] = -2.0e19;
            g_u[composition_constraint_row()] = composition_sum_upper_;
        }
        if (has_source_domain_constraint()) {
            g_l[source_domain_constraint_row()] = -2.0e19;
            g_u[source_domain_constraint_row()] = source_domain_upper_;
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
        Ipopt::Index m,
        bool init_lambda,
        Ipopt::Number*
    ) override {
        if (n != static_cast<Ipopt::Index>(initial_.size()) || m != constraint_count() || !init_x
            || init_z || init_lambda) {
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
        if (!evaluation_domain_valid(n, x)) {
            return false;
        }
        try {
            objective = evaluate(n, x).objective;
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
        if (!evaluation_domain_valid(n, x)) {
            return false;
        }
        try {
            const Held2StateEvaluation evaluation = evaluate(n, x);
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
        if (m != constraint_count()) {
            return false;
        }
        if (has_composition_constraint()) {
            constraints[composition_constraint_row()] =
                std::accumulate(x, x + n - 1, 0.0);
        }
        if (has_source_domain_constraint()) {
            constraints[source_domain_constraint_row()] = 0.0;
            for (Ipopt::Index index = 0; index < n - 1; ++index) {
                constraints[source_domain_constraint_row()] +=
                    source_domain_coefficients_[static_cast<std::size_t>(index)]
                    * x[index];
            }
        }
        return true;
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
        Ipopt::Index expected_nonzeros = has_composition_constraint() ? n - 1 : 0;
        if (has_source_domain_constraint()) {
            expected_nonzeros += static_cast<Ipopt::Index>(std::count_if(
                source_domain_coefficients_.begin(),
                source_domain_coefficients_.end(),
                [](double value) { return value != 0.0; }
            ));
        }
        if (m != constraint_count() || nonzero_count != expected_nonzeros) {
            return false;
        }
        Ipopt::Index position = 0;
        if (has_composition_constraint()) {
            for (Ipopt::Index index = 0; index < n - 1; ++index) {
                if (values == nullptr) {
                    rows[position] = composition_constraint_row();
                    columns[position] = index;
                } else {
                    values[position] = 1.0;
                }
                ++position;
            }
        }
        if (has_source_domain_constraint()) {
            for (Ipopt::Index index = 0; index < n - 1; ++index) {
                const double coefficient =
                    source_domain_coefficients_[static_cast<std::size_t>(index)];
                if (coefficient == 0.0) {
                    continue;
                }
                if (values == nullptr) {
                    rows[position] = source_domain_constraint_row();
                    columns[position] = index;
                } else {
                    values[position] = coefficient;
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
        Ipopt::Index m,
        const Ipopt::Number* constraint_multipliers,
        bool,
        Ipopt::Index nonzero_count,
        Ipopt::Index* rows,
        Ipopt::Index* columns,
        Ipopt::Number* values
    ) override {
        if (m != constraint_count() || nonzero_count != n * (n + 1) / 2) {
            return false;
        }
        Ipopt::Index position = 0;
        if (values == nullptr) {
            for (Ipopt::Index row = 0; row < n; ++row) {
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    rows[position] = row;
                    columns[position] = column;
                    ++position;
                }
            }
            return true;
        }
        if (!evaluation_domain_valid(n, x)) {
            return false;
        }
        try {
            const Held2StateEvaluation evaluation = evaluate(n, x);
            for (Ipopt::Index row = 0; row < n; ++row) {
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    values[position] = objective_factor
                        * evaluation.hessian[static_cast<std::size_t>(row * n + column)];
                    ++position;
                }
            }
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
    }

    void finalize_solution(
        Ipopt::SolverReturn status,
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
        if (m == constraint_count() && x != nullptr) {
            variables_.assign(x, x + n);
            if (z_l != nullptr && z_u != nullptr) {
                lower_bound_multipliers_.assign(z_l, z_l + n);
                upper_bound_multipliers_.assign(z_u, z_u + n);
            }
            if (lambda != nullptr) {
                constraint_multipliers_.assign(lambda, lambda + m);
            }
        }
        solver_converged_ = status == Ipopt::SUCCESS
            || status == Ipopt::STOP_AT_ACCEPTABLE_POINT;
    }

    bool intermediate_callback(
        Ipopt::AlgorithmMode,
        Ipopt::Index iteration,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number step_norm,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Index,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*
    ) override {
        iterations_ = static_cast<int>(iteration);
        final_step_norm_ = step_norm;
        return true;
    }

    [[nodiscard]] Held2SearchRun result() const {
        return {
            solver_converged_,
            callback_error_,
            variables_,
            999,
            iterations_,
            final_step_norm_,
            lower_bound_multipliers_,
            upper_bound_multipliers_,
            constraint_multipliers_,
            {},
        };
    }

private:
    [[nodiscard]] Ipopt::Index constraint_count() const {
        return static_cast<Ipopt::Index>(
            has_composition_constraint() + has_source_domain_constraint()
        );
    }

    [[nodiscard]] bool has_composition_constraint() const {
        return std::isfinite(composition_sum_upper_);
    }

    [[nodiscard]] Ipopt::Index composition_constraint_row() const {
        return 0;
    }

    [[nodiscard]] bool has_source_domain_constraint() const {
        return std::isfinite(source_domain_upper_);
    }

    [[nodiscard]] Ipopt::Index source_domain_constraint_row() const {
        return static_cast<Ipopt::Index>(has_composition_constraint());
    }

    [[nodiscard]] bool source_domain_valid(
        Ipopt::Index n,
        const Ipopt::Number* x
    ) const {
        if (!has_source_domain_constraint()) {
            return true;
        }
        if (source_domain_coefficients_.size()
            != static_cast<std::size_t>(n - 1)) {
            return false;
        }
        double total = 0.0;
        for (Ipopt::Index index = 0; index < n - 1; ++index) {
            total += source_domain_coefficients_[static_cast<std::size_t>(index)]
                * x[index];
        }
        return std::isfinite(total) && total <= source_domain_upper_;
    }

    [[nodiscard]] bool evaluation_domain_valid(
        Ipopt::Index n,
        const Ipopt::Number* x
    ) {
        if (!source_domain_valid(n, x)) {
            return false;
        }
        if (!trial_domain_predicate_) {
            return true;
        }
        try {
            const std::vector<double> independent(x, x + n - 1);
            return trial_domain_predicate_(independent, x[n - 1]);
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
    }

    [[nodiscard]] Held2StateEvaluation evaluate(
        Ipopt::Index n,
        const Ipopt::Number* x
    ) const {
        if (n != static_cast<Ipopt::Index>(initial_.size())) {
            throw std::invalid_argument("HELD2 search coordinate count changed");
        }
        std::vector<double> independent(x, x + n - 1);
        Held2StateEvaluation evaluation = evaluator_(independent, x[n - 1]);
        if (!use_tpd_) {
            return evaluation;
        }
        evaluation.objective -= reference_.objective;
        for (std::size_t index = 0; index < evaluation.gradient.size(); ++index) {
            evaluation.objective -= reference_.gradient[index]
                * (x[index] - reference_variables_[index]);
            evaluation.gradient[index] -= reference_.gradient[index];
        }
        return evaluation;
    }

    Held2StateEvaluator evaluator_;
    std::vector<double> reference_variables_;
    Held2StateEvaluation reference_;
    bool use_tpd_ = false;
    std::vector<double> initial_;
    std::vector<double> lower_;
    std::vector<double> upper_;
    double composition_sum_upper_ = std::numeric_limits<double>::infinity();
    std::vector<double> source_domain_coefficients_;
    double source_domain_upper_ = std::numeric_limits<double>::infinity();
    Held2SearchDomainPredicate trial_domain_predicate_;
    bool solver_converged_ = false;
    std::string callback_error_;
    std::vector<double> variables_;
    int iterations_ = -1;
    double final_step_norm_ = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> lower_bound_multipliers_;
    std::vector<double> upper_bound_multipliers_;
    std::vector<double> constraint_multipliers_;
};

Held2SearchRun solve_held2_search(
    const Held2StateEvaluator& evaluator,
    const std::vector<double>& reference_variables,
    const Held2StateEvaluation& reference,
    bool use_tpd,
    const std::vector<double>& initial,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    double composition_sum_upper = std::numeric_limits<double>::infinity(),
    const std::vector<double>& source_domain_coefficients = {},
    double source_domain_upper = std::numeric_limits<double>::infinity(),
    const Held2SearchDomainPredicate& trial_domain_predicate = {}
) {
    Ipopt::SmartPtr<Held2SearchTnlp> problem = new Held2SearchTnlp(
        evaluator,
        reference_variables,
        reference,
        use_tpd,
        initial,
        lower,
        upper,
        composition_sum_upper,
        source_domain_coefficients,
        source_domain_upper,
        trial_domain_predicate
    );
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application = IpoptApplicationFactory();
    configure_held2_ipopt(application);
    if (application->Initialize() != Ipopt::Solve_Succeeded) {
        return {};
    }
    const Ipopt::ApplicationReturnStatus status = application->OptimizeTNLP(problem);
    Held2SearchRun result = problem->result();
    result.solver_status = static_cast<int>(status);
    return result;
}

Held2SearchRun solve_held2_stage_ii_search(
    const Held2StateEvaluator& evaluator,
    const std::vector<double>& feed,
    const std::vector<double>& multiplier,
    const std::vector<double>& initial,
    const std::vector<double>& physical_lower,
    const std::vector<double>& physical_upper,
    double composition_sum_upper
) {
    const std::size_t dimension = feed.size();
    Held2SearchRun invalid;
    if (dimension == 0 || multiplier.size() != dimension
        || initial.size() != dimension + 1
        || physical_lower.size() != dimension + 1
        || physical_upper.size() != dimension + 1) {
        invalid.callback_error = "HELD2 Stage II simplex search dimensions changed";
        return invalid;
    }
    const std::vector<double> independent_initial(
        initial.begin(),
        initial.end() - 1
    );
    const std::vector<double> independent_lower(
        physical_lower.begin(),
        physical_lower.end() - 1
    );
    const std::vector<double> independent_upper(
        physical_upper.begin(),
        physical_upper.end() - 1
    );
    Held2StageIISimplexChart initial_chart;
    try {
        initial_chart = evaluate_held2_stage_ii_simplex_chart(
            independent_lower,
            independent_upper,
            composition_sum_upper,
            independent_initial,
            true
        );
    } catch (const std::exception& error) {
        invalid.callback_error = error.what();
        return invalid;
    }
    if (initial_chart.singular) {
        invalid.callback_error = "HELD2 Stage II simplex start is singular";
        return invalid;
    }

    std::vector<double> chart_initial = initial_chart.chart_coordinates;
    chart_initial.push_back(initial.back());
    std::vector<double> chart_lower(dimension, 0.0);
    std::vector<double> chart_upper(dimension, 1.0);
    chart_lower.push_back(physical_lower.back());
    chart_upper.push_back(physical_upper.back());
    const Held2StateEvaluator chart_evaluator = [
        evaluator,
        feed,
        multiplier,
        independent_lower,
        independent_upper,
        composition_sum_upper
    ](const std::vector<double>& chart_coordinates, double phase_coordinate) {
        const Held2StageIISimplexChart chart =
            evaluate_held2_stage_ii_simplex_chart(
                independent_lower,
                independent_upper,
                composition_sum_upper,
                chart_coordinates,
                false
            );
        if (chart.singular) {
            throw std::invalid_argument("HELD2 Stage II simplex trial is singular");
        }
        Held2StateEvaluation physical_state = evaluator(
            chart.physical_coordinates,
            phase_coordinate
        );
        if (physical_state.gradient.size() != feed.size() + 1
            || physical_state.hessian.size()
                != (feed.size() + 1) * (feed.size() + 1)) {
            throw std::invalid_argument(
                "HELD2 Stage II physical derivative block is incomplete"
            );
        }
        for (std::size_t coordinate = 0; coordinate < feed.size(); ++coordinate) {
            physical_state.objective += multiplier[coordinate]
                * (feed[coordinate] - chart.physical_coordinates[coordinate]);
            physical_state.gradient[coordinate] -= multiplier[coordinate];
        }
        return transform_held2_stage_ii_simplex_state(physical_state, chart);
    };
    Held2SearchRun run = solve_held2_search(
        chart_evaluator,
        {},
        {},
        false,
        chart_initial,
        chart_lower,
        chart_upper
    );
    if (run.variables.size() != dimension + 1) {
        return run;
    }

    const double terminal_phase_coordinate = run.variables.back();
    const std::vector<double> terminal_chart_coordinates(
        run.variables.begin(),
        run.variables.end() - 1
    );
    Held2StageIISimplexChart terminal_chart;
    try {
        terminal_chart = evaluate_held2_stage_ii_simplex_chart(
            independent_lower,
            independent_upper,
            composition_sum_upper,
            terminal_chart_coordinates,
            false
        );
    } catch (const std::exception& error) {
        run.callback_error = error.what();
        run.solver_converged = false;
        return run;
    }
    if (terminal_chart.singular) {
        run.callback_error = "HELD2 Stage II simplex terminal is singular";
        run.solver_converged = false;
        return run;
    }
    run.variables = terminal_chart.physical_coordinates;
    run.variables.push_back(terminal_phase_coordinate);
    run.coordinate_jacobian = std::move(terminal_chart.jacobian);
    run.constraint_multipliers.clear();
    return run;
}

std::array<double, 2> manufactured_modified_potentials(
    double composition,
    double molar_volume
) {
    const double gibbs = manufactured_gibbs(composition, molar_volume);
    const auto gradient = manufactured_gibbs_gradient(composition, molar_volume);
    const double common_volume_term = -molar_volume * gradient[1];
    return {
        gibbs - composition * gradient[0] + common_volume_term,
        gibbs + (1.0 - composition) * gradient[0] + common_volume_term,
    };
}

double enumerated_manufactured_objective(double feed_composition) {
    std::vector<double> compositions;
    std::vector<double> volumes;
    compositions.reserve(104);
    volumes.reserve(102);
    for (int index = 0; index <= 100; ++index) {
        compositions.push_back(static_cast<double>(index) / 100.0);
        volumes.push_back(0.5 + static_cast<double>(index) / 100.0);
    }
    compositions.insert(compositions.end(), {0.2, 0.5, 0.8});
    volumes.push_back(1.0);
    std::sort(compositions.begin(), compositions.end());
    compositions.erase(std::unique(compositions.begin(), compositions.end()), compositions.end());
    std::sort(volumes.begin(), volumes.end());
    volumes.erase(std::unique(volumes.begin(), volumes.end()), volumes.end());

    std::vector<double> phase_minima;
    phase_minima.reserve(compositions.size());
    for (double composition : compositions) {
        double minimum = std::numeric_limits<double>::infinity();
        for (double volume : volumes) {
            minimum = std::min(minimum, manufactured_gibbs(composition, volume));
        }
        phase_minima.push_back(minimum);
    }

    double best = std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < compositions.size(); ++index) {
        if (std::abs(compositions[index] - feed_composition) <= kCoordinateTolerance) {
            best = std::min(best, phase_minima[index]);
        }
        if (compositions[index] > feed_composition) {
            continue;
        }
        for (std::size_t right = index + 1; right < compositions.size(); ++right) {
            if (compositions[right] < feed_composition) {
                continue;
            }
            const double separation = compositions[right] - compositions[index];
            if (separation <= 0.0) {
                continue;
            }
            const double fraction =
                (compositions[right] - feed_composition) / separation;
            best = std::min(
                best,
                fraction * phase_minima[index]
                    + (1.0 - fraction) * phase_minima[right]
            );
        }
    }
    return best;
}

struct Held2ReducedTransform {
    std::vector<double> modified_fractions;
    std::vector<double> physical_amounts;
    std::vector<double> jacobian;
    double volume = 0.0;
    std::size_t provider_coordinate_count = 0;
    std::size_t reduced_coordinate_count = 0;
};

Held2ReducedTransform make_held2_reduced_transform(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent_modified_fractions,
    double log_volume
) {
    const std::size_t component_count = coordinates.charges.size();
    const std::size_t independent_count = coordinates.independent_indices.size();
    if (independent_modified_fractions.size() != independent_count) {
        throw std::invalid_argument(
            "independent modified composition size does not match the HELD2 chart"
        );
    }
    require_finite_vector(independent_modified_fractions, "independent modified fractions");
    if (!std::isfinite(log_volume)) {
        throw std::invalid_argument("HELD2 log volume must be finite");
    }

    Held2ReducedTransform result;
    result.provider_coordinate_count = component_count + 1;
    result.reduced_coordinate_count = independent_count + 1;
    result.modified_fractions.resize(coordinates.retained_indices.size(), 0.0);
    for (std::size_t independent = 0; independent < independent_count; ++independent) {
        const auto retained = static_cast<std::size_t>(
            std::find(
                coordinates.retained_indices.begin(),
                coordinates.retained_indices.end(),
                coordinates.independent_indices[independent]
            ) - coordinates.retained_indices.begin()
        );
        result.modified_fractions[retained] = independent_modified_fractions[independent];
    }
    result.physical_amounts = held2_lift_independent_fractions(
        coordinates,
        independent_modified_fractions
    );
    const auto dependent_retained = static_cast<std::size_t>(
        std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            coordinates.dependent_index
        ) - coordinates.retained_indices.begin()
    );
    result.modified_fractions[dependent_retained] = 1.0 - std::accumulate(
        independent_modified_fractions.begin(),
        independent_modified_fractions.end(),
        0.0
    );
    result.volume = std::exp(log_volume);
    if (!std::isfinite(result.volume) || result.volume <= 0.0) {
        throw std::invalid_argument("HELD2 phase volume must be finite and positive");
    }

    result.jacobian.assign(
        result.provider_coordinate_count * result.reduced_coordinate_count,
        0.0
    );
    for (std::size_t independent = 0; independent < independent_count; ++independent) {
        for (std::size_t retained = 0; retained < coordinates.retained_indices.size(); ++retained) {
            const std::size_t component = coordinates.retained_indices[retained];
            double modified_derivative = 0.0;
            if (component == coordinates.independent_indices[independent]) {
                modified_derivative = 1.0;
            } else if (component == coordinates.dependent_index) {
                modified_derivative = -1.0;
            }
            result.jacobian[
                component * result.reduced_coordinate_count + independent
            ] = modified_derivative / coordinates.modified_factors[retained];
        }
        double eliminated_derivative = 0.0;
        const double eliminated_charge = coordinates.charges[coordinates.eliminated_index];
        for (std::size_t component : coordinates.retained_indices) {
            eliminated_derivative -= coordinates.charges[component] / eliminated_charge
                * result.jacobian[
                    component * result.reduced_coordinate_count + independent
                ];
        }
        result.jacobian[
            coordinates.eliminated_index * result.reduced_coordinate_count + independent
        ] = eliminated_derivative;
    }
    result.jacobian[
        component_count * result.reduced_coordinate_count + independent_count
    ] = result.volume;
    return result;
}

Held2PackingEvaluation transform_held2_scalar(
    const Held2ReducedTransform& transform,
    double value,
    const std::vector<double>& gradient,
    const std::vector<double>& hessian
) {
    if (!std::isfinite(value)
        || gradient.size() != transform.provider_coordinate_count
        || hessian.size()
            != transform.provider_coordinate_count * transform.provider_coordinate_count) {
        throw std::invalid_argument("HELD2 physical scalar block has invalid dimensions");
    }
    require_finite_vector(gradient, "HELD2 physical scalar gradient");
    require_finite_vector(hessian, "HELD2 physical scalar Hessian");

    Held2PackingEvaluation result;
    result.value = value;
    result.gradient.assign(transform.reduced_coordinate_count, 0.0);
    for (std::size_t reduced = 0; reduced < transform.reduced_coordinate_count; ++reduced) {
        for (std::size_t provider = 0; provider < transform.provider_coordinate_count; ++provider) {
            result.gradient[reduced] += transform.jacobian[
                provider * transform.reduced_coordinate_count + reduced
            ] * gradient[provider];
        }
    }
    result.hessian.assign(
        transform.reduced_coordinate_count * transform.reduced_coordinate_count,
        0.0
    );
    for (std::size_t row = 0; row < transform.reduced_coordinate_count; ++row) {
        for (std::size_t column = 0; column < transform.reduced_coordinate_count; ++column) {
            double transformed = 0.0;
            for (std::size_t left = 0; left < transform.provider_coordinate_count; ++left) {
                for (std::size_t right = 0; right < transform.provider_coordinate_count; ++right) {
                    transformed += transform.jacobian[
                        left * transform.reduced_coordinate_count + row
                    ] * hessian[left * transform.provider_coordinate_count + right]
                        * transform.jacobian[
                            right * transform.reduced_coordinate_count + column
                        ];
                }
            }
            result.hessian[row * transform.reduced_coordinate_count + column] = transformed;
        }
    }
    result.hessian.back() += transform.volume * gradient.back();
    return result;
}

}  // namespace

namespace {

Held2StageIISimplexChart evaluate_held2_stage_ii_simplex_chart(
    const std::vector<double>& independent_lower_bounds,
    const std::vector<double>& independent_upper_bounds,
    double composition_sum_upper,
    const std::vector<double>& coordinates,
    bool inverse
) {
    const std::size_t dimension = independent_lower_bounds.size();
    if (dimension == 0 || coordinates.size() != dimension
        || independent_upper_bounds.size() != dimension
        || !std::isfinite(composition_sum_upper)) {
        throw std::invalid_argument("HELD2 Stage II simplex chart dimensions are invalid");
    }
    require_finite_vector(
        independent_lower_bounds,
        "HELD2 Stage II independent lower bounds"
    );
    require_finite_vector(
        independent_upper_bounds,
        "HELD2 Stage II independent upper bounds"
    );
    require_finite_vector(coordinates, "HELD2 Stage II simplex coordinates");
    const double free_scale = composition_sum_upper
        - std::accumulate(
            independent_lower_bounds.begin(),
            independent_lower_bounds.end(),
            0.0
        );
    if (!(free_scale > 0.0)
        || std::any_of(
            independent_lower_bounds.begin(),
            independent_lower_bounds.end(),
            [](double value) { return value < 0.0; }
        )) {
        throw std::invalid_argument("HELD2 Stage II simplex chart has no feasible interior");
    }
    for (std::size_t index = 0; index < dimension; ++index) {
        const double simplex_maximum = composition_sum_upper
            - std::accumulate(
                independent_lower_bounds.begin(),
                independent_lower_bounds.end(),
                0.0
            )
            + independent_lower_bounds[index];
        if (independent_upper_bounds[index] < simplex_maximum) {
            throw std::invalid_argument(
                "HELD2 Stage II simplex chart requires redundant independent upper bounds"
            );
        }
    }

    Held2StageIISimplexChart result;
    result.chart_coordinates.assign(dimension, 0.0);
    result.physical_coordinates.assign(dimension, 0.0);
    if (inverse) {
        double remaining = free_scale;
        for (std::size_t index = 0; index < dimension; ++index) {
            const double shifted = coordinates[index]
                - independent_lower_bounds[index];
            if (shifted < 0.0 || shifted > remaining || !(remaining > 0.0)) {
                throw std::invalid_argument(
                    "HELD2 Stage II physical start is outside the feasible simplex"
                );
            }
            result.chart_coordinates[index] = shifted / remaining;
            result.physical_coordinates[index] = coordinates[index];
            remaining -= shifted;
        }
    } else {
        result.chart_coordinates = coordinates;
        for (double value : result.chart_coordinates) {
            if (value < 0.0 || value > 1.0) {
                throw std::invalid_argument(
                    "HELD2 Stage II chart coordinate is outside [0, 1]"
                );
            }
        }
        for (std::size_t index = 0; index < dimension; ++index) {
            double product = free_scale * result.chart_coordinates[index];
            for (std::size_t previous = 0; previous < index; ++previous) {
                product *= 1.0 - result.chart_coordinates[previous];
            }
            result.physical_coordinates[index] =
                independent_lower_bounds[index] + product;
        }
    }

    result.jacobian.assign(dimension * dimension, 0.0);
    result.component_hessians.assign(dimension * dimension * dimension, 0.0);
    for (std::size_t component = 0; component < dimension; ++component) {
        for (std::size_t column = 0; column <= component; ++column) {
            double derivative = free_scale;
            if (column != component) {
                derivative *= -result.chart_coordinates[component];
            }
            for (std::size_t previous = 0; previous < component; ++previous) {
                if (previous != column) {
                    derivative *= 1.0 - result.chart_coordinates[previous];
                }
            }
            result.jacobian[component * dimension + column] = derivative;
        }
        if (std::abs(result.jacobian[component * dimension + component])
            <= std::numeric_limits<double>::epsilon()
                * std::max(1.0, free_scale)) {
            result.singular = true;
        }
        for (std::size_t row = 0; row <= component; ++row) {
            for (std::size_t column = 0; column <= component; ++column) {
                if (row == column || (row == component && column == component)) {
                    continue;
                }
                double derivative = free_scale;
                if (row != component && column != component) {
                    derivative *= result.chart_coordinates[component];
                } else {
                    derivative *= -1.0;
                }
                for (std::size_t previous = 0; previous < component; ++previous) {
                    if (previous != row && previous != column) {
                        derivative *= 1.0 - result.chart_coordinates[previous];
                    }
                }
                result.component_hessians[
                    (component * dimension + row) * dimension + column
                ] = derivative;
            }
        }
    }
    return result;
}

Held2StateEvaluation transform_held2_stage_ii_simplex_state(
    const Held2StateEvaluation& physical_state,
    const Held2StageIISimplexChart& chart
) {
    if (chart.singular) {
        throw std::invalid_argument("HELD2 Stage II simplex chart is singular");
    }
    const std::size_t dimension = chart.chart_coordinates.size();
    const std::size_t coordinate_count = dimension + 1;
    if (physical_state.gradient.size() != coordinate_count
        || physical_state.hessian.size() != coordinate_count * coordinate_count
        || chart.jacobian.size() != dimension * dimension
        || chart.component_hessians.size() != dimension * dimension * dimension) {
        throw std::invalid_argument("HELD2 Stage II simplex derivative block is incomplete");
    }

    const auto transform_derivatives = [&chart, dimension, coordinate_count](
                                           const std::vector<double>& gradient,
                                           const std::vector<double>& hessian
                                       ) {
        if (gradient.size() != coordinate_count
            || hessian.size() != coordinate_count * coordinate_count) {
            throw std::invalid_argument(
                "HELD2 Stage II simplex derivative dimensions changed"
            );
        }
        std::pair<std::vector<double>, std::vector<double>> transformed;
        transformed.first.assign(coordinate_count, 0.0);
        transformed.second.assign(coordinate_count * coordinate_count, 0.0);
        for (std::size_t chart_column = 0; chart_column < dimension; ++chart_column) {
            for (std::size_t physical = 0; physical < dimension; ++physical) {
                transformed.first[chart_column] += chart.jacobian[
                    physical * dimension + chart_column
                ] * gradient[physical];
            }
        }
        transformed.first.back() = gradient.back();
        for (std::size_t row = 0; row < dimension; ++row) {
            for (std::size_t column = 0; column < dimension; ++column) {
                double value = 0.0;
                for (std::size_t left = 0; left < dimension; ++left) {
                    value += gradient[left] * chart.component_hessians[
                        (left * dimension + row) * dimension + column
                    ];
                    for (std::size_t right = 0; right < dimension; ++right) {
                        value += chart.jacobian[left * dimension + row]
                            * hessian[left * coordinate_count + right]
                            * chart.jacobian[right * dimension + column];
                    }
                }
                transformed.second[row * coordinate_count + column] = value;
            }
            double cross = 0.0;
            for (std::size_t physical = 0; physical < dimension; ++physical) {
                cross += chart.jacobian[physical * dimension + row]
                    * hessian[physical * coordinate_count + dimension];
            }
            transformed.second[row * coordinate_count + dimension] = cross;
            transformed.second[dimension * coordinate_count + row] = cross;
        }
        transformed.second.back() = hessian.back();
        return transformed;
    };

    Held2StateEvaluation result = physical_state;
    std::tie(result.gradient, result.hessian) = transform_derivatives(
        physical_state.gradient,
        physical_state.hessian
    );
    if (physical_state.has_packing_evaluation) {
        std::tie(result.packing.gradient, result.packing.hessian) =
            transform_derivatives(
                physical_state.packing.gradient,
                physical_state.packing.hessian
            );
    }
    return result;
}

Held2PhysicalKkt evaluate_held2_stage_ii_physical_kkt(
    const std::vector<double>& physical_gradient,
    const std::vector<double>& master_multiplier,
    const std::vector<double>& variables,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    double composition_sum_upper,
    const std::vector<double>& chart_jacobian,
    const std::vector<double>& chart_lower_bound_multipliers,
    const std::vector<double>& chart_upper_bound_multipliers
) {
    const std::size_t dimension = master_multiplier.size();
    if (physical_gradient.size() != dimension + 1
        || variables.size() != dimension + 1
        || lower.size() != dimension + 1
        || upper.size() != dimension + 1
        || chart_jacobian.size() != dimension * dimension
        || chart_lower_bound_multipliers.size() != dimension + 1
        || chart_upper_bound_multipliers.size() != dimension + 1) {
        throw std::invalid_argument("HELD2 Stage II physical KKT dimensions changed");
    }
    require_finite_vector(
        physical_gradient,
        "HELD2 Stage II physical KKT gradient"
    );
    require_finite_vector(
        master_multiplier,
        "HELD2 Stage II physical KKT master multiplier"
    );
    require_finite_vector(variables, "HELD2 Stage II physical KKT variables");
    require_finite_vector(lower, "HELD2 Stage II physical KKT lower bounds");
    require_finite_vector(upper, "HELD2 Stage II physical KKT upper bounds");
    require_finite_vector(
        chart_jacobian,
        "HELD2 Stage II physical KKT chart Jacobian"
    );
    require_finite_vector(
        chart_lower_bound_multipliers,
        "HELD2 Stage II physical KKT chart lower multipliers"
    );
    require_finite_vector(
        chart_upper_bound_multipliers,
        "HELD2 Stage II physical KKT chart upper multipliers"
    );
    if (!std::isfinite(composition_sum_upper)) {
        throw std::invalid_argument(
            "HELD2 Stage II physical KKT composition sum must be finite"
        );
    }
    for (std::size_t index = 0; index < variables.size(); ++index) {
        if (lower[index] > upper[index]
            || variables[index] < lower[index] - kCertificateTolerance
            || variables[index] > upper[index] + kCertificateTolerance) {
            throw std::invalid_argument(
                "HELD2 Stage II physical KKT point is outside its physical domain"
            );
        }
    }
    std::vector<double> gradient(dimension, 0.0);
    for (std::size_t index = 0; index < dimension; ++index) {
        gradient[index] = physical_gradient[index] - master_multiplier[index];
    }
    const double simplex_slack = composition_sum_upper
        - std::accumulate(variables.begin(), variables.end() - 1, 0.0);
    if (simplex_slack < -kCertificateTolerance) {
        throw std::invalid_argument(
            "HELD2 Stage II physical KKT point is outside its physical domain"
        );
    }

    Held2PhysicalKkt result;
    result.dual_signs_valid = std::all_of(
        chart_lower_bound_multipliers.begin(),
        chart_lower_bound_multipliers.end(),
        [](double value) { return value >= 0.0; }
    ) && std::all_of(
        chart_upper_bound_multipliers.begin(),
        chart_upper_bound_multipliers.end(),
        [](double value) { return value >= 0.0; }
    );
    if (!result.dual_signs_valid) {
        return result;
    }

    std::vector<double> chart_bound_contribution(dimension + 1, 0.0);
    for (std::size_t index = 0; index <= dimension; ++index) {
        chart_bound_contribution[index] = chart_upper_bound_multipliers[index]
            - chart_lower_bound_multipliers[index];
    }
    std::vector<double> physical_bound_contribution(dimension + 1, 0.0);
    for (std::size_t column = dimension; column-- > 0;) {
        double value = chart_bound_contribution[column];
        for (std::size_t row = column + 1; row < dimension; ++row) {
            value -= chart_jacobian[row * dimension + column]
                * physical_bound_contribution[row];
        }
        const double diagonal = chart_jacobian[column * dimension + column];
        if (std::abs(diagonal) <= std::numeric_limits<double>::epsilon()) {
            return result;
        }
        physical_bound_contribution[column] = value / diagonal;
    }
    physical_bound_contribution.back() = chart_bound_contribution.back();

    const auto multiplier_cap = [](double slack) {
        return slack == 0.0
            ? std::numeric_limits<double>::infinity()
            : kCertificateTolerance / std::abs(slack);
    };
    double simplex_lower = 0.0;
    double simplex_upper = multiplier_cap(simplex_slack);
    for (std::size_t index = 0; index < dimension; ++index) {
        const double lower_cap = multiplier_cap(variables[index] - lower[index]);
        const double upper_cap = multiplier_cap(upper[index] - variables[index]);
        simplex_lower = std::max(
            simplex_lower,
            physical_bound_contribution[index] - upper_cap
        );
        simplex_upper = std::min(
            simplex_upper,
            physical_bound_contribution[index] + lower_cap
        );
    }
    if (simplex_lower > simplex_upper + kCoordinateTolerance) {
        return result;
    }
    const double simplex_multiplier = std::isfinite(simplex_upper)
        ? simplex_lower + 0.5 * (simplex_upper - simplex_lower)
        : simplex_lower;

    std::vector<double> reconstructed_contribution(dimension + 1, 0.0);
    double complementarity = std::abs(simplex_multiplier * simplex_slack);
    for (std::size_t index = 0; index < dimension; ++index) {
        const double difference = physical_bound_contribution[index]
            - simplex_multiplier;
        const double lower_multiplier = std::max(0.0, -difference);
        const double upper_multiplier = std::max(0.0, difference);
        reconstructed_contribution[index] = -lower_multiplier
            + upper_multiplier + simplex_multiplier;
        complementarity = std::max({
            complementarity,
            std::abs(lower_multiplier * (variables[index] - lower[index])),
            std::abs(upper_multiplier * (upper[index] - variables[index])),
        });
    }
    const double phase_lower_cap = multiplier_cap(variables.back() - lower.back());
    const double phase_upper_cap = multiplier_cap(upper.back() - variables.back());
    if (physical_bound_contribution.back() < -phase_lower_cap - kCoordinateTolerance
        || physical_bound_contribution.back()
            > phase_upper_cap + kCoordinateTolerance) {
        return result;
    }
    const double phase_lower_multiplier = std::max(
        0.0,
        -physical_bound_contribution.back()
    );
    const double phase_upper_multiplier = std::max(
        0.0,
        physical_bound_contribution.back()
    );
    reconstructed_contribution.back() = -phase_lower_multiplier
        + phase_upper_multiplier;
    complementarity = std::max({
        complementarity,
        std::abs(phase_lower_multiplier * (variables.back() - lower.back())),
        std::abs(phase_upper_multiplier * (upper.back() - variables.back())),
    });

    double reconstruction = 0.0;
    for (std::size_t index = 0; index <= dimension; ++index) {
        reconstruction = std::max(
            reconstruction,
            std::abs(
                physical_bound_contribution[index]
                - reconstructed_contribution[index]
            )
        );
    }
    for (std::size_t column = 0; column < dimension; ++column) {
        double chart_reconstruction = 0.0;
        for (std::size_t row = 0; row < dimension; ++row) {
            chart_reconstruction += chart_jacobian[row * dimension + column]
                * reconstructed_contribution[row];
        }
        reconstruction = std::max(
            reconstruction,
            std::abs(chart_bound_contribution[column] - chart_reconstruction)
        );
    }
    reconstruction = std::max(
        reconstruction,
        std::abs(
            chart_bound_contribution.back() - reconstructed_contribution.back()
        )
    );

    double stationarity = 0.0;
    for (std::size_t index = 0; index < dimension; ++index) {
        stationarity = std::max(
            stationarity,
            std::abs(gradient[index] + reconstructed_contribution[index])
        );
    }
    stationarity = std::max(
        stationarity,
        std::abs(physical_gradient.back() + reconstructed_contribution.back())
    );
    result.stationarity_inf_norm = stationarity;
    result.complementarity = complementarity;
    result.reconstruction_inf_norm = reconstruction;
    return result;
}

}  // namespace

std::tuple<
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    double,
    double,
    double,
    bool,
    bool> evaluate_held2_stage_ii_simplex_test_adapter(
    const std::vector<double>& independent_lower_bounds,
    const std::vector<double>& independent_upper_bounds,
    double composition_sum_upper,
    const std::vector<double>& values,
    const std::vector<double>& physical_gradient,
    const std::vector<double>& physical_hessian,
    const std::vector<double>& master_multiplier,
    const std::array<double, 2>& phase_bounds,
    bool inverse,
    bool physical_kkt,
    const std::vector<double>& chart_lower_bound_multipliers,
    const std::vector<double>& chart_upper_bound_multipliers
) {
    if (physical_kkt) {
        if (values.size() != independent_lower_bounds.size() + 1) {
            throw std::invalid_argument(
                "HELD2 Stage II physical KKT dimensions changed"
            );
        }
        std::vector<double> lower = independent_lower_bounds;
        std::vector<double> upper = independent_upper_bounds;
        lower.push_back(phase_bounds[0]);
        upper.push_back(phase_bounds[1]);
        const std::vector<double> independent(
            values.begin(),
            values.end() - 1
        );
        const Held2StageIISimplexChart chart =
            evaluate_held2_stage_ii_simplex_chart(
                independent_lower_bounds,
                independent_upper_bounds,
                composition_sum_upper,
                independent,
                true
            );
        const Held2PhysicalKkt kkt = evaluate_held2_stage_ii_physical_kkt(
            physical_gradient,
            master_multiplier,
            values,
            lower,
            upper,
            composition_sum_upper,
            chart.jacobian,
            chart_lower_bound_multipliers,
            chart_upper_bound_multipliers
        );
        return {{}, {}, {}, {}, {}, {}, kkt.stationarity_inf_norm,
                kkt.complementarity, kkt.reconstruction_inf_norm,
                kkt.dual_signs_valid, false};
    }
    const Held2StageIISimplexChart chart = evaluate_held2_stage_ii_simplex_chart(
        independent_lower_bounds,
        independent_upper_bounds,
        composition_sum_upper,
        values,
        inverse
    );
    std::vector<double> transformed_gradient;
    std::vector<double> transformed_hessian;
    if (!physical_gradient.empty() || !physical_hessian.empty()) {
        Held2StateEvaluation physical_state;
        physical_state.gradient = physical_gradient;
        physical_state.hessian = physical_hessian;
        Held2StateEvaluation transformed = transform_held2_stage_ii_simplex_state(
            physical_state,
            chart
        );
        transformed_gradient = std::move(transformed.gradient);
        transformed_hessian = std::move(transformed.hessian);
    }
    return {
        chart.physical_coordinates,
        chart.chart_coordinates,
        chart.jacobian,
        chart.component_hessians,
        std::move(transformed_gradient),
        std::move(transformed_hessian),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        false,
        chart.singular,
    };
}

Held2StateEvaluation evaluate_held2_manufactured_state(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent_modified_fractions,
    double log_volume
) {
    return evaluate_manufactured_state_impl(
        coordinates,
        independent_modified_fractions,
        log_volume
    );
}

double held2_manufactured_enumerated_objective(double feed_composition) {
    return enumerated_manufactured_objective(feed_composition);
}

Held2Coordinates make_held2_coordinates(const std::vector<double>& charges) {
    if (charges.size() < 3) {
        throw std::invalid_argument("HELD2 requires at least three species");
    }
    require_finite_vector(charges, "charge numbers");
    std::size_t charged_count = 0;
    std::size_t neutral_count = 0;
    double largest_abs_charge = 0.0;
    for (double charge : charges) {
        if (charge == 0.0) {
            ++neutral_count;
        } else {
            ++charged_count;
            largest_abs_charge = std::max(largest_abs_charge, std::abs(charge));
        }
    }
    if (charged_count < 2 || neutral_count < 1) {
        throw std::invalid_argument(
            "HELD2 requires at least two charged species and one molecular species"
        );
    }

    std::size_t eliminated = charges.size();
    for (std::size_t reverse = charges.size(); reverse > 0; --reverse) {
        const std::size_t candidate = reverse - 1;
        if (std::abs(charges[candidate]) != largest_abs_charge) {
            continue;
        }
        bool nonsingular = true;
        for (std::size_t index = 0; index < charges.size(); ++index) {
            if (index == candidate) {
                continue;
            }
            const double factor = 1.0 - charges[index] / charges[candidate];
            if (!(factor > 0.0) || !std::isfinite(factor)) {
                nonsingular = false;
                break;
            }
        }
        if (nonsingular) {
            eliminated = candidate;
            break;
        }
    }
    if (eliminated == charges.size()) {
        throw std::invalid_argument(
            "largest-absolute-charge choice has a singular modified coordinate factor"
        );
    }

    Held2Coordinates result;
    result.charges = charges;
    result.eliminated_index = eliminated;
    for (std::size_t index = 0; index < charges.size(); ++index) {
        if (index == eliminated) {
            continue;
        }
        result.retained_indices.push_back(index);
        result.modified_factors.push_back(1.0 - charges[index] / charges[eliminated]);
        if (charges[index] == 0.0) {
            result.dependent_index = index;
        }
    }
    for (std::size_t retained = 0; retained < result.retained_indices.size(); ++retained) {
        const std::size_t index = result.retained_indices[retained];
        if (index == result.dependent_index) {
            continue;
        }
        result.independent_indices.push_back(index);
        const double factor = result.modified_factors[retained];
        result.independent_lower_bounds.push_back(kHeld2ModifiedLowerScale * factor);
        if (charges[index] == 0.0) {
            result.independent_upper_bounds.push_back(1.0);
            continue;
        }
        double largest_opposite = 0.0;
        for (double other_charge : charges) {
            if (other_charge * charges[index] < 0.0) {
                largest_opposite = std::max(largest_opposite, std::abs(other_charge));
            }
        }
        if (largest_opposite == 0.0) {
            throw std::invalid_argument("charged species has no opposite-sign counterion");
        }
        const double physical_upper = largest_opposite
            / (std::abs(charges[index]) + largest_opposite);
        result.independent_upper_bounds.push_back(factor * physical_upper);
    }
    return result;
}

std::vector<double> held2_transform_physical_fractions(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_fractions
) {
    if (physical_fractions.size() != coordinates.charges.size()) {
        throw std::invalid_argument("physical composition size does not match charges");
    }
    require_finite_vector(physical_fractions, "physical mole fractions");
    if (!std::all_of(physical_fractions.begin(), physical_fractions.end(), [](double value) {
            return value >= 0.0;
        })) {
        throw std::invalid_argument("physical mole fractions must be nonnegative");
    }
    double total = 0.0;
    for (double value : physical_fractions) {
        total += value;
    }
    if (std::abs(total - 1.0) > kCoordinateTolerance) {
        throw std::invalid_argument("physical mole fractions must sum to one");
    }
    if (std::abs(charge_residual(coordinates.charges, physical_fractions))
        > kCoordinateTolerance) {
        throw std::invalid_argument("physical feed must be electroneutral");
    }
    std::vector<double> modified;
    modified.reserve(coordinates.retained_indices.size());
    for (std::size_t retained = 0; retained < coordinates.retained_indices.size(); ++retained) {
        modified.push_back(
            coordinates.modified_factors[retained]
            * physical_fractions[coordinates.retained_indices[retained]]
        );
    }
    double modified_total = 0.0;
    for (double value : modified) {
        modified_total += value;
    }
    if (std::abs(modified_total - 1.0) > kCoordinateTolerance) {
        throw std::invalid_argument("modified mole fractions do not sum to one");
    }
    return modified;
}

std::vector<double> held2_lift_modified_fractions(
    const Held2Coordinates& coordinates,
    const std::vector<double>& modified_fractions
) {
    if (modified_fractions.size() != coordinates.retained_indices.size()) {
        throw std::invalid_argument("modified composition size does not match retained species");
    }
    require_finite_vector(modified_fractions, "modified mole fractions");
    double modified_total = 0.0;
    for (double value : modified_fractions) {
        if (value < 0.0) {
            throw std::invalid_argument("modified mole fractions must be nonnegative");
        }
        modified_total += value;
    }
    if (std::abs(modified_total - 1.0) > kCoordinateTolerance) {
        throw std::invalid_argument("modified mole fractions must sum to one");
    }

    std::vector<double> physical(coordinates.charges.size(), 0.0);
    for (std::size_t retained = 0; retained < coordinates.retained_indices.size(); ++retained) {
        physical[coordinates.retained_indices[retained]] =
            modified_fractions[retained] / coordinates.modified_factors[retained];
    }
    double eliminated = 0.0;
    const double eliminated_charge = coordinates.charges[coordinates.eliminated_index];
    for (std::size_t index : coordinates.retained_indices) {
        eliminated -= coordinates.charges[index] / eliminated_charge * physical[index];
    }
    if (eliminated < -kCoordinateTolerance) {
        throw std::invalid_argument("eliminated ion amount must be nonnegative");
    }
    physical[coordinates.eliminated_index] = std::max(0.0, eliminated);
    return physical;
}

std::vector<double> held2_lift_independent_fractions(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent_modified_fractions
) {
    const std::size_t independent_count = coordinates.independent_indices.size();
    if (independent_modified_fractions.size() != independent_count) {
        throw std::invalid_argument(
            "independent modified composition size does not match the HELD2 chart"
        );
    }
    require_finite_vector(independent_modified_fractions, "independent modified fractions");
    std::vector<double> modified_fractions(coordinates.retained_indices.size(), 0.0);
    double independent_sum = 0.0;
    for (std::size_t independent = 0; independent < independent_count; ++independent) {
        const double value = independent_modified_fractions[independent];
        if (value < coordinates.independent_lower_bounds[independent]
            || value > coordinates.independent_upper_bounds[independent]) {
            throw std::invalid_argument("independent modified fraction is outside its source bound");
        }
        const auto retained = static_cast<std::size_t>(
            std::find(
                coordinates.retained_indices.begin(),
                coordinates.retained_indices.end(),
                coordinates.independent_indices[independent]
            ) - coordinates.retained_indices.begin()
        );
        modified_fractions[retained] = value;
        independent_sum += value;
    }
    if (!(independent_sum < 1.0)) {
        throw std::invalid_argument("independent modified fractions must leave a positive dependent fraction");
    }
    const auto dependent_retained = static_cast<std::size_t>(
        std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            coordinates.dependent_index
        ) - coordinates.retained_indices.begin()
    );
    modified_fractions[dependent_retained] = 1.0 - independent_sum;
    return held2_lift_modified_fractions(coordinates, modified_fractions);
}

double held2_stage_i_total_ion_mole_fraction(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent_modified_fractions
) {
    if (independent_modified_fractions.size()
        != coordinates.independent_indices.size()) {
        throw std::invalid_argument(
            "independent modified composition size does not match the HELD2 chart"
        );
    }
    require_finite_vector(
        independent_modified_fractions,
        "independent modified fractions"
    );
    double total = 0.0;
    for (std::size_t index = 0; index < independent_modified_fractions.size(); ++index) {
        if (coordinates.charges[coordinates.independent_indices[index]] != 0.0) {
            total += independent_modified_fractions[index];
        }
    }
    return total;
}

std::vector<double> held2_transform_modified_potentials(
    const Held2Coordinates& coordinates,
    const std::vector<double>& chemical_potentials
) {
    if (chemical_potentials.size() != coordinates.charges.size()) {
        throw std::invalid_argument("chemical-potential size does not match charges");
    }
    require_finite_vector(chemical_potentials, "chemical potentials");
    const double eliminated_potential =
        chemical_potentials[coordinates.eliminated_index];
    const double eliminated_charge = coordinates.charges[coordinates.eliminated_index];
    std::vector<double> modified;
    modified.reserve(coordinates.retained_indices.size());
    for (std::size_t retained = 0; retained < coordinates.retained_indices.size(); ++retained) {
        const std::size_t index = coordinates.retained_indices[retained];
        const double constrained = chemical_potentials[index]
            - coordinates.charges[index] / eliminated_charge * eliminated_potential;
        modified.push_back(constrained / coordinates.modified_factors[retained]);
    }
    return modified;
}

Held2StateEvaluation evaluate_held2_phase_block(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent_modified_fractions,
    double log_volume,
    double pressure_over_rt,
    double target_pressure_pa,
    const Held2PhysicalPhaseBlock& block
) {
    if (!std::isfinite(log_volume) || !std::isfinite(pressure_over_rt)
        || !std::isfinite(target_pressure_pa) || pressure_over_rt <= 0.0
        || target_pressure_pa <= 0.0) {
        throw std::invalid_argument("HELD2 phase state and pressure scales must be finite and positive");
    }
    if (!std::isfinite(block.helmholtz_over_rt) || !std::isfinite(block.pressure_pa)) {
        throw std::invalid_argument("HELD2 physical phase block scalars must be finite");
    }

    const Held2ReducedTransform transform = make_held2_reduced_transform(
        coordinates,
        independent_modified_fractions,
        log_volume
    );
    std::vector<double> augmented_gradient = block.gradient;
    if (augmented_gradient.empty()) {
        throw std::invalid_argument("HELD2 physical phase gradient must not be empty");
    }
    augmented_gradient.back() += pressure_over_rt;
    const Held2PackingEvaluation scalar = transform_held2_scalar(
        transform,
        block.helmholtz_over_rt + pressure_over_rt * transform.volume,
        augmented_gradient,
        block.hessian
    );
    Held2StateEvaluation result;
    result.modified_fractions = transform.modified_fractions;
    result.physical_amounts = transform.physical_amounts;
    result.volume = transform.volume;
    result.objective = scalar.value;
    result.gradient = scalar.gradient;
    result.hessian = scalar.hessian;
    std::vector<double> physical_potentials(block.gradient.begin(), block.gradient.end() - 1);
    result.modified_potentials = held2_transform_modified_potentials(
        coordinates,
        physical_potentials
    );
    result.pressure_stationarity_relative =
        (block.pressure_pa - target_pressure_pa) / target_pressure_pa;
    result.log_volume_gradient = result.gradient.back();
    result.pressure_stationarity_derivative_log_volume =
        (result.log_volume_gradient - result.hessian.back())
        / (pressure_over_rt * result.volume);
    return result;
}

Held2PackingEvaluation evaluate_held2_packing_block(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent_modified_fractions,
    double log_volume,
    double packing_fraction,
    const std::vector<double>& gradient,
    const std::vector<double>& hessian
) {
    return transform_held2_scalar(
        make_held2_reduced_transform(
            coordinates,
            independent_modified_fractions,
            log_volume
        ),
        packing_fraction,
        gradient,
        hessian
    );
}

Held2StateEvaluation evaluate_held2_log_packing_state(
    const Held2StateEvaluator& log_volume_evaluator,
    const std::vector<double>& independent_modified_fractions,
    double log_packing_fraction,
    const std::array<double, 2>& molar_volume_bounds
) {
    if (!std::isfinite(log_packing_fraction)
        || log_packing_fraction < std::log(kHeld2PackingFractionMinimum)
        || log_packing_fraction > std::log(kHeld2PackingFractionMaximum)
        || !std::isfinite(molar_volume_bounds[0])
        || !std::isfinite(molar_volume_bounds[1])
        || molar_volume_bounds[0] <= 0.0
        || molar_volume_bounds[1] <= molar_volume_bounds[0]) {
        throw std::invalid_argument("HELD2 log-packing coordinate is outside its domain");
    }
    const auto strict_log_volume_bounds = strict_interior_log_bounds(
        molar_volume_bounds
    );
    double left = strict_log_volume_bounds[0];
    double right = strict_log_volume_bounds[1];
    Held2StateEvaluation left_state = log_volume_evaluator(
        independent_modified_fractions, left
    );
    Held2StateEvaluation right_state = log_volume_evaluator(
        independent_modified_fractions, right
    );
    const auto residual = [log_packing_fraction](const Held2StateEvaluation& state) {
        if (!state.has_packing_evaluation || !std::isfinite(state.packing.value)
            || state.packing.value <= 0.0) {
            throw std::invalid_argument("HELD2 Provider packing evaluation is incomplete");
        }
        return std::log(state.packing.value) - log_packing_fraction;
    };
    double left_residual = residual(left_state);
    double right_residual = residual(right_state);
    if (left_residual * right_residual > 0.0) {
        throw std::invalid_argument("HELD2 Provider packing root is not bracketed");
    }

    Held2StateEvaluation physical_state;
    double root_residual = std::numeric_limits<double>::infinity();
    double log_volume = 0.5 * (left + right);
    for (int iteration = 0; iteration < 100; ++iteration) {
        physical_state = log_volume_evaluator(
            independent_modified_fractions, log_volume
        );
        root_residual = residual(physical_state);
        if (std::abs(root_residual) <= 1.0e-13) {
            break;
        }
        if (left_residual * root_residual <= 0.0) {
            right = log_volume;
            right_residual = root_residual;
        } else {
            left = log_volume;
            left_residual = root_residual;
        }
        const double log_packing_volume_derivative =
            physical_state.packing.gradient.back() / physical_state.packing.value;
        const double newton = log_volume
            - root_residual / log_packing_volume_derivative;
        log_volume = std::isfinite(newton) && newton > left && newton < right
            ? newton
            : 0.5 * (left + right);
    }
    if (std::abs(root_residual) > 1.0e-13) {
        throw std::runtime_error("HELD2 Provider packing root refinement was incomplete");
    }

    const std::size_t coordinate_count = physical_state.gradient.size();
    if (coordinate_count == 0
        || physical_state.hessian.size() != coordinate_count * coordinate_count
        || physical_state.packing.gradient.size() != coordinate_count
        || physical_state.packing.hessian.size() != coordinate_count * coordinate_count) {
        throw std::invalid_argument("HELD2 log-packing derivative block is incomplete");
    }
    const double packing = physical_state.packing.value;
    std::vector<double> log_packing_gradient(coordinate_count, 0.0);
    std::vector<double> log_packing_hessian(
        coordinate_count * coordinate_count, 0.0
    );
    for (std::size_t row = 0; row < coordinate_count; ++row) {
        log_packing_gradient[row] = physical_state.packing.gradient[row] / packing;
        for (std::size_t column = 0; column < coordinate_count; ++column) {
            log_packing_hessian[row * coordinate_count + column] =
                physical_state.packing.hessian[row * coordinate_count + column] / packing
                - physical_state.packing.gradient[row]
                    * physical_state.packing.gradient[column] / (packing * packing);
        }
    }
    const std::size_t volume_index = coordinate_count - 1;
    const double log_packing_volume_gradient = log_packing_gradient[volume_index];
    if (!std::isfinite(log_packing_volume_gradient)
        || std::abs(log_packing_volume_gradient)
            <= std::numeric_limits<double>::epsilon()) {
        throw std::invalid_argument("HELD2 log-packing chart is singular");
    }

    std::vector<double> jacobian(coordinate_count * coordinate_count, 0.0);
    for (std::size_t index = 0; index < volume_index; ++index) {
        jacobian[index * coordinate_count + index] = 1.0;
        jacobian[volume_index * coordinate_count + index] =
            -log_packing_gradient[index] / log_packing_volume_gradient;
    }
    jacobian.back() = 1.0 / log_packing_volume_gradient;

    std::vector<double> log_volume_hessian(
        coordinate_count * coordinate_count, 0.0
    );
    for (std::size_t row = 0; row < coordinate_count; ++row) {
        for (std::size_t column = 0; column < coordinate_count; ++column) {
            double contraction = 0.0;
            for (std::size_t left_index = 0; left_index < coordinate_count; ++left_index) {
                for (std::size_t right_index = 0; right_index < coordinate_count; ++right_index) {
                    contraction += jacobian[left_index * coordinate_count + row]
                        * log_packing_hessian[left_index * coordinate_count + right_index]
                        * jacobian[right_index * coordinate_count + column];
                }
            }
            log_volume_hessian[row * coordinate_count + column] =
                -contraction / log_packing_volume_gradient;
        }
    }

    Held2StateEvaluation result = physical_state;
    result.log_volume_gradient = physical_state.gradient.back();
    result.gradient.assign(coordinate_count, 0.0);
    result.hessian.assign(coordinate_count * coordinate_count, 0.0);
    for (std::size_t reduced = 0; reduced < coordinate_count; ++reduced) {
        for (std::size_t physical = 0; physical < coordinate_count; ++physical) {
            result.gradient[reduced] += jacobian[physical * coordinate_count + reduced]
                * physical_state.gradient[physical];
        }
    }
    for (std::size_t row = 0; row < coordinate_count; ++row) {
        for (std::size_t column = 0; column < coordinate_count; ++column) {
            double transformed = physical_state.gradient[volume_index]
                * log_volume_hessian[row * coordinate_count + column];
            for (std::size_t left_index = 0; left_index < coordinate_count; ++left_index) {
                for (std::size_t right_index = 0; right_index < coordinate_count; ++right_index) {
                    transformed += jacobian[left_index * coordinate_count + row]
                        * physical_state.hessian[left_index * coordinate_count + right_index]
                        * jacobian[right_index * coordinate_count + column];
                }
            }
            result.hessian[row * coordinate_count + column] = transformed;
        }
    }
    result.packing.gradient.assign(coordinate_count, 0.0);
    result.packing.gradient.back() = packing;
    result.packing.hessian.assign(coordinate_count * coordinate_count, 0.0);
    result.packing.hessian.back() = packing;
    return result;
}

Held2StageIResult solve_held2_stage_i(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const Held2StateEvaluator& log_volume_evaluator,
    const Held2StateEvaluator& search_evaluator,
    const std::array<double, 2>& molar_volume_bounds,
    const Held2VolumeBoundsEvaluator& volume_bounds_evaluator,
    bool search_uses_log_packing,
    bool volume_domain_search_complete,
    int solve_start_limit,
    double total_ion_mole_fraction_max
) {
    if (!std::isfinite(molar_volume_bounds[0]) || !std::isfinite(molar_volume_bounds[1])
        || molar_volume_bounds[0] <= 0.0
        || molar_volume_bounds[1] <= molar_volume_bounds[0]) {
        throw std::invalid_argument(
            "HELD2 Stage I molar-volume bounds must be finite, positive, and ordered"
        );
    }
    if (!std::isnan(total_ion_mole_fraction_max)
        && (!std::isfinite(total_ion_mole_fraction_max)
            || total_ion_mole_fraction_max < 0.0
            || total_ion_mole_fraction_max > 1.0)) {
        throw std::invalid_argument(
            "HELD2 Stage I total-ion source bound must be NaN or in [0, 1]"
        );
    }
    const std::vector<double> modified_feed = held2_transform_physical_fractions(
        coordinates,
        physical_feed
    );
    std::vector<double> feed_independent;
    feed_independent.reserve(coordinates.independent_indices.size());
    for (std::size_t independent : coordinates.independent_indices) {
        const auto retained = static_cast<std::size_t>(
            std::find(
                coordinates.retained_indices.begin(),
                coordinates.retained_indices.end(),
                independent
            ) - coordinates.retained_indices.begin()
        );
        feed_independent.push_back(modified_feed[retained]);
    }
    if (std::isfinite(total_ion_mole_fraction_max)
        && held2_stage_i_total_ion_mole_fraction(coordinates, feed_independent)
            > total_ion_mole_fraction_max) {
        throw std::invalid_argument(
            "HELD2 feed ion mole fraction exceeds the Provider source domain"
        );
    }
    std::vector<double> source_domain_coefficients(feed_independent.size(), 0.0);
    for (std::size_t index = 0; index < feed_independent.size(); ++index) {
        if (coordinates.charges[coordinates.independent_indices[index]] != 0.0) {
            source_domain_coefficients[index] = 1.0;
        }
    }

    Held2StageIResult result;
    result.volume_domain_search_complete = volume_domain_search_complete;
    result.declared_start_count = 10 * static_cast<int>(coordinates.charges.size());
    result.reference_scan_interval_count = result.declared_start_count;
    result.minimum_tpd = std::numeric_limits<double>::infinity();
    const auto strict_log_volume_bounds = strict_interior_log_bounds(
        molar_volume_bounds
    );
    const double reference_log_lower = strict_log_volume_bounds[0];
    const double reference_log_upper = strict_log_volume_bounds[1];
    const double reference_log_span = reference_log_upper - reference_log_lower;
    std::vector<std::pair<double, Held2StateEvaluation>> reference_scan;
    reference_scan.reserve(
        static_cast<std::size_t>(result.reference_scan_interval_count + 1)
    );
    for (int point = 0; point <= result.reference_scan_interval_count; ++point) {
        const double fraction = static_cast<double>(point)
            / static_cast<double>(result.reference_scan_interval_count);
        const double log_volume = point == 0
            ? reference_log_lower
            : point == result.reference_scan_interval_count
            ? reference_log_upper
            : reference_log_lower + fraction * reference_log_span;
        try {
            reference_scan.emplace_back(
                log_volume,
                log_volume_evaluator(feed_independent, log_volume)
            );
            ++result.reference_scan_point_count;
        } catch (const std::exception&) {
            ++result.reference_evaluation_failure_count;
        }
    }
    if (result.reference_evaluation_failure_count != 0
        || result.reference_scan_point_count
            != result.reference_scan_interval_count + 1) {
        result.outcome = "indeterminate";
        result.reference_failure_reason = "reference_scan_incomplete";
        return result;
    }

    const auto brackets_zero = [](double left, double right) {
        return (left <= 0.0 && right >= 0.0)
            || (left >= 0.0 && right <= 0.0);
    };
    const auto refine_bracket = [
        &log_volume_evaluator,
        &feed_independent,
        &result,
        &brackets_zero
    ](
        double left_log_volume,
        Held2StateEvaluation left_state,
        double right_log_volume,
        Held2StateEvaluation right_state,
        const auto& residual,
        double tolerance,
        bool continue_after_certificate,
        bool allow_representable_bracket,
        double& root_log_volume,
        Held2StateEvaluation& root_state,
        Held2StateEvaluation& classification_state
    ) {
        double left_residual = residual(left_state);
        double right_residual = residual(right_state);
        bool certified = false;
        if (std::abs(left_residual) <= tolerance) {
            root_log_volume = left_log_volume;
            root_state = left_state;
            certified = true;
            if (!continue_after_certificate) {
                classification_state = root_state;
                return true;
            }
        }
        if (!certified && std::abs(right_residual) <= tolerance) {
            root_log_volume = right_log_volume;
            root_state = right_state;
            certified = true;
            if (!continue_after_certificate) {
                classification_state = root_state;
                return true;
            }
        }
        if (!brackets_zero(left_residual, right_residual)) {
            return false;
        }
        for (int iteration = 0; iteration < 100; ++iteration) {
            if (right_log_volume - left_log_volume <= kCoordinateTolerance) {
                if (std::abs(left_residual) <= std::abs(right_residual)) {
                    classification_state = left_state;
                } else {
                    classification_state = right_state;
                }
                if (certified) {
                    return true;
                }
                if (allow_representable_bracket) {
                    root_log_volume = std::abs(left_residual)
                            <= std::abs(right_residual)
                        ? left_log_volume
                        : right_log_volume;
                    root_state = classification_state;
                    return true;
                }
            }
            const double trial_log_volume = 0.5 * (
                left_log_volume + right_log_volume
            );
            if (!(trial_log_volume > left_log_volume)
                || !(trial_log_volume < right_log_volume)) {
                if (std::abs(left_residual) <= std::abs(right_residual)) {
                    classification_state = left_state;
                } else {
                    classification_state = right_state;
                }
                if (certified) {
                    return true;
                }
                if (allow_representable_bracket) {
                    root_log_volume = std::abs(left_residual)
                            <= std::abs(right_residual)
                        ? left_log_volume
                        : right_log_volume;
                    root_state = classification_state;
                    return true;
                }
                ++result.reference_refinement_failure_count;
                return false;
            }
            Held2StateEvaluation trial_state;
            try {
                trial_state = log_volume_evaluator(
                    feed_independent, trial_log_volume
                );
            } catch (const std::exception&) {
                ++result.reference_evaluation_failure_count;
                return false;
            }
            const double root_residual = residual(trial_state);
            if (!certified && std::abs(root_residual) <= tolerance) {
                root_log_volume = trial_log_volume;
                root_state = trial_state;
                certified = true;
                if (!continue_after_certificate) {
                    classification_state = root_state;
                    return true;
                }
            }
            if (brackets_zero(left_residual, root_residual)) {
                right_log_volume = trial_log_volume;
                right_residual = root_residual;
                right_state = std::move(trial_state);
            } else {
                left_log_volume = trial_log_volume;
                left_residual = root_residual;
                left_state = std::move(trial_state);
            }
        }
        ++result.reference_refinement_failure_count;
        return false;
    };
    struct ReferenceCandidate {
        double log_volume = 0.0;
        Held2StateEvaluation state;
        Held2StateEvaluation classification_state;
        bool domain_boundary = false;
    };
    std::vector<ReferenceCandidate> reference_candidates;
    const auto add_reference_candidate = [
        &reference_candidates,
        reference_log_lower,
        reference_log_upper
    ](
        double log_volume,
        Held2StateEvaluation state,
        Held2StateEvaluation classification_state
    ) {
        const bool domain_boundary = log_volume <= reference_log_lower
            || log_volume >= reference_log_upper;
        const auto duplicate = std::find_if(
            reference_candidates.begin(),
            reference_candidates.end(),
            [log_volume](const ReferenceCandidate& candidate) {
                return std::abs(candidate.log_volume - log_volume) < 1.0e-8;
            }
        );
        if (duplicate != reference_candidates.end()) {
            duplicate->domain_boundary = duplicate->domain_boundary
                || domain_boundary;
            return;
        }
        reference_candidates.push_back({
            log_volume,
            std::move(state),
            std::move(classification_state),
            domain_boundary,
        });
    };
    std::vector<double> stationary_points;
    for (std::size_t interval = 0; interval + 1 < reference_scan.size(); ++interval) {
        const double left_log_volume = reference_scan[interval].first;
        const Held2StateEvaluation& left_state = reference_scan[interval].second;
        const double right_log_volume = reference_scan[interval + 1].first;
        const Held2StateEvaluation& right_state = reference_scan[interval + 1].second;
        const auto pressure_residual = [](const Held2StateEvaluation& state) {
            return state.pressure_stationarity_relative;
        };
        if (brackets_zero(
                pressure_residual(left_state), pressure_residual(right_state)
            )) {
            double root_log_volume = 0.0;
            Held2StateEvaluation root_state;
            Held2StateEvaluation classification_state;
            if (refine_bracket(
                    left_log_volume,
                    left_state,
                    right_log_volume,
                    right_state,
                    pressure_residual,
                    kCertificateTolerance,
                    true,
                    false,
                    root_log_volume,
                    root_state,
                    classification_state
                )) {
                add_reference_candidate(
                    root_log_volume,
                    std::move(root_state),
                    std::move(classification_state)
                );
            }
        }

        const auto pressure_derivative = [](const Held2StateEvaluation& state) {
            return state.pressure_stationarity_derivative_log_volume;
        };
        const double left_derivative = pressure_derivative(left_state);
        const double right_derivative = pressure_derivative(right_state);
        if (!std::isfinite(left_derivative) || !std::isfinite(right_derivative)) {
            ++result.reference_refinement_failure_count;
            result.reference_failure_reason =
                "reference_stationary_derivative_incomplete";
            continue;
        }
        if (brackets_zero(left_derivative, right_derivative)) {
            double stationary_log_volume = 0.0;
            Held2StateEvaluation stationary_state;
            Held2StateEvaluation classification_state;
            if (!refine_bracket(
                    left_log_volume,
                    left_state,
                    right_log_volume,
                    right_state,
                    pressure_derivative,
                    kCoordinateTolerance,
                    false,
                    true,
                    stationary_log_volume,
                    stationary_state,
                    classification_state
                )) {
                if (result.reference_failure_reason.empty()) {
                    result.reference_failure_reason =
                        "reference_stationary_refinement_incomplete";
                }
                continue;
            }
            const bool duplicate_stationary = std::any_of(
                stationary_points.begin(),
                stationary_points.end(),
                [stationary_log_volume](double retained) {
                    return std::abs(retained - stationary_log_volume) < 1.0e-8;
                }
            );
            if (!duplicate_stationary) {
                stationary_points.push_back(stationary_log_volume);
                ++result.reference_stationary_point_count;
                if (std::abs(stationary_state.pressure_stationarity_relative)
                    <= kCertificateTolerance) {
                    ++result.reference_tangential_root_count;
                    add_reference_candidate(
                        stationary_log_volume,
                        std::move(stationary_state),
                        std::move(classification_state)
                    );
                }
            }
        }
    }
    std::sort(
        reference_candidates.begin(),
        reference_candidates.end(),
        [](const ReferenceCandidate& left, const ReferenceCandidate& right) {
            return left.log_volume < right.log_volume;
        }
    );
    std::vector<std::size_t> stable_candidates;
    for (std::size_t index = 0; index < reference_candidates.size(); ++index) {
        const ReferenceCandidate& candidate = reference_candidates[index];
        const double curvature = candidate.classification_state.hessian.back();
        const bool marginal = !std::isfinite(curvature)
            || std::abs(curvature) <= kCertificateTolerance;
        const bool mechanically_stable = !marginal && curvature > 0.0;
        if (candidate.domain_boundary) {
            ++result.reference_boundary_root_count;
        }
        if (marginal) {
            ++result.reference_marginal_root_count;
        }
        result.reference_roots.push_back({
            std::log(candidate.classification_state.volume),
            candidate.classification_state.volume,
            candidate.classification_state.objective,
            candidate.classification_state.pressure_stationarity_relative,
            curvature,
            mechanically_stable,
        });
        if (mechanically_stable) {
            stable_candidates.push_back(index);
            ++result.reference_stable_root_count;
        }
    }
    result.reference_root_count = static_cast<int>(result.reference_roots.size());
    Held2StateEvaluation reference;
    std::vector<double> reference_variables;
    if (!stable_candidates.empty()) {
        double lowest_objective = std::numeric_limits<double>::infinity();
        for (std::size_t index : stable_candidates) {
            lowest_objective = std::min(
                lowest_objective,
                reference_candidates[index].classification_state.objective
            );
        }
        std::vector<std::size_t> lowest_candidates;
        for (std::size_t index : stable_candidates) {
            if (std::abs(
                    reference_candidates[index].classification_state.objective
                        - lowest_objective
                ) <= kCertificateTolerance) {
                lowest_candidates.push_back(index);
            }
        }
        if (lowest_candidates.empty()) {
            result.reference_failure_reason = "reference_objective_incomplete";
        } else if (lowest_candidates.size() > 1) {
            result.reference_objective_tie_count = static_cast<int>(
                lowest_candidates.size() - 1
            );
        } else {
            const ReferenceCandidate& selected =
                reference_candidates[lowest_candidates.front()];
            reference = selected.state;
            reference_variables = feed_independent;
            reference_variables.push_back(selected.log_volume);
        }
    }
    if (result.reference_failure_reason.empty()
        && (result.reference_evaluation_failure_count != 0
            || result.reference_refinement_failure_count != 0)) {
        result.reference_failure_reason = "reference_root_refinement_incomplete";
    }
    if (result.reference_failure_reason.empty()
        && result.reference_boundary_root_count != 0) {
        result.reference_failure_reason = "reference_root_at_domain_boundary";
    }
    if (result.reference_failure_reason.empty()
        && result.reference_marginal_root_count != 0) {
        result.reference_failure_reason = "reference_root_marginal";
    }
    if (result.reference_failure_reason.empty()
        && result.reference_objective_tie_count != 0) {
        result.reference_failure_reason = "reference_objective_tie";
    }
    if (result.reference_failure_reason.empty()
        && result.reference_root_count == 0) {
        result.reference_failure_reason = "reference_root_not_found";
    }
    if (result.reference_failure_reason.empty()
        && result.reference_stable_root_count == 0) {
        result.reference_failure_reason = "reference_stable_root_not_found";
    }
    if (result.reference_evaluation_failure_count != 0
        || result.reference_refinement_failure_count != 0
        || result.reference_boundary_root_count != 0
        || result.reference_marginal_root_count != 0
        || result.reference_objective_tie_count != 0
        || result.reference_root_count == 0
        || result.reference_stable_root_count == 0
        || reference_variables.empty()) {
        result.outcome = "indeterminate";
        return result;
    }
    result.reference_modified_fractions = reference.modified_fractions;
    result.reference_volume = reference.volume;
    if (search_uses_log_packing) {
        if (!reference.has_packing_evaluation || reference.packing.value <= 0.0) {
            result.outcome = "indeterminate";
            return result;
        }
        reference_variables = feed_independent;
        reference_variables.push_back(std::log(reference.packing.value));
        try {
            reference = search_evaluator(
                feed_independent, reference_variables.back()
            );
        } catch (const std::exception&) {
            result.outcome = "indeterminate";
            return result;
        }
    }

    std::vector<double> lower = coordinates.independent_lower_bounds;
    std::vector<double> upper = coordinates.independent_upper_bounds;
    for (std::size_t index = 0; index < upper.size(); ++index) {
        upper[index] = std::nextafter(upper[index], lower[index]);
    }
    const auto strict_log_packing_bounds = strict_interior_log_bounds({
        kHeld2PackingFractionMinimum,
        kHeld2PackingFractionMaximum,
    });
    lower.push_back(
        search_uses_log_packing
            ? strict_log_packing_bounds[0]
            : strict_log_volume_bounds[0]
    );
    upper.push_back(
        search_uses_log_packing
            ? strict_log_packing_bounds[1]
            : strict_log_volume_bounds[1]
    );
    std::mt19937 generator(kStageISeed);
    const auto dependent_retained = static_cast<std::size_t>(
        std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            coordinates.dependent_index
        ) - coordinates.retained_indices.begin()
    );
    const double composition_sum_upper =
        1.0
        - kHeld2ModifiedLowerScale * coordinates.modified_factors[dependent_retained];
    std::vector<std::vector<double>> starts;
    starts.reserve(static_cast<std::size_t>(result.declared_start_count));
    if (feed_independent.size() == 1 && molar_volume_bounds[0] <= 1.0
        && molar_volume_bounds[1] >= 1.0) {
        for (double composition : {0.2, 0.8}) {
            const std::vector<double> fixed{composition};
            if (!std::isfinite(total_ion_mole_fraction_max)
                || held2_stage_i_total_ion_mole_fraction(coordinates, fixed)
                    < total_ion_mole_fraction_max) {
                starts.push_back({composition, 0.0});
            }
        }
    }
    int rejected_draws = 0;
    while (starts.size() < static_cast<std::size_t>(result.declared_start_count)) {
        std::vector<double> start;
        start.reserve(feed_independent.size() + 1);
        double sum = 0.0;
        for (std::size_t index = 0; index < feed_independent.size(); ++index) {
            std::uniform_real_distribution<double> distribution(lower[index], upper[index]);
            start.push_back(distribution(generator));
            sum += start.back();
        }
        if (!(sum <= composition_sum_upper) && ++rejected_draws < 10000) {
            continue;
        }
        if (!(sum <= composition_sum_upper)) {
            result.outcome = "indeterminate";
            return result;
        }
        if (std::isfinite(total_ion_mole_fraction_max)
            && !(held2_stage_i_total_ion_mole_fraction(coordinates, start)
                < total_ion_mole_fraction_max)) {
            if (++rejected_draws < 10000) {
                continue;
            }
            result.outcome = "indeterminate";
            return result;
        }
        std::array<double, 2> start_volume_bounds = molar_volume_bounds;
        if (volume_bounds_evaluator) {
            try {
                start_volume_bounds = volume_bounds_evaluator(
                    held2_lift_independent_fractions(coordinates, start)
                );
            } catch (const std::exception&) {
                if (++rejected_draws < 10000) {
                    continue;
                }
                result.outcome = "indeterminate";
                return result;
            }
        }
        const auto start_log_volume_bounds = strict_interior_log_bounds(
            start_volume_bounds
        );
        std::uniform_real_distribution<double> log_volume_distribution(
            start_log_volume_bounds[0], start_log_volume_bounds[1]
        );
        const double start_log_volume = log_volume_distribution(generator);
        if (search_uses_log_packing) {
            try {
                const Held2StateEvaluation start_state = log_volume_evaluator(
                    start, start_log_volume
                );
                if (!start_state.has_packing_evaluation
                    || start_state.packing.value <= 0.0) {
                    throw std::invalid_argument("HELD2 start packing evidence is incomplete");
                }
                start.push_back(std::log(start_state.packing.value));
            } catch (const std::exception&) {
                if (++rejected_draws < 10000) {
                    continue;
                }
                result.outcome = "indeterminate";
                return result;
            }
        } else {
            start.push_back(start_log_volume);
        }
        starts.push_back(std::move(start));
    }
    result.planned_starts = starts;
    Held2SearchDomainPredicate trial_domain_predicate;
    if (!search_uses_log_packing) {
        trial_domain_predicate = [
            coordinates,
            molar_volume_bounds,
            volume_bounds_evaluator
        ](const std::vector<double>& independent, double log_volume) {
            std::array<double, 2> bounds = molar_volume_bounds;
            if (volume_bounds_evaluator) {
                bounds = volume_bounds_evaluator(
                    held2_lift_independent_fractions(coordinates, independent)
                );
            }
            const double volume = std::exp(log_volume);
            return std::isfinite(volume) && volume > bounds[0] && volume < bounds[1];
        };
    }
    const std::size_t attempted_limit = solve_start_limit > 0
        ? std::min(starts.size(), static_cast<std::size_t>(solve_start_limit))
        : starts.size();
    for (std::size_t start_index = 0; start_index < attempted_limit; ++start_index) {
        const std::vector<double>& start = starts[start_index];
        ++result.attempted_start_count;
        const Held2SearchRun run = solve_held2_search(
            search_evaluator,
            reference_variables,
            reference,
            true,
            start,
            lower,
            upper,
            composition_sum_upper,
            source_domain_coefficients,
            total_ion_mole_fraction_max,
            trial_domain_predicate
        );
        if (!run.solver_converged || !run.callback_error.empty()
            || run.variables.size() != feed_independent.size() + 1) {
            if (result.failed_start_index < 0) {
                result.failed_start_index = static_cast<int>(start_index);
                result.failed_start_solver_status = run.solver_status;
                result.failed_start_solver_converged = run.solver_converged;
                result.failed_start_reason = run.callback_error.empty()
                    ? "TPD solve did not return a complete accepted state"
                    : run.callback_error;
                result.failed_start_initial = start;
            }
            ++result.failed_start_count;
            continue;
        }
        ++result.completed_start_count;
        const std::vector<double> composition(
            run.variables.begin(),
            run.variables.end() - 1
        );
        Held2StateEvaluation state;
        try {
            state = search_evaluator(composition, run.variables.back());
        } catch (const std::exception& error) {
            --result.completed_start_count;
            if (result.failed_start_index < 0) {
                result.failed_start_index = static_cast<int>(start_index);
                result.failed_start_solver_status = run.solver_status;
                result.failed_start_solver_converged = run.solver_converged;
                result.failed_start_reason = error.what();
                result.failed_start_initial = start;
            }
            ++result.failed_start_count;
            continue;
        }
        double tpd = state.objective - reference.objective;
        for (std::size_t index = 0; index < state.gradient.size(); ++index) {
            tpd -= reference.gradient[index]
                * (run.variables[index] - reference_variables[index]);
        }
        result.minimum_tpd = std::min(result.minimum_tpd, tpd);
        if (tpd >= kStageITpdThreshold) {
            continue;
        }
        std::array<double, 2> candidate_volume_bounds = molar_volume_bounds;
        if (volume_bounds_evaluator) {
            try {
                candidate_volume_bounds = volume_bounds_evaluator(state.physical_amounts);
            } catch (const std::exception&) {
                ++result.candidate_domain_evaluation_failure_count;
                continue;
            }
        }
        const double lower_relative_distance = std::abs(
            state.volume - candidate_volume_bounds[0]
        ) / std::max(state.volume, candidate_volume_bounds[0]);
        const double upper_relative_distance = std::abs(
            state.volume - candidate_volume_bounds[1]
        ) / std::max(state.volume, candidate_volume_bounds[1]);
        const bool lower_active = lower_relative_distance <= kCertificateTolerance;
        const bool upper_active = upper_relative_distance <= kCertificateTolerance;
        const double volume_gradient = search_uses_log_packing
            ? state.log_volume_gradient
            : state.gradient.back();
        const bool volume_in_domain = state.volume >= candidate_volume_bounds[0]
            && state.volume <= candidate_volume_bounds[1];
        const bool volume_stationary =
            std::abs(state.pressure_stationarity_relative) <= kCertificateTolerance
            || (lower_active && volume_gradient >= -kCertificateTolerance)
            || (upper_active && volume_gradient <= kCertificateTolerance);
        if (!volume_in_domain || !volume_stationary) {
            ++result.candidate_domain_rejection_count;
            continue;
        }
        const bool duplicate = std::any_of(
            result.candidates.begin(),
            result.candidates.end(),
            [&state](const Held2StageICandidate& candidate) {
                return maximum_abs_difference(
                           candidate.modified_fractions,
                           state.modified_fractions
                       ) < 1.0e-7
                    && std::abs(candidate.volume - state.volume) < 1.0e-7;
            }
        );
        if (!duplicate) {
            result.candidates.push_back({
                state.modified_fractions,
                state.volume,
                tpd,
                candidate_volume_bounds,
                state.pressure_stationarity_relative,
                volume_gradient,
                state.has_packing_evaluation ? state.packing.value : 0.0,
                lower_active,
                upper_active,
            });
        }
    }
    std::sort(
        result.candidates.begin(),
        result.candidates.end(),
        [](const Held2StageICandidate& left, const Held2StageICandidate& right) {
            return left.modified_fractions < right.modified_fractions;
        }
    );
    if (!result.candidates.empty()) {
        result.outcome = "negative_tpd";
    } else if (result.failed_start_count != 0
        || result.candidate_domain_evaluation_failure_count != 0
        || !result.volume_domain_search_complete) {
        result.outcome = "indeterminate";
    } else {
        result.outcome = "no_negative_found";
    }
    return result;
}

Held2StageIResult solve_held2_manufactured_stage_i(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::string& reference_scenario
) {
    const Held2Coordinates coordinates = make_held2_coordinates(charges);
    if (coordinates.independent_indices.size() != 1) {
        throw std::invalid_argument(
            "manufactured HELD2 Stage I requires one independent modified composition"
        );
    }
    std::array<double, 2> volume_bounds{0.5, 1.5};
    Held2StateEvaluator evaluator;
    if (reference_scenario == "stage_i") {
        evaluator = [coordinates](
            const std::vector<double>& independent,
            double log_volume
        ) {
            return evaluate_manufactured_state_impl(
                coordinates,
                independent,
                log_volume
            );
        };
    } else {
        if (reference_scenario != "reference_tangential"
            && reference_scenario != "reference_marginal"
            && reference_scenario != "reference_domain_boundary"
            && reference_scenario != "reference_objective_tie") {
            throw std::invalid_argument(
                "unsupported manufactured HELD2 reference scenario"
            );
        }
        volume_bounds = {std::exp(-1.0), std::exp(1.0)};
        const auto strict_bounds = strict_interior_log_bounds(volume_bounds);
        const double boundary_root = strict_bounds[0];
        evaluator = [
            coordinates,
            reference_scenario,
            boundary_root
        ](
            const std::vector<double>& independent,
            double log_volume
        ) {
            Held2StateEvaluation state;
            state.physical_amounts = held2_lift_independent_fractions(
                coordinates, independent
            );
            state.modified_fractions = held2_transform_physical_fractions(
                coordinates, state.physical_amounts
            );
            state.volume = std::exp(log_volume);
            double gradient = 0.0;
            double curvature = 0.0;
            if (reference_scenario == "reference_tangential") {
                const double delta = log_volume - 0.137;
                state.objective = delta * delta * delta / 3.0;
                gradient = delta * delta;
                curvature = 2.0 * delta;
            } else if (reference_scenario == "reference_marginal") {
                const double delta = log_volume - 0.137;
                constexpr double marginal_curvature = 0.5e-8;
                state.objective = 0.5 * marginal_curvature * delta * delta
                    + 0.25 * std::pow(delta, 4);
                gradient = marginal_curvature * delta
                    + delta * delta * delta;
                curvature = marginal_curvature + 3.0 * delta * delta;
            } else if (reference_scenario == "reference_domain_boundary") {
                const double delta = log_volume - boundary_root;
                state.objective = 0.5 * delta * delta;
                gradient = delta;
                curvature = 1.0;
            } else {
                constexpr double root = 0.4;
                const double difference = log_volume * log_volume - root * root;
                state.objective = difference * difference;
                gradient = 4.0 * log_volume * difference;
                curvature = 12.0 * log_volume * log_volume - 4.0 * root * root;
            }
            state.gradient = {0.0, gradient};
            state.hessian = {1.0, 0.0, 0.0, curvature};
            state.modified_potentials.assign(
                coordinates.retained_indices.size(), 0.0
            );
            state.pressure_stationarity_relative = gradient;
            state.pressure_stationarity_derivative_log_volume = curvature;
            state.log_volume_gradient = gradient;
            return state;
        };
    }
    const Held2VolumeBoundsEvaluator volume_bounds_evaluator = [volume_bounds](
        const std::vector<double>&
    ) {
        return volume_bounds;
    };
    return solve_held2_stage_i(
        coordinates,
        physical_feed,
        evaluator,
        evaluator,
        volume_bounds,
        volume_bounds_evaluator,
        false,
        true,
        -1,
        std::numeric_limits<double>::quiet_NaN()
    );
}

namespace {

struct Held2StageIIStep6Cut {
    Held2StateEvaluation state;
    std::vector<double> independent_modified_fractions;
    std::vector<double> fixed_volume_composition_gradient;
    double phase_coordinate = 0.0;
};

[[nodiscard]] std::vector<Held2StageIICandidate>
select_held2_stage_ii_candidates(
    const Held2Coordinates& coordinates,
    const std::vector<double>& feed_independent_modified_fractions,
    double upper_bound,
    const std::vector<double>& multipliers,
    const std::vector<Held2StageIIStep6Cut>& cuts
);

}  // namespace

Held2StageIIResult solve_held2_manufactured_stage_ii(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed
) {
    constexpr int major_iteration_cap = 100;
    constexpr double stage_ii_tolerance = 1.0e-8;
    const Held2Coordinates coordinates = make_held2_coordinates(charges);
    if (coordinates.independent_indices.size() != 1) {
        throw std::invalid_argument(
            "manufactured HELD2 Stage II requires one independent modified composition"
        );
    }
    const std::vector<double> modified_feed = held2_transform_physical_fractions(
        coordinates,
        physical_feed
    );
    const auto independent_retained = static_cast<std::size_t>(
        std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            coordinates.independent_indices.front()
        ) - coordinates.retained_indices.begin()
    );
    const double feed = modified_feed[independent_retained];
    const Held2StateEvaluator evaluator = [coordinates](
        const std::vector<double>& independent,
        double log_volume
    ) {
        return evaluate_manufactured_state_impl(
            coordinates,
            independent,
            log_volume
        );
    };
    std::vector<Held2StateEvaluation> cuts = {
        evaluator({feed}, 0.0),
        evaluator({0.5 * feed}, 0.0),
        evaluator({0.5 * (1.0 + feed)}, 0.0),
    };
    const double feed_gibbs = cuts.front().objective;
    const double lower_composition = coordinates.independent_lower_bounds.front();
    const double upper_composition = std::nextafter(
        coordinates.independent_upper_bounds.front(),
        lower_composition
    );
    const std::vector<double> lower = {lower_composition, std::log(0.5)};
    const std::vector<double> upper = {upper_composition, std::log(1.5)};
    Held2StageIIResult result;
    result.lower_starts_per_iteration = 10 * static_cast<int>(charges.size());

    for (int major = 0; major < major_iteration_cap; ++major) {
        struct Line {
            double intercept = 0.0;
            double slope = 0.0;
        };
        std::vector<Line> lines;
        lines.reserve(cuts.size() + 1);
        for (const Held2StateEvaluation& cut : cuts) {
            lines.push_back({
                cut.objective,
                feed - cut.modified_fractions[independent_retained],
            });
        }
        lines.push_back({feed_gibbs, 0.0});
        std::vector<double> multiplier_candidates = {0.0};
        for (std::size_t left = 0; left < lines.size(); ++left) {
            for (std::size_t right = left + 1; right < lines.size(); ++right) {
                const double slope_delta = lines[left].slope - lines[right].slope;
                if (std::abs(slope_delta) <= std::numeric_limits<double>::epsilon()) {
                    continue;
                }
                multiplier_candidates.push_back(
                    (lines[right].intercept - lines[left].intercept) / slope_delta
                );
            }
        }
        double upper_bound = -std::numeric_limits<double>::infinity();
        double multiplier = 0.0;
        for (double candidate : multiplier_candidates) {
            double envelope = std::numeric_limits<double>::infinity();
            for (const Line& line : lines) {
                envelope = std::min(envelope, line.intercept + line.slope * candidate);
            }
            if (envelope > upper_bound) {
                upper_bound = envelope;
                multiplier = candidate;
            }
        }

        Held2StateEvaluation lower_reference;
        lower_reference.gradient = {multiplier, 0.0};
        lower_reference.hessian.assign(4, 0.0);
        const std::vector<double> lower_reference_variables = {feed, 0.0};
        std::mt19937 generator(kStageISeed + static_cast<unsigned int>(major));
        std::uniform_real_distribution<double> composition_distribution(
            lower_composition,
            upper_composition
        );
        std::uniform_real_distribution<double> log_volume_distribution(
            lower.back(),
            upper.back()
        );
        std::vector<std::vector<double>> starts = {{0.2, 0.0}, {0.8, 0.0}};
        while (starts.size() < static_cast<std::size_t>(result.lower_starts_per_iteration)) {
            starts.push_back({
                composition_distribution(generator),
                log_volume_distribution(generator),
            });
        }
        double lower_bound = std::numeric_limits<double>::infinity();
        Held2StateEvaluation best_state;
        bool lower_failed = false;
        for (std::size_t start_index = 0; start_index < starts.size(); ++start_index) {
            const std::vector<double>& start = starts[start_index];
            ++result.lower_attempted_start_count;
            const Held2SearchRun run = solve_held2_search(
                evaluator,
                lower_reference_variables,
                lower_reference,
                true,
                start,
                lower,
                upper
            );
            if (!run.solver_converged || !run.callback_error.empty()
                || run.variables.size() != 2) {
                lower_failed = true;
                continue;
            }
            Held2StateEvaluation state = evaluator(
                {run.variables.front()},
                run.variables.back()
            );
            const double value = state.objective
                + multiplier * (feed - state.modified_fractions[independent_retained]);
            if (value < lower_bound) {
                lower_bound = value;
                best_state = std::move(state);
            }
        }
        ++result.major_iterations;
        if (lower_failed || !std::isfinite(lower_bound)) {
            result.outcome = "indeterminate";
            result.cut_count = static_cast<int>(cuts.size());
            return result;
        }
        result.bound_history.push_back({
            lower_bound,
            upper_bound,
            multiplier,
            static_cast<int>(cuts.size()),
        });
        const bool duplicate = std::any_of(
            cuts.begin(),
            cuts.end(),
            [&best_state](const Held2StateEvaluation& cut) {
                return maximum_abs_difference(
                           cut.modified_fractions,
                           best_state.modified_fractions
                       ) < 1.0e-7
                    && std::abs(cut.volume - best_state.volume) < 1.0e-7;
            }
        );
        if (!duplicate) {
            cuts.push_back(best_state);
        }
        result.candidates.clear();
        for (const Held2StateEvaluation& cut : cuts) {
            const double lower_value = cut.objective
                + multiplier * (feed - cut.modified_fractions[independent_retained]);
            const double gap = upper_bound - lower_value;
            if (gap >= -stage_ii_tolerance && gap <= stage_ii_tolerance
                && std::abs(cut.gradient.front() - multiplier) <= stage_ii_tolerance
                && std::abs(cut.gradient.back()) <= stage_ii_tolerance) {
                result.candidates.push_back({
                    cut.modified_fractions,
                    {cut.modified_fractions[independent_retained]},
                    cut.volume,
                    std::log(cut.volume),
                    gap,
                });
            }
        }
        if (result.candidates.size() > 1) {
            result.outcome = "candidate_set";
            result.cut_count = static_cast<int>(cuts.size());
            return result;
        }
        if (duplicate) {
            result.outcome = "no_progress";
            result.cut_count = static_cast<int>(cuts.size());
            return result;
        }
    }
    result.outcome = "resource_limit";
    result.cut_count = static_cast<int>(cuts.size());
    return result;
}

namespace {

std::vector<Held2StageIICandidate> select_held2_stage_ii_candidates(
    const Held2Coordinates& coordinates,
    const std::vector<double>& feed_independent_modified_fractions,
    double upper_bound,
    const std::vector<double>& multipliers,
    const std::vector<Held2StageIIStep6Cut>& cuts
) {
    const std::size_t dimension = coordinates.independent_indices.size();
    if (!std::isfinite(upper_bound)
        || feed_independent_modified_fractions.size() != dimension
        || multipliers.size() != dimension
        || coordinates.independent_lower_bounds.size() != dimension) {
        throw std::invalid_argument("HELD2 Stage II Step 6 evidence is incomplete");
    }
    require_finite_vector(
        feed_independent_modified_fractions,
        "HELD2 Stage II Step 6 feed"
    );
    require_finite_vector(multipliers, "HELD2 Stage II Step 6 multipliers");

    const std::size_t component_count_to_check = std::min(
        dimension,
        coordinates.charges.size() - 2
    );
    std::vector<Held2StageIICandidate> result;
    std::vector<std::size_t> selected_cut_indices;
    for (std::size_t cut_index = 0; cut_index < cuts.size(); ++cut_index) {
        const Held2StageIIStep6Cut& cut = cuts[cut_index];
        if (cut.independent_modified_fractions.size() != dimension
            || cut.fixed_volume_composition_gradient.size() != dimension
            || cut.state.gradient.size() != dimension + 1
            || cut.state.modified_fractions.size()
                != coordinates.retained_indices.size()
            || !cut.state.has_packing_evaluation
            || !std::isfinite(cut.state.objective)
            || !std::isfinite(cut.state.volume)
            || !std::isfinite(cut.state.packing.value)) {
            throw std::invalid_argument(
                "HELD2 Stage II Step 6 cut evidence is incomplete"
            );
        }
        require_finite_vector(
            cut.independent_modified_fractions,
            "HELD2 Stage II Step 6 modified fractions"
        );
        require_finite_vector(
            cut.fixed_volume_composition_gradient,
            "HELD2 Stage II Step 6 fixed-volume gradient"
        );
        require_finite_vector(cut.state.gradient, "HELD2 Stage II Step 6 q gradient");

        double lower_value = cut.state.objective;
        for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
            lower_value += multipliers[coordinate]
                * (feed_independent_modified_fractions[coordinate]
                   - cut.independent_modified_fractions[coordinate]);
        }
        const double gap = upper_bound - lower_value;
        if (std::abs(gap) > kStageIIStep6EpsilonB
            || std::abs(cut.state.gradient.back())
                > kStageIIStep6QStationarityTolerance) {
            continue;
        }

        std::vector<std::size_t> checked_component_indices;
        bool fixed_volume_stationary = true;
        for (std::size_t coordinate = 0;
             coordinate < dimension
             && checked_component_indices.size() < component_count_to_check;
             ++coordinate) {
            if (cut.independent_modified_fractions[coordinate]
                <= coordinates.independent_lower_bounds[coordinate]) {
                continue;
            }
            checked_component_indices.push_back(coordinate);
            if (std::abs(
                    cut.fixed_volume_composition_gradient[coordinate]
                    - multipliers[coordinate]
                ) > kStageIIStep6EpsilonLambda * std::abs(multipliers[coordinate])) {
                fixed_volume_stationary = false;
                break;
            }
        }
        if (!fixed_volume_stationary) {
            continue;
        }

        const bool distinct = std::all_of(
            selected_cut_indices.begin(),
            selected_cut_indices.end(),
            [&cuts, &cut](std::size_t selected_index) {
                const Held2StageIIStep6Cut& selected = cuts[selected_index];
                return std::abs(cut.state.packing.value - selected.state.packing.value)
                        >= kStageIIStep6EpsilonEta
                    || maximum_abs_difference(
                           cut.independent_modified_fractions,
                           selected.independent_modified_fractions
                       ) >= kStageIIStep6EpsilonX;
            }
        );
        if (!distinct) {
            continue;
        }
        selected_cut_indices.push_back(cut_index);
        result.push_back({
            cut.state.modified_fractions,
            cut.independent_modified_fractions,
            cut.state.volume,
            cut.phase_coordinate,
            gap,
        });
    }
    return result;
}

}  // namespace

std::vector<Held2StageIICandidate> evaluate_held2_stage_ii_step6_test_adapter(
    const Held2Coordinates& coordinates,
    const std::vector<double>& feed_independent_modified_fractions,
    double upper_bound,
    const std::vector<double>& multipliers,
    const std::vector<Held2StateEvaluation>& states,
    const std::vector<std::vector<double>>& independent_modified_fractions,
    const std::vector<std::vector<double>>& fixed_volume_composition_gradients,
    const std::vector<double>& phase_coordinates
) {
    if (states.size() != independent_modified_fractions.size()
        || states.size() != fixed_volume_composition_gradients.size()
        || states.size() != phase_coordinates.size()) {
        throw std::invalid_argument("HELD2 Stage II Step 6 test evidence is incomplete");
    }
    std::vector<Held2StageIIStep6Cut> cuts;
    cuts.reserve(states.size());
    for (std::size_t index = 0; index < states.size(); ++index) {
        cuts.push_back({
            states[index],
            independent_modified_fractions[index],
            fixed_volume_composition_gradients[index],
            phase_coordinates[index],
        });
    }
    return select_held2_stage_ii_candidates(
        coordinates,
        feed_independent_modified_fractions,
        upper_bound,
        multipliers,
        cuts
    );
}

Held2StageIIResult solve_held2_stage_ii(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const Held2StateEvaluator& search_evaluator,
    const Held2StateEvaluator& fixed_volume_evaluator,
    const Held2StateEvaluation& reference,
    const std::vector<Held2StageICandidate>& stage_i_candidates
) {
    constexpr int major_iteration_cap = 100;
    constexpr double stage_ii_tolerance = 1.0e-8;
    const std::size_t dimension = coordinates.independent_indices.size();
    const std::vector<double> modified_feed = held2_transform_physical_fractions(
        coordinates, physical_feed
    );
    std::vector<std::size_t> retained_positions;
    std::vector<double> feed;
    for (std::size_t component : coordinates.independent_indices) {
        const auto position = std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            component
        );
        if (position == coordinates.retained_indices.end()) {
            throw std::invalid_argument("HELD2 Stage II independent coordinate is not retained");
        }
        retained_positions.push_back(static_cast<std::size_t>(
            position - coordinates.retained_indices.begin()
        ));
        feed.push_back(modified_feed[retained_positions.back()]);
    }
    if (dimension == 0 || reference.packing.value <= 0.0) {
        throw std::invalid_argument("HELD2 Stage II reference evidence is incomplete");
    }
    const auto make_cut = [
        &fixed_volume_evaluator,
        dimension
    ](
        Held2StateEvaluation state,
        std::vector<double> independent,
        double phase_coordinate
    ) {
        if (!(state.volume > 0.0) || !std::isfinite(state.volume)) {
            throw std::invalid_argument(
                "HELD2 Stage II cut has an invalid phase volume"
            );
        }
        const Held2StateEvaluation fixed_volume_state = fixed_volume_evaluator(
            independent,
            std::log(state.volume)
        );
        if (fixed_volume_state.gradient.size() < dimension) {
            throw std::invalid_argument(
                "HELD2 Stage II fixed-volume gradient is incomplete"
            );
        }
        std::vector<double> fixed_volume_composition_gradient(
            fixed_volume_state.gradient.begin(),
            fixed_volume_state.gradient.begin()
                + static_cast<std::ptrdiff_t>(dimension)
        );
        return Held2StageIIStep6Cut{
            std::move(state),
            std::move(independent),
            std::move(fixed_volume_composition_gradient),
            phase_coordinate,
        };
    };
    std::vector<Held2StageIIStep6Cut> cuts;
    cuts.push_back(make_cut(reference, feed, std::log(reference.packing.value)));
    for (const Held2StageICandidate& candidate : stage_i_candidates) {
        std::vector<double> independent;
        for (std::size_t position : retained_positions) {
            independent.push_back(candidate.modified_fractions[position]);
        }
        if (candidate.packing_fraction <= 0.0) {
            continue;
        }
        const double coordinate = std::log(candidate.packing_fraction);
        cuts.push_back(make_cut(
            search_evaluator(independent, coordinate),
            independent,
            coordinate
        ));
    }
    const double composition_sum_upper = 1.0 - kHeld2ModifiedLowerScale
        * coordinates.modified_factors[static_cast<std::size_t>(
            std::find(
                coordinates.retained_indices.begin(),
                coordinates.retained_indices.end(),
                coordinates.dependent_index
            ) - coordinates.retained_indices.begin()
        )];
    for (std::size_t coordinate = 0; cuts.size() < dimension + 1
         && coordinate < dimension; ++coordinate) {
        std::vector<double> point = feed;
        const double span = coordinates.independent_upper_bounds[coordinate]
            - coordinates.independent_lower_bounds[coordinate];
        double delta = 0.05 * span;
        if (std::accumulate(point.begin(), point.end(), 0.0) + delta
            > composition_sum_upper) {
            delta = -delta;
        }
        point[coordinate] = std::clamp(
            point[coordinate] + delta,
            coordinates.independent_lower_bounds[coordinate],
            std::nextafter(
                coordinates.independent_upper_bounds[coordinate],
                coordinates.independent_lower_bounds[coordinate]
            )
        );
        cuts.push_back(make_cut(
            search_evaluator(point, cuts.front().phase_coordinate),
            point,
            cuts.front().phase_coordinate
        ));
    }
    Held2StageIIResult result;
    result.lower_starts_per_iteration = 10 * static_cast<int>(coordinates.charges.size());
    if (cuts.size() < dimension + 1) {
        result.outcome = "indeterminate";
        result.cut_count = static_cast<int>(cuts.size());
        return result;
    }

    const auto solve_dense = [](std::vector<std::vector<double>> matrix,
                                std::vector<double> rhs) {
        const std::size_t size = rhs.size();
        for (std::size_t pivot = 0; pivot < size; ++pivot) {
            std::size_t best = pivot;
            for (std::size_t row = pivot + 1; row < size; ++row) {
                if (std::abs(matrix[row][pivot]) > std::abs(matrix[best][pivot])) {
                    best = row;
                }
            }
            if (std::abs(matrix[best][pivot]) < 1.0e-13) {
                return std::vector<double>{};
            }
            std::swap(matrix[pivot], matrix[best]);
            std::swap(rhs[pivot], rhs[best]);
            const double diagonal = matrix[pivot][pivot];
            for (std::size_t column = pivot; column < size; ++column) {
                matrix[pivot][column] /= diagonal;
            }
            rhs[pivot] /= diagonal;
            for (std::size_t row = 0; row < size; ++row) {
                if (row == pivot) {
                    continue;
                }
                const double factor = matrix[row][pivot];
                for (std::size_t column = pivot; column < size; ++column) {
                    matrix[row][column] -= factor * matrix[pivot][column];
                }
                rhs[row] -= factor * rhs[pivot];
            }
        }
        return rhs;
    };

    for (int major = 0; major < major_iteration_cap; ++major) {
        double upper_bound = -std::numeric_limits<double>::infinity();
        std::vector<double> multiplier(dimension, 0.0);
        std::vector<std::size_t> selected;
        std::function<void(std::size_t)> enumerate = [&](std::size_t next) {
            if (selected.size() == dimension + 1) {
                std::vector<std::vector<double>> matrix(
                    dimension + 1, std::vector<double>(dimension + 1, 0.0)
                );
                std::vector<double> rhs(dimension + 1, 0.0);
                for (std::size_t row = 0; row < selected.size(); ++row) {
                    const Held2StageIIStep6Cut& cut = cuts[selected[row]];
                    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                        matrix[row][coordinate] = -(
                            feed[coordinate]
                            - cut.independent_modified_fractions[coordinate]
                        );
                    }
                    matrix[row][dimension] = 1.0;
                    rhs[row] = cut.state.objective;
                }
                const std::vector<double> solution = solve_dense(matrix, rhs);
                if (solution.empty()) {
                    return;
                }
                const double envelope = solution.back();
                for (const Held2StageIIStep6Cut& cut : cuts) {
                    double line = cut.state.objective;
                    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                        line += solution[coordinate]
                            * (feed[coordinate]
                               - cut.independent_modified_fractions[coordinate]);
                    }
                    if (envelope > line + stage_ii_tolerance) {
                        return;
                    }
                }
                if (envelope > upper_bound) {
                    upper_bound = envelope;
                    std::copy(
                        solution.begin(), solution.end() - 1, multiplier.begin()
                    );
                }
                return;
            }
            for (std::size_t index = next; index < cuts.size(); ++index) {
                selected.push_back(index);
                enumerate(index + 1);
                selected.pop_back();
            }
        };
        enumerate(0);
        if (!std::isfinite(upper_bound)) {
            result.outcome = "indeterminate";
            result.cut_count = static_cast<int>(cuts.size());
            return result;
        }

        std::vector<double> lower = coordinates.independent_lower_bounds;
        std::vector<double> upper = coordinates.independent_upper_bounds;
        for (std::size_t index = 0; index < upper.size(); ++index) {
            upper[index] = std::nextafter(upper[index], lower[index]);
        }
        lower.push_back(std::log(kHeld2PackingFractionMinimum));
        upper.push_back(std::log(kHeld2PackingFractionMaximum));
        std::vector<std::vector<double>> starts;
        for (const Held2StageIIStep6Cut& cut : cuts) {
            std::vector<double> start = cut.independent_modified_fractions;
            start.push_back(cut.phase_coordinate);
            starts.push_back(std::move(start));
        }
        std::mt19937 generator(kStageISeed + static_cast<unsigned int>(major));
        int rejected = 0;
        while (starts.size() < static_cast<std::size_t>(result.lower_starts_per_iteration)) {
            std::vector<double> start;
            double sum = 0.0;
            for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                std::uniform_real_distribution<double> distribution(
                    lower[coordinate], upper[coordinate]
                );
                start.push_back(distribution(generator));
                sum += start.back();
            }
            if (sum > composition_sum_upper && ++rejected < 10000) continue;
            if (sum > composition_sum_upper) {
                result.outcome = "indeterminate";
                result.cut_count = static_cast<int>(cuts.size());
                return result;
            }
            std::uniform_real_distribution<double> packing_distribution(
                lower.back(), upper.back()
            );
            start.push_back(packing_distribution(generator));
            starts.push_back(std::move(start));
        }
        double best_certified_value = std::numeric_limits<double>::infinity();
        Held2StageIIStep6Cut best;
        int certified_start_count = 0;
        for (std::size_t start_index = 0; start_index < starts.size(); ++start_index) {
            const std::vector<double>& start = starts[start_index];
            ++result.lower_attempted_start_count;
            const Held2SearchRun run = solve_held2_stage_ii_search(
                search_evaluator,
                feed,
                multiplier,
                start,
                lower,
                upper,
                composition_sum_upper
            );
            Held2StageIIAttempt attempt;
            attempt.start_index = static_cast<int>(start_index);
            attempt.solver_status = run.solver_status;
            attempt.iterations = run.iterations;
            attempt.solver_converged = run.solver_converged;
            attempt.callback_error = run.callback_error;
            attempt.final_step_norm = run.final_step_norm;
            attempt.lower_value = std::numeric_limits<double>::infinity();
            attempt.projected_kkt_inf_norm = std::numeric_limits<double>::infinity();
            attempt.constraint_violation = std::numeric_limits<double>::infinity();
            attempt.complementarity = std::numeric_limits<double>::infinity();
            attempt.dual_reconstruction_inf_norm =
                std::numeric_limits<double>::infinity();
            attempt.dual_signs_valid = false;
            Held2StateEvaluation state;
            std::vector<double> independent;
            if (run.variables.size() == dimension + 1) {
                try {
                    independent.assign(run.variables.begin(), run.variables.end() - 1);
                    state = search_evaluator(independent, run.variables.back());
                    attempt.provider_terminal_valid = true;
                    attempt.lower_value = state.objective;
                    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                        attempt.lower_value += multiplier[coordinate]
                            * (feed[coordinate] - independent[coordinate]);
                    }
                    attempt.constraint_violation = std::max(
                        0.0,
                        std::accumulate(independent.begin(), independent.end(), 0.0)
                            - composition_sum_upper
                    );
                    for (std::size_t index = 0; index < run.variables.size(); ++index) {
                        attempt.constraint_violation = std::max({
                            attempt.constraint_violation,
                            lower[index] - run.variables[index],
                            run.variables[index] - upper[index],
                        });
                    }
                    const bool packing_domain = state.has_packing_evaluation
                        && state.packing.value >= kHeld2PackingFractionMinimum
                        && state.packing.value <= kHeld2PackingFractionMaximum
                        && std::abs(
                            std::log(state.packing.value) - run.variables.back()
                        ) <= kCertificateTolerance;
                    const bool physical = std::abs(
                        std::accumulate(
                            state.modified_fractions.begin(),
                            state.modified_fractions.end(),
                            0.0
                        ) - 1.0
                    ) <= kCertificateTolerance
                        && std::abs(
                            std::accumulate(
                                state.physical_amounts.begin(),
                                state.physical_amounts.end(),
                                0.0
                            ) - 1.0
                        ) <= kCertificateTolerance
                        && std::abs(charge_residual(
                            coordinates.charges, state.physical_amounts
                        )) <= kCertificateTolerance;
                    const Held2PhysicalKkt kkt =
                        evaluate_held2_stage_ii_physical_kkt(
                            state.gradient,
                            multiplier,
                            run.variables,
                            lower,
                            upper,
                            composition_sum_upper,
                            run.coordinate_jacobian,
                            run.lower_bound_multipliers,
                            run.upper_bound_multipliers
                        );
                    attempt.projected_kkt_inf_norm = kkt.stationarity_inf_norm;
                    attempt.complementarity = kkt.complementarity;
                    attempt.dual_reconstruction_inf_norm =
                        kkt.reconstruction_inf_norm;
                    attempt.dual_signs_valid = kkt.dual_signs_valid;
                    attempt.numerical_certified =
                        std::isfinite(attempt.lower_value)
                        && attempt.projected_kkt_inf_norm <= kStageIIKktTolerance
                        && attempt.constraint_violation <= kCertificateTolerance
                        && attempt.complementarity <= kCertificateTolerance
                        && kkt.reconstruction_inf_norm <= kCoordinateTolerance
                        && kkt.dual_signs_valid
                        && packing_domain
                        && physical;
                } catch (const std::exception& error) {
                    if (attempt.callback_error.empty()) {
                        attempt.callback_error = error.what();
                    }
                }
            }
            attempt.improving = attempt.numerical_certified
                && attempt.lower_value < upper_bound - stage_ii_tolerance;
            result.attempt_log.push_back(attempt);
            if (!attempt.numerical_certified) {
                ++result.lower_failed_start_count;
                if (result.first_failed_lower_start_index < 0) {
                    result.first_failed_lower_start_index =
                        static_cast<int>(start_index);
                    result.first_failed_lower_solver_status = run.solver_status;
                    result.first_failed_lower_reason = attempt.callback_error.empty()
                        ? "Stage II lower terminal failed its numerical certificate"
                        : attempt.callback_error;
                    result.first_failed_lower_initial = start;
                }
                continue;
            }
            ++certified_start_count;
            ++result.lower_completed_start_count;
            if (attempt.lower_value < best_certified_value) {
                best_certified_value = attempt.lower_value;
                best = {
                    std::move(state),
                    std::move(independent),
                    {},
                    run.variables.back(),
                };
            }
        }
        ++result.major_iterations;
        if (!std::isfinite(best_certified_value)) {
            result.outcome = "indeterminate";
            result.cut_count = static_cast<int>(cuts.size());
            return result;
        }
        const Held2StageIILowerDecision decision = decide_held2_stage_ii_lower(
            upper_bound,
            best_certified_value,
            certified_start_count,
            result.lower_starts_per_iteration
        );
        if (result.search_completeness == "not_run") {
            result.search_completeness = decision.search_completeness;
        } else if (decision.search_completeness == "partial") {
            result.search_completeness = "partial";
        }
        result.final_lower_search_complete = decision.lower_bound_certified;
        if (decision.decision == "indeterminate") {
            result.outcome = "indeterminate";
            result.cut_count = static_cast<int>(cuts.size());
            return result;
        }
        best = make_cut(
            std::move(best.state),
            std::move(best.independent_modified_fractions),
            best.phase_coordinate
        );
        const bool duplicate = std::any_of(
            cuts.begin(),
            cuts.end(),
            [&best](const Held2StageIIStep6Cut& cut) {
                return maximum_abs_difference(
                           cut.state.modified_fractions,
                           best.state.modified_fractions
                       ) < 1.0e-7
                    && std::abs(cut.state.volume - best.state.volume) < 1.0e-7;
            }
        );
        const auto select_current_candidate_set = [&]() {
            result.candidates = select_held2_stage_ii_candidates(
                coordinates, feed, upper_bound, multiplier, cuts
            );
            if (result.candidates.size() <= 1) {
                return false;
            }
            result.outcome = "candidate_set";
            result.cut_count = static_cast<int>(cuts.size());
            return true;
        };
        if (decision.decision == "add_improving_cut") {
            if (duplicate) {
                result.outcome = "indeterminate";
                result.cut_count = static_cast<int>(cuts.size());
                return result;
            }
            cuts.push_back(std::move(best));
            ++result.certified_improving_cut_count;
            if (decision.lower_bound_certified) {
                result.bound_history.push_back({
                    best_certified_value,
                    upper_bound,
                    multiplier.empty() ? 0.0 : multiplier.front(),
                    static_cast<int>(cuts.size() - 1),
                });
            }
            if (select_current_candidate_set()) {
                return result;
            }
            continue;
        }
        result.bound_history.push_back({
            best_certified_value,
            upper_bound,
            multiplier.empty() ? 0.0 : multiplier.front(),
            static_cast<int>(cuts.size()),
        });
        if (!duplicate) {
            cuts.push_back(std::move(best));
        }
        if (select_current_candidate_set()) {
            return result;
        }
        if (duplicate) {
            result.outcome = "no_progress";
            result.cut_count = static_cast<int>(cuts.size());
            return result;
        }
    }
    result.outcome = "resource_limit";
    result.cut_count = static_cast<int>(cuts.size());
    return result;
}

Held2StageIILowerDecision decide_held2_stage_ii_lower(
    double upper_bound,
    double best_certified_value,
    int certified_start_count,
    int declared_start_count
) {
    constexpr double stage_ii_tolerance = 1.0e-8;
    if (!std::isfinite(upper_bound) || !std::isfinite(best_certified_value)
        || declared_start_count <= 0 || certified_start_count < 0
        || certified_start_count > declared_start_count) {
        throw std::invalid_argument("HELD2 Stage II lower evidence is invalid");
    }
    Held2StageIILowerDecision result;
    result.search_completeness = certified_start_count == declared_start_count
        ? "complete"
        : "partial";
    result.lower_bound_certified = result.search_completeness == "complete";
    if (best_certified_value < upper_bound - stage_ii_tolerance) {
        result.decision = "add_improving_cut";
    } else if (!result.lower_bound_certified) {
        result.decision = "indeterminate";
    } else {
        result.decision = "evaluate_final_gap";
    }
    return result;
}

Held2ManufacturedEvaluation evaluate_held2_manufactured(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::array<double, 4>& variables,
    const std::vector<double>& chemical_potentials
) {
    if (!std::all_of(variables.begin(), variables.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::invalid_argument("HELD2 reduced vector must contain four finite entries");
    }
    Held2ManufacturedEvaluation result;
    result.coordinates = make_held2_coordinates(charges);
    if (result.coordinates.retained_indices.size() != 2
        || result.coordinates.independent_indices.size() != 1) {
        throw std::invalid_argument(
            "manufactured HELD2 seam requires one independent modified composition"
        );
    }
    result.modified_feed = held2_transform_physical_fractions(
        result.coordinates,
        physical_feed
    );
    result.transformed_modified_potentials = held2_transform_modified_potentials(
        result.coordinates,
        chemical_potentials
    );
    const auto dependent_position = static_cast<std::size_t>(
        std::find(
            result.coordinates.retained_indices.begin(),
            result.coordinates.retained_indices.end(),
            result.coordinates.dependent_index
        ) - result.coordinates.retained_indices.begin()
    );
    const std::size_t independent_position = 1 - dependent_position;
    const double feed = result.modified_feed[independent_position];
    const double u_alpha = variables[0];
    const double u_beta = variables[1];
    const double volume_alpha = variables[2];
    const double volume_beta = variables[3];
    if (!(u_alpha < feed && feed < u_beta)) {
        throw std::invalid_argument("modified compositions must straddle the modified feed");
    }
    const double lower_bound = result.coordinates.independent_lower_bounds[0];
    const double upper_bound = result.coordinates.independent_upper_bounds[0];
    if (u_alpha < lower_bound || u_beta > upper_bound) {
        throw std::invalid_argument(
            "HELD2 modified compositions must satisfy the declared source bounds"
        );
    }
    if (volume_alpha <= 0.0 || volume_beta <= 0.0) {
        throw std::invalid_argument("HELD2 manufactured state is outside its physical domain");
    }
    result.modified_phases[0] = std::vector<double>(2, 0.0);
    result.modified_phases[1] = std::vector<double>(2, 0.0);
    result.modified_phases[0][dependent_position] = 1.0 - u_alpha;
    result.modified_phases[0][independent_position] = u_alpha;
    result.modified_phases[1][dependent_position] = 1.0 - u_beta;
    result.modified_phases[1][independent_position] = u_beta;
    result.physical_phases[0] = held2_lift_modified_fractions(
        result.coordinates,
        result.modified_phases[0]
    );
    result.physical_phases[1] = held2_lift_modified_fractions(
        result.coordinates,
        result.modified_phases[1]
    );

    const double separation = u_beta - u_alpha;
    result.phase_fraction = (u_beta - feed) / separation;
    const double complement = 1.0 - result.phase_fraction;
    const double gibbs_alpha = manufactured_gibbs(u_alpha, volume_alpha);
    const double gibbs_beta = manufactured_gibbs(u_beta, volume_beta);
    result.phase_gibbs_gradients[0] = manufactured_gibbs_gradient(
        u_alpha,
        volume_alpha
    );
    result.phase_gibbs_gradients[1] = manufactured_gibbs_gradient(
        u_beta,
        volume_beta
    );
    result.objective =
        result.phase_fraction * gibbs_alpha + complement * gibbs_beta;
    result.gradient = {
        result.phase_fraction / separation * (gibbs_alpha - gibbs_beta)
            + result.phase_fraction * result.phase_gibbs_gradients[0][0],
        complement / separation * (gibbs_alpha - gibbs_beta)
            + complement * result.phase_gibbs_gradients[1][0],
        result.phase_fraction * result.phase_gibbs_gradients[0][1],
        complement * result.phase_gibbs_gradients[1][1],
    };
    result.phase_modified_potentials[0] = manufactured_modified_potentials(
        u_alpha,
        volume_alpha
    );
    result.phase_modified_potentials[1] = manufactured_modified_potentials(
        u_beta,
        volume_beta
    );
    result.modified_potential_gap = std::max(
        std::abs(
            result.phase_modified_potentials[0][0]
            - result.phase_modified_potentials[1][0]
        ),
        std::abs(
            result.phase_modified_potentials[0][1]
            - result.phase_modified_potentials[1][1]
        )
    );
    result.pressure_stationarity_inf_norm = std::max(
        std::abs(result.phase_gibbs_gradients[0][1]),
        std::abs(result.phase_gibbs_gradients[1][1])
    );

    result.modified_balance.resize(2, 0.0);
    result.ordinary_balance.resize(charges.size(), 0.0);
    for (std::size_t index = 0; index < result.modified_balance.size(); ++index) {
        result.modified_balance[index] =
            result.phase_fraction * result.modified_phases[0][index]
            + complement * result.modified_phases[1][index];
    }
    for (std::size_t index = 0; index < result.ordinary_balance.size(); ++index) {
        result.ordinary_balance[index] =
            result.phase_fraction * result.physical_phases[0][index]
            + complement * result.physical_phases[1][index];
    }
    result.phase_charge_residuals = {
        charge_residual(charges, result.physical_phases[0]),
        charge_residual(charges, result.physical_phases[1]),
    };

    result.certificate.modified_balance_abs = maximum_abs_difference(
        result.modified_balance,
        result.modified_feed
    );
    result.certificate.ordinary_balance_inf_norm = maximum_abs_difference(
        result.ordinary_balance,
        physical_feed
    );
    result.certificate.phase_charge_inf_norm = std::max(
        std::abs(result.phase_charge_residuals[0]),
        std::abs(result.phase_charge_residuals[1])
    );
    result.certificate.modified_potential_gap = result.modified_potential_gap;
    result.certificate.pressure_stationarity_inf_norm =
        result.pressure_stationarity_inf_norm;
    for (double value : result.gradient) {
        result.certificate.reduced_kkt_inf_norm = std::max(
            result.certificate.reduced_kkt_inf_norm,
            std::abs(value)
        );
    }
    result.certificate.enumeration_objective_gap =
        result.objective - enumerated_manufactured_objective(feed);
    result.certificate.independent_modified_composition_count = 1.0;
    const double maximum_metric = std::max({
        result.certificate.modified_balance_abs,
        result.certificate.ordinary_balance_inf_norm,
        result.certificate.phase_charge_inf_norm,
        result.certificate.modified_potential_gap,
        result.certificate.pressure_stationarity_inf_norm,
        result.certificate.reduced_kkt_inf_norm,
        std::abs(result.certificate.enumeration_objective_gap),
    });
    result.certificate.accepted = maximum_metric < kCertificateTolerance;
    result.certificate.independent_evidence =
        std::abs(result.certificate.enumeration_objective_gap) < kCertificateTolerance;
    return result;
}

}  // namespace epcsaft_equilibrium
