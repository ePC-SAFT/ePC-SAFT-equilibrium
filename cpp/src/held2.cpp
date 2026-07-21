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
};

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
        double composition_sum_upper
    )
        : evaluator_(std::move(evaluator)),
          reference_variables_(std::move(reference_variables)),
          reference_(std::move(reference)),
          use_tpd_(use_tpd),
          initial_(std::move(initial)),
          lower_(std::move(lower)),
          upper_(std::move(upper)),
          composition_sum_upper_(composition_sum_upper) {}

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
        const Ipopt::Index expected_nonzeros = has_composition_constraint() ? n - 1 : 0;
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
        };
    }

private:
    [[nodiscard]] Ipopt::Index constraint_count() const {
        return static_cast<Ipopt::Index>(has_composition_constraint());
    }

    [[nodiscard]] bool has_composition_constraint() const {
        return std::isfinite(composition_sum_upper_);
    }

    [[nodiscard]] Ipopt::Index composition_constraint_row() const {
        return 0;
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
    double composition_sum_upper = std::numeric_limits<double>::infinity()
) {
    Ipopt::SmartPtr<Held2SearchTnlp> problem = new Held2SearchTnlp(
        evaluator,
        reference_variables,
        reference,
        use_tpd,
        initial,
        lower,
        upper,
        composition_sum_upper
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
    double left = std::log(molar_volume_bounds[0]);
    double right = std::log(molar_volume_bounds[1]);
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
    int solve_start_limit
) {
    if (!std::isfinite(molar_volume_bounds[0]) || !std::isfinite(molar_volume_bounds[1])
        || molar_volume_bounds[0] <= 0.0
        || molar_volume_bounds[1] <= molar_volume_bounds[0]) {
        throw std::invalid_argument(
            "HELD2 Stage I molar-volume bounds must be finite, positive, and ordered"
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

    Held2StageIResult result;
    result.volume_domain_search_complete = volume_domain_search_complete;
    result.declared_start_count = 10 * static_cast<int>(coordinates.charges.size());
    result.reference_scan_interval_count = result.declared_start_count;
    result.minimum_tpd = std::numeric_limits<double>::infinity();
    const double reference_log_lower = std::log(molar_volume_bounds[0]);
    const double reference_log_upper = std::log(molar_volume_bounds[1]);
    const double reference_log_span = reference_log_upper - reference_log_lower;
    std::vector<std::pair<double, Held2StateEvaluation>> reference_scan;
    reference_scan.reserve(
        static_cast<std::size_t>(result.reference_scan_interval_count + 1)
    );
    for (int point = 0; point <= result.reference_scan_interval_count; ++point) {
        const double fraction = static_cast<double>(point)
            / static_cast<double>(result.reference_scan_interval_count);
        const double log_volume = reference_log_lower + fraction * reference_log_span;
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
        return result;
    }

    Held2StateEvaluation reference;
    std::vector<double> reference_variables;
    double reference_objective = std::numeric_limits<double>::infinity();
    for (std::size_t interval = 0; interval + 1 < reference_scan.size(); ++interval) {
        double left_log_volume = reference_scan[interval].first;
        Held2StateEvaluation left_state = reference_scan[interval].second;
        double right_log_volume = reference_scan[interval + 1].first;
        Held2StateEvaluation right_state = reference_scan[interval + 1].second;
        double left_residual = left_state.pressure_stationarity_relative;
        double right_residual = right_state.pressure_stationarity_relative;
        if (left_residual * right_residual > 0.0) {
            continue;
        }

        double root_log_volume = 0.0;
        Held2StateEvaluation root_state;
        bool refined = false;
        if (std::abs(left_residual) <= kCertificateTolerance) {
            root_log_volume = left_log_volume;
            root_state = std::move(left_state);
            refined = true;
        } else if (std::abs(right_residual) <= kCertificateTolerance) {
            root_log_volume = right_log_volume;
            root_state = std::move(right_state);
            refined = true;
        } else {
            for (int iteration = 0; iteration < 100; ++iteration) {
                root_log_volume = 0.5 * (left_log_volume + right_log_volume);
                try {
                    root_state = log_volume_evaluator(feed_independent, root_log_volume);
                } catch (const std::exception&) {
                    ++result.reference_evaluation_failure_count;
                    break;
                }
                const double root_residual =
                    root_state.pressure_stationarity_relative;
                if (std::abs(root_residual) <= kCertificateTolerance) {
                    refined = true;
                    break;
                }
                if (left_residual * root_residual <= 0.0) {
                    right_log_volume = root_log_volume;
                    right_residual = root_residual;
                } else {
                    left_log_volume = root_log_volume;
                    left_residual = root_residual;
                }
            }
        }
        if (!refined) {
            ++result.reference_refinement_failure_count;
            continue;
        }
        const bool duplicate = std::any_of(
            result.reference_roots.begin(),
            result.reference_roots.end(),
            [root_log_volume](const Held2ReferenceRoot& root) {
                return std::abs(root.log_volume - root_log_volume) < 1.0e-8;
            }
        );
        if (duplicate) {
            continue;
        }
        const double curvature = root_state.hessian.back();
        if (!std::isfinite(curvature) || curvature == 0.0) {
            ++result.reference_refinement_failure_count;
            continue;
        }
        const bool mechanically_stable = curvature > 0.0;
        result.reference_roots.push_back({
            root_log_volume,
            root_state.volume,
            root_state.objective,
            root_state.pressure_stationarity_relative,
            curvature,
            mechanically_stable,
        });
        if (!mechanically_stable) {
            continue;
        }
        ++result.reference_stable_root_count;
        if (root_state.objective < reference_objective) {
            reference_objective = root_state.objective;
            reference = std::move(root_state);
            reference_variables = feed_independent;
            reference_variables.push_back(root_log_volume);
        }
    }
    result.reference_root_count = static_cast<int>(result.reference_roots.size());
    if (result.reference_evaluation_failure_count != 0
        || result.reference_refinement_failure_count != 0
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
    lower.push_back(
        search_uses_log_packing
            ? std::log(kHeld2PackingFractionMinimum)
            : std::log(molar_volume_bounds[0])
    );
    upper.push_back(
        search_uses_log_packing
            ? std::log(kHeld2PackingFractionMaximum)
            : std::log(molar_volume_bounds[1])
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
        starts.push_back({0.2, 0.0});
        starts.push_back({0.8, 0.0});
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
        std::uniform_real_distribution<double> log_volume_distribution(
            std::log(start_volume_bounds[0]),
            std::log(start_volume_bounds[1])
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
            composition_sum_upper
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
    const std::vector<double>& physical_feed
) {
    const Held2Coordinates coordinates = make_held2_coordinates(charges);
    if (coordinates.independent_indices.size() != 1) {
        throw std::invalid_argument(
            "manufactured HELD2 Stage I requires one independent modified composition"
        );
    }
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
    const Held2VolumeBoundsEvaluator volume_bounds_evaluator = [](
        const std::vector<double>&
    ) {
        return std::array<double, 2>{0.5, 1.5};
    };
    return solve_held2_stage_i(
        coordinates,
        physical_feed,
        evaluator,
        evaluator,
        {0.5, 1.5},
        volume_bounds_evaluator,
        false,
        true,
        -1
    );
}

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

Held2StageIIResult solve_held2_stage_ii(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const Held2StateEvaluator& evaluator,
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
    struct Cut {
        Held2StateEvaluation state;
        std::vector<double> independent;
        double phase_coordinate = 0.0;
    };
    std::vector<Cut> cuts;
    cuts.push_back({reference, feed, std::log(reference.packing.value)});
    for (const Held2StageICandidate& candidate : stage_i_candidates) {
        std::vector<double> independent;
        for (std::size_t position : retained_positions) {
            independent.push_back(candidate.modified_fractions[position]);
        }
        if (candidate.packing_fraction <= 0.0) {
            continue;
        }
        const double coordinate = std::log(candidate.packing_fraction);
        cuts.push_back({evaluator(independent, coordinate), independent, coordinate});
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
        cuts.push_back({evaluator(point, cuts.front().phase_coordinate), point,
                        cuts.front().phase_coordinate});
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
                if (row == pivot) continue;
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
                    const Cut& cut = cuts[selected[row]];
                    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                        matrix[row][coordinate] = -(
                            feed[coordinate] - cut.independent[coordinate]
                        );
                    }
                    matrix[row][dimension] = 1.0;
                    rhs[row] = cut.state.objective;
                }
                const std::vector<double> solution = solve_dense(matrix, rhs);
                if (solution.empty()) return;
                const double envelope = solution.back();
                for (const Cut& cut : cuts) {
                    double line = cut.state.objective;
                    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                        line += solution[coordinate]
                            * (feed[coordinate] - cut.independent[coordinate]);
                    }
                    if (envelope > line + stage_ii_tolerance) return;
                }
                if (envelope > upper_bound) {
                    upper_bound = envelope;
                    std::copy(solution.begin(), solution.end() - 1, multiplier.begin());
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

        Held2StateEvaluation lower_reference;
        lower_reference.gradient = multiplier;
        lower_reference.gradient.push_back(0.0);
        lower_reference.hessian.assign((dimension + 1) * (dimension + 1), 0.0);
        std::vector<double> lower_reference_variables = feed;
        lower_reference_variables.push_back(0.0);
        std::vector<double> lower = coordinates.independent_lower_bounds;
        std::vector<double> upper = coordinates.independent_upper_bounds;
        for (std::size_t index = 0; index < upper.size(); ++index) {
            upper[index] = std::nextafter(upper[index], lower[index]);
        }
        lower.push_back(std::log(kHeld2PackingFractionMinimum));
        upper.push_back(std::log(kHeld2PackingFractionMaximum));
        std::vector<std::vector<double>> starts;
        for (const Cut& cut : cuts) {
            std::vector<double> start = cut.independent;
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
        Cut best;
        int certified_start_count = 0;
        for (std::size_t start_index = 0; start_index < starts.size(); ++start_index) {
            const std::vector<double>& start = starts[start_index];
            ++result.lower_attempted_start_count;
            const Held2SearchRun run = solve_held2_search(
                evaluator, lower_reference_variables, lower_reference, true,
                start, lower, upper, composition_sum_upper
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
            Held2StateEvaluation state;
            std::vector<double> independent;
            if (run.variables.size() == dimension + 1) {
                try {
                    independent.assign(run.variables.begin(), run.variables.end() - 1);
                    state = evaluator(independent, run.variables.back());
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
                    if (run.lower_bound_multipliers.size() == dimension + 1
                        && run.upper_bound_multipliers.size() == dimension + 1) {
                        std::vector<double> kkt = state.gradient;
                        for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                            kkt[coordinate] -= multiplier[coordinate];
                            if (!run.constraint_multipliers.empty()) {
                                kkt[coordinate] += run.constraint_multipliers.front();
                            }
                        }
                        attempt.projected_kkt_inf_norm = 0.0;
                        attempt.complementarity = 0.0;
                        for (std::size_t index = 0; index < kkt.size(); ++index) {
                            kkt[index] -= run.lower_bound_multipliers[index];
                            kkt[index] += run.upper_bound_multipliers[index];
                            attempt.projected_kkt_inf_norm = std::max(
                                attempt.projected_kkt_inf_norm,
                                std::abs(kkt[index])
                            );
                            attempt.complementarity = std::max({
                                attempt.complementarity,
                                std::abs(run.lower_bound_multipliers[index]
                                    * (run.variables[index] - lower[index])),
                                std::abs(run.upper_bound_multipliers[index]
                                    * (upper[index] - run.variables[index])),
                            });
                        }
                        if (!run.constraint_multipliers.empty()) {
                            attempt.complementarity = std::max(
                                attempt.complementarity,
                                std::abs(run.constraint_multipliers.front()
                                    * (composition_sum_upper
                                       - std::accumulate(
                                           independent.begin(), independent.end(), 0.0
                                       )))
                            );
                        }
                    }
                    attempt.numerical_certified =
                        std::isfinite(attempt.lower_value)
                        && attempt.projected_kkt_inf_norm <= kStageIIKktTolerance
                        && attempt.constraint_violation <= kCertificateTolerance
                        && attempt.complementarity <= kCertificateTolerance
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
                best = {std::move(state), std::move(independent), run.variables.back()};
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
        const bool duplicate = std::any_of(cuts.begin(), cuts.end(), [&best](const Cut& cut) {
            return maximum_abs_difference(cut.state.modified_fractions, best.state.modified_fractions) < 1.0e-7
                && std::abs(cut.state.volume - best.state.volume) < 1.0e-7;
        });
        if (decision.decision == "add_improving_cut") {
            if (duplicate) {
                result.outcome = "indeterminate";
                result.cut_count = static_cast<int>(cuts.size());
                return result;
            }
            cuts.push_back(best);
            ++result.certified_improving_cut_count;
            if (decision.lower_bound_certified) {
                result.bound_history.push_back({
                    best_certified_value,
                    upper_bound,
                    multiplier.empty() ? 0.0 : multiplier.front(),
                    static_cast<int>(cuts.size() - 1),
                });
            }
            continue;
        }
        if (decision.decision == "indeterminate") {
            result.outcome = "indeterminate";
            result.cut_count = static_cast<int>(cuts.size());
            return result;
        }
        result.bound_history.push_back({
            best_certified_value,
            upper_bound,
            multiplier.empty() ? 0.0 : multiplier.front(),
            static_cast<int>(cuts.size()),
        });
        if (!duplicate) cuts.push_back(best);
        result.candidates.clear();
        for (const Cut& cut : cuts) {
            double value = cut.state.objective;
            for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                value += multiplier[coordinate]
                    * (feed[coordinate] - cut.independent[coordinate]);
            }
            const double gap = upper_bound - value;
            bool stationary = std::abs(cut.state.gradient.back()) <= stage_ii_tolerance;
            for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                stationary = stationary
                    && std::abs(cut.state.gradient[coordinate] - multiplier[coordinate])
                        <= stage_ii_tolerance;
            }
            if (std::abs(gap) <= stage_ii_tolerance && stationary) {
                result.candidates.push_back({
                    cut.state.modified_fractions,
                    cut.independent,
                    cut.state.volume,
                    cut.phase_coordinate,
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
