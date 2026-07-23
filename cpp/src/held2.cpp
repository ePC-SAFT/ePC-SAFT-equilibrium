#include "held2.hpp"
#include "held2_progress.hpp"
#include "held2_stage_ii_basin.hpp"
#include "held2_stage_ii_upper.hpp"
#include "held2_tolerances.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpDenseVector.hpp>
#include <coin/IpIpoptData.hpp>
#include <coin/IpTNLP.hpp>

namespace epcsaft_equilibrium {
namespace {

// Perdomo Eq. (61) uses 1e-10 times the eliminated-ion coordinate factor.
constexpr double kModifiedLowerScale = 1.0e-10;
constexpr int kStageIIpoptIterations = 300;
constexpr int kStageIISobolEvaluationCount = 12;
constexpr int kStageIIDirectEvaluationBudget = 128;
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

double charge_scale(
    const std::vector<double>& charges,
    const std::vector<double>& fractions
) {
    double result = 0.0;
    for (std::size_t index = 0; index < charges.size(); ++index) {
        result += std::abs(charges[index] * fractions[index]);
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
    application->Options()->SetNumericValue("tol", kHeld2IpoptTarget.atol);
    application->Options()->SetNumericValue(
        "acceptable_tol", kHeld2IpoptAcceptable.atol
    );
    application->Options()->SetIntegerValue("acceptable_iter", 0);
    application->Options()->SetStringValue("jacobian_approximation", "exact");
    application->Options()->SetStringValue("hessian_approximation", "exact");
    application->Options()->SetStringValue("nlp_scaling_method", "none");
    application->Options()->SetNumericValue(
        "constr_viol_tol", kHeld2IpoptConstraint.atol
    );
    application->Options()->SetNumericValue("bound_relax_factor", 0.0);
    application->Options()->SetStringValue("honor_original_bounds", "yes");
    application->Options()->SetStringValue("check_derivatives_for_naninf", "yes");
}

struct Held2SearchRun {
    bool solver_converged = false;
    std::string solver_status = "not_run";
    std::string callback_error;
    std::vector<double> variables;
    std::vector<double> lower_bound_multipliers;
    std::vector<double> upper_bound_multipliers;
    std::vector<double> coordinate_jacobian;
    std::vector<Held2StageIILocalEvaluation> evaluation_trace;
};

struct Held2StageIISimplexChart {
    std::vector<double> physical_coordinates;
    std::vector<double> chart_coordinates;
    std::vector<double> jacobian;
    std::vector<double> component_hessians;
    bool singular = false;
    bool normalized_boundary_contact = false;
};

struct Held2SecondOrderScalar {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
};

Held2SecondOrderScalar second_order_constant(double value, std::size_t dimension) {
    return {value, std::vector<double>(dimension, 0.0),
            std::vector<double>(dimension * dimension, 0.0)};
}

Held2SecondOrderScalar second_order_variable(
    double value,
    std::size_t dimension,
    std::size_t index
) {
    Held2SecondOrderScalar result = second_order_constant(value, dimension);
    result.gradient[index] = 1.0;
    return result;
}

Held2SecondOrderScalar second_order_add(
    const Held2SecondOrderScalar& left,
    const Held2SecondOrderScalar& right,
    double right_factor = 1.0
) {
    Held2SecondOrderScalar result = left;
    result.value += right_factor * right.value;
    for (std::size_t index = 0; index < result.gradient.size(); ++index) {
        result.gradient[index] += right_factor * right.gradient[index];
    }
    for (std::size_t index = 0; index < result.hessian.size(); ++index) {
        result.hessian[index] += right_factor * right.hessian[index];
    }
    return result;
}

Held2SecondOrderScalar second_order_multiply(
    const Held2SecondOrderScalar& left,
    const Held2SecondOrderScalar& right
) {
    const std::size_t dimension = left.gradient.size();
    Held2SecondOrderScalar result = second_order_constant(
        left.value * right.value, dimension
    );
    for (std::size_t row = 0; row < dimension; ++row) {
        result.gradient[row] = left.gradient[row] * right.value
            + left.value * right.gradient[row];
        for (std::size_t column = 0; column < dimension; ++column) {
            result.hessian[row * dimension + column] =
                left.hessian[row * dimension + column] * right.value
                + left.gradient[row] * right.gradient[column]
                + left.gradient[column] * right.gradient[row]
                + left.value * right.hessian[row * dimension + column];
        }
    }
    return result;
}

struct Held2PhysicalKkt {
    double primal_inf_norm = std::numeric_limits<double>::infinity();
    double dual_sign_violation_inf_norm = std::numeric_limits<double>::infinity();
    double stationarity_inf_norm = std::numeric_limits<double>::infinity();
    double complementarity = std::numeric_limits<double>::infinity();
    double reconstruction_inf_norm = std::numeric_limits<double>::infinity();
    double reconstruction_scale = std::numeric_limits<double>::infinity();
    bool dual_signs_valid = false;
};

std::optional<std::vector<double>> solve_dense_linear_system(
    const std::vector<double>& matrix,
    const std::vector<double>& right_hand_side,
    std::size_t dimension
) {
    if (dimension == 0 || matrix.size() != dimension * dimension
        || right_hand_side.size() != dimension) {
        return std::nullopt;
    }
    std::vector<double> factors = matrix;
    std::vector<double> solution = right_hand_side;
    double maximum_entry = 0.0;
    for (double value : factors) {
        maximum_entry = std::max(maximum_entry, std::abs(value));
    }
    const double singular_threshold = std::numeric_limits<double>::epsilon()
        * std::max(1.0, maximum_entry) * static_cast<double>(dimension);
    for (std::size_t pivot = 0; pivot < dimension; ++pivot) {
        std::size_t pivot_row = pivot;
        for (std::size_t row = pivot + 1; row < dimension; ++row) {
            if (std::abs(factors[row * dimension + pivot])
                > std::abs(factors[pivot_row * dimension + pivot])) {
                pivot_row = row;
            }
        }
        if (std::abs(factors[pivot_row * dimension + pivot])
            <= singular_threshold) {
            return std::nullopt;
        }
        if (pivot_row != pivot) {
            for (std::size_t column = 0; column < dimension; ++column) {
                std::swap(
                    factors[pivot * dimension + column],
                    factors[pivot_row * dimension + column]
                );
            }
            std::swap(solution[pivot], solution[pivot_row]);
        }
        for (std::size_t row = pivot + 1; row < dimension; ++row) {
            const double multiplier = factors[row * dimension + pivot]
                / factors[pivot * dimension + pivot];
            factors[row * dimension + pivot] = 0.0;
            for (std::size_t column = pivot + 1; column < dimension; ++column) {
                factors[row * dimension + column] -= multiplier
                    * factors[pivot * dimension + column];
            }
            solution[row] -= multiplier * solution[pivot];
        }
    }
    for (std::size_t row = dimension; row-- > 0;) {
        for (std::size_t column = row + 1; column < dimension; ++column) {
            solution[row] -= factors[row * dimension + column] * solution[column];
        }
        solution[row] /= factors[row * dimension + row];
    }
    return solution;
}

double matrix_condition_inf(
    const std::vector<double>& matrix,
    std::size_t dimension
) {
    if (dimension == 0 || matrix.size() != dimension * dimension) {
        return std::numeric_limits<double>::infinity();
    }
    double matrix_norm = 0.0;
    for (std::size_t row = 0; row < dimension; ++row) {
        double row_sum = 0.0;
        for (std::size_t column = 0; column < dimension; ++column) {
            row_sum += std::abs(matrix[row * dimension + column]);
        }
        matrix_norm = std::max(matrix_norm, row_sum);
    }
    std::vector<double> inverse(dimension * dimension, 0.0);
    for (std::size_t column = 0; column < dimension; ++column) {
        std::vector<double> unit(dimension, 0.0);
        unit[column] = 1.0;
        const std::optional<std::vector<double>> inverse_column =
            solve_dense_linear_system(matrix, unit, dimension);
        if (!inverse_column.has_value()) {
            return std::numeric_limits<double>::infinity();
        }
        for (std::size_t row = 0; row < dimension; ++row) {
            inverse[row * dimension + column] = (*inverse_column)[row];
        }
    }
    double inverse_norm = 0.0;
    for (std::size_t row = 0; row < dimension; ++row) {
        double row_sum = 0.0;
        for (std::size_t column = 0; column < dimension; ++column) {
            row_sum += std::abs(inverse[row * dimension + column]);
        }
        inverse_norm = std::max(inverse_norm, row_sum);
    }
    return matrix_norm * inverse_norm;
}

class Held2SearchObjectiveEvaluator {
public:
    Held2SearchObjectiveEvaluator(
        Held2StateEvaluator phase_evaluator,
        std::vector<double> reference_variables,
        Held2StateEvaluation reference,
        bool use_tpd
    )
        : phase_evaluator_(std::move(phase_evaluator)),
          reference_variables_(std::move(reference_variables)),
          reference_(std::move(reference)),
          use_tpd_(use_tpd) {}

    [[nodiscard]] Held2StateEvaluation evaluate(
        const std::vector<double>& variables
    ) const {
        if (variables.empty()
            || variables.size() != reference_variables_.size()) {
            throw std::invalid_argument("HELD2 search coordinate count changed");
        }
        std::vector<double> independent(variables.begin(), variables.end() - 1);
        Held2StateEvaluation evaluation = phase_evaluator_(
            independent,
            variables.back()
        );
        const std::size_t dimension = variables.size();
        if (evaluation.gradient.size() != dimension
            || evaluation.hessian.size() != dimension * dimension) {
            throw std::invalid_argument("HELD2 search derivative dimensions changed");
        }
        if (!use_tpd_) {
            return evaluation;
        }
        if (reference_.gradient.size() != dimension) {
            throw std::invalid_argument("HELD2 reference gradient dimension changed");
        }
        evaluation.objective -= reference_.objective;
        for (std::size_t index = 0; index < dimension; ++index) {
            evaluation.objective -= reference_.gradient[index]
                * (variables[index] - reference_variables_[index]);
            evaluation.gradient[index] -= reference_.gradient[index];
        }
        return evaluation;
    }

private:
    Held2StateEvaluator phase_evaluator_;
    std::vector<double> reference_variables_;
    Held2StateEvaluation reference_;
    bool use_tpd_ = false;
};

class Held2SearchTnlp final : public Ipopt::TNLP {
public:
    Held2SearchTnlp(
        std::function<Held2StateEvaluation(const std::vector<double>&)> objective,
        std::vector<double> initial,
        std::vector<double> lower,
        std::vector<double> upper,
        Held2ProgressObserver* observer,
        std::string progress_stage,
        int progress_major,
        int progress_attempt
    )
        : objective_(std::move(objective)),
          initial_(std::move(initial)),
          lower_(std::move(lower)),
          upper_(std::move(upper)),
          observer_(observer),
          progress_stage_(std::move(progress_stage)),
          progress_major_(progress_major),
          progress_attempt_(progress_attempt) {}

    bool get_nlp_info(
        Ipopt::Index& n,
        Ipopt::Index& m,
        Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag,
        IndexStyleEnum& index_style
    ) override {
        n = static_cast<Ipopt::Index>(initial_.size());
        m = 0;
        nnz_jac_g = 0;
        nnz_h_lag = n * (n + 1) / 2;
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
        if (n != static_cast<Ipopt::Index>(initial_.size()) || m != 0) {
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
        if (n != static_cast<Ipopt::Index>(initial_.size()) || m != 0 || !init_x
            || init_z || init_lambda) {
            return false;
        }
        std::copy(initial_.begin(), initial_.end(), x);
        return true;
    }

    bool eval_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Number& objective
    ) override {
        try {
            objective = evaluate(n, x).objective;
            record_evaluation("eval_f", n, x, new_x, true, "");
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            record_evaluation("eval_f", n, x, new_x, false, error.what());
            return false;
        }
    }

    bool eval_grad_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Number* gradient
    ) override {
        try {
            const Held2StateEvaluation evaluation = evaluate(n, x);
            std::copy(evaluation.gradient.begin(), evaluation.gradient.end(), gradient);
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            record_evaluation("eval_grad_f", n, x, new_x, false, error.what());
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
        bool new_x,
        Ipopt::Number objective_factor,
        Ipopt::Index m,
        const Ipopt::Number*,
        bool,
        Ipopt::Index nonzero_count,
        Ipopt::Index* rows,
        Ipopt::Index* columns,
        Ipopt::Number* values
    ) override {
        if (m != 0 || nonzero_count != n * (n + 1) / 2) {
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
            record_evaluation("eval_h", n, x, new_x, false, error.what());
            return false;
        }
    }

    bool intermediate_callback(
        Ipopt::AlgorithmMode,
        Ipopt::Index iteration,
        Ipopt::Number objective,
        Ipopt::Number primal_residual,
        Ipopt::Number dual_residual,
        Ipopt::Number barrier,
        Ipopt::Number step_norm,
        Ipopt::Number,
        Ipopt::Number dual_step,
        Ipopt::Number primal_step,
        Ipopt::Index line_search_steps,
        const Ipopt::IpoptData* ip_data,
        Ipopt::IpoptCalculatedQuantities*
    ) override {
        Held2ProgressEvent progress;
        progress.kind = Held2ProgressKind::LocalIteration;
        progress.stage = progress_stage_;
        progress.iteration = static_cast<int>(iteration);
        progress.major_iteration = progress_major_;
        progress.attempt = progress_attempt_;
        progress.objective = objective;
        progress.primal_residual = primal_residual;
        progress.dual_residual = dual_residual;
        progress.complementarity = barrier;
        progress.step_norm = step_norm;
        progress.dual_step = dual_step;
        progress.primal_step = primal_step;
        progress.line_search_steps = static_cast<int>(line_search_steps);
        observe_held2(observer_, progress);
        if (ip_data != nullptr && Ipopt::IsValid(ip_data->curr())
            && Ipopt::IsValid(ip_data->curr()->x())) {
            const auto* dense = dynamic_cast<const Ipopt::DenseVector*>(
                Ipopt::GetRawPtr(ip_data->curr()->x())
            );
            if (dense != nullptr
                && dense->Dim() == static_cast<Ipopt::Index>(initial_.size())) {
                record_evaluation(
                    "accepted_iterate",
                    dense->Dim(),
                    dense->Values(),
                    true,
                    true,
                    "",
                    true
                );
            }
        }
        return true;
    }

    void finalize_solution(
        Ipopt::SolverReturn status,
        Ipopt::Index n,
        const Ipopt::Number* x,
        const Ipopt::Number* z_lower,
        const Ipopt::Number* z_upper,
        Ipopt::Index,
        const Ipopt::Number*,
        const Ipopt::Number*,
        Ipopt::Number,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*
    ) override {
        record_evaluation("terminal", n, x, true, true, "", true);
        variables_.assign(x, x + n);
        lower_bound_multipliers_.assign(z_lower, z_lower + n);
        upper_bound_multipliers_.assign(z_upper, z_upper + n);
        solver_converged_ = status == Ipopt::SUCCESS
            || status == Ipopt::STOP_AT_ACCEPTABLE_POINT;
        if (status == Ipopt::SUCCESS) {
            solver_status_ = "success";
        } else if (status == Ipopt::STOP_AT_ACCEPTABLE_POINT) {
            solver_status_ = "acceptable_point";
        } else {
            solver_status_ = "failed_" + std::to_string(static_cast<int>(status));
        }
    }

    [[nodiscard]] Held2SearchRun result() const {
        return {
            solver_converged_,
            solver_status_,
            callback_error_,
            variables_,
            lower_bound_multipliers_,
            upper_bound_multipliers_,
            {},
            evaluation_trace_,
        };
    }

private:
    void record_evaluation(
        const std::string& callback,
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool new_x,
        bool succeeded,
        const std::string& error,
        bool accepted_iterate = false
    ) {
        Held2StageIILocalEvaluation record;
        record.sequence = static_cast<int>(evaluation_trace_.size());
        record.callback = callback;
        record.new_x = new_x;
        record.callback_succeeded = succeeded;
        record.accepted_iterate = accepted_iterate;
        record.error = error;
        if (x != nullptr && n == static_cast<Ipopt::Index>(initial_.size())) {
            record.raw_variables.assign(x, x + n);
            for (Ipopt::Index index = 0; index < n; ++index) {
                const double value = x[index];
                record.maximum_bound_violation = std::max({
                    record.maximum_bound_violation,
                    lower_[static_cast<std::size_t>(index)] - value,
                    value - upper_[static_cast<std::size_t>(index)],
                });
            }
        }
        record.mapping_status = succeeded ? "mapped" : "callback_failed";
        evaluation_trace_.push_back(std::move(record));
    }

    [[nodiscard]] Held2StateEvaluation evaluate(
        Ipopt::Index n,
        const Ipopt::Number* x
    ) const {
        if (n != static_cast<Ipopt::Index>(initial_.size())) {
            throw std::invalid_argument("HELD2 search coordinate count changed");
        }
        const std::vector<double> variables(x, x + n);
        if (cached_evaluation_.has_value() && cached_variables_ == variables) {
            return *cached_evaluation_;
        }
        cached_variables_ = variables;
        cached_evaluation_ = objective_(variables);
        return *cached_evaluation_;
    }

    std::function<Held2StateEvaluation(const std::vector<double>&)> objective_;
    mutable std::vector<double> cached_variables_;
    mutable std::optional<Held2StateEvaluation> cached_evaluation_;
    std::vector<double> initial_;
    std::vector<double> lower_;
    std::vector<double> upper_;
    bool solver_converged_ = false;
    std::string solver_status_ = "not_run";
    std::string callback_error_;
    std::vector<double> variables_;
    std::vector<double> lower_bound_multipliers_;
    std::vector<double> upper_bound_multipliers_;
    std::vector<Held2StageIILocalEvaluation> evaluation_trace_;
    Held2ProgressObserver* observer_ = nullptr;
    std::string progress_stage_;
    int progress_major_ = -1;
    int progress_attempt_ = -1;
};

Held2SearchRun solve_held2_vector_search(
    std::function<Held2StateEvaluation(const std::vector<double>&)> objective,
    const std::vector<double>& initial,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    Held2ProgressObserver* observer = nullptr,
    const std::string& progress_stage = {},
    int progress_major = -1,
    int progress_attempt = -1
) {
    Ipopt::SmartPtr<Held2SearchTnlp> problem = new Held2SearchTnlp(
        std::move(objective),
        initial,
        lower,
        upper,
        observer,
        progress_stage,
        progress_major,
        progress_attempt
    );
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application = IpoptApplicationFactory();
    configure_held2_ipopt(application);
    if (application->Initialize() != Ipopt::Solve_Succeeded) {
        return {};
    }
    application->OptimizeTNLP(problem);
    return problem->result();
}

Held2SearchRun solve_held2_search(
    const Held2StateEvaluator& evaluator,
    const std::vector<double>& reference_variables,
    const Held2StateEvaluation& reference,
    bool use_tpd,
    const std::vector<double>& initial,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    Held2ProgressObserver* observer = nullptr,
    const std::string& progress_stage = {},
    int progress_major = -1,
    int progress_attempt = -1
) {
    Held2SearchObjectiveEvaluator objective(
        evaluator,
        reference_variables,
        reference,
        use_tpd
    );
    return solve_held2_vector_search(
        [objective = std::move(objective)](
            const std::vector<double>& variables
        ) {
            return objective.evaluate(variables);
        },
        initial,
        lower,
        upper,
        observer,
        progress_stage,
        progress_major,
        progress_attempt
    );
}

Held2StageIISimplexChart evaluate_stage_ii_simplex_chart(
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    double sum_upper,
    const std::vector<double>& coordinates,
    bool inverse,
    const std::vector<bool>& charged_coordinates = {},
    double total_ion_mole_fraction_max =
        std::numeric_limits<double>::quiet_NaN()
) {
    const std::size_t dimension = lower.size();
    if (dimension == 0 || upper.size() != dimension
        || coordinates.size() != dimension || !std::isfinite(sum_upper)) {
        throw std::invalid_argument("HELD2 Stage II simplex chart dimensions are invalid");
    }
    require_finite_vector(lower, "HELD2 Stage II independent lower bounds");
    require_finite_vector(upper, "HELD2 Stage II independent upper bounds");
    require_finite_vector(coordinates, "HELD2 Stage II simplex coordinates");
    const double lower_sum = std::accumulate(lower.begin(), lower.end(), 0.0);
    const double free_scale = sum_upper - lower_sum;
    if (!(free_scale > 0.0)) {
        throw std::invalid_argument("HELD2 Stage II simplex chart has no feasible interior");
    }
    for (std::size_t index = 0; index < dimension; ++index) {
        if (lower[index] < 0.0
            || upper[index] < free_scale + lower[index]) {
            throw std::invalid_argument(
                "HELD2 Stage II simplex chart requires redundant upper bounds"
            );
        }
    }

    const bool source_domain_chart =
        std::isfinite(total_ion_mole_fraction_max);
    if (source_domain_chart && charged_coordinates.size() != dimension) {
        throw std::invalid_argument(
            "HELD2 Stage II source-domain chart dimensions are invalid"
        );
    }
    if (!source_domain_chart && !charged_coordinates.empty()) {
        throw std::invalid_argument(
            "HELD2 Stage II source-domain chart requires a finite ion ceiling"
        );
    }
    double charged_lower_sum = 0.0;
    if (source_domain_chart) {
        for (std::size_t index = 0; index < dimension; ++index) {
            if (charged_coordinates[index]) {
                charged_lower_sum += lower[index];
            }
        }
        if (!(total_ion_mole_fraction_max > charged_lower_sum)
            || total_ion_mole_fraction_max > sum_upper) {
            throw std::invalid_argument(
                "HELD2 Stage II Provider ion ceiling is invalid or infeasible"
            );
        }
    }

    Held2StageIISimplexChart result;
    result.chart_coordinates.assign(dimension, 0.0);
    result.physical_coordinates.assign(dimension, 0.0);
    if (source_domain_chart) {
        if (inverse) {
            double remaining_ion = total_ion_mole_fraction_max
                - charged_lower_sum;
            double charged_shift_sum = 0.0;
            for (std::size_t index = 0; index < dimension; ++index) {
                if (!charged_coordinates[index]) {
                    continue;
                }
                const double shifted = coordinates[index] - lower[index];
                if (shifted < 0.0 || shifted > remaining_ion
                    || !(remaining_ion > 0.0)) {
                    throw std::invalid_argument(
                        "HELD2 Stage II physical start is outside the Provider ion domain"
                    );
                }
                result.chart_coordinates[index] = shifted / remaining_ion;
                charged_shift_sum += shifted;
                remaining_ion -= shifted;
            }
            double remaining_composition = free_scale - charged_shift_sum;
            for (std::size_t index = 0; index < dimension; ++index) {
                if (charged_coordinates[index]) {
                    continue;
                }
                const double shifted = coordinates[index] - lower[index];
                if (shifted < 0.0 || shifted > remaining_composition
                    || !(remaining_composition > 0.0)) {
                    throw std::invalid_argument(
                        "HELD2 Stage II physical start is outside the feasible simplex"
                    );
                }
                result.chart_coordinates[index] = shifted / remaining_composition;
                remaining_composition -= shifted;
            }
        } else {
            for (std::size_t index = 0; index < dimension; ++index) {
                const Held2StageIIChartCoordinate normalized =
                    normalize_held2_stage_ii_chart_coordinate(coordinates[index]);
                result.chart_coordinates[index] = normalized.normalized;
                result.normalized_boundary_contact =
                    result.normalized_boundary_contact
                    || normalized.normalized_boundary_contact;
            }
        }

        std::vector<Held2SecondOrderScalar> chart;
        chart.reserve(dimension);
        for (std::size_t index = 0; index < dimension; ++index) {
            chart.push_back(second_order_variable(
                result.chart_coordinates[index], dimension, index
            ));
        }
        std::vector<Held2SecondOrderScalar> shifted(
            dimension, second_order_constant(0.0, dimension)
        );
        Held2SecondOrderScalar remaining_ion = second_order_constant(
            total_ion_mole_fraction_max - charged_lower_sum, dimension
        );
        Held2SecondOrderScalar charged_shift = second_order_constant(0.0, dimension);
        for (std::size_t index = 0; index < dimension; ++index) {
            if (!charged_coordinates[index]) {
                continue;
            }
            shifted[index] = second_order_multiply(remaining_ion, chart[index]);
            charged_shift = second_order_add(charged_shift, shifted[index]);
            remaining_ion = second_order_add(
                remaining_ion, shifted[index], -1.0
            );
        }
        Held2SecondOrderScalar remaining_composition = second_order_add(
            second_order_constant(free_scale, dimension), charged_shift, -1.0
        );
        for (std::size_t index = 0; index < dimension; ++index) {
            if (charged_coordinates[index]) {
                continue;
            }
            shifted[index] = second_order_multiply(
                remaining_composition, chart[index]
            );
            remaining_composition = second_order_add(
                remaining_composition, shifted[index], -1.0
            );
        }

        result.jacobian.assign(dimension * dimension, 0.0);
        result.component_hessians.assign(
            dimension * dimension * dimension, 0.0
        );
        for (std::size_t component = 0; component < dimension; ++component) {
            result.physical_coordinates[component] = lower[component]
                + shifted[component].value;
            for (std::size_t column = 0; column < dimension; ++column) {
                result.jacobian[component * dimension + column] =
                    shifted[component].gradient[column];
                for (std::size_t row = 0; row < dimension; ++row) {
                    result.component_hessians[
                        (component * dimension + row) * dimension + column
                    ] = shifted[component].hessian[row * dimension + column];
                }
            }
            if (std::abs(result.jacobian[component * dimension + component])
                <= std::numeric_limits<double>::epsilon()
                    * std::max(1.0, free_scale)) {
                result.singular = true;
            }
        }
        return result;
    }

    if (inverse) {
        double remaining = free_scale;
        for (std::size_t index = 0; index < dimension; ++index) {
            const double shifted = coordinates[index] - lower[index];
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
        for (std::size_t index = 0; index < dimension; ++index) {
            const Held2StageIIChartCoordinate normalized =
                normalize_held2_stage_ii_chart_coordinate(coordinates[index]);
            result.chart_coordinates[index] = normalized.normalized;
            result.normalized_boundary_contact =
                result.normalized_boundary_contact
                || normalized.normalized_boundary_contact;
        }
        for (std::size_t index = 0; index < dimension; ++index) {
            double product = free_scale * result.chart_coordinates[index];
            for (std::size_t previous = 0; previous < index; ++previous) {
                product *= 1.0 - result.chart_coordinates[previous];
            }
            result.physical_coordinates[index] = lower[index] + product;
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

Held2StateEvaluation transform_stage_ii_reduced_state(
    const Held2StateEvaluation& physical,
    const Held2StageIIPressureRootReduction& reduced,
    const Held2StageIISimplexChart& chart
) {
    const std::size_t dimension = reduced.gradient.size();
    if (reduced.hessian.size() != dimension * dimension
        || chart.jacobian.size() != dimension * dimension
        || chart.component_hessians.size()
            != dimension * dimension * dimension) {
        throw std::invalid_argument(
            "HELD2 Stage II reduced derivative block is incomplete"
        );
    }
    Held2StateEvaluation result = physical;
    result.objective = reduced.objective;
    result.gradient.assign(dimension, 0.0);
    result.hessian.assign(dimension * dimension, 0.0);
    for (std::size_t column = 0; column < dimension; ++column) {
        for (std::size_t physical_index = 0; physical_index < dimension;
             ++physical_index) {
            result.gradient[column] += chart.jacobian[
                physical_index * dimension + column
            ] * reduced.gradient[physical_index];
        }
    }
    for (std::size_t row = 0; row < dimension; ++row) {
        for (std::size_t column = 0; column < dimension; ++column) {
            double value = 0.0;
            for (std::size_t left = 0; left < dimension; ++left) {
                value += reduced.gradient[left] * chart.component_hessians[
                    (left * dimension + row) * dimension + column
                ];
                for (std::size_t right = 0; right < dimension; ++right) {
                    value += chart.jacobian[left * dimension + row]
                        * reduced.hessian[left * dimension + right]
                        * chart.jacobian[right * dimension + column];
                }
            }
            result.hessian[row * dimension + column] = value;
        }
    }
    return result;
}

Held2PressureRoot select_stage_ii_stable_branch(
    const Held2StateEvaluator& evaluator,
    const Held2VolumeBoundsEvaluator& volume_bounds_evaluator,
    const std::vector<double>& independent,
    int stable_branch_index
) {
    if (stable_branch_index < 0) {
        throw std::invalid_argument(
            "HELD2 Stage II stable branch index must be nonnegative"
        );
    }
    const Held2PressureEnvelopeResult envelope =
        evaluate_held2_pressure_envelope(
            independent,
            volume_bounds_evaluator(independent),
            evaluator,
            64,
            8
        );
    if (envelope.outcome != "selected"
        && envelope.failure_reason != "stable_objective_tie") {
        throw std::runtime_error(
            "HELD2 Stage II pressure-root local evaluation failed: "
            + envelope.failure_reason
        );
    }
    int current = 0;
    for (const Held2PressureRoot& root : envelope.roots) {
        if (root.mechanical_class != "strict_stable" || root.boundary) {
            continue;
        }
        if (current++ == stable_branch_index) {
            if (!audit_held2_tolerance(
                    kHeld2RootPressure, root.pressure_residual
                ).passed) {
                throw std::runtime_error(
                    "HELD2 Stage II selected branch failed pressure closure"
                );
            }
            return root;
        }
    }
    throw std::runtime_error(
        "HELD2 Stage II selected stable branch disappeared"
    );
}

Held2PressureRoot polish_stage_ii_stable_branch(
    const Held2StateEvaluator& evaluator,
    const Held2VolumeBoundsEvaluator& volume_bounds_evaluator,
    const std::vector<double>& independent,
    int stable_branch_index,
    double log_volume_hint
) {
    const std::array<double, 2> volume_bounds =
        volume_bounds_evaluator(independent);
    const double lower = std::log(volume_bounds[0]);
    const double upper = std::log(volume_bounds[1]);
    double log_volume = std::clamp(log_volume_hint, lower, upper);
    for (int iteration = 0; iteration < 20; ++iteration) {
        const Held2StateEvaluation state = evaluator(independent, log_volume);
        const double residual = state.pressure_stationarity_relative;
        const double derivative =
            state.pressure_stationarity_derivative_log_volume;
        const double curvature = state.hessian.empty()
            ? std::numeric_limits<double>::quiet_NaN()
            : state.hessian.back();
        if (std::abs(residual) <= kHeld2IpoptTarget.atol
            && derivative < 0.0 && curvature > 0.0
            && audit_held2_tolerance(
                   kHeld2MechanicalMargin,
                   std::min(std::abs(derivative), std::abs(curvature))
               ).passed) {
            Held2PressureRoot root;
            root.log_volume = log_volume;
            root.volume = state.volume;
            root.objective = state.objective;
            root.pressure_residual = residual;
            root.pressure_derivative_log_volume = derivative;
            root.objective_curvature_log_volume = curvature;
            root.mechanical_class = "strict_stable";
            root.origin = "continued_exact_derivative";
            root.state = state;
            return root;
        }
        if (!std::isfinite(residual) || !std::isfinite(derivative)
            || std::abs(derivative)
                <= std::numeric_limits<double>::epsilon()) {
            break;
        }
        const double step = std::clamp(-residual / derivative, -0.75, 0.75);
        const double candidate = std::clamp(
            log_volume + step,
            std::nextafter(lower, upper),
            std::nextafter(upper, lower)
        );
        if (candidate == log_volume) {
            break;
        }
        log_volume = candidate;
    }
    return select_stage_ii_stable_branch(
        evaluator,
        volume_bounds_evaluator,
        independent,
        stable_branch_index
    );
}

Held2SearchRun solve_stage_ii_pressure_root_local(
    const Held2StateEvaluator& evaluator,
    const Held2VolumeBoundsEvaluator& volume_bounds_evaluator,
    const std::vector<double>& feed,
    const std::vector<double>& multipliers,
    const std::vector<double>& independent_initial,
    int stable_branch_index,
    double branch_log_volume_hint,
    const std::vector<double>& independent_lower,
    const std::vector<double>& independent_upper,
    double composition_sum_upper,
    const std::vector<bool>& charged_coordinates,
    double total_ion_mole_fraction_max,
    Held2ProgressObserver* observer = nullptr,
    int major_iteration = -1,
    int attempt_index = -1
) {
    const std::size_t dimension = feed.size();
    Held2SearchRun invalid;
    if (dimension == 0 || multipliers.size() != dimension
        || independent_initial.size() != dimension
        || independent_lower.size() != dimension
        || independent_upper.size() != dimension
        || stable_branch_index < 0 || !std::isfinite(branch_log_volume_hint)) {
        invalid.callback_error =
            "HELD2 Stage II pressure-root search dimensions changed";
        return invalid;
    }
    Held2StageIISimplexChart initial_chart;
    try {
        initial_chart = evaluate_stage_ii_simplex_chart(
            independent_lower,
            independent_upper,
            composition_sum_upper,
            independent_initial,
            true,
            charged_coordinates,
            total_ion_mole_fraction_max
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
    std::vector<double> chart_lower(dimension, 0.0);
    std::vector<double> chart_upper(dimension, 1.0);
    const auto chart_evaluator = [
        evaluator,
        volume_bounds_evaluator,
        feed,
        multipliers,
        stable_branch_index,
        branch_log_volume_hint,
        independent_lower,
        independent_upper,
        composition_sum_upper,
        charged_coordinates,
        total_ion_mole_fraction_max
    ](const std::vector<double>& chart_coordinates) {
        const Held2StageIISimplexChart chart = evaluate_stage_ii_simplex_chart(
            independent_lower,
            independent_upper,
            composition_sum_upper,
            chart_coordinates,
            false,
            charged_coordinates,
            total_ion_mole_fraction_max
        );
        const Held2PressureRoot root = polish_stage_ii_stable_branch(
            evaluator,
            volume_bounds_evaluator,
            chart.physical_coordinates,
            stable_branch_index,
            branch_log_volume_hint
        );
        const Held2StageIIPressureRootReduction reduced =
            reduce_held2_stage_ii_pressure_root(
                chart.physical_coordinates,
                feed,
                multipliers,
                root.state
            );
        return transform_stage_ii_reduced_state(
            root.state,
            reduced,
            chart
        );
    };
    Held2SearchRun run = solve_held2_vector_search(
        chart_evaluator,
        chart_initial,
        chart_lower,
        chart_upper,
        observer,
        "STAGE II LOCAL IPOPT",
        major_iteration,
        attempt_index
    );
    for (Held2StageIILocalEvaluation& evaluation : run.evaluation_trace) {
        if (evaluation.raw_variables.size() != dimension) {
            evaluation.mapping_status = "variables_unavailable";
            continue;
        }
        try {
            const Held2StageIISimplexChart mapped = evaluate_stage_ii_simplex_chart(
                independent_lower,
                independent_upper,
                composition_sum_upper,
                evaluation.raw_variables,
                false,
                charged_coordinates,
                total_ion_mole_fraction_max
            );
            evaluation.physical_variables = mapped.physical_coordinates;
            const Held2PressureRoot root = polish_stage_ii_stable_branch(
                evaluator,
                volume_bounds_evaluator,
                mapped.physical_coordinates,
                stable_branch_index,
                branch_log_volume_hint
            );
            evaluation.physical_variables.push_back(root.log_volume);
            evaluation.mapping_status = mapped.normalized_boundary_contact
                ? "normalized_boundary_contact"
                : "mapped";
        } catch (const std::exception&) {
            evaluation.mapping_status = "chart_out_of_bounds";
        }
    }
    if (run.variables.size() != dimension) {
        return run;
    }
    try {
        const Held2StageIISimplexChart physical = evaluate_stage_ii_simplex_chart(
            independent_lower,
            independent_upper,
            composition_sum_upper,
            run.variables,
            false,
            charged_coordinates,
            total_ion_mole_fraction_max
        );
        if (physical.singular) {
            throw std::invalid_argument("HELD2 Stage II simplex terminal is singular");
        }
        const Held2PressureRoot root = polish_stage_ii_stable_branch(
            evaluator,
            volume_bounds_evaluator,
            physical.physical_coordinates,
            stable_branch_index,
            branch_log_volume_hint
        );
        run.variables = physical.physical_coordinates;
        run.variables.push_back(root.log_volume);
        run.lower_bound_multipliers.push_back(0.0);
        run.upper_bound_multipliers.push_back(0.0);
        run.coordinate_jacobian = physical.jacobian;
    } catch (const std::exception& error) {
        run.callback_error = error.what();
        run.solver_converged = false;
    }
    return run;
}

Held2PhysicalKkt evaluate_stage_ii_physical_kkt(
    const std::vector<double>& physical_gradient,
    const std::vector<double>& master_multipliers,
    const std::vector<double>& variables,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    double composition_sum_upper,
    const std::vector<double>& chart_jacobian,
    const std::vector<double>& chart_lower_multipliers,
    const std::vector<double>& chart_upper_multipliers
) {
    const std::size_t dimension = master_multipliers.size();
    if (physical_gradient.size() != dimension + 1
        || variables.size() != dimension + 1
        || lower.size() != dimension + 1
        || upper.size() != dimension + 1
        || chart_jacobian.size() != dimension * dimension
        || chart_lower_multipliers.size() != dimension + 1
        || chart_upper_multipliers.size() != dimension + 1) {
        throw std::invalid_argument("HELD2 Stage II physical KKT dimensions changed");
    }
    Held2PhysicalKkt result;
    result.dual_sign_violation_inf_norm = 0.0;
    for (double value : chart_lower_multipliers) {
        result.dual_sign_violation_inf_norm = std::max(
            result.dual_sign_violation_inf_norm, std::max(0.0, -value)
        );
    }
    for (double value : chart_upper_multipliers) {
        result.dual_sign_violation_inf_norm = std::max(
            result.dual_sign_violation_inf_norm, std::max(0.0, -value)
        );
    }
    result.dual_signs_valid = audit_held2_tolerance(
        kHeld2Stage2DualSign, result.dual_sign_violation_inf_norm
    ).passed;
    if (!result.dual_signs_valid) {
        return result;
    }
    result.primal_inf_norm = 0.0;
    for (std::size_t index = 0; index < variables.size(); ++index) {
        if (lower[index] > upper[index]) {
            return result;
        }
        result.primal_inf_norm = std::max({
            result.primal_inf_norm,
            std::max(0.0, lower[index] - variables[index]),
            std::max(0.0, variables[index] - upper[index]),
        });
    }
    const double simplex_slack = composition_sum_upper
        - std::accumulate(variables.begin(), variables.end() - 1, 0.0);
    result.primal_inf_norm = std::max(
        result.primal_inf_norm, std::max(0.0, -simplex_slack)
    );
    if (!audit_held2_tolerance(kHeld2Stage2Primal, result.primal_inf_norm).passed) {
        return result;
    }

    std::vector<double> chart_force(dimension + 1, 0.0);
    for (std::size_t index = 0; index <= dimension; ++index) {
        chart_force[index] = chart_upper_multipliers[index]
            - chart_lower_multipliers[index];
    }
    std::vector<double> transposed_jacobian(dimension * dimension, 0.0);
    for (std::size_t row = 0; row < dimension; ++row) {
        for (std::size_t column = 0; column < dimension; ++column) {
            transposed_jacobian[row * dimension + column] =
                chart_jacobian[column * dimension + row];
        }
    }
    const std::optional<std::vector<double>> composition_force =
        solve_dense_linear_system(
            transposed_jacobian,
            std::vector<double>(chart_force.begin(), chart_force.end() - 1),
            dimension
        );
    if (!composition_force.has_value()) {
        return result;
    }
    std::vector<double> physical_force = *composition_force;
    physical_force.push_back(0.0);
    physical_force.back() = chart_force.back();

    const auto multiplier_cap = [](double slack) {
        return slack <= 0.0
            ? std::numeric_limits<double>::infinity()
            : kHeld2Stage2Complementarity.atol / slack;
    };
    double simplex_lower = 0.0;
    double simplex_upper = multiplier_cap(simplex_slack);
    for (std::size_t index = 0; index < dimension; ++index) {
        simplex_lower = std::max(
            simplex_lower,
            physical_force[index]
                - multiplier_cap(upper[index] - variables[index])
        );
        simplex_upper = std::min(
            simplex_upper,
            physical_force[index]
                + multiplier_cap(variables[index] - lower[index])
        );
    }
    if (simplex_lower > simplex_upper + kHeld2Stage2DualPullback.atol) {
        return result;
    }
    const double simplex_multiplier = std::isfinite(simplex_upper)
        ? 0.5 * (simplex_lower + simplex_upper)
        : simplex_lower;
    std::vector<double> reconstructed(dimension + 1, 0.0);
    double complementarity = std::abs(simplex_multiplier * simplex_slack);
    for (std::size_t index = 0; index < dimension; ++index) {
        const double difference = physical_force[index] - simplex_multiplier;
        const double lower_multiplier = std::max(0.0, -difference);
        const double upper_multiplier = std::max(0.0, difference);
        reconstructed[index] = -lower_multiplier + upper_multiplier
            + simplex_multiplier;
        complementarity = std::max({
            complementarity,
            std::abs(lower_multiplier * (variables[index] - lower[index])),
            std::abs(upper_multiplier * (upper[index] - variables[index])),
        });
    }
    const double phase_lower_multiplier = std::max(0.0, -physical_force.back());
    const double phase_upper_multiplier = std::max(0.0, physical_force.back());
    reconstructed.back() = -phase_lower_multiplier + phase_upper_multiplier;
    complementarity = std::max({
        complementarity,
        std::abs(phase_lower_multiplier * (variables.back() - lower.back())),
        std::abs(phase_upper_multiplier * (upper.back() - variables.back())),
    });

    double reconstruction = maximum_abs_difference(physical_force, reconstructed);
    for (std::size_t column = 0; column < dimension; ++column) {
        double pulled = 0.0;
        for (std::size_t row = 0; row < dimension; ++row) {
            pulled += chart_jacobian[row * dimension + column]
                * reconstructed[row];
        }
        reconstruction = std::max(
            reconstruction,
            std::abs(chart_force[column] - pulled)
        );
    }
    reconstruction = std::max(
        reconstruction,
        std::abs(chart_force.back() - reconstructed.back())
    );
    double reconstruction_scale = 0.0;
    for (double value : physical_force) {
        reconstruction_scale = std::max(reconstruction_scale, std::abs(value));
    }
    for (double value : reconstructed) {
        reconstruction_scale = std::max(reconstruction_scale, std::abs(value));
    }
    double stationarity = 0.0;
    for (std::size_t index = 0; index < dimension; ++index) {
        stationarity = std::max(
            stationarity,
            std::abs(
                physical_gradient[index] - master_multipliers[index]
                + reconstructed[index]
            )
        );
    }
    stationarity = std::max(
        stationarity,
        std::abs(physical_gradient.back() + reconstructed.back())
    );
    result.stationarity_inf_norm = stationarity;
    result.complementarity = complementarity;
    result.reconstruction_inf_norm = reconstruction;
    result.reconstruction_scale = reconstruction_scale;
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
        if (std::abs(compositions[index] - feed_composition)
            <= kHeld2BasinDuplicateComposition.atol) {
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

}  // namespace

Held2StageIIPressureRootReduction reduce_held2_stage_ii_pressure_root(
    const std::vector<double>& independent,
    const std::vector<double>& feed,
    const std::vector<double>& multipliers,
    const Held2StateEvaluation& state
) {
    const std::size_t dimension = independent.size();
    const std::size_t full_dimension = dimension + 1;
    if (dimension == 0 || feed.size() != dimension
        || multipliers.size() != dimension
        || state.gradient.size() != full_dimension
        || state.hessian.size() != full_dimension * full_dimension) {
        throw std::invalid_argument(
            "HELD2 Stage II pressure-root reduction dimensions changed"
        );
    }
    require_finite_vector(independent, "HELD2 Stage II independent composition");
    require_finite_vector(feed, "HELD2 Stage II feed composition");
    require_finite_vector(multipliers, "HELD2 Stage II multipliers");
    require_finite_vector(state.gradient, "HELD2 Stage II physical gradient");
    require_finite_vector(state.hessian, "HELD2 Stage II physical Hessian");
    if (!std::isfinite(state.objective)) {
        throw std::invalid_argument(
            "HELD2 Stage II physical objective must be finite"
        );
    }
    const double pressure_curvature = state.hessian.back();
    if (!audit_held2_tolerance(
            kHeld2MechanicalMargin, pressure_curvature
        ).passed) {
        throw std::domain_error(
            "HELD2 Stage II pressure root is not strictly mechanically stable"
        );
    }

    Held2StageIIPressureRootReduction result;
    result.objective = state.objective;
    result.gradient.assign(dimension, 0.0);
    result.hessian.assign(dimension * dimension, 0.0);
    result.pressure_coordinate_gradient = state.gradient.back();
    result.pressure_coordinate_curvature = pressure_curvature;
    for (std::size_t row = 0; row < dimension; ++row) {
        result.objective += multipliers[row]
            * (feed[row] - independent[row]);
        result.gradient[row] = state.gradient[row] - multipliers[row];
        for (std::size_t column = 0; column < dimension; ++column) {
            result.hessian[row * dimension + column] =
                state.hessian[row * full_dimension + column]
                - state.hessian[row * full_dimension + dimension]
                    * state.hessian[dimension * full_dimension + column]
                    / pressure_curvature;
        }
    }
    return result;
}

Held2StageIIStep5Assessment assess_held2_stage_ii_step5(
    double upper_bound,
    double local_value,
    bool local_state_certified
) {
    Held2StageIIStep5Assessment result;
    result.gap = upper_bound - local_value;
    if (!std::isfinite(upper_bound) || !std::isfinite(local_value)) {
        result.reason = "nonfinite_lower_upper_comparison";
        return result;
    }
    if (!local_state_certified) {
        result.reason = "local_state_not_certified";
        return result;
    }
    if (local_value <= upper_bound) {
        result.qualified = true;
        result.reason = "lower_not_above_upper";
        return result;
    }
    if (audit_held2_tolerance(kHeld2Step6Gap, result.gap).passed) {
        result.qualified = true;
        result.reason = "lower_equal_within_step6_gate";
        return result;
    }
    result.reason = "certified_local_above_upper";
    return result;
}

Held2StageIIChartCoordinate normalize_held2_stage_ii_chart_coordinate(
    double coordinate
) {
    if (!std::isfinite(coordinate)) {
        throw std::invalid_argument("HELD2 Stage II chart coordinate is not finite");
    }
    if (coordinate >= 0.0 && coordinate <= 1.0) {
        return {coordinate, coordinate, false};
    }
    if (coordinate < 0.0
        && audit_held2_tolerance(kHeld2ChartContact, coordinate).passed) {
        return {coordinate, 0.0, true};
    }
    if (coordinate > 1.0
        && audit_held2_tolerance(kHeld2ChartContact, coordinate - 1.0).passed) {
        return {coordinate, 1.0, true};
    }
    throw std::invalid_argument("HELD2 Stage II chart coordinate is outside [0, 1]");
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

Held2StateEvaluation evaluate_held2_manufactured_search_objective(
    const Held2Coordinates& coordinates,
    const std::vector<double>& variables,
    const std::vector<double>& reference_variables,
    bool use_tpd
) {
    if (reference_variables.empty()) {
        throw std::invalid_argument("HELD2 search reference variables must not be empty");
    }
    const Held2StateEvaluator phase_evaluator = [&coordinates](
        const std::vector<double>& independent,
        double log_volume
    ) {
        return evaluate_manufactured_state_impl(
            coordinates,
            independent,
            log_volume
        );
    };
    const std::vector<double> reference_independent(
        reference_variables.begin(),
        reference_variables.end() - 1
    );
    Held2SearchObjectiveEvaluator objective(
        phase_evaluator,
        reference_variables,
        phase_evaluator(reference_independent, reference_variables.back()),
        use_tpd
    );
    return objective.evaluate(variables);
}

Held2PressureEnvelopeResult evaluate_held2_pressure_envelope(
    const std::vector<double>& independent_modified_fractions,
    const std::array<double, 2>& molar_volume_bounds,
    const Held2StateEvaluator& evaluator,
    int initial_interval_count,
    int maximum_subdivision_depth
) {
    if (independent_modified_fractions.empty()) {
        throw std::invalid_argument(
            "HELD2 pressure envelope requires independent composition coordinates"
        );
    }
    require_finite_vector(
        independent_modified_fractions,
        "HELD2 pressure-envelope composition"
    );
    if (!std::isfinite(molar_volume_bounds[0])
        || !std::isfinite(molar_volume_bounds[1])
        || molar_volume_bounds[0] <= 0.0
        || molar_volume_bounds[1] <= molar_volume_bounds[0]) {
        throw std::invalid_argument(
            "HELD2 pressure-envelope volume bounds must be finite, positive, and ordered"
        );
    }
    if (initial_interval_count < 1 || maximum_subdivision_depth < 0) {
        throw std::invalid_argument(
            "HELD2 pressure-envelope scan controls must be nonnegative"
        );
    }

    struct EvaluatedPoint {
        Held2PressureScanPoint diagnostic;
        Held2StateEvaluation state;
    };
    struct WorkingInterval {
        double lower = 0.0;
        double upper = 0.0;
        int depth = 0;
    };
    struct Bracket {
        double lower = 0.0;
        double upper = 0.0;
    };

    Held2PressureEnvelopeResult result;
    std::map<double, EvaluatedPoint> cache;
    const auto evaluate_point = [
        &cache,
        &evaluator,
        &independent_modified_fractions,
        &result
    ](double log_volume) -> const EvaluatedPoint& {
        const auto found = cache.find(log_volume);
        if (found != cache.end()) {
            return found->second;
        }
        EvaluatedPoint point;
        point.diagnostic.log_volume = log_volume;
        point.diagnostic.volume = std::exp(log_volume);
        try {
            point.state = evaluator(
                independent_modified_fractions,
                log_volume
            );
            if (!std::isfinite(point.state.volume)
                || point.state.volume <= 0.0
                || !std::isfinite(point.state.objective)
                || !std::isfinite(point.state.pressure_stationarity_relative)
                || !std::isfinite(
                    point.state.pressure_stationarity_derivative_log_volume
                )
                || point.state.hessian.empty()
                || !std::isfinite(point.state.hessian.back())) {
                throw std::invalid_argument(
                    "HELD2 pressure-envelope evaluation returned incomplete evidence"
                );
            }
            point.diagnostic.volume = point.state.volume;
            point.diagnostic.pressure_residual =
                point.state.pressure_stationarity_relative;
            point.diagnostic.pressure_derivative_log_volume =
                point.state.pressure_stationarity_derivative_log_volume;
            point.diagnostic.objective = point.state.objective;
            point.diagnostic.valid = true;
        } catch (const std::exception& error) {
            point.diagnostic.failure = error.what();
            ++result.evaluation_failure_count;
        }
        return cache.emplace(log_volume, std::move(point)).first->second;
    };
    const auto brackets_zero = [](double left, double right) {
        return (left <= 0.0 && right >= 0.0)
            || (left >= 0.0 && right <= 0.0);
    };
    const auto pressure_bracketed = [&brackets_zero](
        const EvaluatedPoint& left,
        const EvaluatedPoint& right
    ) {
        return left.diagnostic.valid && right.diagnostic.valid
            && brackets_zero(
                left.diagnostic.pressure_residual,
                right.diagnostic.pressure_residual
            );
    };
    const auto stationary_bracketed = [&brackets_zero](
        const EvaluatedPoint& left,
        const EvaluatedPoint& right
    ) {
        return left.diagnostic.valid && right.diagnostic.valid
            && brackets_zero(
                left.diagnostic.pressure_derivative_log_volume,
                right.diagnostic.pressure_derivative_log_volume
            );
    };

    const double lower_log_volume = std::log(molar_volume_bounds[0]);
    const double upper_log_volume = std::log(molar_volume_bounds[1]);
    result.lower_log_volume = lower_log_volume;
    result.upper_log_volume = upper_log_volume;
    const double span = upper_log_volume - lower_log_volume;
    std::vector<WorkingInterval> pending;
    pending.reserve(static_cast<std::size_t>(initial_interval_count));
    for (int interval = 0; interval < initial_interval_count; ++interval) {
        const double lower = interval == 0
            ? lower_log_volume
            : lower_log_volume + span * static_cast<double>(interval)
                / static_cast<double>(initial_interval_count);
        const double upper = interval + 1 == initial_interval_count
            ? upper_log_volume
            : lower_log_volume + span * static_cast<double>(interval + 1)
                / static_cast<double>(initial_interval_count);
        pending.push_back({lower, upper, 0});
    }

    std::vector<Bracket> pressure_brackets;
    std::vector<Bracket> stationary_brackets;
    std::size_t cursor = 0;
    while (cursor < pending.size()) {
        const WorkingInterval interval = pending[cursor++];
        const double midpoint = 0.5 * (interval.lower + interval.upper);
        const EvaluatedPoint& left = evaluate_point(interval.lower);
        const EvaluatedPoint& middle = evaluate_point(midpoint);
        const EvaluatedPoint& right = evaluate_point(interval.upper);
        if (!left.diagnostic.valid || !middle.diagnostic.valid
            || !right.diagnostic.valid) {
            result.intervals.push_back({
                interval.lower,
                interval.upper,
                interval.depth,
                "invalid_evaluation",
            });
            continue;
        }

        const bool left_pressure = pressure_bracketed(left, middle);
        const bool right_pressure = pressure_bracketed(middle, right);
        const bool left_stationary = stationary_bracketed(left, middle);
        const bool right_stationary = stationary_bracketed(middle, right);
        const auto joint_root = [](const EvaluatedPoint& point) {
            return audit_held2_tolerance(
                       kHeld2RootPressure,
                       point.diagnostic.pressure_residual
                   ).passed
                && audit_held2_tolerance(
                       kHeld2RootStationary,
                       point.diagnostic.pressure_derivative_log_volume
                   ).passed;
        };
        const bool left_joint_root = joint_root(left) || joint_root(middle);
        const bool right_joint_root = joint_root(middle) || joint_root(right);
        const bool ambiguous =
            (left_pressure && left_stationary && !left_joint_root)
            || (right_pressure && right_stationary && !right_joint_root);
        if (ambiguous && interval.depth < maximum_subdivision_depth) {
            pending.push_back({interval.lower, midpoint, interval.depth + 1});
            pending.push_back({midpoint, interval.upper, interval.depth + 1});
            continue;
        }
        if (ambiguous) {
            ++result.refinement_failure_count;
            result.intervals.push_back({
                interval.lower,
                interval.upper,
                interval.depth,
                "unresolved_multiple_topology",
            });
            continue;
        }
        if (left_pressure) {
            pressure_brackets.push_back({interval.lower, midpoint});
        }
        if (right_pressure) {
            pressure_brackets.push_back({midpoint, interval.upper});
        }
        if (left_stationary) {
            stationary_brackets.push_back({interval.lower, midpoint});
        }
        if (right_stationary) {
            stationary_brackets.push_back({midpoint, interval.upper});
        }
        std::string status = "root_free_finite_scan";
        if ((left_pressure || right_pressure)
            && (left_stationary || right_stationary)) {
            status = "pressure_and_stationary_brackets";
        } else if (left_pressure || right_pressure) {
            status = "pressure_bracket";
        } else if (left_stationary || right_stationary) {
            status = "stationary_bracket";
        }
        result.intervals.push_back({
            interval.lower,
            interval.upper,
            interval.depth,
            std::move(status),
        });
    }

    const auto refine = [
        &evaluate_point,
        &brackets_zero,
        &result
    ](Bracket bracket, bool stationary) -> std::optional<EvaluatedPoint> {
        const auto residual = [stationary](const EvaluatedPoint& point) {
            return stationary
                ? point.diagnostic.pressure_derivative_log_volume
                : point.diagnostic.pressure_residual;
        };
        EvaluatedPoint left = evaluate_point(bracket.lower);
        EvaluatedPoint right = evaluate_point(bracket.upper);
        if (!left.diagnostic.valid || !right.diagnostic.valid
            || !brackets_zero(residual(left), residual(right))) {
            ++result.refinement_failure_count;
            return std::nullopt;
        }
        const double tolerance = stationary
            ? kHeld2RootStationary.atol
            : kHeld2RootPressure.atol;
        EvaluatedPoint best = std::abs(residual(left)) <= std::abs(residual(right))
            ? left
            : right;
        for (int iteration = 0; iteration < 120; ++iteration) {
            if (audit_held2_tolerance(
                    kHeld2RootLogVolumeWidth,
                    bracket.upper - bracket.lower
                ).passed
                && (stationary || std::abs(residual(best)) <= tolerance)) {
                return best;
            }
            const double midpoint = 0.5 * (bracket.lower + bracket.upper);
            if (!(midpoint > bracket.lower) || !(midpoint < bracket.upper)) {
                if (std::abs(residual(best)) <= tolerance) {
                    return best;
                }
                ++result.refinement_failure_count;
                return std::nullopt;
            }
            const EvaluatedPoint middle = evaluate_point(midpoint);
            if (!middle.diagnostic.valid) {
                ++result.refinement_failure_count;
                return std::nullopt;
            }
            if (std::abs(residual(middle)) < std::abs(residual(best))) {
                best = middle;
            }
            if (brackets_zero(residual(left), residual(middle))) {
                bracket.upper = midpoint;
                right = middle;
            } else {
                bracket.lower = midpoint;
                left = middle;
            }
        }
        ++result.refinement_failure_count;
        return std::nullopt;
    };

    std::vector<Held2PressureRoot> candidates;
    const auto make_root = [lower_log_volume, upper_log_volume](
        const EvaluatedPoint& point,
        std::string origin
    ) {
        Held2PressureRoot root;
        root.log_volume = point.diagnostic.log_volume;
        root.volume = point.state.volume;
        root.objective = point.state.objective;
        root.pressure_residual = point.state.pressure_stationarity_relative;
        root.pressure_derivative_log_volume =
            point.state.pressure_stationarity_derivative_log_volume;
        root.objective_curvature_log_volume = point.state.hessian.back();
        root.origin = std::move(origin);
        root.boundary = audit_held2_tolerance(
                            kHeld2RootBoundary,
                            root.log_volume - lower_log_volume
                        ).passed
            || audit_held2_tolerance(
                   kHeld2RootBoundary,
                   upper_log_volume - root.log_volume
               ).passed;
        root.state = point.state;
        return root;
    };
    for (const Bracket& bracket : pressure_brackets) {
        const std::optional<EvaluatedPoint> root = refine(bracket, false);
        if (root.has_value()) {
            candidates.push_back(make_root(*root, "sign_change"));
        }
    }

    std::vector<double> stationary_points;
    for (const Bracket& bracket : stationary_brackets) {
        const std::optional<EvaluatedPoint> stationary = refine(bracket, true);
        if (!stationary.has_value()) {
            continue;
        }
        const bool duplicate = std::any_of(
            stationary_points.begin(),
            stationary_points.end(),
            [&stationary](double retained) {
                return std::abs(
                    retained - stationary->diagnostic.log_volume
                ) <= kHeld2RootDuplicate.atol;
            }
        );
        if (duplicate) {
            continue;
        }
        stationary_points.push_back(stationary->diagnostic.log_volume);
        ++result.stationary_point_count;
        if (audit_held2_tolerance(
                kHeld2RootPressure,
                stationary->diagnostic.pressure_residual
            ).passed) {
            ++result.tangential_root_count;
            candidates.push_back(make_root(*stationary, "tangential"));
        }
    }

    std::sort(
        candidates.begin(),
        candidates.end(),
        [](const Held2PressureRoot& left, const Held2PressureRoot& right) {
            return left.log_volume < right.log_volume;
        }
    );
    for (Held2PressureRoot& candidate : candidates) {
        if (!result.roots.empty()
            && std::abs(
                result.roots.back().log_volume - candidate.log_volume
            ) <= kHeld2RootDuplicate.atol) {
            result.roots.back().boundary = result.roots.back().boundary
                || candidate.boundary;
            if (candidate.origin == "tangential") {
                result.roots.back().origin = candidate.origin;
            }
            ++result.deduplicated_root_count;
            continue;
        }
        result.roots.push_back(std::move(candidate));
    }

    std::vector<int> stable_indices;
    bool unclassified = false;
    for (std::size_t index = 0; index < result.roots.size(); ++index) {
        Held2PressureRoot& root = result.roots[index];
        if (root.boundary) {
            ++result.boundary_root_count;
        }
        const double mechanical_margin = std::min(
            std::abs(root.pressure_derivative_log_volume),
            std::abs(root.objective_curvature_log_volume)
        );
        const bool mechanical_margin_passed = audit_held2_tolerance(
            kHeld2MechanicalMargin,
            mechanical_margin
        ).passed;
        const bool marginal = !mechanical_margin_passed;
        if (marginal) {
            root.mechanical_class = "marginal";
            ++result.marginal_root_count;
        } else if (root.pressure_derivative_log_volume < 0.0
            && root.objective_curvature_log_volume > 0.0) {
            root.mechanical_class = "strict_stable";
            stable_indices.push_back(static_cast<int>(index));
        } else if (root.pressure_derivative_log_volume > 0.0
            && root.objective_curvature_log_volume < 0.0) {
            root.mechanical_class = "unstable";
        } else {
            unclassified = true;
        }
    }

    if (result.evaluation_failure_count != 0) {
        result.failure_reason = "evaluation_failed";
    } else if (result.refinement_failure_count != 0) {
        result.failure_reason = "refinement_failed";
    } else if (result.boundary_root_count != 0) {
        result.failure_reason = "boundary_root";
    } else if (result.marginal_root_count != 0) {
        result.failure_reason = "marginal_root";
    } else if (unclassified) {
        result.failure_reason = "mechanical_classification_failed";
    } else if (result.roots.empty()) {
        result.failure_reason = "root_not_found";
    } else if (stable_indices.empty()) {
        result.failure_reason = "strict_stable_root_not_found";
    } else {
        double lowest_objective = std::numeric_limits<double>::infinity();
        for (int index : stable_indices) {
            lowest_objective = std::min(
                lowest_objective,
                result.roots[static_cast<std::size_t>(index)].objective
            );
        }
        std::vector<int> lowest;
        for (int index : stable_indices) {
            const double objective =
                result.roots[static_cast<std::size_t>(index)].objective;
            const double scale = std::max(std::abs(objective), std::abs(lowest_objective));
            if (audit_held2_tolerance(
                    kHeld2StableObjectiveTie,
                    objective - lowest_objective,
                    scale
                ).passed) {
                lowest.push_back(index);
            }
        }
        if (lowest.size() != 1) {
            result.objective_tie_count = static_cast<int>(lowest.size()) - 1;
            result.failure_reason = "stable_objective_tie";
        } else {
            result.selected_root_index = lowest.front();
            result.outcome = "selected";
        }
    }

    std::sort(
        result.intervals.begin(),
        result.intervals.end(),
        [](const Held2PressureScanInterval& left,
           const Held2PressureScanInterval& right) {
            return left.lower_log_volume < right.lower_log_volume;
        }
    );
    result.scan_points.reserve(cache.size());
    for (const auto& [log_volume, point] : cache) {
        static_cast<void>(log_volume);
        result.scan_points.push_back(point.diagnostic);
    }
    return result;
}

Held2PressureEnvelopeResult evaluate_held2_manufactured_pressure_envelope(
    const std::string& topology,
    double composition,
    int initial_interval_count
) {
    if (!std::isfinite(composition) || composition < 0.0 || composition > 1.0) {
        throw std::invalid_argument(
            "manufactured pressure-envelope composition must be in [0, 1]"
        );
    }
    const Held2StateEvaluator evaluator = [topology, composition](
        const std::vector<double>& independent,
        double log_volume
    ) {
        if (independent.size() != 1 || independent.front() != composition) {
            throw std::invalid_argument(
                "manufactured pressure-envelope composition changed"
            );
        }
        if (topology == "invalid" && log_volume >= 0.1 && log_volume <= 0.2) {
            throw std::domain_error("manufactured invalid Provider interval");
        }

        double phase_gradient = 0.0;
        double curvature = 0.0;
        double objective = 0.0;
        if (topology == "one_root" || topology == "node_root"
            || topology == "boundary" || topology == "invalid") {
            double root = 0.15;
            if (topology == "node_root") {
                root = 0.0;
            } else if (topology == "boundary") {
                root = -1.5;
            } else if (topology == "invalid") {
                root = -0.5;
            }
            phase_gradient = log_volume - root;
            curvature = 1.0;
            objective = 0.5 * phase_gradient * phase_gradient;
        } else {
            double first = -1.0;
            double middle = -0.2;
            double third = 1.0;
            if (topology == "tied") {
                middle = 0.0;
            } else if (topology == "close_roots") {
                first = -0.02;
                middle = 0.0;
                third = 0.02;
            } else if (topology == "branch_switch") {
                middle = composition < 0.5 ? 0.2 : -0.2;
            } else if (topology == "tangential") {
                phase_gradient = log_volume * log_volume
                    * (log_volume - 0.6);
                curvature = 3.0 * log_volume * log_volume
                    - 1.2 * log_volume;
                objective = std::pow(log_volume, 4) / 4.0
                    - 0.2 * std::pow(log_volume, 3);
            } else if (topology != "three_root") {
                throw std::invalid_argument(
                    "unknown manufactured pressure-envelope topology"
                );
            }
            if (topology != "tangential") {
                const double sum = first + middle + third;
                const double pair_sum = first * middle + first * third
                    + middle * third;
                const double product = first * middle * third;
                phase_gradient = (log_volume - first)
                    * (log_volume - middle) * (log_volume - third);
                curvature = 3.0 * log_volume * log_volume
                    - 2.0 * sum * log_volume + pair_sum;
                objective = std::pow(log_volume, 4) / 4.0
                    - sum * std::pow(log_volume, 3) / 3.0
                    + pair_sum * log_volume * log_volume / 2.0
                    - product * log_volume;
            }
        }

        Held2StateEvaluation state;
        state.modified_fractions = {1.0 - composition, composition};
        state.physical_amounts = state.modified_fractions;
        state.volume = std::exp(log_volume);
        state.objective = objective;
        state.gradient = {0.0, phase_gradient};
        state.hessian = {1.0, 0.0, 0.0, curvature};
        state.modified_potentials = {0.0, 0.0};
        state.pressure_stationarity_relative = -phase_gradient;
        state.pressure_stationarity_derivative_log_volume = -curvature;
        return state;
    };
    return evaluate_held2_pressure_envelope(
        {composition},
        {std::exp(-1.5), std::exp(1.5)},
        evaluator,
        initial_interval_count,
        8
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
        result.independent_lower_bounds.push_back(kModifiedLowerScale * factor);
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
    if (!audit_held2_tolerance(kHeld2CompositionSum, total - 1.0).passed) {
        throw std::invalid_argument("physical mole fractions must sum to one");
    }
    if (!audit_held2_tolerance(
            kHeld2ChargeBalance,
            charge_residual(coordinates.charges, physical_fractions),
            charge_scale(coordinates.charges, physical_fractions)
        ).passed) {
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
    if (!audit_held2_tolerance(
            kHeld2CompositionSum, modified_total - 1.0
        ).passed) {
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
    if (!audit_held2_tolerance(
            kHeld2CompositionSum, modified_total - 1.0
        ).passed) {
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
    if (eliminated < 0.0
        && !audit_held2_tolerance(kHeld2ReconstructedIon, eliminated).passed) {
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

std::vector<double> held2_map_unit_cube_to_independent_fractions(
    const Held2Coordinates& coordinates,
    const std::vector<double>& unit_cube_coordinates,
    double total_ion_mole_fraction_max
) {
    const std::size_t dimension = coordinates.independent_indices.size();
    if (dimension == 0 || unit_cube_coordinates.size() != dimension
        || coordinates.independent_lower_bounds.size() != dimension
        || coordinates.independent_upper_bounds.size() != dimension) {
        throw std::invalid_argument("HELD2 simplex chart dimensions are invalid");
    }
    require_finite_vector(unit_cube_coordinates, "HELD2 simplex cube coordinates");
    const auto dependent_retained = static_cast<std::size_t>(
        std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            coordinates.dependent_index
        ) - coordinates.retained_indices.begin()
    );
    const double composition_sum_upper = 1.0 - kModifiedLowerScale
        * coordinates.modified_factors[dependent_retained];
    const double lower_sum = std::accumulate(
        coordinates.independent_lower_bounds.begin(),
        coordinates.independent_lower_bounds.end(),
        0.0
    );
    double remaining_composition = composition_sum_upper - lower_sum;
    if (!(remaining_composition > 0.0)) {
        throw std::invalid_argument("HELD2 simplex chart has no feasible interior");
    }
    double charged_lower_sum = 0.0;
    for (std::size_t index = 0; index < dimension; ++index) {
        if (coordinates.charges[coordinates.independent_indices[index]] != 0.0) {
            charged_lower_sum += coordinates.independent_lower_bounds[index];
        }
    }
    if (!std::isnan(total_ion_mole_fraction_max)
        && (!std::isfinite(total_ion_mole_fraction_max)
            || total_ion_mole_fraction_max <= charged_lower_sum
            || total_ion_mole_fraction_max > 1.0)) {
        throw std::invalid_argument(
            "HELD2 Provider total-ion mole-fraction ceiling is invalid or infeasible"
        );
    }
    double remaining_ion = std::isnan(total_ion_mole_fraction_max)
        ? std::numeric_limits<double>::infinity()
        : total_ion_mole_fraction_max - charged_lower_sum;
    for (std::size_t index = 0; index < dimension; ++index) {
        if (coordinates.independent_lower_bounds[index] < 0.0
            || coordinates.independent_upper_bounds[index]
                < coordinates.independent_lower_bounds[index]) {
            throw std::invalid_argument("HELD2 simplex chart bounds are invalid");
        }
        if (unit_cube_coordinates[index] < 0.0
            || unit_cube_coordinates[index] > 1.0) {
            throw std::invalid_argument("HELD2 simplex cube coordinate is outside [0, 1]");
        }
    }

    std::vector<double> independent(dimension, 0.0);
    for (std::size_t index = 0; index < dimension; ++index) {
        const bool charged =
            coordinates.charges[coordinates.independent_indices[index]] != 0.0;
        double available = charged
            ? std::min(remaining_composition, remaining_ion)
            : remaining_composition;
        available = std::min(
            available,
            coordinates.independent_upper_bounds[index]
                - coordinates.independent_lower_bounds[index]
        );
        const double shifted = available * unit_cube_coordinates[index];
        independent[index] = coordinates.independent_lower_bounds[index] + shifted;
        remaining_composition -= shifted;
        if (charged) {
            remaining_ion -= shifted;
        }
    }
    return independent;
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
    const std::size_t component_count = coordinates.charges.size();
    const std::size_t independent_count = coordinates.independent_indices.size();
    const std::size_t provider_coordinate_count = component_count + 1;
    const std::size_t reduced_coordinate_count = independent_count + 1;
    if (independent_modified_fractions.size() != independent_count) {
        throw std::invalid_argument(
            "independent modified composition size does not match the HELD2 chart"
        );
    }
    require_finite_vector(independent_modified_fractions, "independent modified fractions");
    if (!std::isfinite(log_volume) || !std::isfinite(pressure_over_rt)
        || !std::isfinite(target_pressure_pa) || pressure_over_rt <= 0.0
        || target_pressure_pa <= 0.0) {
        throw std::invalid_argument("HELD2 phase state and pressure scales must be finite and positive");
    }
    if (block.gradient.size() != provider_coordinate_count
        || block.hessian.size() != provider_coordinate_count * provider_coordinate_count) {
        throw std::invalid_argument("HELD2 physical phase block has the wrong tensor dimensions");
    }
    require_finite_vector(block.gradient, "HELD2 physical phase gradient");
    require_finite_vector(block.hessian, "HELD2 physical phase Hessian");
    if (!std::isfinite(block.helmholtz_over_rt) || !std::isfinite(block.pressure_pa)) {
        throw std::invalid_argument("HELD2 physical phase block scalars must be finite");
    }

    Held2StateEvaluation result;
    result.modified_fractions.resize(coordinates.retained_indices.size(), 0.0);
    for (std::size_t independent = 0; independent < independent_count; ++independent) {
        const double value = independent_modified_fractions[independent];
        const auto retained = static_cast<std::size_t>(
            std::find(
                coordinates.retained_indices.begin(),
                coordinates.retained_indices.end(),
                coordinates.independent_indices[independent]
            ) - coordinates.retained_indices.begin()
        );
        result.modified_fractions[retained] = value;
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
    result.modified_fractions[dependent_retained] =
        1.0 - std::accumulate(
            independent_modified_fractions.begin(),
            independent_modified_fractions.end(),
            0.0
        );
    result.volume = std::exp(log_volume);
    if (!std::isfinite(result.volume) || result.volume <= 0.0) {
        throw std::invalid_argument("HELD2 phase volume must be finite and positive");
    }

    std::vector<double> jacobian(
        provider_coordinate_count * reduced_coordinate_count,
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
            jacobian[component * reduced_coordinate_count + independent] =
                modified_derivative / coordinates.modified_factors[retained];
        }
        double eliminated_derivative = 0.0;
        const double eliminated_charge = coordinates.charges[coordinates.eliminated_index];
        for (std::size_t component : coordinates.retained_indices) {
            eliminated_derivative -= coordinates.charges[component] / eliminated_charge
                * jacobian[component * reduced_coordinate_count + independent];
        }
        jacobian[coordinates.eliminated_index * reduced_coordinate_count + independent] =
            eliminated_derivative;
    }
    jacobian[component_count * reduced_coordinate_count + independent_count] = result.volume;

    std::vector<double> augmented_gradient = block.gradient;
    augmented_gradient.back() += pressure_over_rt;
    result.objective = block.helmholtz_over_rt + pressure_over_rt * result.volume;
    result.gradient.assign(reduced_coordinate_count, 0.0);
    for (std::size_t reduced = 0; reduced < reduced_coordinate_count; ++reduced) {
        for (std::size_t provider = 0; provider < provider_coordinate_count; ++provider) {
            result.gradient[reduced] +=
                jacobian[provider * reduced_coordinate_count + reduced]
                * augmented_gradient[provider];
        }
    }
    result.hessian.assign(reduced_coordinate_count * reduced_coordinate_count, 0.0);
    for (std::size_t row = 0; row < reduced_coordinate_count; ++row) {
        for (std::size_t column = 0; column < reduced_coordinate_count; ++column) {
            double value = 0.0;
            for (std::size_t left = 0; left < provider_coordinate_count; ++left) {
                for (std::size_t right = 0; right < provider_coordinate_count; ++right) {
                    value += jacobian[left * reduced_coordinate_count + row]
                        * block.hessian[left * provider_coordinate_count + right]
                        * jacobian[right * reduced_coordinate_count + column];
                }
            }
            result.hessian[row * reduced_coordinate_count + column] = value;
        }
    }
    result.hessian.back() += result.volume * augmented_gradient.back();
    std::vector<double> physical_potentials(block.gradient.begin(), block.gradient.end() - 1);
    result.modified_potentials = held2_transform_modified_potentials(
        coordinates,
        physical_potentials
    );
    result.pressure_stationarity_relative =
        (block.pressure_pa - target_pressure_pa) / target_pressure_pa;
    result.pressure_stationarity_derivative_log_volume =
        (result.gradient.back() - result.hessian.back())
        / (pressure_over_rt * result.volume);
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
    const double feed_independent = modified_feed[independent_retained];
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

    const std::vector<double> reference_initial = {feed_independent, 0.0};
    const Held2StateEvaluation initial_reference = evaluator({feed_independent}, 0.0);
    const Held2SearchRun reference_run = solve_held2_search(
        evaluator,
        reference_initial,
        initial_reference,
        false,
        reference_initial,
        {feed_independent, std::log(0.5)},
        {feed_independent, std::log(1.5)}
    );
    Held2StageIResult result;
    result.declared_start_count = 10 * static_cast<int>(charges.size());
    result.minimum_tpd = std::numeric_limits<double>::infinity();
    if (!reference_run.solver_converged || !reference_run.callback_error.empty()
        || reference_run.variables.size() != 2) {
        result.outcome = "indeterminate";
        return result;
    }
    const Held2StateEvaluation reference = evaluator(
        {reference_run.variables.front()},
        reference_run.variables.back()
    );
    result.reference_modified_fractions = reference.modified_fractions;
    result.reference_volume = reference.volume;

    const double lower_composition = coordinates.independent_lower_bounds.front();
    const double upper_composition = std::nextafter(
        coordinates.independent_upper_bounds.front(),
        lower_composition
    );
    const std::vector<double> lower = {lower_composition, std::log(0.5)};
    const std::vector<double> upper = {upper_composition, std::log(1.5)};
    std::mt19937 generator(kStageISeed);
    std::uniform_real_distribution<double> composition_distribution(
        lower_composition,
        upper_composition
    );
    std::uniform_real_distribution<double> log_volume_distribution(
        lower.back(),
        upper.back()
    );
    std::vector<std::vector<double>> starts;
    starts.reserve(static_cast<std::size_t>(result.declared_start_count));
    starts.push_back({0.2, 0.0});
    starts.push_back({0.8, 0.0});
    while (starts.size() < static_cast<std::size_t>(result.declared_start_count)) {
        starts.push_back({
            composition_distribution(generator),
            log_volume_distribution(generator),
        });
    }
    const std::vector<double> reference_variables = {
        reference_run.variables.front(),
        reference_run.variables.back(),
    };
    for (const std::vector<double>& start : starts) {
        const Held2SearchRun run = solve_held2_search(
            evaluator,
            reference_variables,
            reference,
            true,
            start,
            lower,
            upper
        );
        if (!run.solver_converged || !run.callback_error.empty() || run.variables.size() != 2) {
            ++result.failed_start_count;
            continue;
        }
        ++result.completed_start_count;
        Held2StateEvaluation state = evaluator({run.variables.front()}, run.variables.back());
        double tpd = state.objective - reference.objective;
        for (std::size_t index = 0; index < state.gradient.size(); ++index) {
            tpd -= reference.gradient[index]
                * (run.variables[index] - reference_variables[index]);
        }
        result.minimum_tpd = std::min(result.minimum_tpd, tpd);
        if (!audit_held2_tolerance(kHeld2TpdNegativeMargin, tpd).passed) {
            continue;
        }
        const bool duplicate = std::any_of(
            result.candidates.begin(),
            result.candidates.end(),
            [&state](const Held2StageICandidate& candidate) {
                return maximum_abs_difference(
                           candidate.modified_fractions,
                           state.modified_fractions
                       ) <= kHeld2BasinDuplicateComposition.atol
                    && std::abs(
                           std::log(candidate.volume) - std::log(state.volume)
                       ) <= kHeld2BasinDuplicateLogVolume.atol;
            }
        );
        if (!duplicate) {
            result.candidates.push_back({state.modified_fractions, state.volume, tpd});
        }
    }
    std::sort(
        result.candidates.begin(),
        result.candidates.end(),
        [](const Held2StageICandidate& left, const Held2StageICandidate& right) {
            return left.modified_fractions < right.modified_fractions;
        }
    );
    if (result.failed_start_count != 0) {
        result.outcome = "indeterminate";
    } else if (!result.candidates.empty()) {
        result.outcome = "negative_tpd";
    } else {
        result.outcome = "no_negative_found";
    }
    return result;
}

Held2StageIIResult solve_held2_manufactured_stage_ii_impl(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed
) {
    constexpr int major_iteration_cap = 100;
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
    result.search_strategy = "continuation_sobol_direct_l_ipopt_v1";
    result.global_explorer = "continuation_sobol_direct_l";
    result.lower_starts_per_iteration = 0;
    std::vector<Held2StateEvaluation> traced_basins;
    bool request_direct_escalation = false;

    for (int major = 0; major < major_iteration_cap; ++major) {
        Held2StageIIUpperProblem upper_problem;
        upper_problem.multiplier_lower_bounds = {
            -std::numeric_limits<double>::infinity(),
        };
        upper_problem.multiplier_upper_bounds = {
            std::numeric_limits<double>::infinity(),
        };
        upper_problem.cuts.reserve(cuts.size() + 1);
        for (std::size_t cut_index = 0; cut_index < cuts.size(); ++cut_index) {
            const Held2StateEvaluation& cut = cuts[cut_index];
            upper_problem.cuts.push_back({
                static_cast<int>(cut_index),
                cut.objective,
                {feed - cut.modified_fractions[independent_retained]},
            });
        }
        upper_problem.cuts.push_back({
            static_cast<int>(cuts.size()),
            feed_gibbs,
            {0.0},
        });
        const Held2StageIIUpperResult upper_solve =
            solve_held2_stage_ii_upper_highs(upper_problem);
        if (upper_solve.outcome != "optimal"
            || upper_solve.multipliers.size() != 1) {
            result.outcome = "indeterminate";
            result.cut_count = static_cast<int>(cuts.size());
            return result;
        }
        const double upper_bound = upper_solve.upper_bound;
        const double multiplier = upper_solve.multipliers.front();

        Held2StateEvaluation lower_reference;
        lower_reference.gradient = {multiplier, 0.0};
        lower_reference.hessian.assign(4, 0.0);
        const std::vector<double> lower_reference_variables = {feed, 0.0};
        std::vector<std::vector<double>> starts;
        std::vector<std::string> start_sources;
        {
            std::vector<Held2StageIIBasinSeed> seeds;
            for (const Held2StateEvaluation& basin : traced_basins) {
                seeds.push_back({
                    {basin.modified_fractions[independent_retained]},
                    "continuation",
                });
            }
            for (const Held2StateEvaluation& cut : cuts) {
                seeds.push_back({
                    {cut.modified_fractions[independent_retained]},
                    "cut_state",
                });
            }
            seeds.push_back({{0.2}, "stage_i_witness"});
            seeds.push_back({{0.8}, "stage_i_witness"});
            seeds.push_back({{feed}, "homogeneous_reference"});
            seeds.push_back({
                {lower_composition + 0.05 * (upper_composition - lower_composition)},
                "boundary_aware_seed",
            });
            seeds.push_back({
                {upper_composition - 0.05 * (upper_composition - lower_composition)},
                "boundary_aware_seed",
            });
            const Held2StageIIBasinEvaluator basin_evaluator = [
                &coordinates,
                &evaluator,
                multiplier,
                feed,
                independent_retained
            ](const std::vector<double>& independent) {
                Held2StageIIBasinEvaluation evaluation;
                evaluation.independent_modified_fractions = independent;
                const Held2StateEvaluator pressure_evaluator = [&evaluator](
                    const std::vector<double>& composition,
                    double log_volume
                ) {
                    Held2StateEvaluation state = evaluator(composition, log_volume);
                    state.pressure_stationarity_relative *= -1.0;
                    state.pressure_stationarity_derivative_log_volume *= -1.0;
                    return state;
                };
                evaluation.pressure_envelope = evaluate_held2_pressure_envelope(
                    independent,
                    {0.5, 1.5},
                    pressure_evaluator,
                    64,
                    8
                );
                if (evaluation.pressure_envelope.outcome != "selected"
                    && evaluation.pressure_envelope.failure_reason
                        != "stable_objective_tie") {
                    evaluation.failure_reason =
                        evaluation.pressure_envelope.failure_reason;
                    return evaluation;
                }
                for (const Held2PressureRoot& root :
                     evaluation.pressure_envelope.roots) {
                    if (root.mechanical_class == "strict_stable"
                        && !root.boundary) {
                        evaluation.reduced_lower_value = std::min(
                            evaluation.reduced_lower_value,
                            root.objective + multiplier * (
                                feed
                                - root.state.modified_fractions[
                                    independent_retained
                                ]
                            )
                        );
                    }
                }
                if (!std::isfinite(evaluation.reduced_lower_value)) {
                    evaluation.failure_reason = "no_strict_stable_root";
                    return evaluation;
                }
                evaluation.certified = true;
                return evaluation;
            };
            const Held2StageIIBasinExplorationResult exploration =
                explore_held2_stage_ii_basins(
                    coordinates,
                    seeds,
                    12,
                    request_direct_escalation,
                    24,
                    std::numeric_limits<double>::quiet_NaN(),
                    basin_evaluator
                );
            result.exploration_evaluation_count +=
                exploration.completed_evaluation_count;
            result.exploration_failure_count +=
                exploration.failed_evaluation_count;
            result.exploration_representative_count += static_cast<int>(
                exploration.representatives.size()
            );
            result.duplicate_representative_count +=
                exploration.duplicate_start_count;
            result.direct_escalation_used = result.direct_escalation_used
                || exploration.direct_escalation_used;
            if (exploration.outcome != "representatives_found") {
                result.outcome = "indeterminate_exploration_failure";
                result.cut_count = static_cast<int>(cuts.size());
                return result;
            }
            for (const Held2StageIIPhysicalStart& representative :
                 exploration.representatives) {
                starts.push_back({
                    representative.independent_modified_fractions.front(),
                    representative.log_volume,
                });
                start_sources.push_back(representative.source);
            }
            result.lower_starts_per_iteration = std::max(
                result.lower_starts_per_iteration,
                static_cast<int>(starts.size())
            );
        }
        double lower_bound = std::numeric_limits<double>::infinity();
        Held2StateEvaluation best_state;
        bool lower_failed = false;
        for (std::size_t start_index = 0; start_index < starts.size(); ++start_index) {
            const std::vector<double>& start = starts[start_index];
            Held2StageIIAttempt attempt;
            attempt.attempt_id = static_cast<int>(result.attempt_trace.size());
            attempt.major_iteration = major;
            attempt.start_index = static_cast<int>(start_index);
            attempt.start_source = start_sources[start_index];
            attempt.provider_status = "manufactured_oracle";
            attempt.internal_start = start;
            attempt.same_major_upper_bound = upper_bound;
            attempt.same_major_multipliers = {multiplier};
            const Held2StateEvaluation start_state = evaluator(
                {start.front()},
                start.back()
            );
            attempt.physical_start_modified_fractions =
                start_state.modified_fractions;
            attempt.physical_start_volume = start_state.volume;

            const Held2SearchRun run = solve_held2_search(
                evaluator,
                lower_reference_variables,
                lower_reference,
                true,
                start,
                lower,
                upper
            );
            attempt.solver_status = run.solver_status;
            attempt.solver_converged = run.solver_converged;
            attempt.callback_error = run.callback_error;
            attempt.internal_terminal = run.variables;
            attempt.lower_bound_multipliers = run.lower_bound_multipliers;
            attempt.upper_bound_multipliers = run.upper_bound_multipliers;
            if (!run.solver_converged || !run.callback_error.empty()
                || run.variables.size() != 2) {
                lower_failed = true;
                result.attempt_trace.push_back(std::move(attempt));
                continue;
            }
            Held2StateEvaluation state = evaluator(
                {run.variables.front()},
                run.variables.back()
            );
            const double value = state.objective
                + multiplier * (feed - state.modified_fractions[independent_retained]);
            attempt.terminal_modified_fractions = state.modified_fractions;
            attempt.terminal_volume = state.volume;
            attempt.objective = state.objective;
            attempt.lower_value = value;
            attempt.pressure_residual = state.pressure_stationarity_relative;
            attempt.pressure_passed = audit_held2_tolerance(
                kHeld2RootPressure, state.pressure_stationarity_relative
            ).passed;

            if (run.lower_bound_multipliers.size() == run.variables.size()
                && run.upper_bound_multipliers.size() == run.variables.size()) {
                const std::array<double, 2> lower_gradient = {
                    state.gradient.front() - multiplier,
                    state.gradient.back(),
                };
                double stationarity_inf = 0.0;
                double complementarity_inf = 0.0;
                double dual_sign_violation = 0.0;
                bool dual_signs_valid = true;
                for (std::size_t index = 0; index < run.variables.size(); ++index) {
                    dual_sign_violation = std::max({
                        dual_sign_violation,
                        -run.lower_bound_multipliers[index],
                        -run.upper_bound_multipliers[index],
                    });
                    dual_signs_valid = dual_signs_valid
                        && audit_held2_tolerance(
                            kHeld2Stage2DualSign,
                            std::max({
                                0.0,
                                -run.lower_bound_multipliers[index],
                                -run.upper_bound_multipliers[index],
                            })
                        ).passed;
                    stationarity_inf = std::max(
                        stationarity_inf,
                        std::abs(
                            lower_gradient[index]
                            - run.lower_bound_multipliers[index]
                            + run.upper_bound_multipliers[index]
                        )
                    );
                    complementarity_inf = std::max(
                        complementarity_inf,
                        std::abs(
                            (run.variables[index] - lower[index])
                            * run.lower_bound_multipliers[index]
                        )
                    );
                    complementarity_inf = std::max(
                        complementarity_inf,
                        std::abs(
                            (upper[index] - run.variables[index])
                            * run.upper_bound_multipliers[index]
                        )
                    );
                }
                attempt.chart_kkt_inf_norm = stationarity_inf;
                // The manufactured one-dimensional composition chart is the
                // physical coordinate itself, so the pullback is the identity.
                attempt.physical_kkt_inf_norm = stationarity_inf;
                attempt.complementarity_inf_norm = complementarity_inf;
                attempt.primal_inf_norm = 0.0;
                attempt.dual_sign_violation_inf_norm = dual_sign_violation;
                attempt.dual_pullback_inf_norm = 0.0;
                attempt.dual_pullback_scale = std::max({
                    std::abs(run.lower_bound_multipliers.front()),
                    std::abs(run.upper_bound_multipliers.front()),
                    std::abs(run.lower_bound_multipliers.back()),
                    std::abs(run.upper_bound_multipliers.back()),
                });
                attempt.dual_signs_valid = dual_signs_valid;
                attempt.physical_kkt_passed = attempt.pressure_passed
                    && dual_signs_valid
                    && audit_held2_tolerance(
                           kHeld2Stage2Stationarity, stationarity_inf
                       ).passed
                    && audit_held2_tolerance(
                           kHeld2Stage2Complementarity, complementarity_inf
                       ).passed;
            }
            attempt.cut_eligible = attempt.physical_kkt_passed;
            const double lower_gap = upper_bound - value;
            attempt.step6_gap = lower_gap;
            attempt.fixed_volume_gradient_inf_norm = std::abs(
                state.gradient.front() - multiplier
            );
            attempt.fixed_volume_gradient_scale = std::max(
                std::abs(state.gradient.front()), std::abs(multiplier)
            );
            attempt.step6_eligible = attempt.cut_eligible
                && audit_held2_tolerance(kHeld2Step6Gap, lower_gap).passed
                && audit_held2_tolerance(
                       kHeld2Step6Gradient,
                       attempt.fixed_volume_gradient_inf_norm,
                       attempt.fixed_volume_gradient_scale
                   ).passed;

            const auto basin = std::find_if(
                traced_basins.begin(),
                traced_basins.end(),
                [&state](const Held2StateEvaluation& known) {
                    return maximum_abs_difference(
                               known.modified_fractions,
                               state.modified_fractions
                           ) <= kHeld2BasinDuplicateComposition.atol
                        && std::abs(
                               std::log(known.volume) - std::log(state.volume)
                           ) <= kHeld2BasinDuplicateLogVolume.atol;
                }
            );
            if (basin == traced_basins.end()) {
                attempt.basin_id = static_cast<int>(traced_basins.size());
                traced_basins.push_back(state);
            } else {
                ++result.duplicate_terminal_count;
                attempt.basin_id = static_cast<int>(
                    std::distance(traced_basins.begin(), basin)
                );
            }
            result.attempt_trace.push_back(std::move(attempt));
            if (value < lower_bound) {
                lower_bound = value;
                best_state = state;
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
            true,
            upper_bound,
            {multiplier},
            static_cast<int>(cuts.size()),
            upper_solve.solver,
            upper_solve.solver_version,
            upper_solve.solver_status,
            upper_solve.primal_feasible,
            upper_solve.dual_feasible,
            upper_solve.primal_residual_inf,
            upper_solve.primal_scale,
            upper_solve.dual_residual_inf,
            upper_solve.dual_scale,
            upper_solve.complementarity_inf,
            upper_solve.cut_slacks,
            upper_solve.cut_duals,
            upper_solve.active_cut_ids,
        });
        const bool duplicate = std::any_of(
            cuts.begin(),
            cuts.end(),
            [&best_state](const Held2StateEvaluation& cut) {
                return maximum_abs_difference(
                           cut.modified_fractions,
                           best_state.modified_fractions
                       ) <= kHeld2BasinDuplicateComposition.atol
                    && std::abs(
                           std::log(cut.volume) - std::log(best_state.volume)
                       ) <= kHeld2BasinDuplicateLogVolume.atol;
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
            if (audit_held2_tolerance(kHeld2Step6Gap, gap).passed
                && audit_held2_tolerance(
                       kHeld2Step6Gradient,
                       cut.gradient.front() - multiplier,
                       std::max(std::abs(cut.gradient.front()), std::abs(multiplier))
                   ).passed
                && audit_held2_tolerance(
                       kHeld2RootPressure, cut.gradient.back()
                   ).passed) {
                bool duplicate_candidate = false;
                bool unresolved_identity = false;
                for (const Held2StageIICandidate& known : result.candidates) {
                    const double composition_distance = maximum_abs_difference(
                        known.modified_fractions, cut.modified_fractions
                    );
                    const double log_volume_distance = std::abs(
                        std::log(known.volume) - std::log(cut.volume)
                    );
                    if (composition_distance
                            <= kHeld2BasinDuplicateComposition.atol
                        && log_volume_distance
                            <= kHeld2BasinDuplicateLogVolume.atol) {
                        duplicate_candidate = true;
                        break;
                    }
                    if (!audit_held2_tolerance(
                            kHeld2CandidateDistinctComposition,
                            composition_distance
                        ).passed
                        && !audit_held2_tolerance(
                            kHeld2CandidateDistinctLogVolume,
                            log_volume_distance
                        ).passed) {
                        unresolved_identity = true;
                        break;
                    }
                }
                if (unresolved_identity) {
                    ++result.unresolved_candidate_identity_count;
                } else if (!duplicate_candidate) {
                    result.candidates.push_back({
                        cut.modified_fractions,
                        {cut.modified_fractions[independent_retained]},
                        cut.volume,
                        std::log(cut.volume),
                        gap,
                    });
                }
            }
        }
        if (result.unresolved_candidate_identity_count != 0) {
            result.outcome = "indeterminate_candidate_identity";
            result.cut_count = static_cast<int>(cuts.size());
            result.distinct_basin_count = static_cast<int>(traced_basins.size());
            return result;
        }
        if (result.candidates.size() > 1) {
            result.outcome = "candidate_set";
            result.cut_count = static_cast<int>(cuts.size());
            result.distinct_basin_count = static_cast<int>(traced_basins.size());
            return result;
        }
        if (duplicate) {
            if (!request_direct_escalation) {
                request_direct_escalation = true;
                continue;
            }
            result.outcome = "indeterminate_finite_search_stalled";
            result.cut_count = static_cast<int>(cuts.size());
            result.distinct_basin_count = static_cast<int>(traced_basins.size());
            return result;
        }
    }
    result.outcome = "resource_limit";
    result.cut_count = static_cast<int>(cuts.size());
    result.distinct_basin_count = static_cast<int>(traced_basins.size());
    return result;
}

Held2StageIIResult solve_held2_manufactured_stage_ii(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed
) {
    return solve_held2_manufactured_stage_ii_impl(
        charges,
        physical_feed
    );
}

Held2StageIIResult solve_held2_stage_ii(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const Held2StateEvaluator& evaluator,
    const Held2VolumeBoundsEvaluator& volume_bounds_evaluator,
    const Held2StateEvaluation& reference,
    const std::vector<Held2StageICandidate>& stage_i_candidates,
    double total_ion_mole_fraction_max,
    int major_iteration_cap,
    int local_attempt_cap_per_major,
    Held2ProgressObserver* observer
) {
    if (major_iteration_cap <= 0 || local_attempt_cap_per_major <= 0) {
        throw std::invalid_argument(
            "HELD2 Stage II resource caps must be positive"
        );
    }
    const std::size_t dimension = coordinates.independent_indices.size();
    if (dimension == 0 || reference.gradient.size() != dimension + 1
        || !(reference.volume > 0.0)) {
        throw std::invalid_argument("HELD2 Stage II reference evidence is incomplete");
    }
    const std::vector<double> modified_feed = held2_transform_physical_fractions(
        coordinates,
        physical_feed
    );
    std::vector<std::size_t> retained_positions;
    std::vector<double> feed;
    for (std::size_t component : coordinates.independent_indices) {
        const auto found = std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            component
        );
        if (found == coordinates.retained_indices.end()) {
            throw std::invalid_argument(
                "HELD2 Stage II independent coordinate is not retained"
            );
        }
        retained_positions.push_back(static_cast<std::size_t>(
            found - coordinates.retained_indices.begin()
        ));
        feed.push_back(modified_feed[retained_positions.back()]);
    }
    const auto dependent = std::find(
        coordinates.retained_indices.begin(),
        coordinates.retained_indices.end(),
        coordinates.dependent_index
    );
    if (dependent == coordinates.retained_indices.end()) {
        throw std::invalid_argument("HELD2 Stage II dependent coordinate is not retained");
    }
    const std::size_t dependent_position = static_cast<std::size_t>(
        dependent - coordinates.retained_indices.begin()
    );
    const double composition_sum_upper = 1.0
        - kModifiedLowerScale * coordinates.modified_factors[dependent_position];
    std::vector<double> physical_lower = coordinates.independent_lower_bounds;
    std::vector<double> physical_upper = coordinates.independent_upper_bounds;
    for (std::size_t index = 0; index < dimension; ++index) {
        physical_upper[index] = std::nextafter(
            physical_upper[index],
            physical_lower[index]
        );
    }
    std::vector<bool> source_domain_charged;
    double local_ion_mole_fraction_max =
        std::numeric_limits<double>::quiet_NaN();
    if (!std::isnan(total_ion_mole_fraction_max)) {
        if (!std::isfinite(total_ion_mole_fraction_max)
            || total_ion_mole_fraction_max <= 0.0
            || total_ion_mole_fraction_max > 1.0) {
            throw std::invalid_argument(
                "HELD2 Stage II Provider ion ceiling is invalid"
            );
        }
        source_domain_charged.assign(dimension, false);
        double charged_lower_sum = 0.0;
        bool has_independent_charge = false;
        for (std::size_t index = 0; index < dimension; ++index) {
            const std::size_t component = coordinates.independent_indices[index];
            source_domain_charged[index] = coordinates.charges[component] != 0.0;
            if (source_domain_charged[index]) {
                // In Perdomo's electroneutral modified coordinates, the sum of
                // retained charged fractions is the explicit total-ion fraction.
                has_independent_charge = true;
                charged_lower_sum += physical_lower[index];
            }
        }
        if (has_independent_charge
            && total_ion_mole_fraction_max < composition_sum_upper) {
            if (!(total_ion_mole_fraction_max > charged_lower_sum)) {
                throw std::invalid_argument(
                    "HELD2 Stage II Provider ion ceiling excludes the feasible interior"
                );
            }
            local_ion_mole_fraction_max = total_ion_mole_fraction_max;
        } else {
            source_domain_charged.clear();
        }
    }
    const auto independent_from_state = [
        &retained_positions
    ](const Held2StateEvaluation& state) {
        std::vector<double> independent;
        independent.reserve(retained_positions.size());
        for (std::size_t retained : retained_positions) {
            independent.push_back(state.modified_fractions[retained]);
        }
        return independent;
    };
    struct Cut {
        Held2StateEvaluation state;
        std::vector<double> independent;
        std::vector<double> fixed_volume_gradient;
        bool candidate_phase_eligible = false;
    };
    const auto make_cut = [&evaluator, dimension](
        const std::vector<double>& independent,
        double volume,
        bool candidate_phase_eligible = false
    ) {
        if (!(volume > 0.0) || !std::isfinite(volume)) {
            throw std::invalid_argument("HELD2 Stage II cut volume is invalid");
        }
        Held2StateEvaluation state = evaluator(independent, std::log(volume));
        if (state.gradient.size() != dimension + 1) {
            throw std::invalid_argument(
                "HELD2 Stage II cut gradient is incomplete"
            );
        }
        return Cut{
            state,
            independent,
            std::vector<double>(state.gradient.begin(), state.gradient.end() - 1),
            candidate_phase_eligible,
        };
    };
    const auto make_pressure_cut = [
        &evaluator,
        &volume_bounds_evaluator,
        &make_cut
    ](const std::vector<double>& independent) {
        const Held2PressureEnvelopeResult envelope =
            evaluate_held2_pressure_envelope(
                independent,
                volume_bounds_evaluator(independent),
                evaluator,
                64,
                8
            );
        if (envelope.outcome != "selected" || envelope.selected_root_index < 0) {
            throw std::runtime_error(
                "HELD2 Stage II initial cut pressure envelope failed: "
                + envelope.failure_reason
            );
        }
        return make_cut(
            independent,
            envelope.roots[static_cast<std::size_t>(
                envelope.selected_root_index
            )].volume
        );
    };

    std::vector<Cut> cuts;
    cuts.push_back(make_cut(feed, reference.volume));
    for (const Held2StageICandidate& candidate : stage_i_candidates) {
        if (candidate.modified_fractions.size() != modified_feed.size()) {
            continue;
        }
        Held2StateEvaluation state;
        state.modified_fractions = candidate.modified_fractions;
        cuts.push_back(make_cut(independent_from_state(state), candidate.volume));
    }
    for (std::size_t coordinate = 0;
         cuts.size() < dimension + 1 && coordinate < dimension;
         ++coordinate) {
        std::vector<double> point = feed;
        const double lower_sum = std::accumulate(
            physical_lower.begin(), physical_lower.end(), 0.0
        );
        const double free = composition_sum_upper - lower_sum;
        const double delta = std::min(
            0.05 * free,
            0.5 * (composition_sum_upper
                - std::accumulate(point.begin(), point.end(), 0.0))
        );
        double domain_safe_delta = delta;
        if (!source_domain_charged.empty()
            && source_domain_charged[coordinate]) {
            double charged_sum = 0.0;
            for (std::size_t index = 0; index < dimension; ++index) {
                if (source_domain_charged[index]) {
                    charged_sum += point[index];
                }
            }
            domain_safe_delta = std::min(
                domain_safe_delta,
                0.5 * (local_ion_mole_fraction_max - charged_sum)
            );
        }
        if (!(domain_safe_delta > 0.0)) {
            continue;
        }
        point[coordinate] += domain_safe_delta;
        cuts.push_back(make_pressure_cut(point));
    }

    Held2StageIIResult result;
    result.search_strategy = "continuation_sobol_direct_l_ipopt_v1";
    result.global_explorer = "continuation_sobol_direct_l";
    result.local_attempt_cap_per_major = local_attempt_cap_per_major;
    if (cuts.size() < dimension + 1) {
        result.outcome = "indeterminate_initial_cut_set";
        result.cut_count = static_cast<int>(cuts.size());
        return result;
    }
    const double feed_gibbs = reference.objective;
    std::vector<Held2StateEvaluation> traced_basins;
    bool request_direct_escalation = false;

    for (int major = 0; major < major_iteration_cap; ++major) {
        Held2StageIIUpperProblem upper_problem;
        upper_problem.multiplier_lower_bounds.assign(
            dimension, -std::numeric_limits<double>::infinity()
        );
        upper_problem.multiplier_upper_bounds.assign(
            dimension, std::numeric_limits<double>::infinity()
        );
        for (std::size_t cut_index = 0; cut_index < cuts.size(); ++cut_index) {
            std::vector<double> slopes(dimension, 0.0);
            for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                slopes[coordinate] = feed[coordinate]
                    - cuts[cut_index].independent[coordinate];
            }
            upper_problem.cuts.push_back({
                static_cast<int>(cut_index),
                cuts[cut_index].state.objective,
                std::move(slopes),
            });
        }
        upper_problem.cuts.push_back({
            static_cast<int>(cuts.size()),
            feed_gibbs,
            std::vector<double>(dimension, 0.0),
        });
        const Held2StageIIUpperResult upper_solve =
            solve_held2_stage_ii_upper_highs(upper_problem);
        Held2ProgressEvent upper_progress;
        upper_progress.kind = Held2ProgressKind::StageIIUpper;
        upper_progress.major_iteration = major;
        upper_progress.upper_bound = upper_solve.upper_bound;
        upper_progress.primal_residual = upper_solve.primal_residual_inf;
        upper_progress.dual_residual = upper_solve.dual_residual_inf;
        upper_progress.complementarity = upper_solve.complementarity_inf;
        upper_progress.count = static_cast<int>(
            upper_solve.active_cut_ids.size()
        );
        upper_progress.status = upper_solve.outcome;
        observe_held2(observer, upper_progress);
        if (upper_solve.outcome != "optimal"
            || upper_solve.multipliers.size() != dimension
            || !upper_solve.primal_feasible || !upper_solve.dual_feasible) {
            result.outcome = "indeterminate_upper_lp";
            result.cut_count = static_cast<int>(cuts.size());
            return result;
        }
        const double upper_bound = upper_solve.upper_bound;
        const std::vector<double>& multipliers = upper_solve.multipliers;
        result.bound_history.push_back({
            0.0,
            false,
            upper_bound,
            multipliers,
            static_cast<int>(cuts.size()),
            upper_solve.solver,
            upper_solve.solver_version,
            upper_solve.solver_status,
            upper_solve.primal_feasible,
            upper_solve.dual_feasible,
            upper_solve.primal_residual_inf,
            upper_solve.primal_scale,
            upper_solve.dual_residual_inf,
            upper_solve.dual_scale,
            upper_solve.complementarity_inf,
            upper_solve.cut_slacks,
            upper_solve.cut_duals,
            upper_solve.active_cut_ids,
        });

        std::vector<Held2StageIIBasinSeed> seeds;
        for (const Held2StateEvaluation& basin : traced_basins) {
            seeds.push_back({independent_from_state(basin), "continuation"});
        }
        for (const Cut& cut : cuts) {
            seeds.push_back({cut.independent, "cut_state"});
        }
        for (const Held2StageICandidate& candidate : stage_i_candidates) {
            if (candidate.modified_fractions.size() == modified_feed.size()) {
                Held2StateEvaluation state;
                state.modified_fractions = candidate.modified_fractions;
                seeds.push_back({independent_from_state(state), "stage_i_witness"});
            }
        }
        seeds.push_back({feed, "homogeneous_reference"});
        for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
            std::vector<double> low = feed;
            low[coordinate] = physical_lower[coordinate]
                + 0.05 * (feed[coordinate] - physical_lower[coordinate]);
            seeds.push_back({std::move(low), "boundary_aware_seed"});
        }
        for (std::size_t vertex = 0; vertex < dimension; ++vertex) {
            std::vector<double> cube(dimension, 0.05);
            cube[vertex] = 0.95;
            seeds.push_back({
                held2_map_unit_cube_to_independent_fractions(
                    coordinates,
                    cube,
                    total_ion_mole_fraction_max
                ),
                "simplex_vertex_seed",
            });
        }
        const Held2StageIIBasinEvaluator basin_evaluator = [
            &coordinates,
            &evaluator,
            &volume_bounds_evaluator,
            &feed,
            &multipliers,
            &retained_positions
        ](const std::vector<double>& independent) {
            Held2StageIIBasinEvaluation evaluation;
            evaluation.independent_modified_fractions = independent;
            const std::array<double, 2> bounds = volume_bounds_evaluator(independent);
            evaluation.pressure_envelope = evaluate_held2_pressure_envelope(
                independent,
                bounds,
                evaluator,
                64,
                8
            );
            if (evaluation.pressure_envelope.outcome != "selected"
                && evaluation.pressure_envelope.failure_reason
                    != "stable_objective_tie") {
                evaluation.failure_reason =
                    evaluation.pressure_envelope.failure_reason;
                return evaluation;
            }
            for (const Held2PressureRoot& root : evaluation.pressure_envelope.roots) {
                if (root.mechanical_class != "strict_stable" || root.boundary) {
                    continue;
                }
                double value = root.objective;
                for (std::size_t coordinate = 0;
                     coordinate < multipliers.size(); ++coordinate) {
                    value += multipliers[coordinate] * (
                        feed[coordinate]
                        - root.state.modified_fractions[
                            retained_positions[coordinate]
                        ]
                    );
                }
                evaluation.reduced_lower_value = std::min(
                    evaluation.reduced_lower_value,
                    value
                );
            }
            if (!std::isfinite(evaluation.reduced_lower_value)) {
                evaluation.failure_reason = "no_strict_stable_root";
                return evaluation;
            }
            evaluation.certified = true;
            return evaluation;
        };
        const Held2StageIIBasinExplorationResult exploration =
            explore_held2_stage_ii_basins(
                coordinates,
                seeds,
                kStageIISobolEvaluationCount,
                request_direct_escalation,
                kStageIIDirectEvaluationBudget,
                total_ion_mole_fraction_max,
                basin_evaluator
            );
        result.exploration_evaluation_count += exploration.completed_evaluation_count;
        result.exploration_failure_count += exploration.failed_evaluation_count;
        result.exploration_representative_count += static_cast<int>(
            exploration.representatives.size()
        );
        result.duplicate_representative_count += exploration.duplicate_start_count;
        result.direct_escalation_used = result.direct_escalation_used
            || exploration.direct_escalation_used;
        if (exploration.outcome != "representatives_found") {
            result.outcome = "indeterminate_exploration_failure";
            result.cut_count = static_cast<int>(cuts.size());
            return result;
        }

        double lower_bound = std::numeric_limits<double>::infinity();
        Cut best;
        const std::size_t attempt_count = std::min(
            exploration.representatives.size(),
            static_cast<std::size_t>(local_attempt_cap_per_major)
        );
        const bool representative_cap_truncated =
            attempt_count < exploration.representatives.size();
        bool step5_condition_satisfied = false;
        const std::size_t attempt_trace_start = result.attempt_trace.size();
        std::vector<Cut> same_major_step6_cuts;
        std::vector<std::size_t> representative_order(
            exploration.representatives.size()
        );
        std::iota(
            representative_order.begin(),
            representative_order.end(),
            std::size_t{0}
        );
        if (representative_order.size() > 1) {
            const Held2StageIIPhysicalStart& first =
                exploration.representatives.front();
            std::stable_sort(
                representative_order.begin() + 1,
                representative_order.end(),
                [&exploration, &first](std::size_t left, std::size_t right) {
                    const Held2StageIIPhysicalStart& left_start =
                        exploration.representatives[left];
                    const Held2StageIIPhysicalStart& right_start =
                        exploration.representatives[right];
                    const bool left_same_branch =
                        left_start.stable_branch_index
                        == first.stable_branch_index;
                    const bool right_same_branch =
                        right_start.stable_branch_index
                        == first.stable_branch_index;
                    if (left_same_branch != right_same_branch) {
                        return left_same_branch;
                    }
                    const double left_distance = maximum_abs_difference(
                        left_start.independent_modified_fractions,
                        first.independent_modified_fractions
                    );
                    const double right_distance = maximum_abs_difference(
                        right_start.independent_modified_fractions,
                        first.independent_modified_fractions
                    );
                    if (left_distance != right_distance) {
                        return left_distance > right_distance;
                    }
                    return left_start.reduced_lower_value
                        < right_start.reduced_lower_value;
                }
            );
        }
        for (std::size_t start_index = 0;
             start_index < attempt_count; ++start_index) {
            const std::size_t representative_index =
                representative_order[start_index];
            const Held2StageIIPhysicalStart& representative =
                exploration.representatives[representative_index];
            Held2StageIIAttempt attempt;
            attempt.attempt_id = static_cast<int>(result.attempt_trace.size());
            attempt.major_iteration = major;
            attempt.start_index = static_cast<int>(representative_index);
            attempt.start_source = representative.source;
            attempt.internal_start = representative.independent_modified_fractions;
            attempt.internal_start.push_back(representative.log_volume);
            attempt.same_major_upper_bound = upper_bound;
            attempt.same_major_multipliers = multipliers;
            try {
                const Held2StateEvaluation start_state = evaluator(
                    representative.independent_modified_fractions,
                    representative.log_volume
                );
                attempt.physical_start_modified_fractions =
                    start_state.modified_fractions;
                attempt.physical_start_volume = start_state.volume;
                attempt.provider_status = "provider_exact";
                const Held2SearchRun run = solve_stage_ii_pressure_root_local(
                    evaluator,
                    volume_bounds_evaluator,
                    feed,
                    multipliers,
                    representative.independent_modified_fractions,
                    representative.stable_branch_index,
                    representative.log_volume,
                    physical_lower,
                    physical_upper,
                    composition_sum_upper,
                    source_domain_charged,
                    local_ion_mole_fraction_max,
                    observer,
                    major,
                    attempt.attempt_id
                );
                attempt.solver_status = run.solver_status;
                attempt.solver_converged = run.solver_converged;
                attempt.callback_error = run.callback_error;
                attempt.internal_terminal = run.variables;
                attempt.lower_bound_multipliers = run.lower_bound_multipliers;
                attempt.upper_bound_multipliers = run.upper_bound_multipliers;
                attempt.evaluation_trace = run.evaluation_trace;
                for (const Held2StageIILocalEvaluation& evaluation : run.evaluation_trace) {
                    if (!evaluation.physical_variables.empty()) {
                        attempt.last_valid_physical_variables =
                            evaluation.physical_variables;
                    }
                }
                attempt.provider_status = "provider_exact";
                if (!run.solver_converged || !run.callback_error.empty()
                    || run.variables.size() != dimension + 1) {
                    result.attempt_trace.push_back(std::move(attempt));
                    continue;
                }
                const std::vector<double> independent(
                    run.variables.begin(), run.variables.end() - 1
                );
                const Held2StateEvaluation state = evaluator(
                    independent,
                    run.variables.back()
                );
                const std::array<double, 2> terminal_volume_bounds =
                    volume_bounds_evaluator(independent);
                std::vector<double> lower = physical_lower;
                std::vector<double> upper = physical_upper;
                lower.push_back(std::log(terminal_volume_bounds[0]));
                upper.push_back(std::log(terminal_volume_bounds[1]));
                double value = state.objective;
                for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                    value += multipliers[coordinate]
                        * (feed[coordinate] - independent[coordinate]);
                }
                attempt.terminal_modified_fractions = state.modified_fractions;
                attempt.terminal_volume = state.volume;
                attempt.objective = state.objective;
                attempt.lower_value = value;
                attempt.pressure_residual = state.pressure_stationarity_relative;
                attempt.pressure_passed = audit_held2_tolerance(
                    kHeld2RootPressure, attempt.pressure_residual
                ).passed;
                attempt.chart_jacobian_condition = matrix_condition_inf(
                    run.coordinate_jacobian, dimension
                );
                attempt.chart_kkt_inf_norm = 0.0;
                for (std::size_t column = 0; column < dimension; ++column) {
                    double residual = run.upper_bound_multipliers[column]
                        - run.lower_bound_multipliers[column];
                    for (std::size_t row = 0; row < dimension; ++row) {
                        residual += run.coordinate_jacobian[
                            row * dimension + column
                        ] * (state.gradient[row] - multipliers[row]);
                    }
                    attempt.chart_kkt_inf_norm = std::max(
                        attempt.chart_kkt_inf_norm, std::abs(residual)
                    );
                }
                attempt.chart_kkt_inf_norm = std::max(
                    attempt.chart_kkt_inf_norm,
                    std::abs(
                        state.gradient.back()
                        - run.lower_bound_multipliers.back()
                        + run.upper_bound_multipliers.back()
                    )
                );
                const Held2PhysicalKkt kkt = evaluate_stage_ii_physical_kkt(
                    state.gradient,
                    multipliers,
                    run.variables,
                    lower,
                    upper,
                    composition_sum_upper,
                    run.coordinate_jacobian,
                    run.lower_bound_multipliers,
                    run.upper_bound_multipliers
                );
                attempt.physical_kkt_inf_norm = kkt.stationarity_inf_norm;
                attempt.complementarity_inf_norm = kkt.complementarity;
                attempt.dual_pullback_inf_norm = kkt.reconstruction_inf_norm;
                attempt.dual_pullback_scale = kkt.reconstruction_scale;
                attempt.primal_inf_norm = kkt.primal_inf_norm;
                attempt.dual_sign_violation_inf_norm =
                    kkt.dual_sign_violation_inf_norm;
                attempt.dual_signs_valid = kkt.dual_signs_valid;
                attempt.physical_kkt_passed = attempt.pressure_passed
                    && std::isfinite(kkt.primal_inf_norm)
                    && std::isfinite(kkt.stationarity_inf_norm)
                    && std::isfinite(kkt.complementarity)
                    && std::isfinite(kkt.reconstruction_inf_norm)
                    && std::isfinite(kkt.reconstruction_scale)
                    && audit_held2_tolerance(
                           kHeld2Stage2Primal, kkt.primal_inf_norm
                       ).passed
                    && kkt.dual_signs_valid
                    && audit_held2_tolerance(
                           kHeld2Stage2Stationarity,
                           kkt.stationarity_inf_norm
                       ).passed
                    && audit_held2_tolerance(
                           kHeld2Stage2Complementarity,
                           kkt.complementarity
                       ).passed
                    && audit_held2_tolerance(
                           kHeld2Stage2DualPullback,
                           kkt.reconstruction_inf_norm,
                           kkt.reconstruction_scale
                       ).passed;
                const Held2StageIIStep5Assessment step5 =
                    assess_held2_stage_ii_step5(
                        upper_bound,
                        value,
                        attempt.physical_kkt_passed
                    );
                attempt.step5_qualified = step5.qualified;
                attempt.step5_reason = step5.reason;
                attempt.cut_eligible = attempt.step5_qualified;
                Held2ProgressEvent certificate_progress;
                certificate_progress.kind = Held2ProgressKind::Certificate;
                certificate_progress.stage = "STAGE II LOCAL KKT";
                certificate_progress.major_iteration = major;
                certificate_progress.attempt = attempt.attempt_id;
                certificate_progress.status = attempt.physical_kkt_passed
                    ? "kkt_passed"
                    : "kkt_failed";
                certificate_progress.primal_residual = kkt.primal_inf_norm;
                certificate_progress.dual_residual = kkt.stationarity_inf_norm;
                certificate_progress.complementarity = kkt.complementarity;
                observe_held2(observer, certificate_progress);
                const double gap = upper_bound - value;
                attempt.step6_gap = gap;
                attempt.step6_gap_passed = audit_held2_tolerance(
                    kHeld2Step6Gap, gap
                ).passed;
                bool fixed_volume_stationary = attempt.cut_eligible;
                if (attempt.cut_eligible) {
                    for (std::size_t coordinate = 0; coordinate < dimension;
                         ++coordinate) {
                        if (independent[coordinate]
                            <= physical_lower[coordinate]
                                + kHeld2BoundActivity.atol) {
                            continue;
                        }
                        const double gradient_scale = std::max(
                            std::abs(state.gradient[coordinate]),
                            std::abs(multipliers[coordinate])
                        );
                        attempt.fixed_volume_gradient_inf_norm = std::max(
                            attempt.fixed_volume_gradient_inf_norm,
                            std::abs(
                                state.gradient[coordinate]
                                    - multipliers[coordinate]
                            )
                        );
                        attempt.fixed_volume_gradient_scale = std::max(
                            attempt.fixed_volume_gradient_scale, gradient_scale
                        );
                        fixed_volume_stationary = fixed_volume_stationary
                            && audit_held2_tolerance(
                                   kHeld2Step6Gradient,
                                   state.gradient[coordinate]
                                       - multipliers[coordinate],
                                   gradient_scale
                               ).passed;
                    }
                }
                attempt.step6_gradient_passed = fixed_volume_stationary;
                attempt.step6_eligible = attempt.step5_qualified
                    && attempt.step6_gap_passed
                    && fixed_volume_stationary;
                if (attempt.step6_eligible) {
                    attempt.step6_rejection_reason = "eligible";
                } else if (!attempt.step5_qualified) {
                    attempt.step6_rejection_reason =
                        "step5_lower_not_qualified";
                } else if (!attempt.step6_gap_passed) {
                    attempt.step6_rejection_reason =
                        "same_major_gap_failed";
                } else {
                    attempt.step6_rejection_reason =
                        "fixed_volume_gradient_failed";
                }
                Held2ProgressEvent step6_progress;
                step6_progress.kind = Held2ProgressKind::Certificate;
                step6_progress.stage = "STAGE II STEP 6";
                step6_progress.major_iteration = major;
                step6_progress.attempt = attempt.attempt_id;
                step6_progress.status = attempt.step6_eligible
                    ? "eligible"
                    : "ineligible";
                step6_progress.reason = attempt.step6_eligible
                    ? "same_major_gap_and_fixed_volume_gradient"
                    : attempt.step6_rejection_reason;
                step6_progress.primal_residual = gap;
                step6_progress.dual_residual =
                    attempt.fixed_volume_gradient_inf_norm;
                observe_held2(observer, step6_progress);
                const auto basin = std::find_if(
                    traced_basins.begin(),
                    traced_basins.end(),
                    [&state](const Held2StateEvaluation& known) {
                        return maximum_abs_difference(
                                   known.modified_fractions,
                                   state.modified_fractions
                               ) <= kHeld2BasinDuplicateComposition.atol
                            && std::abs(
                                   std::log(known.volume)
                                       - std::log(state.volume)
                               ) <= kHeld2BasinDuplicateLogVolume.atol;
                    }
                );
                if (basin == traced_basins.end()) {
                    attempt.basin_id = static_cast<int>(traced_basins.size());
                    traced_basins.push_back(state);
                } else {
                    attempt.basin_id = static_cast<int>(
                        basin - traced_basins.begin()
                    );
                    ++result.duplicate_terminal_count;
                }
                if (attempt.cut_eligible && value < lower_bound) {
                    lower_bound = value;
                    best = make_cut(independent, state.volume, true);
                }
                if (attempt.step6_eligible) {
                    const bool duplicate_step6_state = std::any_of(
                        same_major_step6_cuts.begin(),
                        same_major_step6_cuts.end(),
                        [&state](const Cut& known) {
                            return maximum_abs_difference(
                                       known.state.modified_fractions,
                                       state.modified_fractions
                                   )
                                    <= kHeld2BasinDuplicateComposition.atol
                                && std::abs(
                                       std::log(known.state.volume)
                                           - std::log(state.volume)
                                   )
                                    <= kHeld2BasinDuplicateLogVolume.atol;
                        }
                    );
                    if (!duplicate_step6_state) {
                        same_major_step6_cuts.push_back(
                            make_cut(independent, state.volume, true)
                        );
                    }
                }
                const bool qualified = attempt.step5_qualified;
                const bool step6_eligible = attempt.step6_eligible;
                result.attempt_trace.push_back(std::move(attempt));
                if (qualified && !step6_eligible
                    && same_major_step6_cuts.empty()) {
                    step5_condition_satisfied = true;
                    break;
                }
                if (step6_eligible) {
                    step5_condition_satisfied = true;
                    if (same_major_step6_cuts.size() >= 2) {
                        break;
                    }
                }
            } catch (const std::exception& error) {
                attempt.callback_error = error.what();
                attempt.provider_status = "provider_failed";
                result.attempt_trace.push_back(std::move(attempt));
            }
        }
        ++result.major_iterations;
        result.lower_starts_per_iteration = std::max(
            result.lower_starts_per_iteration,
            static_cast<int>(
                result.attempt_trace.size() - attempt_trace_start
            )
        );
        if (!std::isfinite(lower_bound)) {
            result.local_attempts_truncated =
                result.local_attempts_truncated
                || representative_cap_truncated;
            result.outcome = representative_cap_truncated
                ? "resource_limit"
                : "indeterminate_step5_lower_not_below_upper";
            result.cut_count = static_cast<int>(cuts.size());
            result.distinct_basin_count = static_cast<int>(traced_basins.size());
            return result;
        }
        if (!step5_condition_satisfied) {
            result.outcome = "indeterminate_step5_completion";
            result.cut_count = static_cast<int>(cuts.size());
            result.distinct_basin_count = static_cast<int>(traced_basins.size());
            return result;
        }
        if (same_major_step6_cuts.size() == 1
            && representative_cap_truncated) {
            result.local_attempts_truncated = true;
        }
        result.bound_history.back().lower_bound = lower_bound;
        result.bound_history.back().lower_bound_available = true;
        upper_progress.lower_bound = lower_bound;
        upper_progress.gap = upper_bound - lower_bound;
        upper_progress.status = "lower_bound_certified";
        observe_held2(observer, upper_progress);
        const bool duplicate_cut = std::any_of(
            cuts.begin(),
            cuts.end(),
            [&best](const Cut& known) {
                return maximum_abs_difference(
                           known.state.modified_fractions,
                           best.state.modified_fractions
                       ) <= kHeld2BasinDuplicateComposition.atol
                    && std::abs(
                           std::log(known.state.volume)
                               - std::log(best.state.volume)
                       ) <= kHeld2BasinDuplicateLogVolume.atol;
            }
        );
        bool new_cut_added = false;
        if (!duplicate_cut) {
            cuts.push_back(std::move(best));
            new_cut_added = true;
        } else {
            const auto known = std::find_if(
                cuts.begin(),
                cuts.end(),
                [&best](const Cut& cut) {
                    return maximum_abs_difference(
                               cut.state.modified_fractions,
                               best.state.modified_fractions
                           ) <= kHeld2BasinDuplicateComposition.atol
                        && std::abs(
                               std::log(cut.state.volume)
                                   - std::log(best.state.volume)
                           ) <= kHeld2BasinDuplicateLogVolume.atol;
                }
            );
            if (known != cuts.end()) {
                known->candidate_phase_eligible = true;
            }
        }
        for (Cut& candidate_cut : same_major_step6_cuts) {
            const bool duplicate_candidate_cut = std::any_of(
                cuts.begin(),
                cuts.end(),
                [&candidate_cut](const Cut& known) {
                    return maximum_abs_difference(
                               known.state.modified_fractions,
                               candidate_cut.state.modified_fractions
                           ) <= kHeld2BasinDuplicateComposition.atol
                        && std::abs(
                               std::log(known.state.volume)
                                   - std::log(candidate_cut.state.volume)
                           ) <= kHeld2BasinDuplicateLogVolume.atol;
                }
            );
            if (!duplicate_candidate_cut) {
                cuts.push_back(std::move(candidate_cut));
                new_cut_added = true;
            } else {
                const auto known = std::find_if(
                    cuts.begin(),
                    cuts.end(),
                    [&candidate_cut](const Cut& cut) {
                        return maximum_abs_difference(
                                   cut.state.modified_fractions,
                                   candidate_cut.state.modified_fractions
                               ) <= kHeld2BasinDuplicateComposition.atol
                            && std::abs(
                                   std::log(cut.state.volume)
                                       - std::log(candidate_cut.state.volume)
                               ) <= kHeld2BasinDuplicateLogVolume.atol;
                    }
                );
                if (known != cuts.end()) {
                    known->candidate_phase_eligible = true;
                }
            }
        }
        for (const Held2StageIIPhysicalStart& representative :
             exploration.representatives) {
            Cut exploration_cut = make_cut(
                representative.independent_modified_fractions,
                representative.volume
            );
            const bool duplicate_exploration_cut = std::any_of(
                cuts.begin(),
                cuts.end(),
                [&exploration_cut](const Cut& known) {
                    return maximum_abs_difference(
                               known.state.modified_fractions,
                               exploration_cut.state.modified_fractions
                           ) <= kHeld2BasinDuplicateComposition.atol
                        && std::abs(
                               std::log(known.state.volume)
                                   - std::log(exploration_cut.state.volume)
                           ) <= kHeld2BasinDuplicateLogVolume.atol;
                }
            );
            if (!duplicate_exploration_cut) {
                cuts.push_back(std::move(exploration_cut));
                new_cut_added = true;
            }
        }
        result.candidates.clear();
        for (const Cut& cut : cuts) {
            if (!cut.candidate_phase_eligible) {
                continue;
            }
            double value = cut.state.objective;
            for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                value += multipliers[coordinate]
                    * (feed[coordinate] - cut.independent[coordinate]);
            }
            const double gap = upper_bound - value;
            bool fixed_volume_stationary = true;
            for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                if (cut.independent[coordinate]
                    <= physical_lower[coordinate] + kHeld2BoundActivity.atol) {
                    continue;
                }
                const double gradient_scale = std::max(
                    std::abs(cut.fixed_volume_gradient[coordinate]),
                    std::abs(multipliers[coordinate])
                );
                fixed_volume_stationary = fixed_volume_stationary
                    && audit_held2_tolerance(
                           kHeld2Step6Gradient,
                           cut.fixed_volume_gradient[coordinate]
                               - multipliers[coordinate],
                           gradient_scale
                       ).passed;
            }
            if (!audit_held2_tolerance(kHeld2Step6Gap, gap).passed
                || !audit_held2_tolerance(
                        kHeld2RootPressure,
                        cut.state.pressure_stationarity_relative
                    ).passed
                || !fixed_volume_stationary) {
                continue;
            }
            bool distinct = true;
            bool unresolved_identity = false;
            for (const Held2StageIICandidate& known : result.candidates) {
                const double composition_distance = maximum_abs_difference(
                    known.modified_fractions, cut.state.modified_fractions
                );
                const double log_volume_distance = std::abs(
                    known.phase_coordinate - std::log(cut.state.volume)
                );
                const bool duplicate =
                    composition_distance <= kHeld2BasinDuplicateComposition.atol
                    && log_volume_distance <= kHeld2BasinDuplicateLogVolume.atol;
                if (duplicate) {
                    distinct = false;
                    break;
                }
                const bool confidently_distinct =
                    audit_held2_tolerance(
                        kHeld2CandidateDistinctComposition,
                        composition_distance
                    ).passed
                    || audit_held2_tolerance(
                        kHeld2CandidateDistinctLogVolume,
                        log_volume_distance
                    ).passed;
                if (!confidently_distinct) {
                    distinct = false;
                    unresolved_identity = true;
                    break;
                }
            }
            if (unresolved_identity) {
                ++result.unresolved_candidate_identity_count;
            }
            if (distinct) {
                result.candidates.push_back({
                    cut.state.modified_fractions,
                    cut.independent,
                    cut.state.volume,
                    std::log(cut.state.volume),
                    gap,
                });
            }
        }
        result.cut_count = static_cast<int>(cuts.size());
        result.distinct_basin_count = static_cast<int>(traced_basins.size());
        if (result.unresolved_candidate_identity_count != 0) {
            result.outcome = "indeterminate_candidate_identity";
            return result;
        }
        if (result.candidates.size() > 1) {
            result.outcome = "candidate_set";
            return result;
        }
        if (!new_cut_added) {
            if (!request_direct_escalation) {
                request_direct_escalation = true;
                continue;
            }
            result.outcome = "indeterminate_finite_search_stalled";
            return result;
        }
    }
    result.outcome = "resource_limit";
    result.cut_count = static_cast<int>(cuts.size());
    result.distinct_basin_count = static_cast<int>(traced_basins.size());
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
    const double phase_charge_scale = std::max(
        charge_scale(charges, result.physical_phases[0]),
        charge_scale(charges, result.physical_phases[1])
    );
    double modified_potential_scale = 0.0;
    for (const auto& phase_potentials : result.phase_modified_potentials) {
        for (double potential : phase_potentials) {
            modified_potential_scale = std::max(
                modified_potential_scale, std::abs(potential)
            );
        }
    }
    result.certificate.accepted =
        audit_held2_tolerance(
            kHeld2Stage3ModifiedBalance,
            result.certificate.modified_balance_abs
        ).passed
        && audit_held2_tolerance(
            kHeld2Stage3ExplicitBalance,
            result.certificate.ordinary_balance_inf_norm
        ).passed
        && audit_held2_tolerance(
            kHeld2Stage3Charge,
            result.certificate.phase_charge_inf_norm,
            phase_charge_scale
        ).passed
        && audit_held2_tolerance(
            kHeld2Stage3Potential,
            result.certificate.modified_potential_gap,
            modified_potential_scale
        ).passed
        && audit_held2_tolerance(
            kHeld2Stage3Pressure,
            result.certificate.pressure_stationarity_inf_norm
        ).passed
        && audit_held2_tolerance(
            kHeld2Stage3Stationarity,
            result.certificate.reduced_kkt_inf_norm
        ).passed
        && audit_held2_tolerance(
            kHeld2Stage3FreeEnergyGap,
            result.certificate.enumeration_objective_gap
        ).passed;
    result.certificate.independent_evidence = audit_held2_tolerance(
        kHeld2Stage3FreeEnergyGap,
        result.certificate.enumeration_objective_gap
    ).passed;
    return result;
}

}  // namespace epcsaft_equilibrium
