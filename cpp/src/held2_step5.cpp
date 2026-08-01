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

inline constexpr double kStep5DiluteFaceMuInitial = 1.0e-16;
inline constexpr double kStep5DiluteFaceBoundPush = 1.0e-14;

struct Step5Assessment {
    bool qualified = false;
    double gap = 0.0;
    std::string reason;
};

bool finite_values(const std::vector<double>& values) {
    return std::all_of(
        values.begin(), values.end(),
        [](double value) { return std::isfinite(value); }
    );
}

bool finite_step5_terminal(const Held2StateEvaluation& terminal) {
    return std::isfinite(terminal.volume)
        && terminal.volume > 0.0
        && std::isfinite(terminal.objective)
        && finite_values(terminal.gradient)
        && finite_values(terminal.hessian)
        && finite_values(terminal.modified_fractions)
        && finite_values(terminal.physical_amounts)
        && finite_values(terminal.modified_potentials)
        && finite_values(terminal.chemical_potentials_over_rt)
        && std::isfinite(terminal.pressure_stationarity_relative)
        && std::isfinite(
            terminal.pressure_stationarity_derivative_log_volume
        );
}

struct Step5Start {
    std::vector<double> independent;
    double log_volume = 0.0;
    const char* family = "random_interior";
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

bool representation_equivalent(
    const Held2MPoint& left,
    const Held2MPoint& right
) {
    if (left.independent_modified_fractions.size()
        != right.independent_modified_fractions.size()
        || !std::isfinite(left.volume) || left.volume <= 0.0
        || !std::isfinite(right.volume) || right.volume <= 0.0
        || !audit_held2_tolerance(
            kHeld2MRepresentationEquivalent,
            std::log(left.volume) - std::log(right.volume)
        ).passed) {
        return false;
    }
    return std::equal(
        left.independent_modified_fractions.begin(),
        left.independent_modified_fractions.end(),
        right.independent_modified_fractions.begin(),
        [](double left_value, double right_value) {
            return audit_held2_tolerance(
                kHeld2MRepresentationEquivalent,
                left_value - right_value
            ).passed;
        }
    );
}

Step5Start make_step5_start(
    const Held2Coordinates& coordinates,
    const Held2PersistentState& state,
    const std::vector<double>& random_point,
    const std::array<double, 2>& volume_bounds,
    double total_ion_mole_fraction_max,
    int attempt,
    int random_start_cap
) {
    const std::size_t dimension = state.feed.size();
    const std::vector<double> random_independent =
        held2_map_unit_cube_to_independent_fractions(
            coordinates,
            std::vector<double>(
                random_point.begin(),
                random_point.begin()
                    + static_cast<std::ptrdiff_t>(dimension)
            ),
            total_ion_mole_fraction_max
        );
    const double lower_log_volume = std::log(volume_bounds[0]);
    const double upper_log_volume = std::log(volume_bounds[1]);
    const double random_log_volume = lower_log_volume
        + random_point.back()
            * (upper_log_volume - lower_log_volume);
    if (!state.step5_requires_new_member
        || attempt < random_start_cap || state.M.empty()) {
        return {random_independent, random_log_volume, "random_interior"};
    }
    const std::uint64_t structured_ordinal = static_cast<std::uint64_t>(
        attempt - random_start_cap
    );
    const std::uint64_t boundary_start_count =
        4 * static_cast<std::uint64_t>(dimension);
    if (structured_ordinal >= boundary_start_count) {
        const std::uint64_t retained_ordinal =
            structured_ordinal
            - boundary_start_count;
        const Held2MPoint& retained = state.M[
            static_cast<std::size_t>(retained_ordinal)
                % state.M.size()
        ];
        std::vector<double> shifted(dimension, 0.0);
        for (std::size_t index = 0; index < dimension; ++index) {
            shifted[index] = 0.99
                    * retained.independent_modified_fractions[index]
                + 0.01 * random_independent[index];
        }
        const double retained_log_volume = std::clamp(
            std::log(retained.volume), lower_log_volume, upper_log_volume
        );
        return {
            std::move(shifted),
            0.99 * retained_log_volume + 0.01 * random_log_volume,
            "shifted_retained",
        };
    }
    std::vector<double> boundary_cube(dimension, 0.5);
    const std::uint64_t boundary_ordinal = structured_ordinal;
    boundary_cube[
        static_cast<std::size_t>(boundary_ordinal) % dimension
    ] = ((boundary_ordinal / dimension) % 2 == 0) ? 0.0 : 1.0;
    const std::vector<double> boundary =
        held2_map_unit_cube_to_independent_fractions(
            coordinates, boundary_cube, total_ion_mole_fraction_max
        );
    std::vector<double> boundary_aware(dimension, 0.0);
    for (std::size_t index = 0; index < dimension; ++index) {
        boundary_aware[index] = 0.999 * boundary[index]
            + 0.001 * random_independent[index];
    }
    const double volume_fraction =
        ((boundary_ordinal / (2 * dimension)) % 2 == 0) ? 0.05 : 0.95;
    return {
        std::move(boundary_aware),
        lower_log_volume
            + volume_fraction * (upper_log_volume - lower_log_volume),
        "boundary_aware",
    };
}

bool same_coordinates(
    const Held2Coordinates& left,
    const Held2Coordinates& right
) {
    if (left.charges != right.charges
        || left.eliminated_index != right.eliminated_index
        || left.dependent_index != right.dependent_index
        || left.paper_to_provider_indices
            != right.paper_to_provider_indices
        || left.provider_to_paper_indices
            != right.provider_to_paper_indices
        || left.compact_to_paper_indices
            != right.compact_to_paper_indices
        || left.retained_indices != right.retained_indices
        || left.independent_indices != right.independent_indices
        || left.modified_factors != right.modified_factors
        || left.independent_lower_bounds
            != right.independent_lower_bounds
        || left.independent_upper_bounds
            != right.independent_upper_bounds
        || left.polytope_constraints.size()
            != right.polytope_constraints.size()) {
        return false;
    }
    for (std::size_t index = 0;
         index < left.polytope_constraints.size();
         ++index) {
        const Held2PolytopeConstraint& left_constraint =
            left.polytope_constraints[index];
        const Held2PolytopeConstraint& right_constraint =
            right.polytope_constraints[index];
        if (left_constraint.name != right_constraint.name
            || left_constraint.coefficients
                != right_constraint.coefficients
            || left_constraint.upper_bound
                != right_constraint.upper_bound) {
            return false;
        }
    }
    return true;
}

double maximum_abs(const std::vector<double>& values) {
    double result = 0.0;
    for (double value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

std::optional<std::vector<double>> solve_dense_system(
    std::vector<double> matrix,
    std::vector<double> right_hand_side
) {
    const std::size_t dimension = right_hand_side.size();
    if (matrix.size() != dimension * dimension) {
        return std::nullopt;
    }
    for (std::size_t column = 0; column < dimension; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < dimension; ++row) {
            if (std::abs(matrix[row * dimension + column])
                > std::abs(matrix[pivot * dimension + column])) {
                pivot = row;
            }
        }
        if (!std::isfinite(matrix[pivot * dimension + column])
            || std::abs(matrix[pivot * dimension + column]) < 1.0e-14) {
            return std::nullopt;
        }
        for (std::size_t local = column; local < dimension; ++local) {
            std::swap(
                matrix[column * dimension + local],
                matrix[pivot * dimension + local]
            );
        }
        std::swap(right_hand_side[column], right_hand_side[pivot]);
        for (std::size_t row = column + 1; row < dimension; ++row) {
            const double factor =
                matrix[row * dimension + column]
                / matrix[column * dimension + column];
            for (std::size_t local = column; local < dimension; ++local) {
                matrix[row * dimension + local] -=
                    factor * matrix[column * dimension + local];
            }
            right_hand_side[row] -= factor * right_hand_side[column];
        }
    }
    std::vector<double> solution(dimension, 0.0);
    for (std::size_t reverse = 0; reverse < dimension; ++reverse) {
        const std::size_t row = dimension - reverse - 1;
        double value = right_hand_side[row];
        for (std::size_t column = row + 1;
             column < dimension;
             ++column) {
            value -= matrix[row * dimension + column] * solution[column];
        }
        solution[row] = value / matrix[row * dimension + row];
    }
    if (!std::all_of(solution.begin(), solution.end(), [](double value) {
            return std::isfinite(value);
        })) {
        return std::nullopt;
    }
    return solution;
}

struct EqualityGaugeFit {
    double residual_inf = std::numeric_limits<double>::infinity();
    double normalization_multiplier =
        std::numeric_limits<double>::quiet_NaN();
    double charge_multiplier = std::numeric_limits<double>::quiet_NaN();
};

EqualityGaugeFit fit_physical_equality_gauges(
    const std::vector<double>& covector,
    const std::vector<double>& charges
) {
    EqualityGaugeFit result;
    if (covector.size() != charges.size() || covector.empty()
        || !std::all_of(covector.begin(), covector.end(), [](double value) {
            return std::isfinite(value);
        })
        || !std::all_of(charges.begin(), charges.end(), [](double value) {
            return std::isfinite(value);
        })) {
        return result;
    }
    const double count = static_cast<double>(covector.size());
    const double charge_sum =
        std::accumulate(charges.begin(), charges.end(), 0.0);
    const double charge_square_sum = std::inner_product(
        charges.begin(), charges.end(), charges.begin(), 0.0
    );
    const double covector_sum =
        std::accumulate(covector.begin(), covector.end(), 0.0);
    const double charge_covector_sum = std::inner_product(
        charges.begin(), charges.end(), covector.begin(), 0.0
    );
    const double determinant =
        count * charge_square_sum - charge_sum * charge_sum;
    if (!std::isfinite(determinant)
        || determinant <= kHeld2Stage2KktRankPivot.atol) {
        return result;
    }
    result.normalization_multiplier = (
        -covector_sum * charge_square_sum
        + charge_covector_sum * charge_sum
    ) / determinant;
    result.charge_multiplier = (
        -count * charge_covector_sum
        + charge_sum * covector_sum
    ) / determinant;
    result.residual_inf = 0.0;
    for (std::size_t index = 0; index < covector.size(); ++index) {
        result.residual_inf = std::max(
            result.residual_inf,
            std::abs(
                covector[index]
                + result.normalization_multiplier
                + result.charge_multiplier * charges[index]
            )
        );
    }
    return result;
}

struct Step5LocalRun {
    bool converged = false;
    std::string status = "not_run";
    std::vector<double> variables;
    std::vector<double> lower_bound_multipliers;
    std::vector<double> upper_bound_multipliers;
    std::vector<double> constraint_multipliers;
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
        variables_.assign(x, x + n);
        lower_bound_multipliers_.assign(z_lower, z_lower + n);
        upper_bound_multipliers_.assign(z_upper, z_upper + n);
        constraint_multipliers_.assign(lambda, lambda + m);
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
            lower_bound_multipliers_,
            upper_bound_multipliers_,
            constraint_multipliers_,
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
    std::vector<double> lower_bound_multipliers_;
    std::vector<double> upper_bound_multipliers_;
    std::vector<double> constraint_multipliers_;
    int iterations_ = 0;
};

Held2Step5KktCertificate certify_step5_terminal(
    const Held2Coordinates& coordinates,
    const Held2Step4Result& step4,
    const Held2PersistentState& state,
    std::uint64_t start_ordinal,
    const Step5LocalRun& run,
    const std::vector<double>& audited_solver_variables,
    const Held2StateEvaluation& terminal,
    const std::vector<double>& composition_lower,
    const std::vector<double>& composition_upper,
    const std::array<double, 2>& volume_bounds
) {
    const std::size_t dimension = composition_lower.size();
    std::vector<double> variables(
        audited_solver_variables.begin(),
        audited_solver_variables.end() - 1
    );
    variables.push_back(terminal.volume);
    std::vector<double> lower = composition_lower;
    std::vector<double> upper = composition_upper;
    lower.push_back(volume_bounds[0]);
    upper.push_back(volume_bounds[1]);
    std::vector<std::vector<double>> constraints;
    std::vector<double> constraint_upper;
    std::vector<std::size_t> constraint_indices;
    struct CanonicalBoundConstraint {
        std::size_t constraint_index;
        std::size_t coordinate;
        double coefficient;
    };
    std::vector<CanonicalBoundConstraint> canonical_bound_constraints;
    constraints.reserve(coordinates.polytope_constraints.size());
    constraint_upper.reserve(coordinates.polytope_constraints.size());
    constraint_indices.reserve(coordinates.polytope_constraints.size());
    canonical_bound_constraints.reserve(
        coordinates.polytope_constraints.size()
    );
    const auto canonical_variable_bound = [
        &composition_lower, &composition_upper, dimension
    ](const Held2PolytopeConstraint& constraint)
        -> std::optional<std::pair<std::size_t, double>> {
        if (constraint.coefficients.size() != dimension) {
            return std::nullopt;
        }
        std::optional<std::size_t> candidate;
        for (std::size_t coordinate = 0;
             coordinate < dimension;
             ++coordinate) {
            if (std::abs(constraint.coefficients[coordinate])
                <= kHeld2PolytopeFeasibility.atol) {
                continue;
            }
            if (candidate) {
                return std::nullopt;
            }
            candidate = coordinate;
        }
        if (!candidate) {
            return std::nullopt;
        }
        const double coefficient =
            constraint.coefficients[*candidate];
        const double implied_bound =
            constraint.upper_bound / coefficient;
        const double variable_bound = coefficient > 0.0
            ? composition_upper[*candidate]
            : composition_lower[*candidate];
        if (!std::isfinite(implied_bound)
            || !std::isfinite(variable_bound)
            || !audit_held2_tolerance(
                kHeld2BoundActivity,
                implied_bound - variable_bound
            ).passed) {
            return std::nullopt;
        }
        return std::pair{*candidate, coefficient};
    };
    for (std::size_t index = 0;
         index < coordinates.polytope_constraints.size();
         ++index) {
        const Held2PolytopeConstraint& constraint =
            coordinates.polytope_constraints[index];
        if (const auto canonical = canonical_variable_bound(constraint)) {
            canonical_bound_constraints.push_back({
                index, canonical->first, canonical->second,
            });
            continue;
        }
        std::vector<double> coefficients = constraint.coefficients;
        coefficients.push_back(0.0);
        constraints.push_back(std::move(coefficients));
        constraint_upper.push_back(constraint.upper_bound);
        constraint_indices.push_back(index);
    }

    std::vector<double> lower_multipliers =
        run.lower_bound_multipliers;
    std::vector<double> upper_multipliers =
        run.upper_bound_multipliers;
    const bool volume_dual_convertible =
        terminal.volume > 0.0 && std::isfinite(terminal.volume)
        && lower_multipliers.size() == dimension + 1
        && upper_multipliers.size() == dimension + 1;
    if (volume_dual_convertible) {
        lower_multipliers.back() /= terminal.volume;
        upper_multipliers.back() /= terminal.volume;
    }
    if (lower_multipliers.size() == variables.size()
        && upper_multipliers.size() == variables.size()) {
        for (std::size_t index = 0; index < variables.size(); ++index) {
            if (!audit_held2_tolerance(
                    kHeld2BoundActivity,
                    variables[index] - lower[index]
                ).passed
                && std::isfinite(lower_multipliers[index])
                && (
                    std::abs(run.lower_bound_multipliers[index])
                        <= kHeld2Stage2KktDualSign.atol
                    || audit_held2_tolerance(
                        kHeld2Stage2KktComplementarity,
                        lower_multipliers[index]
                            * (variables[index] - lower[index])
                    ).passed
                )) {
                lower_multipliers[index] = 0.0;
            }
            if (!audit_held2_tolerance(
                    kHeld2BoundActivity,
                    upper[index] - variables[index]
                ).passed
                && std::isfinite(upper_multipliers[index])
                && (
                    std::abs(run.upper_bound_multipliers[index])
                        <= kHeld2Stage2KktDualSign.atol
                    || audit_held2_tolerance(
                        kHeld2Stage2KktComplementarity,
                        upper_multipliers[index]
                            * (upper[index] - variables[index])
                    ).passed
                )) {
                upper_multipliers[index] = 0.0;
            }
        }
    }
    std::vector<double> constraint_multipliers;
    if (run.constraint_multipliers.size()
        == coordinates.polytope_constraints.size()) {
        constraint_multipliers.reserve(constraint_indices.size());
        for (std::size_t index : constraint_indices) {
            constraint_multipliers.push_back(
                run.constraint_multipliers[index]
            );
        }
    } else {
        constraint_multipliers = run.constraint_multipliers;
    }
    if (run.constraint_multipliers.size()
            == coordinates.polytope_constraints.size()
        && lower_multipliers.size() == variables.size()
        && upper_multipliers.size() == variables.size()) {
        for (const auto& canonical : canonical_bound_constraints) {
            double multiplier =
                run.constraint_multipliers[
                    canonical.constraint_index
                ];
            const Held2PolytopeConstraint& constraint =
                coordinates.polytope_constraints[
                    canonical.constraint_index
                ];
            const double slack = constraint.upper_bound
                - canonical.coefficient
                    * variables[canonical.coordinate];
            if (!audit_held2_tolerance(
                    kHeld2BoundActivity, slack
                ).passed
                && std::isfinite(multiplier)
                && audit_held2_tolerance(
                    kHeld2Stage2KktComplementarity,
                    multiplier * slack
                ).passed) {
                multiplier = 0.0;
            }
            if (canonical.coefficient > 0.0) {
                upper_multipliers[canonical.coordinate] +=
                    multiplier * canonical.coefficient;
            } else {
                lower_multipliers[canonical.coordinate] -=
                    multiplier * canonical.coefficient;
            }
        }
    }
    if (constraint_multipliers.size() == constraints.size()) {
        for (std::size_t row = 0; row < constraints.size(); ++row) {
            const double value = std::inner_product(
                constraints[row].begin(), constraints[row].end(),
                variables.begin(), 0.0
            );
            if (!audit_held2_tolerance(
                    kHeld2BoundActivity,
                    constraint_upper[row] - value
                ).passed
                && std::isfinite(constraint_multipliers[row])
                && audit_held2_tolerance(
                    kHeld2Stage2KktComplementarity,
                    constraint_multipliers[row]
                        * (constraint_upper[row] - value)
                ).passed) {
                constraint_multipliers[row] = 0.0;
            }
        }
    }

    std::vector<double> objective_gradient(dimension + 1, 0.0);
    double pullback_residual = std::numeric_limits<double>::infinity();
    double pullback_scale = 0.0;
    try {
        const std::vector<double> modified_potentials =
            held2_transform_modified_potentials(
                coordinates, terminal.chemical_potentials_over_rt
            );
        const auto dependent = std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            coordinates.dependent_index
        );
        if (dependent == coordinates.retained_indices.end()
            || terminal.gradient.size() != dimension + 1
            || state.multipliers.size() != dimension) {
            throw std::invalid_argument(
                "HELD2 Step-5 pullback dimensions changed"
            );
        }
        const std::size_t dependent_position =
            static_cast<std::size_t>(
                dependent - coordinates.retained_indices.begin()
            );
        pullback_residual = 0.0;
        for (std::size_t coordinate = 0;
             coordinate < dimension;
             ++coordinate) {
            const auto independent = std::find(
                coordinates.retained_indices.begin(),
                coordinates.retained_indices.end(),
                coordinates.independent_indices[coordinate]
            );
            if (independent == coordinates.retained_indices.end()) {
                throw std::invalid_argument(
                    "HELD2 Step-5 independent chart changed"
                );
            }
            const std::size_t independent_position =
                static_cast<std::size_t>(
                    independent - coordinates.retained_indices.begin()
                );
            const double reconstructed =
                modified_potentials[independent_position]
                - modified_potentials[dependent_position];
            pullback_residual = std::max(
                pullback_residual,
                std::abs(reconstructed - terminal.gradient[coordinate])
            );
            pullback_scale = std::max(
                {pullback_scale, std::abs(reconstructed),
                 std::abs(terminal.gradient[coordinate])}
            );
            objective_gradient[coordinate] =
                reconstructed - state.multipliers[coordinate];
        }
        objective_gradient.back() =
            terminal.gradient.back() / terminal.volume;
    } catch (...) {
        objective_gradient.clear();
    }
    const bool pressure_branch_valid =
        volume_dual_convertible
        && std::isfinite(volume_bounds[0])
        && std::isfinite(volume_bounds[1])
        && volume_bounds[0] > 0.0
        && volume_bounds[0] < volume_bounds[1]
        && audited_solver_variables.size() == dimension + 1
        && audit_held2_tolerance(
            kHeld2JointVolumeConsistency,
            audited_solver_variables.back()
                - std::log(terminal.volume)
        ).passed;
    bool step4_certificate_valid = false;
    if (step4.certificate) {
        try {
            step4_certificate_valid =
                step4.certificate->primal_feasible
                && step4.certificate->dual_feasible
                && audit_held2_tolerance(
                    kHeld2LpPrimal,
                    step4.certificate->primal_residual_inf
                ).passed
                && audit_held2_tolerance(
                    kHeld2LpDual,
                    step4.certificate->dual_residual_inf
                ).passed
                && audit_held2_tolerance(
                    kHeld2LpComplementarity,
                    step4.certificate->complementarity_inf
                ).passed;
        } catch (...) {
            step4_certificate_valid = false;
        }
    }
    bool cut_snapshot_valid = false;
    bool active_cut_ids_valid = false;
    try {
        if (step4.cut_snapshot.size() == state.M.size() + 1
            && step4.upper_bound && step4.multipliers
            && step4.multipliers->size() == dimension) {
            cut_snapshot_valid = true;
            std::vector<int> reconstructed_active_ids;
            for (std::size_t index = 0;
                 index < step4.cut_snapshot.size();
                 ++index) {
                const Held2Step4CutEvidence& cut =
                    step4.cut_snapshot[index];
                int expected_id = -1;
                double expected_intercept =
                    state.feed_reduced_gibbs;
                std::vector<double> expected_slopes(dimension, 0.0);
                if (index < state.M.size()) {
                    const Held2MPoint& point = state.M[index];
                    expected_id =
                        static_cast<int>(point.insertion_id);
                    expected_intercept = point.reduced_gibbs;
                    if (point.independent_modified_fractions.size()
                        != dimension) {
                        cut_snapshot_valid = false;
                        break;
                    }
                    for (std::size_t coordinate = 0;
                         coordinate < dimension;
                         ++coordinate) {
                        expected_slopes[coordinate] =
                            state.feed[coordinate]
                            - point.independent_modified_fractions[
                                coordinate
                            ];
                    }
                }
                if (cut.id != expected_id
                    || cut.intercept != expected_intercept
                    || cut.slopes != expected_slopes) {
                    cut_snapshot_valid = false;
                    break;
                }
                double value = cut.intercept;
                for (std::size_t coordinate = 0;
                     coordinate < dimension;
                     ++coordinate) {
                    value += cut.slopes[coordinate]
                        * (*step4.multipliers)[coordinate];
                }
                if (audit_held2_tolerance(
                        kHeld2LpActiveCut,
                        value - *step4.upper_bound
                    ).passed) {
                    reconstructed_active_ids.push_back(cut.id);
                }
            }
            active_cut_ids_valid = cut_snapshot_valid
                && !reconstructed_active_ids.empty()
                && reconstructed_active_ids
                    == step4.active_cut_ids;
        }
    } catch (...) {
        cut_snapshot_valid = false;
        active_cut_ids_valid = false;
    }
    const bool step4_binding_valid =
        step4.upper_bound && step4.multipliers
        && state.coordinates && step4.coordinate_snapshot
        && same_coordinates(coordinates, *state.coordinates)
        && same_coordinates(
            coordinates, *step4.coordinate_snapshot
        )
        && step4_certificate_valid && cut_snapshot_valid
        && active_cut_ids_valid
        && step4.upper_solve_count == state.upper_solve_count
        && *step4.upper_bound == state.upper_bound
        && *step4.multipliers == state.multipliers;
    Held2Step5KktCertificate certificate = audit_held2_step5_kkt(
        variables,
        lower,
        upper,
        objective_gradient,
        constraints,
        constraint_upper,
        lower_multipliers,
        upper_multipliers,
        constraint_multipliers,
        pullback_residual,
        pullback_scale,
        terminal.pressure_stationarity_relative,
        step4.major_iteration == state.major_iteration,
        step4_binding_valid,
        pressure_branch_valid
    );
    EqualityGaugeFit physical_stationarity;
    double physical_composition_residual =
        std::numeric_limits<double>::infinity();
    try {
        if (terminal.physical_amounts.size()
                != coordinates.charges.size()
            || terminal.chemical_potentials_over_rt.size()
                != coordinates.charges.size()
            || lower_multipliers.size() < dimension
            || upper_multipliers.size() < dimension
            || constraint_multipliers.size() != constraints.size()) {
            throw std::invalid_argument(
                "HELD2 Step-5 physical evidence dimensions changed"
            );
        }
        const double normalization = std::accumulate(
            terminal.physical_amounts.begin(),
            terminal.physical_amounts.end(),
            0.0
        );
        const double charge = std::inner_product(
            terminal.physical_amounts.begin(),
            terminal.physical_amounts.end(),
            coordinates.charges.begin(),
            0.0
        );
        physical_composition_residual = std::max(
            std::abs(normalization - 1.0), std::abs(charge)
        );
        const std::vector<double> audited_composition(
            audited_solver_variables.begin(),
            audited_solver_variables.begin()
                + static_cast<std::ptrdiff_t>(dimension)
        );
        const std::vector<double> reconstructed_physical =
            held2_lift_independent_fractions(
                coordinates, audited_composition
            );
        if (reconstructed_physical.size()
            != terminal.physical_amounts.size()) {
            throw std::invalid_argument(
                "HELD2 Step-5 reconstructed composition changed"
            );
        }
        for (double amount : terminal.physical_amounts) {
            if (!std::isfinite(amount)) {
                throw std::invalid_argument(
                    "HELD2 Step-5 physical composition is nonfinite"
                );
            }
            physical_composition_residual = std::max(
                physical_composition_residual, -amount
            );
        }
        for (std::size_t index = 0;
             index < reconstructed_physical.size();
             ++index) {
            physical_composition_residual = std::max(
                physical_composition_residual,
                std::abs(
                    terminal.physical_amounts[index]
                    - reconstructed_physical[index]
                )
            );
        }

        std::vector<double> physical_covector =
            terminal.chemical_potentials_over_rt;
        std::vector<double> reduced_dual_force(dimension, 0.0);
        for (std::size_t coordinate = 0;
             coordinate < dimension;
             ++coordinate) {
            reduced_dual_force[coordinate] =
                -lower_multipliers[coordinate]
                + upper_multipliers[coordinate];
        }
        for (std::size_t row = 0; row < constraints.size(); ++row) {
            for (std::size_t coordinate = 0;
                 coordinate < dimension;
                 ++coordinate) {
                reduced_dual_force[coordinate] +=
                    constraint_multipliers[row]
                    * constraints[row][coordinate];
            }
        }
        for (std::size_t coordinate = 0;
             coordinate < dimension;
             ++coordinate) {
            const std::size_t component =
                coordinates.independent_indices[coordinate];
            const auto retained = std::find(
                coordinates.retained_indices.begin(),
                coordinates.retained_indices.end(),
                component
            );
            if (retained == coordinates.retained_indices.end()) {
                throw std::invalid_argument(
                    "HELD2 Step-5 physical pullback chart changed"
                );
            }
            const std::size_t retained_position =
                static_cast<std::size_t>(
                    retained - coordinates.retained_indices.begin()
                );
            physical_covector[component] +=
                coordinates.modified_factors[retained_position]
                * (
                    -state.multipliers[coordinate]
                    + reduced_dual_force[coordinate]
                );
        }
        physical_stationarity = fit_physical_equality_gauges(
            physical_covector, coordinates.charges
        );
    } catch (...) {
        physical_composition_residual =
            std::numeric_limits<double>::infinity();
    }
    certificate.physical_composition_residual_inf =
        physical_composition_residual;
    certificate.physical_stationarity_residual_inf =
        physical_stationarity.residual_inf;
    certificate.normalization_multiplier =
        physical_stationarity.normalization_multiplier;
    certificate.charge_multiplier =
        physical_stationarity.charge_multiplier;
    if (certificate.accepted) {
        if (!std::isfinite(physical_composition_residual)
            || !audit_held2_tolerance(
                kHeld2CompositionSum,
                physical_composition_residual
            ).passed) {
            certificate.reason = "physical_composition_failed";
            certificate.accepted = false;
        } else if (!std::isfinite(physical_stationarity.residual_inf)
            || !audit_held2_tolerance(
                kHeld2Stage2KktStationarity,
                physical_stationarity.residual_inf
            ).passed) {
            certificate.reason = "physical_stationarity_failed";
            certificate.accepted = false;
        }
    }
    certificate.major_iteration = state.major_iteration;
    certificate.start_ordinal = start_ordinal;
    certificate.step4_upper_bound = state.upper_bound;
    certificate.step4_multipliers = state.multipliers;
    certificate.step4_active_cut_ids = step4.active_cut_ids;
    certificate.solver_variables = run.variables;
    certificate.audited_variables = audited_solver_variables;
    certificate.solver_lower_bound_multipliers =
        run.lower_bound_multipliers;
    certificate.solver_upper_bound_multipliers =
        run.upper_bound_multipliers;
    certificate.solver_constraint_multipliers =
        run.constraint_multipliers;
    certificate.lower_bound_multipliers = lower_multipliers;
    certificate.upper_bound_multipliers = upper_multipliers;
    certificate.constraint_multipliers = constraint_multipliers;
    certificate.physical_volume_lower = volume_bounds[0];
    certificate.physical_volume_upper = volume_bounds[1];
    certificate.pressure_derivative_log_volume =
        terminal.pressure_stationarity_derivative_log_volume;
    return certificate;
}

Step5LocalRun solve_step5_local(
    Held2StateEvaluator objective,
    const std::vector<Held2PolytopeConstraint>& constraints,
    const std::vector<double>& initial,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    bool near_bound_restart = false
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
    application->Options()->SetStringValue(
        "nlp_scaling_method", "gradient-based"
    );
    application->Options()->SetNumericValue("bound_relax_factor", 0.0);
    if (near_bound_restart) {
        application->Options()->SetStringValue("mu_strategy", "monotone");
        application->Options()->SetNumericValue(
            "mu_init", kStep5DiluteFaceMuInitial
        );
        application->Options()->SetNumericValue(
            "bound_push", kStep5DiluteFaceBoundPush
        );
        application->Options()->SetNumericValue(
            "bound_frac", kStep5DiluteFaceBoundPush
        );
    }
    if (application->Initialize() != Ipopt::Solve_Succeeded) {
        return {};
    }
    application->OptimizeTNLP(problem);
    return problem->result();
}

}  // namespace

bool retain_held2_m_point(
    Held2PersistentState& state,
    Held2MPoint& point
) {
    const auto equivalent = std::find_if(
        state.M.begin(),
        state.M.end(),
        [&point](const Held2MPoint& member) {
            return representation_equivalent(member, point);
        }
    );
    if (equivalent != state.M.end()) {
        point.insertion_id = equivalent->insertion_id;
        return false;
    }
    point.insertion_id = static_cast<std::uint64_t>(state.M.size());
    state.M.push_back(point);
    return true;
}

Held2Step5Result run_held2_step5(
    const Held2Step1Result& step1,
    const Held2Step4Result& step4,
    Held2PersistentState& state,
    const Held2StateEvaluator& evaluator,
    const Held2ResourceProfile& resources,
    Held2ProgressObserver* observer
) {
    Held2Step5Result result;
    result.timing.invocation_count = 1;
    if (step4.status != "complete" || !step4.certificate
        || !step4.certificate->primal_feasible
        || !step4.certificate->dual_feasible
        || step4.major_iteration != state.major_iteration
        || step4.upper_solve_count != state.upper_solve_count
        || !step1.coordinates || !state.coordinates
        || !step4.coordinate_snapshot
        || !same_coordinates(
            *step1.coordinates, *state.coordinates
        )
        || !same_coordinates(
            *step1.coordinates, *step4.coordinate_snapshot
        )
        || !step1.volume_bounds || !evaluator
        || resources.step5_start_cap < 1
        || resources.step5_recovery_start_cap
            < resources.step5_start_cap) {
        result.reason = "invalid_step5_input";
        return result;
    }
    const Held2Coordinates& coordinates = *step1.coordinates;
    const std::size_t dimension = state.feed.size();
    const std::size_t recovery_capacity = static_cast<std::size_t>(
        resources.step5_recovery_start_cap
    );
    const std::size_t nominal_capacity = static_cast<std::size_t>(
        resources.step5_start_cap
    );
    const std::size_t structured_capacity = std::min(
        4 * dimension + state.M.size(),
        recovery_capacity - nominal_capacity
    );
    const int attempt_cap = state.step5_requires_new_member
        ? resources.step5_recovery_start_cap
        : resources.step5_start_cap;
    const int random_start_cap = state.step5_requires_new_member
        ? static_cast<int>(recovery_capacity - structured_capacity)
        : resources.step5_start_cap;
    std::vector<double> lower = coordinates.independent_lower_bounds;
    std::vector<double> upper = coordinates.independent_upper_bounds;
    const std::size_t begin =
        static_cast<std::size_t>(state.next_start_ordinal);
    std::mt19937_64 random(0);
    random.discard(begin * (dimension + 1));
    std::uniform_real_distribution<double> unit;
    std::vector<std::vector<double>> starts(
        static_cast<std::size_t>(attempt_cap),
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
            for (int index = 0; index < attempt_cap; ++index) {
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
    double best = std::numeric_limits<double>::infinity();
    Held2MPoint best_point;
    std::optional<Held2MPoint> selected_new_point;
    for (int attempt = 0;
         attempt < attempt_cap;
        ++attempt) {
        const std::vector<double>& point =
            starts[static_cast<std::size_t>(attempt)];
        const std::uint64_t ordinal = state.next_start_ordinal++;
        ++result.starts_consumed;
        Step5Start start = make_step5_start(
            coordinates,
            state,
            point,
            *state.step5_volume_bounds,
            step1.total_ion_mole_fraction_max,
            attempt,
            random_start_cap
        );
        Held2LocalCertificate certificate;
        certificate.start_ordinal = ordinal;
        certificate.start_family = start.family;
        std::array<double, 2> start_volume_bounds;
        try {
            start_volume_bounds = counted_volume_bounds(start.independent);
        } catch (...) {
            certificate.solver_status = "start_volume_bounds_failed";
            result.attempts.push_back(certificate);
            continue;
        }
        (*state.step5_volume_bounds)[0] = std::min(
            (*state.step5_volume_bounds)[0], start_volume_bounds[0]
        );
        (*state.step5_volume_bounds)[1] = std::max(
            (*state.step5_volume_bounds)[1], start_volume_bounds[1]
        );
        const std::array<double, 2> solve_volume_bounds =
            *state.step5_volume_bounds;
        start.log_volume = std::clamp(
            start.log_volume,
            std::log(start_volume_bounds[0]),
            std::log(start_volume_bounds[1])
        );
        std::vector<double> initial = start.independent;
        initial.push_back(start.log_volume);
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
        Step5LocalRun run = solve_step5_local(
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
        if (!run.converged || run.variables.size() != dimension + 1
            || !std::all_of(
                run.variables.begin(), run.variables.end(),
                [](double value) { return std::isfinite(value); }
            )) {
            if (run.converged) {
                certificate.kkt = Held2Step5KktCertificate{};
                certificate.kkt->reason = "nonfinite_evidence";
            }
            result.attempts.push_back(certificate);
            continue;
        }
        std::vector<double> audited_variables = run.variables;
        std::vector<double> independent(
            run.variables.begin(), run.variables.end() - 1
        );
        Held2StateEvaluation terminal;
        try {
            terminal = counted_evaluator(
                independent, run.variables.back()
            );
        } catch (...) {
            certificate.kkt = Held2Step5KktCertificate{};
            certificate.kkt->reason = "terminal_evaluation_failed";
            result.attempts.push_back(certificate);
            continue;
        }
        double value = 0.0;
        const auto assess_terminal = [&]() {
            value = terminal.objective;
            for (std::size_t coordinate = 0;
                 coordinate < dimension;
                 ++coordinate) {
                value += state.multipliers[coordinate]
                    * (state.feed[coordinate] - independent[coordinate]);
            }
            certificate.local_value = value;
            std::array<double, 2> terminal_volume_bounds;
            try {
                terminal_volume_bounds =
                    counted_volume_bounds(independent);
            } catch (...) {
                certificate.finite_and_in_domain = false;
                certificate.accepted = false;
                certificate.kkt = Held2Step5KktCertificate{};
                certificate.kkt->reason =
                    "terminal_volume_bounds_failed";
                return;
            }
            certificate.finite_and_in_domain =
                std::isfinite(value)
                && finite_step5_terminal(terminal)
                && terminal.volume >= terminal_volume_bounds[0]
                && terminal.volume <= terminal_volume_bounds[1];
            certificate.kkt = certify_step5_terminal(
                coordinates,
                step4,
                state,
                ordinal,
                run,
                audited_variables,
                terminal,
                lower,
                upper,
                terminal_volume_bounds
            );
            certificate.accepted = certificate.finite_and_in_domain
                && certificate.kkt->accepted;
        };
        assess_terminal();
        if (!certificate.accepted && certificate.finite_and_in_domain
            && certificate.kkt
            && certificate.kkt->reason == "stationarity_failed"
            && assess_step5(state.upper_bound, value, true).qualified
            && run.lower_bound_multipliers.size() == dimension + 1) {
            std::vector<double> face_initial = run.variables;
            bool near_dilute_face = false;
            const double face_radius = std::sqrt(
                kHeld2Stage2KktComplementarity.atol
            );
            Held2DiluteFaceRestartEvidence face_evidence;
            face_evidence.face_radius = face_radius;
            face_evidence.bound_activity_atol =
                kHeld2BoundActivity.atol;
            face_evidence.dual_sign_atol =
                kHeld2Stage2KktDualSign.atol;
            face_evidence.complementarity_atol =
                kHeld2Stage2KktComplementarity.atol;
            for (std::size_t coordinate = 0;
                 coordinate < dimension;
                 ++coordinate) {
                const double distance =
                    run.variables[coordinate] - solver_lower[coordinate];
                const double multiplier =
                    run.lower_bound_multipliers[coordinate];
                if (std::isfinite(distance)
                    && std::isfinite(multiplier)
                    && distance > kHeld2BoundActivity.atol
                    && distance <= face_radius
                    && multiplier > kHeld2Stage2KktDualSign.atol
                    && multiplier * distance
                        <= kHeld2Stage2KktComplementarity.atol) {
                    face_initial[coordinate] = solver_lower[coordinate];
                    near_dilute_face = true;
                    face_evidence.coordinate_indices.push_back(
                        coordinate
                    );
                    face_evidence.provider_component_indices.push_back(
                        coordinates.independent_indices[coordinate]
                    );
                    face_evidence.lower_bound_distances.push_back(
                        distance
                    );
                    face_evidence.lower_bound_multipliers.push_back(
                        multiplier
                    );
                    face_evidence.complementarity_products.push_back(
                        multiplier * distance
                    );
                }
            }
            if (near_dilute_face) {
                face_evidence.attempted = true;
                Step5LocalRun face_run = solve_step5_local(
                    objective,
                    coordinates.polytope_constraints,
                    face_initial,
                    solver_lower,
                    solver_upper,
                    true
                );
                ++result.timing.optimizer_solves;
                result.timing.optimizer_iterations +=
                    static_cast<std::uint64_t>(face_run.iterations);
                face_evidence.solver_status = face_run.status;
                if (face_run.converged
                    && face_run.variables.size() == dimension + 1
                    && std::all_of(
                        face_run.variables.begin(),
                        face_run.variables.end(),
                        [](double candidate) {
                            return std::isfinite(candidate);
                        }
                    )) {
                    try {
                        std::vector<double> face_independent(
                            face_run.variables.begin(),
                            face_run.variables.end() - 1
                        );
                        Held2StateEvaluation face_terminal =
                            counted_evaluator(
                                face_independent,
                                face_run.variables.back()
                            );
                        run = std::move(face_run);
                        audited_variables = run.variables;
                        independent = std::move(face_independent);
                        terminal = std::move(face_terminal);
                        certificate.solver_status =
                            run.status + "_dilute_face";
                        assess_terminal();
                        face_evidence.accepted = certificate.accepted;
                    } catch (...) {
                        // Preserve the original rejected certificate.
                    }
                }
                certificate.dilute_face_restart =
                    std::move(face_evidence);
            }
        }
        if (!certificate.accepted && certificate.kkt
            && certificate.kkt->reason == "pressure_stationarity_failed"
            && assess_step5(state.upper_bound, value, true).qualified) {
            try {
                std::vector<double> polished = run.variables;
                for (int polish = 0; polish < 12; ++polish) {
                    const std::vector<double> composition(
                        polished.begin(), polished.end() - 1
                    );
                    const Held2StateEvaluation polished_state =
                        objective(composition, polished.back());
                    const double derivative =
                        polished_state
                            .pressure_stationarity_derivative_log_volume;
                    if (!audit_held2_tolerance(
                            kHeld2RootPressure,
                            polished_state.pressure_stationarity_relative
                        ).passed) {
                        if (!std::isfinite(derivative)
                            || std::abs(derivative) < 1.0e-12) {
                            break;
                        }
                        const double correction = std::clamp(
                            -polished_state.pressure_stationarity_relative
                                / derivative,
                            -1.0,
                            1.0
                        );
                        polished.back() = std::clamp(
                            polished.back() + correction,
                            solver_lower.back(),
                            solver_upper.back()
                        );
                    }
                    const Held2StateEvaluation composition_state =
                        objective(
                            std::vector<double>(
                                polished.begin(), polished.end() - 1
                            ),
                            polished.back()
                        );
                    std::vector<double> composition_gradient(
                        composition_state.gradient.begin(),
                        composition_state.gradient.begin()
                            + static_cast<std::ptrdiff_t>(dimension)
                    );
                    if (audit_held2_tolerance(
                            kHeld2Stage2KktStationarity,
                            maximum_abs(composition_gradient)
                        ).passed
                        && audit_held2_tolerance(
                            kHeld2RootPressure,
                            composition_state.pressure_stationarity_relative
                        ).passed) {
                        break;
                    }
                    std::vector<double> composition_hessian(
                        dimension * dimension, 0.0
                    );
                    for (std::size_t row = 0; row < dimension; ++row) {
                        for (std::size_t column = 0;
                             column < dimension;
                             ++column) {
                            composition_hessian[row * dimension + column] =
                                composition_state.hessian[
                                    row * (dimension + 1) + column
                                ];
                        }
                    }
                    for (double& value : composition_gradient) {
                        value = -value;
                    }
                    const auto direction = solve_dense_system(
                        std::move(composition_hessian),
                        std::move(composition_gradient)
                    );
                    if (!direction) {
                        break;
                    }
                    double step_scale = 1.0;
                    const double maximum_step = maximum_abs(*direction);
                    if (maximum_step > 0.1) {
                        step_scale = 0.1 / maximum_step;
                    }
                    bool advanced = false;
                    for (int line_search = 0;
                         line_search < 16;
                         ++line_search) {
                        std::vector<double> candidate(
                            polished.begin(), polished.end() - 1
                        );
                        for (std::size_t coordinate = 0;
                             coordinate < dimension;
                             ++coordinate) {
                            candidate[coordinate] +=
                                step_scale * (*direction)[coordinate];
                        }
                        bool feasible = true;
                        for (std::size_t coordinate = 0;
                             coordinate < dimension;
                             ++coordinate) {
                            feasible = feasible
                                && candidate[coordinate] >= lower[coordinate]
                                && candidate[coordinate] <= upper[coordinate];
                        }
                        for (const Held2PolytopeConstraint& constraint :
                             coordinates.polytope_constraints) {
                            feasible = feasible
                                && std::inner_product(
                                    constraint.coefficients.begin(),
                                    constraint.coefficients.end(),
                                    candidate.begin(),
                                    0.0
                                ) <= constraint.upper_bound;
                        }
                        if (feasible) {
                            std::copy(
                                candidate.begin(), candidate.end(),
                                polished.begin()
                            );
                            advanced = true;
                            break;
                        }
                        step_scale *= 0.5;
                    }
                    if (!advanced) {
                        break;
                    }
                }
                for (int pressure_polish = 0;
                     pressure_polish < 8;
                     ++pressure_polish) {
                    const std::vector<double> composition(
                        polished.begin(), polished.end() - 1
                    );
                    const Held2StateEvaluation pressure_state =
                        counted_evaluator(composition, polished.back());
                    if (audit_held2_tolerance(
                            kHeld2RootPressure,
                            pressure_state.pressure_stationarity_relative
                        ).passed) {
                        break;
                    }
                    const double derivative =
                        pressure_state
                            .pressure_stationarity_derivative_log_volume;
                    if (!std::isfinite(derivative)
                        || std::abs(derivative) < 1.0e-12) {
                        break;
                    }
                    polished.back() = std::clamp(
                        polished.back()
                            - pressure_state.pressure_stationarity_relative
                                / derivative,
                        solver_lower.back(),
                        solver_upper.back()
                    );
                }
                const std::vector<double> polished_composition(
                    polished.begin(), polished.end() - 1
                );
                const Held2StateEvaluation polished_terminal =
                    counted_evaluator(polished_composition, polished.back());
                if (audit_held2_tolerance(
                        kHeld2RootPressure,
                        polished_terminal.pressure_stationarity_relative
                    ).passed) {
                    audited_variables = std::move(polished);
                    certificate.solver_status =
                        run.status + "_pressure_polished";
                    independent = polished_composition;
                    terminal = polished_terminal;
                    assess_terminal();
                }
            } catch (...) {
                // Preserve the original rejected certificate.
            }
        }
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
            ).qualified) {
            continue;
        }
        Held2MPoint qualified_point{
            static_cast<std::uint64_t>(state.M.size()),
            independent,
            terminal.volume,
            std::numeric_limits<double>::quiet_NaN(),
            terminal.objective,
            terminal.gradient,
            "step5_local",
        };
        if (value < best) {
            best = value;
            best_point = qualified_point;
        }
        const bool changes_persistent_set = std::none_of(
            state.M.begin(),
            state.M.end(),
            [&qualified_point](const Held2MPoint& member) {
                return representation_equivalent(member, qualified_point);
            }
        );
        if (changes_persistent_set) {
            selected_new_point = std::move(qualified_point);
            state.step5_requires_new_member = false;
            break;
        }
        if (!state.step5_requires_new_member) {
            break;
        }
    }
    if (!std::isfinite(best)) {
        result.reason = "step5_start_budget_exhausted";
        return result;
    }
    if (state.step5_requires_new_member) {
        best_point.origin = "step5_equivalent_member";
        result.lower_value = best;
        result.terminal = best_point;
        result.reason = "step5_recovery_exhausted";
        result.timing.terminal_status = result.status;
        result.timing.terminal_reason = result.reason;
        return result;
    }
    Held2MPoint terminal_point = selected_new_point
        ? std::move(*selected_new_point) : best_point;
    state.lower_value = best;
    result.lower_value = best;
    if (!retain_held2_m_point(state, terminal_point)) {
        terminal_point.origin = "step5_equivalent_member";
        result.reason = "equivalent_member";
    } else {
        result.reason = "step5_complete";
    }
    result.terminal = terminal_point;
    result.status = "complete";
    result.timing.terminal_status = result.status;
    result.timing.terminal_reason = result.reason;
    result.timing.next_step = 6;
    return result;
}

}  // namespace epcsaft_equilibrium
