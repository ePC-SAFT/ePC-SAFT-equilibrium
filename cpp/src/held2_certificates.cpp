#include "held2_certificates.hpp"

#include "held2_tolerances.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace epcsaft_equilibrium {
namespace {

bool finite_values(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

bool valid_bounds(
    const std::vector<double>& lower,
    const std::vector<double>& upper
) {
    if (lower.size() != upper.size()) {
        return false;
    }
    for (std::size_t index = 0; index < lower.size(); ++index) {
        if (std::isnan(lower[index]) || std::isnan(upper[index])
            || lower[index] == std::numeric_limits<double>::infinity()
            || upper[index] == -std::numeric_limits<double>::infinity()
            || lower[index] > upper[index]) {
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

std::size_t matrix_rank(std::vector<std::vector<double>> matrix) {
    if (matrix.empty()) {
        return 0;
    }
    const std::size_t rows = matrix.size();
    const std::size_t columns = matrix.front().size();
    double scale = 0.0;
    for (const auto& row : matrix) {
        scale = std::max(scale, maximum_abs(row));
    }
    std::size_t rank = 0;
    for (std::size_t column = 0;
         column < columns && rank < rows;
         ++column) {
        std::size_t pivot = rank;
        for (std::size_t row = rank + 1; row < rows; ++row) {
            if (std::abs(matrix[row][column])
                > std::abs(matrix[pivot][column])) {
                pivot = row;
            }
        }
        const double pivot_value = std::abs(matrix[pivot][column]);
        if (!audit_held2_tolerance(
                kHeld2Stage2KktRankPivot,
                pivot_value,
                scale
            ).passed) {
            continue;
        }
        std::swap(matrix[rank], matrix[pivot]);
        for (std::size_t row = rank + 1; row < rows; ++row) {
            const double factor =
                matrix[row][column] / matrix[rank][column];
            for (std::size_t local = column; local < columns; ++local) {
                matrix[row][local] -= factor * matrix[rank][local];
            }
        }
        ++rank;
    }
    return rank;
}

}  // namespace

Held2FarkasCertificate audit_held2_farkas_certificate(
    const std::vector<std::vector<double>>& matrix,
    const std::vector<double>& row_lower,
    const std::vector<double>& row_upper,
    const std::vector<double>& column_lower,
    const std::vector<double>& column_upper,
    const std::vector<double>& row_ray
) {
    Held2FarkasCertificate result;
    const std::size_t rows = matrix.size();
    const std::size_t columns = column_lower.size();
    if (rows != row_lower.size() || rows != row_upper.size()
        || rows != row_ray.size() || column_upper.size() != columns
        || std::any_of(
            matrix.begin(), matrix.end(),
            [columns](const auto& row) { return row.size() != columns; }
        )) {
        result.reason = "dimension_mismatch";
        return result;
    }
    if (!valid_bounds(row_lower, row_upper)
        || !valid_bounds(column_lower, column_upper)
        || !finite_values(row_ray)
        || std::any_of(
            matrix.begin(), matrix.end(),
            [](const auto& row) { return !finite_values(row); }
        )) {
        result.reason = "nonfinite_evidence";
        return result;
    }
    result.normalization = maximum_abs(row_ray);
    if (!std::isfinite(result.normalization) || result.normalization == 0.0) {
        result.reason = "zero_ray";
        return result;
    }
    result.row_multipliers = row_ray;
    for (double& value : result.row_multipliers) {
        value /= result.normalization;
    }

    double row_support = 0.0;
    result.row_sign_violation_inf = 0.0;
    std::vector<double> column_dual(columns, 0.0);
    for (std::size_t row = 0; row < rows; ++row) {
        const double multiplier = result.row_multipliers[row];
        const bool has_lower = std::isfinite(row_lower[row]);
        const bool has_upper = std::isfinite(row_upper[row]);
        if (!has_lower && !has_upper) {
            result.row_sign_violation_inf = std::max(
                result.row_sign_violation_inf, std::abs(multiplier)
            );
        } else if (!has_lower && multiplier < 0.0) {
            result.row_sign_violation_inf = std::max(
                result.row_sign_violation_inf, -multiplier
            );
        } else if (!has_upper && multiplier > 0.0) {
            result.row_sign_violation_inf = std::max(
                result.row_sign_violation_inf, multiplier
            );
        }
        if (multiplier > 0.0 && has_upper) {
            row_support += multiplier * row_upper[row];
        } else if (multiplier < 0.0 && has_lower) {
            row_support += multiplier * row_lower[row];
        }
        for (std::size_t column = 0; column < columns; ++column) {
            column_dual[column] += multiplier * matrix[row][column];
        }
    }
    if (!std::isfinite(row_support) || !finite_values(column_dual)) {
        result.reason = "nonfinite_derived_evidence";
        return result;
    }
    if (result.row_sign_violation_inf > 0.0
        || !audit_held2_tolerance(
            kHeld2FarkasRowSign, result.row_sign_violation_inf
        ).passed) {
        result.reason = "row_sign_failed";
        return result;
    }

    double column_minimum = 0.0;
    result.dual_feasibility_violation_inf = 0.0;
    for (std::size_t column = 0; column < columns; ++column) {
        const double dual = column_dual[column];
        const bool has_lower = std::isfinite(column_lower[column]);
        const bool has_upper = std::isfinite(column_upper[column]);
        if (!has_lower && !has_upper) {
            result.dual_feasibility_violation_inf = std::max(
                result.dual_feasibility_violation_inf, std::abs(dual)
            );
        } else if (!has_lower && dual > 0.0) {
            result.dual_feasibility_violation_inf = std::max(
                result.dual_feasibility_violation_inf, dual
            );
        } else if (!has_upper && dual < 0.0) {
            result.dual_feasibility_violation_inf = std::max(
                result.dual_feasibility_violation_inf, -dual
            );
        }
        if (dual >= 0.0 && has_lower) {
            column_minimum += dual * column_lower[column];
        } else if (dual < 0.0 && has_upper) {
            column_minimum += dual * column_upper[column];
        }
    }
    if (!std::isfinite(column_minimum)) {
        result.reason = "nonfinite_derived_evidence";
        return result;
    }
    if (result.dual_feasibility_violation_inf > 0.0
        || !audit_held2_tolerance(
            kHeld2FarkasDual, result.dual_feasibility_violation_inf
        ).passed) {
        result.reason = "dual_feasibility_failed";
        return result;
    }
    result.contradiction_margin = column_minimum - row_support;
    if (!std::isfinite(result.contradiction_margin)) {
        result.reason = "nonfinite_derived_evidence";
        return result;
    }
    result.contradiction_scale = std::max(
        {1.0, std::abs(column_minimum), std::abs(row_support)}
    );
    const Held2ToleranceAudit contradiction = audit_held2_tolerance(
        kHeld2FarkasContradiction,
        result.contradiction_margin,
        result.contradiction_scale
    );
    result.contradiction_threshold = contradiction.threshold;
    if (!contradiction.passed) {
        result.reason = "contradiction_margin_failed";
        return result;
    }
    result.reason = "certified";
    result.accepted = true;
    return result;
}

std::string adjudicate_held2_farkas_status(
    bool solver_infeasible,
    const std::optional<Held2FarkasCertificate>& certificate
) {
    return solver_infeasible && certificate && certificate->accepted
        ? "certified_infeasible" : "indeterminate";
}

Held2Step5KktCertificate audit_held2_step5_kkt(
    const std::vector<double>& variables,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    const std::vector<double>& objective_gradient,
    const std::vector<std::vector<double>>& constraints,
    const std::vector<double>& constraint_upper,
    const std::vector<double>& lower_bound_multipliers,
    const std::vector<double>& upper_bound_multipliers,
    const std::vector<double>& constraint_multipliers,
    double pullback_residual,
    double pullback_scale,
    double pressure_residual,
    bool same_major_iteration,
    bool step4_binding_valid,
    bool pressure_branch_valid
) {
    Held2Step5KktCertificate result;
    const std::size_t dimension = variables.size();
    result.same_major_iteration = same_major_iteration;
    result.step4_binding_valid = step4_binding_valid;
    result.pressure_branch_valid = pressure_branch_valid;
    result.pullback_residual_inf = pullback_residual;
    result.pullback_scale = pullback_scale;
    result.pressure_residual = pressure_residual;
    if (lower.size() != dimension || upper.size() != dimension
        || objective_gradient.size() != dimension
        || lower_bound_multipliers.size() != dimension
        || upper_bound_multipliers.size() != dimension
        || constraints.size() != constraint_upper.size()
        || constraints.size() != constraint_multipliers.size()
        || std::any_of(
            constraints.begin(), constraints.end(),
            [dimension](const auto& row) {
                return row.size() != dimension;
            }
        )) {
        result.reason = "dimension_mismatch";
        return result;
    }
    if (!valid_bounds(lower, upper) || !finite_values(variables)
        || !finite_values(objective_gradient)
        || !finite_values(lower_bound_multipliers)
        || !finite_values(upper_bound_multipliers)
        || !finite_values(constraint_upper)
        || !finite_values(constraint_multipliers)
        || !std::isfinite(pullback_residual)
        || !std::isfinite(pullback_scale) || pullback_scale < 0.0
        || !std::isfinite(pressure_residual)
        || std::any_of(
            constraints.begin(), constraints.end(),
            [](const auto& row) { return !finite_values(row); }
        )) {
        result.reason = "nonfinite_evidence";
        return result;
    }
    if (!same_major_iteration) {
        result.reason = "stale_major_iteration";
        return result;
    }
    if (!step4_binding_valid) {
        result.reason = "step4_binding_failed";
        return result;
    }
    if (!pressure_branch_valid) {
        result.reason = "pressure_branch_failed";
        return result;
    }

    result.primal_residual_inf = 0.0;
    std::vector<double> constraint_values(constraints.size(), 0.0);
    std::vector<std::vector<double>> active_rows;
    for (std::size_t index = 0; index < dimension; ++index) {
        result.primal_residual_inf = std::max({
            result.primal_residual_inf,
            lower[index] - variables[index],
            variables[index] - upper[index],
            0.0,
        });
        if (std::isfinite(lower[index])
            && audit_held2_tolerance(
                kHeld2BoundActivity, variables[index] - lower[index]
            ).passed) {
            std::vector<double> row(dimension, 0.0);
            row[index] = -1.0;
            active_rows.push_back(std::move(row));
        }
        if (std::isfinite(upper[index])
            && audit_held2_tolerance(
                kHeld2BoundActivity, upper[index] - variables[index]
            ).passed) {
            std::vector<double> row(dimension, 0.0);
            row[index] = 1.0;
            active_rows.push_back(std::move(row));
        }
    }
    for (std::size_t row = 0; row < constraints.size(); ++row) {
        constraint_values[row] = std::inner_product(
            constraints[row].begin(), constraints[row].end(),
            variables.begin(), 0.0
        );
        const double slack = constraint_upper[row] - constraint_values[row];
        result.primal_residual_inf = std::max(
            result.primal_residual_inf, -slack
        );
        if (audit_held2_tolerance(kHeld2BoundActivity, slack).passed) {
            active_rows.push_back(constraints[row]);
        }
    }
    if (!audit_held2_tolerance(
            kHeld2Stage2KktPrimal, result.primal_residual_inf
        ).passed) {
        result.reason = "primal_feasibility_failed";
        return result;
    }
    result.active_constraint_count = active_rows.size();
    result.active_jacobian_rank = matrix_rank(active_rows);
    if (result.active_jacobian_rank != result.active_constraint_count) {
        result.reason = "active_jacobian_rank_deficient";
        return result;
    }

    result.dual_sign_violation_inf = 0.0;
    for (double value : lower_bound_multipliers) {
        result.dual_sign_violation_inf = std::max(
            result.dual_sign_violation_inf, -value
        );
    }
    for (double value : upper_bound_multipliers) {
        result.dual_sign_violation_inf = std::max(
            result.dual_sign_violation_inf, -value
        );
    }
    for (double value : constraint_multipliers) {
        result.dual_sign_violation_inf = std::max(
            result.dual_sign_violation_inf, -value
        );
    }
    if (!audit_held2_tolerance(
            kHeld2Stage2KktDualSign, result.dual_sign_violation_inf
        ).passed) {
        result.reason = "dual_sign_failed";
        return result;
    }
    if (!audit_held2_tolerance(
            kHeld2Stage2KktPullback, pullback_residual, pullback_scale
        ).passed) {
        result.reason = "pullback_failed";
        return result;
    }
    if (!audit_held2_tolerance(
            kHeld2RootPressure, pressure_residual
        ).passed) {
        result.reason = "pressure_stationarity_failed";
        return result;
    }

    std::vector<double> stationarity = objective_gradient;
    for (std::size_t index = 0; index < dimension; ++index) {
        stationarity[index] -= lower_bound_multipliers[index];
        stationarity[index] += upper_bound_multipliers[index];
    }
    for (std::size_t row = 0; row < constraints.size(); ++row) {
        for (std::size_t column = 0; column < dimension; ++column) {
            stationarity[column] +=
                constraint_multipliers[row] * constraints[row][column];
        }
    }
    result.stationarity_residual_inf = maximum_abs(stationarity);
    if (!audit_held2_tolerance(
            kHeld2Stage2KktStationarity,
            result.stationarity_residual_inf
        ).passed) {
        result.reason = "stationarity_failed";
        return result;
    }

    result.complementarity_inf = 0.0;
    for (std::size_t index = 0; index < dimension; ++index) {
        if (std::isfinite(lower[index])) {
            result.complementarity_inf = std::max(
                result.complementarity_inf,
                std::abs(
                    lower_bound_multipliers[index]
                    * (variables[index] - lower[index])
                )
            );
        }
        if (std::isfinite(upper[index])) {
            result.complementarity_inf = std::max(
                result.complementarity_inf,
                std::abs(
                    upper_bound_multipliers[index]
                    * (upper[index] - variables[index])
                )
            );
        }
    }
    for (std::size_t row = 0; row < constraints.size(); ++row) {
        result.complementarity_inf = std::max(
            result.complementarity_inf,
            std::abs(
                constraint_multipliers[row]
                * (constraint_upper[row] - constraint_values[row])
            )
        );
    }
    if (!audit_held2_tolerance(
            kHeld2Stage2KktComplementarity,
            result.complementarity_inf
        ).passed) {
        result.reason = "complementarity_failed";
        return result;
    }
    result.reason = "certified";
    result.accepted = true;
    return result;
}

}  // namespace epcsaft_equilibrium
