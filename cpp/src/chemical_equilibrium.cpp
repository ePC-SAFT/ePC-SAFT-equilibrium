#include "chemical_equilibrium.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

namespace epcsaft_equilibrium {
namespace {

constexpr double kResidualMultiplier = 4096.0;

void require_finite(double value, const char* field) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(field) + " must be finite");
    }
}

void require_finite_vector(const std::vector<double>& values, const char* field) {
    if (!std::all_of(values.begin(), values.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::invalid_argument(std::string(field) + " must be finite");
    }
}

double matrix_scale(const DenseMatrix& matrix) {
    double scale = 0.0;
    for (double value : matrix.values) {
        scale = std::max(scale, std::abs(value));
    }
    return scale;
}

double numerical_tolerance(double scale, std::size_t dimension) {
    return kResidualMultiplier * std::numeric_limits<double>::epsilon()
        * std::max(1.0, scale) * static_cast<double>(std::max<std::size_t>(1, dimension));
}

double row_dot(const std::vector<double>& left, const std::vector<double>& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument("row dimensions do not match");
    }
    return std::inner_product(left.begin(), left.end(), right.begin(), 0.0);
}

double row_norm(const std::vector<double>& row) {
    return std::sqrt(row_dot(row, row));
}

std::vector<double> matrix_row(const DenseMatrix& matrix, std::size_t row) {
    std::vector<double> result(matrix.columns, 0.0);
    for (std::size_t column = 0; column < matrix.columns; ++column) {
        result[column] = matrix(row, column);
    }
    return result;
}

DenseMatrix matrix_from_rows(const std::vector<std::vector<double>>& rows) {
    if (rows.empty()) {
        return {};
    }
    const std::size_t columns = rows.front().size();
    DenseMatrix result{rows.size(), columns, {}};
    result.values.reserve(rows.size() * columns);
    for (const std::vector<double>& row : rows) {
        if (row.size() != columns) {
            throw std::invalid_argument("rows must have a common column count");
        }
        result.values.insert(result.values.end(), row.begin(), row.end());
    }
    return result;
}

std::vector<std::vector<double>> orthonormal_rows(
    const std::vector<std::vector<double>>& rows
) {
    std::vector<std::vector<double>> result;
    result.reserve(rows.size());
    for (const std::vector<double>& row : rows) {
        std::vector<double> residual = row;
        for (int pass = 0; pass < 2; ++pass) {
            for (const std::vector<double>& basis : result) {
                const double projection = row_dot(residual, basis);
                for (std::size_t column = 0; column < residual.size(); ++column) {
                    residual[column] -= projection * basis[column];
                }
            }
        }
        const double residual_norm = row_norm(residual);
        if (residual_norm <= numerical_tolerance(1.0, row.size())) {
            throw std::invalid_argument("row basis is numerically deficient");
        }
        for (double& value : residual) {
            value /= residual_norm;
        }
        result.push_back(std::move(residual));
    }
    return result;
}

double scaled_residual_norm(
    const std::vector<double>& row,
    const std::vector<std::vector<double>>& basis
) {
    std::vector<double> residual = row;
    for (int pass = 0; pass < 2; ++pass) {
        for (const std::vector<double>& vector : basis) {
            const double projection = row_dot(residual, vector);
            for (std::size_t column = 0; column < residual.size(); ++column) {
                residual[column] -= projection * vector[column];
            }
        }
    }
    const double scale = row_norm(row);
    return scale == 0.0 ? 0.0 : row_norm(residual) / scale;
}

std::vector<std::size_t> select_independent_rows(const DenseMatrix& matrix) {
    if (matrix.values.size() != matrix.rows * matrix.columns) {
        throw std::invalid_argument("matrix dimensions do not match storage");
    }
    std::vector<std::size_t> selected;
    std::vector<std::vector<double>> basis;
    std::vector<bool> used(matrix.rows, false);
    const double tolerance = numerical_tolerance(
        1.0, std::max(matrix.rows, matrix.columns)
    );
    while (selected.size() < matrix.columns) {
        std::size_t best = matrix.rows;
        double best_score = tolerance;
        for (std::size_t candidate = 0; candidate < matrix.rows; ++candidate) {
            if (used[candidate]) {
                continue;
            }
            const double score = scaled_residual_norm(matrix_row(matrix, candidate), basis);
            if (score > best_score) {
                best = candidate;
                best_score = score;
            }
        }
        if (best == matrix.rows) {
            break;
        }
        used[best] = true;
        selected.push_back(best);
        basis = orthonormal_rows([&] {
            std::vector<std::vector<double>> rows;
            rows.reserve(selected.size());
            for (std::size_t index : selected) {
                rows.push_back(matrix_row(matrix, index));
            }
            return rows;
        }());
    }
    return selected;
}

struct RowQr {
    DenseMatrix upper;
    std::vector<std::vector<double>> orthonormal;
};

RowQr factor_row_basis(const DenseMatrix& basis) {
    RowQr result{
        DenseMatrix{basis.rows, basis.rows, std::vector<double>(basis.rows * basis.rows, 0.0)},
        {},
    };
    result.orthonormal.reserve(basis.rows);
    for (std::size_t column = 0; column < basis.rows; ++column) {
        std::vector<double> residual = matrix_row(basis, column);
        for (int pass = 0; pass < 2; ++pass) {
            for (std::size_t row = 0; row < result.orthonormal.size(); ++row) {
                const double projection = row_dot(residual, result.orthonormal[row]);
                result.upper(row, column) += projection;
                for (std::size_t species = 0; species < residual.size(); ++species) {
                    residual[species] -= projection * result.orthonormal[row][species];
                }
            }
        }
        const double diagonal = row_norm(residual);
        if (diagonal <= numerical_tolerance(1.0, basis.columns)) {
            throw std::invalid_argument("reaction matrix rank is deficient");
        }
        result.upper(column, column) = diagonal;
        for (double& value : residual) {
            value /= diagonal;
        }
        result.orthonormal.push_back(std::move(residual));
    }
    return result;
}

std::vector<double> row_coordinates(
    const std::vector<double>& row,
    const RowQr& factor
) {
    std::vector<double> projected(factor.orthonormal.size(), 0.0);
    for (std::size_t index = 0; index < factor.orthonormal.size(); ++index) {
        projected[index] = row_dot(row, factor.orthonormal[index]);
    }
    std::vector<double> coordinates(projected.size(), 0.0);
    for (std::size_t reverse = projected.size(); reverse > 0; --reverse) {
        const std::size_t row_index = reverse - 1;
        double value = projected[row_index];
        for (std::size_t column = row_index + 1; column < coordinates.size(); ++column) {
            value -= factor.upper(row_index, column) * coordinates[column];
        }
        coordinates[row_index] = value / factor.upper(row_index, row_index);
    }
    return coordinates;
}

struct ReferenceReconstruction {
    std::vector<double> reference;
    std::size_t rank = 0;
    double qr_diagonal_ratio = 0.0;
};

ReferenceReconstruction construct_minimum_norm_reference(
    const DenseMatrix& reaction_matrix,
    const std::vector<double>& ln_k
) {
    if (reaction_matrix.rows == 0 || reaction_matrix.rows > reaction_matrix.columns
        || ln_k.size() != reaction_matrix.rows) {
        throw std::invalid_argument("reference reconstruction linear system is invalid");
    }

    const std::size_t species_count = reaction_matrix.columns;
    const std::size_t reaction_count = reaction_matrix.rows;
    DenseMatrix factor{
        species_count,
        reaction_count,
        std::vector<double>(species_count * reaction_count, 0.0),
    };
    std::vector<double> scaled_right_hand_side(reaction_count, 0.0);
    for (std::size_t reaction = 0; reaction < reaction_count; ++reaction) {
        double row_norm = 0.0;
        for (std::size_t species = 0; species < species_count; ++species) {
            row_norm = std::hypot(row_norm, reaction_matrix(reaction, species));
        }
        if (row_norm == 0.0) {
            throw std::invalid_argument("reaction matrix rank is deficient");
        }
        scaled_right_hand_side[reaction] = -ln_k[reaction] / row_norm;
        for (std::size_t species = 0; species < species_count; ++species) {
            factor(species, reaction) = reaction_matrix(reaction, species) / row_norm;
        }
    }

    std::vector<std::size_t> permutation(reaction_count, 0);
    std::iota(permutation.begin(), permutation.end(), 0);
    std::vector<double> reflector_scales(reaction_count, 0.0);
    std::vector<double> diagonal_magnitudes(reaction_count, 0.0);
    const double rank_tolerance = numerical_tolerance(1.0, species_count);

    for (std::size_t column = 0; column < reaction_count; ++column) {
        std::size_t pivot = column;
        double pivot_norm = 0.0;
        for (std::size_t candidate = column; candidate < reaction_count; ++candidate) {
            double trailing_norm = 0.0;
            for (std::size_t row = column; row < species_count; ++row) {
                trailing_norm = std::hypot(trailing_norm, factor(row, candidate));
            }
            if (trailing_norm > pivot_norm) {
                pivot = candidate;
                pivot_norm = trailing_norm;
            }
        }
        if (pivot_norm <= rank_tolerance) {
            throw std::invalid_argument("reaction matrix rank is deficient");
        }
        if (pivot != column) {
            for (std::size_t row = 0; row < species_count; ++row) {
                std::swap(factor(row, pivot), factor(row, column));
            }
            std::swap(permutation[pivot], permutation[column]);
        }

        double column_norm = 0.0;
        for (std::size_t row = column; row < species_count; ++row) {
            column_norm = std::hypot(column_norm, factor(row, column));
        }
        const double leading_value = factor(column, column);
        const double diagonal = -std::copysign(column_norm, leading_value);
        const double reflector_leading = leading_value - diagonal;
        if (std::abs(diagonal) <= rank_tolerance || reflector_leading == 0.0) {
            throw std::invalid_argument("reaction matrix rank is deficient");
        }
        reflector_scales[column] = (diagonal - leading_value) / diagonal;
        factor(column, column) = diagonal;
        for (std::size_t row = column + 1; row < species_count; ++row) {
            factor(row, column) /= reflector_leading;
        }
        for (std::size_t other = column + 1; other < reaction_count; ++other) {
            double projection = factor(column, other);
            for (std::size_t row = column + 1; row < species_count; ++row) {
                projection += factor(row, column) * factor(row, other);
            }
            projection *= reflector_scales[column];
            factor(column, other) -= projection;
            for (std::size_t row = column + 1; row < species_count; ++row) {
                factor(row, other) -= factor(row, column) * projection;
            }
        }
        diagonal_magnitudes[column] = std::abs(diagonal);
    }

    std::vector<double> projected_reference(reaction_count, 0.0);
    for (std::size_t row = 0; row < reaction_count; ++row) {
        double value = scaled_right_hand_side[permutation[row]];
        for (std::size_t column = 0; column < row; ++column) {
            value -= factor(column, row) * projected_reference[column];
        }
        projected_reference[row] = value / factor(row, row);
    }

    std::vector<double> reference(species_count, 0.0);
    std::copy(projected_reference.begin(), projected_reference.end(), reference.begin());
    for (std::size_t reverse = reaction_count; reverse > 0; --reverse) {
        const std::size_t column = reverse - 1;
        double projection = reference[column];
        for (std::size_t row = column + 1; row < species_count; ++row) {
            projection += factor(row, column) * reference[row];
        }
        projection *= reflector_scales[column];
        reference[column] -= projection;
        for (std::size_t row = column + 1; row < species_count; ++row) {
            reference[row] -= factor(row, column) * projection;
        }
    }

    const auto [minimum_diagonal, maximum_diagonal] = std::minmax_element(
        diagonal_magnitudes.begin(), diagonal_magnitudes.end()
    );
    return {
        std::move(reference),
        reaction_count,
        *minimum_diagonal / *maximum_diagonal,
    };
}

double reference_residual_inf_norm(
    const DenseMatrix& reaction_matrix,
    const std::vector<double>& reference,
    const std::vector<double>& ln_k
) {
    double norm = 0.0;
    for (std::size_t reaction = 0; reaction < reaction_matrix.rows; ++reaction) {
        double value = ln_k[reaction];
        for (std::size_t species = 0; species < reaction_matrix.columns; ++species) {
            value += reaction_matrix(reaction, species) * reference[species];
        }
        norm = std::max(norm, std::abs(value));
    }
    return norm;
}

double charge_reaction_inf_norm(
    const std::vector<int>& charges,
    const DenseMatrix& reaction_matrix
) {
    double norm = 0.0;
    for (std::size_t reaction = 0; reaction < reaction_matrix.rows; ++reaction) {
        double value = 0.0;
        for (std::size_t species = 0; species < reaction_matrix.columns; ++species) {
            value += static_cast<double>(charges[species])
                * reaction_matrix(reaction, species);
        }
        norm = std::max(norm, std::abs(value));
    }
    return norm;
}

void validate_identities(const ReactionSystemInput& input) {
    const std::size_t species_count = input.species_ids.size();
    if (species_count < 2) {
        throw std::invalid_argument("reaction system requires at least two species");
    }
    std::unordered_set<std::string> unique_ids;
    for (const std::string& species_id : input.species_ids) {
        if (species_id.empty() || !unique_ids.insert(species_id).second) {
            throw std::invalid_argument("species identities must be nonempty and unique");
        }
    }
    if (input.charges.size() != species_count) {
        throw std::invalid_argument("species charges do not match species count");
    }
    if (input.provider_fingerprint.empty()) {
        throw std::invalid_argument("Provider fingerprint is incomplete");
    }
}

void validate_reference_records(const ReactionSystemInput& input) {
    if (input.equilibrium_constant_records.size() != input.ln_k.size()) {
        throw std::invalid_argument("equilibrium-constant source records are incomplete");
    }
    for (const EquilibriumConstantRecord& record : input.equilibrium_constant_records) {
        if (record.source_id.empty()) {
            throw std::invalid_argument("equilibrium constant source identity is incomplete");
        }
        if (record.reference_id.empty()) {
            throw std::invalid_argument("equilibrium constant reference identity is incomplete");
        }
        if (record.reaction_orientation != "products_positive") {
            throw std::invalid_argument(
                "equilibrium constant reaction orientation must be products_positive"
            );
        }
        if (record.conversion_id != "already-provider-basis") {
            throw std::invalid_argument(
                "equilibrium constant conversion provenance is not already-provider-basis"
            );
        }
        if (record.reference_id != "provider-helmholtz-coordinate-basis") {
            throw std::invalid_argument(
                "equilibrium constant reference is not the Provider chemical basis"
            );
        }
        if (!record.dimensionless) {
            throw std::invalid_argument("lnK must be explicitly dimensionless");
        }
        require_finite(record.temperature_k, "equilibrium constant temperature");
        require_finite(record.pressure_pa, "equilibrium constant pressure");
        if (record.temperature_k != input.temperature_k) {
            throw std::invalid_argument("equilibrium constant temperature binding does not match");
        }
        if (record.pressure_pa != input.pressure_pa) {
            throw std::invalid_argument("equilibrium constant pressure binding does not match");
        }
    }
}

DenseMatrix selected_rows(
    const DenseMatrix& matrix,
    const std::vector<std::size_t>& row_indices
) {
    std::vector<std::vector<double>> rows;
    rows.reserve(row_indices.size());
    for (std::size_t row : row_indices) {
        rows.push_back(matrix_row(matrix, row));
    }
    return matrix_from_rows(rows);
}

double reaction_conservation_inf_norm(
    const std::vector<double>& molar_masses,
    const DenseMatrix& balance_matrix,
    const DenseMatrix& reaction_matrix
) {
    double norm = 0.0;
    for (std::size_t reaction = 0; reaction < reaction_matrix.rows; ++reaction) {
        double mass_value = 0.0;
        for (std::size_t species = 0; species < reaction_matrix.columns; ++species) {
            mass_value += molar_masses[species] * reaction_matrix(reaction, species);
        }
        norm = std::max(norm, std::abs(mass_value));
        for (std::size_t balance = 0; balance < balance_matrix.rows; ++balance) {
            double balance_value = 0.0;
            for (std::size_t species = 0; species < reaction_matrix.columns; ++species) {
                balance_value += balance_matrix(balance, species)
                    * reaction_matrix(reaction, species);
            }
            norm = std::max(norm, std::abs(balance_value));
        }
    }
    return norm;
}

std::vector<std::vector<double>> seeded_balance_rows(
    const ReactionSystemInput& input,
    bool& charge_is_independent
) {
    std::vector<std::vector<double>> basis_rows;
    std::vector<std::vector<double>> orthonormal;
    const double rank_tolerance = numerical_tolerance(1.0, input.species_ids.size());
    const auto add_seed = [&](const std::vector<double>& row) {
        const double score = scaled_residual_norm(row, orthonormal);
        if (score <= rank_tolerance) {
            return false;
        }
        basis_rows.push_back(row);
        orthonormal = orthonormal_rows(basis_rows);
        return true;
    };

    add_seed(input.molar_masses_kg_per_mol);
    std::vector<double> charge_row(input.species_ids.size(), 0.0);
    for (std::size_t species = 0; species < input.charges.size(); ++species) {
        charge_row[species] = static_cast<double>(input.charges[species]);
    }
    charge_is_independent = !std::all_of(
        charge_row.begin(), charge_row.end(), [](double value) { return value == 0.0; }
    ) && add_seed(charge_row);

    std::vector<std::vector<double>> retained_balance_rows;
    for (std::size_t row = 0; row < input.balance_matrix.rows; ++row) {
        const std::vector<double> candidate = matrix_row(input.balance_matrix, row);
        const double score = scaled_residual_norm(candidate, orthonormal);
        if (score <= rank_tolerance) {
            continue;
        }
        retained_balance_rows.push_back(candidate);
        basis_rows.push_back(candidate);
        orthonormal = orthonormal_rows(basis_rows);
    }

    std::vector<std::vector<double>> result;
    result.reserve(1 + retained_balance_rows.size());
    result.push_back(input.molar_masses_kg_per_mol);
    result.insert(result.end(), retained_balance_rows.begin(), retained_balance_rows.end());
    return result;
}

double matrix_reconstruction_inf_norm(
    const DenseMatrix& transform,
    const DenseMatrix& basis,
    const DenseMatrix& supplied
) {
    double norm = 0.0;
    for (std::size_t row = 0; row < supplied.rows; ++row) {
        for (std::size_t species = 0; species < supplied.columns; ++species) {
            double value = supplied(row, species);
            for (std::size_t basis_row = 0; basis_row < basis.rows; ++basis_row) {
                value -= transform(row, basis_row) * basis(basis_row, species);
            }
            norm = std::max(norm, std::abs(value));
        }
    }
    return norm;
}

double constant_reconstruction_inf_norm(
    const DenseMatrix& transform,
    const std::vector<double>& independent,
    const std::vector<double>& supplied
) {
    double norm = 0.0;
    for (std::size_t row = 0; row < supplied.size(); ++row) {
        double value = -supplied[row];
        for (std::size_t basis = 0; basis < independent.size(); ++basis) {
            value += transform(row, basis) * independent[basis];
        }
        norm = std::max(norm, std::abs(value));
    }
    return norm;
}

}  // namespace

double DenseMatrix::operator()(std::size_t row, std::size_t column) const {
    return values.at(row * columns + column);
}

double& DenseMatrix::operator()(std::size_t row, std::size_t column) {
    return values.at(row * columns + column);
}

std::size_t AmountChart::coordinate_count() const {
    if (!ionic()) {
        return neutral_indices.size();
    }
    return 1 + (cation_indices.size() - 1) + (anion_indices.size() - 1)
        + neutral_indices.size();
}

bool AmountChart::ionic() const {
    return !cation_indices.empty() || !anion_indices.empty();
}

AmountChart make_amount_chart(const std::vector<int>& charges) {
    if (charges.empty()) {
        throw std::invalid_argument("amount chart requires at least one species");
    }
    AmountChart chart;
    chart.charges = charges;
    for (std::size_t species = 0; species < charges.size(); ++species) {
        if (charges[species] > 0) {
            chart.cation_indices.push_back(species);
        } else if (charges[species] < 0) {
            chart.anion_indices.push_back(species);
        } else {
            chart.neutral_indices.push_back(species);
        }
    }
    if (chart.cation_indices.empty() != chart.anion_indices.empty()) {
        throw std::invalid_argument(
            "ionic amount chart requires both cations and anions"
        );
    }
    return chart;
}

namespace {

struct SimplexDerivatives {
    std::vector<double> shares;
    std::vector<double> jacobian;
    std::vector<double> hessians;
};

SimplexDerivatives reference_softmax(
    const std::vector<double>& coordinates,
    std::size_t offset,
    std::size_t category_count
) {
    if (category_count == 0) {
        throw std::invalid_argument("softmax category count must be positive");
    }
    const std::size_t dimension = category_count - 1;
    double maximum = 0.0;
    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
        maximum = std::max(maximum, coordinates[offset + coordinate]);
    }
    std::vector<double> weights(category_count, std::exp(-maximum));
    for (std::size_t category = 0; category < dimension; ++category) {
        weights[category] = std::exp(coordinates[offset + category] - maximum);
    }
    const double total = std::accumulate(weights.begin(), weights.end(), 0.0);
    SimplexDerivatives result;
    result.shares.resize(category_count, 0.0);
    std::transform(weights.begin(), weights.end(), result.shares.begin(), [total](double value) {
        return value / total;
    });
    if (!std::all_of(result.shares.begin(), result.shares.end(), [](double share) {
            return std::isfinite(share) && share > 0.0;
        })) {
        throw std::invalid_argument(
            "amount chart simplex is outside the strictly positive representable domain"
        );
    }
    result.jacobian.assign(category_count * dimension, 0.0);
    result.hessians.assign(category_count * dimension * dimension, 0.0);
    for (std::size_t category = 0; category < category_count; ++category) {
        for (std::size_t first = 0; first < dimension; ++first) {
            const double first_delta = category == first ? 1.0 : 0.0;
            const double first_factor = first_delta - result.shares[first];
            result.jacobian[category * dimension + first] =
                result.shares[category] * first_factor;
            for (std::size_t second = 0; second < dimension; ++second) {
                const double second_delta = category == second ? 1.0 : 0.0;
                const double simplex_delta = first == second ? 1.0 : 0.0;
                result.hessians[
                    category * dimension * dimension + first * dimension + second
                ] = result.shares[category]
                    * (
                        first_factor * (second_delta - result.shares[second])
                        - result.shares[first]
                            * (simplex_delta - result.shares[second])
                    );
            }
        }
    }
    return result;
}

void fill_charged_group(
    const AmountChart& chart,
    const std::vector<std::size_t>& species_indices,
    const SimplexDerivatives& simplex,
    std::size_t coordinate_offset,
    double charge_equivalents,
    AmountChartEvaluation& result
) {
    const std::size_t coordinate_count = chart.coordinate_count();
    const std::size_t simplex_dimension = species_indices.size() - 1;
    for (std::size_t category = 0; category < species_indices.size(); ++category) {
        const std::size_t species = species_indices[category];
        const double charge = std::abs(static_cast<double>(chart.charges[species]));
        const double amount = charge_equivalents * simplex.shares[category] / charge;
        if (!std::isfinite(amount) || amount <= 0.0) {
            throw std::invalid_argument(
                "amount chart is outside the strictly positive representable domain"
            );
        }
        result.amounts[species] = amount;
        result.jacobian[species * coordinate_count] = amount;
        result.amount_hessians[species * coordinate_count * coordinate_count] = amount;
        for (std::size_t first = 0; first < simplex_dimension; ++first) {
            const std::size_t first_coordinate = coordinate_offset + first;
            const double first_derivative = charge_equivalents
                * simplex.jacobian[category * simplex_dimension + first] / charge;
            result.jacobian[species * coordinate_count + first_coordinate] =
                first_derivative;
            result.amount_hessians[
                species * coordinate_count * coordinate_count + first_coordinate
            ] = first_derivative;
            result.amount_hessians[
                species * coordinate_count * coordinate_count
                    + first_coordinate * coordinate_count
            ] = first_derivative;
            for (std::size_t second = 0; second < simplex_dimension; ++second) {
                const std::size_t second_coordinate = coordinate_offset + second;
                result.amount_hessians[
                    species * coordinate_count * coordinate_count
                        + first_coordinate * coordinate_count + second_coordinate
                ] = charge_equivalents
                    * simplex.hessians[
                        category * simplex_dimension * simplex_dimension
                            + first * simplex_dimension + second
                    ] / charge;
            }
        }
    }
}

}  // namespace

AmountChartEvaluation evaluate_amount_chart(
    const AmountChart& chart,
    const std::vector<double>& coordinates
) {
    const std::size_t coordinate_count = chart.coordinate_count();
    if (coordinates.size() != coordinate_count) {
        throw std::invalid_argument("amount chart coordinate count does not match topology");
    }
    require_finite_vector(coordinates, "amount chart coordinates");
    AmountChartEvaluation result;
    result.amounts.assign(chart.charges.size(), 0.0);
    result.jacobian.assign(chart.charges.size() * coordinate_count, 0.0);
    result.amount_hessians.assign(
        chart.charges.size() * coordinate_count * coordinate_count,
        0.0
    );
    if (!chart.ionic()) {
        for (std::size_t species = 0; species < chart.charges.size(); ++species) {
            const double amount = std::exp(coordinates[species]);
            if (!std::isfinite(amount) || amount <= 0.0) {
                throw std::invalid_argument("amount chart produced a non-finite amount");
            }
            result.amounts[species] = amount;
            result.jacobian[species * coordinate_count + species] = amount;
            result.amount_hessians[
                species * coordinate_count * coordinate_count
                    + species * coordinate_count + species
            ] = amount;
        }
    } else {
        const double charge_equivalents = std::exp(coordinates[0]);
        if (!std::isfinite(charge_equivalents) || charge_equivalents <= 0.0) {
            throw std::invalid_argument("amount chart produced non-finite charge equivalents");
        }
        const std::size_t cation_offset = 1;
        const std::size_t anion_offset = cation_offset + chart.cation_indices.size() - 1;
        const std::size_t neutral_offset = anion_offset + chart.anion_indices.size() - 1;
        const SimplexDerivatives cations = reference_softmax(
            coordinates, cation_offset, chart.cation_indices.size()
        );
        const SimplexDerivatives anions = reference_softmax(
            coordinates, anion_offset, chart.anion_indices.size()
        );
        fill_charged_group(
            chart, chart.cation_indices, cations, cation_offset, charge_equivalents, result
        );
        fill_charged_group(
            chart, chart.anion_indices, anions, anion_offset, charge_equivalents, result
        );
        for (std::size_t neutral = 0; neutral < chart.neutral_indices.size(); ++neutral) {
            const std::size_t species = chart.neutral_indices[neutral];
            const std::size_t coordinate = neutral_offset + neutral;
            const double amount = std::exp(coordinates[coordinate]);
            if (!std::isfinite(amount) || amount <= 0.0) {
                throw std::invalid_argument("amount chart produced a non-finite neutral amount");
            }
            result.amounts[species] = amount;
            result.jacobian[species * coordinate_count + coordinate] = amount;
            result.amount_hessians[
                species * coordinate_count * coordinate_count
                    + coordinate * coordinate_count + coordinate
            ] = amount;
        }
    }
    result.minimum_amount = *std::min_element(result.amounts.begin(), result.amounts.end());
    for (std::size_t species = 0; species < chart.charges.size(); ++species) {
        result.charge_residual += static_cast<double>(chart.charges[species])
            * result.amounts[species];
    }
    return result;
}

std::vector<double> invert_amount_chart(
    const AmountChart& chart,
    const std::vector<double>& amounts
) {
    if (amounts.size() != chart.charges.size()
        || !std::all_of(amounts.begin(), amounts.end(), [](double value) {
            return std::isfinite(value) && value > 0.0;
        })) {
        throw std::invalid_argument("amount chart inverse requires finite positive amounts");
    }
    std::vector<double> coordinates(chart.coordinate_count(), 0.0);
    if (!chart.ionic()) {
        std::transform(amounts.begin(), amounts.end(), coordinates.begin(), [](double value) {
            return std::log(value);
        });
        return coordinates;
    }
    double positive_equivalents = 0.0;
    for (std::size_t species : chart.cation_indices) {
        positive_equivalents += static_cast<double>(chart.charges[species]) * amounts[species];
    }
    double negative_equivalents = 0.0;
    for (std::size_t species : chart.anion_indices) {
        negative_equivalents += std::abs(static_cast<double>(chart.charges[species]))
            * amounts[species];
    }
    const double charge_tolerance = numerical_tolerance(
        std::max(positive_equivalents, negative_equivalents), chart.charges.size()
    );
    if (std::abs(positive_equivalents - negative_equivalents) > charge_tolerance) {
        throw std::invalid_argument("amount chart inverse requires an electroneutral state");
    }
    const double charge_equivalents = 0.5 * (positive_equivalents + negative_equivalents);
    coordinates[0] = std::log(charge_equivalents);
    std::size_t offset = 1;
    for (const std::vector<std::size_t>* group : {&chart.cation_indices, &chart.anion_indices}) {
        const std::size_t reference_species = group->back();
        const double reference_share = std::abs(static_cast<double>(chart.charges[reference_species]))
            * amounts[reference_species] / charge_equivalents;
        for (std::size_t category = 0; category + 1 < group->size(); ++category) {
            const std::size_t species = (*group)[category];
            const double share = std::abs(static_cast<double>(chart.charges[species]))
                * amounts[species] / charge_equivalents;
            coordinates[offset++] = std::log(share / reference_share);
        }
    }
    for (std::size_t species : chart.neutral_indices) {
        coordinates[offset++] = std::log(amounts[species]);
    }
    return coordinates;
}

CompiledReactionSystem compile_reaction_system(const ReactionSystemInput& input) {
    validate_identities(input);
    require_finite(input.temperature_k, "temperature");
    require_finite(input.pressure_pa, "pressure");
    if (input.temperature_k <= 0.0 || input.pressure_pa <= 0.0) {
        throw std::invalid_argument("temperature and pressure must be positive");
    }
    const std::size_t species_count = input.species_ids.size();
    if (input.reaction_matrix.rows == 0) {
        throw std::invalid_argument(
            "complete closed system requires independent reactions and lnK"
        );
    }
    if (input.feed_amounts.size() != species_count
        || input.molar_masses_kg_per_mol.size() != species_count
        || input.balance_matrix.columns != species_count
        || input.reaction_matrix.columns != species_count) {
        throw std::invalid_argument("reaction-system dimensions do not match species count");
    }
    require_finite_vector(input.feed_amounts, "feed amounts");
    if (!std::all_of(input.feed_amounts.begin(), input.feed_amounts.end(), [](double value) {
            return value >= 0.0;
        })) {
        throw std::invalid_argument("feed amounts must be nonnegative");
    }
    if (std::accumulate(input.feed_amounts.begin(), input.feed_amounts.end(), 0.0) <= 0.0) {
        throw std::invalid_argument("feed must contain a positive total amount");
    }
    require_finite_vector(input.balance_matrix.values, "balance matrix");
    require_finite_vector(input.reaction_matrix.values, "reaction matrix");
    require_finite_vector(input.ln_k, "lnK");
    require_finite_vector(input.molar_masses_kg_per_mol, "molar masses");
    if (!std::all_of(
            input.molar_masses_kg_per_mol.begin(),
            input.molar_masses_kg_per_mol.end(),
            [](double value) { return value > 0.0; }
        )) {
        throw std::invalid_argument("molar masses must be finite and positive");
    }
    if (input.ln_k.size() != input.reaction_matrix.rows) {
        throw std::invalid_argument("equilibrium constants are incomplete for the reaction matrix");
    }
    validate_reference_records(input);

    const double conservation_norm = reaction_conservation_inf_norm(
        input.molar_masses_kg_per_mol,
        input.balance_matrix,
        input.reaction_matrix
    );
    const double molar_mass_scale = *std::max_element(
        input.molar_masses_kg_per_mol.begin(),
        input.molar_masses_kg_per_mol.end()
    );
    const double conservation_tolerance = numerical_tolerance(
        std::max(molar_mass_scale, matrix_scale(input.balance_matrix))
            * matrix_scale(input.reaction_matrix),
        species_count
    );
    if (conservation_norm > conservation_tolerance) {
        throw std::invalid_argument(
            "reaction stoichiometry does not conserve molar mass and the balance matrix"
        );
    }
    double feed_charge = 0.0;
    for (std::size_t species = 0; species < species_count; ++species) {
        feed_charge += static_cast<double>(input.charges[species]) * input.feed_amounts[species];
    }
    if (std::abs(feed_charge) > numerical_tolerance(
            *std::max_element(input.feed_amounts.begin(), input.feed_amounts.end()), species_count
        )) {
        throw std::invalid_argument("feed must be exactly electroneutral within numerical precision");
    }
    const double reaction_charge_norm = charge_reaction_inf_norm(
        input.charges, input.reaction_matrix
    );
    if (reaction_charge_norm > numerical_tolerance(matrix_scale(input.reaction_matrix), species_count)) {
        throw std::invalid_argument("reaction stoichiometry does not conserve charge");
    }

    const std::vector<std::size_t> reaction_basis_rows = select_independent_rows(
        input.reaction_matrix
    );
    if (reaction_basis_rows.empty()) {
        throw std::invalid_argument("reaction matrix rank is deficient");
    }
    const DenseMatrix reaction_matrix = selected_rows(
        input.reaction_matrix, reaction_basis_rows
    );
    const RowQr reaction_factor = factor_row_basis(reaction_matrix);
    DenseMatrix reaction_transform{
        input.reaction_matrix.rows,
        reaction_matrix.rows,
        std::vector<double>(
            input.reaction_matrix.rows * reaction_matrix.rows, 0.0
        ),
    };
    for (std::size_t supplied = 0; supplied < input.reaction_matrix.rows; ++supplied) {
        const auto basis_position = std::find(
            reaction_basis_rows.begin(), reaction_basis_rows.end(), supplied
        );
        if (basis_position != reaction_basis_rows.end()) {
            reaction_transform(
                supplied,
                static_cast<std::size_t>(
                    std::distance(reaction_basis_rows.begin(), basis_position)
                )
            ) = 1.0;
            continue;
        }
        const std::vector<double> coordinates = row_coordinates(
            matrix_row(input.reaction_matrix, supplied), reaction_factor
        );
        for (std::size_t basis = 0; basis < coordinates.size(); ++basis) {
            reaction_transform(supplied, basis) = coordinates[basis];
        }
    }
    const double reaction_transform_norm = matrix_reconstruction_inf_norm(
        reaction_transform, reaction_matrix, input.reaction_matrix
    );
    if (reaction_transform_norm > numerical_tolerance(
            matrix_scale(input.reaction_matrix), species_count
        )) {
        throw std::invalid_argument("reaction matrix rank reduction is inconsistent");
    }
    std::vector<double> independent_ln_k;
    independent_ln_k.reserve(reaction_basis_rows.size());
    for (std::size_t row : reaction_basis_rows) {
        independent_ln_k.push_back(input.ln_k[row]);
    }
    const double reaction_cycle_norm = constant_reconstruction_inf_norm(
        reaction_transform, independent_ln_k, input.ln_k
    );
    double ln_k_scale = 0.0;
    for (double value : input.ln_k) {
        ln_k_scale = std::max(ln_k_scale, std::abs(value));
    }
    if (reaction_cycle_norm > numerical_tolerance(
            ln_k_scale, std::max(input.reaction_matrix.rows, species_count)
        )) {
        throw std::invalid_argument(
            "reaction constant cycle is inconsistent after Provider-basis conversion"
        );
    }

    bool charge_is_independent = false;
    const DenseMatrix balance_matrix = matrix_from_rows(
        seeded_balance_rows(input, charge_is_independent)
    );
    const std::size_t balance_rank = balance_matrix.rows;
    const std::size_t reaction_rank = reaction_matrix.rows;
    if (balance_rank + (charge_is_independent ? 1 : 0) + reaction_rank
        != species_count) {
        throw std::invalid_argument(
            "complete closed system reaction matrix rank sum does not equal species count"
        );
    }

    ReferenceReconstruction reconstruction = construct_minimum_norm_reference(
        reaction_matrix, independent_ln_k
    );

    const double reference_residual = reference_residual_inf_norm(
        reaction_matrix, reconstruction.reference, independent_ln_k
    );
    if (reference_residual > numerical_tolerance(ln_k_scale, species_count)) {
        throw std::invalid_argument("standard chemical reference reconstruction is inconsistent");
    }

    std::vector<double> balance_totals(balance_matrix.rows, 0.0);
    for (std::size_t row = 0; row < balance_matrix.rows; ++row) {
        for (std::size_t species = 0; species < species_count; ++species) {
            balance_totals[row] += balance_matrix(row, species)
                * input.feed_amounts[species];
        }
    }
    CompiledReactionSystem result;
    result.species_ids = input.species_ids;
    result.charges = input.charges;
    result.molar_masses_kg_per_mol = input.molar_masses_kg_per_mol;
    result.supplied_balance_matrix = input.balance_matrix;
    result.supplied_reaction_matrix = input.reaction_matrix;
    result.supplied_ln_k = input.ln_k;
    result.reaction_transform = std::move(reaction_transform);
    result.reaction_basis_rows = reaction_basis_rows;
    result.reaction_cycle_inf_norm = reaction_cycle_norm;
    result.reaction_transform_inf_norm = reaction_transform_norm;
    result.balance_matrix = balance_matrix;
    result.reaction_matrix = reaction_matrix;
    result.balance_totals = std::move(balance_totals);
    result.feed_amounts = input.feed_amounts;
    result.ln_k = std::move(independent_ln_k);
    result.g_ref = std::move(reconstruction.reference);
    result.provider_fingerprint = input.provider_fingerprint;
    result.balance_rank = balance_rank;
    result.reaction_rank = reaction_rank;
    result.reaction_qr_diagonal_ratio = reconstruction.qr_diagonal_ratio;
    result.reference_reconstruction_inf_norm = reference_residual;
    result.conservation_reaction_inf_norm = conservation_norm;
    result.charge_reaction_inf_norm = reaction_charge_norm;
    return result;
}

}  // namespace epcsaft_equilibrium
