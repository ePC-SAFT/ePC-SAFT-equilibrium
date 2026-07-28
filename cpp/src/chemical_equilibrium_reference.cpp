#include "chemical_equilibrium_reference.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace epcsaft_equilibrium {
namespace {

constexpr double kResidualMultiplier = 4096.0;

double numerical_tolerance(double scale, std::size_t dimension) {
    return kResidualMultiplier * std::numeric_limits<double>::epsilon()
        * std::max(1.0, scale) * static_cast<double>(std::max<std::size_t>(1, dimension));
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

}  // namespace

ReactionReferenceReconstruction reconstruct_reaction_reference(
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
    const double residual = reference_residual_inf_norm(
        reaction_matrix, reference, ln_k
    );
    return {
        std::move(reference),
        reaction_count,
        *minimum_diagonal / *maximum_diagonal,
        residual,
    };
}

}  // namespace epcsaft_equilibrium
