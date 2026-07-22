#include "chemical_equilibrium.hpp"
#include "provider.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
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

std::size_t matrix_rank(DenseMatrix matrix) {
    if (matrix.values.size() != matrix.rows * matrix.columns) {
        throw std::invalid_argument("matrix dimensions do not match storage");
    }
    const double tolerance = numerical_tolerance(
        matrix_scale(matrix), std::max(matrix.rows, matrix.columns)
    );
    std::size_t pivot_row = 0;
    for (std::size_t column = 0; column < matrix.columns && pivot_row < matrix.rows; ++column) {
        std::size_t pivot = pivot_row;
        for (std::size_t row = pivot_row + 1; row < matrix.rows; ++row) {
            if (std::abs(matrix(row, column)) > std::abs(matrix(pivot, column))) {
                pivot = row;
            }
        }
        if (std::abs(matrix(pivot, column)) <= tolerance) {
            continue;
        }
        if (pivot != pivot_row) {
            for (std::size_t index = 0; index < matrix.columns; ++index) {
                std::swap(matrix(pivot, index), matrix(pivot_row, index));
            }
        }
        const double pivot_value = matrix(pivot_row, column);
        for (std::size_t row = pivot_row + 1; row < matrix.rows; ++row) {
            const double factor = matrix(row, column) / pivot_value;
            for (std::size_t index = column; index < matrix.columns; ++index) {
                matrix(row, index) -= factor * matrix(pivot_row, index);
            }
        }
        ++pivot_row;
    }
    return pivot_row;
}

double multiply_inf_norm(const DenseMatrix& left, const DenseMatrix& right_transposed) {
    if (left.columns != right_transposed.columns) {
        throw std::invalid_argument("matrix product dimensions are inconsistent");
    }
    double norm = 0.0;
    for (std::size_t row = 0; row < left.rows; ++row) {
        for (std::size_t reaction = 0; reaction < right_transposed.rows; ++reaction) {
            double value = 0.0;
            for (std::size_t species = 0; species < left.columns; ++species) {
                value += left(row, species) * right_transposed(reaction, species);
            }
            norm = std::max(norm, std::abs(value));
        }
    }
    return norm;
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
    if (input.provider_component_ids != input.expected_provider_component_ids
        || input.provider_component_ids != input.species_ids) {
        throw std::invalid_argument("Provider component order does not match true species order");
    }
    if (input.charges.size() != species_count
        || input.provider_charges != input.expected_provider_charges
        || input.provider_charges != input.charges) {
        throw std::invalid_argument("Provider charges do not match compiled species charges");
    }
    if (input.provider_fingerprint.empty()
        || input.provider_fingerprint != input.expected_provider_fingerprint) {
        throw std::invalid_argument("Provider fingerprint does not match the expected artifact");
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
        if (!record.dimensionless) {
            throw std::invalid_argument("lnK must be explicitly dimensionless");
        }
        require_finite(record.temperature_k, "equilibrium constant temperature");
        require_finite(record.pressure_pa, "equilibrium constant pressure");
        if (record.temperature_k != input.temperature_k) {
            throw std::invalid_argument("equilibrium constant temperature binding does not match");
        }
        if (record.pressure_binding == "fixed") {
            if (record.pressure_pa != input.pressure_pa) {
                throw std::invalid_argument("equilibrium constant pressure binding does not match");
            }
        } else if (record.pressure_binding != "temperature_only") {
            throw std::invalid_argument("equilibrium constant pressure binding is incomplete");
        }
    }
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
    if (input.balance_matrix.rows != input.declared_balance_rank) {
        throw std::invalid_argument(
            "declared balance rank does not match balance matrix rank or independent row count"
        );
    }
    const std::size_t balance_rank = matrix_rank(input.balance_matrix);
    if (balance_rank != input.declared_balance_rank) {
        throw std::invalid_argument("declared balance rank does not match balance matrix rank");
    }
    if (input.ln_k.size() != input.reaction_matrix.rows) {
        throw std::invalid_argument("equilibrium constants are incomplete for the reaction matrix");
    }
    ReferenceReconstruction reconstruction = construct_minimum_norm_reference(
        input.reaction_matrix, input.ln_k
    );
    const std::size_t reaction_rank = reconstruction.rank;
    if (input.complete_closed_system
        && balance_rank + reaction_rank != species_count) {
        throw std::invalid_argument("complete closed system rank sum does not equal species count");
    }

    const double conservation_norm = multiply_inf_norm(
        input.balance_matrix, input.reaction_matrix
    );
    const double conservation_tolerance = numerical_tolerance(
        matrix_scale(input.balance_matrix) * matrix_scale(input.reaction_matrix),
        species_count
    );
    if (conservation_norm > conservation_tolerance) {
        throw std::invalid_argument("reaction stoichiometry does not conserve the balance matrix");
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
    validate_reference_records(input);

    const double reference_residual = reference_residual_inf_norm(
        input.reaction_matrix, reconstruction.reference, input.ln_k
    );
    double ln_k_scale = 0.0;
    for (double value : input.ln_k) {
        ln_k_scale = std::max(ln_k_scale, std::abs(value));
    }
    if (reference_residual > numerical_tolerance(ln_k_scale, species_count)) {
        throw std::invalid_argument("standard chemical reference reconstruction is inconsistent");
    }

    std::vector<double> balance_totals(input.balance_matrix.rows, 0.0);
    for (std::size_t row = 0; row < input.balance_matrix.rows; ++row) {
        for (std::size_t species = 0; species < species_count; ++species) {
            balance_totals[row] += input.balance_matrix(row, species)
                * input.feed_amounts[species];
        }
    }
    return {
        input.species_ids,
        input.charges,
        input.balance_matrix,
        input.reaction_matrix,
        std::move(balance_totals),
        input.feed_amounts,
        input.ln_k,
        std::move(reconstruction.reference),
        input.provider_fingerprint,
        balance_rank,
        reaction_rank,
        reconstruction.qr_diagonal_ratio,
        reference_residual,
        conservation_norm,
        reaction_charge_norm,
    };
}

SourceStandardStateResult transform_source_standard_state(
    const DenseMatrix& reaction_matrix,
    const std::vector<double>& source_ln_k,
    const std::vector<double>& log_activity_scale_factors,
    const std::vector<int>& charges,
    const std::vector<std::string>& component_ids,
    const std::string& provider_fingerprint,
    double temperature_k,
    double pressure_pa,
    const NeutralReferenceEvaluation& reference
) {
    constexpr std::string_view kProviderBasis =
        EPCSAFT_NATIVE_HELMHOLTZ_BASIS_ID_V1;
    if (reaction_matrix.rows == 0 || reaction_matrix.columns == 0
        || reaction_matrix.values.size() != reaction_matrix.rows * reaction_matrix.columns
        || source_ln_k.size() != reaction_matrix.rows
        || log_activity_scale_factors.size() != reaction_matrix.columns
        || charges.size() != reaction_matrix.columns
        || component_ids.size() != reaction_matrix.columns
        || reference.component_count != reaction_matrix.columns
        || reference.neutral_basis_row_count == 0
        || reference.neutral_basis_row_count > reaction_matrix.columns
        || reference.neutral_basis.size()
            != reference.neutral_basis_row_count * reaction_matrix.columns
        || reference.log_fugacity_contractions.size()
            != reference.neutral_basis_row_count
        || reference.reference_composition.size() != reaction_matrix.columns) {
        throw std::invalid_argument("neutral-reference transformation dimensions are inconsistent");
    }
    if (reference.basis_id != kProviderBasis
        || reference.parameter_fingerprint.empty()
        || provider_fingerprint.empty()
        || provider_fingerprint != reference.parameter_fingerprint) {
        throw std::invalid_argument("Provider neutral-reference identity is incompatible");
    }
    if (reference.derivative_availability
        != EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_NONE_V1) {
        throw std::invalid_argument(
            "Provider neutral-reference derivative availability is unsupported"
        );
    }
    if (!std::isfinite(temperature_k) || !std::isfinite(pressure_pa)
        || temperature_k != reference.temperature_k || pressure_pa != reference.pressure_pa) {
        throw std::invalid_argument("source/reference temperature and pressure are not bound");
    }
    std::unordered_set<std::string> unique_ids;
    for (const std::string& id : component_ids) {
        if (id.empty() || !unique_ids.insert(id).second) {
            throw std::invalid_argument("source component identity is incomplete");
        }
    }
    require_finite_vector(source_ln_k, "source lnK");
    require_finite_vector(log_activity_scale_factors, "log activity-scale factors");
    require_finite_vector(reaction_matrix.values, "source reactions");
    require_finite_vector(reference.reference_composition, "Provider reference composition");
    if (reference.reference_amount_mol != 1.0
        || reference.reference_number_density_mol_per_m3 != 1.0
        || !std::isfinite(reference.solvent_molar_mass_kg_per_mol)
        || reference.solvent_molar_mass_kg_per_mol <= 0.0
        || !std::isfinite(reference.reference_molality_mol_per_kg)
        || reference.reference_molality_mol_per_kg <= 0.0
        || !std::isfinite(reference.reference_convergence_error)
        || reference.reference_convergence_error < 0.0
        || reference.reference_convergence_error > 5.0e-5) {
        throw std::invalid_argument("Provider neutral-reference scalar identity is incompatible");
    }
    double reference_composition_sum = 0.0;
    for (double value : reference.reference_composition) {
        if (value < 0.0) {
            throw std::invalid_argument("Provider reference composition is negative");
        }
        reference_composition_sum += value;
    }
    if (std::abs(reference_composition_sum - 1.0)
        > numerical_tolerance(1.0, reference.reference_composition.size())) {
        throw std::invalid_argument("Provider reference composition is not normalized");
    }
    for (std::size_t row = 0; row < reference.neutral_basis_row_count; ++row) {
        for (std::size_t component = 0; component < reaction_matrix.columns; ++component) {
            require_finite(
                reference.neutral_basis[row * reaction_matrix.columns + component],
                "Provider neutral basis"
            );
        }
    }
    require_finite_vector(reference.log_fugacity_contractions, "Provider contractions");
    DenseMatrix neutral_basis{
        reference.neutral_basis_row_count,
        reaction_matrix.columns,
        reference.neutral_basis,
    };
    if (matrix_rank(neutral_basis) != neutral_basis.rows) {
        throw std::invalid_argument("Provider neutral-reference basis is rank deficient");
    }
    const double charge_tolerance = numerical_tolerance(
        matrix_scale(reaction_matrix), reaction_matrix.columns
    );
    for (std::size_t row = 0; row < reaction_matrix.rows; ++row) {
        double reaction_charge = 0.0;
        for (std::size_t component = 0; component < reaction_matrix.columns; ++component) {
            reaction_charge += static_cast<double>(charges[component])
                * reaction_matrix(row, component);
        }
        if (std::abs(reaction_charge) > charge_tolerance) {
            throw std::invalid_argument("source reaction is not charge neutral");
        }
    }
    for (std::size_t row = 0; row < neutral_basis.rows; ++row) {
        double basis_charge = 0.0;
        for (std::size_t component = 0; component < neutral_basis.columns; ++component) {
            basis_charge += static_cast<double>(charges[component])
                * neutral_basis(row, component);
        }
        if (std::abs(basis_charge) > charge_tolerance) {
            throw std::invalid_argument("Provider neutral basis row is not charge neutral");
        }
    }

    const std::size_t basis_count = neutral_basis.rows;
    DenseMatrix gram{basis_count, basis_count, std::vector<double>(basis_count * basis_count, 0.0)};
    for (std::size_t first = 0; first < basis_count; ++first) {
        for (std::size_t second = 0; second < basis_count; ++second) {
            for (std::size_t component = 0; component < neutral_basis.columns; ++component) {
                gram(first, second) += neutral_basis(first, component)
                    * neutral_basis(second, component);
            }
        }
    }
    auto solve_full_column_system = [&](const std::vector<double>& values) {
        std::vector<double> matrix = gram.values;
        std::vector<double> rhs(basis_count, 0.0);
        for (std::size_t row = 0; row < basis_count; ++row) {
            for (std::size_t component = 0; component < neutral_basis.columns; ++component) {
                rhs[row] += neutral_basis(row, component) * values[component];
            }
        }
        for (std::size_t pivot = 0; pivot < basis_count; ++pivot) {
            std::size_t selected = pivot;
            for (std::size_t row = pivot + 1; row < basis_count; ++row) {
                if (std::abs(matrix[row * basis_count + pivot])
                    > std::abs(matrix[selected * basis_count + pivot])) {
                    selected = row;
                }
            }
            if (std::abs(matrix[selected * basis_count + pivot])
                <= numerical_tolerance(1.0, basis_count)) {
                throw std::invalid_argument("Provider neutral-reference basis solve is rank deficient");
            }
            if (selected != pivot) {
                for (std::size_t column = pivot; column < basis_count; ++column) {
                    std::swap(
                        matrix[selected * basis_count + column],
                        matrix[pivot * basis_count + column]
                    );
                }
                std::swap(rhs[selected], rhs[pivot]);
            }
            for (std::size_t row = pivot + 1; row < basis_count; ++row) {
                const double factor = matrix[row * basis_count + pivot]
                    / matrix[pivot * basis_count + pivot];
                for (std::size_t column = pivot; column < basis_count; ++column) {
                    matrix[row * basis_count + column] -= factor
                        * matrix[pivot * basis_count + column];
                }
                rhs[row] -= factor * rhs[pivot];
            }
        }
        std::vector<double> result(basis_count, 0.0);
        for (std::size_t reverse = basis_count; reverse > 0; --reverse) {
            const std::size_t row = reverse - 1;
            double value = rhs[row];
            for (std::size_t column = row + 1; column < basis_count; ++column) {
                value -= matrix[row * basis_count + column] * result[column];
            }
            result[row] = value / matrix[row * basis_count + row];
        }
        return result;
    };

    SourceStandardStateResult result;
    result.standard_offsets.assign(reaction_matrix.rows, 0.0);
    result.ln_k_provider_basis.assign(reaction_matrix.rows, 0.0);
    result.reaction_to_neutral_basis = DenseMatrix{
        reaction_matrix.rows,
        basis_count,
        std::vector<double>(reaction_matrix.rows * basis_count, 0.0),
    };
    for (std::size_t reaction = 0; reaction < reaction_matrix.rows; ++reaction) {
        std::vector<double> nu(reaction_matrix.columns, 0.0);
        for (std::size_t component = 0; component < reaction_matrix.columns; ++component) {
            nu[component] = reaction_matrix(reaction, component);
        }
        const std::vector<double> coefficients = solve_full_column_system(nu);
        double residual = 0.0;
        double scale_offset = 0.0;
        double contraction_offset = 0.0;
        for (std::size_t component = 0; component < reaction_matrix.columns; ++component) {
            scale_offset += nu[component] * log_activity_scale_factors[component];
        }
        for (std::size_t basis = 0; basis < basis_count; ++basis) {
            result.reaction_to_neutral_basis(reaction, basis) = coefficients[basis];
            contraction_offset += coefficients[basis]
                * reference.log_fugacity_contractions[basis];
        }
        for (std::size_t component = 0; component < reaction_matrix.columns; ++component) {
            double reconstructed = 0.0;
            for (std::size_t basis = 0; basis < basis_count; ++basis) {
                reconstructed += coefficients[basis] * neutral_basis(basis, component);
            }
            residual = std::max(residual, std::abs(reconstructed - nu[component]));
        }
        result.representation_residual_inf_norm = std::max(
            result.representation_residual_inf_norm, residual
        );
        result.standard_offsets[reaction] = scale_offset + contraction_offset;
        result.ln_k_provider_basis[reaction] = source_ln_k[reaction]
            + result.standard_offsets[reaction];
    }
    result.derivative_availability = reference.derivative_availability;
    result.basis_id = reference.basis_id;
    result.parameter_fingerprint = reference.parameter_fingerprint;
    return result;
}

}  // namespace epcsaft_equilibrium
