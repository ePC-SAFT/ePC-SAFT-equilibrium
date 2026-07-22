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

std::vector<double> solve_square(DenseMatrix matrix, std::vector<double> right_hand_side) {
    if (matrix.rows == 0 || matrix.rows != matrix.columns
        || matrix.values.size() != matrix.rows * matrix.columns
        || right_hand_side.size() != matrix.rows) {
        throw std::invalid_argument("reference reconstruction linear system is invalid");
    }
    const double tolerance = numerical_tolerance(matrix_scale(matrix), matrix.rows);
    for (std::size_t column = 0; column < matrix.columns; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < matrix.rows; ++row) {
            if (std::abs(matrix(row, column)) > std::abs(matrix(pivot, column))) {
                pivot = row;
            }
        }
        if (std::abs(matrix(pivot, column)) <= tolerance) {
            throw std::invalid_argument("reaction matrix rank is deficient");
        }
        if (pivot != column) {
            for (std::size_t index = column; index < matrix.columns; ++index) {
                std::swap(matrix(pivot, index), matrix(column, index));
            }
            std::swap(right_hand_side[pivot], right_hand_side[column]);
        }
        const double pivot_value = matrix(column, column);
        for (std::size_t row = column + 1; row < matrix.rows; ++row) {
            const double factor = matrix(row, column) / pivot_value;
            for (std::size_t index = column; index < matrix.columns; ++index) {
                matrix(row, index) -= factor * matrix(column, index);
            }
            right_hand_side[row] -= factor * right_hand_side[column];
        }
    }
    std::vector<double> solution(matrix.rows, 0.0);
    for (std::size_t reverse = matrix.rows; reverse > 0; --reverse) {
        const std::size_t row = reverse - 1;
        double value = right_hand_side[row];
        for (std::size_t column = row + 1; column < matrix.columns; ++column) {
            value -= matrix(row, column) * solution[column];
        }
        solution[row] = value / matrix(row, row);
    }
    return solution;
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

std::vector<double> construct_minimum_norm_reference(
    const DenseMatrix& reaction_matrix,
    const std::vector<double>& ln_k
) {
    DenseMatrix gram{
        reaction_matrix.rows,
        reaction_matrix.rows,
        std::vector<double>(reaction_matrix.rows * reaction_matrix.rows, 0.0),
    };
    for (std::size_t row = 0; row < reaction_matrix.rows; ++row) {
        for (std::size_t column = 0; column < reaction_matrix.rows; ++column) {
            for (std::size_t species = 0; species < reaction_matrix.columns; ++species) {
                gram(row, column) += reaction_matrix(row, species)
                    * reaction_matrix(column, species);
            }
        }
    }
    std::vector<double> right_hand_side(ln_k.size(), 0.0);
    std::transform(ln_k.begin(), ln_k.end(), right_hand_side.begin(), [](double value) {
        return -value;
    });
    const std::vector<double> dual = solve_square(std::move(gram), std::move(right_hand_side));
    std::vector<double> reference(reaction_matrix.columns, 0.0);
    for (std::size_t species = 0; species < reaction_matrix.columns; ++species) {
        for (std::size_t reaction = 0; reaction < reaction_matrix.rows; ++reaction) {
            reference[species] += reaction_matrix(reaction, species) * dual[reaction];
        }
    }
    return reference;
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
    if (!input.has_standard_state_transformation) {
        if (!input.provider_basis_id.empty()
            || input.standard_state_transformation_residual != 0.0) {
            throw std::invalid_argument(
                "standard-state transformation evidence is inconsistent"
            );
        }
        return;
    }
    if (input.provider_basis_id.empty()) {
        throw std::invalid_argument("standard-state Provider basis identity is incomplete");
    }
    require_finite(
        input.standard_state_transformation_residual,
        "standard-state transformation residual"
    );
    double ln_k_scale = 0.0;
    for (double value : input.ln_k) {
        ln_k_scale = std::max(ln_k_scale, std::abs(value));
    }
    const double residual_tolerance = numerical_tolerance(
        ln_k_scale, input.species_ids.size()
    );
    if (std::abs(input.standard_state_transformation_residual) > residual_tolerance) {
        throw std::invalid_argument("standard-state transformation residual is inconsistent");
    }
    for (const EquilibriumConstantRecord& record : input.equilibrium_constant_records) {
        if (record.reference_id != input.provider_basis_id) {
            throw std::invalid_argument(
                "equilibrium constant Provider basis identity does not match transformation"
            );
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
    const std::size_t reaction_rank = matrix_rank(input.reaction_matrix);
    if (reaction_rank != input.reaction_matrix.rows) {
        throw std::invalid_argument("reaction matrix rank is deficient");
    }
    if (input.ln_k.size() != input.reaction_matrix.rows) {
        throw std::invalid_argument("equilibrium constants are incomplete for the reaction matrix");
    }
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

    std::vector<double> reference = construct_minimum_norm_reference(
        input.reaction_matrix, input.ln_k
    );
    const double reference_residual = reference_residual_inf_norm(
        input.reaction_matrix, reference, input.ln_k
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
        std::move(reference),
        input.provider_fingerprint,
        input.has_standard_state_transformation,
        input.provider_basis_id,
        input.standard_state_transformation_residual,
        balance_rank,
        reaction_rank,
        reference_residual,
        conservation_norm,
        reaction_charge_norm,
    };
}

MixedStandardStateResult transform_water_self_ionization_standard_state(
    double ln_kw_mixed_standard,
    const StandardReferenceEvaluation& reference
) {
    constexpr std::string_view kProviderBasis =
        "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3";
    if (!std::isfinite(ln_kw_mixed_standard)) {
        throw std::invalid_argument("mixed-standard lnKw must be finite");
    }
    if (reference.basis_id != kProviderBasis) {
        throw std::invalid_argument("Provider standard-reference basis identity mismatch");
    }
    if (reference.parameter_fingerprint.empty()) {
        throw std::invalid_argument("Provider standard-reference fingerprint is incomplete");
    }
    for (double value : {
             reference.formula_unit_log_fugacity,
             reference.pure_solvent_log_fugacity,
             reference.solvent_molar_mass_kg_per_mol,
             reference.reference_molality_mol_per_kg,
             reference.convergence_error,
             reference.pure_solvent_molar_volume_m3_per_mol,
         }) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Provider standard-reference scalar is not finite");
        }
    }
    if (reference.solvent_molar_mass_kg_per_mol <= 0.0
        || reference.solvent_molar_mass_kg_per_mol >= 1.0) {
        throw std::invalid_argument("Provider solvent molar mass must be in kg/mol");
    }
    if (reference.reference_molality_mol_per_kg <= 0.0) {
        throw std::invalid_argument("Provider terminal molality must be positive");
    }
    if (reference.convergence_error < 0.0) {
        throw std::invalid_argument("Provider reference convergence error must be nonnegative");
    }
    if (reference.pure_solvent_molar_volume_m3_per_mol <= 0.0) {
        throw std::invalid_argument("Provider pure-solvent molar volume must be positive");
    }

    // This is the IAPWS thermodynamic standard molality, not Provider's
    // terminal infinite-dilution evaluation coordinate.
    constexpr double kStandardMolalityMolPerKg = 1.0;
    const double delta = 2.0 * std::log(
        reference.solvent_molar_mass_kg_per_mol * kStandardMolalityMolPerKg
    ) + reference.formula_unit_log_fugacity
        - 2.0 * reference.pure_solvent_log_fugacity;
    return {delta, ln_kw_mixed_standard + delta};
}

}  // namespace epcsaft_equilibrium
