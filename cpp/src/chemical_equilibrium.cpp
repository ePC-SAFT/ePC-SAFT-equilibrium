#include "chemical_equilibrium.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace epcsaft_equilibrium {
namespace {

constexpr double kResidualMultiplier = 4096.0;

void require_finite_vector(const std::vector<double>& values, const char* field) {
    if (!std::all_of(values.begin(), values.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::invalid_argument(std::string(field) + " must be finite");
    }
}

double numerical_tolerance(double scale, std::size_t dimension) {
    return kResidualMultiplier * std::numeric_limits<double>::epsilon()
        * std::max(1.0, scale) * static_cast<double>(std::max<std::size_t>(1, dimension));
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

}  // namespace epcsaft_equilibrium
