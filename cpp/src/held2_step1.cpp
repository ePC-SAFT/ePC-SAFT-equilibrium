#include "held2.hpp"
#include "held2_tolerances.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace epcsaft_equilibrium {
namespace {

void require_finite_vector(
    const std::vector<double>& values,
    const char* name
) {
    if (!std::all_of(values.begin(), values.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
}

double charge_residual(
    const std::vector<double>& charges,
    const std::vector<double>& fractions
) {
    return std::inner_product(
        charges.begin(), charges.end(), fractions.begin(), 0.0
    );
}

double charge_scale(
    const std::vector<double>& charges,
    const std::vector<double>& fractions
) {
    double scale = 0.0;
    for (std::size_t index = 0; index < charges.size(); ++index) {
        scale += std::abs(charges[index] * fractions[index]);
    }
    return scale;
}

}  // namespace

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
        result.modified_factors.push_back(
            1.0 - charges[index] / charges[eliminated]
        );
        if (charges[index] == 0.0) {
            result.dependent_index = index;
        }
    }
    for (std::size_t retained = 0;
         retained < result.retained_indices.size();
         ++retained) {
        const std::size_t index = result.retained_indices[retained];
        if (index == result.dependent_index) {
            continue;
        }
        result.independent_indices.push_back(index);
        const double factor = result.modified_factors[retained];
        result.independent_lower_bounds.push_back(
            kHeld2ModifiedLowerScale * factor
        );
        if (charges[index] == 0.0) {
            result.independent_upper_bounds.push_back(1.0);
            continue;
        }
        double largest_opposite = 0.0;
        for (double other_charge : charges) {
            if (other_charge * charges[index] < 0.0) {
                largest_opposite = std::max(
                    largest_opposite, std::abs(other_charge)
                );
            }
        }
        if (largest_opposite == 0.0) {
            throw std::invalid_argument(
                "charged species has no opposite-sign counterion"
            );
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
        throw std::invalid_argument(
            "physical composition size does not match charges"
        );
    }
    require_finite_vector(physical_fractions, "physical mole fractions");
    if (!std::all_of(
            physical_fractions.begin(),
            physical_fractions.end(),
            [](double value) { return value >= 0.0; }
        )) {
        throw std::invalid_argument(
            "physical mole fractions must be nonnegative"
        );
    }
    const double total = std::accumulate(
        physical_fractions.begin(), physical_fractions.end(), 0.0
    );
    if (!audit_held2_tolerance(kHeld2CompositionSum, total - 1.0).passed) {
        throw std::invalid_argument(
            "physical mole fractions must sum to one"
        );
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
    for (std::size_t retained = 0;
         retained < coordinates.retained_indices.size();
         ++retained) {
        modified.push_back(
            coordinates.modified_factors[retained]
            * physical_fractions[coordinates.retained_indices[retained]]
        );
    }
    const double modified_total = std::accumulate(
        modified.begin(), modified.end(), 0.0
    );
    if (!audit_held2_tolerance(
            kHeld2CompositionSum, modified_total - 1.0
        ).passed) {
        throw std::invalid_argument(
            "modified mole fractions do not sum to one"
        );
    }
    return modified;
}

std::vector<double> held2_lift_modified_fractions(
    const Held2Coordinates& coordinates,
    const std::vector<double>& modified_fractions
) {
    if (modified_fractions.size() != coordinates.retained_indices.size()) {
        throw std::invalid_argument(
            "modified composition size does not match retained species"
        );
    }
    require_finite_vector(modified_fractions, "modified mole fractions");
    double modified_total = 0.0;
    for (double value : modified_fractions) {
        if (value < 0.0) {
            throw std::invalid_argument(
                "modified mole fractions must be nonnegative"
            );
        }
        modified_total += value;
    }
    if (!audit_held2_tolerance(
            kHeld2CompositionSum, modified_total - 1.0
        ).passed) {
        throw std::invalid_argument(
            "modified mole fractions must sum to one"
        );
    }

    std::vector<double> physical(coordinates.charges.size(), 0.0);
    for (std::size_t retained = 0;
         retained < coordinates.retained_indices.size();
         ++retained) {
        physical[coordinates.retained_indices[retained]] =
            modified_fractions[retained]
            / coordinates.modified_factors[retained];
    }
    double eliminated = 0.0;
    const double eliminated_charge =
        coordinates.charges[coordinates.eliminated_index];
    for (std::size_t index : coordinates.retained_indices) {
        eliminated -= coordinates.charges[index] / eliminated_charge
            * physical[index];
    }
    if (eliminated < 0.0
        && !audit_held2_tolerance(
            kHeld2ReconstructedIon, eliminated
        ).passed) {
        throw std::invalid_argument(
            "eliminated ion amount must be nonnegative"
        );
    }
    physical[coordinates.eliminated_index] = std::max(0.0, eliminated);
    return physical;
}

std::vector<double> held2_lift_independent_fractions(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent_modified_fractions
) {
    const std::size_t independent_count =
        coordinates.independent_indices.size();
    if (independent_modified_fractions.size() != independent_count) {
        throw std::invalid_argument(
            "independent modified composition size does not match the HELD2 chart"
        );
    }
    require_finite_vector(
        independent_modified_fractions,
        "independent modified fractions"
    );
    std::vector<double> modified_fractions(
        coordinates.retained_indices.size(), 0.0
    );
    double independent_sum = 0.0;
    for (std::size_t independent = 0;
         independent < independent_count;
         ++independent) {
        const double value = independent_modified_fractions[independent];
        if (value < coordinates.independent_lower_bounds[independent]
            || value > coordinates.independent_upper_bounds[independent]) {
            throw std::invalid_argument(
                "independent modified fraction is outside its source bound"
            );
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
        throw std::invalid_argument(
            "independent modified fractions must leave a positive dependent fraction"
        );
    }
    const auto dependent_retained = static_cast<std::size_t>(
        std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            coordinates.dependent_index
        ) - coordinates.retained_indices.begin()
    );
    modified_fractions[dependent_retained] = 1.0 - independent_sum;
    return held2_lift_modified_fractions(
        coordinates, modified_fractions
    );
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
        throw std::invalid_argument(
            "HELD2 simplex chart dimensions are invalid"
        );
    }
    require_finite_vector(
        unit_cube_coordinates, "HELD2 simplex cube coordinates"
    );
    const auto dependent_retained = static_cast<std::size_t>(
        std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            coordinates.dependent_index
        ) - coordinates.retained_indices.begin()
    );
    const double composition_sum_upper = 1.0
        - kHeld2ModifiedLowerScale
            * coordinates.modified_factors[dependent_retained];
    const double lower_sum = std::accumulate(
        coordinates.independent_lower_bounds.begin(),
        coordinates.independent_lower_bounds.end(),
        0.0
    );
    double remaining_composition = composition_sum_upper - lower_sum;
    if (!(remaining_composition > 0.0)) {
        throw std::invalid_argument(
            "HELD2 simplex chart has no feasible interior"
        );
    }
    double charged_lower_sum = 0.0;
    for (std::size_t index = 0; index < dimension; ++index) {
        if (coordinates.charges[
                coordinates.independent_indices[index]
            ] != 0.0) {
            charged_lower_sum +=
                coordinates.independent_lower_bounds[index];
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
            throw std::invalid_argument(
                "HELD2 simplex chart bounds are invalid"
            );
        }
        if (unit_cube_coordinates[index] < 0.0
            || unit_cube_coordinates[index] > 1.0) {
            throw std::invalid_argument(
                "HELD2 simplex cube coordinate is outside [0, 1]"
            );
        }
    }

    std::vector<double> independent(dimension, 0.0);
    for (std::size_t index = 0; index < dimension; ++index) {
        const bool charged =
            coordinates.charges[
                coordinates.independent_indices[index]
            ] != 0.0;
        double available = charged
            ? std::min(remaining_composition, remaining_ion)
            : remaining_composition;
        available = std::min(
            available,
            coordinates.independent_upper_bounds[index]
                - coordinates.independent_lower_bounds[index]
        );
        const double shifted = available * unit_cube_coordinates[index];
        independent[index] =
            coordinates.independent_lower_bounds[index] + shifted;
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
        throw std::invalid_argument(
            "chemical-potential size does not match charges"
        );
    }
    require_finite_vector(chemical_potentials, "chemical potentials");
    const double eliminated_potential =
        chemical_potentials[coordinates.eliminated_index];
    const double eliminated_charge =
        coordinates.charges[coordinates.eliminated_index];
    std::vector<double> modified;
    modified.reserve(coordinates.retained_indices.size());
    for (std::size_t retained = 0;
         retained < coordinates.retained_indices.size();
         ++retained) {
        const std::size_t index =
            coordinates.retained_indices[retained];
        const double constrained = chemical_potentials[index]
            - coordinates.charges[index] / eliminated_charge
                * eliminated_potential;
        modified.push_back(
            constrained / coordinates.modified_factors[retained]
        );
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
    const std::size_t independent_count =
        coordinates.independent_indices.size();
    const std::size_t provider_coordinate_count = component_count + 1;
    const std::size_t reduced_coordinate_count = independent_count + 1;
    if (independent_modified_fractions.size() != independent_count) {
        throw std::invalid_argument(
            "independent modified composition size does not match the HELD2 chart"
        );
    }
    require_finite_vector(
        independent_modified_fractions,
        "independent modified fractions"
    );
    if (!std::isfinite(log_volume) || !std::isfinite(pressure_over_rt)
        || !std::isfinite(target_pressure_pa) || pressure_over_rt <= 0.0
        || target_pressure_pa <= 0.0) {
        throw std::invalid_argument(
            "HELD2 phase state and pressure scales must be finite and positive"
        );
    }
    if (block.gradient.size() != provider_coordinate_count
        || block.hessian.size()
            != provider_coordinate_count * provider_coordinate_count) {
        throw std::invalid_argument(
            "HELD2 physical phase block has the wrong tensor dimensions"
        );
    }
    require_finite_vector(
        block.gradient, "HELD2 physical phase gradient"
    );
    require_finite_vector(
        block.hessian, "HELD2 physical phase Hessian"
    );
    if (!std::isfinite(block.helmholtz_over_rt)
        || !std::isfinite(block.pressure_pa)) {
        throw std::invalid_argument(
            "HELD2 physical phase block scalars must be finite"
        );
    }

    Held2StateEvaluation result;
    result.modified_fractions.resize(
        coordinates.retained_indices.size(), 0.0
    );
    for (std::size_t independent = 0;
         independent < independent_count;
         ++independent) {
        const double value =
            independent_modified_fractions[independent];
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
        coordinates, independent_modified_fractions
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
        throw std::invalid_argument(
            "HELD2 phase volume must be finite and positive"
        );
    }

    std::vector<double> jacobian(
        provider_coordinate_count * reduced_coordinate_count, 0.0
    );
    for (std::size_t independent = 0;
         independent < independent_count;
         ++independent) {
        for (std::size_t retained = 0;
             retained < coordinates.retained_indices.size();
             ++retained) {
            const std::size_t component =
                coordinates.retained_indices[retained];
            double modified_derivative = 0.0;
            if (component
                == coordinates.independent_indices[independent]) {
                modified_derivative = 1.0;
            } else if (component == coordinates.dependent_index) {
                modified_derivative = -1.0;
            }
            jacobian[
                component * reduced_coordinate_count + independent
            ] = modified_derivative
                / coordinates.modified_factors[retained];
        }
        double eliminated_derivative = 0.0;
        const double eliminated_charge =
            coordinates.charges[coordinates.eliminated_index];
        for (std::size_t component : coordinates.retained_indices) {
            eliminated_derivative -=
                coordinates.charges[component] / eliminated_charge
                * jacobian[
                    component * reduced_coordinate_count + independent
                ];
        }
        jacobian[
            coordinates.eliminated_index * reduced_coordinate_count
                + independent
        ] = eliminated_derivative;
    }
    jacobian[
        component_count * reduced_coordinate_count + independent_count
    ] = result.volume;

    std::vector<double> augmented_gradient = block.gradient;
    augmented_gradient.back() += pressure_over_rt;
    result.objective =
        block.helmholtz_over_rt + pressure_over_rt * result.volume;
    result.gradient.assign(reduced_coordinate_count, 0.0);
    for (std::size_t reduced = 0;
         reduced < reduced_coordinate_count;
         ++reduced) {
        for (std::size_t provider = 0;
             provider < provider_coordinate_count;
             ++provider) {
            result.gradient[reduced] +=
                jacobian[
                    provider * reduced_coordinate_count + reduced
                ] * augmented_gradient[provider];
        }
    }
    result.hessian.assign(
        reduced_coordinate_count * reduced_coordinate_count, 0.0
    );
    for (std::size_t row = 0;
         row < reduced_coordinate_count;
         ++row) {
        for (std::size_t column = 0;
             column < reduced_coordinate_count;
             ++column) {
            double value = 0.0;
            for (std::size_t left = 0;
                 left < provider_coordinate_count;
                 ++left) {
                for (std::size_t right = 0;
                     right < provider_coordinate_count;
                     ++right) {
                    value += jacobian[
                        left * reduced_coordinate_count + row
                    ] * block.hessian[
                        left * provider_coordinate_count + right
                    ] * jacobian[
                        right * reduced_coordinate_count + column
                    ];
                }
            }
            result.hessian[
                row * reduced_coordinate_count + column
            ] = value;
        }
    }
    result.hessian.back() +=
        result.volume * augmented_gradient.back();
    const std::vector<double> physical_potentials(
        block.gradient.begin(), block.gradient.end() - 1
    );
    result.modified_potentials =
        held2_transform_modified_potentials(
            coordinates, physical_potentials
        );
    result.pressure_stationarity_relative =
        (block.pressure_pa - target_pressure_pa) / target_pressure_pa;
    result.pressure_stationarity_derivative_log_volume =
        (result.gradient.back() - result.hessian.back())
        / (pressure_over_rt * result.volume);
    return result;
}

}  // namespace epcsaft_equilibrium
