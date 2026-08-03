#include "held2_step1.hpp"
#include "held2_tolerances.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
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

std::size_t retained_position(
    const Held2Coordinates& coordinates,
    std::size_t provider_index
) {
    const auto position = std::find(
        coordinates.retained_indices.begin(),
        coordinates.retained_indices.end(),
        provider_index
    );
    if (position == coordinates.retained_indices.end()) {
        throw std::invalid_argument("species is absent from HELD2 coordinates");
    }
    return static_cast<std::size_t>(
        position - coordinates.retained_indices.begin()
    );
}

void require_complete_polytope(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent,
    Held2CompositionDomain domain = Held2CompositionDomain::FiniteSearch
) {
    for (std::size_t constraint_index = 0;
         constraint_index < coordinates.polytope_constraints.size();
         ++constraint_index) {
        const Held2PolytopeConstraint& constraint =
            coordinates.polytope_constraints[constraint_index];
        if (constraint.coefficients.size() != independent.size()) {
            throw std::invalid_argument(
                "HELD2 Step-1 polytope dimensions are invalid"
            );
        }
        const double value = std::inner_product(
            constraint.coefficients.begin(),
            constraint.coefficients.end(),
            independent.begin(),
            0.0
        );
        const double violation = std::max(
            0.0, value - constraint.upper_bound
        );
        const std::size_t compact = constraint_index / 2;
        const bool charged_coordinate = constraint_index < 2 * independent.size()
            && coordinates.charges[
                coordinates.independent_indices[compact]
            ] != 0.0;
        if (domain == Held2CompositionDomain::TraceRefinement
            && constraint_index % 2 == 0 && charged_coordinate
            && independent[compact] <= 0.0) {
            throw std::invalid_argument(
                "charged trace composition must be positive"
            );
        }
        const bool charged_trace_lower_bound =
            domain == Held2CompositionDomain::TraceRefinement
            && constraint_index % 2 == 0 && charged_coordinate
            && independent[compact] > 0.0;
        if (!std::isfinite(value) || (
            !audit_held2_tolerance(
                kHeld2PolytopeFeasibility, violation
            ).passed && !charged_trace_lower_bound
        )) {
            throw std::invalid_argument(
                "independent modified composition violates "
                + constraint.name
            );
        }
    }
}

std::vector<double> independent_from_modified(
    const Held2Coordinates& coordinates,
    const std::vector<double>& modified
);

std::vector<double> physical_total_ion_coefficients(
    const Held2Coordinates& coordinates
) {
    std::vector<double> coefficients(
        coordinates.independent_indices.size(), 0.0
    );
    const double eliminated_charge =
        coordinates.charges[coordinates.eliminated_index];
    for (std::size_t compact = 0; compact < coefficients.size(); ++compact) {
        const std::size_t provider = coordinates.independent_indices[compact];
        const double factor = coordinates.modified_factors[
            retained_position(coordinates, provider)
        ];
        if (coordinates.charges[provider] != 0.0) {
            coefficients[compact] += 1.0 / factor;
        }
        coefficients[compact] -= coordinates.charges[provider]
            / (eliminated_charge * factor);
    }
    return coefficients;
}

std::vector<double> polytope_anchor(
    const Held2Coordinates& coordinates,
    double total_ion_mole_fraction_max
) {
    std::size_t positive_count = 0;
    std::size_t negative_count = 0;
    std::size_t neutral_count = 0;
    for (double charge : coordinates.charges) {
        positive_count += charge > 0.0 ? 1U : 0U;
        negative_count += charge < 0.0 ? 1U : 0U;
        neutral_count += charge == 0.0 ? 1U : 0U;
    }
    const double ion_total = std::isfinite(total_ion_mole_fraction_max)
        ? 0.5 * total_ion_mole_fraction_max
        : 0.5;
    double charge_weight_sum = 0.0;
    for (double charge : coordinates.charges) {
        if (charge > 0.0) {
            charge_weight_sum += 1.0
                / (static_cast<double>(positive_count) * charge);
        } else if (charge < 0.0) {
            charge_weight_sum += 1.0
                / (static_cast<double>(negative_count) * -charge);
        }
    }
    const double charge_amount = ion_total / charge_weight_sum;
    std::vector<double> physical(coordinates.charges.size(), 0.0);
    for (std::size_t index = 0; index < physical.size(); ++index) {
        const double charge = coordinates.charges[index];
        if (charge > 0.0) {
            physical[index] = charge_amount
                / (static_cast<double>(positive_count) * charge);
        } else if (charge < 0.0) {
            physical[index] = charge_amount
                / (static_cast<double>(negative_count) * -charge);
        } else {
            physical[index] =
                (1.0 - ion_total) / static_cast<double>(neutral_count);
        }
    }
    return independent_from_modified(
        coordinates,
        held2_transform_physical_fractions(coordinates, physical)
    );
}

std::array<double, 2> checked_volume_bounds(
    const Held2PhysicalVolumeBoundsEvaluator& evaluator,
    const std::vector<double>& physical
) {
    const std::array<double, 2> bounds = evaluator(physical);
    if (
        !std::isfinite(bounds[0]) || !std::isfinite(bounds[1])
        || bounds[0] <= 0.0 || bounds[1] <= bounds[0]
    ) {
        throw std::invalid_argument("empty_physical_volume_domain");
    }
    return bounds;
}

std::vector<double> independent_from_modified(
    const Held2Coordinates& coordinates,
    const std::vector<double>& modified
) {
    std::vector<double> independent;
    independent.reserve(coordinates.independent_indices.size());
    for (std::size_t provider_index : coordinates.independent_indices) {
        independent.push_back(
            modified[retained_position(coordinates, provider_index)]
        );
    }
    return independent;
}

Held2Step1Result finish_step1(
    Held2Step1Result result,
    const char* reason,
    int next_step,
    const std::chrono::steady_clock::time_point& wall_start,
    std::clock_t cpu_start
) {
    const char* status = next_step == 0 ? "indeterminate" : "complete";
    result.status = status;
    result.reason = reason;
    result.timing.invocation_count = 1;
    result.timing.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start
    ).count();
    result.timing.cpu_seconds = static_cast<double>(
        std::clock() - cpu_start
    ) / static_cast<double>(CLOCKS_PER_SEC);
    result.timing.terminal_status = status;
    result.timing.terminal_reason = reason;
    result.timing.next_step = next_step;
    return result;
}

}  // namespace

Held2Coordinates make_held2_coordinates(
    const std::vector<double>& charges,
    const std::vector<std::string>& component_ids
) {
    if (charges.size() < 3) {
        throw std::invalid_argument("HELD2 requires at least three species");
    }
    if (!component_ids.empty() && component_ids.size() != charges.size()) {
        throw std::invalid_argument(
            "component IDs and charge numbers must have equal size"
        );
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
    result.dependent_index = charges.size();

    for (std::size_t index = 0; index < charges.size(); ++index) {
        if (charges[index] != 0.0 && index != eliminated) {
            result.paper_to_provider_indices.push_back(index);
        }
    }
    result.paper_to_provider_indices.push_back(eliminated);
    for (std::size_t index = 0; index < charges.size(); ++index) {
        if (charges[index] == 0.0) {
            result.dependent_index = index;
        }
    }
    for (std::size_t index = 0; index < charges.size(); ++index) {
        if (charges[index] == 0.0 && index != result.dependent_index) {
            result.paper_to_provider_indices.push_back(index);
        }
    }
    result.paper_to_provider_indices.push_back(result.dependent_index);
    result.provider_to_paper_indices.assign(charges.size(), charges.size());
    for (std::size_t paper = 0;
         paper < result.paper_to_provider_indices.size();
         ++paper) {
        result.provider_to_paper_indices[
            result.paper_to_provider_indices[paper]
        ] = paper;
    }

    for (std::size_t paper = 0;
         paper < result.paper_to_provider_indices.size();
         ++paper) {
        const std::size_t index = result.paper_to_provider_indices[paper];
        if (index == eliminated) {
            continue;
        }
        result.retained_indices.push_back(index);
        result.modified_factors.push_back(
            1.0 - charges[index] / charges[eliminated]
        );
        if (index == result.dependent_index) {
            continue;
        }
        result.compact_to_paper_indices.push_back(paper);
        result.independent_indices.push_back(index);
        const double factor = result.modified_factors.back();
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

    const std::size_t dimension = result.independent_indices.size();
    for (std::size_t compact = 0; compact < dimension; ++compact) {
        std::vector<double> lower(dimension, 0.0);
        lower[compact] = -1.0;
        const std::string coordinate_name = component_ids.empty()
            ? std::to_string(result.independent_indices[compact])
            : component_ids[result.independent_indices[compact]];
        result.polytope_constraints.push_back({
            "lower_bound:" + coordinate_name,
            std::move(lower),
            -result.independent_lower_bounds[compact],
        });
        std::vector<double> upper(dimension, 0.0);
        upper[compact] = 1.0;
        result.polytope_constraints.push_back({
            "upper_bound:" + coordinate_name,
            std::move(upper),
            result.independent_upper_bounds[compact],
        });
    }
    const double dependent_lower = kHeld2ModifiedLowerScale
        * result.modified_factors[
            retained_position(result, result.dependent_index)
        ];
    result.polytope_constraints.push_back({
        "closure_species_lower_bound",
        std::vector<double>(dimension, 1.0),
        1.0 - dependent_lower,
    });
    std::vector<double> eliminated_nonnegative(dimension, 0.0);
    for (std::size_t compact = 0; compact < dimension; ++compact) {
        const std::size_t provider = result.independent_indices[compact];
        const double factor = result.modified_factors[
            retained_position(result, provider)
        ];
        eliminated_nonnegative[compact] =
            charges[provider] / (charges[eliminated] * factor);
    }
    result.polytope_constraints.push_back({
        "eliminated_ion_nonnegative",
        std::move(eliminated_nonnegative),
        0.0,
    });
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

static std::vector<double> held2_lift_modified_fractions(
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
    const std::vector<double>& independent_modified_fractions,
    Held2CompositionDomain domain
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
    require_complete_polytope(
        coordinates, independent_modified_fractions, domain
    );
    std::vector<double> modified_fractions(
        coordinates.retained_indices.size(), 0.0
    );
    double independent_sum = 0.0;
    for (std::size_t independent = 0;
         independent < independent_count;
         ++independent) {
        const double value = independent_modified_fractions[independent];
        const std::size_t retained = retained_position(
            coordinates, coordinates.independent_indices[independent]
        );
        modified_fractions[retained] = value;
        independent_sum += value;
    }
    const std::size_t dependent_retained = retained_position(
        coordinates, coordinates.dependent_index
    );
    modified_fractions[dependent_retained] = 1.0 - independent_sum;
    return held2_lift_modified_fractions(
        coordinates, modified_fractions
    );
}

Held2Step1Result run_held2_step1(
    const std::vector<std::string>& component_ids,
    const std::vector<double>& charges,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& physical_feed,
    const Held2PhysicalVolumeBoundsEvaluator& physical_volume_bounds,
    std::string provider_fingerprint,
    double total_ion_mole_fraction_max
) {
    const auto wall_start = std::chrono::steady_clock::now();
    const std::clock_t cpu_start = std::clock();
    Held2Step1Result result;
    result.temperature_k = temperature_k;
    result.pressure_pa = pressure_pa;
    result.total_ion_mole_fraction_max =
        total_ion_mole_fraction_max;
    result.provider_fingerprint = std::move(provider_fingerprint);
    const auto finish = [&](const char* reason, int next_step = 0) {
        return finish_step1(
            std::move(result), reason, next_step, wall_start, cpu_start
        );
    };
    if (!std::isfinite(temperature_k) || temperature_k <= 0.0) {
        return finish("invalid_temperature");
    }
    if (!std::isfinite(pressure_pa) || pressure_pa <= 0.0) {
        return finish("invalid_pressure");
    }
    if (result.provider_fingerprint.empty()) {
        return finish("invalid_provider_fingerprint");
    }
    if (
        component_ids.size() != charges.size()
        || physical_feed.size() != charges.size()
        || component_ids.empty()
    ) {
        return finish("invalid_component_metadata");
    }
    for (std::size_t index = 0; index < component_ids.size(); ++index) {
        if (
            component_ids[index].empty()
            || std::find(
                component_ids.begin(),
                component_ids.begin() + static_cast<std::ptrdiff_t>(index),
                component_ids[index]
            ) != component_ids.begin() + static_cast<std::ptrdiff_t>(index)
        ) {
            return finish("invalid_component_metadata");
        }
    }

    Held2Coordinates coordinates;
    try {
        coordinates = make_held2_coordinates(charges, component_ids);
        if (std::isfinite(total_ion_mole_fraction_max)) {
            if (total_ion_mole_fraction_max <= 0.0
                || total_ion_mole_fraction_max > 1.0) {
                return finish("invalid_provider_ion_domain");
            }
            coordinates.polytope_constraints.push_back({
                "provider_total_ion_mole_fraction_max",
                physical_total_ion_coefficients(coordinates),
                total_ion_mole_fraction_max,
            });
            require_complete_polytope(
                coordinates,
                polytope_anchor(
                    coordinates, total_ion_mole_fraction_max
                )
            );
        }
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        const char* reason = message.find("singular") != std::string::npos
            ? "unsupported_singular_charge_transformation"
            : "unsupported_charge_topology";
        return finish(reason);
    }

    std::vector<double> independent_feed;
    std::vector<double> lifted_feed;
    try {
        const std::vector<double> modified =
            held2_transform_physical_fractions(coordinates, physical_feed);
        independent_feed = independent_from_modified(coordinates, modified);
        lifted_feed = held2_lift_independent_fractions(
            coordinates, independent_feed
        );
    } catch (const std::invalid_argument&) {
        return finish("invalid_feed");
    }
    if (!physical_volume_bounds) {
        return finish("missing_physical_volume_bounds");
    }
    std::array<double, 2> feed_volume_bounds;
    try {
        ++result.timing.provider_evaluations;
        feed_volume_bounds = checked_volume_bounds(
            physical_volume_bounds, lifted_feed
        );
    } catch (const std::invalid_argument& error) {
        const std::string reason = error.what();
        return finish(
            reason == "empty_physical_volume_domain"
            ? "empty_physical_volume_domain"
            : "provider_volume_domain_failure"
        );
    } catch (...) {
        return finish("provider_volume_domain_failure");
    }

    result.coordinates = coordinates;
    result.independent_feed = independent_feed;
    result.feed_volume_bounds = feed_volume_bounds;
    result.volume_bounds = Held2VolumeBoundsEvaluator(
        [coordinates, physical_volume_bounds](
            const std::vector<double>& independent
        ) {
            return checked_volume_bounds(
                physical_volume_bounds,
                held2_lift_independent_fractions(coordinates, independent)
            );
        }
    );
    return finish("step1_complete", 2);
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
    if (!std::isnan(total_ion_mole_fraction_max)
        && (!std::isfinite(total_ion_mole_fraction_max)
            || total_ion_mole_fraction_max <= 0.0
            || total_ion_mole_fraction_max > 1.0)) {
        throw std::invalid_argument(
            "HELD2 Provider total-ion mole-fraction ceiling is invalid or infeasible"
        );
    }
    for (std::size_t index = 0; index < dimension; ++index) {
        if (unit_cube_coordinates[index] < 0.0
            || unit_cube_coordinates[index] > 1.0) {
            throw std::invalid_argument(
                "HELD2 simplex cube coordinate is outside [0, 1]"
            );
        }
    }

    std::vector<Held2PolytopeConstraint> constraints =
        coordinates.polytope_constraints;
    if (std::isfinite(total_ion_mole_fraction_max)) {
        constraints.push_back({
            "provider_total_ion_mole_fraction_max",
            physical_total_ion_coefficients(coordinates),
            total_ion_mole_fraction_max,
        });
    }
    const std::vector<double> anchor = polytope_anchor(
        coordinates, total_ion_mole_fraction_max
    );
    require_complete_polytope(coordinates, anchor);
    std::vector<double> direction(dimension, 0.0);
    double radial_fraction = 0.0;
    for (std::size_t index = 0; index < dimension; ++index) {
        direction[index] = 2.0 * unit_cube_coordinates[index] - 1.0;
        radial_fraction = std::max(
            radial_fraction, std::abs(direction[index])
        );
    }
    if (radial_fraction == 0.0) {
        return anchor;
    }
    for (double& value : direction) {
        value /= radial_fraction;
    }
    double boundary_distance = std::numeric_limits<double>::infinity();
    for (const Held2PolytopeConstraint& constraint : constraints) {
        const double anchor_value = std::inner_product(
            constraint.coefficients.begin(),
            constraint.coefficients.end(),
            anchor.begin(),
            0.0
        );
        const double direction_value = std::inner_product(
            constraint.coefficients.begin(),
            constraint.coefficients.end(),
            direction.begin(),
            0.0
        );
        if (direction_value > 0.0) {
            boundary_distance = std::min(
                boundary_distance,
                (constraint.upper_bound - anchor_value) / direction_value
            );
        }
    }
    if (!std::isfinite(boundary_distance) || boundary_distance < 0.0) {
        throw std::invalid_argument(
            "HELD2 complete polytope has no radial boundary"
        );
    }
    std::vector<double> independent = anchor;
    for (std::size_t index = 0; index < dimension; ++index) {
        independent[index] +=
            radial_fraction * boundary_distance * direction[index];
    }
    require_complete_polytope(coordinates, independent);
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
    const Held2PhysicalPhaseBlock& block,
    Held2CompositionDomain domain
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
        const std::size_t retained = retained_position(
            coordinates, coordinates.independent_indices[independent]
        );
        result.modified_fractions[retained] = value;
    }
    result.physical_amounts = held2_lift_independent_fractions(
        coordinates, independent_modified_fractions, domain
    );
    const std::size_t dependent_retained = retained_position(
        coordinates, coordinates.dependent_index
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
    result.helmholtz_over_rt_reference_amount =
        block.helmholtz_over_rt;
    result.pressure_pa = block.pressure_pa;
    result.chemical_potentials_over_rt = std::move(physical_potentials);
    return result;
}

}  // namespace epcsaft_equilibrium
