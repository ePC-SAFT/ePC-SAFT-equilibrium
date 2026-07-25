#include "held2_step1.hpp"

#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace epcsaft_equilibrium::test {
namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_close(
    const std::vector<double>& actual,
    const std::vector<double>& expected,
    double tolerance,
    const char* message
) {
    require(actual.size() == expected.size(), message);
    for (std::size_t index = 0; index < actual.size(); ++index) {
        require(
            std::isfinite(actual[index])
                && std::abs(actual[index] - expected[index]) <= tolerance,
            message
        );
    }
}

Held2PhysicalVolumeBoundsEvaluator volume_bounds() {
    return [](const std::vector<double>& physical) {
        return std::array<double, 2>{
            1.0e-5 + 1.0e-6 * physical.at(2), 1.0e-3,
        };
    };
}

Held2Step1Result step1(
    const std::vector<double>& charges,
    const std::vector<double>& feed,
    const Held2PhysicalVolumeBoundsEvaluator& bounds = volume_bounds(),
    double temperature = 298.15,
    double pressure = 100000.0
) {
    std::vector<std::string> ids(charges.size());
    for (std::size_t index = 0; index < ids.size(); ++index) {
        ids[index] = "component-" + std::to_string(index);
    }
    return run_held2_step1(
        ids, charges, temperature, pressure, feed, bounds
    );
}

void expect_reason(const Held2Step1Result& result, const char* reason) {
    require(
        result.status == "indeterminate" && result.reason == reason
            && result.timing.next_step == 0,
        reason
    );
}

void check_coordinates() {
    const std::vector<double> charges{0.0, 1.0, 0.0, -1.0};
    const std::vector<double> feed{0.70, 0.05, 0.20, 0.05};
    const Held2Step1Result result = step1(charges, feed);
    require(
        result.status == "complete" && result.reason == "step1_complete"
            && result.timing.invocation_count == 1
            && result.timing.provider_evaluations == 1
            && result.timing.optimizer_solves == 0
            && result.timing.next_step == 2,
        "LiCl Step 1 evidence changed"
    );
    const Held2Coordinates& coordinates = *result.coordinates;
    require(
        coordinates.eliminated_index == 3
            && coordinates.dependent_index == 2
            && coordinates.paper_to_provider_indices
                == std::vector<std::size_t>({1, 3, 0, 2})
            && coordinates.provider_to_paper_indices
                == std::vector<std::size_t>({2, 0, 3, 1})
            && coordinates.compact_to_paper_indices
                == std::vector<std::size_t>({0, 2})
            && coordinates.independent_indices
                == std::vector<std::size_t>({1, 0}),
        "LiCl paper coordinate contract changed"
    );
    require_close(
        coordinates.independent_lower_bounds,
        {2.0e-10, 1.0e-10},
        1.0e-20,
        "LiCl Eq. (61) bounds changed"
    );
    require_close(
        coordinates.independent_upper_bounds,
        {1.0, 1.0},
        1.0e-15,
        "LiCl corrected Eqs. (59)-(60) bounds changed"
    );
    require_close(
        *result.independent_feed, {0.10, 0.70}, 1.0e-15,
        "LiCl transformed feed changed"
    );
    require_close(
        held2_lift_independent_fractions(
            coordinates, *result.independent_feed
        ),
        feed,
        1.0e-15,
        "LiCl inverse lift changed"
    );
    const std::array<double, 2> bounds =
        (*result.volume_bounds)(*result.independent_feed);
    require_close(
        std::vector<double>(bounds.begin(), bounds.end()),
        {1.02e-5, 1.0e-3},
        1.0e-18,
        "Provider volume bounds changed"
    );

    const std::vector<double> potentials{3.0, -2.0, 1.0, 4.0};
    std::vector<double> shifted = potentials;
    for (std::size_t index = 0; index < shifted.size(); ++index) {
        shifted[index] += 17.0 * charges[index];
    }
    require_close(
        held2_transform_modified_potentials(coordinates, potentials),
        held2_transform_modified_potentials(coordinates, shifted),
        1.0e-12,
        "Galvani gauge changed modified potentials"
    );

    const Held2Step1Result calcium = step1(
        {0.0, -1.0, 0.0, 2.0}, {0.65, 0.10, 0.20, 0.05}
    );
    require_close(
        calcium.coordinates->independent_upper_bounds,
        {1.0, 1.0},
        1.0e-15,
        "CaCl2 corrected bounds changed"
    );
    require_close(
        held2_lift_independent_fractions(
            *calcium.coordinates, *calcium.independent_feed
        ),
        {0.65, 0.10, 0.20, 0.05},
        1.0e-15,
        "CaCl2 inverse lift changed"
    );
}

void check_polytope() {
    const Held2Coordinates mixed = *step1(
        {0.0, 1.0, -1.0, 2.0, 0.0},
        {0.55, 0.05, 0.20, 0.075, 0.125}
    ).coordinates;
    bool rejected = false;
    try {
        static_cast<void>(held2_lift_independent_fractions(
            mixed, {0.20, 0.10, 0.30}
        ));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "polytope admitted a negative eliminated ion");
    for (const std::vector<double>& cube :
         std::vector<std::vector<double>>{
             {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
             {1.0, 1.0, 1.0}, {0.25, 0.75, 0.40},
         }) {
        const std::vector<double> physical =
            held2_lift_independent_fractions(
                mixed,
                held2_map_unit_cube_to_independent_fractions(
                    mixed, cube, 0.60
                )
            );
        double charge = 0.0;
        double ions = 0.0;
        for (std::size_t index = 0; index < physical.size(); ++index) {
            charge += mixed.charges[index] * physical[index];
            ions += mixed.charges[index] == 0.0 ? 0.0 : physical[index];
        }
        require(
            std::abs(charge) <= 1.0e-12 && ions <= 0.60 + 1.0e-12,
            "polytope map violated charge or Provider ceiling"
        );
    }
}

void check_failures() {
    expect_reason(
        step1({1.0, 1.0, -1.0, -1.0, 0.0},
              {0.10, 0.10, 0.10, 0.10, 0.60}),
        "unsupported_singular_charge_transformation"
    );
    expect_reason(
        step1({0.0, 1.0, -1.0}, {0.80, 0.10, 0.10}, volume_bounds(), 0.0),
        "invalid_temperature"
    );
    expect_reason(
        step1({0.0, 1.0, -1.0}, {0.80, 0.10, 0.10},
              volume_bounds(), 298.15, -1.0),
        "invalid_pressure"
    );
    for (const std::vector<double>& feed :
         std::vector<std::vector<double>>{
             {0.69, 0.05, 0.20, 0.05}, {0.69, 0.06, 0.20, 0.05},
         }) {
        expect_reason(
            step1({0.0, 1.0, 0.0, -1.0}, feed), "invalid_feed"
        );
    }
    expect_reason(
        step1({0.0, 1.0, -1.0}, {0.80, 0.10, 0.10},
              [](const std::vector<double>&) {
                  return std::array<double, 2>{1.0e-3, 1.0e-3};
              }),
        "empty_physical_volume_domain"
    );
    const std::shared_ptr<int> calls = std::make_shared<int>(0);
    const Held2Step1Result deferred = step1(
        {0.0, 1.0, -1.0},
        {0.80, 0.10, 0.10},
        [calls](const std::vector<double>&) {
            return ++*calls == 1
                ? std::array<double, 2>{1.0e-5, 1.0e-3}
                : std::array<double, 2>{1.0e-3, 1.0e-3};
        }
    );
    bool rejected = false;
    try {
        static_cast<void>((*deferred.volume_bounds)(
            *deferred.independent_feed
        ));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "deferred Provider domain failure was ignored");
}

}  // namespace

void run_held2_step1_checks() {
    check_coordinates();
    check_polytope();
    check_failures();
}

}  // namespace epcsaft_equilibrium::test
