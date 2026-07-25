#include "held2_step_checks.hpp"

#include "held2_step1.hpp"

#include <algorithm>
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
    double actual,
    double expected,
    double tolerance,
    const char* message
) {
    require(
        std::isfinite(actual)
            && std::abs(actual - expected) <= tolerance,
        message
    );
}

void require_vector_close(
    const std::vector<double>& actual,
    const std::vector<double>& expected,
    double tolerance,
    const char* message
) {
    require(actual.size() == expected.size(), message);
    for (std::size_t index = 0; index < actual.size(); ++index) {
        require_close(actual[index], expected[index], tolerance, message);
    }
}

Held2PhysicalVolumeBoundsEvaluator volume_bounds(
    std::vector<double>* observed_physical = nullptr
) {
    return [observed_physical](const std::vector<double>& physical) {
        if (observed_physical != nullptr) {
            *observed_physical = physical;
        }
        return std::array<double, 2>{
            1.0e-5 + 1.0e-6 * physical.at(2),
            1.0e-3,
        };
    };
}

void check_licl_interleaved_order() {
    const std::vector<std::string> components{
        "water", "lithium-cation", "1-butanol", "chloride-anion",
    };
    const std::vector<double> charges{0.0, 1.0, 0.0, -1.0};
    const std::vector<double> feed{0.70, 0.05, 0.20, 0.05};
    std::vector<double> observed_physical;

    const Held2Step1Result result = run_held2_step1(
        components,
        charges,
        298.15,
        100000.0,
        feed,
        volume_bounds(&observed_physical)
    );

    require(result.status == "complete", "LiCl Step 1 did not complete");
    require(result.reason == "step1_complete", "LiCl Step 1 reason changed");
    require(result.coordinates.has_value(), "LiCl coordinates are absent");
    require(result.independent_feed.has_value(), "LiCl feed is absent");
    require(result.volume_bounds.has_value(), "LiCl volume domain is absent");
    require(result.timing.invocation_count == 1, "Step-1 invocation count changed");
    require(
        result.timing.provider_evaluations == 1,
        "Step 1 did not count its Provider domain evaluation"
    );
    require(result.timing.optimizer_solves == 0, "Step 1 ran an optimizer");
    require(result.timing.optimizer_iterations == 0, "Step 1 ran optimizer iterations");
    require(result.timing.next_step == 2, "Step 1 did not transition to Step 2");
    require(
        result.timing.terminal_status == result.status
            && result.timing.terminal_reason == result.reason,
        "Step-1 timing evidence diverged from its result"
    );
    const Held2Coordinates& coordinates = *result.coordinates;
    require(coordinates.eliminated_index == 3, "LiCl eliminated ion changed");
    require(coordinates.dependent_index == 2, "LiCl closure species changed");
    require(
        coordinates.paper_to_provider_indices
            == std::vector<std::size_t>({1, 3, 0, 2}),
        "LiCl paper-to-Provider permutation is incorrect"
    );
    require(
        coordinates.provider_to_paper_indices
            == std::vector<std::size_t>({2, 0, 3, 1}),
        "LiCl Provider-to-paper permutation is incorrect"
    );
    require(
        coordinates.compact_to_paper_indices
            == std::vector<std::size_t>({0, 2}),
        "LiCl compact-to-paper map is incorrect"
    );
    require(
        coordinates.independent_indices
            == std::vector<std::size_t>({1, 0}),
        "LiCl compact coordinates do not follow paper order"
    );
    require_vector_close(
        coordinates.independent_lower_bounds,
        {2.0e-10, 1.0e-10},
        1.0e-20,
        "LiCl lower bounds do not match Eq. (61)"
    );
    require_vector_close(
        coordinates.independent_upper_bounds,
        {1.0, 1.0},
        1.0e-15,
        "LiCl upper bounds do not match corrected Eqs. (59)-(60)"
    );
    require(
        coordinates.polytope_constraints.size() == 6,
        "LiCl complete Step-1 polytope is incomplete"
    );
    require_vector_close(
        *result.independent_feed,
        {0.10, 0.70},
        1.0e-15,
        "LiCl transformed feed is incorrect"
    );
    require_vector_close(
        held2_lift_independent_fractions(
            coordinates, *result.independent_feed
        ),
        feed,
        1.0e-15,
        "LiCl inverse lift does not round trip"
    );
    const std::array<double, 2> bounds =
        (*result.volume_bounds)(*result.independent_feed);
    require_close(bounds[0], 1.02e-5, 1.0e-18, "lower volume bound changed");
    require_close(bounds[1], 1.0e-3, 1.0e-18, "upper volume bound changed");
    require_vector_close(
        observed_physical,
        feed,
        1.0e-15,
        "Provider volume bounds did not receive the physical composition"
    );

    const std::vector<double> potentials{3.0, -2.0, 1.0, 4.0};
    std::vector<double> shifted = potentials;
    for (std::size_t index = 0; index < shifted.size(); ++index) {
        shifted[index] += 17.0 * charges[index];
    }
    require_vector_close(
        held2_transform_modified_potentials(coordinates, potentials),
        held2_transform_modified_potentials(coordinates, shifted),
        1.0e-12,
        "modified potentials are not Galvani-gauge invariant"
    );
}

void check_cacl2_asymmetric_bounds() {
    const Held2Step1Result result = run_held2_step1(
        {"water", "chloride-anion", "1-butanol", "calcium-cation"},
        {0.0, -1.0, 0.0, 2.0},
        298.15,
        100000.0,
        {0.65, 0.10, 0.20, 0.05},
        volume_bounds()
    );

    require(result.status == "complete", "CaCl2 Step 1 did not complete");
    require(result.coordinates.has_value(), "CaCl2 coordinates are absent");
    const Held2Coordinates& coordinates = *result.coordinates;
    require(coordinates.eliminated_index == 3, "CaCl2 eliminated ion changed");
    require_vector_close(
        coordinates.independent_upper_bounds,
        {1.0, 1.0},
        1.0e-15,
        "CaCl2 corrected upper bounds are incorrect"
    );
    require_vector_close(
        held2_lift_independent_fractions(
            coordinates, *result.independent_feed
        ),
        {0.65, 0.10, 0.20, 0.05},
        1.0e-15,
        "CaCl2 inverse lift does not round trip"
    );
}

void check_polytope_mapper_preserves_provider_ceiling() {
    const Held2Step1Result result = run_held2_step1(
        {"water", "cation", "anion"},
        {0.0, 1.0, -1.0},
        298.15,
        100000.0,
        {0.80, 0.10, 0.10},
        volume_bounds()
    );
    const Held2Coordinates& coordinates = *result.coordinates;
    require_vector_close(
        held2_map_unit_cube_to_independent_fractions(
            coordinates, {0.0}, 0.38
        ),
        {2.0e-10},
        1.0e-12,
        "radial polytope map changed its lower endpoint"
    );
    require_vector_close(
        held2_map_unit_cube_to_independent_fractions(
            coordinates, {0.5}, 0.38
        ),
        {0.19},
        1.0e-12,
        "radial polytope map changed its interior anchor"
    );
    require_vector_close(
        held2_map_unit_cube_to_independent_fractions(
            coordinates, {1.0}, 0.38
        ),
        {0.38},
        1.0e-12,
        "radial polytope map violated the Provider ion ceiling"
    );
    static_cast<void>(held2_map_unit_cube_to_independent_fractions(
        coordinates, {0.5}, 1.0
    ));
}

void check_coupled_polytope_constraint() {
    const Held2Step1Result result = run_held2_step1(
        {"water", "monovalent-cation", "anion", "divalent-cation", "solvent"},
        {0.0, 1.0, -1.0, 2.0, 0.0},
        298.15,
        100000.0,
        {0.55, 0.05, 0.20, 0.075, 0.125},
        volume_bounds()
    );
    require(result.status == "complete", "mixed-charge Step 1 did not complete");

    bool rejected = false;
    try {
        static_cast<void>(held2_lift_independent_fractions(
            *result.coordinates,
            {0.20, 0.10, 0.30}
        ));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(
        rejected,
        "complete polytope admitted a negative eliminated-ion amount"
    );
    for (const std::vector<double>& cube : std::vector<std::vector<double>>{
             {1.0, 0.0, 0.0},
             {0.0, 1.0, 0.0},
             {1.0, 1.0, 1.0},
             {0.25, 0.75, 0.40},
         }) {
        const std::vector<double> independent =
            held2_map_unit_cube_to_independent_fractions(
                *result.coordinates, cube, 0.60
            );
        const std::vector<double> physical =
            held2_lift_independent_fractions(
                *result.coordinates, independent
            );
        double charge = 0.0;
        double total_ion = 0.0;
        for (std::size_t index = 0; index < physical.size(); ++index) {
            charge += (*result.coordinates).charges[index] * physical[index];
            if ((*result.coordinates).charges[index] != 0.0) {
                total_ion += physical[index];
            }
        }
        require_close(
            charge, 0.0, 1.0e-12,
            "complete-polytope mapper violated electroneutrality"
        );
        require(
            total_ion <= 0.60 + 1.0e-12,
            "complete-polytope mapper violated the Provider ion ceiling"
        );
    }
}

void check_fail_closed_inputs() {
    const Held2Step1Result singular = run_held2_step1(
        {"cation-a", "cation-b", "anion-a", "anion-b", "water"},
        {1.0, 1.0, -1.0, -1.0, 0.0},
        298.15,
        100000.0,
        {0.10, 0.10, 0.10, 0.10, 0.60},
        volume_bounds()
    );
    require(
        singular.reason == "unsupported_singular_charge_transformation",
        "singular charge topology was not typed"
    );
    require(!singular.coordinates.has_value(), "singular topology fabricated coordinates");

    const Held2Step1Result invalid_temperature = run_held2_step1(
        {"water", "cation", "anion"},
        {0.0, 1.0, -1.0},
        0.0,
        100000.0,
        {0.80, 0.10, 0.10},
        volume_bounds()
    );
    require(
        invalid_temperature.reason == "invalid_temperature",
        "invalid temperature was not typed"
    );
    require(
        invalid_temperature.timing.provider_evaluations == 0
            && invalid_temperature.timing.next_step == 0,
        "invalid temperature fabricated work or a transition"
    );
    const Held2Step1Result invalid_pressure = run_held2_step1(
        {"water", "cation", "anion"},
        {0.0, 1.0, -1.0},
        298.15,
        -1.0,
        {0.80, 0.10, 0.10},
        volume_bounds()
    );
    require(
        invalid_pressure.reason == "invalid_pressure",
        "invalid pressure was not typed"
    );

    const Held2Step1Result empty_domain = run_held2_step1(
        {"water", "cation", "anion"},
        {0.0, 1.0, -1.0},
        298.15,
        100000.0,
        {0.80, 0.10, 0.10},
        [](const std::vector<double>&) {
            return std::array<double, 2>{1.0e-3, 1.0e-3};
        }
    );
    require(
        empty_domain.reason == "empty_physical_volume_domain",
        "empty Provider volume domain was not typed"
    );
    require(
        !empty_domain.volume_bounds.has_value(),
        "empty Provider domain fabricated an evaluator"
    );
    const Held2Step1Result provider_failure = run_held2_step1(
        {"water", "cation", "anion"},
        {0.0, 1.0, -1.0},
        298.15,
        100000.0,
        {0.80, 0.10, 0.10},
        [](const std::vector<double>&) -> std::array<double, 2> {
            throw std::runtime_error("manufactured Provider failure");
        }
    );
    require(
        provider_failure.reason == "provider_volume_domain_failure",
        "Provider callback exception was not typed"
    );

    const Held2Step1Result nonnormalized = run_held2_step1(
        {"water", "cation", "1-butanol", "anion"},
        {0.0, 1.0, 0.0, -1.0},
        298.15,
        100000.0,
        {0.69, 0.05, 0.20, 0.05},
        volume_bounds()
    );
    require(
        nonnormalized.reason == "invalid_feed",
        "nonnormalized feed was not rejected"
    );
    const Held2Step1Result charged = run_held2_step1(
        {"water", "cation", "1-butanol", "anion"},
        {0.0, 1.0, 0.0, -1.0},
        298.15,
        100000.0,
        {0.69, 0.06, 0.20, 0.05},
        volume_bounds()
    );
    require(
        charged.reason == "invalid_feed",
        "nonelectroneutral feed was not rejected"
    );

    const std::shared_ptr<int> calls = std::make_shared<int>(0);
    const Held2Step1Result deferred_failure = run_held2_step1(
        {"water", "cation", "anion"},
        {0.0, 1.0, -1.0},
        298.15,
        100000.0,
        {0.80, 0.10, 0.10},
        [calls](const std::vector<double>&) {
            ++*calls;
            return *calls == 1
                ? std::array<double, 2>{1.0e-5, 1.0e-3}
                : std::array<double, 2>{1.0e-3, 1.0e-3};
        }
    );
    require(
        deferred_failure.status == "complete",
        "valid feed did not retain its transformed volume evaluator"
    );
    bool later_domain_rejected = false;
    try {
        static_cast<void>((*deferred_failure.volume_bounds)(
            *deferred_failure.independent_feed
        ));
    } catch (const std::invalid_argument&) {
        later_domain_rejected = true;
    }
    require(
        later_domain_rejected,
        "later empty Provider volume domain was not rejected"
    );
}

}  // namespace

void run_held2_step1_checks() {
    check_licl_interleaved_order();
    check_cacl2_asymmetric_bounds();
    check_polytope_mapper_preserves_provider_ceiling();
    check_coupled_polytope_constraint();
    check_fail_closed_inputs();
}

}  // namespace epcsaft_equilibrium::test
