#include "held2_step2.hpp"
#include "held2_step4.hpp"

#include <array>
#include <cmath>
#include <functional>
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

Held2StateEvaluation quadratic_state(
    const std::vector<double>& independent,
    double log_volume
) {
    const double composition = independent.front();
    const double shifted = composition - 0.4;
    Held2StateEvaluation state;
    state.modified_fractions = {1.0 - composition, composition};
    state.physical_amounts = state.modified_fractions;
    state.volume = std::exp(log_volume);
    state.objective = 5.0 * shifted + 2.0 * shifted * shifted
        + 3.0 * shifted * log_volume + 4.0 * log_volume * log_volume;
    state.gradient = {
        5.0 + 4.0 * shifted + 3.0 * log_volume,
        3.0 * shifted + 8.0 * log_volume,
    };
    state.hessian = {4.0, 3.0, 3.0, 8.0};
    state.pressure_stationarity_relative = -state.gradient.back();
    state.pressure_stationarity_derivative_log_volume = -8.0;
    return state;
}

Held2StateEvaluation search_state(
    const std::vector<double>& independent,
    double log_volume,
    bool negative
) {
    const double composition = independent.front();
    const double delta = composition - 0.5;
    double objective = delta * delta;
    double gradient = 2.0 * delta;
    double curvature = 2.0;
    if (negative) {
        const double scaled = (composition - 0.75) / 0.08;
        const double well = std::exp(-scaled * scaled);
        objective -= 0.10 * well;
        gradient += 0.2 * scaled * well / 0.08;
        curvature += 0.2 * well * (1.0 - 2.0 * scaled * scaled)
            / (0.08 * 0.08);
    }
    Held2StateEvaluation state;
    state.modified_fractions = {1.0 - composition, composition};
    state.physical_amounts = state.modified_fractions;
    state.volume = std::exp(log_volume);
    state.objective = objective + 0.5 * log_volume * log_volume;
    state.gradient = {gradient, log_volume};
    state.hessian = {curvature, 0.0, 0.0, 1.0};
    state.pressure_stationarity_relative = -log_volume;
    state.pressure_stationarity_derivative_log_volume = -1.0;
    return state;
}

Held2StateEvaluation tied_state(
    const std::vector<double>& independent,
    double log_volume
) {
    Held2StateEvaluation state = search_state(
        independent, log_volume, false
    );
    state.objective += 0.25 * std::pow(log_volume, 4)
        - 0.5 * log_volume * log_volume;
    state.gradient.back() = log_volume * (log_volume * log_volume - 1.0);
    state.hessian.back() = 3.0 * log_volume * log_volume - 1.0;
    state.pressure_stationarity_relative = -state.gradient.back();
    state.pressure_stationarity_derivative_log_volume =
        -state.hessian.back();
    return state;
}

Held2Step3Result appendix_c_result() {
    const Held2Step1Result prepared = step1(
        {0.0, 1.0, -1.0},
        {0.50, 0.25, 0.25},
        [](const std::vector<double>&) {
            return std::array<double, 2>{
                std::exp(-1.5), std::exp(1.5),
            };
        }
    );
    Held2Step2Result step2;
    step2.outcome = Held2Step2Outcome::NegativeWitness;
    step2.reference = search_state({0.5}, 0.0, false);
    return run_held2_step3(
        prepared,
        step2,
        [](const std::vector<double>& independent) {
            return evaluate_held2_pressure_envelope(
                independent,
                {std::exp(-1.5), std::exp(1.5)},
                [](const auto& composition, double log_volume) {
                    return search_state(composition, log_volume, false);
                },
                64
            );
        }
    );
}

}  // namespace

void run_held2_step1_checks() {
    check_coordinates();
    check_polytope();
    check_failures();
}

std::string run_held2_step2_checks(Held2ProgressObserver* observer) {
    const Held2StateEvaluation reference = quadratic_state({0.4}, 0.0);
    const Held2TpdEvaluation tpd = evaluate_held2_tpd(
        reference, {0.4}, quadratic_state({0.6}, 0.1), {0.6}
    );
    require(
        std::abs(evaluate_held2_tpd(
            reference, {0.4}, reference, {0.4}
        ).value) <= 1.0e-12,
        "TPD is not zero at the reference"
    );
    require_close(
        tpd.gradient, {1.1, 1.4}, 1.0e-12,
        "corrected Eq. (62) gradient changed"
    );
    require_close(
        tpd.hessian, {4.0, 3.0, 3.0, 8.0}, 1.0e-12,
        "corrected Eq. (62) Hessian changed"
    );
    constexpr double step = 1.0e-6;
    const auto evaluate = [&](double composition, double log_volume) {
        return evaluate_held2_tpd(
            reference,
            {0.4},
            quadratic_state({composition}, log_volume),
            {composition}
        );
    };
    require_close(
        {
            (evaluate(0.6 + step, 0.1).value
                - evaluate(0.6 - step, 0.1).value) / (2.0 * step),
            (evaluate(0.6, 0.1 + step).value
                - evaluate(0.6, 0.1 - step).value) / (2.0 * step),
        },
        tpd.gradient,
        1.0e-9,
        "TPD centered gradient check failed"
    );
    const Held2TpdEvaluation plus = evaluate(
        0.6 + 0.6 * step, 0.1 - 0.8 * step
    );
    const Held2TpdEvaluation minus = evaluate(
        0.6 - 0.6 * step, 0.1 + 0.8 * step
    );
    require_close(
        {
            (plus.gradient[0] - minus.gradient[0]) / (2.0 * step),
            (plus.gradient[1] - minus.gradient[1]) / (2.0 * step),
        },
        {0.0, -4.6},
        1.0e-9,
        "TPD Hessian-vector check failed"
    );

    const Held2Step1Result prepared = step1(
        {0.0, 1.0, -1.0},
        {0.50, 0.25, 0.25},
        [](const std::vector<double>&) {
            return std::array<double, 2>{
                std::exp(-1.5), std::exp(1.5),
            };
        }
    );
    const Held2Step2Result negative = run_held2_step2(
        prepared,
        [](const auto& composition, double log_volume) {
            return search_state(composition, log_volume, true);
        },
        200,
        observer
    );
    require(
        negative.outcome == Held2Step2Outcome::NegativeWitness
            && negative.negative_witness.has_value()
            && negative.negative_witness->tpd < -1.0e-8,
        "Step 2 missed a strict negative TPD witness"
    );
    const Held2Step2Result nonnegative = run_held2_step2(
        prepared,
        [](const auto& composition, double log_volume) {
            return search_state(composition, log_volume, false);
        },
        80
    );
    require(
        nonnegative.outcome
                == Held2Step2Outcome::NoNegativeWitnessDetected
            && nonnegative.globality_certificate == "not_guaranteed"
            && nonnegative.reference_envelope->selected_root_index >= 0,
        "finite Step-2 search semantics changed"
    );
    require(
        run_held2_step2(prepared, tied_state, 80).outcome
            == Held2Step2Outcome::Indeterminate,
        "tied reference roots did not fail closed"
    );
    const std::string summary = negative.reason + '|'
        + nonnegative.reason + "|stable_objective_tie";
    return std::to_string(std::hash<std::string>{}(summary));
}

void run_held2_step3_checks() {
    const Held2Step3Result result = appendix_c_result();
    require(
        result.status == "complete" && result.state->M.size() == 3,
        "Appendix C did not create 1 + 2(C - 2) points"
    );
    require_close(
        {
            result.state->M[0].independent_modified_fractions[0],
            result.state->M[1].independent_modified_fractions[0],
            result.state->M[2].independent_modified_fractions[0],
        },
        {0.5, 0.25, 0.75},
        1.0e-12,
        "corrected Appendix-C bracketing changed"
    );
    require(
        result.state->M[0].insertion_id == 0
            && result.state->M[1].insertion_id == 1
            && result.state->M[2].insertion_id == 2,
        "Appendix-C insertion order changed"
    );
}

void run_held2_step4_checks() {
    Held2PersistentState state = std::move(*appendix_c_result().state);
    const Held2Step4Result result = run_held2_step4(state);
    require(
        result.status == "complete" && result.certificate->primal_feasible
            && result.certificate->dual_feasible
            && state.upper_solve_count == 1,
        "Step-4 LP was not certified exactly once"
    );
    require_close(
        {*result.upper_bound},
        {0.0},
        1.0e-10,
        "Step-4 analytic envelope changed"
    );
    require(
        std::abs(result.multipliers->front()) <= 0.25 + 1.0e-10,
        "Step-4 multiplier left the analytic optimum face"
    );
}

}  // namespace epcsaft_equilibrium::test
