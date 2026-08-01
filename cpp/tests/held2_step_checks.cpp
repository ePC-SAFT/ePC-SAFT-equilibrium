#include "held2_algorithm.hpp"
#include "held2_step2.hpp"
#include "held2_step4.hpp"
#include "held2_step5.hpp"
#include "held2_step7.hpp"
#include "held2_step8.hpp"
#include "held2_step9.hpp"
#include "held2_step10.hpp"
#include "held2_tolerances.hpp"
#include "flash.hpp"
#include "result_json.hpp"

#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
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

double manufactured_gibbs(double composition, double volume) {
    const double shifted = composition - 0.5;
    constexpr double inner_squared = 0.15 * 0.15;
    constexpr double outer_squared = 0.30 * 0.30;
    const double composition_energy = 1.0e4
        * (std::pow(shifted, 6) / 6.0
           - (inner_squared + outer_squared) * std::pow(shifted, 4) / 4.0
           + inner_squared * outer_squared * shifted * shifted / 2.0);
    const double volume_delta = volume - 1.2;
    return composition_energy + 2.5 * volume_delta * volume_delta + volume;
}

std::array<double, 2> manufactured_gibbs_gradient(
    double composition,
    double volume
) {
    const double shifted = composition - 0.5;
    constexpr double inner_squared = 0.15 * 0.15;
    constexpr double outer_squared = 0.30 * 0.30;
    return {
        1.0e4
            * (std::pow(shifted, 5)
               - (inner_squared + outer_squared) * std::pow(shifted, 3)
               + inner_squared * outer_squared * shifted),
        5.0 * (volume - 1.2) + 1.0,
    };
}

Held2StateEvaluation evaluate_manufactured_state_impl(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent,
    double log_volume
) {
    if (independent.size() != 1 || !std::isfinite(log_volume)) {
        throw std::invalid_argument("manufactured HELD2 phase state has the wrong size");
    }
    Held2StateEvaluation state;
    state.physical_amounts = held2_lift_independent_fractions(coordinates, independent);
    state.modified_fractions = held2_transform_physical_fractions(
        coordinates, state.physical_amounts
    );
    state.volume = std::exp(log_volume);
    if (!std::isfinite(state.volume) || state.volume <= 0.0) {
        throw std::invalid_argument("manufactured HELD2 phase volume must be positive");
    }
    const double composition = independent.front();
    const auto gradient = manufactured_gibbs_gradient(composition, state.volume);
    state.objective = manufactured_gibbs(composition, state.volume);
    state.gradient = {gradient[0], state.volume * gradient[1]};
    const double shifted = composition - 0.5;
    constexpr double inner_squared = 0.15 * 0.15;
    constexpr double outer_squared = 0.30 * 0.30;
    state.hessian = {
        1.0e4
            * (5.0 * std::pow(shifted, 4)
               - 3.0 * (inner_squared + outer_squared) * shifted * shifted
               + inner_squared * outer_squared),
        0.0,
        0.0,
        5.0 * state.volume * state.volume + state.volume * gradient[1],
    };
    const double common_volume_term = -state.volume * gradient[1];
    state.modified_potentials = {
        state.objective - composition * gradient[0] + common_volume_term,
        state.objective + (1.0 - composition) * gradient[0] + common_volume_term,
    };
    state.chemical_potentials_over_rt = {
        state.modified_potentials[0],
        2.0 * state.modified_potentials[1],
        0.0,
    };
    state.pressure_stationarity_relative = gradient[1];
    state.pressure_stationarity_derivative_log_volume = 5.0 * state.volume;
    return state;
}

Held2Step1Result step1(
    const std::vector<double>& charges,
    const std::vector<double>& feed,
    const Held2PhysicalVolumeBoundsEvaluator& bounds = volume_bounds(),
    double temperature = 298.15,
    double pressure = 100000.0,
    double total_ion_mole_fraction_max =
        std::numeric_limits<double>::quiet_NaN()
) {
    std::vector<std::string> ids(charges.size());
    for (std::size_t index = 0; index < ids.size(); ++index) {
        ids[index] = "component-" + std::to_string(index);
    }
    return run_held2_step1(
        ids,
        charges,
        temperature,
        pressure,
        feed,
        bounds,
        total_ion_mole_fraction_max
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
    const Held2Step1Result prepared = step1(
        {0.0, 1.0, -1.0, 2.0, 0.0},
        {0.55, 0.05, 0.20, 0.075, 0.125}
    );
    const Held2Coordinates& mixed = *prepared.coordinates;
    bool rejected = false;
    try {
        static_cast<void>(held2_lift_independent_fractions(
            mixed, {0.20, 0.10, 0.30}
        ));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "polytope admitted a negative eliminated ion");
    std::vector<double> charged_trace = *prepared.independent_feed;
    charged_trace.front() = 1.0e-40;
    static_cast<void>(held2_lift_independent_fractions(
        mixed,
        charged_trace,
        Held2CompositionDomain::TraceRefinement
    ));
    charged_trace.front() = 0.0;
    rejected = false;
    try {
        static_cast<void>(held2_lift_independent_fractions(
            mixed,
            charged_trace,
            Held2CompositionDomain::TraceRefinement
        ));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "trace domain admitted a zero charged fraction");
    std::vector<double> neutral_trace = *prepared.independent_feed;
    neutral_trace.back() = 1.0e-40;
    rejected = false;
    try {
        static_cast<void>(held2_lift_independent_fractions(
            mixed,
            neutral_trace,
            Held2CompositionDomain::TraceRefinement
        ));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(
        rejected,
        "trace domain bypassed a neutral finite-search lower bound"
    );
    std::vector<double> invalid_trace = *prepared.independent_feed;
    invalid_trace.front() = 2.0;
    rejected = false;
    try {
        static_cast<void>(held2_lift_independent_fractions(
            mixed,
            invalid_trace,
            Held2CompositionDomain::TraceRefinement
        ));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(
        rejected,
        "trace domain bypassed a non-lower-bound polytope constraint"
    );
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
    state.physical_amounts = {
        1.0 - composition, 0.5 * composition, 0.5 * composition,
    };
    state.volume = std::exp(log_volume);
    state.objective = objective + 0.5 * log_volume * log_volume;
    state.gradient = {gradient, log_volume};
    state.hessian = {curvature, 0.0, 0.0, 1.0};
    state.chemical_potentials_over_rt = {0.0, 2.0 * gradient, 0.0};
    state.modified_potentials = {0.0, gradient};
    state.pressure_stationarity_relative = -log_volume;
    state.pressure_stationarity_derivative_log_volume = -1.0;
    return state;
}

Held2StateEvaluation multi_root_state(
    const std::vector<double>& independent,
    double log_volume,
    double tilt
) {
    Held2StateEvaluation state = search_state(
        independent, log_volume, false
    );
    state.objective += 0.25 * std::pow(
        log_volume * log_volume - 1.0, 2
    ) - 0.5 * log_volume * log_volume + tilt * log_volume;
    state.gradient.back() =
        log_volume * (log_volume * log_volume - 1.0) + tilt;
    state.hessian.back() = 3.0 * log_volume * log_volume - 1.0;
    state.pressure_stationarity_relative = -state.gradient.back();
    state.pressure_stationarity_derivative_log_volume =
        -state.hessian.back();
    return state;
}

Held2StateEvaluation tied_state(
    const std::vector<double>& independent,
    double log_volume
) {
    return multi_root_state(independent, log_volume, 0.0);
}

Held2StateEvaluation narrow_liquid_state(
    const std::vector<double>& independent,
    double log_volume
) {
    Held2StateEvaluation state = search_state(
        independent, log_volume + 10.0, true
    );
    state.volume = std::exp(log_volume);
    return state;
}

Held2StateEvaluation multiple_root_state(
    const std::vector<double>& independent,
    double log_volume
) {
    return multi_root_state(independent, log_volume, 0.05);
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
    const Held2StateEvaluation witness =
        search_state({0.75}, 0.0, true);
    step2.negative_witness = Held2StageICandidate{
        held2_transform_physical_fractions(
            *prepared.coordinates,
            held2_lift_independent_fractions(
                *prepared.coordinates, {0.75}
            )
        ),
        witness.volume,
        -0.1,
        witness.objective,
        witness.gradient,
    };
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

Held2Step3Result electrolyte_appendix_c_result() {
    const Held2Step1Result prepared = step1(
        {1.0, -1.0, 0.0, 0.0},
        {
            0.11052410423969840,
            0.11052410423969840,
            0.59895663485604511,
            0.17999515666455806,
        },
        volume_bounds(),
        298.15,
        100000.0,
        0.6
    );
    Held2Step2Result step2;
    step2.outcome = Held2Step2Outcome::NegativeWitness;
    step2.reference = search_state({0.2210482084793968, 0.5989566348560451}, 0.0, false);
    const std::vector<double> witness_independent{0.20, 0.55};
    const Held2StateEvaluation witness =
        search_state(witness_independent, 0.0, true);
    step2.negative_witness = Held2StageICandidate{
        held2_transform_physical_fractions(
            *prepared.coordinates,
            held2_lift_independent_fractions(
                *prepared.coordinates, witness_independent
            )
        ),
        witness.volume,
        -0.1,
        witness.objective,
        witness.gradient,
    };
    return run_held2_step3(
        prepared,
        step2,
        [](const std::vector<double>& independent) {
            Held2PressureEnvelopeResult envelope;
            envelope.outcome = "selected";
            envelope.selected_root_index = 0;
            Held2PressureRoot root;
            root.volume = 1.0;
            root.objective = independent.front();
            envelope.roots.push_back(std::move(root));
            return envelope;
        }
    );
}

std::pair<Held2Step1Result, Held2Step6Result> stage_iii_fixture() {
    Held2Step1Result prepared = step1(
        {0.0, 1.0, -1.0},
        {0.50, 0.25, 0.25},
        [](const std::vector<double>&) {
            return std::array<double, 2>{0.5, 1.5};
        }
    );
    Held2Step6Result candidates;
    candidates.status = "complete";
    candidates.candidates = {
        {7, {0.20}, 1.0, 0.20, 0.0, {}, "manufactured"},
        {9, {0.80}, 1.0, 0.80, 0.0, {}, "manufactured"},
    };
    return {std::move(prepared), std::move(candidates)};
}

Held2Step8Result manufactured_step8(
    const Held2Step1Result& prepared,
    const Held2Step6Result& candidates
) {
    return run_held2_step8(
        prepared,
        candidates,
        [coordinates = *prepared.coordinates](
            const auto& composition, double log_volume
        ) {
            return evaluate_manufactured_state_impl(
                coordinates, composition, log_volume
            );
        },
        [](const auto& composition, double) {
            return composition.front();
        }
    );
}

void run_held2_step1_checks() {
    check_coordinates();
    check_polytope();
    check_failures();
}

void run_held2_step2_checks() {
    const Held2StateEvaluation reference = quadratic_state({0.4}, 0.0);
    const Held2StateEvaluation trial = quadratic_state({0.6}, 0.1);
    const Held2TpdEvaluation tpd = evaluate_held2_tpd(
        reference, {0.4}, trial, {0.6}
    );
    const double grepe_reduced_cost = trial.objective
        - (reference.objective - 0.4 * reference.gradient.front())
        - 0.6 * reference.gradient.front();
    require(
        std::abs(tpd.value - grepe_reduced_cost) <= 1.0e-12,
        "HELD2 TPD no longer equals the nonreactive GREPE reduced cost"
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
        200
    );
    require(
        negative.outcome == Held2Step2Outcome::NegativeWitness
            && negative.negative_witness.has_value()
            && negative.negative_witness->tpd < -1.0e-8
            && negative.timing.optimizer_solves == 1,
        "Step 2 missed a strict negative TPD witness"
    );
    const Held2StateEvaluation transient_reference =
        search_state({0.5}, 0.0, true);
    std::optional<std::pair<double, double>> transient_negative;
    const Held2Step2Result recertification_failure = run_held2_step2(
        prepared,
        [&transient_negative, &transient_reference](
            const auto& composition,
            double log_volume
        ) {
            if (transient_negative
                && composition.front() == transient_negative->first
                && log_volume == transient_negative->second) {
                throw std::runtime_error(
                    "manufactured transient negative state"
                );
            }
            Held2StateEvaluation state =
                search_state(composition, log_volume, true);
            if (evaluate_held2_tpd(
                    transient_reference,
                    {0.5},
                    state,
                    composition
                ).value < -1.0e-8) {
                transient_negative = {
                    composition.front(), log_volume,
                };
            }
            return state;
        },
        200
    );
    require(
        recertification_failure.outcome
            == Held2Step2Outcome::Indeterminate,
        "Step 2 accepted a negative state that failed recertification"
    );
    const Held2Step1Result narrow_prepared = step1(
        {0.0, 1.0, -1.0},
        {0.50, 0.25, 0.25},
        [](const std::vector<double>&) {
            return std::array<double, 2>{
                std::exp(-12.0), std::exp(3.0),
            };
        }
    );
    require(
        run_held2_step2(
            narrow_prepared, narrow_liquid_state, 6500
        ).outcome == Held2Step2Outcome::NegativeWitness,
        "Step 2 missed a narrow joint composition-volume TPD basin"
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
    const Held2Step2Result multiple = run_held2_step2(
        prepared, multiple_root_state, 80
    );
    require(
        multiple.outcome
                == Held2Step2Outcome::NoNegativeWitnessDetected
            && multiple.reference_envelope
            && multiple.reference_envelope->roots.size() == 3
            && multiple.reference_envelope->selected_root_index == 0
            && multiple.reference_envelope->roots[0].mechanical_class
                == "strict_stable"
            && multiple.reference_envelope->roots[1].mechanical_class
                == "unstable"
            && multiple.reference_envelope->roots[2].mechanical_class
                == "strict_stable"
            && multiple.reference_envelope->roots[0].objective
                < multiple.reference_envelope->roots[2].objective,
        "multiple-root homogeneous reference selection changed"
    );
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
        {0.5, 0.2500000001, 0.74999999995},
        1.0e-12,
        "polytope-aware Appendix-C bracketing changed"
    );
    require(
        result.state->M[0].insertion_id == 0
            && result.state->M[1].insertion_id == 1
            && result.state->M[2].insertion_id == 2,
        "Appendix-C insertion order changed"
    );
    const Held2Step3Result electrolyte = electrolyte_appendix_c_result();
    require(
        electrolyte.status == "complete"
            && electrolyte.state->M.size() == 6
            && electrolyte.state->M.back().origin
                == "stage_i_negative_witness",
        "Appendix C left the four-component electroneutral simplex"
    );
    for (std::size_t coordinate = 0; coordinate < 2; ++coordinate) {
        require(
            electrolyte.state->M[1 + 2 * coordinate]
                    .independent_modified_fractions[coordinate]
                < electrolyte.state->feed[coordinate]
                && electrolyte.state->M[2 + 2 * coordinate]
                        .independent_modified_fractions[coordinate]
                    > electrolyte.state->feed[coordinate],
            "Appendix-C cuts do not bracket an electrolyte feed"
        );
    }
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

    Held2PersistentState master;
    master.feed = {0.5};
    master.feed_reduced_gibbs = 10.0;
    master.M = {
        {1, {0.2}, 1.0, 0.0, 1.0, {}, "master"},
        {2, {0.8}, 1.0, 0.0, 2.0, {}, "master"},
    };
    const Held2Step4Result dual = run_held2_step4(master);
    constexpr double slope = 5.0 / 3.0;
    require(dual.status == "complete", "restricted-master dual solve failed");
    const double intercept = *dual.upper_bound - 0.5 * slope;
    require_close(
        {
            *dual.upper_bound, dual.multipliers->front(),
            1.0 - intercept - 0.2 * slope,
            2.0 - intercept - 0.8 * slope,
        },
        {1.5, slope, 0.0, 0.0}, 1.0e-12,
        "HELD2 upper envelope no longer equals the restricted phase master"
    );
    master.M.push_back({3, {0.5}, 1.0, 0.0, 1.4, {}, "pricing"});
    const Held2Step4Result expanded = run_held2_step4(master);
    require(expanded.status == "complete", "expanded restricted-master dual failed");
    require_close(
        {1.4 - intercept - 0.5 * slope, *expanded.upper_bound},
        {-0.1, 1.4}, 1.0e-12,
        "negative reduced-cost column did not improve the restricted master"
    );
}

void run_held2_step5_checks() {
    require(
        kHeld2Stage2KktStationarity.atol == 1.0e-6,
        "Step-5 discovery stationarity tolerance changed"
    );
    Held2PersistentState identity_state;
    identity_state.M.push_back({
        42, {0.25}, 1.0, std::numeric_limits<double>::quiet_NaN(),
        0.0, {}, "appendix_c",
    });
    Held2MPoint numerical_copy{
        0, {0.25 + 1.0e-9}, std::exp(1.0e-9),
        std::numeric_limits<double>::quiet_NaN(), 0.0, {}, "step5",
    };
    require(
        !retain_held2_m_point(identity_state, numerical_copy)
            && numerical_copy.insertion_id == 42
            && identity_state.M.size() == 1,
        "Step-5 M identity retained a numerical representation copy"
    );
    Held2MPoint volume_distinct{
        0, {0.25}, std::exp(1.0e-5),
        std::numeric_limits<double>::quiet_NaN(), 0.0, {}, "step5",
    };
    require(
        retain_held2_m_point(identity_state, volume_distinct)
            && identity_state.M.size() == 2,
        "Step-5 M identity collapsed a distinct molar volume"
    );

    std::size_t volume_evaluations = 0;
    const Held2Step1Result prepared = step1(
        {0.0, 1.0, -1.0},
        {0.50, 0.25, 0.25},
        [&volume_evaluations](const std::vector<double>&) {
            ++volume_evaluations;
            return std::array<double, 2>{
                std::exp(-1.5), std::exp(1.5),
            };
        }
    );
    volume_evaluations = 0;
    Held2PersistentState state = std::move(*appendix_c_result().state);
    const Held2Step4Result step4 = run_held2_step4(state);
    std::size_t state_evaluations = 0;
    const Held2Step5Result result = run_held2_step5(
        prepared,
        step4,
        state,
        [&state_evaluations](
            const auto& composition, double log_volume
        ) {
            ++state_evaluations;
            return search_state(composition, log_volume, false);
        },
        {0, 8, 10}
    );
    std::string step5_failure =
        "Step-5 certified persistent multistart changed";
    if (!result.attempts.empty() && result.attempts.back().kkt) {
        step5_failure += " (last KKT: "
            + result.attempts.back().kkt->reason + ")";
    }
    require(
        result.status == "complete" && result.lower_value <= step4.upper_bound
            && !result.attempts.empty()
            && result.attempts.back().accepted
            && result.attempts.front().start_family == "random_interior"
            && result.attempts.back().kkt
            && result.attempts.back().kkt->accepted
            && result.attempts.back().kkt->same_major_iteration
            && result.attempts.back().kkt->step4_binding_valid
            && state.next_start_ordinal == result.starts_consumed,
        step5_failure.c_str()
    );
    require(
        result.timing.provider_evaluations
            == state_evaluations + volume_evaluations,
        "Step-5 Provider evaluation accounting changed"
    );
    Held2PersistentState recovery_state = state;
    recovery_state.step5_requires_new_member = true;
    const Held2Step4Result recovery_step4 = run_held2_step4(recovery_state);
    const Held2MPoint recovery_target = recovery_state.M.front();
    const double recovery_multiplier = recovery_state.multipliers.front();
    const double recovery_offset = recovery_state.upper_bound
        - 100.0 * (1.0 + std::abs(recovery_multiplier));
    Held2ResourceProfile recovery_resources;
    recovery_resources.step5_start_cap = 8;
    recovery_resources.step5_recovery_start_cap = 144;
    const Held2Step5Result exhausted_recovery = run_held2_step5(
        prepared,
        recovery_step4,
        recovery_state,
        [recovery_target, recovery_multiplier, recovery_offset](
            const auto& composition, double log_volume
        ) {
            const double delta = composition.front()
                - recovery_target.independent_modified_fractions.front();
            const double log_volume_delta = log_volume
                - std::log(recovery_target.volume);
            const double composition_gradient =
                recovery_multiplier + 2.0 * delta;
            Held2StateEvaluation evaluated;
            evaluated.modified_fractions = {
                1.0 - composition.front(), composition.front()
            };
            evaluated.physical_amounts = {
                1.0 - composition.front(),
                0.5 * composition.front(),
                0.5 * composition.front(),
            };
            evaluated.volume = std::exp(log_volume);
            evaluated.objective = recovery_offset
                + recovery_multiplier * delta + delta * delta
                + 0.5 * log_volume_delta * log_volume_delta;
            evaluated.gradient = {
                composition_gradient, log_volume_delta
            };
            evaluated.hessian = {2.0, 0.0, 0.0, 1.0};
            evaluated.chemical_potentials_over_rt = {
                0.0, 2.0 * composition_gradient, 0.0
            };
            evaluated.modified_potentials = {
                0.0, composition_gradient
            };
            evaluated.pressure_stationarity_relative = -log_volume_delta;
            evaluated.pressure_stationarity_derivative_log_volume = -1.0;
            return evaluated;
        },
        recovery_resources
    );
    const std::string recovery_failure =
        "Step-5 repeated an exhausted unchanged-hull recovery pass (status="
        + exhausted_recovery.status + ", reason="
        + exhausted_recovery.reason + ", starts="
        + std::to_string(exhausted_recovery.starts_consumed) + ")";
    require(
        exhausted_recovery.status == "indeterminate"
            && exhausted_recovery.reason == "step5_recovery_exhausted"
            && exhausted_recovery.starts_consumed > 1
            && exhausted_recovery.attempts[0].start_family
                == "random_interior"
            && exhausted_recovery.attempts[1].start_family
                == "boundary_aware"
            && std::any_of(
                exhausted_recovery.attempts.begin(),
                exhausted_recovery.attempts.end(),
                [](const Held2LocalCertificate& attempt) {
                    return attempt.start_family
                        == "space_filling_interior";
                }
            )
            && exhausted_recovery.attempts.back().start_family
                == "random_interior"
            && exhausted_recovery.lower_value
            && *exhausted_recovery.lower_value == *std::min_element(
                exhausted_recovery.attempts.begin(),
                exhausted_recovery.attempts.end(),
                [](const Held2LocalCertificate& left,
                   const Held2LocalCertificate& right) {
                    return left.local_value.value_or(
                        std::numeric_limits<double>::infinity()
                    ) < right.local_value.value_or(
                        std::numeric_limits<double>::infinity()
                    );
                }
            )->local_value
            && recovery_state.step5_requires_new_member,
        recovery_failure.c_str()
    );
    Held2PersistentState missing_certificate_state =
        std::move(*appendix_c_result().state);
    Held2Step4Result missing_certificate_step4 =
        run_held2_step4(missing_certificate_state);
    missing_certificate_step4.certificate.reset();
    const Held2Step5Result missing_certificate = run_held2_step5(
        prepared,
        missing_certificate_step4,
        missing_certificate_state,
        [](const auto& composition, double log_volume) {
            return search_state(composition, log_volume, false);
        },
        {0, 8, 10}
    );
    require(
        missing_certificate.status == "indeterminate"
            && missing_certificate.reason == "invalid_step5_input"
            && missing_certificate.attempts.empty(),
        "Step-5 accepted a missing Step-4 certificate"
    );

    Held2PersistentState corrupt_cut_state =
        std::move(*appendix_c_result().state);
    Held2Step4Result corrupt_cut_step4 =
        run_held2_step4(corrupt_cut_state);
    corrupt_cut_step4.active_cut_ids.push_back(
        corrupt_cut_step4.active_cut_ids.front()
    );
    const Held2Step5Result corrupt_cut = run_held2_step5(
        prepared,
        corrupt_cut_step4,
        corrupt_cut_state,
        [](const auto& composition, double log_volume) {
            return search_state(composition, log_volume, false);
        },
        {0, 8, 10}
    );
    require(
        corrupt_cut.status == "indeterminate"
            && !corrupt_cut.attempts.empty()
            && std::none_of(
                corrupt_cut.attempts.begin(),
                corrupt_cut.attempts.end(),
                [](const Held2LocalCertificate& attempt) {
                    return attempt.accepted;
                }
            ),
        "Step-5 accepted corrupted Step-4 active-cut identity"
    );

    Held2PersistentState mutated_cut_state =
        std::move(*appendix_c_result().state);
    const Held2Step4Result mutated_cut_step4 =
        run_held2_step4(mutated_cut_state);
    mutated_cut_state.M.front().reduced_gibbs += 0.25;
    const Held2Step5Result mutated_cut = run_held2_step5(
        prepared,
        mutated_cut_step4,
        mutated_cut_state,
        [](const auto& composition, double log_volume) {
            return search_state(composition, log_volume, false);
        },
        {0, 8, 10}
    );
    require(
        mutated_cut.status == "indeterminate"
            && !mutated_cut.attempts.empty()
            && std::none_of(
                mutated_cut.attempts.begin(),
                mutated_cut.attempts.end(),
                [](const Held2LocalCertificate& attempt) {
                    return attempt.accepted;
                }
            ),
        "Step-5 accepted mutated Step-4 cut data"
    );

    Held2PersistentState invalid_physical_state =
        std::move(*appendix_c_result().state);
    const Held2Step4Result invalid_physical_step4 =
        run_held2_step4(invalid_physical_state);
    const Held2Step5Result invalid_physical = run_held2_step5(
        prepared,
        invalid_physical_step4,
        invalid_physical_state,
        [](const auto& composition, double log_volume) {
            Held2StateEvaluation evaluated =
                search_state(composition, log_volume, false);
            const double shift =
                0.1 * evaluated.physical_amounts[1];
            evaluated.physical_amounts[0] += 2.0 * shift;
            evaluated.physical_amounts[1] -= shift;
            evaluated.physical_amounts[2] -= shift;
            return evaluated;
        },
        {0, 8, 10}
    );
    require(
        invalid_physical.status == "indeterminate"
            && !invalid_physical.attempts.empty()
            && std::none_of(
                invalid_physical.attempts.begin(),
                invalid_physical.attempts.end(),
                [](const Held2LocalCertificate& attempt) {
                    return attempt.accepted;
                }
            ),
        "Step-5 accepted invalid physical composition evidence"
    );

    const double local_volume_target =
        1.0 - state.M.back().independent_modified_fractions.front();
    const Held2Step1Result local_volume_prepared = step1(
        {0.0, 1.0, -1.0},
        {0.50, 0.25, 0.25},
        [local_volume_target](const std::vector<double>& composition) {
            return std::abs(
                composition.front() - local_volume_target
            ) < 1.0e-6
                ? std::array<double, 2>{2.0, 3.0}
                : std::array<double, 2>{
                    std::exp(-1.5), std::exp(1.5),
                };
        }
    );
    Held2PersistentState local_volume_state =
        std::move(*appendix_c_result().state);
    const Held2Step4Result local_volume_step4 =
        run_held2_step4(local_volume_state);
    const Held2Step5Result local_volume = run_held2_step5(
        local_volume_prepared,
        local_volume_step4,
        local_volume_state,
        [](const auto& composition, double log_volume) {
            return search_state(composition, log_volume, false);
        },
        {0, 8, 10}
    );
    require(
        local_volume.status == "indeterminate"
            && !local_volume.attempts.empty()
            && std::none_of(
                local_volume.attempts.begin(),
                local_volume.attempts.end(),
                [](const Held2LocalCertificate& attempt) {
                    return attempt.accepted;
                }
            ),
        "Step-5 certified against sampled rather than terminal volume bounds"
    );

    Held2PersistentState gauged_state =
        std::move(*appendix_c_result().state);
    const Held2Step4Result gauged_step4 =
        run_held2_step4(gauged_state);
    const Held2Step5Result gauged = run_held2_step5(
        prepared,
        gauged_step4,
        gauged_state,
        [](const auto& composition, double log_volume) {
            Held2StateEvaluation evaluated =
                search_state(composition, log_volume, false);
            const std::vector<double> charges{0.0, 1.0, -1.0};
            for (std::size_t index = 0;
                 index < charges.size();
                 ++index) {
                evaluated.chemical_potentials_over_rt[index] +=
                    19.0 * charges[index];
            }
            return evaluated;
        },
        {0, 8, 10}
    );
    require(
        gauged.status == "complete"
            && gauged.attempts.back().kkt
            && gauged.attempts.back().kkt->accepted
            && gauged.attempts.back().kkt
                    ->physical_stationarity_residual_inf
                <= kHeld2Stage2KktStationarity.atol,
        "Step-5 physical KKT certificate changed under a Galvani gauge"
    );

    Held2PersistentState polished_state =
        std::move(*appendix_c_result().state);
    const Held2Step4Result polished_step4 =
        run_held2_step4(polished_state);
    const Held2Step5Result polished = run_held2_step5(
        prepared,
        polished_step4,
        polished_state,
        [](const auto& composition, double log_volume) {
            Held2StateEvaluation evaluated =
                search_state(composition, log_volume, false);
            evaluated.pressure_stationarity_relative += 5.0e-8;
            return evaluated;
        },
        {0, 8, 10}
    );
    const auto polished_attempt = std::find_if(
        polished.attempts.begin(),
        polished.attempts.end(),
        [](const Held2LocalCertificate& attempt) {
            return attempt.accepted
                && attempt.solver_status.find("_pressure_polished")
                    != std::string::npos;
        }
    );
    require(
        polished.status == "complete"
            && polished_attempt != polished.attempts.end()
            && polished_attempt->kkt
            && polished_attempt->kkt->solver_variables
                != polished_attempt->kkt->audited_variables,
        "Step-5 pressure polish obscured raw solver-variable provenance"
    );

    const Held2Step1Result reordered_prepared = step1(
        {1.0, -1.0, 0.0},
        {0.25, 0.25, 0.50},
        [](const std::vector<double>&) {
            return std::array<double, 2>{
                std::exp(-1.5), std::exp(1.5),
            };
        }
    );
    const Held2StateEvaluator reordered_evaluator =
        [](const auto& composition, double log_volume) {
            Held2StateEvaluation evaluated =
                search_state(composition, log_volume, false);
            evaluated.physical_amounts = {
                evaluated.physical_amounts[1],
                evaluated.physical_amounts[2],
                evaluated.physical_amounts[0],
            };
            evaluated.chemical_potentials_over_rt = {
                evaluated.chemical_potentials_over_rt[1],
                evaluated.chemical_potentials_over_rt[2],
                evaluated.chemical_potentials_over_rt[0],
            };
            return evaluated;
        };
    Held2PersistentState mismatched_chart_state =
        std::move(*appendix_c_result().state);
    const Held2Step4Result mismatched_chart_step4 =
        run_held2_step4(mismatched_chart_state);
    const Held2Step5Result mismatched_chart = run_held2_step5(
        reordered_prepared,
        mismatched_chart_step4,
        mismatched_chart_state,
        reordered_evaluator,
        {0, 8, 10}
    );
    require(
        mismatched_chart.status == "indeterminate"
            && mismatched_chart.reason == "invalid_step5_input"
            && mismatched_chart.attempts.empty(),
        "Step-5 accepted a mismatched Step-1 coordinate chart"
    );
    Held2PersistentState reordered_state =
        std::move(*appendix_c_result().state);
    reordered_state.coordinates = *reordered_prepared.coordinates;
    const Held2Step4Result reordered_step4 =
        run_held2_step4(reordered_state);
    const Held2Step5Result reordered = run_held2_step5(
        reordered_prepared,
        reordered_step4,
        reordered_state,
        reordered_evaluator,
        {0, 8, 10}
    );
    require(
        reordered.status == "complete"
            && reordered.attempts.back().kkt
            && reordered.attempts.back().kkt->accepted
            && reordered.attempts.back().kkt
                    ->physical_stationarity_residual_inf
                <= kHeld2Stage2KktStationarity.atol,
        "Step-5 physical KKT certificate changed under species reordering"
    );

    Held2PersistentState bound_state;
    bound_state.coordinates = *prepared.coordinates;
    bound_state.feed = {0.5};
    bound_state.feed_reduced_gibbs = 10.0;
    bound_state.M = {
        {1, {0.2}, 1.0, 0.0, 1.0, {}, "master"},
        {2, {0.8}, 1.0, 0.0, 2.0, {}, "master"},
    };
    const Held2Step4Result bound_step4 =
        run_held2_step4(bound_state);
    const Held2Step5Result active_bound = run_held2_step5(
        prepared,
        bound_step4,
        bound_state,
        [](const auto& composition, double log_volume) {
            return search_state(composition, log_volume, false);
        },
        {0, 8, 10}
    );
    std::string active_bound_failure =
        "Step-5 active composition bound was double-counted";
    if (!active_bound.attempts.empty()
        && active_bound.attempts.back().kkt) {
        active_bound_failure += " (last KKT: "
            + active_bound.attempts.back().kkt->reason
            + ", active="
            + std::to_string(
                active_bound.attempts.back().kkt
                    ->active_constraint_count
            )
            + ", rank="
            + std::to_string(
                active_bound.attempts.back().kkt
                    ->active_jacobian_rank
            )
            + ")";
    }
    require(
        active_bound.status == "complete"
            && active_bound.attempts.back().kkt
            && active_bound.attempts.back().kkt->accepted
            && active_bound.attempts.back().kkt
                    ->active_constraint_count
                == active_bound.attempts.back().kkt
                    ->active_jacobian_rank
            && active_bound.attempts.back().kkt
                    ->active_constraint_count
                >= 1,
        active_bound_failure.c_str()
    );

    Held2PersistentState dilute_state;
    dilute_state.coordinates = *prepared.coordinates;
    dilute_state.feed = {0.5};
    dilute_state.feed_reduced_gibbs = 10.0;
    dilute_state.M = {
        {1, {0.2}, 1.0, 0.0, 0.0, {}, "master"},
        {2, {0.8}, 1.0, 0.0, 0.0, {}, "master"},
    };
    const Held2Step4Result dilute_step4 = run_held2_step4(dilute_state);
    const double dilute_multiplier = dilute_state.multipliers.front();
    const double dilute_lower =
        prepared.coordinates->independent_lower_bounds.front();
    const double dilute_offset = dilute_state.upper_bound - 1.0;
    const Held2Step5Result dilute = run_held2_step5(
        prepared,
        dilute_step4,
        dilute_state,
        [dilute_multiplier, dilute_lower, dilute_offset](
            const auto& composition, double log_volume
        ) {
            constexpr double slope = 2.25e-5;
            const double delta = composition.front() - dilute_lower;
            const double composition_gradient =
                dilute_multiplier + slope + delta;
            Held2StateEvaluation evaluated;
            evaluated.modified_fractions = {
                1.0 - composition.front(), composition.front()
            };
            evaluated.physical_amounts = {
                1.0 - composition.front(),
                0.5 * composition.front(),
                0.5 * composition.front(),
            };
            evaluated.volume = std::exp(log_volume);
            evaluated.objective = dilute_offset
                + dilute_multiplier * delta + slope * delta
                + 0.5 * delta * delta
                + 0.5 * log_volume * log_volume;
            evaluated.gradient = {
                composition_gradient, log_volume
            };
            evaluated.hessian = {1.0, 0.0, 0.0, 1.0};
            evaluated.chemical_potentials_over_rt = {
                0.0, 2.0 * composition_gradient, 0.0
            };
            evaluated.modified_potentials = {
                0.0, composition_gradient
            };
            evaluated.pressure_stationarity_relative = -log_volume;
            evaluated.pressure_stationarity_derivative_log_volume = -1.0;
            return evaluated;
        },
        {0, 8, 10}
    );
    const std::string dilute_failure =
        "Step-5 did not resolve a dilute positive finite-search face"
        " (status=" + dilute.status + ", reason=" + dilute.reason
        + ", starts=" + std::to_string(dilute.starts_consumed)
        + (dilute.attempts.empty() || !dilute.attempts.front().kkt
            ? ")"
            : ", kkt=" + dilute.attempts.front().kkt->reason
                + ", x=" + std::to_string(
                    dilute.attempts.front().kkt
                        ->audited_variables.front()
                )
                + ", x_nano=" + std::to_string(
                    1.0e9 * dilute.attempts.front().kkt
                        ->audited_variables.front()
                )
                + ", residual=" + std::to_string(
                    dilute.attempts.front().kkt
                        ->stationarity_residual_inf
                )
                + ", solver_zl=" + std::to_string(
                    dilute.attempts.front().kkt
                        ->solver_lower_bound_multipliers.front()
                )
                + ", audited_zl=" + std::to_string(
                    dilute.attempts.front().kkt
                        ->lower_bound_multipliers.front()
                ) + ")");
    require(
        dilute.status == "complete"
            && dilute.starts_consumed == 1
            && dilute.timing.optimizer_solves == 2
            && dilute.attempts.size() == 1
            && dilute.attempts.front().accepted
            && dilute.attempts.front().kkt
            && dilute.attempts.front().kkt->accepted
            && dilute.attempts.front().dilute_face_restart
            && dilute.attempts.front().dilute_face_restart->attempted
            && dilute.attempts.front().dilute_face_restart->accepted
            && dilute.attempts.front().dilute_face_restart
                ->coordinate_indices == std::vector<std::size_t>{0}
            && dilute.attempts.front().dilute_face_restart
                ->provider_component_indices
                == std::vector<std::size_t>{
                    prepared.coordinates->independent_indices.front()
                }
            && dilute.attempts.front().dilute_face_restart
                ->lower_bound_distances.size() == 1
            && dilute.attempts.front().dilute_face_restart
                ->complementarity_products.front()
                <= kHeld2Stage2KktComplementarity.atol
            && dilute.attempts.front().kkt->audited_variables.front()
                > 0.0
            && dilute.attempts.front().kkt->audited_variables.front()
                <= dilute_lower + kHeld2BoundActivity.atol,
        dilute_failure.c_str()
    );
}

void run_held2_step6_checks() {
    const Held2Step1Result prepared = step1(
        {0.0, 1.0, -1.0}, {0.50, 0.25, 0.25},
        [](const std::vector<double>&) {
            return std::array<double, 2>{
                std::exp(-1.5), std::exp(1.5),
            };
        }
    );
    Held2PersistentState state = std::move(*appendix_c_result().state);
    const Held2Step4Result step4 = run_held2_step4(state);
    std::size_t packing_evaluations = 0;
    const Held2Step6Result result = run_held2_step6(
        prepared,
        step4,
        state,
        [&packing_evaluations](const auto& composition, double) {
            ++packing_evaluations;
            return composition.front();
        }
    );
    require(
        result.status == "complete"
            && result.decisions.size() == state.M.size()
            && result.timing.next_step == 7,
        "Step-6 full-M search changed"
    );
    require(
        result.timing.provider_evaluations == packing_evaluations,
        "Step-6 Provider evaluation accounting changed"
    );
    require(
        kHeld2PaperStep6Gap.atol == 5.0e-2
            && audit_held2_tolerance(
                kHeld2PaperStep6Gap, 4.0e-2
            ).passed
            && !audit_held2_tolerance(
                kHeld2PaperStep6Gap, 6.0e-2
            ).passed,
        "Step-6 finite-search gap policy changed"
    );
}

void run_held2_step7_checks() {
    Held2PersistentState state = std::move(*appendix_c_result().state);
    Held2Step5Result step5;
    step5.status = "complete";
    Held2Step6Result step6;
    step6.status = "complete";
    const Held2Step7Result next = run_held2_step7(
        state, step5, step6, {0, 20, 10}
    );
    require(
        next.status == "complete" && next.next_step == 4
            && state.major_iteration == 1,
        "Step 7 did not advance exactly one major iteration"
    );
    const Held2Step7Result again = run_held2_step7(
        state, step5, step6, {0, 20, 10}
    );
    require(
        again.status == "complete" && again.next_step == 4
            && state.major_iteration == 2,
        "Step 7 did not advance exactly one major iteration"
    );
}

void run_held2_step8_checks() {
    Held2Step5Result equivalent_step5;
    equivalent_step5.status = "complete";
    equivalent_step5.reason = "equivalent_member";
    Held2Step6Result insufficient_step6;
    insufficient_step6.status = "complete";
    insufficient_step6.candidates.push_back(
        {7, {0.20}, 1.0, 0.20, 0.0, {}, "manufactured"}
    );
    require(
        held2_insufficient_candidate_recovery_required(
            equivalent_step5.reason, insufficient_step6.candidates.size()
        ),
        "Step-6 equivalent-member stall did not activate Step-5 recovery"
    );
    insufficient_step6.candidates.push_back(
        {9, {0.80}, 1.0, 0.80, 0.0, {}, "manufactured"}
    );
    require(
        !held2_insufficient_candidate_recovery_required(
            equivalent_step5.reason, insufficient_step6.candidates.size()
        ),
        "Step-5 recovery activated for a complete Step-6 candidate set"
    );
    Held2PersistentState feedback_state;
    Held2Step8Result feedback;
    feedback.active_phases.push_back({
        7,
        0.5,
        {0.20},
        {0.80, 0.10, 0.10},
        1.0,
        0.2,
        -1.0,
        100000.0,
        {-1.0, -1.0, -1.0},
        -0.9,
        {0.0},
    });
    require(
        retain_held2_step8_feedback(feedback_state, feedback)
            && !retain_held2_step8_feedback(feedback_state, feedback)
            && feedback_state.M.size() == 1,
        "Step-8 feedback did not distinguish progress from an unchanged set"
    );

    for (std::string_view reason : {
             "problem_67_not_converged",
             "stage_iii_solver_not_converged",
             "stage_iii_active_set_resolve_failed",
             "stage_iii_active_set_balance_failed",
         }) {
        require(
            held2_step8_failure_returns_to_stage_ii(reason),
            "Step-8 rejected acceleration failure became terminal"
        );
    }
    require(
        !held2_step8_failure_returns_to_stage_ii(
            "problem_67_feasibility_indeterminate"
        )
            && !held2_step8_failure_returns_to_stage_ii(
                "stage_iii_pressure_refinement_failed"
            )
            && !held2_step8_failure_returns_to_stage_ii(
                "stage_iii_solve_budget_exhausted"
            ),
        "Step-8 malformed or uncertified evidence entered Stage II"
    );
    auto [prepared, candidates] = stage_iii_fixture();
    const Held2Step8Result result =
        manufactured_step8(prepared, candidates);
    require(
        result.outcome == Held2Step8Outcome::CertifiedFeasible
            && result.nlp->accepted
            && !result.warm_start_used
            && !result.cold_fallback_used
            && result.timing.optimizer_solves == 1
            && result.provider_state_evaluations
                + result.provider_volume_bound_evaluations
                + result.provider_packing_evaluations
                == result.timing.provider_evaluations
            && result.phase_coalescences.empty()
            && result.active_phases.size() == 2
            && result.active_phases[0].stable_id == 7
            && result.active_phases[1].stable_id == 9,
        "Step-8 Eq. (67) solve changed"
    );

    auto provider_queries =
        std::make_shared<std::vector<std::vector<double>>>();
    auto volume_bound_queries =
        std::make_shared<std::vector<std::vector<double>>>();
    Held2Step1Result cache_prepared = step1(
        {0.0, 1.0, -1.0},
        {0.50, 0.25, 0.25},
        [volume_bound_queries](const std::vector<double>& composition) {
            if (std::find(
                    volume_bound_queries->begin(),
                    volume_bound_queries->end(),
                    composition
                ) != volume_bound_queries->end()) {
                throw std::runtime_error(
                    "manufactured Provider received duplicate volume bounds"
                );
            }
            volume_bound_queries->push_back(composition);
            return std::array<double, 2>{0.5, 1.5};
        }
    );
    volume_bound_queries->clear();
    const Held2StateEvaluator duplicate_rejecting_provider =
        [coordinates = *cache_prepared.coordinates, provider_queries](
            const auto& composition, double log_volume
        ) {
            std::vector<double> key = composition;
            key.push_back(log_volume);
            if (std::find(
                    provider_queries->begin(), provider_queries->end(), key
                ) != provider_queries->end()) {
                throw std::runtime_error(
                    "manufactured Provider received a duplicate state"
                );
            }
            provider_queries->push_back(std::move(key));
            return evaluate_manufactured_state_impl(
                coordinates, composition, log_volume
            );
        };
    const Held2PackingFractionEvaluator cache_packing =
        [](const auto& composition, double) {
            return composition.front();
        };
    const int provider_context = 0;
    Held2AlgorithmCache provider_cache;
    const auto run_cached_problem = [
        &cache_prepared,
        &duplicate_rejecting_provider,
        &cache_packing,
        &provider_context,
        &provider_cache
    ](
        const Held2Step6Result& selected,
        double radius = kHeld2Problem67InitialRadius
    ) {
        return run_held2_step8(
            cache_prepared,
            selected,
            duplicate_rejecting_provider,
            cache_packing,
            nullptr,
            radius,
            &provider_cache,
            &provider_context
        );
    };
    const Held2Step8Result deduplicated_provider =
        run_cached_problem(candidates);
    require(
        deduplicated_provider.outcome
                == Held2Step8Outcome::CertifiedFeasible
            && deduplicated_provider.nlp
            && deduplicated_provider.nlp->accepted,
        "Step-8 repeated an exact Provider state evaluation"
    );
    Held2Step6Result intervening_candidates = candidates;
    intervening_candidates.candidates.push_back(
        {11, {0.55}, 1.0, 0.55, 0.0, {}, "manufactured"}
    );
    const Held2Step8Result intervening_problem =
        run_cached_problem(intervening_candidates);
    require(
        intervening_problem.outcome == Held2Step8Outcome::CertifiedFeasible
            && intervening_problem.timing.terminal_reason
                != "cached_problem_67",
        "Step-8 cache hid a changed Problem (67)"
    );
    const Held2Step8Result expanded_problem = run_cached_problem(
        candidates, kHeld2Problem67ExpandedRadius
    );
    require(
        expanded_problem.outcome == Held2Step8Outcome::CertifiedFeasible
            && expanded_problem.timing.optimizer_solves > 0
            && expanded_problem.timing.terminal_reason
                != "cached_problem_67",
        "Step-8 cache hid a changed Problem (67) radius"
    );
    const std::vector<double> original_feed = *cache_prepared.independent_feed;
    cache_prepared.independent_feed = std::vector<double>{0.40};
    const Held2Step8Result changed_feed_problem = run_held2_step8(
        cache_prepared,
        candidates,
        duplicate_rejecting_provider,
        cache_packing,
        nullptr,
        kHeld2Problem67InitialRadius,
        &provider_cache,
        &provider_context
    );
    cache_prepared.independent_feed = original_feed;
    require(
        changed_feed_problem.timing.optimizer_solves > 0
            && changed_feed_problem.timing.terminal_reason
                != "cached_problem_67",
        "Step-8 cache hid a changed feed"
    );
    const Held2Step8Result reused_provider =
        run_cached_problem(candidates);
    require(
        reused_provider.outcome == Held2Step8Outcome::CertifiedFeasible
            && reused_provider.provider_state_evaluations == 0
            && reused_provider.timing.optimizer_solves == 0
            && reused_provider.timing.terminal_reason
                == "cached_problem_67",
        "Step-8 re-solved an already-adjudicated Problem (67)"
    );

    auto context_queries = std::make_shared<std::uint64_t>(0);
    Held2Step1Result context_prepared = cache_prepared;
    context_prepared.volume_bounds = Held2VolumeBoundsEvaluator(
        [](const std::vector<double>&) {
            return std::array<double, 2>{0.5, 1.5};
        }
    );
    Held2StateEvaluator mutable_provider =
        [coordinates = *cache_prepared.coordinates, context_queries](
            const auto& composition, double log_volume
        ) {
            ++*context_queries;
            return evaluate_manufactured_state_impl(
                coordinates, composition, log_volume
            );
        };
    Held2AlgorithmCache context_cache;
    const int original_context = 0;
    const Held2Step8Result original_context_result = run_held2_step8(
        context_prepared,
        candidates,
        mutable_provider,
        cache_packing,
        nullptr,
        kHeld2Problem67InitialRadius,
        &context_cache,
        &original_context
    );
    const std::uint64_t first_context_queries = *context_queries;
    mutable_provider =
        [coordinates = *cache_prepared.coordinates, context_queries](
            const auto& composition, double log_volume
        ) {
            ++*context_queries;
            return evaluate_manufactured_state_impl(
                coordinates, composition, log_volume
            );
        };
    const int changed_context = 0;
    const Held2Step8Result changed_context_result = run_held2_step8(
        context_prepared,
        candidates,
        mutable_provider,
        cache_packing,
        nullptr,
        kHeld2Problem67InitialRadius,
        &context_cache,
        &changed_context
    );
    require(
        original_context_result.outcome
                == Held2Step8Outcome::CertifiedFeasible
            && changed_context_result.outcome
                == Held2Step8Outcome::CertifiedFeasible
            && changed_context_result.timing.optimizer_solves > 0
            && changed_context_result.timing.terminal_reason
                != "cached_problem_67"
            && *context_queries > first_context_queries,
        "Step-8 cache context did not invalidate an in-place Provider change"
    );

    Held2Step6Result strict_two_phase_duplicate = candidates;
    strict_two_phase_duplicate.candidates[1]
        .independent_modified_fractions = {0.200004};
    const Held2Step1Result strict_duplicate_prepared = step1(
        {0.0, 1.0, -1.0},
        held2_lift_independent_fractions(
            *prepared.coordinates, {0.200002}
        ),
        [](const std::vector<double>&) {
            return std::array<double, 2>{0.5, 1.5};
        }
    );
    const Held2Step8Result collapsed_two_phase =
        manufactured_step8(
            strict_duplicate_prepared, strict_two_phase_duplicate
        );
    require(
        collapsed_two_phase.outcome
                == Held2Step8Outcome::InsufficientCandidates
            && collapsed_two_phase.phase_coalescences.size() == 1
            && collapsed_two_phase.phase_coalescences.front().tolerance
                == kHeld2PhaseMerge.atol
            && !collapsed_two_phase.phase_coalescences.front()
                .reduced_solve_accepted,
        "Step-8 two-phase collapse lost strict identity evidence"
    );

    Held2Step6Result numerical_duplicates = candidates;
    numerical_duplicates.candidates.push_back(
        {10, {0.200004}, 1.0, 0.200004, 0.0, {}, "manufactured"}
    );
    const Held2Step8Result duplicate_reduced =
        manufactured_step8(prepared, numerical_duplicates);
    require(
        duplicate_reduced.outcome == Held2Step8Outcome::CertifiedFeasible
            && duplicate_reduced.active_phases.size() == 2
            && duplicate_reduced.timing.optimizer_solves >= 2
            && duplicate_reduced.phase_coalescences.size() == 1
            && duplicate_reduced.phase_coalescences.front()
                .tolerance == kHeld2PaperPhaseCoalescence.atol
            && duplicate_reduced.phase_coalescences.front()
                .reduced_solve_accepted
            && duplicate_reduced.ordinary_balance_inf
                <= kHeld2Stage3ExplicitBalance.atol,
        "Step-8 did not re-solve a numerical duplicate phase set"
    );

    auto recentered_bound_queries =
        std::make_shared<std::vector<double>>();
    const Held2Step1Result composition_bounded_prepared = step1(
        {0.0, 1.0, -1.0},
        {0.50, 0.25, 0.25},
        [recentered_bound_queries](const std::vector<double>& composition) {
            recentered_bound_queries->push_back(composition.front());
            return std::array<double, 2>{0.5, 1.5};
        }
    );
    const Held2Step8Result composition_bounded_reduced =
        manufactured_step8(
            composition_bounded_prepared, numerical_duplicates
        );
    const bool queried_weighted_center = std::any_of(
        recentered_bound_queries->begin(),
        recentered_bound_queries->end(),
        [](double composition) {
            return composition > 0.20 && composition < 0.200004;
        }
    );
    require(
        composition_bounded_reduced.outcome
                == Held2Step8Outcome::CertifiedFeasible
            && queried_weighted_center,
        "Step-8 did not refresh composition-dependent bounds at the"
        " weighted coalescence center"
    );

    auto nonretained_query_count = std::make_shared<int>(0);
    const Held2Step1Result refreshed_nonretained_prepared = step1(
        {0.0, 1.0, -1.0},
        {0.50, 0.25, 0.25},
        [nonretained_query_count](const std::vector<double>& composition) {
            if (composition.front() > 0.5) {
                ++*nonretained_query_count;
                if (*nonretained_query_count > 1) {
                    return std::array<double, 2>{0.5, 1.4};
                }
            }
            return std::array<double, 2>{0.5, 1.5};
        }
    );
    const Held2Step8Result refreshed_nonretained =
        manufactured_step8(
            refreshed_nonretained_prepared, numerical_duplicates
        );
    const std::string nonretained_failure =
        "Step-8 did not refresh and apply bounds for a recentered"
        " non-retained phase (outcome="
        + std::to_string(static_cast<int>(refreshed_nonretained.outcome))
        + ", queries=" + std::to_string(*nonretained_query_count)
        + ", reason=" + refreshed_nonretained.reason + ")";
    require(
        refreshed_nonretained.outcome
                == Held2Step8Outcome::CertifiedFeasible
            && *nonretained_query_count >= 2,
        nonretained_failure.c_str()
    );

    const Held2Step1Result rejecting_bounds_prepared = step1(
        {0.0, 1.0, -1.0},
        {0.50, 0.25, 0.25},
        [](const std::vector<double>& composition) {
            if (composition.front() > 0.20
                && composition.front() < 0.200004) {
                throw std::domain_error("manufactured bound failure");
            }
            return std::array<double, 2>{0.5, 1.5};
        }
    );
    const Held2Step8Result rejected_recenter = manufactured_step8(
        rejecting_bounds_prepared, numerical_duplicates
    );
    require(
        rejected_recenter.outcome == Held2Step8Outcome::Indeterminate
            && rejected_recenter.reason
                == "stage_iii_recenter_volume_bounds_failed"
            && rejected_recenter.phase_coalescences.size() == 1
            && !rejected_recenter.phase_coalescences.front()
                .reduced_solve_accepted,
        "Step-8 did not fail closed when recentered Provider bounds failed"
    );

    const Held2Step1Result excluding_bounds_prepared = step1(
        {0.0, 1.0, -1.0},
        {0.50, 0.25, 0.25},
        [](const std::vector<double>& composition) {
            if (composition.front() > 0.20
                && composition.front() < 0.200004) {
                return std::array<double, 2>{1.1, 1.5};
            }
            return std::array<double, 2>{0.5, 1.5};
        }
    );
    const Held2Step8Result excluded_weighted_center =
        manufactured_step8(
            excluding_bounds_prepared, numerical_duplicates
        );
    require(
        excluded_weighted_center.outcome
                == Held2Step8Outcome::Indeterminate
            && excluded_weighted_center.reason
                == "stage_iii_recenter_volume_out_of_bounds"
            && excluded_weighted_center.phase_coalescences.size() == 1
            && !excluded_weighted_center.phase_coalescences.front()
                .reduced_solve_accepted,
        "Step-8 silently changed an out-of-domain weighted center"
    );

    std::vector<Held2StageIICandidate> budget_candidates;
    std::vector<std::array<double, 2>> budget_bounds;
    for (const Held2MPoint& point : numerical_duplicates.candidates) {
        const auto volume_interval =
            (*prepared.volume_bounds)(point.independent_modified_fractions);
        budget_candidates.push_back({
            point.independent_modified_fractions,
            point.volume,
            std::log(point.volume),
        });
        budget_bounds.push_back({
            std::log(volume_interval[0]), std::log(volume_interval[1]),
        });
    }
    const Held2Problem67Result budget_exhausted = solve_held2_problem67(
        *prepared.coordinates,
        held2_lift_independent_fractions(
            *prepared.coordinates, *prepared.independent_feed
        ),
        budget_candidates,
        [coordinates = *prepared.coordinates](
            const auto& composition, double log_volume
        ) {
            return evaluate_manufactured_state_impl(
                coordinates, composition, log_volume
            );
        },
        budget_bounds,
        {},
        1
    );
    require(
        budget_exhausted.failure_reason
                == "stage_iii_solve_budget_exhausted"
            && budget_exhausted.phase_coalescences.size() == 1
            && !budget_exhausted.phase_coalescences.front()
                .reduced_solve_accepted,
        "Step-8 budget exhaustion lost coalescence evidence"
    );

    require(
        kHeld2PaperPhaseCoalescence.atol
            == kHeld2PaperPotentialRatio.atol,
        "Step-8 phase coalescence is inconsistent with paper convergence"
    );

    candidates.candidates.push_back(
        {10, {0.50}, 1.0, 0.50, 0.0, {}, "manufactured"}
    );
    const Held2Step8Result expanded =
        manufactured_step8(prepared, candidates);
    require(
        expanded.outcome == Held2Step8Outcome::CertifiedFeasible
            && expanded.timing.optimizer_solves == 2
            && expanded.problem_candidate_ids
                == std::vector<std::uint64_t>{7, 9, 10}
            && expanded.candidate_ids
                == std::vector<std::uint64_t>{7, 9}
            && expanded.continuation_variables.size()
                == expanded.candidate_ids.size() * 3
            && expanded.ordinary_balance_inf
                <= kHeld2Stage3ExplicitBalance.atol,
        "Step-8 inactive-phase re-solve lost its active-set state"
    );
    Held2Step6Result batch_retirement_candidates = candidates;
    batch_retirement_candidates.candidates.push_back(
        {11, {0.55}, 1.0, 0.55, 0.0, {}, "manufactured"}
    );
    const Held2Step8Result batch_retired =
        manufactured_step8(prepared, batch_retirement_candidates);
    require(
        batch_retired.outcome == Held2Step8Outcome::CertifiedFeasible
            && batch_retired.timing.optimizer_solves == 2
            && batch_retired.candidate_ids
                == std::vector<std::uint64_t>{7, 9},
        "Step-8 replayed the reduced NLP for every same-certificate retirement"
    );
    Held2Step8Result stable_previous = expanded;
    for (Held2Phase& phase : stable_previous.active_phases) {
        phase.independent_modified_fractions = {
            phase.stable_id == 7 ? 0.20 : 0.80
        };
        phase.volume = 1.0;
    }
    std::uint64_t replay_evaluations = 0;
    const Held2Step8Result replayed = run_held2_step8(
        prepared,
        candidates,
        [coordinates = *prepared.coordinates, &replay_evaluations](
            const auto& composition, double log_volume
        ) {
            ++replay_evaluations;
            return evaluate_manufactured_state_impl(
                coordinates, composition, log_volume
            );
        },
        [&replay_evaluations](const auto& composition, double) {
            ++replay_evaluations;
            return composition.front();
        },
        &stable_previous
    );
    require(
        replayed.outcome == Held2Step8Outcome::CertifiedFeasible
            && replayed.neighborhood_radius
                == stable_previous.neighborhood_radius
            && replayed.timing.terminal_reason
                != "unchanged_problem_67"
            && replay_evaluations > 0,
        "Step-8 reused previous evidence without an exact cache key"
    );
    replay_evaluations = 0;
    const Held2Step8Result expanded_replay = run_held2_step8(
        prepared,
        candidates,
        [coordinates = *prepared.coordinates, &replay_evaluations](
            const auto& composition, double log_volume
        ) {
            ++replay_evaluations;
            return evaluate_manufactured_state_impl(
                coordinates, composition, log_volume
            );
        },
        [&replay_evaluations](const auto& composition, double) {
            ++replay_evaluations;
            return composition.front();
        },
        &stable_previous,
        kHeld2Problem67ExpandedRadius
    );
    require(
        expanded_replay.neighborhood_radius
                == kHeld2Problem67ExpandedRadius
            && expanded_replay.timing.terminal_reason
                != "unchanged_problem_67"
            && replay_evaluations > 0,
        "Step-8 memoization hid an expanded candidate neighborhood"
    );
    Held2Step8Result recentered_previous = stable_previous;
    recentered_previous.active_phases.front()
        .independent_modified_fractions.front() =
            0.20 + kHeld2Problem67InitialRadius;
    replay_evaluations = 0;
    const Held2Step8Result recentered = run_held2_step8(
        prepared,
        candidates,
        [coordinates = *prepared.coordinates, &replay_evaluations](
            const auto& composition, double log_volume
        ) {
            ++replay_evaluations;
            return evaluate_manufactured_state_impl(
                coordinates, composition, log_volume
            );
        },
        [&replay_evaluations](const auto& composition, double) {
            ++replay_evaluations;
            return composition.front();
        },
        &recentered_previous
    );
    require(
        recentered.timing.terminal_reason != "unchanged_problem_67"
            && replay_evaluations > 0,
        "Step-8 memoization ignored a recentered candidate neighborhood"
    );
    Held2Step8Result failed_reduction = expanded;
    failed_reduction.outcome = Held2Step8Outcome::Indeterminate;
    failed_reduction.reason = "problem_67_not_converged";
    failed_reduction.nlp->accepted = false;
    failed_reduction.active_phases.clear();
    failed_reduction.total_reduced_gibbs.reset();
    failed_reduction.candidate_variables = {0.20, 1.0, 0.80, 1.0};
    replay_evaluations = 0;
    const Held2Step8Result failed_replay = run_held2_step8(
        prepared,
        candidates,
        [coordinates = *prepared.coordinates, &replay_evaluations](
            const auto& composition, double log_volume
        ) {
            ++replay_evaluations;
            return evaluate_manufactured_state_impl(
                coordinates, composition, log_volume
            );
        },
        [&replay_evaluations](const auto& composition, double) {
            ++replay_evaluations;
            return composition.front();
        },
        &failed_reduction
    );
    require(
        failed_replay.outcome == Held2Step8Outcome::CertifiedFeasible
            && failed_replay.problem_candidate_ids
                == std::vector<std::uint64_t>{7, 9}
            && failed_replay.timing.terminal_reason
                != "unchanged_problem_67"
            && failed_replay.timing.terminal_status == "complete"
            && replay_evaluations > 0,
        "Step-8 reused indeterminate evidence without adjudication"
    );
    candidates.candidates.push_back(
        {11, {0.55}, 1.0, 0.55, 0.0, {}, "manufactured"}
    );
    const Held2Step8Result continued = run_held2_step8(
        prepared,
        candidates,
        [coordinates = *prepared.coordinates](
            const auto& composition, double log_volume
        ) {
            return evaluate_manufactured_state_impl(
                coordinates, composition, log_volume
            );
        },
        [](const auto& composition, double) {
            return composition.front();
        },
        &expanded
    );
    require(
        continued.outcome == Held2Step8Outcome::CertifiedFeasible
            && continued.warm_start_used
            && continued.problem_candidate_ids
                == std::vector<std::uint64_t>{7, 9, 11}
            && std::find(
                continued.candidate_ids.begin(),
                continued.candidate_ids.end(),
                10
            ) == continued.candidate_ids.end(),
        "Step-8 continuation retained a certified inactive candidate"
    );

    candidates.candidates.resize(2);
    candidates.candidates[0].independent_modified_fractions = {0.10};
    candidates.candidates[1].independent_modified_fractions = {0.20};
    const Held2Step8Result infeasible =
        manufactured_step8(prepared, candidates);
    require(
        infeasible.outcome == Held2Step8Outcome::CertifiedInfeasible
            && infeasible.farkas
            && infeasible.farkas->accepted
            && infeasible.farkas
                ->solver_ray_recovered_without_presolve
            && !infeasible.nlp
            && infeasible.farkas->contradiction_margin
                > infeasible.farkas->contradiction_threshold,
        "Step-8 Eq. (67) infeasibility did not return through Step 7"
    );
    const Held2StateEvaluator infeasible_evaluator =
        [coordinates = *prepared.coordinates](
            const auto& composition, double log_volume
        ) {
            return evaluate_manufactured_state_impl(
                coordinates, composition, log_volume
            );
        };
    const Held2PackingFractionEvaluator infeasible_packing =
        [](const auto& composition, double) {
            return composition.front();
        };
    const int infeasible_context = 0;
    Held2AlgorithmCache infeasible_cache;
    const Held2Step8Result retained_infeasible = run_held2_step8(
        prepared,
        candidates,
        infeasible_evaluator,
        infeasible_packing,
        nullptr,
        kHeld2Problem67InitialRadius,
        &infeasible_cache,
        &infeasible_context
    );
    Held2Step6Result intervening_infeasible_candidates = candidates;
    intervening_infeasible_candidates.candidates.push_back(
        {12, {0.80}, 1.0, 0.80, 0.0, {}, "manufactured"}
    );
    const Held2Step8Result intervening_infeasible = run_held2_step8(
        prepared,
        intervening_infeasible_candidates,
        infeasible_evaluator,
        infeasible_packing,
        nullptr,
        kHeld2Problem67InitialRadius,
        &infeasible_cache,
        &infeasible_context
    );
    Held2Step8Result changed_warm_start = result;
    changed_warm_start.active_phases[0]
        .independent_modified_fractions = {0.10};
    changed_warm_start.active_phases[1]
        .independent_modified_fractions = {0.20};
    changed_warm_start.active_phases[0].volume = 1.0;
    changed_warm_start.active_phases[1].volume = 1.0;
    changed_warm_start.continuation_variables = {
        0.6, 0.11, 0.0,
        0.4, 0.19, 0.0,
    };
    const Held2Step8Result replayed_infeasible = run_held2_step8(
        prepared,
        candidates,
        infeasible_evaluator,
        infeasible_packing,
        &changed_warm_start,
        kHeld2Problem67InitialRadius,
        &infeasible_cache,
        &infeasible_context
    );
    require(
        retained_infeasible.outcome
                == Held2Step8Outcome::CertifiedInfeasible
            && intervening_infeasible.timing.terminal_reason
                != "cached_problem_67"
            && replayed_infeasible.outcome
                == Held2Step8Outcome::CertifiedInfeasible
            && replayed_infeasible.farkas
            && replayed_infeasible.farkas->accepted
            && replayed_infeasible.timing.optimizer_solves == 0
            && replayed_infeasible.timing.terminal_reason
                == "cached_problem_67",
        "Step-8 did not retain exact certified-infeasible evidence"
    );
    Held2AlgorithmResult serialized_infeasible;
    serialized_infeasible.step8_history.push_back(infeasible);
    const std::string serialized_farkas = flash_result_to_json({
        {298.15, 100000.0, {0.5, 0.25, 0.25}},
        "sha256:manufactured",
        serialized_infeasible,
    });
    require(
        serialized_farkas.find(
            "\"solver_ray_recovered_without_presolve\":true"
        ) != std::string::npos
            && serialized_farkas.find(
                "\"stage3_ipopt_target\":"
            ) != std::string::npos,
        "Step-8 Farkas or effective solver tolerance was not serialized"
    );

    candidates.candidates[0].independent_modified_fractions = {0.4998};
    candidates.candidates[1].independent_modified_fractions = {0.5002};
    require(
        manufactured_step8(prepared, candidates).outcome
            == Held2Step8Outcome::InsufficientCandidates,
        "Step-8 collapsed phase set did not return through Step 7"
    );
}

void run_held2_step9_checks() {
    auto [prepared, candidates] = stage_iii_fixture();
    const Held2Step8Result step8 =
        manufactured_step8(prepared, candidates);
    Held2Step4Result step4;
    step4.status = "complete";
    step4.upper_bound = step8.total_reduced_gibbs;
    std::size_t provider_evaluations = 0;
    const auto evaluator = [
        coordinates = *prepared.coordinates,
        &provider_evaluations
    ](
        const auto& composition, double log_volume
    ) {
        ++provider_evaluations;
        return evaluate_manufactured_state_impl(
            coordinates, composition, log_volume
        );
    };
    const Held2Step9Result converged =
        run_held2_step9(step4, step8, evaluator);
    require(
        converged.outcome == Held2Step9Outcome::Converged
            && converged.next_action == Held2Step9Action::RunStep10
            && converged.physical->accepted
            && !converged.potential_comparisons.empty(),
        "Step-9 Eqs. (68)-(69) convergence changed"
    );
    require(
        converged.timing.provider_evaluations == provider_evaluations,
        "Step-9 Provider evaluation accounting changed"
    );
    *step4.upper_bound -= 1.0e-9;
    require(
        run_held2_step9(step4, step8, evaluator).outcome
            == Held2Step9Outcome::Converged,
        "Step-9 rejected a roundoff-scale negative Eq. (68) gap"
    );
    *step4.upper_bound = *step8.total_reduced_gibbs + 5.0e-5;
    require(
        run_held2_step9(step4, step8, evaluator).outcome
            == Held2Step9Outcome::Converged,
        "Step-9 rejected an Eq. (68) gap within the paper tolerance"
    );
    *step4.upper_bound = *step8.total_reduced_gibbs + 2.0e-4;
    const Held2Step9Result free_energy_failure =
        run_held2_step9(step4, step8, evaluator);
    require(
        free_energy_failure.outcome
                == Held2Step9Outcome::PaperConvergenceFailed
            && free_energy_failure.next_action
                == Held2Step9Action::ReturnStageII,
        "Step-9 Eq. (68) failure did not return to Step 4"
    );
    *step4.upper_bound = *step8.total_reduced_gibbs;
    const Held2Step9Result zero_denominator = run_held2_step9(
        step4,
        step8,
        [evaluator](const auto& composition, double log_volume) {
            Held2StateEvaluation state = evaluator(
                composition, log_volume
            );
            state.modified_potentials = {
                0.0, composition.front() < 0.5 ? 0.0 : 1.0,
            };
            return state;
        }
    );
    require(
        zero_denominator.outcome
                == Held2Step9Outcome::PaperConvergenceFailed
            && zero_denominator.next_action
                == Held2Step9Action::RunStep10
            && !zero_denominator.potential_comparisons.back().passed
            && std::isinf(
                zero_denominator.potential_comparisons.back().ratio
            ),
        "Step-9 exact-zero denominator rule changed"
    );
}

void run_held2_step10_checks() {
    auto [prepared, candidates] = stage_iii_fixture();
    Held2Step8Result step8 = manufactured_step8(prepared, candidates);
    Held2Step4Result step4;
    step4.status = "complete";
    step4.upper_bound = step8.total_reduced_gibbs;
    const auto evaluator = [coordinates = *prepared.coordinates](
        const auto& composition, double log_volume
    ) {
        return evaluate_manufactured_state_impl(
            coordinates, composition, log_volume
        );
    };
    const Held2Step9Result step9 =
        run_held2_step9(step4, step8, evaluator);
    const Held2Step10Result no_trace = run_held2_step10(
        prepared, step8, step9, evaluator
    );
    require(
        no_trace.status == "complete"
            && no_trace.reason == "trace_refinement_not_required"
            && no_trace.trace_component_indices.empty()
            && no_trace.refinements.empty()
            && no_trace.final_certificate->accepted,
        "Step-10 no-trace path changed"
    );
    Held2Step9Result rejected_kkt = step9;
    rejected_kkt.physical->accepted = false;
    require(
        run_held2_step10(
            prepared, step8, rejected_kkt, evaluator
        ).reason == "trace_final_certificate_failed",
        "Step-10 accepted a rejected Stage-III KKT certificate"
    );
    const double lower =
        prepared.coordinates->independent_lower_bounds.front();
    Held2Step8Result trace_step8 = step8;
    Held2Phase& trace_phase = trace_step8.active_phases.front();
    trace_phase.independent_modified_fractions.front() = lower;
    trace_phase.physical_fractions_provider_order =
        held2_lift_independent_fractions(
            *prepared.coordinates,
            trace_phase.independent_modified_fractions
        );
    Held2Phase& reference_phase = trace_step8.active_phases.back();
    reference_phase.independent_modified_fractions.front() = 1.0 - lower;
    reference_phase.physical_fractions_provider_order =
        held2_lift_independent_fractions(
            *prepared.coordinates,
            reference_phase.independent_modified_fractions
        );
    const std::size_t provider =
        prepared.coordinates->independent_indices.front();
    const std::size_t retained = static_cast<std::size_t>(std::find(
        prepared.coordinates->retained_indices.begin(),
        prepared.coordinates->retained_indices.end(),
        provider
    ) - prepared.coordinates->retained_indices.begin());
    Held2Step9Result trace_failure = step9;
    trace_failure.outcome = Held2Step9Outcome::PaperConvergenceFailed;
    trace_failure.next_action = Held2Step9Action::RunStep10;
    trace_failure.reason = "paper_potential_convergence_failed";
    const Held2Step10Result refined = run_held2_step10(
        prepared,
        trace_step8,
        trace_failure,
        [evaluator, coordinates = *prepared.coordinates, provider, retained](
            const auto& composition,
            double log_volume
        ) {
            const std::vector<double> physical =
                held2_lift_independent_fractions(
                    coordinates,
                    composition,
                    Held2CompositionDomain::TraceRefinement
                );
            std::vector<double> bounded = composition;
            for (std::size_t index = 0; index < bounded.size(); ++index) {
                bounded[index] = std::max(
                    bounded[index],
                    coordinates.independent_lower_bounds[index]
                );
            }
            Held2StateEvaluation state = evaluator(
                bounded, log_volume
            );
            state.modified_fractions =
                held2_transform_physical_fractions(coordinates, physical);
            state.physical_amounts = physical;
            const double fraction = physical[provider];
            state.modified_potentials.assign(
                coordinates.retained_indices.size(), 1.0
            );
            if (fraction <= 5.0e-10) {
                state.modified_potentials[retained] +=
                    std::log10(fraction / 1.0e-12);
            }
            return state;
        }
    );
    require(
        refined.status == "complete"
            && refined.next_action == Held2Step10Action::Accept
            && refined.reason == "trace_refinement_complete"
            && refined.trace_component_indices
                == std::vector<std::size_t>{provider}
            && refined.refinements.size() == 1
            && std::abs(
                refined.refinements.front().refined_mole_fraction - 1.0e-12
            ) <= 1.0e-18
            && refined.final_certificate->accepted,
        "Step-10 bounded trace root changed"
    );
    const Held2Step10Result non_trace = run_held2_step10(
            prepared, manufactured_step8(prepared, candidates),
            trace_failure, evaluator
        );
    require(
        non_trace.reason == "trace_refinement_not_applicable"
            && non_trace.next_action
                == Held2Step10Action::ReturnStageII,
        "Step-10 accepted a non-trace Step-9 potential failure"
    );
    for (Held2Phase& phase : step8.active_phases) {
        phase.independent_modified_fractions.front() = lower;
        phase.physical_fractions_provider_order =
            held2_lift_independent_fractions(
                *prepared.coordinates,
                phase.independent_modified_fractions
            );
    }
    require(
        run_held2_step10(prepared, step8, step9, evaluator).reason
            == "trace_reference_absent",
        "Step-10 absent trace reference did not fail closed"
    );
}

void run_workflow_check() {
    const std::vector<double> charges{0.0, 1.0, -1.0};
    const Held2Coordinates coordinates = make_held2_coordinates(charges);
    auto provider_queries =
        std::make_shared<std::vector<std::vector<double>>>();
    auto volume_bound_queries = std::make_shared<std::uint64_t>(0);
    auto packing_queries = std::make_shared<std::uint64_t>(0);
    const Held2AlgorithmResult result = run_held2_algorithm(
        {
            {"neutral", "cation", "anion"},
            charges,
            [coordinates, provider_queries](
                const auto& composition, double log_volume
            ) {
                std::vector<double> key = composition;
                key.push_back(log_volume);
                if (std::find(
                        provider_queries->begin(),
                        provider_queries->end(),
                        key
                    ) != provider_queries->end()) {
                    throw std::runtime_error(
                        "HELD2 repeated an exact Provider state across stages"
                    );
                }
                provider_queries->push_back(std::move(key));
                Held2StateEvaluation state = evaluate_manufactured_state_impl(
                    coordinates, composition, log_volume
                );
                state.pressure_stationarity_relative *= -1.0;
                state.pressure_stationarity_derivative_log_volume *= -1.0;
                return state;
            },
            [volume_bound_queries](const auto&) {
                ++*volume_bound_queries;
                return std::array<double, 2>{0.5, 1.5};
            },
            [packing_queries](const auto& composition, double) {
                ++*packing_queries;
                return composition.front();
            },
            std::numeric_limits<double>::quiet_NaN(),
            {},
        },
        {298.15, 100000.0, {0.5, 0.25, 0.25}},
        {200, 20, 10}
    );
    std::string workflow_failure =
        "manufactured HELD2 Steps 1-10 failed at "
        + result.failure_stage + ": " + result.failure_reason;
    if (!result.step5_history.empty()
        && !result.step5_history.back().attempts.empty()
        && result.step5_history.back().attempts.back().kkt) {
        workflow_failure += " (last KKT: "
            + result.step5_history.back().attempts.back().kkt->reason + ")";
    }
    require(
        result.outcome == "physical_equilibrium_accepted"
            && !provider_queries->empty()
            && result.step10
            && result.step10->status == "complete"
            && result.phases.size() == 2
            && result.step10->final_certificate->accepted
            && result.globality_certificate == "not_guaranteed",
        workflow_failure.c_str()
    );
    const std::uint64_t reported_provider_work = std::accumulate(
        result.step_timings.begin(), result.step_timings.end(),
        std::uint64_t{0},
        [](std::uint64_t total, const Held2StepTiming& timing) {
            return total + timing.provider_evaluations;
        }
    );
    require(
        reported_provider_work
            == provider_queries->size()
                + *volume_bound_queries + *packing_queries,
        "HELD2 diagnostics counted cache hits as Provider work"
    );
}

}  // namespace
}  // namespace epcsaft_equilibrium::test

int main() {
    using namespace epcsaft_equilibrium::test;
    run_held2_step1_checks();
    run_held2_step2_checks();
    run_held2_step3_checks();
    run_held2_step4_checks();
    run_held2_step5_checks();
    run_held2_step6_checks();
    run_held2_step7_checks();
    run_held2_step8_checks();
    run_held2_step9_checks();
    run_held2_step10_checks();
    run_workflow_check();
    return 0;
}
