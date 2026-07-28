#include "held2_step10.hpp"
#include "held2_tolerances.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace epcsaft_equilibrium {
namespace {

constexpr double kTraceMaximum = 5.0e-10;
constexpr double kLogTraceMinimum = -300.0;
constexpr int kBisectionIterations = 100;

std::size_t retained_position(
    const Held2Coordinates& coordinates,
    std::size_t provider
) {
    return static_cast<std::size_t>(std::find(
        coordinates.retained_indices.begin(),
        coordinates.retained_indices.end(),
        provider
    ) - coordinates.retained_indices.begin());
}

double charge_residual(
    const std::vector<double>& charges,
    const std::vector<double>& fractions
) {
    double residual = 0.0;
    for (std::size_t index = 0; index < charges.size(); ++index) {
        residual += charges[index] * fractions[index];
    }
    return residual;
}

double maximum_abs_difference(
    const std::vector<double>& left,
    const std::vector<double>& right
) {
    double result = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        result = std::max(result, std::abs(left[index] - right[index]));
    }
    return result;
}

}  // namespace

Held2Step10Result run_held2_step10(
    const Held2Step1Result& step1,
    const Held2Step8Result& step8,
    const Held2Step9Result& step9,
    const Held2StateEvaluator& evaluator,
    Held2ProgressObserver* observer
) {
    Held2Step10Result result;
    result.timing.invocation_count = 1;
    const bool trace_potential_failure =
        step9.outcome == Held2Step9Outcome::PaperConvergenceFailed;
    if (!step1.coordinates || !step1.independent_feed
        || step8.outcome != Held2Step8Outcome::CertifiedFeasible
        || step9.next_action != Held2Step9Action::RunStep10
        || !step9.physical || !evaluator) {
        result.reason = "invalid_step10_input";
        return result;
    }
    const Held2Coordinates& coordinates = *step1.coordinates;
    result.phases = step8.active_phases;
    if (trace_potential_failure
        && std::any_of(
            step9.potential_comparisons.begin(),
            step9.potential_comparisons.end(),
            [&](const Held2PotentialComparison& comparison) {
                return !comparison.passed
                    && (comparison.component_index
                            >= coordinates.retained_indices.size()
                        || coordinates.charges[
                            coordinates.retained_indices[
                                comparison.component_index
                            ]
                        ] == 0.0);
            }
        )) {
        result.next_action = Held2Step10Action::ReturnStageII;
        result.reason = "trace_refinement_not_applicable";
        return result;
    }
    bool trace_found = false;
    for (Held2Phase& phase : result.phases) {
        const std::vector<double> modified =
            held2_transform_physical_fractions(
                coordinates, phase.physical_fractions_provider_order
            );
        for (std::size_t coordinate = 0;
             coordinate < coordinates.independent_indices.size();
             ++coordinate) {
            const std::size_t provider =
                coordinates.independent_indices[coordinate];
            const std::size_t retained =
                retained_position(coordinates, provider);
            const double lower = kHeld2ModifiedLowerScale
                * coordinates.modified_factors[retained];
            if (modified[retained] > 10.0 * lower) continue;
            trace_found = true;
            const auto reference = std::max_element(
                result.phases.begin(),
                result.phases.end(),
                [&](const Held2Phase& left, const Held2Phase& right) {
                    const double left_amount =
                        left.physical_fractions_provider_order[provider];
                    const double right_amount =
                        right.physical_fractions_provider_order[provider];
                    return left_amount == right_amount
                        ? left.stable_id > right.stable_id
                        : left_amount < right_amount;
                }
            );
            if (reference == result.phases.end()
                || reference->physical_fractions_provider_order[provider]
                    <= kTraceMaximum) {
                result.reason = "trace_reference_absent";
                return result;
            }
            const Held2StateEvaluation reference_state = evaluator(
                reference->independent_modified_fractions,
                std::log(reference->volume)
            );
            ++result.timing.provider_evaluations;
            if (retained >= reference_state.modified_potentials.size()) {
                result.reason = "trace_potential_evidence_missing";
                return result;
            }
            const double target =
                reference_state.modified_potentials[retained];
            const auto residual = [&](double log10_fraction) {
                std::vector<double> trial =
                    phase.independent_modified_fractions;
                trial[coordinate] = coordinates.modified_factors[retained]
                    * std::pow(10.0, log10_fraction);
                const Held2StateEvaluation state = evaluator(
                    trial, std::log(phase.volume)
                );
                ++result.timing.provider_evaluations;
                if (retained >= state.modified_potentials.size()) {
                    throw std::invalid_argument(
                        "trace_potential_evidence_missing"
                    );
                }
                return state.modified_potentials[retained] - target;
            };
            double lower_log = kLogTraceMinimum;
            double upper_log = std::log10(kTraceMaximum);
            double lower_residual;
            double upper_residual;
            try {
                lower_residual = residual(lower_log);
                upper_residual = residual(upper_log);
            } catch (const std::exception&) {
                result.reason = "trace_provider_evaluation_failed";
                return result;
            }
            if (!std::isfinite(lower_residual)
                || !std::isfinite(upper_residual)
                || lower_residual * upper_residual > 0.0) {
                result.reason = "trace_residual_not_bracketed";
                return result;
            }
            double root = std::abs(lower_residual) <= std::abs(upper_residual)
                ? lower_log : upper_log;
            double root_residual = std::abs(lower_residual)
                <= std::abs(upper_residual)
                ? lower_residual : upper_residual;
            for (int iteration = 0;
                 iteration < kBisectionIterations
                    && std::abs(root_residual) > 1.0e-8;
                 ++iteration) {
                const double middle = 0.5 * (lower_log + upper_log);
                const double middle_residual = residual(middle);
                if (!std::isfinite(middle_residual)) {
                    result.reason = "trace_provider_evaluation_failed";
                    return result;
                }
                root = middle;
                root_residual = middle_residual;
                if ((lower_residual <= 0.0 && middle_residual >= 0.0)
                    || (lower_residual >= 0.0 && middle_residual <= 0.0)) {
                    upper_log = middle;
                } else {
                    lower_log = middle;
                    lower_residual = middle_residual;
                }
            }
            if (std::abs(root_residual) > 1.0e-8) {
                result.reason = "trace_residual_unresolved";
                return result;
            }
            const double initial =
                phase.physical_fractions_provider_order[provider];
            phase.independent_modified_fractions[coordinate] =
                coordinates.modified_factors[retained] * std::pow(10.0, root);
            try {
                phase.physical_fractions_provider_order =
                    held2_lift_independent_fractions(
                        coordinates,
                        phase.independent_modified_fractions,
                        Held2CompositionDomain::TraceRefinement
                    );
            } catch (const std::invalid_argument&) {
                result.reason = "trace_reconstruction_failed";
                return result;
            }
            result.refinements.push_back({
                phase.stable_id,
                provider,
                reference->stable_id,
                initial,
                phase.physical_fractions_provider_order[provider],
                root_residual,
                "refined",
            });
        }
    }
    if (!trace_found && trace_potential_failure) {
        result.next_action = Held2Step10Action::ReturnStageII;
        result.reason = "trace_refinement_not_applicable";
        return result;
    }

    std::vector<double> material_balance(coordinates.charges.size(), 0.0);
    std::vector<double> modified_balance(
        coordinates.retained_indices.size(), 0.0
    );
    double charge = 0.0;
    double pressure = 0.0;
    double total_free_energy = 0.0;
    for (Held2Phase& phase : result.phases) {
        for (std::size_t component = 0;
             component < material_balance.size();
             ++component) {
            material_balance[component] += phase.phase_fraction
                * phase.physical_fractions_provider_order[component];
        }
        const std::vector<double> modified =
            held2_transform_physical_fractions(
                coordinates, phase.physical_fractions_provider_order
            );
        for (std::size_t component = 0;
             component < modified_balance.size();
             ++component) {
            modified_balance[component] +=
                phase.phase_fraction * modified[component];
        }
        charge = std::max(
            charge,
            std::abs(charge_residual(
                coordinates.charges,
                phase.physical_fractions_provider_order
            ))
        );
        const Held2StateEvaluation state = evaluator(
            phase.independent_modified_fractions,
            std::log(phase.volume)
        );
        ++result.timing.provider_evaluations;
        phase.helmholtz_over_rt_reference_amount =
            state.helmholtz_over_rt_reference_amount;
        phase.pressure_pa = state.pressure_pa;
        phase.chemical_potentials_over_rt =
            state.chemical_potentials_over_rt;
        total_free_energy += phase.phase_fraction * state.objective;
        pressure = std::max(
            pressure, std::abs(state.pressure_stationarity_relative)
        );
    }
    result.total_free_energy_over_rt = total_free_energy;
    const std::vector<double> feed = held2_lift_independent_fractions(
        coordinates, *step1.independent_feed
    );
    const std::vector<double> modified_feed =
        held2_transform_physical_fractions(coordinates, feed);
    Held2PhysicalCertificate certificate = *step9.physical;
    certificate.ordinary_balance_inf =
        maximum_abs_difference(material_balance, feed);
    certificate.modified_balance_inf =
        maximum_abs_difference(modified_balance, modified_feed);
    certificate.electroneutrality_inf = charge;
    certificate.pressure_residual_inf = pressure;
    certificate.accepted = step9.physical->accepted
        && audit_held2_tolerance(
        kHeld2Stage3ModifiedBalance, certificate.modified_balance_inf
    ).passed && audit_held2_tolerance(
        kHeld2Stage3ExplicitBalance, certificate.ordinary_balance_inf
    ).passed && audit_held2_tolerance(
        kHeld2Stage3Charge, certificate.electroneutrality_inf
    ).passed && audit_held2_tolerance(
        kHeld2Stage3Pressure, certificate.pressure_residual_inf
    ).passed && audit_held2_tolerance(
        kHeld2Stage3Stationarity, certificate.kkt_residual_inf
    ).passed;
    result.final_certificate = certificate;
    if (!certificate.accepted) {
        result.reason = "trace_final_certificate_failed";
        return result;
    }
    Held2Step8Result refined = step8;
    refined.active_phases = result.phases;
    refined.ordinary_balance_inf = certificate.ordinary_balance_inf;
    refined.electroneutrality_inf = certificate.electroneutrality_inf;
    refined.pressure_residual_inf = certificate.pressure_residual_inf;
    Held2Step4Result bound;
    bound.status = "complete";
    if (!step9.free_energy_gap || !step8.total_reduced_gibbs) {
        result.reason = "trace_free_energy_evidence_missing";
        return result;
    }
    bound.upper_bound =
        *step8.total_reduced_gibbs + *step9.free_energy_gap;
    const Held2Step9Result recertified =
        run_held2_step9(bound, refined, evaluator);
    result.timing.provider_evaluations += result.phases.size();
    if (recertified.outcome != Held2Step9Outcome::Converged) {
        result.final_certificate->accepted = false;
        result.reason = "trace_potential_recertification_failed";
        return result;
    }
    result.status = "complete";
    result.next_action = Held2Step10Action::Accept;
    result.reason = trace_found ? "trace_refinement_complete"
                                : "trace_refinement_not_required";
    result.timing.terminal_status = result.status;
    result.timing.terminal_reason = result.reason;
    result.timing.next_step = 0;
    Held2ProgressEvent progress;
    progress.kind = Held2ProgressKind::Certificate;
    progress.stage = "STEP 10";
    progress.status = "passed";
    progress.count = static_cast<int>(result.refinements.size());
    observe_held2(observer, progress);
    return result;
}

}  // namespace epcsaft_equilibrium
