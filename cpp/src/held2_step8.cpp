#include "held2_step8.hpp"
#include "held2_tolerances.hpp"

#include <Highs.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

namespace epcsaft_equilibrium {
namespace {

constexpr double kCandidateRadius = 1.0e-3;

struct FeasibilityRun {
    Held2FeasibilityCertificate certificate;
    std::vector<double> initial_variables;
};

FeasibilityRun solve_feasibility(
    const Held2Coordinates& coordinates,
    const std::vector<double>& feed,
    const std::vector<Held2MPoint>& candidates
) {
    const std::size_t phases = candidates.size();
    const std::size_t dimension = feed.size();
    const std::size_t columns = phases * (dimension + 1);
    std::vector<std::vector<double>> rows;
    std::vector<double> rhs;
    std::vector<bool> equality;
    const auto add_row = [&](std::vector<double> row, double bound, bool equal) {
        rows.push_back(std::move(row));
        rhs.push_back(bound);
        equality.push_back(equal);
    };
    std::vector<double> sum_phi(columns);
    for (std::size_t phase = 0; phase < phases; ++phase) {
        sum_phi[phase * (dimension + 1)] = 1.0;
    }
    add_row(std::move(sum_phi), 1.0, true);
    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
        std::vector<double> balance(columns);
        for (std::size_t phase = 0; phase < phases; ++phase) {
            balance[phase * (dimension + 1) + 1 + coordinate] = 1.0;
        }
        add_row(std::move(balance), feed[coordinate], true);
    }
    for (std::size_t phase = 0; phase < phases; ++phase) {
        const std::size_t offset = phase * (dimension + 1);
        std::vector<double> nonnegative(columns);
        nonnegative[offset] = -1.0;
        add_row(std::move(nonnegative), 0.0, false);
        for (const Held2PolytopeConstraint& constraint :
             coordinates.polytope_constraints) {
            std::vector<double> row(columns);
            row[offset] = -constraint.upper_bound;
            for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
                row[offset + 1 + coordinate] =
                    constraint.coefficients[coordinate];
            }
            add_row(std::move(row), 0.0, false);
        }
        for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
            const double center =
                candidates[phase].independent_modified_fractions[coordinate];
            const double lower = std::max(
                coordinates.independent_lower_bounds[coordinate],
                center - kCandidateRadius
            );
            const double upper = std::min(
                coordinates.independent_upper_bounds[coordinate],
                center + kCandidateRadius
            );
            std::vector<double> upper_row(columns);
            upper_row[offset] = -upper;
            upper_row[offset + 1 + coordinate] = 1.0;
            add_row(std::move(upper_row), 0.0, false);
            std::vector<double> lower_row(columns);
            lower_row[offset] = lower;
            lower_row[offset + 1 + coordinate] = -1.0;
            add_row(std::move(lower_row), 0.0, false);
        }
    }

    HighsModel model;
    model.lp_.num_col_ = static_cast<HighsInt>(columns);
    model.lp_.num_row_ = static_cast<HighsInt>(rows.size());
    model.lp_.col_cost_.assign(columns, 0.0);
    model.lp_.col_lower_.assign(columns, -kHighsInf);
    model.lp_.col_upper_.assign(columns, kHighsInf);
    model.lp_.row_lower_.reserve(rows.size());
    model.lp_.row_upper_.reserve(rows.size());
    for (std::size_t row = 0; row < rows.size(); ++row) {
        model.lp_.row_lower_.push_back(equality[row] ? rhs[row] : -kHighsInf);
        model.lp_.row_upper_.push_back(rhs[row]);
    }
    model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
    model.lp_.a_matrix_.start_ = {0};
    for (std::size_t column = 0; column < columns; ++column) {
        for (std::size_t row = 0; row < rows.size(); ++row) {
            if (rows[row][column] != 0.0) {
                model.lp_.a_matrix_.index_.push_back(
                    static_cast<HighsInt>(row)
                );
                model.lp_.a_matrix_.value_.push_back(rows[row][column]);
            }
        }
        model.lp_.a_matrix_.start_.push_back(
            static_cast<HighsInt>(model.lp_.a_matrix_.index_.size())
        );
    }

    FeasibilityRun result;
    Highs highs;
    if (highs.setOptionValue("output_flag", false) == HighsStatus::kError
        || highs.setOptionValue("threads", 1) == HighsStatus::kError
        || highs.passModel(model) == HighsStatus::kError
        || highs.run() == HighsStatus::kError) {
        result.certificate.solver_status = "setup_or_solve_failed";
        return result;
    }
    const HighsModelStatus status = highs.getModelStatus();
    result.certificate.solver_status = highs.modelStatusToString(status);
    if (status == HighsModelStatus::kOptimal) {
        const HighsSolution& solution = highs.getSolution();
        double residual = 0.0;
        for (std::size_t row = 0; row < rows.size(); ++row) {
            double value = 0.0;
            for (std::size_t column = 0; column < columns; ++column) {
                value += rows[row][column] * solution.col_value[column];
            }
            residual = std::max(
                residual,
                equality[row] ? std::abs(value - rhs[row])
                              : std::max(0.0, value - rhs[row])
            );
        }
        result.certificate.primal_residual_inf = residual;
        result.certificate.feasible = audit_held2_tolerance(
            kHeld2LpPrimal, residual
        ).passed;
        if (result.certificate.feasible) {
            result.initial_variables.reserve(phases * (dimension + 2));
            for (std::size_t phase = 0; phase < phases; ++phase) {
                const std::size_t offset = phase * (dimension + 1);
                const double fraction = solution.col_value[offset];
                result.initial_variables.push_back(fraction);
                for (std::size_t coordinate = 0;
                     coordinate < dimension;
                     ++coordinate) {
                    result.initial_variables.push_back(
                        fraction > kHeld2LpPrimal.atol
                            ? solution.col_value[offset + 1 + coordinate]
                                / fraction
                            : candidates[phase]
                                .independent_modified_fractions[coordinate]
                    );
                }
                result.initial_variables.push_back(
                    std::log(candidates[phase].volume)
                );
            }
        }
        return result;
    }
    if (status != HighsModelStatus::kInfeasible) {
        return result;
    }
    std::vector<double> ray(rows.size());
    bool available = false;
    if (highs.getDualRay(available, ray.data()) == HighsStatus::kError
        || !available) {
        return result;
    }
    for (double orientation : {1.0, -1.0}) {
        double scale = 0.0;
        for (double value : ray) scale = std::max(scale, std::abs(value));
        if (scale == 0.0) continue;
        double residual = 0.0;
        double contradiction = 0.0;
        for (std::size_t row = 0; row < rows.size(); ++row) {
            const double multiplier = orientation * ray[row] / scale;
            if (!equality[row]) residual = std::max(residual, -multiplier);
            contradiction += multiplier * rhs[row];
        }
        for (std::size_t column = 0; column < columns; ++column) {
            double stationarity = 0.0;
            for (std::size_t row = 0; row < rows.size(); ++row) {
                stationarity += orientation * ray[row] / scale
                    * rows[row][column];
            }
            residual = std::max(residual, std::abs(stationarity));
        }
        residual = std::max(residual, std::max(0.0, contradiction));
        if (contradiction < -1.0e-9 && residual <= 1.0e-7) {
            result.certificate.infeasible = true;
            result.certificate.farkas_certificate_valid = true;
            result.certificate.certificate_residual_inf = residual;
            break;
        }
    }
    return result;
}

std::vector<double> independent(
    const Held2Coordinates& coordinates,
    const std::vector<double>& modified
) {
    std::vector<double> result;
    for (std::size_t provider : coordinates.independent_indices) {
        const auto position = std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            provider
        );
        result.push_back(modified[static_cast<std::size_t>(
            position - coordinates.retained_indices.begin()
        )]);
    }
    return result;
}

std::uint64_t nearest_id(
    const std::vector<Held2MPoint>& candidates,
    const std::vector<double>& composition
) {
    return std::min_element(
        candidates.begin(),
        candidates.end(),
        [&](const Held2MPoint& left, const Held2MPoint& right) {
            const auto distance = [&](const Held2MPoint& point) {
                double value = 0.0;
                for (std::size_t index = 0; index < composition.size(); ++index) {
                    value = std::max(value, std::abs(
                        point.independent_modified_fractions[index]
                        - composition[index]
                    ));
                }
                return value;
            };
            const double left_distance = distance(left);
            const double right_distance = distance(right);
            return left_distance == right_distance
                ? left.insertion_id < right.insertion_id
                : left_distance < right_distance;
        }
    )->insertion_id;
}

}  // namespace

Held2Step8Result run_held2_step8(
    const Held2Step1Result& step1,
    const Held2Step6Result& step6,
    const Held2StateEvaluator& evaluator,
    const Held2PackingFractionEvaluator& packing_fraction,
    Held2ProgressObserver* observer
) {
    Held2Step8Result result;
    result.timing.invocation_count = 1;
    if (!step1.coordinates || !step1.independent_feed || !step1.volume_bounds
        || step6.status != "complete" || step6.candidates.size() < 2
        || !evaluator || !packing_fraction) {
        result.reason = "invalid_step8_input";
        return result;
    }
    const FeasibilityRun feasible = solve_feasibility(
        *step1.coordinates, *step1.independent_feed, step6.candidates
    );
    result.feasibility = feasible.certificate;
    result.timing.optimizer_solves = 1;
    if (feasible.certificate.infeasible) {
        result.outcome = Held2Step8Outcome::CertifiedInfeasible;
        result.reason = "candidate_neighborhoods_infeasible";
        result.timing.terminal_status = "complete";
        result.timing.terminal_reason = result.reason;
        result.timing.next_step = 7;
        return result;
    }
    if (!feasible.certificate.feasible) {
        result.reason = "candidate_feasibility_uncertified";
        return result;
    }

    std::vector<Held2StageIICandidate> candidates;
    std::vector<std::array<double, 2>> bounds;
    std::uint64_t provider_evaluations = 0;
    for (const Held2MPoint& point : step6.candidates) {
        const std::array<double, 2> physical_bounds =
            (*step1.volume_bounds)(point.independent_modified_fractions);
        ++provider_evaluations;
        candidates.push_back({
            {},
            point.independent_modified_fractions,
            point.volume,
            std::log(point.volume),
            point.reduced_gibbs,
        });
        bounds.push_back({
            std::log(physical_bounds[0]), std::log(physical_bounds[1]),
        });
    }
    const auto feasibility_solve_count = std::make_shared<std::uint64_t>(1);
    const Held2StageIIIInitializer initializer =
        [coordinates = *step1.coordinates,
         feed = *step1.independent_feed,
         initial = feasible.initial_variables,
         first = true,
         feasibility_solve_count](
            const std::vector<Held2StageIICandidate>& candidates
        ) mutable {
            if (first) {
                first = false;
                return Held2StageIIIFeasibilityStart{
                    "feasible", std::move(initial)
                };
            }
            std::vector<Held2MPoint> points;
            points.reserve(candidates.size());
            for (std::size_t index = 0; index < candidates.size(); ++index) {
                points.push_back({
                    index,
                    candidates[index].independent_modified_fractions,
                    candidates[index].volume,
                    std::numeric_limits<double>::quiet_NaN(),
                    0.0,
                    "step8_active_set",
                });
            }
            ++*feasibility_solve_count;
            const FeasibilityRun run =
                solve_feasibility(coordinates, feed, points);
            return Held2StageIIIFeasibilityStart{
                run.certificate.feasible
                    ? "feasible"
                    : run.certificate.infeasible
                        ? "infeasible" : "uncertified",
                run.initial_variables,
            };
        };
    const Held2StateEvaluator counted_evaluator =
        [&evaluator, &provider_evaluations](
            const std::vector<double>& composition,
            double log_volume
        ) {
            ++provider_evaluations;
            return evaluator(composition, log_volume);
        };
    const Held2PackingFractionEvaluator counted_packing =
        [&packing_fraction, &provider_evaluations](
            const std::vector<double>& composition,
            double volume
        ) {
            ++provider_evaluations;
            return packing_fraction(composition, volume);
        };
    const Held2StageIIIResult solved = solve_held2_stage_iii(
        *step1.coordinates,
        held2_lift_independent_fractions(
            *step1.coordinates, *step1.independent_feed
        ),
        candidates,
        counted_evaluator,
        bounds,
        std::numeric_limits<double>::quiet_NaN(),
        "unavailable",
        observer,
        initializer,
        counted_packing
    );
    result.timing.provider_evaluations = provider_evaluations;
    result.timing.optimizer_solves = *feasibility_solve_count
        + static_cast<std::uint64_t>(solved.stage_iii_solve_count);
    result.timing.optimizer_iterations =
        static_cast<std::uint64_t>(solved.optimizer_iteration_count);
    const double primal = std::max(
        solved.modified_balance_inf_norm, solved.ordinary_balance_inf_norm
    );
    const bool accepted =
        solved.numerical_status == "converged"
        && solved.phases.size() >= 2
        && audit_held2_tolerance(kHeld2Stage3ModifiedBalance, primal).passed
        && audit_held2_tolerance(
            kHeld2Stage3Charge,
            solved.phase_charge_inf_norm,
            solved.phase_charge_scale
        ).passed
        && audit_held2_tolerance(
            kHeld2Stage3Pressure, solved.pressure_stationarity_inf_norm
        ).passed
        && audit_held2_tolerance(
            kHeld2Stage3Stationarity, solved.kkt_stationarity_inf_norm
        ).passed
        && audit_held2_tolerance(
            kHeld2Stage3DualSign, solved.dual_sign_violation_inf_norm
        ).passed
        && audit_held2_tolerance(
            kHeld2Stage3Complementarity,
            solved.bound_complementarity_inf_norm
        ).passed;
    result.nlp = Held2NlpCertificate{
        solved.solver_status,
        primal,
        solved.kkt_stationarity_inf_norm,
        solved.dual_sign_violation_inf_norm,
        solved.bound_complementarity_inf_norm,
        accepted,
    };
    result.ordinary_balance_inf = solved.ordinary_balance_inf_norm;
    result.electroneutrality_inf = solved.phase_charge_inf_norm;
    result.electroneutrality_scale = solved.phase_charge_scale;
    result.pressure_residual_inf = solved.pressure_stationarity_inf_norm;
    for (const Held2StageIIILifecycleStep& decision : solved.lifecycle) {
        std::uint64_t id = 0;
        if (!decision.candidate_independent_modified_fractions.empty()) {
            id = nearest_id(
                step6.candidates,
                decision.candidate_independent_modified_fractions
            );
        }
        result.lifecycle.push_back({
            id,
            decision.action,
            decision.decision_reason,
            decision.action == "accept_active_set",
        });
    }
    if (!accepted) {
        if (solved.failure_reason == "collapsed_phase_set") {
            result.outcome = Held2Step8Outcome::InsufficientCandidates;
            result.reason = solved.failure_reason;
            result.timing.terminal_status = "complete";
            result.timing.terminal_reason = result.reason;
            result.timing.next_step = 7;
            return result;
        }
        result.reason = solved.failure_reason.empty()
            ? "stage8_nlp_uncertified" : solved.failure_reason;
        return result;
    }
    for (const Held2StageIIIPhase& phase : solved.phases) {
        const std::vector<double> composition =
            independent(*step1.coordinates, phase.modified_fractions);
        result.active_phases.push_back({
            nearest_id(step6.candidates, composition),
            phase.phase_fraction,
            composition,
            phase.physical_fractions,
            phase.volume,
            counted_packing(composition, phase.volume),
        });
    }
    result.timing.provider_evaluations = provider_evaluations;
    std::sort(
        result.active_phases.begin(),
        result.active_phases.end(),
        [](const Held2Phase& left, const Held2Phase& right) {
            return left.stable_id < right.stable_id;
        }
    );
    result.outcome = Held2Step8Outcome::CertifiedFeasible;
    result.reason = "step8_complete";
    result.total_reduced_gibbs = solved.objective;
    result.timing.terminal_status = "complete";
    result.timing.terminal_reason = result.reason;
    result.timing.next_step = 9;
    return result;
}

}  // namespace epcsaft_equilibrium
