#include "held2_stage_ii_basin.hpp"
#include "held2_tolerances.hpp"

#include <nlopt.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace epcsaft_equilibrium {
namespace {

std::string nlopt_version_string() {
    int major = 0;
    int minor = 0;
    int bugfix = 0;
    nlopt_version(&major, &minor, &bugfix);
    return std::to_string(major) + "." + std::to_string(minor) + "."
        + std::to_string(bugfix);
}

double maximum_abs_difference(
    const std::vector<double>& left,
    const std::vector<double>& right
) {
    if (left.size() != right.size()) {
        return std::numeric_limits<double>::infinity();
    }
    double maximum = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        maximum = std::max(maximum, std::abs(left[index] - right[index]));
    }
    return maximum;
}

bool same_composition(
    const std::vector<double>& left,
    const std::vector<double>& right
) {
    return maximum_abs_difference(left, right)
        <= kHeld2BasinDuplicateComposition.atol;
}

bool same_physical_start(
    const Held2StageIIPhysicalStart& left,
    const Held2StageIIPhysicalStart& right
) {
    return same_composition(
               left.independent_modified_fractions,
               right.independent_modified_fractions
           )
        && std::abs(left.log_volume - right.log_volume)
            <= kHeld2BasinDuplicateLogVolume.atol;
}

Held2StageIIBasinEvaluation evaluate_fail_closed(
    const Held2StageIIBasinEvaluator& evaluator,
    const std::vector<double>& independent
) {
    try {
        return evaluator(independent);
    } catch (const std::exception& error) {
        Held2StageIIBasinEvaluation failed;
        failed.independent_modified_fractions = independent;
        failed.failure_reason = std::string("envelope_evaluator_exception: ")
            + error.what();
        return failed;
    }
}

void retain_evaluation(
    Held2StageIIBasinExplorationResult& result,
    Held2StageIIBasinEvaluation evaluation,
    const std::string& source
) {
    if (!evaluation.certified
        || !std::isfinite(evaluation.reduced_lower_value)) {
        ++result.failed_evaluation_count;
        result.evaluations.push_back(std::move(evaluation));
        return;
    }
    ++result.completed_evaluation_count;
    double minimum_stable_root_objective =
        std::numeric_limits<double>::infinity();
    for (const Held2PressureRoot& root : evaluation.pressure_envelope.roots) {
        if (root.mechanical_class == "strict_stable" && !root.boundary) {
            minimum_stable_root_objective = std::min(
                minimum_stable_root_objective,
                root.objective
            );
        }
    }
    const double composition_offset = evaluation.reduced_lower_value
        - minimum_stable_root_objective;
    for (const Held2PressureRoot& root : evaluation.pressure_envelope.roots) {
        if (root.mechanical_class != "strict_stable" || root.boundary) {
            continue;
        }
        Held2StageIIPhysicalStart start;
        start.independent_modified_fractions =
            evaluation.independent_modified_fractions;
        start.log_volume = root.log_volume;
        start.volume = root.volume;
        start.reduced_lower_value = root.objective + composition_offset;
        start.source = source;
        start.root_origin = root.origin;
        start.root_completeness =
            evaluation.pressure_envelope.root_completeness;
        const bool duplicate = std::any_of(
            result.representatives.begin(),
            result.representatives.end(),
            [&start](const Held2StageIIPhysicalStart& known) {
                return same_physical_start(start, known);
            }
        );
        if (duplicate) {
            ++result.duplicate_start_count;
        } else {
            result.representatives.push_back(std::move(start));
        }
    }
    result.evaluations.push_back(std::move(evaluation));
}

struct DirectContext {
    const Held2Coordinates* coordinates = nullptr;
    const Held2StageIIBasinEvaluator* evaluator = nullptr;
    Held2StageIIBasinExplorationResult* result = nullptr;
    nlopt::opt* optimizer = nullptr;
    bool stop_requested = false;
};

double direct_objective(
    const std::vector<double>& cube,
    std::vector<double>& gradient,
    void* opaque
) {
    auto& context = *static_cast<DirectContext*>(opaque);
    if (context.stop_requested) {
        return 0.0;
    }
    if (!gradient.empty()) {
        throw std::invalid_argument(
            "HELD2 Stage-II DIRECT-L requested an unexpected gradient"
        );
    }
    const std::vector<double> independent =
        held2_map_unit_cube_to_independent_fractions(
            *context.coordinates,
            cube,
            std::numeric_limits<double>::quiet_NaN()
        );
    Held2StageIIBasinEvaluation evaluation = evaluate_fail_closed(
        *context.evaluator,
        independent
    );
    const double objective = evaluation.reduced_lower_value;
    retain_evaluation(*context.result, std::move(evaluation), "direct_l");
    if (context.result->failed_evaluation_count != 0) {
        context.stop_requested = true;
        context.result->termination_reason =
            "required_envelope_evaluation_failed";
        context.optimizer->force_stop();
        return 0.0;
    }
    return objective;
}

struct SobolParameters {
    int degree;
    std::uint32_t coefficient;
    std::array<std::uint32_t, 5> initial;
};

constexpr std::array<SobolParameters, 9> kSobolParameters = {{
    {1, 0, {1, 0, 0, 0, 0}},
    {2, 1, {1, 3, 0, 0, 0}},
    {3, 1, {1, 3, 1, 0, 0}},
    {3, 2, {1, 1, 1, 0, 0}},
    {4, 1, {1, 3, 5, 13, 0}},
    {4, 4, {1, 1, 5, 5, 0}},
    {5, 2, {1, 3, 3, 9, 7}},
    {5, 4, {1, 1, 5, 11, 27}},
    {5, 7, {1, 1, 7, 13, 3}},
}};

std::array<std::uint32_t, 32> sobol_directions(std::size_t dimension) {
    std::array<std::uint32_t, 32> directions{};
    if (dimension == 0) {
        for (std::size_t bit = 0; bit < directions.size(); ++bit) {
            directions[bit] = std::uint32_t{1} << (31U - bit);
        }
        return directions;
    }
    if (dimension > kSobolParameters.size()) {
        throw std::invalid_argument(
            "HELD2 Sobol explorer supports at most ten composition coordinates"
        );
    }
    const SobolParameters& parameters = kSobolParameters[dimension - 1];
    for (int bit = 0; bit < parameters.degree; ++bit) {
        directions[static_cast<std::size_t>(bit)] =
            parameters.initial[static_cast<std::size_t>(bit)]
            << (31 - bit);
    }
    for (int bit = parameters.degree; bit < 32; ++bit) {
        std::uint32_t value = directions[static_cast<std::size_t>(
            bit - parameters.degree
        )];
        value ^= value >> parameters.degree;
        for (int offset = 1; offset < parameters.degree; ++offset) {
            const int coefficient_bit = parameters.degree - 1 - offset;
            if ((parameters.coefficient >> coefficient_bit) & 1U) {
                value ^= directions[static_cast<std::size_t>(bit - offset)];
            }
        }
        directions[static_cast<std::size_t>(bit)] = value;
    }
    return directions;
}

}  // namespace

std::vector<std::vector<double>> held2_sobol_points(
    std::size_t dimension,
    int count
) {
    if (dimension == 0 || count < 0 || dimension > 10) {
        throw std::invalid_argument("HELD2 Sobol policy is invalid");
    }
    std::vector<std::array<std::uint32_t, 32>> directions;
    directions.reserve(dimension);
    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
        directions.push_back(sobol_directions(coordinate));
    }
    std::vector<std::vector<double>> points;
    points.reserve(static_cast<std::size_t>(count));
    constexpr double scale = 1.0 / 4294967296.0;
    for (int index = 1; index <= count; ++index) {
        const std::uint32_t gray = static_cast<std::uint32_t>(index)
            ^ (static_cast<std::uint32_t>(index) >> 1U);
        std::vector<double> point(dimension, 0.0);
        for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
            std::uint32_t value = 0;
            for (std::size_t bit = 0; bit < 32; ++bit) {
                if ((gray >> bit) & 1U) {
                    value ^= directions[coordinate][bit];
                }
            }
            point[coordinate] = static_cast<double>(value) * scale;
        }
        points.push_back(std::move(point));
    }
    return points;
}

Held2StageIIBasinExplorationResult explore_held2_stage_ii_basins(
    const Held2Coordinates& coordinates,
    const std::vector<Held2StageIIBasinSeed>& seeds,
    int sobol_count,
    bool use_direct_escalation,
    int direct_evaluation_budget,
    const Held2StageIIBasinEvaluator& evaluator
) {
    if (coordinates.independent_indices.empty() || sobol_count < 0
        || direct_evaluation_budget < 0
        || (use_direct_escalation && direct_evaluation_budget == 0)) {
        throw std::invalid_argument("HELD2 Stage-II basin policy is invalid");
    }
    Held2StageIIBasinExplorationResult result;
    result.direct_solver_version = nlopt_version_string();
    result.declared_sobol_count = sobol_count;
    result.declared_direct_budget = use_direct_escalation
        ? direct_evaluation_budget
        : 0;

    std::vector<Held2StageIIBasinSeed> unique_seeds;
    for (const Held2StageIIBasinSeed& seed : seeds) {
        const bool duplicate = std::any_of(
            unique_seeds.begin(),
            unique_seeds.end(),
            [&seed](const Held2StageIIBasinSeed& known) {
                return same_composition(
                    seed.independent_modified_fractions,
                    known.independent_modified_fractions
                );
            }
        );
        if (duplicate) {
            ++result.duplicate_start_count;
        } else {
            unique_seeds.push_back(seed);
        }
    }
    for (const Held2StageIIBasinSeed& seed : unique_seeds) {
        retain_evaluation(
            result,
            evaluate_fail_closed(
                evaluator,
                seed.independent_modified_fractions
            ),
            seed.source
        );
        if (result.failed_evaluation_count != 0) {
            result.termination_reason = "required_envelope_evaluation_failed";
            return result;
        }
    }
    for (const std::vector<double>& cube : held2_sobol_points(
             coordinates.independent_indices.size(),
             sobol_count
         )) {
        const std::vector<double> independent =
            held2_map_unit_cube_to_independent_fractions(
                coordinates,
                cube,
                std::numeric_limits<double>::quiet_NaN()
            );
        retain_evaluation(
            result,
            evaluate_fail_closed(evaluator, independent),
            "sobol"
        );
        if (result.failed_evaluation_count != 0) {
            result.termination_reason = "required_envelope_evaluation_failed";
            return result;
        }
    }

    if (use_direct_escalation) {
        result.direct_escalation_used = true;
        nlopt::opt optimizer(
            nlopt::GN_DIRECT_L,
            coordinates.independent_indices.size()
        );
        DirectContext context{
            &coordinates,
            &evaluator,
            &result,
            &optimizer,
            false,
        };
        optimizer.set_lower_bounds(std::vector<double>(
            coordinates.independent_indices.size(),
            0.0
        ));
        optimizer.set_upper_bounds(std::vector<double>(
            coordinates.independent_indices.size(),
            1.0
        ));
        optimizer.set_maxeval(direct_evaluation_budget);
        optimizer.set_min_objective(direct_objective, &context);
        std::vector<double> initial(
            coordinates.independent_indices.size(),
            0.5
        );
        double minimum = std::numeric_limits<double>::infinity();
        try {
            const nlopt::result status = optimizer.optimize(initial, minimum);
            if (status != nlopt::MAXEVAL_REACHED) {
                result.termination_reason = "unexpected_direct_termination";
                return result;
            }
        } catch (const nlopt::forced_stop&) {
            if (result.failed_evaluation_count != 0) {
                return result;
            }
            result.termination_reason = "unexpected_direct_forced_stop";
            return result;
        } catch (const std::exception& error) {
            result.termination_reason = std::string("direct_solver_failure: ")
                + error.what();
            return result;
        }
    }

    if (result.representatives.empty()) {
        result.termination_reason = "no_physical_basin_representatives";
        return result;
    }
    std::sort(
        result.representatives.begin(),
        result.representatives.end(),
        [](const Held2StageIIPhysicalStart& left,
           const Held2StageIIPhysicalStart& right) {
            if (left.reduced_lower_value != right.reduced_lower_value) {
                return left.reduced_lower_value < right.reduced_lower_value;
            }
            if (left.independent_modified_fractions
                != right.independent_modified_fractions) {
                return left.independent_modified_fractions
                    < right.independent_modified_fractions;
            }
            return left.log_volume < right.log_volume;
        }
    );
    result.outcome = "representatives_found";
    result.termination_reason = use_direct_escalation
        ? "declared_direct_budget_exhausted"
        : "deterministic_exploration_completed";
    return result;
}

}  // namespace epcsaft_equilibrium
