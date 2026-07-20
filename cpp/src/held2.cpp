#include "held2.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>

namespace epcsaft_equilibrium {
namespace {

constexpr double kCoordinateTolerance = 1.0e-12;
// Perdomo Eq. (61) uses 1e-10 times the eliminated-ion coordinate factor.
constexpr double kModifiedLowerScale = 1.0e-10;
// The approved manufactured oracle uses a strict 1e-8 direct certificate.
constexpr double kCertificateTolerance = 1.0e-8;
constexpr double kStageITpdThreshold = -1.0e-8;
constexpr int kStageIIpoptIterations = 300;
constexpr unsigned int kStageISeed = 2025;

using Held2StateEvaluator = std::function<Held2StateEvaluation(
    const std::vector<double>&,
    double
)>;

void require_finite_vector(const std::vector<double>& values, const char* name) {
    if (!std::all_of(values.begin(), values.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
}

double maximum_abs_difference(
    const std::vector<double>& left,
    const std::vector<double>& right
) {
    if (left.size() != right.size()) {
        throw std::invalid_argument("vectors must have matching sizes");
    }
    double result = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        result = std::max(result, std::abs(left[index] - right[index]));
    }
    return result;
}

double charge_residual(
    const std::vector<double>& charges,
    const std::vector<double>& fractions
) {
    double result = 0.0;
    for (std::size_t index = 0; index < charges.size(); ++index) {
        result += charges[index] * fractions[index];
    }
    return result;
}

double manufactured_helmholtz(double composition, double molar_volume) {
    const double shifted = composition - 0.5;
    constexpr double inner_squared = 0.15 * 0.15;
    constexpr double outer_squared = 0.30 * 0.30;
    const double composition_energy = 1.0e4
        * (std::pow(shifted, 6) / 6.0
           - (inner_squared + outer_squared) * std::pow(shifted, 4) / 4.0
           + inner_squared * outer_squared * shifted * shifted / 2.0);
    const double volume_delta = molar_volume - 1.2;
    return composition_energy + 0.5 * 5.0 * volume_delta * volume_delta;
}

std::array<double, 2> manufactured_gibbs_gradient(
    double composition,
    double molar_volume
) {
    const double shifted = composition - 0.5;
    constexpr double inner_squared = 0.15 * 0.15;
    constexpr double outer_squared = 0.30 * 0.30;
    return {
        1.0e4
            * (std::pow(shifted, 5)
               - (inner_squared + outer_squared) * std::pow(shifted, 3)
               + inner_squared * outer_squared * shifted),
        5.0 * (molar_volume - 1.2) + 1.0,
    };
}

double manufactured_gibbs(double composition, double molar_volume) {
    return manufactured_helmholtz(composition, molar_volume) + molar_volume;
}

double manufactured_composition_hessian(double composition) {
    const double shifted = composition - 0.5;
    constexpr double inner_squared = 0.15 * 0.15;
    constexpr double outer_squared = 0.30 * 0.30;
    return 1.0e4
        * (5.0 * std::pow(shifted, 4)
           - 3.0 * (inner_squared + outer_squared) * shifted * shifted
           + inner_squared * outer_squared);
}

std::array<double, 2> manufactured_modified_potentials(
    double composition,
    double molar_volume
);

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
        coordinates,
        state.physical_amounts
    );
    state.volume = std::exp(log_volume);
    if (!std::isfinite(state.volume) || state.volume <= 0.0) {
        throw std::invalid_argument("manufactured HELD2 phase volume must be positive");
    }
    const double composition = independent.front();
    const auto phase_gradient = manufactured_gibbs_gradient(composition, state.volume);
    state.objective = manufactured_gibbs(composition, state.volume);
    state.gradient = {phase_gradient[0], state.volume * phase_gradient[1]};
    state.hessian = {
        manufactured_composition_hessian(composition),
        0.0,
        0.0,
        5.0 * state.volume * state.volume + state.volume * phase_gradient[1],
    };
    const auto potentials = manufactured_modified_potentials(composition, state.volume);
    state.modified_potentials.assign(potentials.begin(), potentials.end());
    state.pressure_stationarity_relative = phase_gradient[1];
    return state;
}

void configure_held2_ipopt(const Ipopt::SmartPtr<Ipopt::IpoptApplication>& application) {
    application->Options()->SetStringValue("option_file_name", "");
    application->Options()->SetIntegerValue("print_level", 0);
    application->Options()->SetStringValue("sb", "yes");
    application->Options()->SetIntegerValue("max_iter", kStageIIpoptIterations);
    application->Options()->SetNumericValue("tol", 1.0e-10);
    application->Options()->SetNumericValue("acceptable_tol", 1.0e-9);
    application->Options()->SetIntegerValue("acceptable_iter", 0);
    application->Options()->SetStringValue("jacobian_approximation", "exact");
    application->Options()->SetStringValue("hessian_approximation", "exact");
    application->Options()->SetStringValue("nlp_scaling_method", "none");
    application->Options()->SetNumericValue("bound_relax_factor", 0.0);
    application->Options()->SetStringValue("honor_original_bounds", "yes");
    application->Options()->SetStringValue("check_derivatives_for_naninf", "yes");
}

struct Held2SearchRun {
    bool solver_converged = false;
    std::string callback_error;
    std::vector<double> variables;
};

class Held2SearchTnlp final : public Ipopt::TNLP {
public:
    Held2SearchTnlp(
        Held2StateEvaluator evaluator,
        std::vector<double> reference_variables,
        Held2StateEvaluation reference,
        bool use_tpd,
        std::vector<double> initial,
        std::vector<double> lower,
        std::vector<double> upper
    )
        : evaluator_(std::move(evaluator)),
          reference_variables_(std::move(reference_variables)),
          reference_(std::move(reference)),
          use_tpd_(use_tpd),
          initial_(std::move(initial)),
          lower_(std::move(lower)),
          upper_(std::move(upper)) {}

    bool get_nlp_info(
        Ipopt::Index& n,
        Ipopt::Index& m,
        Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag,
        IndexStyleEnum& index_style
    ) override {
        n = static_cast<Ipopt::Index>(initial_.size());
        m = 0;
        nnz_jac_g = 0;
        nnz_h_lag = n * (n + 1) / 2;
        index_style = TNLP::C_STYLE;
        return true;
    }

    bool get_bounds_info(
        Ipopt::Index n,
        Ipopt::Number* x_l,
        Ipopt::Number* x_u,
        Ipopt::Index m,
        Ipopt::Number*,
        Ipopt::Number*
    ) override {
        if (n != static_cast<Ipopt::Index>(initial_.size()) || m != 0) {
            return false;
        }
        std::copy(lower_.begin(), lower_.end(), x_l);
        std::copy(upper_.begin(), upper_.end(), x_u);
        return true;
    }

    bool get_starting_point(
        Ipopt::Index n,
        bool init_x,
        Ipopt::Number* x,
        bool init_z,
        Ipopt::Number*,
        Ipopt::Number*,
        Ipopt::Index m,
        bool init_lambda,
        Ipopt::Number*
    ) override {
        if (n != static_cast<Ipopt::Index>(initial_.size()) || m != 0 || !init_x
            || init_z || init_lambda) {
            return false;
        }
        std::copy(initial_.begin(), initial_.end(), x);
        return true;
    }

    bool eval_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number& objective
    ) override {
        try {
            objective = evaluate(n, x).objective;
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
    }

    bool eval_grad_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number* gradient
    ) override {
        try {
            const Held2StateEvaluation evaluation = evaluate(n, x);
            std::copy(evaluation.gradient.begin(), evaluation.gradient.end(), gradient);
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
    }

    bool eval_g(
        Ipopt::Index,
        const Ipopt::Number*,
        bool,
        Ipopt::Index m,
        Ipopt::Number*
    ) override {
        return m == 0;
    }

    bool eval_jac_g(
        Ipopt::Index,
        const Ipopt::Number*,
        bool,
        Ipopt::Index m,
        Ipopt::Index nonzero_count,
        Ipopt::Index*,
        Ipopt::Index*,
        Ipopt::Number*
    ) override {
        return m == 0 && nonzero_count == 0;
    }

    bool eval_h(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number objective_factor,
        Ipopt::Index m,
        const Ipopt::Number*,
        bool,
        Ipopt::Index nonzero_count,
        Ipopt::Index* rows,
        Ipopt::Index* columns,
        Ipopt::Number* values
    ) override {
        if (m != 0 || nonzero_count != n * (n + 1) / 2) {
            return false;
        }
        Ipopt::Index position = 0;
        if (values == nullptr) {
            for (Ipopt::Index row = 0; row < n; ++row) {
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    rows[position] = row;
                    columns[position] = column;
                    ++position;
                }
            }
            return true;
        }
        try {
            const Held2StateEvaluation evaluation = evaluate(n, x);
            for (Ipopt::Index row = 0; row < n; ++row) {
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    values[position] = objective_factor
                        * evaluation.hessian[static_cast<std::size_t>(row * n + column)];
                    ++position;
                }
            }
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            return false;
        }
    }

    void finalize_solution(
        Ipopt::SolverReturn status,
        Ipopt::Index n,
        const Ipopt::Number* x,
        const Ipopt::Number*,
        const Ipopt::Number*,
        Ipopt::Index,
        const Ipopt::Number*,
        const Ipopt::Number*,
        Ipopt::Number,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*
    ) override {
        variables_.assign(x, x + n);
        solver_converged_ = status == Ipopt::SUCCESS
            || status == Ipopt::STOP_AT_ACCEPTABLE_POINT;
    }

    [[nodiscard]] Held2SearchRun result() const {
        return {solver_converged_, callback_error_, variables_};
    }

private:
    [[nodiscard]] Held2StateEvaluation evaluate(
        Ipopt::Index n,
        const Ipopt::Number* x
    ) const {
        if (n != static_cast<Ipopt::Index>(initial_.size())) {
            throw std::invalid_argument("HELD2 search coordinate count changed");
        }
        std::vector<double> independent(x, x + n - 1);
        Held2StateEvaluation evaluation = evaluator_(independent, x[n - 1]);
        if (!use_tpd_) {
            return evaluation;
        }
        evaluation.objective -= reference_.objective;
        for (std::size_t index = 0; index < evaluation.gradient.size(); ++index) {
            evaluation.objective -= reference_.gradient[index]
                * (x[index] - reference_variables_[index]);
            evaluation.gradient[index] -= reference_.gradient[index];
        }
        return evaluation;
    }

    Held2StateEvaluator evaluator_;
    std::vector<double> reference_variables_;
    Held2StateEvaluation reference_;
    bool use_tpd_ = false;
    std::vector<double> initial_;
    std::vector<double> lower_;
    std::vector<double> upper_;
    bool solver_converged_ = false;
    std::string callback_error_;
    std::vector<double> variables_;
};

Held2SearchRun solve_held2_search(
    const Held2StateEvaluator& evaluator,
    const std::vector<double>& reference_variables,
    const Held2StateEvaluation& reference,
    bool use_tpd,
    const std::vector<double>& initial,
    const std::vector<double>& lower,
    const std::vector<double>& upper
) {
    Ipopt::SmartPtr<Held2SearchTnlp> problem = new Held2SearchTnlp(
        evaluator,
        reference_variables,
        reference,
        use_tpd,
        initial,
        lower,
        upper
    );
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application = IpoptApplicationFactory();
    configure_held2_ipopt(application);
    if (application->Initialize() != Ipopt::Solve_Succeeded) {
        return {};
    }
    application->OptimizeTNLP(problem);
    return problem->result();
}

std::array<double, 2> manufactured_modified_potentials(
    double composition,
    double molar_volume
) {
    const double gibbs = manufactured_gibbs(composition, molar_volume);
    const auto gradient = manufactured_gibbs_gradient(composition, molar_volume);
    const double common_volume_term = -molar_volume * gradient[1];
    return {
        gibbs - composition * gradient[0] + common_volume_term,
        gibbs + (1.0 - composition) * gradient[0] + common_volume_term,
    };
}

double enumerated_manufactured_objective(double feed_composition) {
    std::vector<double> compositions;
    std::vector<double> volumes;
    compositions.reserve(104);
    volumes.reserve(102);
    for (int index = 0; index <= 100; ++index) {
        compositions.push_back(static_cast<double>(index) / 100.0);
        volumes.push_back(0.5 + static_cast<double>(index) / 100.0);
    }
    compositions.insert(compositions.end(), {0.2, 0.5, 0.8});
    volumes.push_back(1.0);
    std::sort(compositions.begin(), compositions.end());
    compositions.erase(std::unique(compositions.begin(), compositions.end()), compositions.end());
    std::sort(volumes.begin(), volumes.end());
    volumes.erase(std::unique(volumes.begin(), volumes.end()), volumes.end());

    std::vector<double> phase_minima;
    phase_minima.reserve(compositions.size());
    for (double composition : compositions) {
        double minimum = std::numeric_limits<double>::infinity();
        for (double volume : volumes) {
            minimum = std::min(minimum, manufactured_gibbs(composition, volume));
        }
        phase_minima.push_back(minimum);
    }

    double best = std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < compositions.size(); ++index) {
        if (std::abs(compositions[index] - feed_composition) <= kCoordinateTolerance) {
            best = std::min(best, phase_minima[index]);
        }
        if (compositions[index] > feed_composition) {
            continue;
        }
        for (std::size_t right = index + 1; right < compositions.size(); ++right) {
            if (compositions[right] < feed_composition) {
                continue;
            }
            const double separation = compositions[right] - compositions[index];
            if (separation <= 0.0) {
                continue;
            }
            const double fraction =
                (compositions[right] - feed_composition) / separation;
            best = std::min(
                best,
                fraction * phase_minima[index]
                    + (1.0 - fraction) * phase_minima[right]
            );
        }
    }
    return best;
}

}  // namespace

Held2StateEvaluation evaluate_held2_manufactured_state(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent_modified_fractions,
    double log_volume
) {
    return evaluate_manufactured_state_impl(
        coordinates,
        independent_modified_fractions,
        log_volume
    );
}

double held2_manufactured_enumerated_objective(double feed_composition) {
    return enumerated_manufactured_objective(feed_composition);
}

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
        result.modified_factors.push_back(1.0 - charges[index] / charges[eliminated]);
        if (charges[index] == 0.0) {
            result.dependent_index = index;
        }
    }
    for (std::size_t retained = 0; retained < result.retained_indices.size(); ++retained) {
        const std::size_t index = result.retained_indices[retained];
        if (index == result.dependent_index) {
            continue;
        }
        result.independent_indices.push_back(index);
        const double factor = result.modified_factors[retained];
        result.independent_lower_bounds.push_back(kModifiedLowerScale * factor);
        if (charges[index] == 0.0) {
            result.independent_upper_bounds.push_back(1.0);
            continue;
        }
        double largest_opposite = 0.0;
        for (double other_charge : charges) {
            if (other_charge * charges[index] < 0.0) {
                largest_opposite = std::max(largest_opposite, std::abs(other_charge));
            }
        }
        if (largest_opposite == 0.0) {
            throw std::invalid_argument("charged species has no opposite-sign counterion");
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
        throw std::invalid_argument("physical composition size does not match charges");
    }
    require_finite_vector(physical_fractions, "physical mole fractions");
    if (!std::all_of(physical_fractions.begin(), physical_fractions.end(), [](double value) {
            return value >= 0.0;
        })) {
        throw std::invalid_argument("physical mole fractions must be nonnegative");
    }
    double total = 0.0;
    for (double value : physical_fractions) {
        total += value;
    }
    if (std::abs(total - 1.0) > kCoordinateTolerance) {
        throw std::invalid_argument("physical mole fractions must sum to one");
    }
    if (std::abs(charge_residual(coordinates.charges, physical_fractions))
        > kCoordinateTolerance) {
        throw std::invalid_argument("physical feed must be electroneutral");
    }
    std::vector<double> modified;
    modified.reserve(coordinates.retained_indices.size());
    for (std::size_t retained = 0; retained < coordinates.retained_indices.size(); ++retained) {
        modified.push_back(
            coordinates.modified_factors[retained]
            * physical_fractions[coordinates.retained_indices[retained]]
        );
    }
    double modified_total = 0.0;
    for (double value : modified) {
        modified_total += value;
    }
    if (std::abs(modified_total - 1.0) > kCoordinateTolerance) {
        throw std::invalid_argument("modified mole fractions do not sum to one");
    }
    return modified;
}

std::vector<double> held2_lift_modified_fractions(
    const Held2Coordinates& coordinates,
    const std::vector<double>& modified_fractions
) {
    if (modified_fractions.size() != coordinates.retained_indices.size()) {
        throw std::invalid_argument("modified composition size does not match retained species");
    }
    require_finite_vector(modified_fractions, "modified mole fractions");
    double modified_total = 0.0;
    for (double value : modified_fractions) {
        if (value < 0.0) {
            throw std::invalid_argument("modified mole fractions must be nonnegative");
        }
        modified_total += value;
    }
    if (std::abs(modified_total - 1.0) > kCoordinateTolerance) {
        throw std::invalid_argument("modified mole fractions must sum to one");
    }

    std::vector<double> physical(coordinates.charges.size(), 0.0);
    for (std::size_t retained = 0; retained < coordinates.retained_indices.size(); ++retained) {
        physical[coordinates.retained_indices[retained]] =
            modified_fractions[retained] / coordinates.modified_factors[retained];
    }
    double eliminated = 0.0;
    const double eliminated_charge = coordinates.charges[coordinates.eliminated_index];
    for (std::size_t index : coordinates.retained_indices) {
        eliminated -= coordinates.charges[index] / eliminated_charge * physical[index];
    }
    if (eliminated < -kCoordinateTolerance) {
        throw std::invalid_argument("eliminated ion amount must be nonnegative");
    }
    physical[coordinates.eliminated_index] = std::max(0.0, eliminated);
    return physical;
}

std::vector<double> held2_lift_independent_fractions(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent_modified_fractions
) {
    const std::size_t independent_count = coordinates.independent_indices.size();
    if (independent_modified_fractions.size() != independent_count) {
        throw std::invalid_argument(
            "independent modified composition size does not match the HELD2 chart"
        );
    }
    require_finite_vector(independent_modified_fractions, "independent modified fractions");
    std::vector<double> modified_fractions(coordinates.retained_indices.size(), 0.0);
    double independent_sum = 0.0;
    for (std::size_t independent = 0; independent < independent_count; ++independent) {
        const double value = independent_modified_fractions[independent];
        if (value < coordinates.independent_lower_bounds[independent]
            || value > coordinates.independent_upper_bounds[independent]) {
            throw std::invalid_argument("independent modified fraction is outside its source bound");
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
        throw std::invalid_argument("independent modified fractions must leave a positive dependent fraction");
    }
    const auto dependent_retained = static_cast<std::size_t>(
        std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            coordinates.dependent_index
        ) - coordinates.retained_indices.begin()
    );
    modified_fractions[dependent_retained] = 1.0 - independent_sum;
    return held2_lift_modified_fractions(coordinates, modified_fractions);
}

std::vector<double> held2_transform_modified_potentials(
    const Held2Coordinates& coordinates,
    const std::vector<double>& chemical_potentials
) {
    if (chemical_potentials.size() != coordinates.charges.size()) {
        throw std::invalid_argument("chemical-potential size does not match charges");
    }
    require_finite_vector(chemical_potentials, "chemical potentials");
    const double eliminated_potential =
        chemical_potentials[coordinates.eliminated_index];
    const double eliminated_charge = coordinates.charges[coordinates.eliminated_index];
    std::vector<double> modified;
    modified.reserve(coordinates.retained_indices.size());
    for (std::size_t retained = 0; retained < coordinates.retained_indices.size(); ++retained) {
        const std::size_t index = coordinates.retained_indices[retained];
        const double constrained = chemical_potentials[index]
            - coordinates.charges[index] / eliminated_charge * eliminated_potential;
        modified.push_back(constrained / coordinates.modified_factors[retained]);
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
    const std::size_t independent_count = coordinates.independent_indices.size();
    const std::size_t provider_coordinate_count = component_count + 1;
    const std::size_t reduced_coordinate_count = independent_count + 1;
    if (independent_modified_fractions.size() != independent_count) {
        throw std::invalid_argument(
            "independent modified composition size does not match the HELD2 chart"
        );
    }
    require_finite_vector(independent_modified_fractions, "independent modified fractions");
    if (!std::isfinite(log_volume) || !std::isfinite(pressure_over_rt)
        || !std::isfinite(target_pressure_pa) || pressure_over_rt <= 0.0
        || target_pressure_pa <= 0.0) {
        throw std::invalid_argument("HELD2 phase state and pressure scales must be finite and positive");
    }
    if (block.gradient.size() != provider_coordinate_count
        || block.hessian.size() != provider_coordinate_count * provider_coordinate_count) {
        throw std::invalid_argument("HELD2 physical phase block has the wrong tensor dimensions");
    }
    require_finite_vector(block.gradient, "HELD2 physical phase gradient");
    require_finite_vector(block.hessian, "HELD2 physical phase Hessian");
    if (!std::isfinite(block.helmholtz_over_rt) || !std::isfinite(block.pressure_pa)) {
        throw std::invalid_argument("HELD2 physical phase block scalars must be finite");
    }

    Held2StateEvaluation result;
    result.modified_fractions.resize(coordinates.retained_indices.size(), 0.0);
    for (std::size_t independent = 0; independent < independent_count; ++independent) {
        const double value = independent_modified_fractions[independent];
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
        coordinates,
        independent_modified_fractions
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
        throw std::invalid_argument("HELD2 phase volume must be finite and positive");
    }

    std::vector<double> jacobian(
        provider_coordinate_count * reduced_coordinate_count,
        0.0
    );
    for (std::size_t independent = 0; independent < independent_count; ++independent) {
        for (std::size_t retained = 0; retained < coordinates.retained_indices.size(); ++retained) {
            const std::size_t component = coordinates.retained_indices[retained];
            double modified_derivative = 0.0;
            if (component == coordinates.independent_indices[independent]) {
                modified_derivative = 1.0;
            } else if (component == coordinates.dependent_index) {
                modified_derivative = -1.0;
            }
            jacobian[component * reduced_coordinate_count + independent] =
                modified_derivative / coordinates.modified_factors[retained];
        }
        double eliminated_derivative = 0.0;
        const double eliminated_charge = coordinates.charges[coordinates.eliminated_index];
        for (std::size_t component : coordinates.retained_indices) {
            eliminated_derivative -= coordinates.charges[component] / eliminated_charge
                * jacobian[component * reduced_coordinate_count + independent];
        }
        jacobian[coordinates.eliminated_index * reduced_coordinate_count + independent] =
            eliminated_derivative;
    }
    jacobian[component_count * reduced_coordinate_count + independent_count] = result.volume;

    std::vector<double> augmented_gradient = block.gradient;
    augmented_gradient.back() += pressure_over_rt;
    result.objective = block.helmholtz_over_rt + pressure_over_rt * result.volume;
    result.gradient.assign(reduced_coordinate_count, 0.0);
    for (std::size_t reduced = 0; reduced < reduced_coordinate_count; ++reduced) {
        for (std::size_t provider = 0; provider < provider_coordinate_count; ++provider) {
            result.gradient[reduced] +=
                jacobian[provider * reduced_coordinate_count + reduced]
                * augmented_gradient[provider];
        }
    }
    result.hessian.assign(reduced_coordinate_count * reduced_coordinate_count, 0.0);
    for (std::size_t row = 0; row < reduced_coordinate_count; ++row) {
        for (std::size_t column = 0; column < reduced_coordinate_count; ++column) {
            double value = 0.0;
            for (std::size_t left = 0; left < provider_coordinate_count; ++left) {
                for (std::size_t right = 0; right < provider_coordinate_count; ++right) {
                    value += jacobian[left * reduced_coordinate_count + row]
                        * block.hessian[left * provider_coordinate_count + right]
                        * jacobian[right * reduced_coordinate_count + column];
                }
            }
            result.hessian[row * reduced_coordinate_count + column] = value;
        }
    }
    result.hessian.back() += result.volume * augmented_gradient.back();
    std::vector<double> physical_potentials(block.gradient.begin(), block.gradient.end() - 1);
    result.modified_potentials = held2_transform_modified_potentials(
        coordinates,
        physical_potentials
    );
    result.pressure_stationarity_relative =
        (block.pressure_pa - target_pressure_pa) / target_pressure_pa;
    return result;
}

Held2StageIResult solve_held2_manufactured_stage_i(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed
) {
    const Held2Coordinates coordinates = make_held2_coordinates(charges);
    if (coordinates.independent_indices.size() != 1) {
        throw std::invalid_argument(
            "manufactured HELD2 Stage I requires one independent modified composition"
        );
    }
    const std::vector<double> modified_feed = held2_transform_physical_fractions(
        coordinates,
        physical_feed
    );
    const auto independent_retained = static_cast<std::size_t>(
        std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            coordinates.independent_indices.front()
        ) - coordinates.retained_indices.begin()
    );
    const double feed_independent = modified_feed[independent_retained];
    const Held2StateEvaluator evaluator = [coordinates](
        const std::vector<double>& independent,
        double log_volume
    ) {
        return evaluate_manufactured_state_impl(
            coordinates,
            independent,
            log_volume
        );
    };

    const std::vector<double> reference_initial = {feed_independent, 0.0};
    const Held2StateEvaluation initial_reference = evaluator({feed_independent}, 0.0);
    const Held2SearchRun reference_run = solve_held2_search(
        evaluator,
        reference_initial,
        initial_reference,
        false,
        reference_initial,
        {feed_independent, std::log(0.5)},
        {feed_independent, std::log(1.5)}
    );
    Held2StageIResult result;
    result.declared_start_count = 10 * static_cast<int>(charges.size());
    result.minimum_tpd = std::numeric_limits<double>::infinity();
    if (!reference_run.solver_converged || !reference_run.callback_error.empty()
        || reference_run.variables.size() != 2) {
        result.outcome = "indeterminate";
        return result;
    }
    const Held2StateEvaluation reference = evaluator(
        {reference_run.variables.front()},
        reference_run.variables.back()
    );
    result.reference_modified_fractions = reference.modified_fractions;
    result.reference_volume = reference.volume;

    const double lower_composition = coordinates.independent_lower_bounds.front();
    const double upper_composition = std::nextafter(
        coordinates.independent_upper_bounds.front(),
        lower_composition
    );
    const std::vector<double> lower = {lower_composition, std::log(0.5)};
    const std::vector<double> upper = {upper_composition, std::log(1.5)};
    std::mt19937 generator(kStageISeed);
    std::uniform_real_distribution<double> composition_distribution(
        lower_composition,
        upper_composition
    );
    std::uniform_real_distribution<double> log_volume_distribution(
        lower.back(),
        upper.back()
    );
    std::vector<std::vector<double>> starts;
    starts.reserve(static_cast<std::size_t>(result.declared_start_count));
    starts.push_back({0.2, 0.0});
    starts.push_back({0.8, 0.0});
    while (starts.size() < static_cast<std::size_t>(result.declared_start_count)) {
        starts.push_back({
            composition_distribution(generator),
            log_volume_distribution(generator),
        });
    }
    const std::vector<double> reference_variables = {
        reference_run.variables.front(),
        reference_run.variables.back(),
    };
    for (const std::vector<double>& start : starts) {
        const Held2SearchRun run = solve_held2_search(
            evaluator,
            reference_variables,
            reference,
            true,
            start,
            lower,
            upper
        );
        if (!run.solver_converged || !run.callback_error.empty() || run.variables.size() != 2) {
            ++result.failed_start_count;
            continue;
        }
        ++result.completed_start_count;
        Held2StateEvaluation state = evaluator({run.variables.front()}, run.variables.back());
        double tpd = state.objective - reference.objective;
        for (std::size_t index = 0; index < state.gradient.size(); ++index) {
            tpd -= reference.gradient[index]
                * (run.variables[index] - reference_variables[index]);
        }
        result.minimum_tpd = std::min(result.minimum_tpd, tpd);
        if (tpd >= kStageITpdThreshold) {
            continue;
        }
        const bool duplicate = std::any_of(
            result.candidates.begin(),
            result.candidates.end(),
            [&state](const Held2StageICandidate& candidate) {
                return maximum_abs_difference(
                           candidate.modified_fractions,
                           state.modified_fractions
                       ) < 1.0e-7
                    && std::abs(candidate.volume - state.volume) < 1.0e-7;
            }
        );
        if (!duplicate) {
            result.candidates.push_back({state.modified_fractions, state.volume, tpd});
        }
    }
    std::sort(
        result.candidates.begin(),
        result.candidates.end(),
        [](const Held2StageICandidate& left, const Held2StageICandidate& right) {
            return left.modified_fractions < right.modified_fractions;
        }
    );
    if (result.failed_start_count != 0) {
        result.outcome = "indeterminate";
    } else if (!result.candidates.empty()) {
        result.outcome = "negative_tpd";
    } else {
        result.outcome = "no_negative_found";
    }
    return result;
}

Held2StageIIResult solve_held2_manufactured_stage_ii(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed
) {
    constexpr int major_iteration_cap = 100;
    constexpr double stage_ii_tolerance = 1.0e-8;
    const Held2Coordinates coordinates = make_held2_coordinates(charges);
    if (coordinates.independent_indices.size() != 1) {
        throw std::invalid_argument(
            "manufactured HELD2 Stage II requires one independent modified composition"
        );
    }
    const std::vector<double> modified_feed = held2_transform_physical_fractions(
        coordinates,
        physical_feed
    );
    const auto independent_retained = static_cast<std::size_t>(
        std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            coordinates.independent_indices.front()
        ) - coordinates.retained_indices.begin()
    );
    const double feed = modified_feed[independent_retained];
    const Held2StateEvaluator evaluator = [coordinates](
        const std::vector<double>& independent,
        double log_volume
    ) {
        return evaluate_manufactured_state_impl(
            coordinates,
            independent,
            log_volume
        );
    };
    std::vector<Held2StateEvaluation> cuts = {
        evaluator({feed}, 0.0),
        evaluator({0.5 * feed}, 0.0),
        evaluator({0.5 * (1.0 + feed)}, 0.0),
    };
    const double feed_gibbs = cuts.front().objective;
    const double lower_composition = coordinates.independent_lower_bounds.front();
    const double upper_composition = std::nextafter(
        coordinates.independent_upper_bounds.front(),
        lower_composition
    );
    const std::vector<double> lower = {lower_composition, std::log(0.5)};
    const std::vector<double> upper = {upper_composition, std::log(1.5)};
    Held2StageIIResult result;
    result.lower_starts_per_iteration = 10 * static_cast<int>(charges.size());

    for (int major = 0; major < major_iteration_cap; ++major) {
        struct Line {
            double intercept = 0.0;
            double slope = 0.0;
        };
        std::vector<Line> lines;
        lines.reserve(cuts.size() + 1);
        for (const Held2StateEvaluation& cut : cuts) {
            lines.push_back({
                cut.objective,
                feed - cut.modified_fractions[independent_retained],
            });
        }
        lines.push_back({feed_gibbs, 0.0});
        std::vector<double> multiplier_candidates = {0.0};
        for (std::size_t left = 0; left < lines.size(); ++left) {
            for (std::size_t right = left + 1; right < lines.size(); ++right) {
                const double slope_delta = lines[left].slope - lines[right].slope;
                if (std::abs(slope_delta) <= std::numeric_limits<double>::epsilon()) {
                    continue;
                }
                multiplier_candidates.push_back(
                    (lines[right].intercept - lines[left].intercept) / slope_delta
                );
            }
        }
        double upper_bound = -std::numeric_limits<double>::infinity();
        double multiplier = 0.0;
        for (double candidate : multiplier_candidates) {
            double envelope = std::numeric_limits<double>::infinity();
            for (const Line& line : lines) {
                envelope = std::min(envelope, line.intercept + line.slope * candidate);
            }
            if (envelope > upper_bound) {
                upper_bound = envelope;
                multiplier = candidate;
            }
        }

        Held2StateEvaluation lower_reference;
        lower_reference.gradient = {multiplier, 0.0};
        lower_reference.hessian.assign(4, 0.0);
        const std::vector<double> lower_reference_variables = {feed, 0.0};
        std::mt19937 generator(kStageISeed + static_cast<unsigned int>(major));
        std::uniform_real_distribution<double> composition_distribution(
            lower_composition,
            upper_composition
        );
        std::uniform_real_distribution<double> log_volume_distribution(
            lower.back(),
            upper.back()
        );
        std::vector<std::vector<double>> starts = {{0.2, 0.0}, {0.8, 0.0}};
        while (starts.size() < static_cast<std::size_t>(result.lower_starts_per_iteration)) {
            starts.push_back({
                composition_distribution(generator),
                log_volume_distribution(generator),
            });
        }
        double lower_bound = std::numeric_limits<double>::infinity();
        Held2StateEvaluation best_state;
        bool lower_failed = false;
        for (const std::vector<double>& start : starts) {
            const Held2SearchRun run = solve_held2_search(
                evaluator,
                lower_reference_variables,
                lower_reference,
                true,
                start,
                lower,
                upper
            );
            if (!run.solver_converged || !run.callback_error.empty()
                || run.variables.size() != 2) {
                lower_failed = true;
                continue;
            }
            Held2StateEvaluation state = evaluator(
                {run.variables.front()},
                run.variables.back()
            );
            const double value = state.objective
                + multiplier * (feed - state.modified_fractions[independent_retained]);
            if (value < lower_bound) {
                lower_bound = value;
                best_state = std::move(state);
            }
        }
        ++result.major_iterations;
        if (lower_failed || !std::isfinite(lower_bound)) {
            result.outcome = "indeterminate";
            result.cut_count = static_cast<int>(cuts.size());
            return result;
        }
        result.bound_history.push_back({
            lower_bound,
            upper_bound,
            multiplier,
            static_cast<int>(cuts.size()),
        });
        const bool duplicate = std::any_of(
            cuts.begin(),
            cuts.end(),
            [&best_state](const Held2StateEvaluation& cut) {
                return maximum_abs_difference(
                           cut.modified_fractions,
                           best_state.modified_fractions
                       ) < 1.0e-7
                    && std::abs(cut.volume - best_state.volume) < 1.0e-7;
            }
        );
        if (!duplicate) {
            cuts.push_back(best_state);
        }
        result.candidates.clear();
        for (const Held2StateEvaluation& cut : cuts) {
            const double lower_value = cut.objective
                + multiplier * (feed - cut.modified_fractions[independent_retained]);
            const double gap = upper_bound - lower_value;
            if (gap >= -stage_ii_tolerance && gap <= stage_ii_tolerance
                && std::abs(cut.gradient.front() - multiplier) <= stage_ii_tolerance
                && std::abs(cut.gradient.back()) <= stage_ii_tolerance) {
                result.candidates.push_back({
                    cut.modified_fractions,
                    cut.volume,
                    gap,
                });
            }
        }
        if (result.candidates.size() > 1) {
            result.outcome = "candidate_set";
            result.cut_count = static_cast<int>(cuts.size());
            return result;
        }
        if (duplicate) {
            result.outcome = "no_progress";
            result.cut_count = static_cast<int>(cuts.size());
            return result;
        }
    }
    result.outcome = "resource_limit";
    result.cut_count = static_cast<int>(cuts.size());
    return result;
}

Held2ManufacturedEvaluation evaluate_held2_manufactured(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::array<double, 4>& variables,
    const std::vector<double>& chemical_potentials
) {
    if (!std::all_of(variables.begin(), variables.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::invalid_argument("HELD2 reduced vector must contain four finite entries");
    }
    Held2ManufacturedEvaluation result;
    result.coordinates = make_held2_coordinates(charges);
    if (result.coordinates.retained_indices.size() != 2
        || result.coordinates.independent_indices.size() != 1) {
        throw std::invalid_argument(
            "manufactured HELD2 seam requires one independent modified composition"
        );
    }
    result.modified_feed = held2_transform_physical_fractions(
        result.coordinates,
        physical_feed
    );
    result.transformed_modified_potentials = held2_transform_modified_potentials(
        result.coordinates,
        chemical_potentials
    );
    const auto dependent_position = static_cast<std::size_t>(
        std::find(
            result.coordinates.retained_indices.begin(),
            result.coordinates.retained_indices.end(),
            result.coordinates.dependent_index
        ) - result.coordinates.retained_indices.begin()
    );
    const std::size_t independent_position = 1 - dependent_position;
    const double feed = result.modified_feed[independent_position];
    const double u_alpha = variables[0];
    const double u_beta = variables[1];
    const double volume_alpha = variables[2];
    const double volume_beta = variables[3];
    if (!(u_alpha < feed && feed < u_beta)) {
        throw std::invalid_argument("modified compositions must straddle the modified feed");
    }
    const double lower_bound = result.coordinates.independent_lower_bounds[0];
    const double upper_bound = result.coordinates.independent_upper_bounds[0];
    if (u_alpha < lower_bound || u_beta > upper_bound) {
        throw std::invalid_argument(
            "HELD2 modified compositions must satisfy the declared source bounds"
        );
    }
    if (volume_alpha <= 0.0 || volume_beta <= 0.0) {
        throw std::invalid_argument("HELD2 manufactured state is outside its physical domain");
    }
    result.modified_phases[0] = std::vector<double>(2, 0.0);
    result.modified_phases[1] = std::vector<double>(2, 0.0);
    result.modified_phases[0][dependent_position] = 1.0 - u_alpha;
    result.modified_phases[0][independent_position] = u_alpha;
    result.modified_phases[1][dependent_position] = 1.0 - u_beta;
    result.modified_phases[1][independent_position] = u_beta;
    result.physical_phases[0] = held2_lift_modified_fractions(
        result.coordinates,
        result.modified_phases[0]
    );
    result.physical_phases[1] = held2_lift_modified_fractions(
        result.coordinates,
        result.modified_phases[1]
    );

    const double separation = u_beta - u_alpha;
    result.phase_fraction = (u_beta - feed) / separation;
    const double complement = 1.0 - result.phase_fraction;
    const double gibbs_alpha = manufactured_gibbs(u_alpha, volume_alpha);
    const double gibbs_beta = manufactured_gibbs(u_beta, volume_beta);
    result.phase_gibbs_gradients[0] = manufactured_gibbs_gradient(
        u_alpha,
        volume_alpha
    );
    result.phase_gibbs_gradients[1] = manufactured_gibbs_gradient(
        u_beta,
        volume_beta
    );
    result.objective =
        result.phase_fraction * gibbs_alpha + complement * gibbs_beta;
    result.gradient = {
        result.phase_fraction / separation * (gibbs_alpha - gibbs_beta)
            + result.phase_fraction * result.phase_gibbs_gradients[0][0],
        complement / separation * (gibbs_alpha - gibbs_beta)
            + complement * result.phase_gibbs_gradients[1][0],
        result.phase_fraction * result.phase_gibbs_gradients[0][1],
        complement * result.phase_gibbs_gradients[1][1],
    };
    result.phase_modified_potentials[0] = manufactured_modified_potentials(
        u_alpha,
        volume_alpha
    );
    result.phase_modified_potentials[1] = manufactured_modified_potentials(
        u_beta,
        volume_beta
    );
    result.modified_potential_gap = std::max(
        std::abs(
            result.phase_modified_potentials[0][0]
            - result.phase_modified_potentials[1][0]
        ),
        std::abs(
            result.phase_modified_potentials[0][1]
            - result.phase_modified_potentials[1][1]
        )
    );
    result.pressure_stationarity_inf_norm = std::max(
        std::abs(result.phase_gibbs_gradients[0][1]),
        std::abs(result.phase_gibbs_gradients[1][1])
    );

    result.modified_balance.resize(2, 0.0);
    result.ordinary_balance.resize(charges.size(), 0.0);
    for (std::size_t index = 0; index < result.modified_balance.size(); ++index) {
        result.modified_balance[index] =
            result.phase_fraction * result.modified_phases[0][index]
            + complement * result.modified_phases[1][index];
    }
    for (std::size_t index = 0; index < result.ordinary_balance.size(); ++index) {
        result.ordinary_balance[index] =
            result.phase_fraction * result.physical_phases[0][index]
            + complement * result.physical_phases[1][index];
    }
    result.phase_charge_residuals = {
        charge_residual(charges, result.physical_phases[0]),
        charge_residual(charges, result.physical_phases[1]),
    };

    result.certificate.modified_balance_abs = maximum_abs_difference(
        result.modified_balance,
        result.modified_feed
    );
    result.certificate.ordinary_balance_inf_norm = maximum_abs_difference(
        result.ordinary_balance,
        physical_feed
    );
    result.certificate.phase_charge_inf_norm = std::max(
        std::abs(result.phase_charge_residuals[0]),
        std::abs(result.phase_charge_residuals[1])
    );
    result.certificate.modified_potential_gap = result.modified_potential_gap;
    result.certificate.pressure_stationarity_inf_norm =
        result.pressure_stationarity_inf_norm;
    for (double value : result.gradient) {
        result.certificate.reduced_kkt_inf_norm = std::max(
            result.certificate.reduced_kkt_inf_norm,
            std::abs(value)
        );
    }
    result.certificate.enumeration_objective_gap =
        result.objective - enumerated_manufactured_objective(feed);
    result.certificate.independent_modified_composition_count = 1.0;
    const double maximum_metric = std::max({
        result.certificate.modified_balance_abs,
        result.certificate.ordinary_balance_inf_norm,
        result.certificate.phase_charge_inf_norm,
        result.certificate.modified_potential_gap,
        result.certificate.pressure_stationarity_inf_norm,
        result.certificate.reduced_kkt_inf_norm,
        std::abs(result.certificate.enumeration_objective_gap),
    });
    result.certificate.accepted = maximum_metric < kCertificateTolerance;
    result.certificate.independent_evidence =
        std::abs(result.certificate.enumeration_objective_gap) < kCertificateTolerance;
    return result;
}

}  // namespace epcsaft_equilibrium
