#include "chemical_equilibrium.hpp"
#include "provider.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Highs.h>
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>
#include <coin/IpTNLP.hpp>

namespace epcsaft_equilibrium {
namespace {

constexpr double kGasConstantJPerMolK = 8.31446261815324;
constexpr double kInfinity = 1.0e19;
constexpr double kBalanceTolerance = 1.0e-9;
constexpr double kPressureTolerance = 1.0e-8;
constexpr double kAffinityTolerance = 1.0e-7;
constexpr double kKktTolerance = 1.0e-7;
constexpr double kProviderDomainTolerance = 1.0e-12;
constexpr double kSensitivityConditionNumberMax = 1.0e6;

std::string ipopt_status_name(Ipopt::ApplicationReturnStatus status) {
    switch (status) {
        case Ipopt::Solve_Succeeded:
            return "solve_succeeded";
        case Ipopt::Solved_To_Acceptable_Level:
            return "solved_to_acceptable_level";
        case Ipopt::Infeasible_Problem_Detected:
            return "infeasible_problem_detected";
        case Ipopt::Search_Direction_Becomes_Too_Small:
            return "search_direction_too_small";
        case Ipopt::Diverging_Iterates:
            return "diverging_iterates";
        case Ipopt::Feasible_Point_Found:
            return "feasible_point_found";
        case Ipopt::Maximum_Iterations_Exceeded:
            return "maximum_iterations_exceeded";
        case Ipopt::Restoration_Failed:
            return "restoration_failed";
        case Ipopt::Error_In_Step_Computation:
            return "error_in_step_computation";
        case Ipopt::Not_Enough_Degrees_Of_Freedom:
            return "not_enough_degrees_of_freedom";
        case Ipopt::Invalid_Problem_Definition:
            return "invalid_problem_definition";
        case Ipopt::Invalid_Option:
            return "invalid_option";
        case Ipopt::Invalid_Number_Detected:
            return "invalid_number_detected";
        case Ipopt::Unrecoverable_Exception:
            return "unrecoverable_exception";
        case Ipopt::NonIpopt_Exception_Thrown:
            return "non_ipopt_exception_thrown";
        case Ipopt::Insufficient_Memory:
            return "insufficient_memory";
        case Ipopt::Internal_Error:
            return "internal_error";
        default:
            return "ipopt_status_" + std::to_string(static_cast<int>(status));
    }
}

void configure_ipopt(
    const Ipopt::SmartPtr<Ipopt::IpoptApplication>& application,
    int max_iterations
) {
    application->Options()->SetStringValue("option_file_name", "");
    application->Options()->SetIntegerValue("print_level", 0);
    application->Options()->SetStringValue("sb", "yes");
    application->Options()->SetIntegerValue("max_iter", max_iterations);
    application->Options()->SetNumericValue("tol", 1.0e-10);
    application->Options()->SetNumericValue("acceptable_tol", 1.0e-9);
    application->Options()->SetIntegerValue("acceptable_iter", 0);
    application->Options()->SetNumericValue("constr_viol_tol", 1.0e-10);
    application->Options()->SetStringValue("jacobian_approximation", "exact");
    application->Options()->SetStringValue("hessian_approximation", "exact");
    application->Options()->SetStringValue("nlp_scaling_method", "none");
    application->Options()->SetNumericValue("bound_relax_factor", 0.0);
    application->Options()->SetStringValue("honor_original_bounds", "yes");
    application->Options()->SetStringValue("check_derivatives_for_naninf", "yes");
}

double vector_inf_norm(const std::vector<double>& values) {
    double norm = 0.0;
    for (double value : values) {
        norm = std::max(norm, std::abs(value));
    }
    return norm;
}

double matrix_row_l2_norm(const DenseMatrix& matrix, std::size_t row) {
    double norm = 0.0;
    for (std::size_t column = 0; column < matrix.columns; ++column) {
        norm = std::hypot(norm, matrix(row, column));
    }
    if (norm == 0.0) {
        throw std::invalid_argument("matrix row must not be zero");
    }
    return norm;
}

std::vector<double> matrix_vector(
    const DenseMatrix& matrix,
    const std::vector<double>& vector
) {
    if (matrix.columns != vector.size()) {
        throw std::invalid_argument("matrix-vector dimensions are inconsistent");
    }
    std::vector<double> result(matrix.rows, 0.0);
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        for (std::size_t column = 0; column < matrix.columns; ++column) {
            result[row] += matrix(row, column) * vector[column];
        }
    }
    return result;
}

double reaction_residual_inf_norm(
    const DenseMatrix& reactions,
    const std::vector<double>& potentials
) {
    double result = 0.0;
    for (std::size_t reaction = 0; reaction < reactions.rows; ++reaction) {
        double value = 0.0;
        for (std::size_t species = 0; species < reactions.columns; ++species) {
            value += reactions(reaction, species) * potentials[species];
        }
        result = std::max(
            result,
            std::abs(value) / matrix_row_l2_norm(reactions, reaction)
        );
    }
    return result;
}

bool add_independent_row(
    const std::vector<double>& row,
    std::vector<std::vector<double>>& orthonormal_basis
) {
    std::vector<double> residual = row;
    const double original_norm = std::sqrt(std::inner_product(
        row.begin(), row.end(), row.begin(), 0.0
    ));
    if (original_norm == 0.0) {
        return false;
    }
    for (const std::vector<double>& basis : orthonormal_basis) {
        const double projection = std::inner_product(
            residual.begin(), residual.end(), basis.begin(), 0.0
        );
        for (std::size_t column = 0; column < residual.size(); ++column) {
            residual[column] -= projection * basis[column];
        }
    }
    const double residual_norm = std::sqrt(std::inner_product(
        residual.begin(), residual.end(), residual.begin(), 0.0
    ));
    if (residual_norm <= 4096.0 * std::numeric_limits<double>::epsilon()
            * std::max(1.0, original_norm) * static_cast<double>(row.size())) {
        return false;
    }
    for (double& value : residual) {
        value /= residual_norm;
    }
    orthonormal_basis.push_back(std::move(residual));
    return true;
}

struct ConstraintRows {
    DenseMatrix matrix;
    std::vector<double> totals;
};

ConstraintRows independent_max_min_rows(
    const DenseMatrix& balance_matrix,
    const std::vector<double>& feed,
    const std::vector<int>& charges
) {
    ConstraintRows result;
    result.matrix.columns = feed.size();
    std::vector<std::vector<double>> basis;
    for (std::size_t row = 0; row < balance_matrix.rows; ++row) {
        std::vector<double> values(feed.size(), 0.0);
        for (std::size_t species = 0; species < feed.size(); ++species) {
            values[species] = balance_matrix(row, species);
        }
        if (add_independent_row(values, basis)) {
            result.matrix.values.insert(
                result.matrix.values.end(), values.begin(), values.end()
            );
            result.totals.push_back(std::inner_product(
                values.begin(), values.end(), feed.begin(), 0.0
            ));
        }
    }
    std::vector<double> charge_row(charges.begin(), charges.end());
    if (add_independent_row(charge_row, basis)) {
        result.matrix.values.insert(
            result.matrix.values.end(), charge_row.begin(), charge_row.end()
        );
        result.totals.push_back(0.0);
    }
    result.matrix.rows = result.totals.size();
    return result;
}

struct LinearInitialization {
    std::string status;
    std::vector<double> values;
};

LinearInitialization solve_initialization_lp(
    const ConstraintRows& constraints,
    const std::vector<int>& charges,
    double ion_fraction_max
) {
    const std::size_t species_count = charges.size();
    const std::size_t minimum_column = species_count;
    const bool ion_limited = std::isfinite(ion_fraction_max);
    const std::size_t row_count =
        constraints.matrix.rows + species_count + (ion_limited ? 1 : 0);
    HighsModel model;
    model.lp_.num_col_ = static_cast<HighsInt>(species_count + 1);
    model.lp_.num_row_ = static_cast<HighsInt>(row_count);
    model.lp_.sense_ = ObjSense::kMaximize;
    model.lp_.col_cost_.assign(species_count + 1, 0.0);
    model.lp_.col_cost_[minimum_column] = 1.0;
    model.lp_.col_lower_.assign(species_count + 1, 0.0);
    model.lp_.col_upper_.assign(species_count + 1, kHighsInf);
    model.lp_.row_lower_.reserve(row_count);
    model.lp_.row_upper_.reserve(row_count);
    for (double total : constraints.totals) {
        model.lp_.row_lower_.push_back(total);
        model.lp_.row_upper_.push_back(total);
    }
    model.lp_.row_lower_.insert(
        model.lp_.row_lower_.end(), species_count, 0.0
    );
    model.lp_.row_upper_.insert(
        model.lp_.row_upper_.end(), species_count, kHighsInf
    );
    if (ion_limited) {
        model.lp_.row_lower_.push_back(-kHighsInf);
        model.lp_.row_upper_.push_back(0.0);
    }
    model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
    model.lp_.a_matrix_.start_ = {0};
    for (std::size_t species = 0; species < species_count; ++species) {
        for (std::size_t row = 0; row < constraints.matrix.rows; ++row) {
            const double value = constraints.matrix(row, species);
            if (value != 0.0) {
                model.lp_.a_matrix_.index_.push_back(static_cast<HighsInt>(row));
                model.lp_.a_matrix_.value_.push_back(value);
            }
        }
        model.lp_.a_matrix_.index_.push_back(static_cast<HighsInt>(
            constraints.matrix.rows + species
        ));
        model.lp_.a_matrix_.value_.push_back(1.0);
        if (ion_limited) {
            model.lp_.a_matrix_.index_.push_back(static_cast<HighsInt>(
                row_count - 1
            ));
            model.lp_.a_matrix_.value_.push_back(
                charges[species] == 0
                    ? -ion_fraction_max
                    : 1.0 - ion_fraction_max
            );
        }
        model.lp_.a_matrix_.start_.push_back(static_cast<HighsInt>(
            model.lp_.a_matrix_.index_.size()
        ));
    }
    for (std::size_t species = 0; species < species_count; ++species) {
        model.lp_.a_matrix_.index_.push_back(static_cast<HighsInt>(
            constraints.matrix.rows + species
        ));
        model.lp_.a_matrix_.value_.push_back(-1.0);
    }
    model.lp_.a_matrix_.start_.push_back(static_cast<HighsInt>(
        model.lp_.a_matrix_.index_.size()
    ));
    Highs highs;
    if (highs.setOptionValue("output_flag", false) == HighsStatus::kError
        || highs.setOptionValue("threads", 1) == HighsStatus::kError
        || highs.setOptionValue("solver", std::string("simplex"))
            == HighsStatus::kError
        || highs.passModel(model) == HighsStatus::kError
        || highs.run() == HighsStatus::kError) {
        return {"setup_failed", {}};
    }
    if (highs.getModelStatus() != HighsModelStatus::kOptimal
        || highs.getInfo().primal_solution_status != kSolutionStatusFeasible) {
        return {highs.modelStatusToString(highs.getModelStatus()), {}};
    }
    return {"solve_succeeded", highs.getSolution().col_value};
}

struct PhysicalPhaseEvaluation {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    double pressure_pa = 0.0;
};

struct PhysicalScalarEvaluation {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
};

struct PhaseBlockEvaluation {
    PhysicalPhaseEvaluation mechanical;
    PhysicalScalarEvaluation packing;
    bool has_packing = false;
    std::vector<std::string> active_parameter_names;
    std::vector<double> state_parameter_derivatives;
    std::vector<double> pressure_parameter_derivatives_pa;
    std::vector<double> chemical_potential_parameter_derivatives_over_rt;
};

using PhaseEvaluator = std::function<PhaseBlockEvaluation(
    double,
    const std::vector<double>&,
    double
)>;

struct ReactionDomain {
    bool enforce_packing = false;
    double packing_min = 0.0;
    double packing_max = 0.0;
    double total_ion_fraction_max = std::numeric_limits<double>::quiet_NaN();
};

PhysicalPhaseEvaluation evaluate_ideal_phase(
    double temperature_k,
    const std::vector<double>& amounts,
    double volume
) {
    if (!std::isfinite(volume) || volume <= 0.0
        || !std::all_of(amounts.begin(), amounts.end(), [](double amount) {
            return std::isfinite(amount) && amount > 0.0;
        })) {
        throw std::domain_error("ideal phase requires positive finite amounts and volume");
    }
    const std::size_t coordinate_count = amounts.size() + 1;
    PhysicalPhaseEvaluation result;
    result.gradient.assign(coordinate_count, 0.0);
    result.hessian.assign(coordinate_count * coordinate_count, 0.0);
    const double total = std::accumulate(amounts.begin(), amounts.end(), 0.0);
    for (std::size_t species = 0; species < amounts.size(); ++species) {
        const double logarithm = std::log(amounts[species] / volume);
        result.value += amounts[species] * (logarithm - 1.0);
        result.gradient[species] = logarithm;
        result.hessian[species * coordinate_count + species] = 1.0 / amounts[species];
        result.hessian[species * coordinate_count + amounts.size()] = -1.0 / volume;
        result.hessian[amounts.size() * coordinate_count + species] = -1.0 / volume;
    }
    result.gradient.back() = -total / volume;
    result.hessian.back() = total / (volume * volume);
    result.pressure_pa = kGasConstantJPerMolK * temperature_k * total / volume;
    return result;
}

PhaseBlockEvaluation evaluate_ideal_phase_block(
    double temperature_k,
    const std::vector<double>& amounts,
    double volume
) {
    PhaseBlockEvaluation result;
    result.mechanical = evaluate_ideal_phase(temperature_k, amounts, volume);
    return result;
}

PhaseEvaluator ideal_phase_evaluator() {
    return evaluate_ideal_phase_block;
}

std::size_t reaction_constraint_count(
    const ConstraintRows& balances,
    const ReactionDomain& domain
) {
    return balances.matrix.rows + (domain.enforce_packing ? 1 : 0)
        + (std::isfinite(domain.total_ion_fraction_max) ? 1 : 0);
}

struct ReactionNlpEvaluation {
    double objective = 0.0;
    std::vector<double> gradient;
    std::vector<double> constraints;
    std::vector<double> jacobian;
    std::vector<double> lagrangian_hessian;
    AmountChartEvaluation amount_chart;
    PhaseBlockEvaluation phase;
    std::vector<double> constraint_lower;
    std::vector<double> constraint_upper;
    double volume = 0.0;
};

ReactionNlpEvaluation evaluate_reaction_nlp(
    const AmountChart& chart,
    const ConstraintRows& balances,
    const std::vector<double>& g_ref,
    double temperature_k,
    double pressure_pa,
    const PhaseEvaluator& phase_evaluator,
    const ReactionDomain& domain,
    const std::vector<double>& variables,
    const std::vector<double>& multipliers
) {
    const std::size_t amount_dimension = chart.coordinate_count();
    const std::size_t variable_count = amount_dimension + 1;
    const std::size_t constraint_count = reaction_constraint_count(balances, domain);
    if (variables.size() != variable_count || multipliers.size() != constraint_count) {
        throw std::invalid_argument("reaction NLP dimensions are inconsistent");
    }
    std::vector<double> amount_coordinates(
        variables.begin(), variables.begin() + static_cast<std::ptrdiff_t>(amount_dimension)
    );
    ReactionNlpEvaluation result;
    result.amount_chart = evaluate_amount_chart(chart, amount_coordinates);
    result.volume = std::exp(variables.back());
    if (!std::isfinite(result.volume) || result.volume <= 0.0) {
        throw std::domain_error("reaction NLP volume is invalid");
    }
    result.phase = phase_evaluator(
        temperature_k, result.amount_chart.amounts, result.volume
    );
    if (result.phase.mechanical.gradient.size() != g_ref.size() + 1
        || result.phase.mechanical.hessian.size() != (g_ref.size() + 1) * (g_ref.size() + 1)) {
        throw std::invalid_argument("phase block derivative dimensions are inconsistent");
    }
    const double pressure_over_rt = pressure_pa / (kGasConstantJPerMolK * temperature_k);
    result.objective = result.phase.mechanical.value + pressure_over_rt * result.volume;
    for (std::size_t species = 0; species < g_ref.size(); ++species) {
        result.objective += g_ref[species] * result.amount_chart.amounts[species];
    }

    std::vector<double> physical_gradient = result.phase.mechanical.gradient;
    for (std::size_t species = 0; species < g_ref.size(); ++species) {
        physical_gradient[species] += g_ref[species];
    }
    physical_gradient.back() += pressure_over_rt;
    result.gradient.assign(variable_count, 0.0);
    for (std::size_t reduced = 0; reduced < amount_dimension; ++reduced) {
        for (std::size_t species = 0; species < g_ref.size(); ++species) {
            result.gradient[reduced] += result.amount_chart.jacobian[
                species * amount_dimension + reduced
            ] * physical_gradient[species];
        }
    }
    result.gradient.back() = result.volume * physical_gradient.back();

    result.constraints.assign(constraint_count, 0.0);
    result.constraint_lower.assign(constraint_count, 0.0);
    result.constraint_upper.assign(constraint_count, 0.0);
    std::vector<double> physical_constraint_gradients(
        constraint_count * (g_ref.size() + 1), 0.0
    );
    std::vector<double> physical_constraint_hessians(
        constraint_count * (g_ref.size() + 1) * (g_ref.size() + 1), 0.0
    );
    for (std::size_t row = 0; row < balances.matrix.rows; ++row) {
        result.constraints[row] = -balances.totals[row];
        for (std::size_t species = 0; species < g_ref.size(); ++species) {
            result.constraints[row] += balances.matrix(row, species)
                * result.amount_chart.amounts[species];
            physical_constraint_gradients[row * (g_ref.size() + 1) + species] =
                balances.matrix(row, species);
        }
    }

    std::size_t constraint = balances.matrix.rows;
    if (domain.enforce_packing) {
        if (!result.phase.has_packing
            || result.phase.packing.gradient.size() != g_ref.size() + 1
            || result.phase.packing.hessian.size()
                != (g_ref.size() + 1) * (g_ref.size() + 1)) {
            throw std::invalid_argument("phase block is missing exact packing derivatives");
        }
        result.constraints[constraint] = result.phase.packing.value;
        result.constraint_lower[constraint] = domain.packing_min;
        result.constraint_upper[constraint] = domain.packing_max;
        std::copy(
            result.phase.packing.gradient.begin(),
            result.phase.packing.gradient.end(),
            physical_constraint_gradients.begin()
                + static_cast<std::ptrdiff_t>(constraint * (g_ref.size() + 1))
        );
        std::copy(
            result.phase.packing.hessian.begin(),
            result.phase.packing.hessian.end(),
            physical_constraint_hessians.begin()
                + static_cast<std::ptrdiff_t>(
                    constraint * (g_ref.size() + 1) * (g_ref.size() + 1)
                )
        );
        ++constraint;
    }
    if (std::isfinite(domain.total_ion_fraction_max)) {
        result.constraint_lower[constraint] = -kInfinity;
        result.constraint_upper[constraint] = 0.0;
        for (std::size_t species = 0; species < g_ref.size(); ++species) {
            const double coefficient = chart.charges[species] == 0
                ? -domain.total_ion_fraction_max
                : 1.0 - domain.total_ion_fraction_max;
            result.constraints[constraint] += coefficient
                * result.amount_chart.amounts[species];
            physical_constraint_gradients[constraint * (g_ref.size() + 1) + species] =
                coefficient;
        }
    }

    result.jacobian.assign(constraint_count * variable_count, 0.0);

    std::vector<double> physical_jacobian((g_ref.size() + 1) * variable_count, 0.0);
    for (std::size_t species = 0; species < g_ref.size(); ++species) {
        for (std::size_t reduced = 0; reduced < amount_dimension; ++reduced) {
            physical_jacobian[species * variable_count + reduced] =
                result.amount_chart.jacobian[species * amount_dimension + reduced];
        }
    }
    physical_jacobian[g_ref.size() * variable_count + amount_dimension] = result.volume;
    for (std::size_t row = 0; row < constraint_count; ++row) {
        for (std::size_t reduced = 0; reduced < amount_dimension; ++reduced) {
            for (std::size_t species = 0; species < g_ref.size(); ++species) {
                result.jacobian[row * variable_count + reduced] +=
                    physical_constraint_gradients[row * (g_ref.size() + 1) + species]
                    * result.amount_chart.jacobian[species * amount_dimension + reduced];
            }
        }
        result.jacobian[row * variable_count + amount_dimension] = result.volume
            * physical_constraint_gradients[row * (g_ref.size() + 1) + g_ref.size()];
    }

    std::vector<double> physical_lagrangian_gradient = physical_gradient;
    std::vector<double> physical_lagrangian_hessian = result.phase.mechanical.hessian;
    for (std::size_t constraint_row = 0; constraint_row < constraint_count; ++constraint_row) {
        for (std::size_t physical = 0; physical < g_ref.size() + 1; ++physical) {
            physical_lagrangian_gradient[physical] += multipliers[constraint_row]
                * physical_constraint_gradients[
                    constraint_row * (g_ref.size() + 1) + physical
                ];
            for (std::size_t right = 0; right < g_ref.size() + 1; ++right) {
                physical_lagrangian_hessian[physical * (g_ref.size() + 1) + right] +=
                    multipliers[constraint_row]
                    * physical_constraint_hessians[
                        constraint_row * (g_ref.size() + 1) * (g_ref.size() + 1)
                            + physical * (g_ref.size() + 1) + right
                    ];
            }
        }
    }
    result.lagrangian_hessian.assign(variable_count * variable_count, 0.0);
    for (std::size_t row = 0; row < variable_count; ++row) {
        for (std::size_t column = 0; column < variable_count; ++column) {
            for (std::size_t left = 0; left < g_ref.size() + 1; ++left) {
                for (std::size_t right = 0; right < g_ref.size() + 1; ++right) {
                    result.lagrangian_hessian[row * variable_count + column] +=
                        physical_jacobian[left * variable_count + row]
                        * physical_lagrangian_hessian[left * (g_ref.size() + 1) + right]
                        * physical_jacobian[right * variable_count + column];
                }
            }
        }
    }
    for (std::size_t species = 0; species < g_ref.size(); ++species) {
        const double component_weight = physical_lagrangian_gradient[species];
        for (std::size_t row = 0; row < amount_dimension; ++row) {
            for (std::size_t column = 0; column < amount_dimension; ++column) {
                result.lagrangian_hessian[row * variable_count + column] +=
                    component_weight * result.amount_chart.amount_hessians[
                        species * amount_dimension * amount_dimension
                            + row * amount_dimension + column
                    ];
            }
        }
    }
    result.lagrangian_hessian.back() += physical_lagrangian_gradient.back() * result.volume;
    return result;
}

std::vector<double> chemical_potentials(
    const ReactionNlpEvaluation& evaluation,
    const std::vector<double>& g_ref
) {
    std::vector<double> result(
        evaluation.phase.mechanical.gradient.begin(),
        evaluation.phase.mechanical.gradient.begin()
            + static_cast<std::ptrdiff_t>(g_ref.size())
    );
    for (std::size_t species = 0; species < result.size(); ++species) {
        result[species] += g_ref[species];
    }
    return result;
}

struct KktPolishEvaluation {
    std::vector<double> constraints;
    std::vector<double> jacobian;
    bool domain_feasible = true;
};

std::size_t kkt_polish_constraint_count(
    const ConstraintRows& balances,
    const DenseMatrix& reactions
) {
    return balances.matrix.rows + reactions.rows + 1;
}

KktPolishEvaluation evaluate_kkt_polish(
    const AmountChart& chart,
    const ConstraintRows& balances,
    const DenseMatrix& reactions,
    const std::vector<double>& g_ref,
    double temperature_k,
    double pressure_pa,
    const PhaseEvaluator& phase_evaluator,
    const ReactionDomain& domain,
    const std::vector<double>& variables
) {
    if (reactions.columns != g_ref.size()) {
        throw std::invalid_argument("reaction polish dimensions are inconsistent");
    }
    const ReactionNlpEvaluation base = evaluate_reaction_nlp(
        chart,
        balances,
        g_ref,
        temperature_k,
        pressure_pa,
        phase_evaluator,
        domain,
        variables,
        std::vector<double>(reaction_constraint_count(balances, domain), 0.0)
    );
    const std::size_t species_count = g_ref.size();
    const std::size_t variable_count = variables.size();
    const std::size_t amount_dimension = chart.coordinate_count();
    const std::size_t physical_count = species_count + 1;
    KktPolishEvaluation result;
    const std::size_t constraint_count = kkt_polish_constraint_count(balances, reactions);
    result.constraints.assign(constraint_count, 0.0);
    result.jacobian.assign(constraint_count * variable_count, 0.0);
    for (std::size_t row = balances.matrix.rows; row < base.constraints.size(); ++row) {
        const double scale = std::max(1.0, std::abs(base.constraints[row]));
        result.domain_feasible = result.domain_feasible
            && base.constraints[row]
                >= base.constraint_lower[row] - kProviderDomainTolerance * scale
            && base.constraints[row]
                <= base.constraint_upper[row] + kProviderDomainTolerance * scale;
    }

    for (std::size_t row = 0; row < balances.matrix.rows; ++row) {
        const double scale = std::max(1.0, std::abs(balances.totals[row]));
        result.constraints[row] = base.constraints[row] / scale;
        for (std::size_t variable = 0; variable < variable_count; ++variable) {
            result.jacobian[row * variable_count + variable] =
                base.jacobian[row * variable_count + variable] / scale;
        }
    }

    std::vector<double> physical_jacobian(physical_count * variable_count, 0.0);
    for (std::size_t species = 0; species < species_count; ++species) {
        for (std::size_t reduced = 0; reduced < amount_dimension; ++reduced) {
            physical_jacobian[species * variable_count + reduced] =
                base.amount_chart.jacobian[species * amount_dimension + reduced];
        }
    }
    physical_jacobian[species_count * variable_count + amount_dimension] =
        base.volume;

    const std::vector<double> potentials = chemical_potentials(base, g_ref);
    std::size_t constraint = balances.matrix.rows;
    for (std::size_t reaction = 0; reaction < reactions.rows; ++reaction) {
        const double reaction_scale = matrix_row_l2_norm(reactions, reaction);
        for (std::size_t species = 0; species < species_count; ++species) {
            result.constraints[constraint] += reactions(reaction, species)
                * potentials[species] / reaction_scale;
            for (std::size_t variable = 0; variable < variable_count; ++variable) {
                for (std::size_t physical = 0; physical < physical_count; ++physical) {
                    result.jacobian[constraint * variable_count + variable] +=
                        reactions(reaction, species)
                        * base.phase.mechanical.hessian[
                            species * physical_count + physical
                        ]
                        * physical_jacobian[physical * variable_count + variable]
                        / reaction_scale;
                }
            }
        }
        ++constraint;
    }

    result.constraints[constraint] =
        (base.phase.mechanical.pressure_pa - pressure_pa) / pressure_pa;
    for (std::size_t variable = 0; variable < variable_count; ++variable) {
        for (std::size_t physical = 0; physical < physical_count; ++physical) {
            result.jacobian[constraint * variable_count + variable] -=
                kGasConstantJPerMolK * temperature_k
                * base.phase.mechanical.hessian[
                    species_count * physical_count + physical
                ]
                * physical_jacobian[physical * variable_count + variable]
                / pressure_pa;
        }
    }
    return result;
}

bool solve_square_system(
    std::vector<double> matrix,
    std::vector<double>& right_hand_side
) {
    const std::size_t dimension = right_hand_side.size();
    if (matrix.size() != dimension * dimension) {
        return false;
    }
    double matrix_scale = 0.0;
    for (double value : matrix) {
        matrix_scale = std::max(matrix_scale, std::abs(value));
    }
    const double tolerance = 4096.0 * std::numeric_limits<double>::epsilon()
        * std::max(1.0, matrix_scale) * static_cast<double>(dimension);
    for (std::size_t column = 0; column < dimension; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < dimension; ++row) {
            if (std::abs(matrix[row * dimension + column])
                > std::abs(matrix[pivot * dimension + column])) {
                pivot = row;
            }
        }
        if (std::abs(matrix[pivot * dimension + column]) <= tolerance) {
            return false;
        }
        if (pivot != column) {
            for (std::size_t entry = column; entry < dimension; ++entry) {
                std::swap(
                    matrix[pivot * dimension + entry],
                    matrix[column * dimension + entry]
                );
            }
            std::swap(right_hand_side[pivot], right_hand_side[column]);
        }
        for (std::size_t row = column + 1; row < dimension; ++row) {
            const double factor = matrix[row * dimension + column]
                / matrix[column * dimension + column];
            for (std::size_t entry = column + 1; entry < dimension; ++entry) {
                matrix[row * dimension + entry] -=
                    factor * matrix[column * dimension + entry];
            }
            right_hand_side[row] -= factor * right_hand_side[column];
        }
    }
    for (std::size_t row = dimension; row-- > 0;) {
        for (std::size_t column = row + 1; column < dimension; ++column) {
            right_hand_side[row] -=
                matrix[row * dimension + column] * right_hand_side[column];
        }
        right_hand_side[row] /= matrix[row * dimension + row];
    }
    return std::all_of(
        right_hand_side.begin(),
        right_hand_side.end(),
        [](double value) { return std::isfinite(value); }
    );
}

struct SquareSystemAnalysis {
    std::size_t rank = 0;
    double condition_number_inf = std::numeric_limits<double>::infinity();
};

SquareSystemAnalysis analyze_square_system(const std::vector<double>& matrix) {
    const std::size_t dimension = static_cast<std::size_t>(
        std::sqrt(static_cast<double>(matrix.size()))
    );
    if (dimension * dimension != matrix.size() || dimension == 0) {
        return {};
    }
    double matrix_norm_inf = 0.0;
    double matrix_scale = 0.0;
    for (std::size_t row = 0; row < dimension; ++row) {
        double row_sum = 0.0;
        for (std::size_t column = 0; column < dimension; ++column) {
            const double value = std::abs(matrix[row * dimension + column]);
            row_sum += value;
            matrix_scale = std::max(matrix_scale, value);
        }
        matrix_norm_inf = std::max(matrix_norm_inf, row_sum);
    }
    const double tolerance = 4096.0 * std::numeric_limits<double>::epsilon()
        * std::max(1.0, matrix_scale) * static_cast<double>(dimension);
    std::vector<double> reduced = matrix;
    std::size_t rank = 0;
    for (std::size_t column = 0; column < dimension; ++column) {
        std::size_t pivot = rank;
        for (std::size_t row = rank; row < dimension; ++row) {
            if (std::abs(reduced[row * dimension + column])
                > std::abs(reduced[pivot * dimension + column])) {
                pivot = row;
            }
        }
        if (std::abs(reduced[pivot * dimension + column]) <= tolerance) {
            continue;
        }
        if (pivot != rank) {
            for (std::size_t entry = column; entry < dimension; ++entry) {
                std::swap(
                    reduced[pivot * dimension + entry],
                    reduced[rank * dimension + entry]
                );
            }
        }
        for (std::size_t row = rank + 1; row < dimension; ++row) {
            const double factor = reduced[row * dimension + column]
                / reduced[rank * dimension + column];
            for (std::size_t entry = column + 1; entry < dimension; ++entry) {
                reduced[row * dimension + entry] -=
                    factor * reduced[rank * dimension + entry];
            }
        }
        ++rank;
    }
    SquareSystemAnalysis result{rank, std::numeric_limits<double>::infinity()};
    if (rank != dimension) {
        return result;
    }
    std::vector<double> inverse_row_sums(dimension, 0.0);
    for (std::size_t column = 0; column < dimension; ++column) {
        std::vector<double> inverse_column(dimension, 0.0);
        inverse_column[column] = 1.0;
        if (!solve_square_system(matrix, inverse_column)) {
            return result;
        }
        for (std::size_t row = 0; row < dimension; ++row) {
            inverse_row_sums[row] += std::abs(inverse_column[row]);
        }
    }
    result.condition_number_inf = matrix_norm_inf * *std::max_element(
        inverse_row_sums.begin(), inverse_row_sums.end()
    );
    return result;
}

ChemicalSensitivityResult evaluate_implicit_sensitivities(
    const AmountChart& chart,
    const ConstraintRows& balances,
    const DenseMatrix& reactions,
    const std::vector<double>& g_ref,
    double temperature_k,
    double pressure_pa,
    const PhaseEvaluator& phase_evaluator,
    const ReactionDomain& domain,
    const std::vector<double>& variables,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    const std::vector<std::size_t>& active_constraint_bounds,
    double trace_floor,
    const std::string& parameter_fingerprint,
    const std::vector<double>& ln_k_pressure_derivatives_per_pa,
    const std::vector<double>& ln_k_parameter_derivatives
) {
    constexpr double kInactiveMargin = 1.0e-7;
    ChemicalSensitivityResult result;
    result.parameter_fingerprint = parameter_fingerprint;
    result.chart_topology = chart.ionic()
        ? "electroneutral_log_total+softmax_shares+neutral_log_amounts+log_volume"
        : "neutral_log_amounts[" + std::to_string(chart.coordinate_count())
            + "]+log_volume";
    for (std::size_t variable = 0; variable < variables.size(); ++variable) {
        if (variables[variable] - lower[variable] <= kInactiveMargin) {
            result.active_lower_bounds.push_back(variable);
        }
        if (upper[variable] - variables[variable] <= kInactiveMargin) {
            result.active_upper_bounds.push_back(variable);
        }
    }
    result.active_constraint_bounds = active_constraint_bounds;
    KktPolishEvaluation evaluation;
    try {
        evaluation = evaluate_kkt_polish(
            chart,
            balances,
            reactions,
            g_ref,
            temperature_k,
            pressure_pa,
            phase_evaluator,
            domain,
            variables
        );
    } catch (const std::exception&) {
        result.failure_reason = "kkt_evaluation_failed";
        return result;
    }
    result.kkt_dimension = variables.size();
    const SquareSystemAnalysis analysis = analyze_square_system(evaluation.jacobian);
    result.kkt_rank = analysis.rank;
    result.condition_number_inf = analysis.condition_number_inf;
    if (!result.active_lower_bounds.empty()
        || !result.active_upper_bounds.empty()
        || !result.active_constraint_bounds.empty()) {
        result.failure_reason = "active_set_change_not_differentiable";
        return result;
    }
    if (!evaluation.domain_feasible) {
        result.failure_reason = "provider_domain_boundary_active";
        return result;
    }
    if (analysis.rank != result.kkt_dimension) {
        result.failure_reason = "singular_kkt_jacobian";
        return result;
    }
    if (!std::isfinite(analysis.condition_number_inf)
        || analysis.condition_number_inf > kSensitivityConditionNumberMax) {
        result.failure_reason = "ill_conditioned_kkt_jacobian";
        return result;
    }
    if (!ln_k_pressure_derivatives_per_pa.empty()
        && ln_k_pressure_derivatives_per_pa.size() != reactions.rows) {
        result.failure_reason =
            "incomplete_transformed_reference_pressure_derivatives";
        return result;
    }
    const std::size_t amount_dimension = chart.coordinate_count();
    const AmountChartEvaluation amount_chart = evaluate_amount_chart(
        chart,
        std::vector<double>(
            variables.begin(),
            variables.begin() + static_cast<std::ptrdiff_t>(amount_dimension)
        )
    );
    const double volume = std::exp(variables.back());
    PhaseBlockEvaluation phase;
    try {
        phase = phase_evaluator(
            temperature_k, amount_chart.amounts, volume
        );
    } catch (const std::exception&) {
        result.failure_reason = "provider_parameter_partial_evaluation_failed";
        return result;
    }
    const std::size_t active_parameter_count =
        phase.active_parameter_names.size();
    if (active_parameter_count != 0
        && (phase.state_parameter_derivatives.size()
                != (g_ref.size() + 1) * active_parameter_count
            || phase.pressure_parameter_derivatives_pa.size()
                != active_parameter_count
            || phase.chemical_potential_parameter_derivatives_over_rt.size()
                != g_ref.size() * active_parameter_count
            || ln_k_parameter_derivatives.size()
                != reactions.rows * active_parameter_count)) {
        result.failure_reason =
            "incomplete_provider_or_reference_parameter_derivatives";
        return result;
    }
    if (active_parameter_count == 0 && !ln_k_parameter_derivatives.empty()) {
        result.failure_reason = "provider_parameter_order_mismatch";
        return result;
    }

    for (std::size_t row = 0; row < balances.matrix.rows; ++row) {
        result.parameter_order.push_back("balance_total[" + std::to_string(row) + "]");
    }
    for (std::size_t row = 0; row < reactions.rows; ++row) {
        result.parameter_order.push_back(
            "ln_k_provider_basis[" + std::to_string(row) + "]"
        );
    }
    result.parameter_order.emplace_back("pressure_pa");
    result.parameter_order.insert(
        result.parameter_order.end(),
        phase.active_parameter_names.begin(),
        phase.active_parameter_names.end()
    );
    const std::size_t parameter_count = result.parameter_order.size();
    result.amount_derivatives.assign(
        parameter_count * g_ref.size(), 0.0
    );
    result.volume_derivatives.assign(parameter_count, 0.0);
    for (std::size_t species = 0;
         species < amount_chart.amounts.size();
         ++species) {
        if (amount_chart.amounts[species] <= trace_floor) {
            result.active_trace_species.push_back(species);
        }
    }
    if (!result.active_trace_species.empty()) {
        result.failure_reason = "active_set_change_not_differentiable";
        result.parameter_order.clear();
        result.amount_derivatives.clear();
        result.volume_derivatives.clear();
        return result;
    }
    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        std::vector<double> right_hand_side(result.kkt_dimension, 0.0);
        if (parameter < balances.matrix.rows) {
            const double total = balances.totals[parameter];
            const double scale = std::max(1.0, std::abs(total));
            right_hand_side[parameter] = 1.0 / scale;
            if (std::abs(total) > 1.0) {
                right_hand_side[parameter] +=
                    evaluation.constraints[parameter] * std::copysign(1.0, total)
                    / scale;
            }
        } else if (parameter < balances.matrix.rows + reactions.rows) {
            const std::size_t reaction = parameter - balances.matrix.rows;
            right_hand_side[parameter] =
                1.0 / matrix_row_l2_norm(reactions, reaction);
        } else if (parameter
            == balances.matrix.rows + reactions.rows) {
            for (std::size_t reaction = 0; reaction < reactions.rows; ++reaction) {
                right_hand_side[balances.matrix.rows + reaction] =
                    ln_k_pressure_derivatives_per_pa.empty()
                    ? 0.0
                    : ln_k_pressure_derivatives_per_pa[reaction]
                        / matrix_row_l2_norm(reactions, reaction);
            }
            right_hand_side.back() =
                (1.0 + evaluation.constraints.back()) / pressure_pa;
        } else {
            const std::size_t active = parameter
                - balances.matrix.rows - reactions.rows - 1;
            for (std::size_t reaction = 0; reaction < reactions.rows; ++reaction) {
                double mechanical_derivative = 0.0;
                for (std::size_t species = 0; species < g_ref.size(); ++species) {
                    mechanical_derivative += reactions(reaction, species)
                        * phase.chemical_potential_parameter_derivatives_over_rt[
                            species * active_parameter_count + active
                        ];
                }
                right_hand_side[balances.matrix.rows + reaction] = (
                    ln_k_parameter_derivatives[
                        reaction * active_parameter_count + active
                    ] - mechanical_derivative
                ) / matrix_row_l2_norm(reactions, reaction);
            }
            right_hand_side.back() =
                -phase.pressure_parameter_derivatives_pa[active] / pressure_pa;
        }
        if (!solve_square_system(evaluation.jacobian, right_hand_side)) {
            result.failure_reason = "singular_kkt_jacobian";
            result.amount_derivatives.clear();
            result.volume_derivatives.clear();
            return result;
        }
        for (std::size_t species = 0; species < g_ref.size(); ++species) {
            for (std::size_t coordinate = 0;
                 coordinate < amount_dimension;
                 ++coordinate) {
                result.amount_derivatives[parameter * g_ref.size() + species] +=
                    amount_chart.jacobian[
                        species * amount_dimension + coordinate
                    ] * right_hand_side[coordinate];
            }
        }
        result.volume_derivatives[parameter] =
            volume * right_hand_side.back();
    }
    result.status = "available";
    result.failure_reason.clear();
    return result;
}

bool polish_interior_kkt(
    const AmountChart& chart,
    const ConstraintRows& balances,
    const DenseMatrix& reactions,
    const std::vector<double>& g_ref,
    double temperature_k,
    double pressure_pa,
    const PhaseEvaluator& phase_evaluator,
    const ReactionDomain& domain,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    int max_iterations,
    std::vector<double>& variables
) {
    const std::size_t dimension = variables.size();
    if (kkt_polish_constraint_count(balances, reactions) != dimension) {
        return false;
    }
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        const KktPolishEvaluation evaluation = evaluate_kkt_polish(
            chart,
            balances,
            reactions,
            g_ref,
            temperature_k,
            pressure_pa,
            phase_evaluator,
            domain,
            variables
        );
        const double residual = vector_inf_norm(evaluation.constraints);
        if (residual <= kAffinityTolerance) {
            return evaluation.domain_feasible;
        }
        std::vector<double> step = evaluation.constraints;
        for (double& value : step) {
            value = -value;
        }
        if (!solve_square_system(evaluation.jacobian, step)) {
            return false;
        }
        double alpha = 1.0;
        for (std::size_t index = 0; index < dimension; ++index) {
            if (step[index] > 0.0) {
                alpha = std::min(
                    alpha, 0.99 * (upper[index] - variables[index]) / step[index]
                );
            } else if (step[index] < 0.0) {
                alpha = std::min(
                    alpha, 0.99 * (lower[index] - variables[index]) / step[index]
                );
            }
        }
        bool accepted_step = false;
        for (int backtrack = 0; backtrack < 60 && alpha > 1.0e-12; ++backtrack) {
            std::vector<double> candidate = variables;
            for (std::size_t index = 0; index < dimension; ++index) {
                candidate[index] += alpha * step[index];
            }
            try {
                const KktPolishEvaluation trial = evaluate_kkt_polish(
                    chart,
                    balances,
                    reactions,
                    g_ref,
                    temperature_k,
                    pressure_pa,
                    phase_evaluator,
                    domain,
                    candidate
                );
                if (trial.domain_feasible
                    && vector_inf_norm(trial.constraints) < residual) {
                    variables = std::move(candidate);
                    accepted_step = true;
                    break;
                }
            } catch (const std::exception&) {
            }
            alpha *= 0.5;
        }
        if (!accepted_step) {
            return false;
        }
    }
    return false;
}

class ReactionTnlp final : public Ipopt::TNLP {
public:
    ReactionTnlp(
        AmountChart chart,
        ConstraintRows balances,
        std::vector<double> g_ref,
        double temperature_k,
        double pressure_pa,
        PhaseEvaluator phase_evaluator,
        ReactionDomain domain,
        std::vector<double> initial,
        std::vector<double> lower,
        std::vector<double> upper
    )
        : chart_(std::move(chart)),
          balances_(std::move(balances)),
          g_ref_(std::move(g_ref)),
          temperature_k_(temperature_k),
          pressure_pa_(pressure_pa),
          phase_evaluator_(std::move(phase_evaluator)),
          domain_(domain),
          initial_(std::move(initial)),
          lower_(std::move(lower)),
          upper_(std::move(upper)),
          solution_(initial_),
          constraint_multipliers_(reaction_constraint_count(balances_, domain_), 0.0),
          lower_multipliers_(initial_.size(), 0.0),
          upper_multipliers_(initial_.size(), 0.0) {}

    bool get_nlp_info(
        Ipopt::Index& n,
        Ipopt::Index& m,
        Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag,
        IndexStyleEnum& index_style
    ) override {
        n = static_cast<Ipopt::Index>(initial_.size());
        m = static_cast<Ipopt::Index>(constraint_multipliers_.size());
        nnz_jac_g = static_cast<Ipopt::Index>(initial_.size() * constraint_multipliers_.size());
        nnz_h_lag = static_cast<Ipopt::Index>(initial_.size() * (initial_.size() + 1) / 2);
        index_style = TNLP::C_STYLE;
        return true;
    }

    bool get_bounds_info(
        Ipopt::Index n,
        Ipopt::Number* x_l,
        Ipopt::Number* x_u,
        Ipopt::Index m,
        Ipopt::Number* g_l,
        Ipopt::Number* g_u
    ) override {
        if (n != static_cast<Ipopt::Index>(initial_.size())
            || m != static_cast<Ipopt::Index>(constraint_multipliers_.size())
            || x_l == nullptr
            || x_u == nullptr
            || g_l == nullptr
            || g_u == nullptr) {
            return false;
        }
        std::copy(lower_.begin(), lower_.end(), x_l);
        std::copy(upper_.begin(), upper_.end(), x_u);
        std::fill(g_l, g_l + m, 0.0);
        std::fill(g_u, g_u + m, 0.0);
        std::size_t constraint = balances_.matrix.rows;
        if (domain_.enforce_packing) {
            g_l[constraint] = domain_.packing_min;
            g_u[constraint] = domain_.packing_max;
            ++constraint;
        }
        if (std::isfinite(domain_.total_ion_fraction_max)) {
            g_l[constraint] = -kInfinity;
            g_u[constraint] = 0.0;
        }
        return true;
    }

    bool get_starting_point(
        Ipopt::Index n,
        bool init_x,
        Ipopt::Number* x,
        bool,
        Ipopt::Number*,
        Ipopt::Number*,
        Ipopt::Index m,
        bool,
        Ipopt::Number*
    ) override {
        if (n != static_cast<Ipopt::Index>(initial_.size())
            || m != static_cast<Ipopt::Index>(constraint_multipliers_.size())
            || !init_x
            || x == nullptr) {
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
        if (n != static_cast<Ipopt::Index>(initial_.size()) || x == nullptr) {
            return false;
        }
        try {
            objective = evaluate(x, std::vector<double>(constraint_multipliers_.size(), 0.0)).objective;
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            objective = 0.0;
            return true;
        }
    }

    bool eval_grad_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number* gradient
    ) override {
        if (n != static_cast<Ipopt::Index>(initial_.size())
            || x == nullptr
            || gradient == nullptr) {
            return false;
        }
        try {
            const ReactionNlpEvaluation evaluation = evaluate(
                x, std::vector<double>(constraint_multipliers_.size(), 0.0)
            );
            std::copy(evaluation.gradient.begin(), evaluation.gradient.end(), gradient);
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            std::fill(gradient, gradient + n, 0.0);
            return true;
        }
    }

    bool eval_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Index m,
        Ipopt::Number* constraints
    ) override {
        if (n != static_cast<Ipopt::Index>(initial_.size())
            || m != static_cast<Ipopt::Index>(constraint_multipliers_.size())
            || x == nullptr
            || constraints == nullptr) {
            return false;
        }
        try {
            const ReactionNlpEvaluation evaluation = evaluate(
                x, std::vector<double>(constraint_multipliers_.size(), 0.0)
            );
            std::copy(evaluation.constraints.begin(), evaluation.constraints.end(), constraints);
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            std::fill(constraints, constraints + m, 0.0);
            return true;
        }
    }

    bool eval_jac_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Index m,
        Ipopt::Index nonzero_count,
        Ipopt::Index* rows,
        Ipopt::Index* columns,
        Ipopt::Number* values
    ) override {
        if (n != static_cast<Ipopt::Index>(initial_.size())
            || m != static_cast<Ipopt::Index>(constraint_multipliers_.size())
            || nonzero_count != n * m) {
            return false;
        }
        if (values == nullptr) {
            if (rows == nullptr || columns == nullptr) {
                return false;
            }
            for (Ipopt::Index row = 0; row < m; ++row) {
                for (Ipopt::Index column = 0; column < n; ++column) {
                    rows[row * n + column] = row;
                    columns[row * n + column] = column;
                }
            }
            return true;
        }
        if (x == nullptr) {
            return false;
        }
        try {
            const ReactionNlpEvaluation evaluation = evaluate(
                x, std::vector<double>(constraint_multipliers_.size(), 0.0)
            );
            std::copy(evaluation.jacobian.begin(), evaluation.jacobian.end(), values);
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            std::fill(values, values + nonzero_count, 0.0);
            return true;
        }
    }

    bool eval_h(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number objective_factor,
        Ipopt::Index m,
        const Ipopt::Number* multipliers,
        bool,
        Ipopt::Index nonzero_count,
        Ipopt::Index* rows,
        Ipopt::Index* columns,
        Ipopt::Number* values
    ) override {
        if (n != static_cast<Ipopt::Index>(initial_.size())
            || m != static_cast<Ipopt::Index>(constraint_multipliers_.size())
            || nonzero_count != n * (n + 1) / 2) {
            return false;
        }
        if (values == nullptr) {
            if (rows == nullptr || columns == nullptr) {
                return false;
            }
            Ipopt::Index entry = 0;
            for (Ipopt::Index row = 0; row < n; ++row) {
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    rows[entry] = row;
                    columns[entry] = column;
                    ++entry;
                }
            }
            return true;
        }
        if (x == nullptr || multipliers == nullptr) {
            return false;
        }
        try {
            std::vector<double> lambda(multipliers, multipliers + m);
            const ReactionNlpEvaluation objective = evaluate(
                x, std::vector<double>(static_cast<std::size_t>(m), 0.0)
            );
            const ReactionNlpEvaluation lagrangian = evaluate(x, lambda);
            Ipopt::Index entry = 0;
            for (Ipopt::Index row = 0; row < n; ++row) {
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    const std::size_t index = static_cast<std::size_t>(row * n + column);
                    values[entry++] = objective_factor * objective.lagrangian_hessian[index]
                        + lagrangian.lagrangian_hessian[index]
                        - objective.lagrangian_hessian[index];
                }
            }
            return true;
        } catch (const std::exception& error) {
            callback_error_ = error.what();
            std::fill(values, values + nonzero_count, 0.0);
            return true;
        }
    }

    bool intermediate_callback(
        Ipopt::AlgorithmMode,
        Ipopt::Index,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Index,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*
    ) override {
        // Ipopt 3.11 may read callback buffers after a false return. Every
        // caught Provider failure writes finite containment values instead;
        // this callback stops the solve and the recorded error rejects the
        // result, so the containment values can never become scientific data.
        return callback_error_.empty();
    }

    void finalize_solution(
        Ipopt::SolverReturn,
        Ipopt::Index n,
        const Ipopt::Number* x,
        const Ipopt::Number* z_l,
        const Ipopt::Number* z_u,
        Ipopt::Index m,
        const Ipopt::Number*,
        const Ipopt::Number* lambda,
        Ipopt::Number,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*
    ) override {
        if (n == static_cast<Ipopt::Index>(solution_.size())
            && x != nullptr
            && z_l != nullptr
            && z_u != nullptr) {
            std::copy(x, x + n, solution_.begin());
            std::copy(z_l, z_l + n, lower_multipliers_.begin());
            std::copy(z_u, z_u + n, upper_multipliers_.begin());
        }
        if (m == static_cast<Ipopt::Index>(constraint_multipliers_.size())
            && lambda != nullptr) {
            std::copy(lambda, lambda + m, constraint_multipliers_.begin());
        }
    }

    [[nodiscard]] const std::vector<double>& solution() const { return solution_; }
    [[nodiscard]] const std::vector<double>& constraint_multipliers() const {
        return constraint_multipliers_;
    }
    [[nodiscard]] const std::vector<double>& lower_multipliers() const {
        return lower_multipliers_;
    }
    [[nodiscard]] const std::vector<double>& upper_multipliers() const {
        return upper_multipliers_;
    }
    [[nodiscard]] const std::string& callback_error() const { return callback_error_; }

private:
    ReactionNlpEvaluation evaluate(
        const Ipopt::Number* values,
        const std::vector<double>& multipliers
    ) const {
        return evaluate_reaction_nlp(
            chart_,
            balances_,
            g_ref_,
            temperature_k_,
            pressure_pa_,
            phase_evaluator_,
            domain_,
            std::vector<double>(values, values + initial_.size()),
            multipliers
        );
    }

    AmountChart chart_;
    ConstraintRows balances_;
    std::vector<double> g_ref_;
    double temperature_k_;
    double pressure_pa_;
    PhaseEvaluator phase_evaluator_;
    ReactionDomain domain_;
    std::vector<double> initial_;
    std::vector<double> lower_;
    std::vector<double> upper_;
    std::vector<double> solution_;
    std::vector<double> constraint_multipliers_;
    std::vector<double> lower_multipliers_;
    std::vector<double> upper_multipliers_;
    std::string callback_error_;
};

std::vector<std::vector<double>> nullspace_basis(
    const std::vector<double>& matrix,
    std::size_t rows,
    std::size_t columns
) {
    std::vector<double> rref = matrix;
    std::vector<std::size_t> pivot_columns;
    std::size_t pivot_row = 0;
    for (std::size_t column = 0; column < columns && pivot_row < rows; ++column) {
        std::size_t pivot = pivot_row;
        for (std::size_t row = pivot_row + 1; row < rows; ++row) {
            if (std::abs(rref[row * columns + column])
                > std::abs(rref[pivot * columns + column])) {
                pivot = row;
            }
        }
        if (std::abs(rref[pivot * columns + column]) <= 1.0e-10) {
            continue;
        }
        for (std::size_t entry = 0; entry < columns; ++entry) {
            std::swap(rref[pivot * columns + entry], rref[pivot_row * columns + entry]);
        }
        const double pivot_value = rref[pivot_row * columns + column];
        for (std::size_t entry = 0; entry < columns; ++entry) {
            rref[pivot_row * columns + entry] /= pivot_value;
        }
        for (std::size_t row = 0; row < rows; ++row) {
            if (row == pivot_row) {
                continue;
            }
            const double factor = rref[row * columns + column];
            for (std::size_t entry = 0; entry < columns; ++entry) {
                rref[row * columns + entry] -= factor
                    * rref[pivot_row * columns + entry];
            }
        }
        pivot_columns.push_back(column);
        ++pivot_row;
    }
    std::vector<bool> is_pivot(columns, false);
    for (std::size_t column : pivot_columns) {
        is_pivot[column] = true;
    }
    std::vector<std::vector<double>> basis;
    for (std::size_t free_column = 0; free_column < columns; ++free_column) {
        if (is_pivot[free_column]) {
            continue;
        }
        std::vector<double> vector(columns, 0.0);
        vector[free_column] = 1.0;
        for (std::size_t row = 0; row < pivot_columns.size(); ++row) {
            vector[pivot_columns[row]] = -rref[row * columns + free_column];
        }
        for (const std::vector<double>& prior : basis) {
            const double projection = std::inner_product(
                vector.begin(), vector.end(), prior.begin(), 0.0
            );
            for (std::size_t entry = 0; entry < columns; ++entry) {
                vector[entry] -= projection * prior[entry];
            }
        }
        const double norm = std::sqrt(std::inner_product(
            vector.begin(), vector.end(), vector.begin(), 0.0
        ));
        for (double& value : vector) {
            value /= norm;
        }
        basis.push_back(std::move(vector));
    }
    return basis;
}

std::vector<double> recompute_equality_multipliers(
    const ConstraintRows& balances,
    const std::vector<double>& potentials
) {
    const std::size_t equality_count = balances.matrix.rows;
    if (balances.matrix.columns != potentials.size()) {
        throw std::invalid_argument("equality multiplier dimensions are inconsistent");
    }
    std::vector<double> gram(equality_count * equality_count, 0.0);
    std::vector<double> right_hand_side(equality_count, 0.0);
    for (std::size_t row = 0; row < equality_count; ++row) {
        for (std::size_t other = 0; other < equality_count; ++other) {
            for (std::size_t species = 0; species < potentials.size(); ++species) {
                gram[row * equality_count + other] +=
                    balances.matrix(row, species)
                    * balances.matrix(other, species);
            }
        }
        for (std::size_t species = 0; species < potentials.size(); ++species) {
            right_hand_side[row] -=
                balances.matrix(row, species) * potentials[species];
        }
    }
    if (!solve_square_system(std::move(gram), right_hand_side)) {
        throw std::domain_error("equality multiplier system is singular");
    }
    return right_hand_side;
}

bool reduced_hessian_positive(
    const ReactionNlpEvaluation& evaluation,
    std::size_t constraint_count
) {
    const std::size_t dimension = evaluation.gradient.size();
    const std::vector<std::vector<double>> basis = nullspace_basis(
        evaluation.jacobian, constraint_count, dimension
    );
    if (basis.empty()) {
        return true;
    }
    std::vector<double> reduced(basis.size() * basis.size(), 0.0);
    for (std::size_t row = 0; row < basis.size(); ++row) {
        for (std::size_t column = 0; column < basis.size(); ++column) {
            for (std::size_t left = 0; left < dimension; ++left) {
                for (std::size_t right = 0; right < dimension; ++right) {
                    reduced[row * basis.size() + column] += basis[row][left]
                        * evaluation.lagrangian_hessian[left * dimension + right]
                        * basis[column][right];
                }
            }
        }
    }
    std::vector<double> diagonal_scales(basis.size(), 0.0);
    for (std::size_t index = 0; index < basis.size(); ++index) {
        const double diagonal = reduced[index * basis.size() + index];
        if (!std::isfinite(diagonal) || diagonal <= 0.0) {
            return false;
        }
        diagonal_scales[index] = std::sqrt(diagonal);
    }
    for (std::size_t row = 0; row < basis.size(); ++row) {
        for (std::size_t column = 0; column < basis.size(); ++column) {
            reduced[row * basis.size() + column] /=
                diagonal_scales[row] * diagonal_scales[column];
        }
    }
    for (std::size_t column = 0; column < basis.size(); ++column) {
        double diagonal = reduced[column * basis.size() + column];
        for (std::size_t prior = 0; prior < column; ++prior) {
            const double value = reduced[column * basis.size() + prior];
            diagonal -= value * value;
        }
        if (diagonal <= 1.0e-10) {
            return false;
        }
        reduced[column * basis.size() + column] = std::sqrt(diagonal);
        for (std::size_t row = column + 1; row < basis.size(); ++row) {
            double value = reduced[row * basis.size() + column];
            for (std::size_t prior = 0; prior < column; ++prior) {
                value -= reduced[row * basis.size() + prior]
                    * reduced[column * basis.size() + prior];
            }
            reduced[row * basis.size() + column] =
                value / reduced[column * basis.size() + column];
        }
    }
    return true;
}

}  // namespace

MaxMinInitializationResult max_min_initialization(
    const DenseMatrix& balance_matrix,
    const std::vector<double>& feed_amounts,
    const std::vector<int>& charges,
    double trace_floor,
    double total_ion_fraction_max
) {
    if (balance_matrix.columns != feed_amounts.size()
        || charges.size() != feed_amounts.size()
        || balance_matrix.values.size() != balance_matrix.rows * balance_matrix.columns
        || !std::isfinite(trace_floor) || trace_floor <= 0.0
        || (std::isfinite(total_ion_fraction_max)
            && (total_ion_fraction_max < 0.0 || total_ion_fraction_max > 1.0))) {
        throw std::invalid_argument("max-min initialization input is invalid");
    }
    const ConstraintRows constraints = independent_max_min_rows(
        balance_matrix, feed_amounts, charges
    );
    MaxMinInitializationResult result;
    const LinearInitialization initialization = solve_initialization_lp(
        constraints, charges, total_ion_fraction_max
    );
    result.solver_status = initialization.status;
    if (initialization.values.size() != feed_amounts.size() + 1) {
        return result;
    }
    const std::vector<double>& solution = initialization.values;
    result.amounts.assign(solution.begin(), solution.begin() + static_cast<std::ptrdiff_t>(feed_amounts.size()));
    const double max_min_amount = solution.back();
    const std::vector<double> equality = matrix_vector(constraints.matrix, result.amounts);
    std::vector<double> residual(equality.size(), 0.0);
    for (std::size_t row = 0; row < equality.size(); ++row) {
        residual[row] = equality[row] - constraints.totals[row];
    }
    const double equality_inf_norm = vector_inf_norm(residual);
    const bool finite_solution = std::isfinite(max_min_amount)
        && max_min_amount < 0.5 * kInfinity
        && std::all_of(result.amounts.begin(), result.amounts.end(), [](double amount) {
            return std::isfinite(amount) && amount >= 0.0;
        });
    const double recomputed_minimum = result.amounts.empty()
        ? 0.0
        : *std::min_element(result.amounts.begin(), result.amounts.end());
    bool max_min_inequalities_feasible = finite_solution;
    for (double amount : result.amounts) {
        max_min_inequalities_feasible = max_min_inequalities_feasible
            && amount + kBalanceTolerance * std::max(1.0, std::abs(amount))
                >= max_min_amount;
    }
    bool source_domain_feasible = true;
    if (std::isfinite(total_ion_fraction_max)) {
        double ionic = 0.0;
        double total = 0.0;
        for (std::size_t species = 0; species < result.amounts.size(); ++species) {
            total += result.amounts[species];
            if (charges[species] != 0) {
                ionic += result.amounts[species];
            }
        }
        source_domain_feasible = total > 0.0
            && ionic / total <= total_ion_fraction_max + kProviderDomainTolerance;
    }
    result.amount_upper_bounds.assign(feed_amounts.size(), 0.0);
    const double bounding_total = constraints.totals.front();
    bool finite_amount_bounds = bounding_total > 0.0;
    for (std::size_t species = 0; species < feed_amounts.size(); ++species) {
        const double coefficient = constraints.matrix(0, species);
        finite_amount_bounds = finite_amount_bounds && coefficient > 0.0;
        result.amount_upper_bounds[species] = std::nextafter(
            2.0 * bounding_total / coefficient + 1.0e-12,
            std::numeric_limits<double>::infinity()
        );
    }
    result.strict_positive_feasible = result.solver_status == "solve_succeeded"
        && equality_inf_norm <= kBalanceTolerance
        && recomputed_minimum > trace_floor
        && max_min_amount > trace_floor
        && max_min_inequalities_feasible
        && source_domain_feasible
        && finite_amount_bounds;
    return result;
}

ChemicalSolveResult solve_reaction(
    const CompiledReactionSystem& system,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& gauge_coefficients,
    double trace_floor,
    int max_iterations,
    const MaxMinInitializationResult& initialization,
    const PhaseEvaluator& phase_evaluator,
    const ReactionDomain& domain,
    double initial_volume,
    const std::vector<double>& ln_k_pressure_derivatives_per_pa,
    const std::vector<double>& ln_k_parameter_derivatives,
    const std::array<double, 2>& volume_bounds
) {
    if (!std::isfinite(temperature_k) || temperature_k <= 0.0
        || !std::isfinite(pressure_pa) || pressure_pa <= 0.0
        || !std::isfinite(trace_floor) || trace_floor <= 0.0) {
        throw std::invalid_argument("reaction solve scales are invalid");
    }
    ChemicalSolveResult result;
    result.solver_status = initialization.solver_status;
    if (!initialization.strict_positive_feasible) {
        result.numerical_status = "failed";
        result.physical_status = "not_adjudicated";
        result.trace_status = "at_or_below_floor";
        return result;
    }
    std::vector<double> g_ref = system.g_ref;
    if (!gauge_coefficients.empty()) {
        if (gauge_coefficients.size() != system.balance_matrix.rows) {
            throw std::invalid_argument("gauge coefficient count does not match balances");
        }
        for (std::size_t species = 0; species < g_ref.size(); ++species) {
            for (std::size_t row = 0; row < system.balance_matrix.rows; ++row) {
                g_ref[species] += system.balance_matrix(row, species)
                    * gauge_coefficients[row];
            }
        }
    }
    const AmountChart chart = make_amount_chart(system.charges);
    const std::vector<double>& initial_amounts = initialization.amounts;
    if (initial_amounts.size() != system.species_count
        || !std::all_of(initial_amounts.begin(), initial_amounts.end(), [](double amount) {
            return std::isfinite(amount) && amount > 0.0;
        })) {
        throw std::invalid_argument("reaction starting amounts are invalid");
    }
    std::vector<double> initial = invert_amount_chart(chart, initial_amounts);
    if (!std::isfinite(initial_volume) || initial_volume <= 0.0) {
        const double pressure_over_rt = pressure_pa / (kGasConstantJPerMolK * temperature_k);
        initial_volume = std::accumulate(
            initialization.amounts.begin(), initialization.amounts.end(), 0.0
        ) / pressure_over_rt;
    }
    initial.push_back(std::log(initial_volume));
    std::vector<double> lower(initial.size(), -40.0);
    std::vector<double> upper(initial.size(), 40.0);
    if (!chart.ionic()) {
        std::fill(lower.begin(), lower.end() - 1, std::log(0.1 * trace_floor));
        for (std::size_t species = 0; species < system.species_count; ++species) {
            upper[species] = std::log(initialization.amount_upper_bounds[species]);
        }
    } else {
        lower[0] = std::log(0.1 * trace_floor);
        double charge_equivalent_upper = 0.0;
        for (std::size_t species : chart.cation_indices) {
            charge_equivalent_upper += static_cast<double>(system.charges[species])
                * initialization.amount_upper_bounds[species];
        }
        upper[0] = std::log(charge_equivalent_upper);
        if (!std::isfinite(charge_equivalent_upper) || charge_equivalent_upper <= 0.0) {
            throw std::invalid_argument("ionic charge-equivalent upper bound is invalid");
        }
        const std::size_t neutral_offset = 1 + chart.cation_indices.size() - 1
            + chart.anion_indices.size() - 1;
        // The charged-group coordinates are log(share/reference_share).
        // Static +/-40 bounds are artificial floors when trace_floor is below
        // exp(-40) mol. Since each share is at most one and charge_equivalents
        // is bounded above, these conservative bounds keep any species above
        // trace_floor away from the coordinate bounds while retaining the
        // neutral coordinates' ten-fold numerical margin.
        const double log_trace_share_floor = std::log(0.1 * trace_floor)
            - std::log(charge_equivalent_upper);
        const auto share_upper = [&](std::size_t reference_species) {
            return std::log(10.0 * charge_equivalent_upper)
                - std::log(trace_floor)
                - std::log(static_cast<double>(std::abs(system.charges[reference_species])));
        };
        std::size_t cation_offset = 1;
        for (std::size_t category = 0; category + 1 < chart.cation_indices.size(); ++category) {
            const std::size_t species = chart.cation_indices[category];
            lower[cation_offset + category] = log_trace_share_floor
                + std::log(static_cast<double>(std::abs(system.charges[species])));
            upper[cation_offset + category] = share_upper(chart.cation_indices.back());
        }
        const std::size_t anion_offset = cation_offset + chart.cation_indices.size() - 1;
        for (std::size_t category = 0; category + 1 < chart.anion_indices.size(); ++category) {
            const std::size_t species = chart.anion_indices[category];
            lower[anion_offset + category] = log_trace_share_floor
                + std::log(static_cast<double>(std::abs(system.charges[species])));
            upper[anion_offset + category] = share_upper(chart.anion_indices.back());
        }
        for (std::size_t neutral = 0; neutral < chart.neutral_indices.size(); ++neutral) {
            lower[neutral_offset + neutral] = std::log(0.1 * trace_floor);
            upper[neutral_offset + neutral] = std::log(
                initialization.amount_upper_bounds[chart.neutral_indices[neutral]]
            );
        }
    }
    if (std::isfinite(volume_bounds[0])
        && std::isfinite(volume_bounds[1])) {
        if (volume_bounds[0] <= 0.0
            || volume_bounds[1] <= volume_bounds[0]
            || initial_volume <= volume_bounds[0]
            || initial_volume >= volume_bounds[1]) {
            throw std::invalid_argument(
                "reaction volume bounds or initial volume are incompatible"
            );
        }
        lower.back() = std::log(volume_bounds[0]);
        upper.back() = std::log(volume_bounds[1]);
    } else {
        lower.back() = std::log(initial_volume) - 30.0;
        upper.back() = std::log(initial_volume) + 30.0;
    }
    const ConstraintRows balances{system.balance_matrix, system.balance_totals};
    auto* raw_problem = new ReactionTnlp(
        chart,
        balances,
        g_ref,
        temperature_k,
        pressure_pa,
        phase_evaluator,
        domain,
        initial,
        lower,
        upper
    );
    Ipopt::SmartPtr<Ipopt::TNLP> problem = raw_problem;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application = IpoptApplicationFactory();
    if (max_iterations < 0) {
        throw std::invalid_argument("reaction solver iteration limit must be nonnegative");
    }
    configure_ipopt(application, max_iterations);
    const Ipopt::ApplicationReturnStatus initialize_status = application->Initialize();
    if (initialize_status != Ipopt::Solve_Succeeded) {
        result.solver_status = "initialization_" + ipopt_status_name(initialize_status);
        result.numerical_status = "failed";
        return result;
    }
    const Ipopt::ApplicationReturnStatus status = application->OptimizeTNLP(problem);
    std::vector<double> variables = raw_problem->solution();
    bool polished = false;
    if (status == Ipopt::Solve_Succeeded && max_iterations > 0) {
        try {
            const ReactionNlpEvaluation preliminary = evaluate_reaction_nlp(
                chart,
                balances,
                g_ref,
                temperature_k,
                pressure_pa,
                phase_evaluator,
                domain,
                variables,
                raw_problem->constraint_multipliers()
            );
            const std::vector<double> preliminary_potentials =
                chemical_potentials(preliminary, g_ref);
            const double preliminary_affinity = reaction_residual_inf_norm(
                system.reaction_matrix, preliminary_potentials
            );
            const double preliminary_balance = vector_inf_norm(std::vector<double>(
                preliminary.constraints.begin(),
                preliminary.constraints.begin()
                    + static_cast<std::ptrdiff_t>(balances.matrix.rows)
            ));
            const double preliminary_pressure = std::abs(
                preliminary.phase.mechanical.pressure_pa - pressure_pa
            ) / pressure_pa;
            if (preliminary_affinity > kAffinityTolerance
                && preliminary_balance <= kBalanceTolerance
                && preliminary_pressure <= kPressureTolerance
                && preliminary.amount_chart.minimum_amount > trace_floor) {
                std::vector<double> candidate = variables;
                if (polish_interior_kkt(
                    chart,
                    balances,
                    system.reaction_matrix,
                    g_ref,
                    temperature_k,
                    pressure_pa,
                    phase_evaluator,
                    domain,
                    lower,
                    upper,
                    max_iterations,
                    candidate
                )) {
                    const ReactionNlpEvaluation polished_evaluation =
                        evaluate_reaction_nlp(
                            chart,
                            balances,
                            g_ref,
                            temperature_k,
                            pressure_pa,
                            phase_evaluator,
                            domain,
                            candidate,
                            std::vector<double>(
                                reaction_constraint_count(balances, domain), 0.0
                            )
                        );
                    const std::vector<double> polished_potentials =
                        chemical_potentials(polished_evaluation, g_ref);
                    const double polished_affinity = reaction_residual_inf_norm(
                        system.reaction_matrix, polished_potentials
                    );
                    if (polished_affinity <= kAffinityTolerance) {
                        variables = std::move(candidate);
                        polished = true;
                    }
                }
            }
        } catch (const std::exception&) {
            // The original minimization remains authoritative and fails its
            // physical certificates when an interior polish is unavailable.
        }
    }
    result.solver_status = ipopt_status_name(status);
    result.callback_error = raw_problem->callback_error();
    ReactionNlpEvaluation evaluation;
    try {
        evaluation = evaluate_reaction_nlp(
            chart,
            balances,
            g_ref,
            temperature_k,
            pressure_pa,
            phase_evaluator,
            domain,
            variables,
            polished
                ? std::vector<double>(
                    reaction_constraint_count(balances, domain), 0.0
                )
                : raw_problem->constraint_multipliers()
        );
    } catch (const std::exception& error) {
        result.callback_error = error.what();
        result.numerical_status = "failed";
        return result;
    }
    result.amounts = evaluation.amount_chart.amounts;
    result.volume_m3 = evaluation.volume;
    result.balance_inf_norm = vector_inf_norm(std::vector<double>(
        evaluation.constraints.begin(),
        evaluation.constraints.begin() + static_cast<std::ptrdiff_t>(balances.matrix.rows)
    ));
    result.charge_inf_norm = std::abs(evaluation.amount_chart.charge_residual);
    result.pressure_relative_residual = std::abs(
        evaluation.phase.mechanical.pressure_pa - pressure_pa
    ) / pressure_pa;
    const std::vector<double> potentials = chemical_potentials(evaluation, g_ref);
    result.reaction_affinity_inf_norm = reaction_residual_inf_norm(
        system.reaction_matrix, potentials
    );
    if (domain.enforce_packing) {
        result.packing_fraction = evaluation.phase.packing.value;
        result.provider_domain_status = result.packing_fraction >= domain.packing_min
                && result.packing_fraction <= domain.packing_max
            ? "passed"
            : "failed";
    } else {
        result.provider_domain_status = "not_applicable";
    }
    if (std::isfinite(domain.total_ion_fraction_max)) {
        double ionic = 0.0;
        double total = 0.0;
        for (std::size_t species = 0; species < result.amounts.size(); ++species) {
            total += result.amounts[species];
            if (system.charges[species] != 0) {
                ionic += result.amounts[species];
            }
        }
        if (total <= 0.0 || ionic / total > domain.total_ion_fraction_max + 1.0e-12) {
            result.provider_domain_status = "failed";
        }
    }
    result.trace_status = evaluation.amount_chart.minimum_amount > trace_floor
        ? "interior"
        : "at_or_below_floor";
    double complementarity = 0.0;
    bool sensitivity_interior = true;
    std::vector<std::size_t> active_constraint_bounds;
    constexpr double kInactiveMargin = 1.0e-7;
    for (std::size_t variable = 0; variable < variables.size(); ++variable) {
        if (!polished) {
            complementarity = std::max(
                complementarity,
                std::abs(raw_problem->lower_multipliers()[variable]
                    * (variables[variable] - lower[variable]))
            );
            complementarity = std::max(
                complementarity,
                std::abs(raw_problem->upper_multipliers()[variable]
                    * (upper[variable] - variables[variable]))
            );
        }
        sensitivity_interior = sensitivity_interior
            && variables[variable] - lower[variable] > kInactiveMargin
            && upper[variable] - variables[variable] > kInactiveMargin;
    }
    for (std::size_t row = balances.matrix.rows; row < evaluation.constraints.size(); ++row) {
        const double value = evaluation.constraints[row];
        const double multiplier = polished
            ? 0.0
            : raw_problem->constraint_multipliers()[row];
        const bool has_lower = evaluation.constraint_lower[row] > -0.5 * kInfinity;
        const bool has_upper = evaluation.constraint_upper[row] < 0.5 * kInfinity;
        if (multiplier < 0.0) {
            complementarity = std::max(
                complementarity,
                has_lower
                    ? std::abs(multiplier * (value - evaluation.constraint_lower[row]))
                    : std::abs(multiplier)
            );
        } else if (multiplier > 0.0) {
            complementarity = std::max(
                complementarity,
                has_upper
                    ? std::abs(multiplier * (evaluation.constraint_upper[row] - value))
                    : std::abs(multiplier)
            );
        }
        const double scale = std::max(1.0, std::abs(value));
        sensitivity_interior = sensitivity_interior
            && (!has_lower
                || value - evaluation.constraint_lower[row] > kInactiveMargin * scale)
            && (!has_upper
                || evaluation.constraint_upper[row] - value > kInactiveMargin * scale);
        if ((has_lower
                && value - evaluation.constraint_lower[row]
                    <= kInactiveMargin * scale)
            || (has_upper
                && evaluation.constraint_upper[row] - value
                    <= kInactiveMargin * scale)) {
            active_constraint_bounds.push_back(row);
        }
    }
    result.complementarity_inf_norm = complementarity;
    std::vector<double> equality_multipliers(
        reaction_constraint_count(balances, domain), 0.0
    );
    const ConstraintRows physical_balances = independent_max_min_rows(
        system.balance_matrix, system.feed_amounts, system.charges
    );
    std::vector<double> physical_multipliers;
    ReactionNlpEvaluation equality_evaluation;
    try {
        const std::vector<double> recomputed =
            recompute_equality_multipliers(balances, potentials);
        physical_multipliers =
            recompute_equality_multipliers(physical_balances, potentials);
        std::copy(
            recomputed.begin(), recomputed.end(), equality_multipliers.begin()
        );
        equality_evaluation = evaluate_reaction_nlp(
            chart,
            balances,
            g_ref,
            temperature_k,
            pressure_pa,
            phase_evaluator,
            domain,
            variables,
            equality_multipliers
        );
    } catch (const std::exception& error) {
        result.callback_error = error.what();
        result.numerical_status = "failed";
        return result;
    }
    if (sensitivity_interior) {
        result.local_minimum_status = reduced_hessian_positive(
            equality_evaluation, balances.matrix.rows
        ) ? "passed" : "failed";
    } else {
        result.local_minimum_status = "not_adjudicated";
    }

    std::vector<double> physical_stationarity = potentials;
    for (std::size_t species = 0; species < physical_stationarity.size(); ++species) {
        for (std::size_t row = 0; row < physical_balances.matrix.rows; ++row) {
            physical_stationarity[species] +=
                physical_balances.matrix(row, species) * physical_multipliers[row];
        }
    }
    result.kkt_stationarity_inf_norm = vector_inf_norm(physical_stationarity);
    result.numerical_status = status == Ipopt::Solve_Succeeded
            && result.balance_inf_norm <= kBalanceTolerance
            && result.kkt_stationarity_inf_norm <= kKktTolerance
            && result.complementarity_inf_norm <= kKktTolerance
        ? "passed"
        : "failed";
    result.physical_status = result.balance_inf_norm <= kBalanceTolerance
            && result.charge_inf_norm <= kBalanceTolerance
            && result.pressure_relative_residual <= kPressureTolerance
            && result.reaction_affinity_inf_norm <= kAffinityTolerance
            && result.trace_status == "interior"
            && result.provider_domain_status != "failed"
        ? "passed"
        : "failed";
    result.accepted = result.solver_status == "solve_succeeded"
        && result.callback_error.empty()
        && result.numerical_status == "passed"
        && result.physical_status == "passed"
        && result.local_minimum_status == "passed";
    if (status == Ipopt::Solve_Succeeded && result.callback_error.empty()) {
        result.sensitivities = evaluate_implicit_sensitivities(
            chart,
            balances,
            system.reaction_matrix,
            g_ref,
            temperature_k,
            pressure_pa,
            phase_evaluator,
            domain,
            variables,
            lower,
            upper,
            active_constraint_bounds,
            trace_floor,
            system.provider_fingerprint,
            ln_k_pressure_derivatives_per_pa,
            ln_k_parameter_derivatives
        );
        if (!result.accepted && result.sensitivities.status == "available") {
            result.sensitivities.status = "unavailable";
            result.sensitivities.failure_reason = "primal_solution_not_certified";
            result.sensitivities.parameter_order.clear();
            result.sensitivities.amount_derivatives.clear();
            result.sensitivities.volume_derivatives.clear();
        }
    }
    return result;
}

namespace {

void attach_support_evidence(
    const CompiledReactionSystem& system,
    ChemicalSolveResult& result
) {
    result.structural_zero_species_indices = system.removed_species_indices;
}

void expand_original_amounts_and_residuals(
    const CompiledReactionSystem& system,
    ChemicalSolveResult& result
) {
    if (!system.removed_species_indices.empty()
        && result.amounts.size() == system.retained_species_indices.size()) {
        std::vector<double> original_amounts(
            system.original_species_count, 0.0
        );
        for (std::size_t retained = 0;
             retained < system.retained_species_indices.size();
             ++retained) {
            original_amounts[system.retained_species_indices[retained]] =
                result.amounts[retained];
        }
        result.amounts = std::move(original_amounts);
    }
    if (result.amounts.size() != system.original_species_count) {
        return;
    }
    for (std::size_t row = 0; row < system.supplied_balance_matrix.rows; ++row) {
        double residual = 0.0;
        for (std::size_t species = 0;
             species < system.original_species_count;
             ++species) {
            residual += system.supplied_balance_matrix(row, species)
                * (
                    result.amounts[species]
                    - system.original_feed_amounts[species]
                );
        }
        result.balance_inf_norm = std::max(
            result.balance_inf_norm, std::abs(residual)
        );
    }
    double mass_residual = 0.0;
    double charge_residual = 0.0;
    for (std::size_t species = 0;
         species < system.original_species_count;
         ++species) {
        mass_residual += system.original_molar_masses_kg_per_mol[species]
            * (
                result.amounts[species]
                - system.original_feed_amounts[species]
            );
        charge_residual += static_cast<double>(system.original_charges[species])
            * result.amounts[species];
    }
    result.balance_inf_norm = std::max(
        result.balance_inf_norm, std::abs(mass_residual)
    );
    result.charge_inf_norm = std::max(
        result.charge_inf_norm, std::abs(charge_residual)
    );
}

ChemicalSolveResult finalize_chemical_result(
    const CompiledReactionSystem& system,
    ChemicalSolveResult result,
    bool structural_face_supported
) {
    attach_support_evidence(system, result);
    if (!system.removed_species_indices.empty()) {
        result.sensitivities.status = "unavailable";
        result.sensitivities.failure_reason =
            "active_set_change_not_differentiable";
        result.sensitivities.parameter_order.clear();
        result.sensitivities.amount_derivatives.clear();
        result.sensitivities.volume_derivatives.clear();
        result.sensitivities.active_trace_species =
            system.removed_species_indices;
        result.sensitivities.chart_topology = "structural_face";
    }
    if (!system.removed_species_indices.empty() && !structural_face_supported) {
        result.accepted = false;
        result.solver_status = "boundary_direction_unresolved";
        result.chemical_certification_level = "BOUNDARY_DIRECTION_UNRESOLVED";
        result.boundary_status = "boundary_direction_unresolved";
        result.amounts.clear();
        return result;
    }
    expand_original_amounts_and_residuals(system, result);
    if (result.accepted
        && (
            result.balance_inf_norm > kBalanceTolerance
            || result.charge_inf_norm > kBalanceTolerance
        )) {
        result.accepted = false;
        result.numerical_status = "failed";
        result.physical_status = "failed";
    }
    if (result.accepted) {
        result.chemical_certification_level = "LOCAL_EQUILIBRIUM";
        result.boundary_status = system.removed_species_indices.empty()
            ? "strict_interior"
            : "structural_face";
    } else if (!system.removed_species_indices.empty()) {
        result.boundary_status = "structural_face";
    }
    return result;
}

}  // namespace

ChemicalSolveResult solve_manufactured_ideal_reaction(
    const CompiledReactionSystem& system,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& gauge_coefficients,
    double trace_floor,
    int max_iterations
) {
    const MaxMinInitializationResult initialization = max_min_initialization(
        system.balance_matrix,
        system.feed_amounts,
        system.charges,
        trace_floor,
        std::numeric_limits<double>::quiet_NaN()
    );
    ChemicalSolveResult result = solve_reaction(
        system,
        temperature_k,
        pressure_pa,
        gauge_coefficients,
        trace_floor,
        max_iterations,
        initialization,
        ideal_phase_evaluator(),
        ReactionDomain{},
        std::numeric_limits<double>::quiet_NaN(),
        {},
        {},
        {
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
        }
    );
    return finalize_chemical_result(system, std::move(result), true);
}

ProviderPhaseBlockEvidence evaluate_provider_phase_block(
    const ProviderContext& provider,
    double temperature_k,
    const std::vector<double>& amounts,
    double volume_m3
) {
    const MixturePhaseEvaluation phase = provider.evaluate_electrolyte(
        temperature_k, amounts, volume_m3
    );
    const PackingFractionEvaluation packing = provider.evaluate_packing_fraction(
        temperature_k, amounts, volume_m3
    );
    return {
        phase.value,
        phase.gradient,
        phase.hessian,
        phase.pressure_pa,
        packing.value,
        packing.gradient,
        packing.hessian,
    };
}

ChemicalSolveResult solve_provider_reaction(
    const CompiledReactionSystem& system,
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    double packing_fraction_min,
    double packing_fraction_max,
    double total_ion_fraction_max,
    double trace_floor,
    const std::vector<double>& ln_k_pressure_derivatives_per_pa,
    const std::vector<double>& ln_k_parameter_derivatives,
    const ProviderActiveParameterSet* active_parameters,
    int max_iterations
) {
    if (!std::isfinite(packing_fraction_min)
        || !std::isfinite(packing_fraction_max)
        || packing_fraction_min <= 0.0
        || packing_fraction_max <= packing_fraction_min) {
        throw std::invalid_argument("packing-fraction bounds must be finite, positive, and ordered");
    }
    if (!system.removed_species_indices.empty()) {
        return finalize_chemical_result(system, ChemicalSolveResult{}, false);
    }
    const MaxMinInitializationResult initialization = max_min_initialization(
        system.balance_matrix,
        system.feed_amounts,
        system.charges,
        trace_floor,
        total_ion_fraction_max
    );
    if (!initialization.strict_positive_feasible) {
        ChemicalSolveResult result;
        result.solver_status = initialization.solver_status;
        result.numerical_status = "failed";
        result.trace_status = "at_or_below_floor";
        return finalize_chemical_result(system, std::move(result), false);
    }
    const std::vector<double>& starting_amounts = initialization.amounts;
    const double total = std::accumulate(
        starting_amounts.begin(), starting_amounts.end(), 0.0
    );
    std::vector<double> mole_fractions(starting_amounts.size(), 0.0);
    for (std::size_t species = 0; species < mole_fractions.size(); ++species) {
        mole_fractions[species] = starting_amounts[species] / total;
    }
    const std::array<double, 2> seed_molar_bounds = provider.evaluate_molar_volume_bounds(
        temperature_k, mole_fractions, packing_fraction_min, packing_fraction_max
    );
    const bool parameterized = active_parameters != nullptr
        && !active_parameters->parameters.empty();
    std::array<double, 2> molar_bounds = seed_molar_bounds;
    if (parameterized) {
        const double seed_volume = std::sqrt(
            total * seed_molar_bounds[0] * total * seed_molar_bounds[1]
        );
        molar_bounds = provider.evaluate_reacting_phase_parameters(
            temperature_k,
            starting_amounts,
            seed_volume,
            packing_fraction_min,
            packing_fraction_max,
            *active_parameters
        ).molar_volume_bounds_m3_per_mol;
    }
    double lower_volume = std::nextafter(
        total * molar_bounds[0], std::numeric_limits<double>::infinity()
    );
    double upper_volume = std::nextafter(total * molar_bounds[1], 0.0);
    const std::array<double, 2> admitted_volume_bounds{
        lower_volume, upper_volume
    };
    const auto pressure_residual = [&](double volume) {
        if (parameterized) {
            return provider.evaluate_reacting_phase_parameters(
                temperature_k,
                starting_amounts,
                volume,
                packing_fraction_min,
                packing_fraction_max,
                *active_parameters
            ).phase.pressure_pa - pressure_pa;
        }
        return provider.evaluate_electrolyte(
            temperature_k, starting_amounts, volume
        ).pressure_pa - pressure_pa;
    };
    double initial_volume = std::sqrt(lower_volume * upper_volume);
    double lower_residual = pressure_residual(lower_volume);
    const double upper_residual = pressure_residual(upper_volume);
    if (std::signbit(lower_residual) != std::signbit(upper_residual)) {
        for (int iteration = 0; iteration < 100; ++iteration) {
            const double midpoint = std::sqrt(lower_volume * upper_volume);
            const double midpoint_residual = pressure_residual(midpoint);
            initial_volume = midpoint;
            if (std::abs(midpoint_residual) <= 1.0e-10 * pressure_pa) {
                break;
            }
            if (std::signbit(midpoint_residual) == std::signbit(lower_residual)) {
                lower_volume = midpoint;
                lower_residual = midpoint_residual;
            } else {
                upper_volume = midpoint;
            }
        }
    }
    const PhaseEvaluator phase_evaluator = [
        &provider,
        active_parameters,
        parameterized,
        packing_fraction_min,
        packing_fraction_max
    ](
        double temperature,
        const std::vector<double>& amounts,
        double volume
    ) {
        if (parameterized) {
            const ParameterizedPhaseEvaluation block =
                provider.evaluate_reacting_phase_parameters(
                    temperature,
                    amounts,
                    volume,
                    packing_fraction_min,
                    packing_fraction_max,
                    *active_parameters
                );
            PhaseBlockEvaluation result;
            result.mechanical = {
                block.phase.value,
                block.phase.gradient,
                block.phase.hessian,
                block.phase.pressure_pa,
            };
            result.packing = {
                block.packing.value,
                block.packing.gradient,
                block.packing.hessian,
            };
            result.has_packing = true;
            result.state_parameter_derivatives =
                block.state_parameter_derivatives;
            result.pressure_parameter_derivatives_pa =
                block.pressure_parameter_derivatives_pa;
            result.chemical_potential_parameter_derivatives_over_rt =
                block.chemical_potential_parameter_derivatives_over_rt;
            result.active_parameter_names.reserve(
                active_parameters->parameters.size()
            );
            for (const ProviderActiveParameter& parameter
                 : active_parameters->parameters) {
                std::string name = "provider_parameter["
                    + parameter.family + ";" + parameter.identity + ";";
                for (std::size_t index = 0;
                     index < parameter.component_ids.size();
                     ++index) {
                    if (index != 0) {
                        name += ",";
                    }
                    name += parameter.component_ids[index];
                }
                result.active_parameter_names.push_back(name + "]");
            }
            return result;
        }
        const ProviderPhaseBlockEvidence block = evaluate_provider_phase_block(
            provider, temperature, amounts, volume
        );
        PhaseBlockEvaluation result;
        result.mechanical = {
            block.value,
            block.gradient,
            block.hessian,
            block.pressure_pa,
        };
        result.packing = {
            block.packing_fraction,
            block.packing_gradient,
            block.packing_hessian,
        };
        result.has_packing = true;
        return result;
    };
    const ReactionDomain domain{
        true,
        packing_fraction_min,
        packing_fraction_max,
        total_ion_fraction_max,
    };
    return finalize_chemical_result(
        system,
        solve_reaction(
            system,
            temperature_k,
            pressure_pa,
            {},
            trace_floor,
            max_iterations,
            initialization,
            phase_evaluator,
            domain,
            initial_volume,
            ln_k_pressure_derivatives_per_pa,
            ln_k_parameter_derivatives,
            admitted_volume_bounds
        ),
        false
    );
}

ManufacturedNlpEvaluation evaluate_manufactured_reaction_nlp(
    const CompiledReactionSystem& system,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& gauge_coefficients,
    const std::vector<double>& variables,
    const std::vector<double>& constraint_multipliers
) {
    std::vector<double> g_ref = system.g_ref;
    if (!gauge_coefficients.empty()) {
        if (gauge_coefficients.size() != system.balance_matrix.rows) {
            throw std::invalid_argument("gauge coefficient count does not match balances");
        }
        for (std::size_t species = 0; species < g_ref.size(); ++species) {
            for (std::size_t row = 0; row < system.balance_matrix.rows; ++row) {
                g_ref[species] += system.balance_matrix(row, species)
                    * gauge_coefficients[row];
            }
        }
    }
    const AmountChart chart = make_amount_chart(system.charges);
    const ConstraintRows balances{system.balance_matrix, system.balance_totals};
    const ReactionNlpEvaluation evaluation = evaluate_reaction_nlp(
        chart,
        balances,
        g_ref,
        temperature_k,
        pressure_pa,
        ideal_phase_evaluator(),
        ReactionDomain{},
        variables,
        constraint_multipliers
    );
    ManufacturedNlpEvaluation result;
    result.objective = evaluation.objective;
    result.objective_gradient = evaluation.gradient;
    result.constraints = evaluation.constraints;
    result.constraint_jacobian = evaluation.jacobian;
    result.lagrangian_gradient = evaluation.gradient;
    for (std::size_t variable = 0; variable < result.lagrangian_gradient.size(); ++variable) {
        for (std::size_t row = 0; row < balances.matrix.rows; ++row) {
            result.lagrangian_gradient[variable] += evaluation.jacobian[
                row * result.lagrangian_gradient.size() + variable
            ] * constraint_multipliers[row];
        }
    }
    result.lagrangian_hessian = evaluation.lagrangian_hessian;
    result.amounts = evaluation.amount_chart.amounts;
    result.volume_m3 = evaluation.volume;
    return result;
}

}  // namespace epcsaft_equilibrium
