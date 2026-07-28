#include "held2.hpp"
#include "held2_tolerances.hpp"

#include <Highs.h>
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft_equilibrium {
namespace {

std::vector<std::size_t> independent_positions(
    const Held2Coordinates& coordinates
) {
    std::vector<std::size_t> positions;
    for (std::size_t component : coordinates.independent_indices) {
        const auto found = std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            component
        );
        if (found == coordinates.retained_indices.end()) {
            throw std::invalid_argument(
                "HELD2 independent coordinate is not retained"
            );
        }
        positions.push_back(static_cast<std::size_t>(
            found - coordinates.retained_indices.begin()
        ));
    }
    return positions;
}

double dependent_upper(const Held2Coordinates& coordinates) {
    const auto found = std::find(
        coordinates.retained_indices.begin(),
        coordinates.retained_indices.end(),
        coordinates.dependent_index
    );
    if (found == coordinates.retained_indices.end()) {
        throw std::invalid_argument(
            "HELD2 dependent coordinate is not retained"
        );
    }
    const std::size_t position = static_cast<std::size_t>(
        found - coordinates.retained_indices.begin()
    );
    return 1.0
        - kHeld2ModifiedLowerScale * coordinates.modified_factors[position];
}

double maximum_abs(const std::vector<double>& values) {
    double maximum = 0.0;
    for (double value : values) {
        maximum = std::max(maximum, std::abs(value));
    }
    return maximum;
}

double maximum_difference(
    const std::vector<double>& left,
    const std::vector<double>& right
) {
    if (left.size() != right.size()) {
        return std::numeric_limits<double>::infinity();
    }
    double maximum = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        maximum = std::max(
            maximum, std::abs(left[index] - right[index])
        );
    }
    return maximum;
}

struct Problem67 {
    const Held2Coordinates& coordinates;
    const std::vector<double>& feed;
    const Held2StateEvaluator& evaluate;
    const Held2StateValueEvaluator& evaluate_value;
    std::size_t phase_count;
    std::size_t dimension;
    std::size_t block_size;
    int objective_evaluations = 0;
    std::string callback_error;
};

double problem67_value(
    const Problem67& problem,
    const std::vector<double>& variables
) {
    double objective = 0.0;
    for (std::size_t phase = 0; phase < problem.phase_count; ++phase) {
        const std::size_t offset = phase * problem.block_size;
        const std::vector<double> composition(
            variables.begin() + static_cast<std::ptrdiff_t>(offset + 1),
            variables.begin() + static_cast<std::ptrdiff_t>(
                offset + 1 + problem.dimension
            )
        );
        objective += variables[offset] * problem.evaluate_value(
            composition, variables[offset + 1 + problem.dimension]
        );
    }
    return objective;
}

double problem67_objective(
    const std::vector<double>& variables,
    std::vector<double>& gradient,
    void* opaque
) {
    auto& problem = *static_cast<Problem67*>(opaque);
    ++problem.objective_evaluations;
    if (!gradient.empty()) {
        std::fill(gradient.begin(), gradient.end(), 0.0);
    }
    double objective = 0.0;
    for (std::size_t phase = 0; phase < problem.phase_count; ++phase) {
        const std::size_t offset = phase * problem.block_size;
        const double fraction = variables[offset];
        const std::vector<double> composition(
            variables.begin() + static_cast<std::ptrdiff_t>(offset + 1),
            variables.begin() + static_cast<std::ptrdiff_t>(
                offset + 1 + problem.dimension
            )
        );
        Held2StateEvaluation state;
        try {
            state = problem.evaluate(
                composition, variables[offset + 1 + problem.dimension]
            );
        } catch (const std::exception& error) {
            problem.callback_error = error.what();
            throw;
        }
        if (state.gradient.size() != problem.dimension + 1) {
            throw std::invalid_argument(
                "HELD2 Problem (67) gradient dimensions changed"
            );
        }
        objective += fraction * state.objective;
        if (gradient.empty()) {
            continue;
        }
        gradient[offset] = state.objective;
        for (std::size_t local = 0;
             local < problem.dimension + 1;
             ++local) {
            gradient[offset + 1 + local] =
                fraction * state.gradient[local];
        }
    }
    return objective;
}

struct Problem67Constraint {
    Problem67* problem;
    int kind;
    std::size_t index;
    std::size_t phase = 0;
};

double problem67_constraint(
    const std::vector<double>& variables,
    std::vector<double>& gradient,
    void* opaque
) {
    const auto& constraint =
        *static_cast<Problem67Constraint*>(opaque);
    const Problem67& problem = *constraint.problem;
    if (!gradient.empty()) {
        std::fill(gradient.begin(), gradient.end(), 0.0);
    }
    double value = constraint.kind == 0
        ? -1.0
        : constraint.kind == 1
            ? -problem.feed[constraint.index]
            : -problem.coordinates
                .polytope_constraints[constraint.index].upper_bound;
    for (std::size_t phase = 0; phase < problem.phase_count; ++phase) {
        const std::size_t offset = phase * problem.block_size;
        const double fraction = variables[offset];
        if (constraint.kind == 0) {
            value += fraction;
            if (!gradient.empty()) {
                gradient[offset] = 1.0;
            }
        } else if (constraint.kind == 1) {
            const std::size_t variable = offset + 1 + constraint.index;
            value += fraction * variables[variable];
            if (!gradient.empty()) {
                gradient[offset] = variables[variable];
                gradient[variable] = fraction;
            }
        } else if (phase == constraint.phase) {
            const Held2PolytopeConstraint& polytope =
                problem.coordinates.polytope_constraints[constraint.index];
            for (std::size_t coordinate = 0;
                 coordinate < problem.dimension;
                 ++coordinate) {
                value += polytope.coefficients[coordinate]
                    * variables[offset + 1 + coordinate];
                if (!gradient.empty()) {
                    gradient[offset + 1 + coordinate] =
                        polytope.coefficients[coordinate];
                }
            }
        }
    }
    return value;
}

Held2Problem67Evaluation evaluate_problem67(
    const Held2Coordinates& coordinates,
    const std::vector<double>& feed,
    const Held2StateEvaluator& evaluator,
    std::size_t phase_count,
    const std::vector<double>& variables,
    const std::vector<double>& multipliers
);

class Problem67Tnlp final : public Ipopt::TNLP {
public:
    Problem67Tnlp(
        Problem67& problem,
        std::vector<double> lower,
        std::vector<double> upper,
        std::vector<double> initial
    )
        : problem_(problem),
          lower_(std::move(lower)),
          upper_(std::move(upper)),
          initial_(std::move(initial)) {
        constraints_.push_back({&problem_, 0, 0});
        for (std::size_t coordinate = 0;
             coordinate < problem_.dimension;
             ++coordinate) {
            constraints_.push_back({&problem_, 1, coordinate});
        }
        for (std::size_t phase = 0;
             phase < problem_.phase_count;
             ++phase) {
            for (std::size_t constraint = 0;
                 constraint
                    < problem_.coordinates.polytope_constraints.size();
                 ++constraint) {
                constraints_.push_back(
                    {&problem_, 2, constraint, phase}
                );
            }
        }
    }

    bool get_nlp_info(
        Ipopt::Index& n,
        Ipopt::Index& m,
        Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag,
        IndexStyleEnum& index_style
    ) override {
        n = static_cast<Ipopt::Index>(initial_.size());
        m = static_cast<Ipopt::Index>(constraints_.size());
        nnz_jac_g = n * m;
        nnz_h_lag = n * (n + 1) / 2;
        index_style = C_STYLE;
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
            || m != static_cast<Ipopt::Index>(constraints_.size())) {
            return false;
        }
        std::copy(lower_.begin(), lower_.end(), x_l);
        std::copy(upper_.begin(), upper_.end(), x_u);
        for (std::size_t index = 0; index < constraints_.size(); ++index) {
            g_l[index] = constraints_[index].kind < 2 ? 0.0 : -1.0e19;
            g_u[index] = 0.0;
        }
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
        if (n != static_cast<Ipopt::Index>(initial_.size())
            || m != static_cast<Ipopt::Index>(constraints_.size())
            || !init_x || init_z || init_lambda) {
            return false;
        }
        std::copy(initial_.begin(), initial_.end(), x);
        return true;
    }

    bool eval_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number& value
    ) override {
        if (!domain_valid(n, x)) {
            return false;
        }
        try {
            if (problem_.evaluate_value) {
                value = problem67_value(
                    problem_, std::vector<double>(x, x + n)
                );
            } else {
                evaluate(n, x);
                value = cached_evaluation_.objective;
            }
            return true;
        } catch (...) {
            return false;
        }
    }

    bool eval_grad_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number* gradient
    ) override {
        if (!domain_valid(n, x)) {
            return false;
        }
        try {
            evaluate(n, x);
            std::copy(
                cached_evaluation_.objective_gradient.begin(),
                cached_evaluation_.objective_gradient.end(),
                gradient
            );
            return true;
        } catch (...) {
            return false;
        }
    }

    bool eval_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Index m,
        Ipopt::Number* values
    ) override {
        if (n != static_cast<Ipopt::Index>(initial_.size())
            || m != static_cast<Ipopt::Index>(constraints_.size())) {
            return false;
        }
        const std::vector<double> variables(x, x + n);
        std::vector<double> no_gradient;
        for (std::size_t index = 0; index < constraints_.size(); ++index) {
            values[index] = problem67_constraint(
                variables, no_gradient, &constraints_[index]
            );
        }
        return true;
    }

    bool eval_jac_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Index m,
        Ipopt::Index nnz,
        Ipopt::Index* rows,
        Ipopt::Index* columns,
        Ipopt::Number* values
    ) override {
        if (m != static_cast<Ipopt::Index>(constraints_.size())
            || nnz != n * m) {
            return false;
        }
        Ipopt::Index position = 0;
        if (values == nullptr) {
            for (Ipopt::Index row = 0; row < m; ++row) {
                for (Ipopt::Index column = 0; column < n; ++column) {
                    rows[position] = row;
                    columns[position++] = column;
                }
            }
            return true;
        }
        const std::vector<double> variables(x, x + n);
        for (Problem67Constraint& constraint : constraints_) {
            std::vector<double> gradient(static_cast<std::size_t>(n));
            static_cast<void>(problem67_constraint(
                variables, gradient, &constraint
            ));
            for (double value : gradient) {
                values[position++] = value;
            }
        }
        return true;
    }

    bool eval_h(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number objective_factor,
        Ipopt::Index m,
        const Ipopt::Number* lambda,
        bool,
        Ipopt::Index nnz,
        Ipopt::Index* rows,
        Ipopt::Index* columns,
        Ipopt::Number* values
    ) override {
        if (m != static_cast<Ipopt::Index>(constraints_.size())
            || nnz != n * (n + 1) / 2) {
            return false;
        }
        Ipopt::Index position = 0;
        if (values == nullptr) {
            for (Ipopt::Index row = 0; row < n; ++row) {
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    rows[position] = row;
                    columns[position++] = column;
                }
            }
            return true;
        }
        if (!domain_valid(n, x) || lambda == nullptr) {
            return false;
        }
        try {
            evaluate(n, x);
            Held2Problem67Evaluation evaluation = cached_evaluation_;
            for (double& value : evaluation.lagrangian_hessian) {
                value *= objective_factor;
            }
            for (std::size_t phase = 0;
                 phase < problem_.phase_count;
                 ++phase) {
                const std::size_t offset = phase * problem_.block_size;
                for (std::size_t coordinate = 0;
                     coordinate < problem_.dimension;
                     ++coordinate) {
                    const std::size_t variable = offset + 1 + coordinate;
                    evaluation.lagrangian_hessian[
                        variable * static_cast<std::size_t>(n) + offset
                    ] += lambda[coordinate + 1];
                    evaluation.lagrangian_hessian[
                        offset * static_cast<std::size_t>(n) + variable
                    ] += lambda[coordinate + 1];
                }
            }
            for (Ipopt::Index row = 0; row < n; ++row) {
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    values[position++] = evaluation.lagrangian_hessian[
                        static_cast<std::size_t>(row * n + column)
                    ];
                }
            }
            return true;
        } catch (...) {
            return false;
        }
    }

    bool intermediate_callback(
        Ipopt::AlgorithmMode,
        Ipopt::Index iteration,
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
        iterations_ = std::max(iterations_, static_cast<int>(iteration) + 1);
        return true;
    }

    void finalize_solution(
        Ipopt::SolverReturn status,
        Ipopt::Index n,
        const Ipopt::Number* x,
        const Ipopt::Number* z_lower,
        const Ipopt::Number* z_upper,
        Ipopt::Index m,
        const Ipopt::Number*,
        const Ipopt::Number* lambda,
        Ipopt::Number objective,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*
    ) override {
        converged_ = status == Ipopt::SUCCESS
            || status == Ipopt::STOP_AT_ACCEPTABLE_POINT;
        objective_ = objective;
        variables_.assign(x, x + n);
        lower_multipliers_.assign(z_lower, z_lower + n);
        upper_multipliers_.assign(z_upper, z_upper + n);
        constraint_multipliers_.assign(lambda, lambda + m);
    }

    bool converged() const { return converged_; }
    int iterations() const { return iterations_; }
    double objective() const { return objective_; }
    const std::vector<double>& variables() const { return variables_; }
    const std::vector<double>& lower_multipliers() const {
        return lower_multipliers_;
    }
    const std::vector<double>& upper_multipliers() const {
        return upper_multipliers_;
    }
    const std::vector<double>& constraint_multipliers() const {
        return constraint_multipliers_;
    }

private:
    bool domain_valid(Ipopt::Index n, const Ipopt::Number* x) const {
        if (n != static_cast<Ipopt::Index>(initial_.size())) {
            return false;
        }
        for (std::size_t phase = 0; phase < problem_.phase_count; ++phase) {
            const std::size_t offset = phase * problem_.block_size;
            for (const Held2PolytopeConstraint& constraint :
                 problem_.coordinates.polytope_constraints) {
                double value = 0.0;
                for (std::size_t coordinate = 0;
                     coordinate < problem_.dimension;
                     ++coordinate) {
                    value += constraint.coefficients[coordinate]
                        * x[offset + 1 + coordinate];
                }
                if (value > constraint.upper_bound) {
                    return false;
                }
            }
        }
        return true;
    }

    void evaluate(Ipopt::Index n, const Ipopt::Number* x) {
        const std::vector<double> variables(x, x + n);
        if (variables == cached_variables_) {
            return;
        }
        cached_variables_ = variables;
        cached_evaluation_ = evaluate_problem67(
            problem_.coordinates,
            problem_.feed,
            problem_.evaluate,
            problem_.phase_count,
            cached_variables_,
            std::vector<double>(
                problem_.dimension + 1 + problem_.phase_count, 0.0
            )
        );
    }

    Problem67& problem_;
    std::vector<double> lower_;
    std::vector<double> upper_;
    std::vector<double> initial_;
    std::vector<Problem67Constraint> constraints_;
    std::vector<double> cached_variables_;
    Held2Problem67Evaluation cached_evaluation_;
    bool converged_ = false;
    int iterations_ = 0;
    double objective_ = 0.0;
    std::vector<double> variables_;
    std::vector<double> lower_multipliers_;
    std::vector<double> upper_multipliers_;
    std::vector<double> constraint_multipliers_;
};

std::string ipopt_status(Ipopt::ApplicationReturnStatus status) {
    switch (status) {
        case Ipopt::Solve_Succeeded: return "solve_succeeded";
        case Ipopt::Solved_To_Acceptable_Level:
            return "solved_to_acceptable_level";
        case Ipopt::Infeasible_Problem_Detected:
            return "infeasible_problem_detected";
        case Ipopt::Maximum_Iterations_Exceeded:
            return "maximum_iterations_exceeded";
        default:
            return "ipopt_status_" + std::to_string(
                static_cast<int>(status)
            );
    }
}

std::string problem67_feasibility(
    const Held2Coordinates& coordinates,
    const std::vector<double>& feed,
    const std::vector<Held2StageIICandidate>& candidates
) {
    const std::size_t phase_count = candidates.size();
    const std::size_t dimension = feed.size();
    const std::size_t fraction_count = phase_count;
    HighsModel model;
    model.lp_.num_col_ = static_cast<HighsInt>(
        fraction_count + phase_count * dimension
    );
    model.lp_.num_row_ = static_cast<HighsInt>(
        1 + dimension + 2 * phase_count * dimension
        + phase_count * coordinates.polytope_constraints.size()
    );
    model.lp_.col_cost_.assign(
        static_cast<std::size_t>(model.lp_.num_col_), 0.0
    );
    model.lp_.col_lower_.assign(
        static_cast<std::size_t>(model.lp_.num_col_), 0.0
    );
    model.lp_.col_upper_.assign(
        static_cast<std::size_t>(model.lp_.num_col_), 1.0
    );
    model.lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
    model.lp_.a_matrix_.start_ = {0};
    const auto row = [&model](double lower, double upper) {
        model.lp_.row_lower_.push_back(lower);
        model.lp_.row_upper_.push_back(upper);
        model.lp_.a_matrix_.start_.push_back(
            static_cast<HighsInt>(model.lp_.a_matrix_.index_.size())
        );
    };
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        model.lp_.a_matrix_.index_.push_back(
            static_cast<HighsInt>(phase)
        );
        model.lp_.a_matrix_.value_.push_back(1.0);
    }
    row(1.0, 1.0);
    for (std::size_t coordinate = 0;
         coordinate < dimension;
         ++coordinate) {
        for (std::size_t phase = 0; phase < phase_count; ++phase) {
            model.lp_.a_matrix_.index_.push_back(static_cast<HighsInt>(
                fraction_count + phase * dimension + coordinate
            ));
            model.lp_.a_matrix_.value_.push_back(1.0);
        }
        row(feed[coordinate], feed[coordinate]);
    }
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        for (std::size_t coordinate = 0;
             coordinate < dimension;
             ++coordinate) {
            const double lower = std::max(
                coordinates.independent_lower_bounds[coordinate],
                candidates[phase]
                    .independent_modified_fractions[coordinate]
                    - kHeld2Problem67Radius
            );
            const double upper = std::min(
                coordinates.independent_upper_bounds[coordinate],
                candidates[phase]
                    .independent_modified_fractions[coordinate]
                    + kHeld2Problem67Radius
            );
            const HighsInt fraction = static_cast<HighsInt>(phase);
            const HighsInt weighted = static_cast<HighsInt>(
                fraction_count + phase * dimension + coordinate
            );
            model.lp_.a_matrix_.index_.push_back(fraction);
            model.lp_.a_matrix_.value_.push_back(-lower);
            model.lp_.a_matrix_.index_.push_back(weighted);
            model.lp_.a_matrix_.value_.push_back(1.0);
            row(0.0, kHighsInf);
            model.lp_.a_matrix_.index_.push_back(fraction);
            model.lp_.a_matrix_.value_.push_back(-upper);
            model.lp_.a_matrix_.index_.push_back(weighted);
            model.lp_.a_matrix_.value_.push_back(1.0);
            row(-kHighsInf, 0.0);
        }
    }
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        for (const Held2PolytopeConstraint& constraint :
             coordinates.polytope_constraints) {
            model.lp_.a_matrix_.index_.push_back(
                static_cast<HighsInt>(phase)
            );
            model.lp_.a_matrix_.value_.push_back(
                -constraint.upper_bound
            );
            for (std::size_t coordinate = 0;
                 coordinate < dimension;
                 ++coordinate) {
                model.lp_.a_matrix_.index_.push_back(
                    static_cast<HighsInt>(
                        fraction_count + phase * dimension + coordinate
                    )
                );
                model.lp_.a_matrix_.value_.push_back(
                    constraint.coefficients[coordinate]
                );
            }
            row(-kHighsInf, 0.0);
        }
    }
    Highs highs;
    if (highs.setOptionValue("output_flag", false) == HighsStatus::kError
        || highs.setOptionValue("threads", 1) == HighsStatus::kError
        || highs.setOptionValue(
            "primal_feasibility_tolerance", kHeld2IpoptConstraint.atol
        ) == HighsStatus::kError
        || highs.passModel(model) == HighsStatus::kError
        || highs.run() == HighsStatus::kError) {
        return "indeterminate";
    }
    if (highs.getModelStatus() == HighsModelStatus::kInfeasible) {
        return "infeasible";
    }
    return highs.getModelStatus() == HighsModelStatus::kOptimal
        && highs.getInfo().primal_solution_status
            == kSolutionStatusFeasible
        ? "feasible" : "indeterminate";
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

bool recover_phase_fractions(
    const std::vector<double>& feed,
    std::vector<Held2Problem67Phase>& phases
) {
    const HighsInt count = static_cast<HighsInt>(phases.size());
    HighsModel model;
    model.lp_.num_col_ = count;
    model.lp_.num_row_ = static_cast<HighsInt>(feed.size() + 1);
    model.lp_.col_lower_.assign(phases.size(), 0.0);
    model.lp_.col_upper_.assign(phases.size(), 1.0);
    model.lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
    model.lp_.a_matrix_.start_ = {0};
    model.hessian_.dim_ = count;
    model.hessian_.start_ = {0};
    for (HighsInt phase = 0; phase < count; ++phase) {
        const double fraction =
            phases[static_cast<std::size_t>(phase)].phase_fraction;
        model.lp_.col_cost_.push_back(-fraction);
        model.hessian_.index_.push_back(phase);
        model.hessian_.value_.push_back(1.0);
        model.hessian_.start_.push_back(phase + 1);
        model.lp_.a_matrix_.index_.push_back(phase);
        model.lp_.a_matrix_.value_.push_back(1.0);
    }
    model.lp_.row_lower_.push_back(1.0);
    model.lp_.row_upper_.push_back(1.0);
    model.lp_.a_matrix_.start_.push_back(count);
    for (std::size_t coordinate = 0; coordinate < feed.size(); ++coordinate) {
        for (HighsInt phase = 0; phase < count; ++phase) {
            model.lp_.a_matrix_.index_.push_back(phase);
            model.lp_.a_matrix_.value_.push_back(
                phases[static_cast<std::size_t>(phase)]
                    .physical_fractions[coordinate]
            );
        }
        model.lp_.row_lower_.push_back(feed[coordinate]);
        model.lp_.row_upper_.push_back(feed[coordinate]);
        model.lp_.a_matrix_.start_.push_back(
            static_cast<HighsInt>(model.lp_.a_matrix_.index_.size())
        );
    }
    Highs highs;
    if (highs.setOptionValue("output_flag", false) == HighsStatus::kError
        || highs.setOptionValue("threads", 1) == HighsStatus::kError
        || highs.passModel(model) == HighsStatus::kError
        || highs.run() == HighsStatus::kError
        || highs.getModelStatus() != HighsModelStatus::kOptimal) {
        return false;
    }
    const HighsSolution& solution = highs.getSolution();
    if (!solution.value_valid
        || solution.col_value.size() != phases.size()) {
        return false;
    }
    for (std::size_t phase = 0; phase < phases.size(); ++phase) {
        phases[phase].phase_fraction = solution.col_value[phase];
    }
    return true;
}

Held2Problem67Evaluation evaluate_problem67(
    const Held2Coordinates& coordinates,
    const std::vector<double>& feed,
    const Held2StateEvaluator& evaluator,
    std::size_t phase_count,
    const std::vector<double>& variables,
    const std::vector<double>& multipliers
) {
    const std::size_t dimension = feed.size();
    const std::size_t block_size = dimension + 2;
    const std::size_t variable_count = phase_count * block_size;
    const std::size_t constraint_count = dimension + 1 + phase_count;
    if (variables.size() != variable_count
        || multipliers.size() != constraint_count) {
        throw std::invalid_argument(
            "HELD2 Problem (67) dimensions changed"
        );
    }
    Held2Problem67Evaluation result;
    result.objective_gradient.assign(variable_count, 0.0);
    result.constraints.assign(constraint_count, 0.0);
    result.constraint_jacobian.assign(
        constraint_count * variable_count, 0.0
    );
    result.lagrangian_hessian.assign(
        variable_count * variable_count, 0.0
    );
    result.constraints[0] = -1.0;
    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
        result.constraints[coordinate + 1] = -feed[coordinate];
    }
    const double sum_upper = dependent_upper(coordinates);
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const std::size_t offset = phase * block_size;
        const double fraction = variables[offset];
        const std::vector<double> composition(
            variables.begin() + static_cast<std::ptrdiff_t>(offset + 1),
            variables.begin() + static_cast<std::ptrdiff_t>(
                offset + 1 + dimension
            )
        );
        const Held2StateEvaluation state = evaluator(
            composition, variables[offset + 1 + dimension]
        );
        if (state.gradient.size() != dimension + 1
            || state.hessian.size() != (dimension + 1) * (dimension + 1)) {
            throw std::invalid_argument(
                "HELD2 Problem (67) derivative dimensions changed"
            );
        }
        result.objective += fraction * state.objective;
        result.objective_gradient[offset] = state.objective;
        result.constraints[0] += fraction;
        result.constraint_jacobian[offset] = 1.0;
        result.constraints[dimension + 1 + phase] = -sum_upper;
        for (std::size_t local = 0; local < dimension + 1; ++local) {
            const std::size_t variable = offset + 1 + local;
            result.objective_gradient[variable] =
                fraction * state.gradient[local];
            const double cross = state.gradient[local]
                + (local < dimension ? multipliers[local + 1] : 0.0);
            result.lagrangian_hessian[
                variable * variable_count + offset
            ] += cross;
            result.lagrangian_hessian[
                offset * variable_count + variable
            ] += cross;
            for (std::size_t other = 0;
                 other < dimension + 1;
                 ++other) {
                result.lagrangian_hessian[
                    variable * variable_count + offset + 1 + other
                ] += fraction
                    * state.hessian[local * (dimension + 1) + other];
            }
            if (local == dimension) {
                continue;
            }
            result.constraints[local + 1] +=
                fraction * variables[variable];
            result.constraints[dimension + 1 + phase] +=
                variables[variable];
            result.constraint_jacobian[
                (local + 1) * variable_count + offset
            ] = variables[variable];
            result.constraint_jacobian[
                (local + 1) * variable_count + variable
            ] = fraction;
            result.constraint_jacobian[
                (dimension + 1 + phase) * variable_count + variable
            ] = 1.0;
        }
    }
    result.lagrangian_gradient = result.objective_gradient;
    for (std::size_t variable = 0;
         variable < variable_count;
         ++variable) {
        for (std::size_t constraint = 0;
             constraint < constraint_count;
             ++constraint) {
            result.lagrangian_gradient[variable] +=
                multipliers[constraint]
                * result.constraint_jacobian[
                    constraint * variable_count + variable
                ];
        }
    }
    return result;
}

}  // namespace

Held2Problem67Result solve_held2_problem67(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_feed,
    const std::vector<Held2StageIICandidate>& candidates,
    const Held2StateEvaluator& evaluator,
    const std::vector<std::array<double, 2>>& phase_coordinate_bounds,
    double free_energy_upper_bound,
    const std::string& free_energy_gap_provenance,
    std::vector<double> variables,
    const Held2StateValueEvaluator& value_evaluator
) {
    Held2Problem67Result result;
    result.input_candidate_count = static_cast<int>(candidates.size());
    result.stage_iii_solve_count = 1;
    const std::size_t dimension = coordinates.independent_indices.size();
    const std::size_t block_size = dimension + 2;
    if (candidates.size() < 2
        || phase_coordinate_bounds.size() != candidates.size()) {
        result.failure_reason = "candidate_set_incomplete";
        return result;
    }

    const std::vector<double> modified_feed =
        held2_transform_physical_fractions(coordinates, physical_feed);
    std::vector<double> feed;
    for (std::size_t position : independent_positions(coordinates)) {
        feed.push_back(modified_feed[position]);
    }
    if (feed.size() != dimension) {
        result.failure_reason = "feed_dimension_changed";
        return result;
    }
    const std::string feasibility = problem67_feasibility(
        coordinates, feed, candidates
    );
    if (feasibility == "infeasible") {
        result.solver_status = "infeasible_problem_detected";
        result.failure_reason = "problem_67_infeasible";
        return result;
    }
    if (feasibility != "feasible") {
        result.solver_status = "feasibility_indeterminate";
        result.failure_reason = "problem_67_feasibility_indeterminate";
        return result;
    }
    const std::size_t variable_count = candidates.size() * block_size;
    if (!variables.empty() && variables.size() != variable_count) {
        result.failure_reason = "initial_dimension_changed";
        return result;
    }
    const bool supplied_initial = !variables.empty();
    if (!supplied_initial) {
        variables.resize(variable_count);
    }
    std::vector<double> lower(variable_count);
    std::vector<double> upper(variable_count);
    for (std::size_t phase = 0; phase < candidates.size(); ++phase) {
        if (candidates[phase].independent_modified_fractions.size()
            != dimension) {
            result.failure_reason = "candidate_dimension_changed";
            return result;
        }
        const std::size_t offset = phase * block_size;
        lower[offset] = 0.0;
        upper[offset] = 1.0;
        if (!supplied_initial) {
            variables[offset] =
                1.0 / static_cast<double>(candidates.size());
        }
        for (std::size_t coordinate = 0;
             coordinate < dimension;
             ++coordinate) {
            const std::size_t variable = offset + 1 + coordinate;
            lower[variable] = std::max(
                coordinates.independent_lower_bounds[coordinate],
                candidates[phase]
                    .independent_modified_fractions[coordinate]
                    - kHeld2Problem67Radius
            );
            upper[variable] = std::min(
                coordinates.independent_upper_bounds[coordinate],
                candidates[phase]
                    .independent_modified_fractions[coordinate]
                    + kHeld2Problem67Radius
            );
            if (!supplied_initial) {
                variables[variable] = candidates[phase]
                    .independent_modified_fractions[coordinate];
            }
        }
        const std::size_t volume = offset + 1 + dimension;
        lower[volume] = phase_coordinate_bounds[phase][0];
        upper[volume] = phase_coordinate_bounds[phase][1];
        if (!supplied_initial) {
            variables[volume] = candidates[phase].phase_coordinate;
        }
    }

    Problem67 problem{
        coordinates,
        feed,
        evaluator,
        value_evaluator,
        candidates.size(),
        dimension,
        block_size,
        0,
        {},
    };
    std::vector<Problem67Constraint> constraints;
    constraints.reserve(
        dimension + 1
        + candidates.size() * coordinates.polytope_constraints.size()
    );
    constraints.push_back({&problem, 0, 0});
    for (std::size_t coordinate = 0;
         coordinate < dimension;
         ++coordinate) {
        constraints.push_back({&problem, 1, coordinate});
    }
    for (std::size_t phase = 0; phase < candidates.size(); ++phase) {
        for (std::size_t constraint = 0;
             constraint < coordinates.polytope_constraints.size();
             ++constraint) {
            constraints.push_back({&problem, 2, constraint, phase});
        }
    }

    auto* raw = new Problem67Tnlp(
        problem, lower, upper, variables
    );
    Ipopt::SmartPtr<Ipopt::TNLP> nlp_problem = raw;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application =
        IpoptApplicationFactory();
    application->Options()->SetStringValue("option_file_name", "");
    application->Options()->SetIntegerValue("print_level", 0);
    application->Options()->SetStringValue("sb", "yes");
    application->Options()->SetIntegerValue("max_iter", 300);
    application->Options()->SetNumericValue(
        "tol", kHeld2IpoptTarget.atol
    );
    application->Options()->SetNumericValue(
        "constr_viol_tol", kHeld2IpoptConstraint.atol
    );
    application->Options()->SetStringValue(
        "jacobian_approximation", "exact"
    );
    application->Options()->SetStringValue(
        "hessian_approximation", "exact"
    );
    application->Options()->SetStringValue("mu_strategy", "adaptive");
    application->Options()->SetNumericValue("bound_relax_factor", 0.0);
    if (application->Initialize() != Ipopt::Solve_Succeeded) {
        result.solver_status = "initialization_failed";
        result.failure_reason = "stage_iii_solver_initialization_failed";
        return result;
    }
    const Ipopt::ApplicationReturnStatus status =
        application->OptimizeTNLP(nlp_problem);
    result.solver_status = ipopt_status(status);
    result.optimizer_iteration_count = raw->iterations();
    if (!raw->converged()
        && raw->variables().size() != variable_count) {
        result.numerical_status = "not_converged";
        result.failure_reason = "stage_iii_solver_not_converged";
        return result;
    }
    variables = raw->variables();
    const double objective = raw->objective();
    const std::vector<double>& z_lower = raw->lower_multipliers();
    const std::vector<double>& z_upper = raw->upper_multipliers();
    const std::vector<double>& lambda = raw->constraint_multipliers();
    for (std::size_t phase = 0; phase < candidates.size(); ++phase) {
        const std::size_t offset = phase * block_size;
        if (z_lower[offset] > kHeld2PhaseRetirementMargin.atol) {
            continue;
        }
        const std::vector<double> independent(
            variables.begin() + static_cast<std::ptrdiff_t>(offset + 1),
            variables.begin() + static_cast<std::ptrdiff_t>(
                offset + 1 + dimension
            )
        );
        double& log_volume = variables[offset + 1 + dimension];
        for (int iteration = 0; iteration < 8; ++iteration) {
            const Held2StateEvaluation state = evaluator(
                independent, log_volume
            );
            if (audit_held2_tolerance(
                    kHeld2Stage3Pressure,
                    state.pressure_stationarity_relative
                ).passed) {
                break;
            }
            const double derivative =
                state.pressure_stationarity_derivative_log_volume;
            if (!std::isfinite(derivative) || derivative == 0.0) {
                result.failure_reason =
                    "stage_iii_pressure_refinement_failed";
                return result;
            }
            const double next = std::clamp(
                log_volume
                    - state.pressure_stationarity_relative / derivative,
                phase_coordinate_bounds[phase][0],
                phase_coordinate_bounds[phase][1]
            );
            if (next == log_volume) {
                result.failure_reason =
                    "stage_iii_pressure_refinement_failed";
                return result;
            }
            log_volume = next;
            ++result.pressure_polish_iteration_count;
        }
        if (!audit_held2_tolerance(
                kHeld2Stage3Pressure,
                evaluator(independent, log_volume)
                    .pressure_stationarity_relative
            ).passed) {
            result.failure_reason =
                "stage_iii_pressure_refinement_failed";
            return result;
        }
    }
    result.pressure_polish_status = "complete";
    result.solution_variables = variables;

    std::vector<double> equality_residuals;
    equality_residuals.reserve(dimension + 1);
    std::vector<double> no_gradient;
    for (std::size_t index = 0; index < dimension + 1; ++index) {
        equality_residuals.push_back(problem67_constraint(
            variables, no_gradient, &constraints[index]
        ));
    }
    result.modified_balance_inf_norm = maximum_abs(equality_residuals);
    if (!audit_held2_tolerance(
            kHeld2Stage3ModifiedBalance,
            result.modified_balance_inf_norm
        ).passed) {
        result.numerical_status = "not_converged";
        result.failure_reason = "stage_iii_material_balance_failed";
        return result;
    }

    std::vector<double> stationarity(variable_count);
    static_cast<void>(problem67_objective(
        variables, stationarity, &problem
    ));
    std::vector<double> constraint_values(constraints.size());
    for (std::size_t row = 0; row < constraints.size(); ++row) {
        std::vector<double> gradient(variable_count);
        constraint_values[row] = problem67_constraint(
            variables, gradient, &constraints[row]
        );
        for (std::size_t column = 0; column < variable_count; ++column) {
            stationarity[column] += lambda[row] * gradient[column];
        }
    }
    for (std::size_t column = 0; column < variable_count; ++column) {
        stationarity[column] += z_upper[column] - z_lower[column];
        result.dual_sign_violation_inf_norm = std::max({
            result.dual_sign_violation_inf_norm,
            -z_lower[column],
            -z_upper[column],
        });
        result.bound_complementarity_inf_norm = std::max({
            result.bound_complementarity_inf_norm,
            std::abs((variables[column] - lower[column]) * z_lower[column]),
            std::abs((upper[column] - variables[column]) * z_upper[column]),
        });
    }
    for (std::size_t row = dimension + 1;
         row < constraints.size();
         ++row) {
        result.dual_sign_violation_inf_norm = std::max(
            result.dual_sign_violation_inf_norm, -lambda[row]
        );
        result.bound_complementarity_inf_norm = std::max(
            result.bound_complementarity_inf_norm,
            std::abs(lambda[row] * constraint_values[row])
        );
    }
    result.kkt_stationarity_inf_norm = maximum_abs(stationarity);

    result.objective = objective;
    result.free_energy_upper_bound = free_energy_upper_bound;
    result.free_energy_gap_provenance = free_energy_gap_provenance;
    result.free_energy_gap_available =
        std::isfinite(free_energy_upper_bound)
        && !free_energy_gap_provenance.empty()
        && free_energy_gap_provenance != "unavailable";
    if (result.free_energy_gap_available) {
        result.free_energy_gap = free_energy_upper_bound - objective;
    }
    result.minimum_phase_fraction =
        std::numeric_limits<double>::infinity();
    for (std::size_t phase = 0; phase < candidates.size(); ++phase) {
        const std::size_t offset = phase * block_size;
        const double fraction = variables[offset];
        if (z_lower[offset] > kHeld2PhaseRetirementMargin.atol) {
            ++result.retired_inactive_count;
            continue;
        }
        const std::vector<double> independent(
            variables.begin() + static_cast<std::ptrdiff_t>(offset + 1),
            variables.begin() + static_cast<std::ptrdiff_t>(
                offset + 1 + dimension
            )
        );
        const Held2StateEvaluation state = evaluator(
            independent, variables[offset + 1 + dimension]
        );
        Held2Problem67Phase current{
            fraction,
            state.modified_fractions,
            state.physical_amounts,
            state.volume,
        };
        auto duplicate = std::find_if(
            result.phases.begin(),
            result.phases.end(),
            [&](const Held2Problem67Phase& known) {
                return std::max(
                    maximum_difference(
                        known.physical_fractions,
                        current.physical_fractions
                    ),
                    std::abs(
                        std::log(known.volume)
                        - std::log(current.volume)
                    )
                ) <= kHeld2PhaseMerge.atol;
            }
        );
        if (duplicate == result.phases.end()) {
            result.phases.push_back(std::move(current));
        } else {
            const double combined =
                duplicate->phase_fraction + current.phase_fraction;
            const double retained_weight =
                duplicate->phase_fraction / combined;
            const double current_weight =
                current.phase_fraction / combined;
            for (std::size_t index = 0;
                 index < duplicate->modified_fractions.size();
                 ++index) {
                duplicate->modified_fractions[index] =
                    retained_weight * duplicate->modified_fractions[index]
                    + current_weight * current.modified_fractions[index];
            }
            for (std::size_t index = 0;
                 index < duplicate->physical_fractions.size();
                 ++index) {
                duplicate->physical_fractions[index] =
                    retained_weight * duplicate->physical_fractions[index]
                    + current_weight * current.physical_fractions[index];
            }
            duplicate->volume =
                retained_weight * duplicate->volume
                + current_weight * current.volume;
            duplicate->phase_fraction = combined;
            ++result.retired_duplicate_count;
        }
    }
    if (result.phases.size() < 2) {
        result.failure_reason = "collapsed_phase_set";
        return result;
    }
    if (result.retired_inactive_count > 0) {
        ++result.stage_iii_solve_count;
        if (!recover_phase_fractions(
                physical_feed, result.phases
            )) {
            result.failure_reason =
                "stage_iii_phase_fraction_recovery_failed";
            return result;
        }
        result.phases.erase(
            std::remove_if(
                result.phases.begin(),
                result.phases.end(),
            [](const Held2Problem67Phase& phase) {
                return phase.phase_fraction <= 0.0;
            }),
            result.phases.end()
        );
        if (result.phases.size() < 2) {
            result.failure_reason = "collapsed_phase_set";
            return result;
        }
    }

    std::vector<double> modified_balance(modified_feed.size());
    std::vector<double> ordinary_balance(physical_feed.size());
    result.objective = 0.0;
    result.pressure_stationarity_inf_norm = 0.0;
    result.minimum_phase_distance =
        std::numeric_limits<double>::infinity();
    for (const Held2Problem67Phase& phase : result.phases) {
        result.minimum_phase_fraction = std::min(
            result.minimum_phase_fraction, phase.phase_fraction
        );
        for (std::size_t index = 0;
             index < modified_balance.size();
             ++index) {
            modified_balance[index] +=
                phase.phase_fraction * phase.modified_fractions[index];
        }
        for (std::size_t index = 0;
             index < ordinary_balance.size();
             ++index) {
            ordinary_balance[index] +=
                phase.phase_fraction * phase.physical_fractions[index];
        }
        const double charge = charge_residual(
            coordinates.charges, phase.physical_fractions
        );
        result.phase_charge_inf_norm = std::max(
            result.phase_charge_inf_norm, std::abs(charge)
        );
        double charge_scale = 0.0;
        for (std::size_t index = 0;
             index < coordinates.charges.size();
             ++index) {
            charge_scale += std::abs(
                coordinates.charges[index]
                    * phase.physical_fractions[index]
            );
        }
        result.phase_charge_scale = std::max(
            result.phase_charge_scale, charge_scale
        );
        std::vector<double> independent;
        for (std::size_t position : independent_positions(coordinates)) {
            independent.push_back(phase.modified_fractions[position]);
        }
        const Held2StateEvaluation state = evaluator(
            independent, std::log(phase.volume)
        );
        result.objective += phase.phase_fraction * state.objective;
        result.pressure_stationarity_inf_norm = std::max(
            result.pressure_stationarity_inf_norm,
            std::abs(state.pressure_stationarity_relative)
        );
    }
    result.modified_balance_inf_norm = maximum_difference(
        modified_balance, modified_feed
    );
    result.ordinary_balance_inf_norm = maximum_difference(
        ordinary_balance, physical_feed
    );
    if (!audit_held2_tolerance(
            kHeld2Stage3ModifiedBalance,
            result.modified_balance_inf_norm
        ).passed
        || !audit_held2_tolerance(
            kHeld2Stage3ExplicitBalance,
            result.ordinary_balance_inf_norm
        ).passed) {
        result.failure_reason = "stage_iii_active_set_balance_failed";
        return result;
    }
    for (std::size_t left = 0; left < result.phases.size(); ++left) {
        for (std::size_t right = left + 1;
             right < result.phases.size();
             ++right) {
            result.minimum_phase_distance = std::min(
                result.minimum_phase_distance,
                std::max(
                    maximum_difference(
                        result.phases[left].physical_fractions,
                        result.phases[right].physical_fractions
                    ),
                    std::abs(
                        std::log(result.phases[left].volume)
                        - std::log(result.phases[right].volume)
                    )
                )
            );
        }
    }
    result.numerical_status = "converged";
    result.physical_status = "accepted";
    result.feedback = "none";
    result.phase_identity_status = "paper_duplicate_removal_complete";
    result.kkt_evidence_available = true;
    result.physical_evidence_available = true;
    result.phase_identity_evidence_available = true;
    return result;
}

}  // namespace epcsaft_equilibrium
