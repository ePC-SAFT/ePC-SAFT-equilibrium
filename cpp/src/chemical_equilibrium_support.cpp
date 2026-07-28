#include "chemical_equilibrium.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/rational.hpp>

#include "Highs.h"

namespace epcsaft_equilibrium {
namespace {

using ExactInteger = boost::multiprecision::cpp_int;
using ExactRational = boost::rational<ExactInteger>;
using ExactDecimal = boost::multiprecision::cpp_dec_float_50;

constexpr std::size_t kMaximumBasisCandidates = 100000;

struct ExactEqualitySystem {
    std::vector<std::vector<ExactRational>> matrix;
    std::vector<ExactRational> right_hand_side;
    std::size_t columns = 0;
};

struct HighsCandidate {
    bool optimal = false;
    std::string status = "setup_failed";
    std::vector<double> values;
    std::vector<double> row_duals;
    std::vector<double> column_duals;
    std::vector<bool> basic_columns;
    double objective = 0.0;
};

ExactRational exact_binary_rational(double value) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument("exact support validation requires finite inputs");
    }
    if (value == 0.0) {
        return ExactRational(0);
    }
    std::uint64_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value));
    std::memcpy(&bits, &value, sizeof(value));
    const bool negative = (bits >> 63U) != 0U;
    const std::uint64_t exponent_bits = (bits >> 52U) & 0x7ffU;
    const std::uint64_t fraction_bits = bits & ((std::uint64_t{1} << 52U) - 1U);
    ExactInteger significand;
    int binary_exponent = 0;
    if (exponent_bits == 0U) {
        significand = fraction_bits;
        binary_exponent = -1022 - 52;
    } else {
        significand = (std::uint64_t{1} << 52U) | fraction_bits;
        binary_exponent = static_cast<int>(exponent_bits) - 1023 - 52;
    }
    if (negative) {
        significand = -significand;
    }
    if (binary_exponent >= 0) {
        significand <<= binary_exponent;
        return ExactRational(significand);
    }
    ExactInteger denominator = 1;
    denominator <<= -binary_exponent;
    return ExactRational(significand, denominator);
}

double to_double(const ExactRational& value) {
    return (
        ExactDecimal(value.numerator()) / ExactDecimal(value.denominator())
    ).convert_to<double>();
}

bool solve_square_system(
    std::vector<std::vector<ExactRational>> matrix,
    std::vector<ExactRational> right_hand_side,
    std::vector<ExactRational>& solution
) {
    const std::size_t dimension = matrix.size();
    if (right_hand_side.size() != dimension
        || !std::all_of(matrix.begin(), matrix.end(), [dimension](const auto& row) {
            return row.size() == dimension;
        })) {
        return false;
    }
    for (std::size_t column = 0; column < dimension; ++column) {
        std::size_t pivot = column;
        while (pivot < dimension && matrix[pivot][column] == 0) {
            ++pivot;
        }
        if (pivot == dimension) {
            return false;
        }
        if (pivot != column) {
            std::swap(matrix[pivot], matrix[column]);
            std::swap(right_hand_side[pivot], right_hand_side[column]);
        }
        const ExactRational diagonal = matrix[column][column];
        for (std::size_t index = column; index < dimension; ++index) {
            matrix[column][index] /= diagonal;
        }
        right_hand_side[column] /= diagonal;
        for (std::size_t row = 0; row < dimension; ++row) {
            if (row == column || matrix[row][column] == 0) {
                continue;
            }
            const ExactRational factor = matrix[row][column];
            for (std::size_t index = column; index < dimension; ++index) {
                matrix[row][index] -= factor * matrix[column][index];
            }
            right_hand_side[row] -= factor * right_hand_side[column];
        }
    }
    solution = std::move(right_hand_side);
    return true;
}

ExactEqualitySystem independent_exact_equalities(
    const DenseMatrix& balance_matrix,
    const std::vector<double>& balance_totals,
    const std::vector<int>& charges
) {
    const std::size_t columns = balance_matrix.columns;
    std::vector<std::vector<ExactRational>> augmented;
    augmented.reserve(balance_matrix.rows + 1);
    for (std::size_t row = 0; row < balance_matrix.rows; ++row) {
        std::vector<ExactRational> values(columns + 1, ExactRational(0));
        for (std::size_t column = 0; column < columns; ++column) {
            values[column] = exact_binary_rational(balance_matrix(row, column));
        }
        values.back() = exact_binary_rational(balance_totals[row]);
        augmented.push_back(std::move(values));
    }
    if (std::any_of(charges.begin(), charges.end(), [](int charge) {
            return charge != 0;
        })) {
        std::vector<ExactRational> charge_row(columns + 1, ExactRational(0));
        for (std::size_t column = 0; column < columns; ++column) {
            charge_row[column] = ExactRational(charges[column]);
        }
        augmented.push_back(std::move(charge_row));
    }

    std::size_t rank = 0;
    for (std::size_t column = 0; column < columns && rank < augmented.size(); ++column) {
        std::size_t pivot = rank;
        while (pivot < augmented.size() && augmented[pivot][column] == 0) {
            ++pivot;
        }
        if (pivot == augmented.size()) {
            continue;
        }
        if (pivot != rank) {
            std::swap(augmented[pivot], augmented[rank]);
        }
        const ExactRational diagonal = augmented[rank][column];
        for (std::size_t index = column; index <= columns; ++index) {
            augmented[rank][index] /= diagonal;
        }
        for (std::size_t row = 0; row < augmented.size(); ++row) {
            if (row == rank || augmented[row][column] == 0) {
                continue;
            }
            const ExactRational factor = augmented[row][column];
            for (std::size_t index = column; index <= columns; ++index) {
                augmented[row][index] -= factor * augmented[rank][index];
            }
        }
        ++rank;
    }
    for (std::size_t row = rank; row < augmented.size(); ++row) {
        const bool zero_coefficients = std::all_of(
            augmented[row].begin(),
            augmented[row].begin() + static_cast<std::ptrdiff_t>(columns),
            [](const ExactRational& value) { return value == 0; }
        );
        if (zero_coefficients && augmented[row].back() != 0) {
            throw std::invalid_argument("homogeneous support equalities are inconsistent");
        }
    }
    if (rank == 0) {
        throw std::invalid_argument("homogeneous support requires a conserved total");
    }

    ExactEqualitySystem result;
    result.columns = columns;
    result.matrix.assign(rank, std::vector<ExactRational>(columns, ExactRational(0)));
    result.right_hand_side.resize(rank, ExactRational(0));
    for (std::size_t row = 0; row < rank; ++row) {
        std::copy_n(
            augmented[row].begin(),
            static_cast<std::ptrdiff_t>(columns),
            result.matrix[row].begin()
        );
        result.right_hand_side[row] = augmented[row].back();
    }
    return result;
}

std::vector<std::size_t> ordered_columns(const HighsCandidate& candidate) {
    std::vector<std::size_t> result(candidate.values.size(), 0);
    std::iota(result.begin(), result.end(), 0);
    std::stable_sort(result.begin(), result.end(), [&candidate](
        std::size_t left, std::size_t right
    ) {
        const bool left_basic = left < candidate.basic_columns.size()
            && candidate.basic_columns[left];
        const bool right_basic = right < candidate.basic_columns.size()
            && candidate.basic_columns[right];
        if (left_basic != right_basic) {
            return left_basic;
        }
        const bool left_positive = candidate.values[left] > 0.0;
        const bool right_positive = candidate.values[right] > 0.0;
        if (left_positive != right_positive) {
            return left_positive;
        }
        if (candidate.values[left] != candidate.values[right]) {
            return candidate.values[left] > candidate.values[right];
        }
        return left < right;
    });
    return result;
}

bool visit_combinations(
    const std::vector<std::size_t>& pool,
    std::size_t choose,
    std::size_t offset,
    std::vector<std::size_t>& selected,
    std::size_t& visited,
    const std::function<bool(const std::vector<std::size_t>&)>& visitor
) {
    if (choose == 0) {
        if (++visited > kMaximumBasisCandidates) {
            return true;
        }
        return visitor(selected);
    }
    if (pool.size() - offset < choose) {
        return false;
    }
    for (std::size_t index = offset; index + choose <= pool.size(); ++index) {
        selected.push_back(pool[index]);
        if (visit_combinations(
                pool, choose - 1, index + 1, selected, visited, visitor
            )) {
            return true;
        }
        selected.pop_back();
    }
    return false;
}

std::optional<std::vector<ExactRational>> exact_primal_witness(
    const ExactEqualitySystem& system,
    const HighsCandidate& candidate,
    std::optional<std::size_t> required_positive
) {
    const std::size_t rank = system.matrix.size();
    if (rank > system.columns || candidate.values.size() != system.columns) {
        return std::nullopt;
    }
    std::vector<std::size_t> pool = ordered_columns(candidate);
    std::vector<std::size_t> selected;
    if (required_positive.has_value()) {
        selected.push_back(*required_positive);
        pool.erase(std::remove(pool.begin(), pool.end(), *required_positive), pool.end());
    }
    const std::size_t required_count = selected.size();
    if (required_count > rank) {
        return std::nullopt;
    }
    std::optional<std::vector<ExactRational>> result;
    std::size_t visited = 0;
    visit_combinations(
        pool,
        rank - required_count,
        0,
        selected,
        visited,
        [&](const std::vector<std::size_t>& basis_columns) {
            std::vector<std::vector<ExactRational>> basis(
                rank, std::vector<ExactRational>(rank, ExactRational(0))
            );
            for (std::size_t row = 0; row < rank; ++row) {
                for (std::size_t column = 0; column < rank; ++column) {
                    basis[row][column] =
                        system.matrix[row][basis_columns[column]];
                }
            }
            std::vector<ExactRational> basis_values;
            if (!solve_square_system(
                    std::move(basis), system.right_hand_side, basis_values
                )
                || std::any_of(
                    basis_values.begin(), basis_values.end(),
                    [](const ExactRational& value) { return value < 0; }
                )) {
                return false;
            }
            std::vector<ExactRational> witness(
                system.columns, ExactRational(0)
            );
            for (std::size_t index = 0; index < rank; ++index) {
                witness[basis_columns[index]] = basis_values[index];
            }
            if (required_positive.has_value()
                && witness[*required_positive] <= 0) {
                return false;
            }
            for (std::size_t row = 0; row < rank; ++row) {
                ExactRational residual = -system.right_hand_side[row];
                for (std::size_t column = 0; column < system.columns; ++column) {
                    residual += system.matrix[row][column] * witness[column];
                }
                if (residual != 0) {
                    return false;
                }
            }
            result = std::move(witness);
            return true;
        }
    );
    return result;
}

std::optional<std::vector<ExactRational>> exact_dual_zero_certificate(
    const ExactEqualitySystem& system,
    const std::vector<ExactRational>& objective,
    const HighsCandidate& candidate
) {
    const std::size_t rank = system.matrix.size();
    if (rank > system.columns || objective.size() != system.columns) {
        return std::nullopt;
    }
    const std::vector<std::size_t> pool = ordered_columns(candidate);
    std::vector<std::size_t> selected;
    std::optional<std::vector<ExactRational>> result;
    std::size_t visited = 0;
    visit_combinations(
        pool,
        rank,
        0,
        selected,
        visited,
        [&](const std::vector<std::size_t>& basis_columns) {
            std::vector<std::vector<ExactRational>> transpose_basis(
                rank, std::vector<ExactRational>(rank, ExactRational(0))
            );
            std::vector<ExactRational> active_objective(rank, ExactRational(0));
            for (std::size_t equation = 0; equation < rank; ++equation) {
                const std::size_t species = basis_columns[equation];
                active_objective[equation] = objective[species];
                for (std::size_t multiplier = 0; multiplier < rank; ++multiplier) {
                    transpose_basis[equation][multiplier] =
                        system.matrix[multiplier][species];
                }
            }
            std::vector<ExactRational> multipliers;
            if (!solve_square_system(
                    std::move(transpose_basis),
                    std::move(active_objective),
                    multipliers
                )) {
                return false;
            }
            for (std::size_t species = 0; species < system.columns; ++species) {
                ExactRational reduced = -objective[species];
                for (std::size_t row = 0; row < rank; ++row) {
                    reduced += system.matrix[row][species] * multipliers[row];
                }
                if (reduced < 0) {
                    return false;
                }
            }
            ExactRational dual_objective(0);
            for (std::size_t row = 0; row < rank; ++row) {
                dual_objective += system.right_hand_side[row] * multipliers[row];
            }
            if (dual_objective != 0) {
                return false;
            }
            result = std::move(multipliers);
            return true;
        }
    );
    return result;
}

std::optional<std::vector<ExactRational>> exact_row_coordinates(
    const ExactEqualitySystem& system,
    const std::vector<ExactRational>& row
) {
    const std::size_t rank = system.matrix.size();
    if (row.size() != system.columns || rank > system.columns) {
        return std::nullopt;
    }
    HighsCandidate ordering;
    ordering.values.assign(system.columns, 0.0);
    ordering.basic_columns.assign(system.columns, false);
    const std::vector<std::size_t> pool = ordered_columns(ordering);
    std::vector<std::size_t> selected;
    std::optional<std::vector<ExactRational>> result;
    std::size_t visited = 0;
    visit_combinations(
        pool,
        rank,
        0,
        selected,
        visited,
        [&](const std::vector<std::size_t>& basis_columns) {
            std::vector<std::vector<ExactRational>> transpose_basis(
                rank, std::vector<ExactRational>(rank, ExactRational(0))
            );
            std::vector<ExactRational> active_values(rank, ExactRational(0));
            for (std::size_t equation = 0; equation < rank; ++equation) {
                active_values[equation] = row[basis_columns[equation]];
                for (std::size_t multiplier = 0; multiplier < rank; ++multiplier) {
                    transpose_basis[equation][multiplier] =
                        system.matrix[multiplier][basis_columns[equation]];
                }
            }
            std::vector<ExactRational> coordinates;
            if (!solve_square_system(
                    std::move(transpose_basis),
                    std::move(active_values),
                    coordinates
                )) {
                return false;
            }
            for (std::size_t column = 0; column < system.columns; ++column) {
                ExactRational residual = -row[column];
                for (std::size_t multiplier = 0; multiplier < rank; ++multiplier) {
                    residual +=
                        coordinates[multiplier] * system.matrix[multiplier][column];
                }
                if (residual != 0) {
                    return false;
                }
            }
            result = std::move(coordinates);
            return true;
        }
    );
    return result;
}

HighsCandidate solve_candidate_lp(
    const ExactEqualitySystem& system,
    const std::vector<double>& objective
) {
    HighsCandidate result;
    if (objective.size() != system.columns) {
        return result;
    }
    HighsModel model;
    model.lp_.num_col_ = static_cast<HighsInt>(system.columns);
    model.lp_.num_row_ = static_cast<HighsInt>(system.matrix.size());
    model.lp_.sense_ = ObjSense::kMaximize;
    model.lp_.col_cost_ = objective;
    model.lp_.col_lower_.assign(system.columns, 0.0);
    model.lp_.col_upper_.assign(system.columns, kHighsInf);
    model.lp_.row_lower_.reserve(system.matrix.size());
    model.lp_.row_upper_.reserve(system.matrix.size());
    for (const ExactRational& value : system.right_hand_side) {
        const double converted = to_double(value);
        model.lp_.row_lower_.push_back(converted);
        model.lp_.row_upper_.push_back(converted);
    }
    model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
    model.lp_.a_matrix_.start_ = {0};
    for (std::size_t column = 0; column < system.columns; ++column) {
        for (std::size_t row = 0; row < system.matrix.size(); ++row) {
            if (system.matrix[row][column] == 0) {
                continue;
            }
            model.lp_.a_matrix_.index_.push_back(static_cast<HighsInt>(row));
            model.lp_.a_matrix_.value_.push_back(
                to_double(system.matrix[row][column])
            );
        }
        model.lp_.a_matrix_.start_.push_back(
            static_cast<HighsInt>(model.lp_.a_matrix_.index_.size())
        );
    }
    Highs highs;
    if (highs.setOptionValue("output_flag", false) == HighsStatus::kError
        || highs.setOptionValue("threads", 1) == HighsStatus::kError
        || highs.setOptionValue("solver", std::string("simplex"))
            == HighsStatus::kError
        || highs.passModel(model) == HighsStatus::kError
        || highs.run() == HighsStatus::kError) {
        return result;
    }
    const HighsModelStatus model_status = highs.getModelStatus();
    result.status = highs.modelStatusToString(model_status);
    if (model_status != HighsModelStatus::kOptimal
        || highs.getInfo().primal_solution_status != kSolutionStatusFeasible) {
        return result;
    }
    const HighsSolution& solution = highs.getSolution();
    if (solution.col_value.size() != system.columns
        || solution.row_dual.size() != system.matrix.size()
        || solution.col_dual.size() != system.columns) {
        result.status = "optimal_without_complete_solution";
        return result;
    }
    result.optimal = true;
    result.status = "optimal";
    result.values = solution.col_value;
    result.row_duals = solution.row_dual;
    result.column_duals = solution.col_dual;
    result.basic_columns.assign(system.columns, false);
    const HighsBasis& basis = highs.getBasis();
    if (basis.valid && basis.col_status.size() == system.columns) {
        for (std::size_t column = 0; column < system.columns; ++column) {
            result.basic_columns[column] =
                basis.col_status[column] == HighsBasisStatus::kBasic;
        }
    }
    result.objective = std::inner_product(
        objective.begin(), objective.end(), result.values.begin(), 0.0
    );
    if (std::abs(result.objective) <= 1.0e-12) {
        result.objective = 0.0;
    }
    return result;
}

std::vector<double> to_doubles(const std::vector<ExactRational>& values) {
    std::vector<double> result;
    result.reserve(values.size());
    std::transform(values.begin(), values.end(), std::back_inserter(result), to_double);
    return result;
}

}  // namespace

HomogeneousSupportAnalysis analyze_homogeneous_support(
    const DenseMatrix& balance_matrix,
    const std::vector<double>& balance_totals,
    const std::vector<int>& charges,
    const std::vector<double>& molar_masses_kg_per_mol
) {
    const std::size_t species_count = balance_matrix.columns;
    if (species_count == 0
        || balance_matrix.values.size() != balance_matrix.rows * species_count
        || balance_totals.size() != balance_matrix.rows
        || charges.size() != species_count
        || molar_masses_kg_per_mol.size() != species_count
        || !std::all_of(
            balance_matrix.values.begin(),
            balance_matrix.values.end(),
            [](double value) { return std::isfinite(value); }
        )
        || !std::all_of(
            balance_totals.begin(),
            balance_totals.end(),
            [](double value) { return std::isfinite(value); }
        )
        || !std::all_of(
            molar_masses_kg_per_mol.begin(),
            molar_masses_kg_per_mol.end(),
            [](double value) { return std::isfinite(value) && value > 0.0; }
        )) {
        throw std::invalid_argument("homogeneous support input is invalid");
    }

    const ExactEqualitySystem exact = independent_exact_equalities(
        balance_matrix, balance_totals, charges
    );
    std::vector<ExactRational> exact_masses;
    exact_masses.reserve(species_count);
    std::transform(
        molar_masses_kg_per_mol.begin(),
        molar_masses_kg_per_mol.end(),
        std::back_inserter(exact_masses),
        exact_binary_rational
    );
    const auto mass_coordinates = exact_row_coordinates(exact, exact_masses);
    if (!mass_coordinates.has_value()) {
        throw std::invalid_argument(
            "molar mass must be in the span of homogeneous support balances"
        );
    }
    ExactRational total_mass(0);
    for (std::size_t row = 0; row < exact.matrix.size(); ++row) {
        total_mass += (*mass_coordinates)[row] * exact.right_hand_side[row];
    }
    if (total_mass <= 0) {
        throw std::invalid_argument("homogeneous support total mass must be positive");
    }

    HomogeneousSupportAnalysis result;
    result.species.resize(species_count);
    const HighsCandidate phase1 = solve_candidate_lp(
        exact, std::vector<double>(species_count, 0.0)
    );
    result.phase1_status = phase1.status;
    if (!phase1.optimal) {
        result.validation_status = "phase1_unresolved";
        return result;
    }
    const auto exact_phase1 = exact_primal_witness(exact, phase1, std::nullopt);
    if (!exact_phase1.has_value()) {
        result.validation_status = "exact_phase1_unresolved";
        return result;
    }

    std::vector<std::vector<ExactRational>> accessible_witnesses;
    accessible_witnesses.reserve(species_count);
    bool complete = true;
    for (std::size_t species = 0; species < species_count; ++species) {
        std::vector<ExactRational> exact_objective(
            species_count, ExactRational(0)
        );
        exact_objective[species] = exact_masses[species] / total_mass;
        std::vector<double> objective(species_count, 0.0);
        objective[species] = to_double(exact_objective[species]);
        const HighsCandidate candidate = solve_candidate_lp(exact, objective);
        SpeciesSupportEvidence& evidence = result.species[species];
        evidence.candidate_maximum_mass_fraction = candidate.objective;
        if (!candidate.optimal) {
            complete = false;
            continue;
        }
        const auto primal = exact_primal_witness(exact, candidate, species);
        if (primal.has_value()) {
            evidence.classification = "proved_accessible";
            evidence.primal_validated = true;
            evidence.witness_amounts = to_doubles(*primal);
            accessible_witnesses.push_back(*primal);
            continue;
        }
        const auto dual = exact_dual_zero_certificate(
            exact, exact_objective, candidate
        );
        if (dual.has_value()) {
            evidence.classification = "proved_structural_zero";
            evidence.dual_validated = true;
            evidence.dual_multipliers = to_doubles(*dual);
            continue;
        }
        complete = false;
    }

    if (!accessible_witnesses.empty()) {
        std::vector<ExactRational> average(species_count, ExactRational(0));
        for (const auto& witness : accessible_witnesses) {
            for (std::size_t species = 0; species < species_count; ++species) {
                average[species] += witness[species];
            }
        }
        const ExactRational divisor(accessible_witnesses.size());
        for (ExactRational& value : average) {
            value /= divisor;
        }
        for (std::size_t species = 0; species < species_count; ++species) {
            if (result.species[species].classification == "proved_accessible"
                && average[species] <= 0) {
                complete = false;
            }
        }
        result.witness_average_amounts = to_doubles(average);
    }
    result.validation_status = complete
        ? "exact_certificates_complete"
        : "exact_certificates_incomplete";
    if (complete) {
        result.equality_inf_norm = 0.0;
    } else if (!result.witness_average_amounts.empty()) {
        for (std::size_t row = 0; row < balance_matrix.rows; ++row) {
            double residual = -balance_totals[row];
            for (std::size_t species = 0; species < species_count; ++species) {
                residual += balance_matrix(row, species)
                    * result.witness_average_amounts[species];
            }
            result.equality_inf_norm = std::max(
                result.equality_inf_norm, std::abs(residual)
            );
        }
        double charge_residual = 0.0;
        for (std::size_t species = 0; species < species_count; ++species) {
            charge_residual += static_cast<double>(charges[species])
                * result.witness_average_amounts[species];
        }
        result.equality_inf_norm = std::max(
            result.equality_inf_norm, std::abs(charge_residual)
        );
    } else {
        result.equality_inf_norm = std::numeric_limits<double>::infinity();
    }
    return result;
}

}  // namespace epcsaft_equilibrium
