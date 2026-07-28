#include "chemical_equilibrium.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <numeric>
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
    std::vector<double> values;
    std::vector<bool> basic_columns;
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

bool exact_positive_witness(
    const ExactEqualitySystem& system,
    const HighsCandidate& candidate,
    std::size_t target_species
) {
    const std::size_t rank = system.matrix.size();
    if (rank == 0 || rank > system.columns
        || target_species >= system.columns
        || candidate.values.size() != system.columns) {
        return false;
    }
    std::vector<std::size_t> pool = ordered_columns(candidate);
    pool.erase(std::remove(pool.begin(), pool.end(), target_species), pool.end());
    std::vector<std::size_t> selected{target_species};
    bool validated = false;
    std::size_t visited = 0;
    visit_combinations(
        pool,
        rank - 1,
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
                )
                || basis_values.front() <= 0) {
                return false;
            }
            validated = true;
            return true;
        }
    );
    return validated;
}

bool exact_zero_certificate(
    const ExactEqualitySystem& system,
    std::size_t target_species,
    const HighsCandidate& candidate
) {
    const std::size_t rank = system.matrix.size();
    if (rank > system.columns || target_species >= system.columns) {
        return false;
    }
    const std::vector<std::size_t> pool = ordered_columns(candidate);
    std::vector<std::size_t> selected;
    bool validated = false;
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
                active_objective[equation] =
                    species == target_species ? ExactRational(1) : ExactRational(0);
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
                ExactRational reduced =
                    species == target_species ? ExactRational(-1) : ExactRational(0);
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
            validated = true;
            return true;
        }
    );
    return validated;
}

HighsCandidate solve_candidate_lp(
    const ExactEqualitySystem& system,
    std::size_t target_species
) {
    HighsCandidate result;
    if (target_species >= system.columns) {
        return result;
    }
    HighsModel model;
    model.lp_.num_col_ = static_cast<HighsInt>(system.columns);
    model.lp_.num_row_ = static_cast<HighsInt>(system.matrix.size());
    model.lp_.sense_ = ObjSense::kMaximize;
    model.lp_.col_cost_.assign(system.columns, 0.0);
    model.lp_.col_cost_[target_species] = 1.0;
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
    if (model_status != HighsModelStatus::kOptimal
        || highs.getInfo().primal_solution_status != kSolutionStatusFeasible) {
        return result;
    }
    const HighsSolution& solution = highs.getSolution();
    if (solution.col_value.size() != system.columns) {
        return result;
    }
    result.optimal = true;
    result.values = solution.col_value;
    result.basic_columns.assign(system.columns, false);
    const HighsBasis& basis = highs.getBasis();
    if (basis.valid && basis.col_status.size() == system.columns) {
        for (std::size_t column = 0; column < system.columns; ++column) {
            result.basic_columns[column] =
                basis.col_status[column] == HighsBasisStatus::kBasic;
        }
    }
    return result;
}

}  // namespace

HomogeneousSupportAnalysis analyze_homogeneous_support(
    const DenseMatrix& balance_matrix,
    const std::vector<double>& balance_totals,
    const std::vector<int>& charges
) {
    const std::size_t species_count = balance_matrix.columns;
    if (species_count == 0
        || balance_matrix.values.size() != balance_matrix.rows * species_count
        || balance_totals.size() != balance_matrix.rows
        || charges.size() != species_count
        || !std::all_of(
            balance_matrix.values.begin(),
            balance_matrix.values.end(),
            [](double value) { return std::isfinite(value); }
        )
        || !std::all_of(
            balance_totals.begin(),
            balance_totals.end(),
            [](double value) { return std::isfinite(value); }
        )) {
        throw std::invalid_argument("homogeneous support input is invalid");
    }

    const ExactEqualitySystem exact = independent_exact_equalities(
        balance_matrix, balance_totals, charges
    );

    HomogeneousSupportAnalysis result;
    result.classifications.assign(species_count, "unresolved");
    bool complete = true;
    bool has_accessible_species = false;
    for (std::size_t species = 0; species < species_count; ++species) {
        const HighsCandidate candidate = solve_candidate_lp(exact, species);
        if (!candidate.optimal) {
            complete = false;
            continue;
        }
        if (exact_positive_witness(exact, candidate, species)) {
            result.classifications[species] = "proved_accessible";
            has_accessible_species = true;
            continue;
        }
        if (exact_zero_certificate(exact, species, candidate)) {
            result.classifications[species] = "proved_structural_zero";
            continue;
        }
        complete = false;
    }
    result.exact_certificates_complete = complete && has_accessible_species;
    return result;
}

}  // namespace epcsaft_equilibrium
