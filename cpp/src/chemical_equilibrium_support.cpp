#include "chemical_equilibrium.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

namespace epcsaft_equilibrium {
namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;

constexpr std::size_t kMaximumBasisCandidates = 100000;

struct ExactSystem {
    std::vector<std::vector<Rational>> matrix;
    std::vector<Rational> totals;
    std::size_t columns = 0;
};

Rational exact(double value) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument("exact support validation requires finite inputs");
    }
    if (value == 0.0) {
        return Rational(0);
    }
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(value));
    Integer numerator =
        (bits & ((std::uint64_t{1} << 52U) - 1U))
        | (std::uint64_t{1} << 52U);
    const std::uint64_t exponent = (bits >> 52U) & 0x7ffU;
    int shift = static_cast<int>(exponent) - 1023 - 52;
    if (exponent == 0U) {
        numerator = bits & ((std::uint64_t{1} << 52U) - 1U);
        shift = -1022 - 52;
    }
    if ((bits >> 63U) != 0U) {
        numerator = -numerator;
    }
    if (shift >= 0) {
        return Rational(numerator << shift);
    }
    return Rational(numerator, Integer(1) << -shift);
}

ExactSystem independent_equalities(
    const DenseMatrix& balances,
    const std::vector<double>& totals,
    const std::vector<int>& charges
) {
    const std::size_t columns = balances.columns;
    std::vector<std::vector<Rational>> rows;
    rows.reserve(balances.rows + 1);
    for (std::size_t row = 0; row < balances.rows; ++row) {
        rows.emplace_back(columns + 1, Rational(0));
        for (std::size_t column = 0; column < columns; ++column) {
            rows.back()[column] = exact(balances(row, column));
        }
        rows.back().back() = exact(totals[row]);
    }
    if (std::any_of(charges.begin(), charges.end(), [](int value) {
            return value != 0;
        })) {
        rows.emplace_back(columns + 1, Rational(0));
        for (std::size_t column = 0; column < columns; ++column) {
            rows.back()[column] = Rational(charges[column]);
        }
    }

    std::size_t rank = 0;
    for (std::size_t column = 0; column < columns && rank < rows.size(); ++column) {
        std::size_t pivot = rank;
        while (pivot < rows.size() && rows[pivot][column] == 0) {
            ++pivot;
        }
        if (pivot == rows.size()) {
            continue;
        }
        std::swap(rows[rank], rows[pivot]);
        const Rational diagonal = rows[rank][column];
        for (std::size_t index = column; index <= columns; ++index) {
            rows[rank][index] /= diagonal;
        }
        for (std::size_t row = 0; row < rows.size(); ++row) {
            if (row == rank || rows[row][column] == 0) {
                continue;
            }
            const Rational factor = rows[row][column];
            for (std::size_t index = column; index <= columns; ++index) {
                rows[row][index] -= factor * rows[rank][index];
            }
        }
        ++rank;
    }
    for (std::size_t row = rank; row < rows.size(); ++row) {
        if (std::all_of(
                rows[row].begin(),
                rows[row].begin() + static_cast<std::ptrdiff_t>(columns),
                [](const Rational& value) { return value == 0; }
            )
            && rows[row].back() != 0) {
            throw std::invalid_argument("homogeneous support equalities are inconsistent");
        }
    }
    ExactSystem result;
    result.columns = columns;
    result.matrix.assign(rank, std::vector<Rational>(columns, Rational(0)));
    result.totals.resize(rank);
    for (std::size_t row = 0; row < rank; ++row) {
        std::copy_n(rows[row].begin(), columns, result.matrix[row].begin());
        result.totals[row] = rows[row].back();
    }
    return result;
}

bool solve_basis(
    std::vector<std::vector<Rational>> matrix,
    std::vector<Rational> values,
    std::vector<Rational>& solution
) {
    for (std::size_t column = 0; column < matrix.size(); ++column) {
        std::size_t pivot = column;
        while (pivot < matrix.size() && matrix[pivot][column] == 0) {
            ++pivot;
        }
        if (pivot == matrix.size()) {
            return false;
        }
        std::swap(matrix[pivot], matrix[column]);
        std::swap(values[pivot], values[column]);
        const Rational diagonal = matrix[column][column];
        for (std::size_t index = column; index < matrix.size(); ++index) {
            matrix[column][index] /= diagonal;
        }
        values[column] /= diagonal;
        for (std::size_t row = 0; row < matrix.size(); ++row) {
            if (row == column || matrix[row][column] == 0) {
                continue;
            }
            const Rational factor = matrix[row][column];
            for (std::size_t index = column; index < matrix.size(); ++index) {
                matrix[row][index] -= factor * matrix[column][index];
            }
            values[row] -= factor * values[column];
        }
    }
    solution = std::move(values);
    return true;
}

bool visit_bases(
    std::size_t columns,
    std::size_t choose,
    std::size_t offset,
    std::vector<std::size_t>& selected,
    const std::function<bool(const std::vector<std::size_t>&)>& visitor
) {
    if (choose == 0) {
        return visitor(selected);
    }
    for (std::size_t column = offset; column + choose <= columns; ++column) {
        selected.push_back(column);
        if (visit_bases(columns, choose - 1, column + 1, selected, visitor)) {
            return true;
        }
        selected.pop_back();
    }
    return false;
}

}  // namespace

std::vector<std::size_t> homogeneous_structural_zeros(
    const DenseMatrix& balance_matrix,
    const std::vector<double>& balance_totals,
    const std::vector<int>& charges
) {
    const std::size_t species_count = balance_matrix.columns;
    if (species_count == 0 || balance_matrix.rows == 0
        || balance_matrix.values.size() != balance_matrix.rows * species_count
        || balance_totals.size() != balance_matrix.rows
        || charges.size() != species_count
        || !std::all_of(
            balance_matrix.values.begin(),
            balance_matrix.values.end(),
            [](double value) { return std::isfinite(value); }
        )
        || !std::all_of(
            balance_matrix.values.begin(),
            balance_matrix.values.begin() + static_cast<std::ptrdiff_t>(species_count),
            [](double value) { return value > 0.0; }
        )) {
        throw std::invalid_argument(
            "homogeneous support requires finite balances with a positive bounding row"
        );
    }
    const ExactSystem system = independent_equalities(
        balance_matrix, balance_totals, charges
    );
    std::vector<bool> accessible(species_count, false);
    std::vector<std::size_t> basis;
    std::size_t visited = 0;
    bool feasible = false;
    visit_bases(
        species_count,
        system.matrix.size(),
        0,
        basis,
        [&](const std::vector<std::size_t>& columns) {
            if (++visited > kMaximumBasisCandidates) {
                throw std::runtime_error(
                    "exact homogeneous support enumeration limit exceeded"
                );
            }
            std::vector<std::vector<Rational>> square(
                columns.size(), std::vector<Rational>(columns.size(), Rational(0))
            );
            for (std::size_t row = 0; row < columns.size(); ++row) {
                for (std::size_t column = 0; column < columns.size(); ++column) {
                    square[row][column] = system.matrix[row][columns[column]];
                }
            }
            std::vector<Rational> values;
            if (!solve_basis(std::move(square), system.totals, values)
                || std::any_of(values.begin(), values.end(), [](const Rational& value) {
                    return value < 0;
                })) {
                return false;
            }
            feasible = true;
            for (std::size_t column = 0; column < columns.size(); ++column) {
                accessible[columns[column]] =
                    accessible[columns[column]] || values[column] > 0;
            }
            return std::all_of(accessible.begin(), accessible.end(), [](bool value) {
                return value;
            });
        }
    );
    if (!feasible) {
        throw std::invalid_argument(
            "homogeneous support has no nonnegative feasible state"
        );
    }
    std::vector<std::size_t> result;
    for (std::size_t species = 0; species < species_count; ++species) {
        if (!accessible[species]) {
            result.push_back(species);
        }
    }
    return result;
}

}  // namespace epcsaft_equilibrium
