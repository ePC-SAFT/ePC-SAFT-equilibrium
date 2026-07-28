#include "chemical_equilibrium.hpp"
#include "chemical_equilibrium_reference.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

namespace epcsaft_equilibrium {
namespace {

constexpr double kResidualMultiplier = 4096.0;

void require_finite(double value, const char* field) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(field) + " must be finite");
    }
}

void require_finite_vector(const std::vector<double>& values, const char* field) {
    if (!std::all_of(values.begin(), values.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::invalid_argument(std::string(field) + " must be finite");
    }
}

double matrix_scale(const DenseMatrix& matrix) {
    double scale = 0.0;
    for (double value : matrix.values) {
        scale = std::max(scale, std::abs(value));
    }
    return scale;
}

double numerical_tolerance(double scale, std::size_t dimension) {
    return kResidualMultiplier * std::numeric_limits<double>::epsilon()
        * std::max(1.0, scale) * static_cast<double>(std::max<std::size_t>(1, dimension));
}

double row_dot(const std::vector<double>& left, const std::vector<double>& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument("row dimensions do not match");
    }
    return std::inner_product(left.begin(), left.end(), right.begin(), 0.0);
}

double row_norm(const std::vector<double>& row) {
    return std::sqrt(row_dot(row, row));
}

std::vector<double> matrix_row(const DenseMatrix& matrix, std::size_t row) {
    std::vector<double> result(matrix.columns, 0.0);
    for (std::size_t column = 0; column < matrix.columns; ++column) {
        result[column] = matrix(row, column);
    }
    return result;
}

DenseMatrix matrix_from_rows(const std::vector<std::vector<double>>& rows) {
    if (rows.empty()) {
        return {};
    }
    const std::size_t columns = rows.front().size();
    DenseMatrix result{rows.size(), columns, {}};
    result.values.reserve(rows.size() * columns);
    for (const std::vector<double>& row : rows) {
        if (row.size() != columns) {
            throw std::invalid_argument("rows must have a common column count");
        }
        result.values.insert(result.values.end(), row.begin(), row.end());
    }
    return result;
}

std::vector<std::vector<double>> orthonormal_rows(
    const std::vector<std::vector<double>>& rows
) {
    std::vector<std::vector<double>> result;
    result.reserve(rows.size());
    for (const std::vector<double>& row : rows) {
        std::vector<double> residual = row;
        for (int pass = 0; pass < 2; ++pass) {
            for (const std::vector<double>& basis : result) {
                const double projection = row_dot(residual, basis);
                for (std::size_t column = 0; column < residual.size(); ++column) {
                    residual[column] -= projection * basis[column];
                }
            }
        }
        const double residual_norm = row_norm(residual);
        if (residual_norm <= numerical_tolerance(1.0, row.size())) {
            throw std::invalid_argument("row basis is numerically deficient");
        }
        for (double& value : residual) {
            value /= residual_norm;
        }
        result.push_back(std::move(residual));
    }
    return result;
}

double scaled_residual_norm(
    const std::vector<double>& row,
    const std::vector<std::vector<double>>& basis
) {
    std::vector<double> residual = row;
    for (int pass = 0; pass < 2; ++pass) {
        for (const std::vector<double>& vector : basis) {
            const double projection = row_dot(residual, vector);
            for (std::size_t column = 0; column < residual.size(); ++column) {
                residual[column] -= projection * vector[column];
            }
        }
    }
    const double scale = row_norm(row);
    return scale == 0.0 ? 0.0 : row_norm(residual) / scale;
}

std::vector<std::size_t> select_independent_rows(const DenseMatrix& matrix) {
    if (matrix.values.size() != matrix.rows * matrix.columns) {
        throw std::invalid_argument("matrix dimensions do not match storage");
    }
    std::vector<std::size_t> selected;
    std::vector<std::vector<double>> basis;
    std::vector<bool> used(matrix.rows, false);
    const double tolerance = numerical_tolerance(
        1.0, std::max(matrix.rows, matrix.columns)
    );
    while (selected.size() < matrix.columns) {
        std::size_t best = matrix.rows;
        double best_score = tolerance;
        for (std::size_t candidate = 0; candidate < matrix.rows; ++candidate) {
            if (used[candidate]) {
                continue;
            }
            const double score = scaled_residual_norm(matrix_row(matrix, candidate), basis);
            if (score > best_score) {
                best = candidate;
                best_score = score;
            }
        }
        if (best == matrix.rows) {
            break;
        }
        used[best] = true;
        selected.push_back(best);
        basis = orthonormal_rows([&] {
            std::vector<std::vector<double>> rows;
            rows.reserve(selected.size());
            for (std::size_t index : selected) {
                rows.push_back(matrix_row(matrix, index));
            }
            return rows;
        }());
    }
    return selected;
}

struct RowQr {
    DenseMatrix upper;
    std::vector<std::vector<double>> orthonormal;
};

RowQr factor_row_basis(const DenseMatrix& basis) {
    RowQr result{
        DenseMatrix{basis.rows, basis.rows, std::vector<double>(basis.rows * basis.rows, 0.0)},
        {},
    };
    result.orthonormal.reserve(basis.rows);
    for (std::size_t column = 0; column < basis.rows; ++column) {
        std::vector<double> residual = matrix_row(basis, column);
        for (int pass = 0; pass < 2; ++pass) {
            for (std::size_t row = 0; row < result.orthonormal.size(); ++row) {
                const double projection = row_dot(residual, result.orthonormal[row]);
                result.upper(row, column) += projection;
                for (std::size_t species = 0; species < residual.size(); ++species) {
                    residual[species] -= projection * result.orthonormal[row][species];
                }
            }
        }
        const double diagonal = row_norm(residual);
        if (diagonal <= numerical_tolerance(1.0, basis.columns)) {
            throw std::invalid_argument("reaction matrix rank is deficient");
        }
        result.upper(column, column) = diagonal;
        for (double& value : residual) {
            value /= diagonal;
        }
        result.orthonormal.push_back(std::move(residual));
    }
    return result;
}

std::vector<double> row_coordinates(
    const std::vector<double>& row,
    const RowQr& factor
) {
    std::vector<double> projected(factor.orthonormal.size(), 0.0);
    for (std::size_t index = 0; index < factor.orthonormal.size(); ++index) {
        projected[index] = row_dot(row, factor.orthonormal[index]);
    }
    std::vector<double> coordinates(projected.size(), 0.0);
    for (std::size_t reverse = projected.size(); reverse > 0; --reverse) {
        const std::size_t row_index = reverse - 1;
        double value = projected[row_index];
        for (std::size_t column = row_index + 1; column < coordinates.size(); ++column) {
            value -= factor.upper(row_index, column) * coordinates[column];
        }
        coordinates[row_index] = value / factor.upper(row_index, row_index);
    }
    return coordinates;
}

double charge_reaction_inf_norm(
    const std::vector<int>& charges,
    const DenseMatrix& reaction_matrix
) {
    double norm = 0.0;
    for (std::size_t reaction = 0; reaction < reaction_matrix.rows; ++reaction) {
        double value = 0.0;
        for (std::size_t species = 0; species < reaction_matrix.columns; ++species) {
            value += static_cast<double>(charges[species])
                * reaction_matrix(reaction, species);
        }
        norm = std::max(norm, std::abs(value));
    }
    return norm;
}

void validate_identities(const ReactionSystemInput& input) {
    const std::size_t species_count = input.species_ids.size();
    if (species_count < 2) {
        throw std::invalid_argument("reaction system requires at least two species");
    }
    std::unordered_set<std::string> unique_ids;
    for (const std::string& species_id : input.species_ids) {
        if (species_id.empty() || !unique_ids.insert(species_id).second) {
            throw std::invalid_argument("species identities must be nonempty and unique");
        }
    }
    if (input.charges.size() != species_count) {
        throw std::invalid_argument("species charges do not match species count");
    }
    if (input.provider_fingerprint.empty()) {
        throw std::invalid_argument("Provider fingerprint is incomplete");
    }
}

void validate_reference_records(const ReactionSystemInput& input) {
    if (input.equilibrium_constant_records.size() != input.ln_k.size()) {
        throw std::invalid_argument("equilibrium-constant source records are incomplete");
    }
    for (const EquilibriumConstantRecord& record : input.equilibrium_constant_records) {
        if (record.source_id.empty()) {
            throw std::invalid_argument("equilibrium constant source identity is incomplete");
        }
        if (record.reference_id.empty()) {
            throw std::invalid_argument("equilibrium constant reference identity is incomplete");
        }
        if (record.reaction_orientation != "products_positive") {
            throw std::invalid_argument(
                "equilibrium constant reaction orientation must be products_positive"
            );
        }
        if (record.conversion_id != "already-provider-basis") {
            throw std::invalid_argument(
                "equilibrium constant conversion provenance is not already-provider-basis"
            );
        }
        if (record.reference_id != "provider-helmholtz-coordinate-basis") {
            throw std::invalid_argument(
                "equilibrium constant reference is not the Provider chemical basis"
            );
        }
        if (!record.dimensionless) {
            throw std::invalid_argument("lnK must be explicitly dimensionless");
        }
        require_finite(record.temperature_k, "equilibrium constant temperature");
        require_finite(record.pressure_pa, "equilibrium constant pressure");
        if (record.temperature_k != input.temperature_k) {
            throw std::invalid_argument("equilibrium constant temperature binding does not match");
        }
        if (record.pressure_pa != input.pressure_pa) {
            throw std::invalid_argument("equilibrium constant pressure binding does not match");
        }
    }
}

DenseMatrix selected_rows(
    const DenseMatrix& matrix,
    const std::vector<std::size_t>& row_indices
) {
    std::vector<std::vector<double>> rows;
    rows.reserve(row_indices.size());
    for (std::size_t row : row_indices) {
        rows.push_back(matrix_row(matrix, row));
    }
    return matrix_from_rows(rows);
}

double reaction_conservation_inf_norm(
    const std::vector<double>& molar_masses,
    const DenseMatrix& balance_matrix,
    const DenseMatrix& reaction_matrix
) {
    double norm = 0.0;
    for (std::size_t reaction = 0; reaction < reaction_matrix.rows; ++reaction) {
        double mass_value = 0.0;
        for (std::size_t species = 0; species < reaction_matrix.columns; ++species) {
            mass_value += molar_masses[species] * reaction_matrix(reaction, species);
        }
        norm = std::max(norm, std::abs(mass_value));
        for (std::size_t balance = 0; balance < balance_matrix.rows; ++balance) {
            double balance_value = 0.0;
            for (std::size_t species = 0; species < reaction_matrix.columns; ++species) {
                balance_value += balance_matrix(balance, species)
                    * reaction_matrix(reaction, species);
            }
            norm = std::max(norm, std::abs(balance_value));
        }
    }
    return norm;
}

std::vector<std::vector<double>> seeded_balance_rows(
    const std::vector<double>& molar_masses,
    const std::vector<int>& charges,
    const DenseMatrix& supplied_balance_matrix,
    bool& charge_is_independent
) {
    std::vector<std::vector<double>> basis_rows;
    std::vector<std::vector<double>> orthonormal;
    const double rank_tolerance = numerical_tolerance(1.0, molar_masses.size());
    const auto add_seed = [&](const std::vector<double>& row) {
        const double score = scaled_residual_norm(row, orthonormal);
        if (score <= rank_tolerance) {
            return false;
        }
        basis_rows.push_back(row);
        orthonormal = orthonormal_rows(basis_rows);
        return true;
    };

    add_seed(molar_masses);
    std::vector<double> charge_row(molar_masses.size(), 0.0);
    for (std::size_t species = 0; species < charges.size(); ++species) {
        charge_row[species] = static_cast<double>(charges[species]);
    }
    charge_is_independent = !std::all_of(
        charge_row.begin(), charge_row.end(), [](double value) { return value == 0.0; }
    ) && add_seed(charge_row);

    std::vector<std::vector<double>> retained_balance_rows;
    for (std::size_t row = 0; row < supplied_balance_matrix.rows; ++row) {
        const std::vector<double> candidate = matrix_row(
            supplied_balance_matrix, row
        );
        const double score = scaled_residual_norm(candidate, orthonormal);
        if (score <= rank_tolerance) {
            continue;
        }
        retained_balance_rows.push_back(candidate);
        basis_rows.push_back(candidate);
        orthonormal = orthonormal_rows(basis_rows);
    }

    std::vector<std::vector<double>> result;
    result.reserve(1 + retained_balance_rows.size());
    result.push_back(molar_masses);
    result.insert(result.end(), retained_balance_rows.begin(), retained_balance_rows.end());
    return result;
}

double matrix_reconstruction_inf_norm(
    const DenseMatrix& transform,
    const DenseMatrix& basis,
    const DenseMatrix& supplied
) {
    double norm = 0.0;
    for (std::size_t row = 0; row < supplied.rows; ++row) {
        for (std::size_t species = 0; species < supplied.columns; ++species) {
            double value = supplied(row, species);
            for (std::size_t basis_row = 0; basis_row < basis.rows; ++basis_row) {
                value -= transform(row, basis_row) * basis(basis_row, species);
            }
            norm = std::max(norm, std::abs(value));
        }
    }
    return norm;
}

double constant_reconstruction_inf_norm(
    const DenseMatrix& transform,
    const std::vector<double>& independent,
    const std::vector<double>& supplied
) {
    double norm = 0.0;
    for (std::size_t row = 0; row < supplied.size(); ++row) {
        double value = -supplied[row];
        for (std::size_t basis = 0; basis < independent.size(); ++basis) {
            value += transform(row, basis) * independent[basis];
        }
        norm = std::max(norm, std::abs(value));
    }
    return norm;
}

DenseMatrix identity_matrix(std::size_t dimension) {
    DenseMatrix result{
        dimension,
        dimension,
        std::vector<double>(dimension * dimension, 0.0),
    };
    for (std::size_t index = 0; index < dimension; ++index) {
        result(index, index) = 1.0;
    }
    return result;
}

DenseMatrix restrict_columns(
    const DenseMatrix& matrix,
    const std::vector<std::size_t>& retained_columns
) {
    DenseMatrix result{
        matrix.rows,
        retained_columns.size(),
        std::vector<double>(matrix.rows * retained_columns.size(), 0.0),
    };
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        for (std::size_t column = 0; column < retained_columns.size(); ++column) {
            result(row, column) = matrix(row, retained_columns[column]);
        }
    }
    return result;
}

DenseMatrix multiply_matrices(const DenseMatrix& left, const DenseMatrix& right) {
    if (left.columns != right.rows) {
        throw std::invalid_argument("accessible reaction matrix dimensions are inconsistent");
    }
    DenseMatrix result{
        left.rows,
        right.columns,
        std::vector<double>(left.rows * right.columns, 0.0),
    };
    for (std::size_t row = 0; row < left.rows; ++row) {
        for (std::size_t inner = 0; inner < left.columns; ++inner) {
            for (std::size_t column = 0; column < right.columns; ++column) {
                result(row, column) += left(row, inner) * right(inner, column);
            }
        }
    }
    return result;
}

std::vector<double> multiply_matrix_vector(
    const DenseMatrix& matrix,
    const std::vector<double>& vector
) {
    if (matrix.columns != vector.size()) {
        throw std::invalid_argument("accessible reaction vector dimensions are inconsistent");
    }
    std::vector<double> result(matrix.rows, 0.0);
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        for (std::size_t column = 0; column < matrix.columns; ++column) {
            result[row] += matrix(row, column) * vector[column];
        }
    }
    return result;
}

DenseMatrix nullspace_basis(const DenseMatrix& matrix) {
    DenseMatrix reduced = matrix;
    const double tolerance = numerical_tolerance(
        matrix_scale(matrix), std::max(matrix.rows, matrix.columns)
    );
    std::vector<std::size_t> pivot_columns;
    std::size_t pivot_row = 0;
    for (std::size_t column = 0;
         column < reduced.columns && pivot_row < reduced.rows;
         ++column) {
        std::size_t pivot = pivot_row;
        for (std::size_t row = pivot_row + 1; row < reduced.rows; ++row) {
            if (std::abs(reduced(row, column)) > std::abs(reduced(pivot, column))) {
                pivot = row;
            }
        }
        if (std::abs(reduced(pivot, column)) <= tolerance) {
            continue;
        }
        if (pivot != pivot_row) {
            for (std::size_t index = 0; index < reduced.columns; ++index) {
                std::swap(reduced(pivot, index), reduced(pivot_row, index));
            }
        }
        const double diagonal = reduced(pivot_row, column);
        for (std::size_t index = column; index < reduced.columns; ++index) {
            reduced(pivot_row, index) /= diagonal;
        }
        for (std::size_t row = 0; row < reduced.rows; ++row) {
            if (row == pivot_row) {
                continue;
            }
            const double factor = reduced(row, column);
            for (std::size_t index = column; index < reduced.columns; ++index) {
                reduced(row, index) -= factor * reduced(pivot_row, index);
            }
        }
        pivot_columns.push_back(column);
        ++pivot_row;
    }

    std::vector<bool> is_pivot(reduced.columns, false);
    for (std::size_t column : pivot_columns) {
        is_pivot[column] = true;
    }
    std::vector<std::vector<double>> rows;
    for (std::size_t free_column = 0; free_column < reduced.columns; ++free_column) {
        if (is_pivot[free_column]) {
            continue;
        }
        std::vector<double> vector(reduced.columns, 0.0);
        vector[free_column] = 1.0;
        for (std::size_t row = 0; row < pivot_columns.size(); ++row) {
            vector[pivot_columns[row]] = -reduced(row, free_column);
        }
        rows.push_back(std::move(vector));
    }
    return matrix_from_rows(rows);
}

double removed_species_residual(
    const DenseMatrix& full_reaction_matrix,
    const std::vector<std::size_t>& removed_species
) {
    double result = 0.0;
    for (std::size_t reaction = 0; reaction < full_reaction_matrix.rows; ++reaction) {
        for (std::size_t species : removed_species) {
            result = std::max(
                result, std::abs(full_reaction_matrix(reaction, species))
            );
        }
    }
    return result;
}

CompiledReactionSystem compile_accessible_face(CompiledReactionSystem system) {
    system.original_species_ids = system.species_ids;
    system.original_charges = system.charges;
    system.original_molar_masses_kg_per_mol = system.molar_masses_kg_per_mol;
    system.original_feed_amounts = system.feed_amounts;
    system.support = analyze_homogeneous_support(
        system.balance_matrix,
        system.balance_totals,
        system.charges
    );
    if (system.support.classifications.size() != system.species_ids.size()) {
        throw std::runtime_error("homogeneous support evidence is incomplete");
    }
    for (std::size_t species = 0; species < system.species_ids.size(); ++species) {
        if (system.support.classifications[species] == "proved_structural_zero") {
            system.removed_species_indices.push_back(species);
        } else {
            system.retained_species_indices.push_back(species);
        }
    }
    if (system.removed_species_indices.empty()) {
        system.accessible_reaction_transform = identity_matrix(
            system.reaction_matrix.rows
        );
        return system;
    }

    DenseMatrix removed_coefficients{
        system.removed_species_indices.size(),
        system.reaction_matrix.rows,
        std::vector<double>(
            system.removed_species_indices.size() * system.reaction_matrix.rows,
            0.0
        ),
    };
    for (std::size_t removed = 0;
         removed < system.removed_species_indices.size();
         ++removed) {
        for (std::size_t reaction = 0;
             reaction < system.reaction_matrix.rows;
             ++reaction) {
            removed_coefficients(removed, reaction) = system.reaction_matrix(
                reaction, system.removed_species_indices[removed]
            );
        }
    }
    DenseMatrix accessible_transform = nullspace_basis(removed_coefficients);
    if (accessible_transform.rows == 0) {
        throw std::invalid_argument(
            "accessible chemical face has no retained reaction direction"
        );
    }
    const DenseMatrix accessible_full_reactions = multiply_matrices(
        accessible_transform, system.reaction_matrix
    );
    const double removed_residual = removed_species_residual(
        accessible_full_reactions, system.removed_species_indices
    );
    if (removed_residual > numerical_tolerance(
            matrix_scale(system.reaction_matrix), system.species_ids.size()
        )) {
        throw std::invalid_argument(
            "accessible reaction transform does not cancel removed species"
        );
    }
    DenseMatrix accessible_reactions = restrict_columns(
        accessible_full_reactions, system.retained_species_indices
    );
    std::vector<double> accessible_ln_k = multiply_matrix_vector(
        accessible_transform, system.ln_k
    );

    std::vector<std::string> species_ids;
    std::vector<int> charges;
    std::vector<double> molar_masses;
    std::vector<double> feed_amounts;
    species_ids.reserve(system.retained_species_indices.size());
    charges.reserve(system.retained_species_indices.size());
    molar_masses.reserve(system.retained_species_indices.size());
    feed_amounts.reserve(system.retained_species_indices.size());
    for (std::size_t original : system.retained_species_indices) {
        species_ids.push_back(system.species_ids[original]);
        charges.push_back(system.charges[original]);
        molar_masses.push_back(system.molar_masses_kg_per_mol[original]);
        feed_amounts.push_back(system.feed_amounts[original]);
    }
    const DenseMatrix restricted_balances = restrict_columns(
        system.balance_matrix, system.retained_species_indices
    );
    bool charge_is_independent = false;
    const DenseMatrix balance_matrix = matrix_from_rows(
        seeded_balance_rows(
            molar_masses, charges, restricted_balances, charge_is_independent
        )
    );
    if (balance_matrix.rows + (charge_is_independent ? 1 : 0)
            + accessible_reactions.rows
        != species_ids.size()) {
        throw std::invalid_argument(
            "accessible closed system reaction matrix rank sum does not equal species count"
        );
    }
    const std::vector<std::size_t> independent_accessible_rows =
        select_independent_rows(accessible_reactions);
    if (independent_accessible_rows.size() != accessible_reactions.rows) {
        throw std::invalid_argument("accessible reaction matrix rank is deficient");
    }

    std::vector<double> balance_totals(balance_matrix.rows, 0.0);
    for (std::size_t row = 0; row < balance_matrix.rows; ++row) {
        for (std::size_t species = 0; species < species_ids.size(); ++species) {
            balance_totals[row] +=
                balance_matrix(row, species) * feed_amounts[species];
        }
    }
    const double conservation_norm = reaction_conservation_inf_norm(
        molar_masses, balance_matrix, accessible_reactions
    );
    if (conservation_norm > numerical_tolerance(
            std::max(
                *std::max_element(molar_masses.begin(), molar_masses.end()),
                matrix_scale(balance_matrix)
            ) * matrix_scale(accessible_reactions),
            species_ids.size()
        )) {
        throw std::invalid_argument(
            "accessible reaction stoichiometry does not conserve balances"
        );
    }
    const double charge_norm = charge_reaction_inf_norm(
        charges, accessible_reactions
    );
    if (charge_norm > numerical_tolerance(
            matrix_scale(accessible_reactions), species_ids.size()
        )) {
        throw std::invalid_argument(
            "accessible reaction stoichiometry does not conserve charge"
        );
    }
    ReactionReferenceReconstruction reconstruction = reconstruct_reaction_reference(
        accessible_reactions, accessible_ln_k
    );
    const double reference_residual = reconstruction.residual_inf_norm;
    double ln_k_scale = 0.0;
    for (double value : accessible_ln_k) {
        ln_k_scale = std::max(ln_k_scale, std::abs(value));
    }
    if (reference_residual > numerical_tolerance(ln_k_scale, species_ids.size())) {
        throw std::invalid_argument(
            "accessible chemical reference reconstruction is inconsistent"
        );
    }

    system.species_ids = std::move(species_ids);
    system.charges = std::move(charges);
    system.molar_masses_kg_per_mol = std::move(molar_masses);
    system.balance_matrix = balance_matrix;
    system.reaction_matrix = std::move(accessible_reactions);
    system.balance_totals = std::move(balance_totals);
    system.feed_amounts = std::move(feed_amounts);
    system.ln_k = std::move(accessible_ln_k);
    system.g_ref = std::move(reconstruction.reference);
    system.accessible_reaction_transform = std::move(accessible_transform);
    system.balance_rank = system.balance_matrix.rows;
    system.reaction_rank = system.reaction_matrix.rows;
    return system;
}

}  // namespace


CompiledReactionSystem compile_reaction_system(const ReactionSystemInput& input) {
    validate_identities(input);
    require_finite(input.temperature_k, "temperature");
    require_finite(input.pressure_pa, "pressure");
    if (input.temperature_k <= 0.0 || input.pressure_pa <= 0.0) {
        throw std::invalid_argument("temperature and pressure must be positive");
    }
    const std::size_t species_count = input.species_ids.size();
    if (input.reaction_matrix.rows == 0) {
        throw std::invalid_argument(
            "complete closed system requires independent reactions and lnK"
        );
    }
    if (input.feed_amounts.size() != species_count
        || input.molar_masses_kg_per_mol.size() != species_count
        || input.balance_matrix.columns != species_count
        || input.reaction_matrix.columns != species_count) {
        throw std::invalid_argument("reaction-system dimensions do not match species count");
    }
    require_finite_vector(input.feed_amounts, "feed amounts");
    if (!std::all_of(input.feed_amounts.begin(), input.feed_amounts.end(), [](double value) {
            return value >= 0.0;
        })) {
        throw std::invalid_argument("feed amounts must be nonnegative");
    }
    if (std::accumulate(input.feed_amounts.begin(), input.feed_amounts.end(), 0.0) <= 0.0) {
        throw std::invalid_argument("feed must contain a positive total amount");
    }
    require_finite_vector(input.balance_matrix.values, "balance matrix");
    require_finite_vector(input.reaction_matrix.values, "reaction matrix");
    require_finite_vector(input.ln_k, "lnK");
    require_finite_vector(input.molar_masses_kg_per_mol, "molar masses");
    if (!std::all_of(
            input.molar_masses_kg_per_mol.begin(),
            input.molar_masses_kg_per_mol.end(),
            [](double value) { return value > 0.0; }
        )) {
        throw std::invalid_argument("molar masses must be finite and positive");
    }
    if (input.ln_k.size() != input.reaction_matrix.rows) {
        throw std::invalid_argument("equilibrium constants are incomplete for the reaction matrix");
    }
    validate_reference_records(input);

    const double conservation_norm = reaction_conservation_inf_norm(
        input.molar_masses_kg_per_mol,
        input.balance_matrix,
        input.reaction_matrix
    );
    const double molar_mass_scale = *std::max_element(
        input.molar_masses_kg_per_mol.begin(),
        input.molar_masses_kg_per_mol.end()
    );
    const double conservation_tolerance = numerical_tolerance(
        std::max(molar_mass_scale, matrix_scale(input.balance_matrix))
            * matrix_scale(input.reaction_matrix),
        species_count
    );
    if (conservation_norm > conservation_tolerance) {
        throw std::invalid_argument(
            "reaction stoichiometry does not conserve molar mass and the balance matrix"
        );
    }
    double feed_charge = 0.0;
    for (std::size_t species = 0; species < species_count; ++species) {
        feed_charge += static_cast<double>(input.charges[species]) * input.feed_amounts[species];
    }
    if (std::abs(feed_charge) > numerical_tolerance(
            *std::max_element(input.feed_amounts.begin(), input.feed_amounts.end()), species_count
        )) {
        throw std::invalid_argument("feed must be exactly electroneutral within numerical precision");
    }
    const double reaction_charge_norm = charge_reaction_inf_norm(
        input.charges, input.reaction_matrix
    );
    if (reaction_charge_norm > numerical_tolerance(matrix_scale(input.reaction_matrix), species_count)) {
        throw std::invalid_argument("reaction stoichiometry does not conserve charge");
    }

    const std::vector<std::size_t> reaction_basis_rows = select_independent_rows(
        input.reaction_matrix
    );
    if (reaction_basis_rows.empty()) {
        throw std::invalid_argument("reaction matrix rank is deficient");
    }
    const DenseMatrix reaction_matrix = selected_rows(
        input.reaction_matrix, reaction_basis_rows
    );
    const RowQr reaction_factor = factor_row_basis(reaction_matrix);
    DenseMatrix reaction_transform{
        input.reaction_matrix.rows,
        reaction_matrix.rows,
        std::vector<double>(
            input.reaction_matrix.rows * reaction_matrix.rows, 0.0
        ),
    };
    for (std::size_t supplied = 0; supplied < input.reaction_matrix.rows; ++supplied) {
        const auto basis_position = std::find(
            reaction_basis_rows.begin(), reaction_basis_rows.end(), supplied
        );
        if (basis_position != reaction_basis_rows.end()) {
            reaction_transform(
                supplied,
                static_cast<std::size_t>(
                    std::distance(reaction_basis_rows.begin(), basis_position)
                )
            ) = 1.0;
            continue;
        }
        const std::vector<double> coordinates = row_coordinates(
            matrix_row(input.reaction_matrix, supplied), reaction_factor
        );
        for (std::size_t basis = 0; basis < coordinates.size(); ++basis) {
            reaction_transform(supplied, basis) = coordinates[basis];
        }
    }
    const double reaction_transform_norm = matrix_reconstruction_inf_norm(
        reaction_transform, reaction_matrix, input.reaction_matrix
    );
    if (reaction_transform_norm > numerical_tolerance(
            matrix_scale(input.reaction_matrix), species_count
        )) {
        throw std::invalid_argument("reaction matrix rank reduction is inconsistent");
    }
    std::vector<double> independent_ln_k;
    independent_ln_k.reserve(reaction_basis_rows.size());
    for (std::size_t row : reaction_basis_rows) {
        independent_ln_k.push_back(input.ln_k[row]);
    }
    const double reaction_cycle_norm = constant_reconstruction_inf_norm(
        reaction_transform, independent_ln_k, input.ln_k
    );
    double ln_k_scale = 0.0;
    for (double value : input.ln_k) {
        ln_k_scale = std::max(ln_k_scale, std::abs(value));
    }
    if (reaction_cycle_norm > numerical_tolerance(
            ln_k_scale, std::max(input.reaction_matrix.rows, species_count)
        )) {
        throw std::invalid_argument(
            "reaction constant cycle is inconsistent after Provider-basis conversion"
        );
    }

    bool charge_is_independent = false;
    const DenseMatrix balance_matrix = matrix_from_rows(
        seeded_balance_rows(
            input.molar_masses_kg_per_mol,
            input.charges,
            input.balance_matrix,
            charge_is_independent
        )
    );
    const std::size_t balance_rank = balance_matrix.rows;
    const std::size_t reaction_rank = reaction_matrix.rows;
    if (balance_rank + (charge_is_independent ? 1 : 0) + reaction_rank
        != species_count) {
        throw std::invalid_argument(
            "complete closed system reaction matrix rank sum does not equal species count"
        );
    }

    ReactionReferenceReconstruction reconstruction = reconstruct_reaction_reference(
        reaction_matrix, independent_ln_k
    );

    const double reference_residual = reconstruction.residual_inf_norm;
    if (reference_residual > numerical_tolerance(ln_k_scale, species_count)) {
        throw std::invalid_argument("standard chemical reference reconstruction is inconsistent");
    }

    std::vector<double> balance_totals(balance_matrix.rows, 0.0);
    for (std::size_t row = 0; row < balance_matrix.rows; ++row) {
        for (std::size_t species = 0; species < species_count; ++species) {
            balance_totals[row] += balance_matrix(row, species)
                * input.feed_amounts[species];
        }
    }
    CompiledReactionSystem result;
    result.species_ids = input.species_ids;
    result.charges = input.charges;
    result.molar_masses_kg_per_mol = input.molar_masses_kg_per_mol;
    result.supplied_balance_matrix = input.balance_matrix;
    result.reaction_transform = std::move(reaction_transform);
    result.reaction_basis_rows = reaction_basis_rows;
    result.balance_matrix = balance_matrix;
    result.reaction_matrix = reaction_matrix;
    result.balance_totals = std::move(balance_totals);
    result.feed_amounts = input.feed_amounts;
    result.ln_k = std::move(independent_ln_k);
    result.g_ref = std::move(reconstruction.reference);
    result.balance_rank = balance_rank;
    result.reaction_rank = reaction_rank;
    return compile_accessible_face(std::move(result));
}

}  // namespace epcsaft_equilibrium
