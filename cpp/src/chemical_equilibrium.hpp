#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace epcsaft_equilibrium {

struct DenseMatrix {
    std::size_t rows = 0;
    std::size_t columns = 0;
    std::vector<double> values;

    [[nodiscard]] double operator()(std::size_t row, std::size_t column) const;
    double& operator()(std::size_t row, std::size_t column);
};

struct EquilibriumConstantRecord {
    std::string source_id;
    std::string reference_id;
    bool dimensionless = false;
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    std::string pressure_binding;
};

struct ReactionSystemInput {
    std::vector<std::string> species_ids;
    std::vector<int> charges;
    std::vector<std::string> provider_component_ids;
    std::vector<int> provider_charges;
    std::string provider_fingerprint;
    std::vector<std::string> expected_provider_component_ids;
    std::vector<int> expected_provider_charges;
    std::string expected_provider_fingerprint;
    DenseMatrix balance_matrix;
    std::size_t declared_balance_rank = 0;
    DenseMatrix reaction_matrix;
    std::vector<double> feed_amounts;
    std::vector<double> ln_k;
    std::vector<EquilibriumConstantRecord> equilibrium_constant_records;
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    bool complete_closed_system = false;
};

struct CompiledReactionSystem {
    std::vector<std::string> species_ids;
    std::vector<int> charges;
    DenseMatrix balance_matrix;
    DenseMatrix reaction_matrix;
    std::vector<double> balance_totals;
    std::vector<double> feed_amounts;
    std::vector<double> ln_k;
    std::vector<double> g_ref;
    std::string provider_fingerprint;
    std::size_t balance_rank = 0;
    std::size_t reaction_rank = 0;
    double reference_reconstruction_inf_norm = 0.0;
    double conservation_reaction_inf_norm = 0.0;
    double charge_reaction_inf_norm = 0.0;
};

[[nodiscard]] CompiledReactionSystem compile_reaction_system(
    const ReactionSystemInput& input
);

}  // namespace epcsaft_equilibrium
