#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace epcsaft_equilibrium {

inline constexpr const char* kHeld2ManufacturedFormulationId =
    "perdomo-held2.modified-mole.manufactured.v1";

struct Held2Coordinates {
    std::vector<double> charges;
    std::size_t eliminated_index = 0;
    std::size_t dependent_index = 0;
    std::vector<std::size_t> retained_indices;
    std::vector<std::size_t> independent_indices;
    std::vector<double> modified_factors;
    std::vector<double> independent_lower_bounds;
    std::vector<double> independent_upper_bounds;
};

[[nodiscard]] Held2Coordinates make_held2_coordinates(
    const std::vector<double>& charges
);

[[nodiscard]] std::vector<double> held2_transform_physical_fractions(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical_fractions
);

[[nodiscard]] std::vector<double> held2_lift_modified_fractions(
    const Held2Coordinates& coordinates,
    const std::vector<double>& modified_fractions
);

[[nodiscard]] std::vector<double> held2_lift_independent_fractions(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent_modified_fractions
);

[[nodiscard]] std::vector<double> held2_transform_modified_potentials(
    const Held2Coordinates& coordinates,
    const std::vector<double>& chemical_potentials
);

struct Held2PhysicalPhaseBlock {
    double helmholtz_over_rt = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    double pressure_pa = 0.0;
};

struct Held2StateEvaluation {
    std::vector<double> modified_fractions;
    std::vector<double> physical_amounts;
    double volume = 0.0;
    double objective = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    std::vector<double> modified_potentials;
    double pressure_stationarity_relative = 0.0;
};

[[nodiscard]] Held2StateEvaluation evaluate_held2_phase_block(
    const Held2Coordinates& coordinates,
    const std::vector<double>& independent_modified_fractions,
    double log_volume,
    double pressure_over_rt,
    double target_pressure_pa,
    const Held2PhysicalPhaseBlock& block
);

struct Held2Certificate {
    bool accepted = false;
    bool independent_evidence = false;
    double modified_balance_abs = 0.0;
    double ordinary_balance_inf_norm = 0.0;
    double phase_charge_inf_norm = 0.0;
    double modified_potential_gap = 0.0;
    double pressure_stationarity_inf_norm = 0.0;
    double reduced_kkt_inf_norm = 0.0;
    double enumeration_objective_gap = 0.0;
    double independent_modified_composition_count = 0.0;
};

struct Held2ManufacturedEvaluation {
    Held2Coordinates coordinates;
    std::vector<double> modified_feed;
    double phase_fraction = 0.0;
    std::array<std::vector<double>, 2> modified_phases;
    std::array<std::vector<double>, 2> physical_phases;
    std::array<double, 2> phase_charge_residuals{};
    std::vector<double> modified_balance;
    std::vector<double> ordinary_balance;
    std::vector<double> transformed_modified_potentials;
    std::array<std::array<double, 2>, 2> phase_gibbs_gradients{};
    std::array<std::array<double, 2>, 2> phase_modified_potentials{};
    double modified_potential_gap = 0.0;
    double pressure_stationarity_inf_norm = 0.0;
    double objective = 0.0;
    std::array<double, 4> gradient{};
    Held2Certificate certificate;
};

[[nodiscard]] Held2ManufacturedEvaluation evaluate_held2_manufactured(
    const std::vector<double>& charges,
    const std::vector<double>& physical_feed,
    const std::array<double, 4>& variables,
    const std::vector<double>& chemical_potentials
);

}  // namespace epcsaft_equilibrium
