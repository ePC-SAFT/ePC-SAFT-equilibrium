#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace epcsaft_equilibrium {

struct Held2FarkasCertificate {
    std::string reason = "not_audited";
    std::vector<double> row_multipliers;
    double normalization = 0.0;
    double row_sign_violation_inf = std::numeric_limits<double>::infinity();
    double dual_feasibility_violation_inf =
        std::numeric_limits<double>::infinity();
    double contradiction_margin = -std::numeric_limits<double>::infinity();
    double contradiction_scale = 0.0;
    double contradiction_threshold = std::numeric_limits<double>::infinity();
    bool solver_ray_recovered_without_presolve = false;
    bool accepted = false;
};

struct Held2Step5KktCertificate {
    std::string reason = "not_audited";
    int major_iteration = -1;
    std::uint64_t start_ordinal = 0;
    double step4_upper_bound = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> step4_multipliers;
    std::vector<int> step4_active_cut_ids;
    std::vector<double> solver_variables;
    std::vector<double> audited_variables;
    std::vector<double> solver_lower_bound_multipliers;
    std::vector<double> solver_upper_bound_multipliers;
    std::vector<double> solver_constraint_multipliers;
    std::vector<double> lower_bound_multipliers;
    std::vector<double> upper_bound_multipliers;
    std::vector<double> constraint_multipliers;
    double physical_volume_lower = std::numeric_limits<double>::quiet_NaN();
    double physical_volume_upper = std::numeric_limits<double>::quiet_NaN();
    double primal_residual_inf = std::numeric_limits<double>::infinity();
    double dual_sign_violation_inf = std::numeric_limits<double>::infinity();
    double pullback_residual_inf = std::numeric_limits<double>::infinity();
    double pullback_scale = 0.0;
    double physical_composition_residual_inf =
        std::numeric_limits<double>::infinity();
    double physical_stationarity_residual_inf =
        std::numeric_limits<double>::infinity();
    double normalization_multiplier =
        std::numeric_limits<double>::quiet_NaN();
    double charge_multiplier = std::numeric_limits<double>::quiet_NaN();
    double pressure_residual = std::numeric_limits<double>::infinity();
    double pressure_derivative_log_volume =
        std::numeric_limits<double>::infinity();
    double stationarity_residual_inf = std::numeric_limits<double>::infinity();
    double complementarity_inf = std::numeric_limits<double>::infinity();
    std::size_t active_constraint_count = 0;
    std::size_t active_jacobian_rank = 0;
    bool same_major_iteration = false;
    bool step4_binding_valid = false;
    bool pressure_branch_valid = false;
    bool accepted = false;
};

[[nodiscard]] Held2FarkasCertificate audit_held2_farkas_certificate(
    const std::vector<std::vector<double>>& matrix,
    const std::vector<double>& row_lower,
    const std::vector<double>& row_upper,
    const std::vector<double>& column_lower,
    const std::vector<double>& column_upper,
    const std::vector<double>& row_ray
);

[[nodiscard]] std::string adjudicate_held2_farkas_status(
    bool solver_infeasible,
    const std::optional<Held2FarkasCertificate>& certificate
);

[[nodiscard]] Held2Step5KktCertificate audit_held2_step5_kkt(
    const std::vector<double>& variables,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    const std::vector<double>& objective_gradient,
    const std::vector<std::vector<double>>& constraints,
    const std::vector<double>& constraint_upper,
    const std::vector<double>& lower_bound_multipliers,
    const std::vector<double>& upper_bound_multipliers,
    const std::vector<double>& constraint_multipliers,
    double pullback_residual,
    double pullback_scale,
    double pressure_residual,
    bool same_major_iteration,
    bool step4_binding_valid,
    bool pressure_branch_valid
);

}  // namespace epcsaft_equilibrium
