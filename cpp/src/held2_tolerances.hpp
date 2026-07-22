#pragma once

#include <array>
#include <cmath>
#include <stdexcept>
#include <string_view>

namespace epcsaft_equilibrium {

enum class Held2ToleranceRelation {
    AbsAtMost,
    AtLeast,
    GreaterThan,
    LessThanNegative,
    SolverTarget,
};

struct Held2Tolerance {
    const char* name;
    const char* category;
    const char* failure_meaning;
    Held2ToleranceRelation relation;
    double atol;
    double rtol;
};

struct Held2ToleranceAudit {
    const Held2Tolerance* tolerance = nullptr;
    double residual = 0.0;
    double scale = 0.0;
    double threshold = 0.0;
    bool passed = false;
};

inline constexpr Held2Tolerance kHeld2ChartContact{
    "chart_contact", "representation", "chart contact exceeds numerical allowance",
    Held2ToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr Held2Tolerance kHeld2CompositionSum{
    "composition_sum", "representation", "composition normalization unresolved",
    Held2ToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr Held2Tolerance kHeld2ChargeBalance{
    "charge_balance", "representation", "scaled electroneutrality unresolved",
    Held2ToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr Held2Tolerance kHeld2ReconstructedIon{
    "reconstructed_ion", "representation", "eliminated-ion reconstruction is negative",
    Held2ToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr Held2Tolerance kHeld2BoundActivity{
    "bound_activity", "representation", "bound activity is numerically unresolved",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2RootPressure{
    "root_pressure", "root", "relative pressure closure failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2RootLogVolumeWidth{
    "root_log_volume_width", "root", "log-volume root bracket is unresolved",
    Held2ToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr Held2Tolerance kHeld2RootStationary{
    "root_stationary", "root", "pressure stationary-point refinement failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr Held2Tolerance kHeld2RootBoundary{
    "root_boundary", "root", "root is indistinguishable from Provider boundary",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2RootDuplicate{
    "root_duplicate", "root", "root identity is numerically unresolved",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2MechanicalMargin{
    "mechanical_margin", "topology", "mechanical stability is marginal",
    Held2ToleranceRelation::GreaterThan, 1.0e-6, 0.0,
};
inline constexpr Held2Tolerance kHeld2StableObjectiveTie{
    "stable_objective_tie", "topology", "stable-root objective ordering is unresolved",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 1.0e-9,
};
inline constexpr Held2Tolerance kHeld2TpdNegativeMargin{
    "tpd_negative_margin", "stage_i", "TPD is not materially negative",
    Held2ToleranceRelation::LessThanNegative, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2LpPrimal{
    "lp_primal", "stage_ii_lp", "upper LP primal feasibility failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-9, 1.0e-8,
};
inline constexpr Held2Tolerance kHeld2LpDual{
    "lp_dual", "stage_ii_lp", "upper LP dual feasibility failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-9, 1.0e-8,
};
inline constexpr Held2Tolerance kHeld2LpComplementarity{
    "lp_complementarity", "stage_ii_lp", "upper LP complementarity failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2LpActiveCut{
    "lp_active_cut", "stage_ii_lp_diagnostic", "cut is not diagnostically active",
    Held2ToleranceRelation::AbsAtMost, 1.0e-7, 0.0,
};
inline constexpr Held2Tolerance kHeld2Stage2Primal{
    "stage2_primal", "stage_ii_kkt", "physical primal feasibility failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2Stage2DualSign{
    "stage2_dual_sign", "stage_ii_kkt", "physical multiplier sign failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr Held2Tolerance kHeld2Stage2DualPullback{
    "stage2_dual_pullback", "stage_ii_kkt", "physical dual reconstruction unresolved",
    Held2ToleranceRelation::AbsAtMost, 1.0e-9, 1.0e-9,
};
inline constexpr Held2Tolerance kHeld2Stage2Stationarity{
    "stage2_stationarity", "stage_ii_kkt", "original-coordinate stationarity failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-7, 0.0,
};
inline constexpr Held2Tolerance kHeld2Stage2Complementarity{
    "stage2_complementarity", "stage_ii_kkt", "physical complementarity failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2Step6Gap{
    "step6_gap", "stage_ii_step6", "same-major upper/lower gap failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2Step6Gradient{
    "step6_gradient", "stage_ii_step6", "fixed-volume gradient agreement failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 1.0e-7,
};
inline constexpr Held2Tolerance kHeld2BasinDuplicateComposition{
    "basin_duplicate_composition", "numerical_identity", "composition is not a numerical duplicate",
    Held2ToleranceRelation::AbsAtMost, 1.0e-7, 0.0,
};
inline constexpr Held2Tolerance kHeld2BasinDuplicateLogVolume{
    "basin_duplicate_log_volume", "numerical_identity", "density branch is not a numerical duplicate",
    Held2ToleranceRelation::AbsAtMost, 1.0e-7, 0.0,
};
inline constexpr Held2Tolerance kHeld2CandidateDistinctComposition{
    "candidate_distinct_composition", "phase_identity", "candidate composition is not confidently distinct",
    Held2ToleranceRelation::GreaterThan, 1.0e-5, 0.0,
};
inline constexpr Held2Tolerance kHeld2CandidateDistinctLogVolume{
    "candidate_distinct_log_volume", "phase_identity", "candidate density is not confidently distinct",
    Held2ToleranceRelation::GreaterThan, 1.0e-5, 0.0,
};
inline constexpr Held2Tolerance kHeld2Stage3ModifiedBalance{
    "stage3_modified_balance", "stage_iii_physical", "modified material balance failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2Stage3ExplicitBalance{
    "stage3_explicit_balance", "stage_iii_physical", "explicit material balance failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2Stage3Charge{
    "stage3_charge", "stage_iii_physical", "scaled phase electroneutrality failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr Held2Tolerance kHeld2Stage3Pressure{
    "stage3_pressure", "stage_iii_physical", "relative phase pressure closure failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2Stage3Potential{
    "stage3_potential", "stage_iii_physical", "modified-potential equality failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 1.0e-7,
};
inline constexpr Held2Tolerance kHeld2Stage3Stationarity{
    "stage3_stationarity", "stage_iii_kkt", "Stage-III stationarity failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-7, 0.0,
};
inline constexpr Held2Tolerance kHeld2Stage3DualSign{
    "stage3_dual_sign", "stage_iii_kkt", "Stage-III multiplier sign failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr Held2Tolerance kHeld2Stage3Complementarity{
    "stage3_complementarity", "stage_iii_kkt", "Stage-III complementarity failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2Stage3FreeEnergyGap{
    "stage3_free_energy_gap", "stage_iii_physical", "reduced free-energy gap failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2PhaseActivity{
    "phase_activity", "phase_identity", "phase amount is not confidently active",
    Held2ToleranceRelation::GreaterThan, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2PhaseRetirementMargin{
    "phase_retirement_margin", "stage_iii_kkt", "phase retirement evidence is marginal",
    Held2ToleranceRelation::GreaterThan, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2PhaseMerge{
    "phase_merge", "phase_identity", "phases are not numerical duplicates",
    Held2ToleranceRelation::AbsAtMost, 1.0e-6, 0.0,
};
inline constexpr Held2Tolerance kHeld2PhaseDistinct{
    "phase_distinct", "phase_identity", "phase identity is not confidently distinct",
    Held2ToleranceRelation::GreaterThan, 1.0e-4, 0.0,
};
inline constexpr Held2Tolerance kHeld2IpoptTarget{
    "ipopt_target", "solver", "Ipopt target tolerance was not met",
    Held2ToleranceRelation::SolverTarget, 1.0e-10, 0.0,
};
inline constexpr Held2Tolerance kHeld2IpoptAcceptable{
    "ipopt_acceptable", "solver", "Ipopt acceptable tolerance was not met",
    Held2ToleranceRelation::SolverTarget, 1.0e-9, 0.0,
};
inline constexpr Held2Tolerance kHeld2IpoptConstraint{
    "ipopt_constraint", "solver", "Ipopt constraint target was not met",
    Held2ToleranceRelation::SolverTarget, 1.0e-10, 0.0,
};

inline constexpr std::array<const Held2Tolerance*, 44> kHeld2ToleranceContract{{
    &kHeld2ChartContact,
    &kHeld2CompositionSum,
    &kHeld2ChargeBalance,
    &kHeld2ReconstructedIon,
    &kHeld2BoundActivity,
    &kHeld2RootPressure,
    &kHeld2RootLogVolumeWidth,
    &kHeld2RootStationary,
    &kHeld2RootBoundary,
    &kHeld2RootDuplicate,
    &kHeld2MechanicalMargin,
    &kHeld2StableObjectiveTie,
    &kHeld2TpdNegativeMargin,
    &kHeld2LpPrimal,
    &kHeld2LpDual,
    &kHeld2LpComplementarity,
    &kHeld2LpActiveCut,
    &kHeld2Stage2Primal,
    &kHeld2Stage2DualSign,
    &kHeld2Stage2DualPullback,
    &kHeld2Stage2Stationarity,
    &kHeld2Stage2Complementarity,
    &kHeld2Step6Gap,
    &kHeld2Step6Gradient,
    &kHeld2BasinDuplicateComposition,
    &kHeld2BasinDuplicateLogVolume,
    &kHeld2CandidateDistinctComposition,
    &kHeld2CandidateDistinctLogVolume,
    &kHeld2Stage3ModifiedBalance,
    &kHeld2Stage3ExplicitBalance,
    &kHeld2Stage3Charge,
    &kHeld2Stage3Pressure,
    &kHeld2Stage3Potential,
    &kHeld2Stage3Stationarity,
    &kHeld2Stage3DualSign,
    &kHeld2Stage3Complementarity,
    &kHeld2Stage3FreeEnergyGap,
    &kHeld2PhaseActivity,
    &kHeld2PhaseRetirementMargin,
    &kHeld2PhaseMerge,
    &kHeld2PhaseDistinct,
    &kHeld2IpoptTarget,
    &kHeld2IpoptAcceptable,
    &kHeld2IpoptConstraint,
}};

[[nodiscard]] inline const char* held2_tolerance_relation_name(
    Held2ToleranceRelation relation
) {
    switch (relation) {
        case Held2ToleranceRelation::AbsAtMost: return "abs_at_most";
        case Held2ToleranceRelation::AtLeast: return "at_least";
        case Held2ToleranceRelation::GreaterThan: return "greater_than";
        case Held2ToleranceRelation::LessThanNegative: return "less_than_negative";
        case Held2ToleranceRelation::SolverTarget: return "solver_target";
    }
    throw std::logic_error("unknown HELD2 tolerance relation");
}

[[nodiscard]] inline double held2_tolerance_threshold(
    const Held2Tolerance& tolerance,
    double scale
) {
    if (!std::isfinite(scale) || scale < 0.0) {
        throw std::invalid_argument("HELD2 tolerance scale must be finite and nonnegative");
    }
    return tolerance.atol + tolerance.rtol * scale;
}

[[nodiscard]] inline Held2ToleranceAudit audit_held2_tolerance(
    const Held2Tolerance& tolerance,
    double residual,
    double scale = 0.0
) {
    if (!std::isfinite(residual)) {
        throw std::invalid_argument("HELD2 tolerance residual must be finite");
    }
    const double threshold = held2_tolerance_threshold(tolerance, scale);
    bool passed = false;
    switch (tolerance.relation) {
        case Held2ToleranceRelation::AbsAtMost:
        case Held2ToleranceRelation::SolverTarget:
            passed = std::abs(residual) <= threshold;
            break;
        case Held2ToleranceRelation::AtLeast:
            passed = residual >= threshold;
            break;
        case Held2ToleranceRelation::GreaterThan:
            passed = residual > threshold;
            break;
        case Held2ToleranceRelation::LessThanNegative:
            passed = residual < -threshold;
            break;
    }
    return {&tolerance, residual, scale, threshold, passed};
}

[[nodiscard]] inline const Held2Tolerance& held2_tolerance_by_name(
    std::string_view name
) {
    for (const Held2Tolerance* tolerance : kHeld2ToleranceContract) {
        if (name == tolerance->name) {
            return *tolerance;
        }
    }
    throw std::invalid_argument("unknown HELD2 tolerance name");
}

}  // namespace epcsaft_equilibrium
