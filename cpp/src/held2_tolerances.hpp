#pragma once

#include <cmath>
#include <stdexcept>

namespace epcsaft_equilibrium {

inline constexpr double kHeld2Problem67Radius = 1.0e-3;

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

inline constexpr Held2Tolerance kHeld2PolytopeFeasibility{
    "polytope_feasibility", "representation",
    "modified composition is outside the complete Step-1 polytope",
    Held2ToleranceRelation::AbsAtMost, 1.0e-12, 0.0,
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
inline constexpr Held2Tolerance kHeld2JointVolumeConsistency{
    "joint_volume_consistency", "representation",
    "Provider state does not match the requested joint-search molar volume",
    Held2ToleranceRelation::AbsAtMost, 1.0e-10, 0.0,
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
inline constexpr Held2Tolerance kHeld2TpdReferenceZero{
    "tpd_reference_zero", "stage_i", "TPD is nonzero at its tangent point",
    Held2ToleranceRelation::AbsAtMost, 1.0e-12, 0.0,
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
inline constexpr Held2Tolerance kHeld2Stage2KktPrimal{
    "stage2_kkt_primal", "stage_ii_kkt",
    "Step-5 original-coordinate primal feasibility failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2Stage2KktDualSign{
    "stage2_kkt_dual_sign", "stage_ii_kkt",
    "Step-5 bound or inequality multiplier sign failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr Held2Tolerance kHeld2Stage2KktPullback{
    "stage2_kkt_pullback", "stage_ii_kkt",
    "Step-5 physical-to-modified dual reconstruction failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-9, 1.0e-9,
};
inline constexpr Held2Tolerance kHeld2Stage2KktStationarity{
    "stage2_kkt_stationarity", "stage_ii_kkt",
    "Step-5 original-coordinate stationarity failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-7, 0.0,
};
inline constexpr Held2Tolerance kHeld2Stage2KktComplementarity{
    "stage2_kkt_complementarity", "stage_ii_kkt",
    "Step-5 bound or inequality complementarity failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2Stage2KktRankPivot{
    "stage2_kkt_rank_pivot", "stage_ii_kkt",
    "Step-5 active-constraint Jacobian is rank deficient",
    Held2ToleranceRelation::GreaterThan, 1.0e-12, 0.0,
};
inline constexpr Held2Tolerance kHeld2FarkasRowSign{
    "farkas_row_sign", "stage_iii_lp_certificate",
    "Step-8 Farkas row multiplier has the wrong bound sign",
    Held2ToleranceRelation::AbsAtMost, 1.0e-10, 0.0,
};
inline constexpr Held2Tolerance kHeld2FarkasDual{
    "farkas_dual", "stage_iii_lp_certificate",
    "Step-8 Farkas column dual is incompatible with variable bounds",
    Held2ToleranceRelation::AbsAtMost, 1.0e-10, 0.0,
};
inline constexpr Held2Tolerance kHeld2FarkasContradiction{
    "farkas_contradiction", "stage_iii_lp_certificate",
    "Step-8 Farkas contradiction margin is not strict",
    Held2ToleranceRelation::GreaterThan, 1.0e-9, 1.0e-9,
};
inline constexpr Held2Tolerance kHeld2Step6Gap{
    "step6_gap", "stage_ii_step6", "same-major upper/lower gap failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2PaperStep6Gap{
    "paper_step6_gap", "stage_ii_paper",
    "Perdomo Eq. (66) upper/lower agreement failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-2, 0.0,
};
inline constexpr Held2Tolerance kHeld2PaperStep6Derivative{
    "paper_step6_derivative", "stage_ii_paper",
    "Perdomo Eq. (66) derivative agreement failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.5,
};
inline constexpr Held2Tolerance kHeld2PaperStep6PackingDistinct{
    "paper_step6_packing_distinct", "phase_identity",
    "Perdomo Eq. (66) packing fractions are not distinct",
    Held2ToleranceRelation::AtLeast, 1.0e-3, 0.0,
};
inline constexpr Held2Tolerance kHeld2PaperStep6CompositionDistinct{
    "paper_step6_composition_distinct", "phase_identity",
    "Perdomo Eq. (66) modified compositions are not distinct",
    Held2ToleranceRelation::AtLeast, 1.0e-3, 0.0,
};
inline constexpr Held2Tolerance kHeld2MRepresentationEquivalent{
    "m_representation_equivalent", "phase_identity",
    "Step-5 terminal differs in modified composition or relative molar volume",
    Held2ToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
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
inline constexpr Held2Tolerance kHeld2PaperFreeEnergyGap{
    "paper_free_energy_gap", "stage_iii_paper",
    "Perdomo Eq. (68) free-energy convergence failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-4, 0.0,
};
inline constexpr Held2Tolerance kHeld2PaperPotentialRatio{
    "paper_potential_ratio", "stage_iii_paper",
    "Perdomo Eq. (69) modified-potential convergence failed",
    Held2ToleranceRelation::AbsAtMost, 1.0e-3, 0.0,
};
inline constexpr Held2Tolerance kHeld2PhaseRetirementMargin{
    "phase_retirement_margin", "stage_iii_kkt", "phase retirement evidence is marginal",
    Held2ToleranceRelation::GreaterThan, 1.0e-8, 0.0,
};
inline constexpr Held2Tolerance kHeld2PhaseMerge{
    "phase_merge", "phase_identity", "phases are not numerical duplicates",
    Held2ToleranceRelation::AbsAtMost, 1.0e-5, 0.0,
};
inline constexpr Held2Tolerance kHeld2IpoptTarget{
    "ipopt_target", "solver", "Ipopt target tolerance was not met",
    Held2ToleranceRelation::SolverTarget, 1.0e-10, 0.0,
};
inline constexpr Held2Tolerance kHeld2IpoptConstraint{
    "ipopt_constraint", "solver", "Ipopt constraint target was not met",
    Held2ToleranceRelation::SolverTarget, 1.0e-10, 0.0,
};
[[nodiscard]] inline Held2ToleranceAudit audit_held2_tolerance(
    const Held2Tolerance& tolerance,
    double residual,
    double scale = 0.0
) {
    if (!std::isfinite(residual) || !std::isfinite(scale) || scale < 0.0) {
        throw std::invalid_argument(
            "HELD2 tolerance residual must be finite and scale must be "
            "finite and nonnegative"
        );
    }
    const double threshold = tolerance.atol + tolerance.rtol * scale;
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

}  // namespace epcsaft_equilibrium
