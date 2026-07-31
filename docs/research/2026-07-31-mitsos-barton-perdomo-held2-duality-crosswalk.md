# Mitsos-Barton and Perdomo HELD2 duality crosswalk

Date: 2026-07-31
Scope: primary-source theorem audit against the current HELD2 evidence contract.

Sources read and visually checked:

- Alexander Mitsos and Paul I. Barton, “A dual extremum principle in
  thermodynamics,” *AIChE Journal* 53(8), 2131-2147 (2007), DOI
  `10.1002/aic.11209`; local source
  `/home/tnnrpolley21/Zotero/storage/LU5Q79R8/Mitsos and Barton - 2007 - A dual extremum principle in thermodynamics.pdf`.
- Felipe A. Perdomo et al., “Phase stability criteria and fluid-phase
  equilibria in strong-electrolyte systems,” *Computers & Chemical
  Engineering* 194, 108977 (2025), DOI
  `10.1016/j.compchemeng.2024.108977`; local source
  `/home/tnnrpolley21/Zotero/storage/GVEMXCZD/Perdomo et al. - 2025 - Phase stability criteria and fluid-phase equilibria in strong-electrolyte systems.pdf`.

All 17 Mitsos-Barton pages and all 21 Perdomo pages were read. The rendered
theorem, algorithm, proof, and appendix pages were checked for equations,
qualifiers, and continuation across page boundaries. `[Fact]` is stated by a
source; `[Inference]` maps those facts to this repository.

## Bottom line

[Fact] Mitsos and Barton provide the missing proof foundation for Perdomo's
Appendix A tangent-plane argument. Their Theorem 5 establishes the stability
criterion, and their Lemma 12 supplies the line-segment argument that Perdomo
adapts to modified electrolyte coordinates. Perdomo Appendix A proves only
Theorem 3.1, not the dual-extremum Theorem 3.2.

[Fact] Mitsos and Barton's stronger dual extremum and finite-convergence
results require an actual dual optimum and global solutions of the nonconvex
inner problem. Perdomo explicitly says its tunnelling implementation does not
guarantee the global TPD minimum and recommends a deterministic global solver
as future strengthening.

[Inference] The paper closes a proof-provenance gap but does not upgrade the
current finite HELD2 controller to a global phase-stability proof. The current
`globality_certificate="not_guaranteed"` and
`root_completeness="not_proven"` contracts are required. Step 4, Step 5,
Step 8, and Step 9 need scope labels that prevent their finite or local
quantities from being read as Mitsos-Barton theorem-level bounds.

## 1. What Mitsos and Barton prove

[Fact] Theorem 5, PDF p. 5 (printed p. 2135), Eqs. (7)-(8), states the Gibbs
tangent-plane stability criterion for continuous Gibbs energy, positive
overall composition, interior phase compositions, nonzero phase amounts, and
material balance. The surrounding text on PDF pp. 5-6 explicitly says the
result does not require differentiability, all species to occur in every
phase, or a finite number of phases.

[Fact] Theorem 6, PDF pp. 6-7 (printed pp. 2136-2137), Eq. (9), is conditional
on a solution of the dual problem and an inner minimizer satisfying
`L(x*, lambda*) = q(lambda*) = G^d`. Under that condition, the resulting
hyperplane supports the complete Gibbs surface, the feed lies in the convex
hull of all common points, and the stable phases are common points. A local
stationary point or one supporting contact is not the theorem's hypothesis.

[Fact] Algorithm 1, PDF p. 7, Eqs. (10)-(11), repeatedly solves an upper
bounding linear semi-infinite relaxation and a nonconvex lower problem to
global optimality. Theorem 17, PDF pp. 16-17 (printed pp. 2146-2147), proves
finite termination from those exact bounding operations and continuity of the
dual function. It does not prove finite termination for a capped local
multistart substitute.

## 2. What Perdomo imports

[Fact] Perdomo eliminates one charged species and uses `C-2` independent
modified mole fractions. Theorem 3.1, PDF p. 5, Eqs. (35)-(40), states the
modified Gibbs tangent-plane criterion. Appendix A, PDF pp. 17-18,
Eqs. (A.1)-(A.9), follows Mitsos and Barton's direct proof and explicitly
adapts their Lemma 12 to a line segment in the modified composition domain.
This supplies the necessary and sufficient stability proof in the reduced
electroneutral coordinates.

[Fact] Perdomo Theorem 3.2, PDF p. 6, Eqs. (43)-(48), states the electrolyte
dual-extremum analogue but says its proof is omitted because it is similar to
Mitsos and Barton. Appendix A does not prove Theorem 3.2. Appendix B,
PDF pp. 18-19, Eqs. (B.1)-(B.6), states the Helmholtz analogue used by HELD2,
but likewise does not reproduce the missing global-inner-minimization proof.

[Fact] Perdomo Section 4.2, PDF p. 8, Eq. (63), calls for a global TPD
minimization and then states that the implemented tunnelling method does not
guarantee identification of the global minimum. Multiple starts only increase
reliability. Section 6, PDF pp. 16-17, proposes a reliable deterministic global
solver as future work.

## 3. Repository consequences

[Inference] A finite, independently recertified negative TPD state remains a
valid one-sided instability witness: finding one feasible state below the feed
tangent is enough to disprove homogeneous stability. Completing a finite
search without finding one is not the converse proof.

[Inference] Step 4's HiGHS result is the optimum of the current finite-cut LP.
It is not, by itself, a certified dual optimum over the complete Gibbs or
Helmholtz surface.

[Inference] Step 5's serialized `lower_value` is the best independently
KKT-certified local Problem-(65) value found by the declared finite start
schedule. It is not a certified Mitsos-Barton dual value `q(lambda)` unless
the inner global minimum has separately been proved.

[Inference] Step 8's Farkas certificate proves infeasibility only for the
current candidate-restricted linear master used to initialize Problem (67).
It does not prove that no physical phase split exists outside those candidate
neighborhoods.

[Inference] Step 9's `free_energy_gap` is the signed difference between the
finite-cut Step-4 value and the restricted-candidate Step-8 total Gibbs value.
It is Perdomo's finite algorithmic convergence test, not the theorem-level
`UBD-LBD` global duality gap of Mitsos-Barton Algorithm 1.

[Inference] The deterministic pressure-envelope implementation detects and
refines roots under finite scan and subdivision controls. Its result must not
be described as exhaustive root enumeration while
`root_completeness="not_proven"` remains the evidence status.

## 4. Applied contract changes

The implementation and specification preserve these independent statements:

1. `negative_witness_found`: a recertified feasible negative TPD state exists;
2. `finite_search_completed_no_negative_witness`: every declared finite search
   completed without such a witness;
3. `global_stability_proved`: reserved for a future route that satisfies the
   global inner-minimization and dual-equality hypotheses; and
4. `not_adjudicated`: the required evidence was not completed.

Current HELD2 results never emit item 3. Structured diagnostics identify the
scope of Step-4 bounds, Step-5 local values, the Step-8 restricted master, and
the Step-9 finite-candidate gap. These labels change no solver decision or
numerical tolerance; they prevent a stronger scientific interpretation than
the primary sources support.
