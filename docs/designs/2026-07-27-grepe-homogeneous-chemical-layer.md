# GREPE Homogeneous Chemical Layer

**Status:** implemented public local-value design

**Date:** 2026-07-27

**Scope:** fixed-temperature, fixed-pressure homogeneous chemical equilibrium

## Decision and boundary

GREPE's reaction-compilation, structural-support, accessible-face, and local
chemical-certification concepts deepen the existing private homogeneous
module. `solve_provider_reaction` remains the sole Provider-backed chemical
optimization owner.

There is no parallel GREPE solver or sequential chemistry-then-flash result.
The public `chemical_equilibrium` operation is a typed cutover to this same
owner. The homogeneous module can emit
`LOCAL_EQUILIBRIUM`; phase pricing and global lower bounds are required before
any `SEARCH_STABLE` or `CERTIFIED_EPS_GLOBAL` claim.

This design supersedes the compiler-input and exact-zero deferrals in
`2026-07-21-private-reacting-phase-kernel.md`. Coupled-work deferrals and
installed Provider ownership in
`2026-07-24-reactive-prework-handoff.md` remain in force.

## Implemented claim

For one declared homogeneous fluid family, the homogeneous module:

1. validates and reduces a possibly redundant reaction set;
2. checks converted equilibrium-constant cycles before reduction;
3. verifies the invariant, charge, and reaction decomposition;
4. classifies species support using independently checked exact certificates;
5. removes only proved structural zeros and recompiles reactions on the
   accessible face;
6. solves the retained strict interior with the existing exact-Hessian Ipopt
   path; and
7. returns local, boundary, or unresolved evidence without a phase-stability
   claim.

The public value cutover adds generic source-reference transport through
Provider-owned neutral-reference metadata. It adds no sensitivity, promotion,
receipt, or authority change.

## Input and reaction compilation

The typed public input is an ordered true-species specification with positive molar
masses. Supplied reaction rows and equilibrium constants are unreduced source
data. Every constant record identifies its immutable source, orientation,
dimensionless converted value, common Provider reference convention, fixed
temperature and pressure, and conversion provenance.

Values either use the declared Provider convention or carry the explicit
source standard-state identity and activity-scale factors consumed by the
generic Provider neutral-reference transform. No chemistry-specific molality,
fugacity, or single-ion convention is inferred.

The caller supplies the conserved-component matrix. The compiler preserves
molar-mass conservation, treats charge as a separate homogeneous
electroneutrality equation, and fails if the rows do not span the invariant
space. It never infers conservation semantics from species names.

For the repository's reaction-row orientation, the retained transform obeys

\[
\nu^{\mathrm{sup}} = R_\nu \nu,\qquad
R_\nu k = k^{\mathrm{sup}}.
\]

Every left-null vector \(a\) of the supplied reaction rows must also satisfy

\[
a^T k^{\mathrm{sup}}=0.
\]

The compiled spaces must satisfy

\[
\operatorname{null}
\begin{bmatrix}
B\\z^T
\end{bmatrix}
=
\operatorname{row}(\nu).
\]

Only evidence required by the compiler tests, public result, or downstream solve is
retained: species maps, ranks, selected basis rows, reaction transforms,
compiled reactions and constants, support classifications, and the
minimum-norm reference vector. Validation-only residuals remain local.

## Exact homogeneous support

For

\[
\mathcal P_{\mathrm{str}}
=
\{n\ge0:Bn=b,\ z^Tn=0\},
\]

one deterministic HiGHS simplex solve maximizes each species amount. Its
floating-point solution and basis are candidate search evidence only.

The validator converts every finite input to its exact binary rational and
independently reconstructs one of:

- an exact feasible basic solution with \(n_i>0\), proving
  `proved_accessible`; or
- an exact dual-feasible upper bound of zero, proving
  `proved_structural_zero`.

Any failed reconstruction is `unresolved` and the species remains retained.
Only an exact zero certificate permits deletion. No positive zero tolerance,
LP transcript, primal vector, dual vector, or witness average becomes
persistent solver state. Numerical trace floors remain interior-solver devices
and never alter structural support.

## Accessible-face recompilation

Let retained species be accessible or unresolved. If \(N_{\mathrm{rem}}\)
contains removed-species columns, accessible reaction coefficients satisfy

\[
N_{\mathrm{rem}}^T\xi=0.
\]

For a null-space basis \(Z_{\mathrm{acc}}\),

\[
\nu^{\mathrm{acc}}
=
Z_{\mathrm{acc}}^T\nu P_{\mathrm{acc}},
\qquad
k^{\mathrm{acc}}
=
Z_{\mathrm{acc}}^Tk.
\]

This preserves reaction combinations whose removed-species coefficients
cancel; deleting species columns alone is invalid. The compiler restricts and
rank-reduces invariants, recomputes totals, verifies the accessible
null-space identity, reconstructs the reference vector, and records the
original-to-accessible maps.

## Solver and certification

The existing amount chart remains the only smooth species-coordinate owner.
The existing max-min initializer, exact-Hessian Ipopt TNLP, generic KKT polish,
Provider phase-block evaluator, and postsolve checks operate on the compiled
accessible system.

The complete manufactured ideal phase may restore certified structural zeros
in original species order. An installed Provider phase cannot omit components
unless its exact SDK supports that reduced topology and common reference
convention. Otherwise it returns `BOUNDARY_DIRECTION_UNRESOLVED` without
calling the Provider phase evaluator or substituting epsilon amounts.

`ChemicalSolveResult` uses these primary levels:

- `FEASIBLE_ONLY`;
- `LOCAL_EQUILIBRIUM`;
- `BOUNDARY_DIRECTION_UNRESOLVED`; and
- `INFEASIBLE_DECLARED_MODEL`.

`LOCAL_EQUILIBRIUM` requires applicable compiler, support, balance, charge,
pressure, affinity, KKT, complementarity, Provider-domain, trace, and local
curvature gates. It also requires a strict interior or an admitted
thermodynamic structural face. Missing exact support evidence adds
`SEARCH_ONLY_LP_SUPPORT` and cannot support the local-equilibrium level.

The module does not emit phase-family or globality levels.

## Ownership and evidence

- `chemical_equilibrium_compiler.cpp` owns reaction/invariant reduction and
  accessible-face recompilation.
- `chemical_equilibrium_support.cpp` owns HiGHS candidate search and exact
  certificate reconstruction.
- `chemical_equilibrium_solver.cpp` owns initialization, Ipopt, KKT polish,
  Provider evaluation, and local certification.
- `chemical_equilibrium_bindings.cpp` exposes one underscored native bridge
  consumed by the typed Python operation, plus derivative evidence seams.

No solver registry, chemistry database, backend selector, or second result
hierarchy is introduced.

Persistent scientific evidence covers:

1. redundant consistent reactions and inconsistent converted cycles;
2. exact structural-zero classification with accessible reaction
   reconstruction and manufactured structural-face restoration;
3. Belov–Aristova trace behavior through the same local solver; and
4. one installed-Provider manufactured thermodynamic-block integration case.

The installed case is nonpredictive and does not establish source-complete
chemistry.

## Deferred

- application chemistry databases;
- Provider physical-boundary certification without Provider-owned oracles;
- phase occurrence matrices, generation, master LPs, and pricing;
- simultaneous reactive multiphase equilibrium and final repricing;
- equilibrium sensitivities requiring absent Provider derivatives; and
- release, promotion, receipt, or authority changes.
