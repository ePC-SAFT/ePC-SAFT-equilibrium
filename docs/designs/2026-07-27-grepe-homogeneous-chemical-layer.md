# GREPE Homogeneous Chemical Layer

**Status:** implemented private non-production design

**Date:** 2026-07-27

**Scope:** private fixed-temperature, fixed-pressure homogeneous chemical
equilibrium

## Decision

Migrate the GREPE chemical compiler, structural-support, accessible-face, and
chemical-certification concepts into the existing private homogeneous
chemical-equilibrium module. Keep `solve_provider_reaction` as the sole
Provider-backed chemical optimization owner.

This design deepens the current private module. It does not create a parallel
GREPE solver, a public chemical-equilibrium route, or a sequential
chemistry-then-flash result.

The implementation may emit a GREPE-compatible local chemical-equilibrium
level. It cannot emit `SEARCH_STABLE` or `CERTIFIED_EPS_GLOBAL` because
phase-family pricing and rigorous global lower bounds remain outside the
homogeneous module.

This design supersedes the compiler-input and exact-zero deferrals in
`2026-07-21-private-reacting-phase-kernel.md`. It does not supersede the
coupled-work deferrals or the installed
Provider ownership rules in
`2026-07-24-reactive-prework-handoff.md`.

## Implementation outcome

The private module now implements the selected design through the existing
compiler and solve owners:

- redundant supplied rows are reduced only after mass, balance, charge,
  provenance, and converted-cycle validation;
- deterministic HiGHS reachability candidates are independently reconstructed
  and checked with exact Boost multiprecision binary rationals;
- only exact dual zero certificates can delete a species;
- accessible reaction combinations are rebuilt before reference
  reconstruction;
- manufactured ideal structural faces restore exact zeros in original order;
  and
- Provider structural faces return `BOUNDARY_DIRECTION_UNRESOLVED` before any
  Provider callback.

Accepted homogeneous interiors and supported manufactured structural faces
emit only `LOCAL_EQUILIBRIUM`. The implementation adds no public export,
source-reference transport, phase-stability level, coupled solve, sensitivity,
promotion, or authority change.

## Scientific claim

For one declared homogeneous fluid family, the private module will:

1. compile a possibly redundant supplied reaction set into one independent
   reaction basis;
2. check converted equilibrium-constant cycles before rank reduction;
3. construct and verify a complete invariant, charge, and reaction
   decomposition;
4. classify homogeneous species support using validated linear certificates;
5. remove only proved structural zeros and recompile the reaction space on the
   accessible face;
6. solve the remaining strictly positive fixed-\(T,P\) equilibrium with the
   existing exact-Hessian Ipopt path; and
7. return explicit compilation, support, local-equilibrium, boundary, and
   unresolved evidence.

The resulting state is a homogeneous local equilibrium. It is suitable as a
candidate or seed for future GREPE phase generation, but it is not a phase
stability result.

## Approaches considered

### 1. Deepen the existing private module

Extend `ReactionSystemInput`, `CompiledReactionSystem`, and
`ChemicalSolveResult` in place. Put structural LP and exact-certificate logic
behind one private support interface. Continue through the existing amount
chart and `solve_provider_reaction`.

This is the selected approach. It preserves one owner and lets existing Belov,
Provider-block, KKT, and trace evidence exercise the migrated compiler.

### 2. Add a parallel GREPE chemical compiler and solver

This would preserve the old private input while adding a second reaction
representation and solve path. It is rejected because the two paths would own
the same reaction reduction, reference reconstruction, amount coordinates, and
certification concepts.

### 3. Implement the complete multiphase GREPE controller now

This would add phase-family occurrence matrices, column generation, pricing,
the simultaneous all-phase NLP, repricing, and global certificates. It is
rejected for this slice because the user requested migration into the current
chemical capability and coupled phase equilibrium remains a separately
governed outcome.

## Input and compilation contract

The private input remains an ordered true-species specification. It adds
positive molar masses and treats the supplied reaction rows and equilibrium
constants as the unreduced source set.

Each equilibrium-constant record must identify:

- its immutable source;
- its reaction orientation;
- its dimensionless converted value;
- the common Provider chemical-reference convention;
- its fixed temperature and pressure; and
- the conversion provenance.

The compiler accepts values only after conversion into the declared Provider
convention. The current slice does not invent a generic molality, fugacity, or
single-ion conversion. A record in another convention fails before cycle
checking. The dormant Provider neutral-reference transport remains deferred
until one source-complete, non-MEA installed-artifact consumer exists.

The caller may supply a formula or conserved-component matrix. The compiler
rank-reduces it while retaining the molar-mass row explicitly and retaining
charge as a separate homogeneous electroneutrality equation. If no supplied
row set spans the full invariant space, compilation fails. The compiler does
not silently invent missing chemical conservation semantics from species
names.

For the repository's reaction-row orientation, the compiler records a matrix
\(R_\nu\) satisfying

\[
\nu^{\mathrm{sup}} = R_\nu \nu
\]

and reconstructs the independent constants from

\[
R_\nu k = k^{\mathrm{sup}}.
\]

Before selecting that basis, every vector \(a\) in the left null space of the
supplied reaction rows must satisfy

\[
a^T k^{\mathrm{sup}}=0.
\]

The compiled invariant and charge rows must satisfy

\[
\operatorname{null}
\begin{bmatrix}
B\\z^T
\end{bmatrix}
=
\operatorname{row}(\nu).
\]

All raw matrices, retained matrices, transformations, ranks, conditioning
diagnostics, cycle residuals, and reconstruction residuals remain in the
compiled evidence.

## Homogeneous structural support

One homogeneous species occurrence exists for every declared species. Define

\[
\mathcal P_{\mathrm{str}}
=
\{n\ge0:Bn=b,\ z^Tn=0\}.
\]

A HiGHS Phase-I candidate establishes numerical feasibility. One reachability
LP then maximizes each scaled species amount over
\(\mathcal P_{\mathrm{str}}\).

Floating-point LP output is search evidence only. Exact deletion requires a
validated certificate:

- an exact primal-feasible witness with positive objective proves
  `proved_accessible`;
- an exact dual-feasible upper bound of zero proves
  `proved_structural_zero`;
- any other outcome is `unresolved` and remains retained.

Validation uses exact rational arithmetic over the supplied finite decimal or
binary-rational values. HiGHS supplies candidate active sets and primal/dual
vectors; the validator reconstructs and checks them independently. Failure to
reconstruct an exact certificate is not an error and cannot delete a species.

The witness-average state must satisfy balances and charge and be positive on
every proved-accessible species. It becomes the preferred max-min
initialization evidence.

No positive zero tolerance is introduced. Numerical trace floors remain solver
devices and never change structural support.

## Accessible-face recompilation

Let retained species be those classified as accessible or unresolved. If
\(N_{\mathrm{rem}}\) contains the removed-species columns of the retained
reaction-row matrix, an accessible reaction combination has coefficients
\(\xi\) satisfying

\[
N_{\mathrm{rem}}^T\xi=0.
\]

For a basis \(Z_{\mathrm{acc}}\) of this null space, the row-oriented
accessible reaction system is

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
cancel. Simple species-column deletion is forbidden.

The compiler then:

1. restricts and rank-reduces the invariant and charge rows;
2. preserves the mass row;
3. recomputes feed totals;
4. verifies the accessible null-space identity;
5. reconstructs the accessible chemical-reference vector; and
6. records original-to-accessible species maps.

## Solver integration

The existing amount chart remains the only smooth species-coordinate owner.
Every retained nonstructural species is strictly positive. The existing
max-min initializer, exact-Hessian Ipopt TNLP, generic KKT polish, Provider
phase-block evaluator, and postsolve checks operate on the compiled accessible
system.

For the manufactured complete ideal phase, proved structural zeros can be
removed before optimization and restored as exact zeros in the result because
the closed ideal term and its accessible-face restriction are known.

For an installed Provider phase, a removed component can be omitted only when
the exact installed SDK explicitly supports that reduced component topology
and common reference convention. Otherwise the solve stops with
`BOUNDARY_DIRECTION_UNRESOLVED`; it must not pass an artificial epsilon amount
to Provider.

The generic Newton KKT polish remains a strictly interior correction. It is
not a boundary solver and cannot change structural classifications.

## Boundary and domain certification

The homogeneous result distinguishes:

- `strict_interior`: every retained species is positive and every Provider
  domain inequality is inactive;
- `structural_face`: exact zeros have validated compiler certificates and the
  admitted thermodynamic implementation supports the reduced face;
- `boundary_direction_unresolved`: an exact-zero or physical-domain boundary
  requires a directional oracle or rigorous separation bound that is not
  available; and
- `rejected_boundary`: an exact feasible activating direction or negative
  bounded directional value disproves local boundary equilibrium.

For a smooth polyhedral face, the module may use a bounded HiGHS separation LP
and exact dual validation. A Provider physical-domain boundary is never
certified without Provider-owned directional feasibility and derivative
oracles.

## Result contract

`ChemicalSolveResult` gains one primary chemical-certification level:

- `FEASIBLE_ONLY`;
- `LOCAL_EQUILIBRIUM`;
- `BOUNDARY_DIRECTION_UNRESOLVED`; or
- `INFEASIBLE_DECLARED_MODEL`.

It also records support qualifiers, including `SEARCH_ONLY_LP_SUPPORT` when LP
support evidence lacks exact validation.

`LOCAL_EQUILIBRIUM` requires all applicable current gates:

- compiler and support validity;
- balance and exact charge;
- pressure;
- reaction affinity;
- original-coordinate KKT stationarity;
- complementarity for active smooth inequalities;
- Provider domain;
- trace classification;
- reduced-Hessian local-minimum evidence; and
- either a strict interior or a fully validated structural face.

The homogeneous module never emits `NEGATIVE_PHASE_FOUND`, `SEARCH_STABLE`,
`CERTIFIED_EPS_GLOBAL`, or `RELATIVE_INTERIOR_UNVERIFIED`; those levels require
phase-family master and pricing information.

## File ownership

- `cpp/src/chemical_equilibrium.hpp` owns the private input, compiled evidence,
  support evidence, and result types.
- `cpp/src/chemical_equilibrium.cpp` owns reaction/invariant compilation,
  accessible-face recompilation, and the amount chart.
- `cpp/src/chemical_equilibrium_support.cpp` owns the homogeneous HiGHS
  reachability problems and exact certificate validation.
- `cpp/src/chemical_equilibrium_solver.cpp` continues to own initialization,
  Ipopt, KKT polish, Provider evaluation, and local certification.
- `cpp/src/chemical_equilibrium_bindings.cpp` keeps the underscored test seam
  and serializes the deeper evidence. No public Python export is added.
- `tests/test_chemical_equilibrium.py` remains the compact persistent
  scientific test surface.

No generic solver registry, reaction database, chemistry-named implementation,
or second result hierarchy is introduced.

## Testing

Persistent evidence is limited to three parameterized scientific families:

1. **Reaction compilation:** redundant and permuted supplied reactions,
   converted-cycle consistency, mass/charge conservation, independent-basis
   reconstruction, and fail-closed inconsistent cycles.
2. **Structural support:** neutral structural zero, jointly accessible ionic
   species, exact zero versus unresolved classification, witness averaging,
   and a reaction combination preserved only by accessible-face
   cancellation.
3. **End-to-end local equilibrium:** the Belov–Aristova trace benchmark and
   one compact manufactured structural-face case, checking composition,
   invariance, balances, pressure, affinities, KKT, local curvature, support
   evidence, and the primary certification level.

The installed Provider manufactured case remains a thermodynamic-block and
domain integration test. It does not become source-complete predictive
chemistry.

Every production behavior is introduced red-before-green. Numerical tolerances
remain tied to the existing physical residual scales or to independently
validated exact certificates; they are not tuned to current solver output.

## Explicitly deferred

- source-standard-state transport without a source-complete non-MEA installed
  consumer;
- application chemistry databases;
- phase-family occurrence matrices;
- phase generation, master LPs, pricing, and final state-matched repricing;
- a simultaneous reactive multiphase NLP;
- rigorous Provider physical-boundary certification without Provider oracles;
- equilibrium sensitivities requiring absent Provider derivatives;
- public exports, release, promotion, receipt, or authority changes.

## Acceptance

The migration is complete for the homogeneous chemical layer when:

1. the compiler accepts redundant consistent reactions and rejects inconsistent
   converted cycles before reduction;
2. the mass, invariant, charge, and reaction spaces are complete and recorded;
3. structural deletion occurs only with exact validated zero certificates;
4. accessible-face reaction combinations are reconstructed correctly;
5. Belov and the retained Provider manufactured interior cases pass through
   the same existing solver owner;
6. structural or Provider boundary cases either carry complete evidence or
   fail closed with the exact unresolved level;
7. no coupled-phase or global-stability claim is emitted; and
8. focused native, static, full-suite, and cleanup verification pass.
