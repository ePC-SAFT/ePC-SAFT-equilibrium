# GREPE Reactive-LLE Implementation Plan

Status: approved documentation plan; implementation not started

Date: 2026-07-29

**Goal:** Add one private, typed `equilibrate` facade that preserves the
existing homogeneous chemical and nonreactive phase-equilibrium routes and
implements full GREPE discovery and final simultaneous optimization for
fixed-\(T,P\), at-most-two-liquid neutral reactive equilibrium.

**Scientific authority:**

- `docs/designs/2026-07-28-grepe-reactive-phase-equilibrium.md`
- `docs/designs/2026-07-24-held2-paper-algorithm.md`
- `docs/designs/2026-07-29-held2-publication-algorithm.md`
- `docs/designs/2026-07-27-grepe-homogeneous-chemical-layer.md`
- `docs/designs/2026-07-29-reactive-lle-benchmark-review.md`

**Architecture:** Existing native HELD2 and homogeneous chemical controllers
remain single owners. A private typed facade delegates limiting problems to
them. The new coupled controller composes their lower-level chemistry,
thermodynamic, and deterministic-search mechanisms but performs one
simultaneous reactive phase calculation. No controller or EOS implementation
is copied.

## Non-negotiable contracts

- The attached checkout and implementation branch own every edit, build,
  review, commit, and PR operation.
- The installed Provider public SDK owns Helmholtz values, nonlinear tensors,
  domains, packing fractions, component identity, and parameter fingerprints.
- The first coupled scope is one Provider liquid model, one or two discovered
  liquid phases, fixed \(T,P\), and the shared compiler-accessible species set.
- Callers supply species feed amounts, source reactions and constants, physical
  conserved-component rows, and immutable provenance. Compilers own every
  reduced basis and coordinate chart.
- Phase amounts can reach zero. Species and phases are never retired because
  they are small.
- Every production derivative is analytic; the final Ipopt problem uses an
  exact Lagrangian Hessian.
- One derived nondimensionalization contract is used by all stages.
- Every finite search retains
  `globality_certificate="not_guaranteed"`.
- Provider failures, uncertified optimizer terminals, stale evidence,
  unresolved topology, and resource exhaustion fail closed.
- Budgets count Provider callbacks and solver work. Wall time is evidence, not
  a correctness assertion.
- `tp_flash` and `chemical_equilibrium` remain unchanged public fallbacks until
  a later promotion decision.

## Task 0: freeze the exact starting behavior

- [ ] Record the Equilibrium commit, installed Provider version/fingerprint,
  compiler configuration, and relevant catalog fingerprints.
- [ ] Run the native manufactured HELD2 Steps 1--10 check.
- [ ] Run Perdomo regression-region tests, including the Table-3 point.
- [ ] Run the complete Perdomo numerical matrix.
- [ ] Run native/Python diagnostic parity.
- [ ] Run the complete chemical compiler, support, solver, sensitivity, and
  public-value suites.
- [ ] Record Provider calls and wall times for the fast Perdomo and installed
  water-ionization sentinels.
- [ ] Store no case-specific values in production code.

Checkpoint: the checkout is clean after the evidence record, and no scientific
behavior has changed.

## Task 1: close HELD2 certificate gaps in isolation

### Step-5 physical-coordinate KKT evidence

- [ ] Reconstruct the full physical-coordinate terminal from Problem (65).
- [ ] Independently audit primal constraints, physical stationarity,
  lower/upper bound multiplier signs, complementarity, and modified-to-physical
  dual pullback.
- [ ] Add manufactured pass, sign-failure, complementarity-failure, and
  reconstruction-failure cases.
- [ ] Make missing evidence indeterminate; do not alter the objective,
  multistart schedule, stop condition, or tolerances incidentally.

### Step-8 Farkas evidence

- [ ] Retain the exact perspective-form feasibility LP.
- [ ] Validate a solver-proposed infeasibility certificate independently
  before returning `certified_infeasible`.
- [ ] Distinguish feasible, validated infeasible, and unresolved solver
  outcomes.
- [ ] Add feasible, analytically infeasible, invalid-certificate, and numerical
  ambiguity tests.

Checkpoint: Task-0 numerical outputs and Provider-call budgets remain in their
frozen regions. This task is its own reviewable change.

## Task 2: extract the shared deterministic joint-search mechanism

- [ ] Identify the smallest HELD2 Step-2 mechanics that are independent of TPD
  and modified electrolyte coordinates.
- [ ] Extract only normalized-box exploration, immutable candidate ordering,
  Provider-call accounting, early-stop plumbing, and deterministic termination.
- [ ] Keep coordinate mapping, physical feasibility, objective evaluation, and
  independent witness recertification in HELD2.
- [ ] Add a formulation-neutral narrow-basin manufactured check.
- [ ] Prove the HELD2 native payload, Perdomo results, and Provider-call budget
  are unchanged.
- [ ] Delete superseded nested or duplicated search plumbing.

Checkpoint: no GREPE production controller exists yet, and HELD2 still has one
scientific Step-2 owner.

## Task 3: add private typed facade and coupled value contracts

- [ ] Define private typed problem variants for:
  - constrained homogeneous reacting equilibrium;
  - nonreactive phase discovery; and
  - reactive phase discovery.
- [ ] Define one discriminated internal result with independent input, solver,
  numerical, physical, search, predictive, and globality axes.
- [ ] Route the first two variants to the existing native owners without
  copying callbacks or result construction.
- [ ] Normalize timing fields and prove exact or tolerance-owned payload parity
  against the existing public operations.
- [ ] Add no public symbol and no compatibility alias.

Checkpoint: the facade adds no new physics and both frozen public operations
remain unchanged.

## Task 4: implement GREPE Stage 0 and scaling

- [ ] Reuse the existing reaction compiler, standard-state conversion,
  conserved-component audit, support classifier, accessible-face compiler, and
  Provider phase-block evaluator.
- [ ] Validate same-Provider liquid-family scope, fixed \(T,P\), ordered
  species, feed amounts, physical conserved rows, reaction cycles, mass,
  charge, null-space identity, fingerprint, and derivative availability.
- [ ] Derive the one amount, volume, objective, balance, and multiplier scaling
  contract.
- [ ] Reject unresolved exact structural faces before phase optimization.
- [ ] Add reaction-basis, species-order, feed-scale, and reference-gauge
  invariance tests.

Checkpoint: Stage 0 either returns one immutable compiled problem or a typed
failure; it performs no phase search.

## Task 5: implement GREPE Stages 1--4

### Stage 1: physical seeds

- [ ] Enumerate distinct mechanically stable feed-density roots.
- [ ] Run the existing homogeneous reacting owner from each applicable root
  seed.
- [ ] Retain only independently certified, physically distinct local states.

### Stage 2: feasibility master

- [ ] Solve the restricted Phase-I master with exact balance columns.
- [ ] Independently audit its primal/dual evidence.
- [ ] When needed, run bounded feasibility pricing through the shared search
  mechanism.
- [ ] Treat unsuccessful finite feasibility search as unresolved, not model
  infeasibility.

### Stage 3: energy master and reactive pricing

- [ ] Solve and independently audit the restricted energy master.
- [ ] Search jointly over invariant composition coordinates and normalized
  log volume using its conserved multipliers.
- [ ] Count every bounds query and phase evaluation as Provider work.
- [ ] Append only independently reevaluated feasible strict-negative columns
  under immutable insertion IDs.

### Stage 4: topology

- [ ] Keep persistent identity, candidate distinctness, and duplicate merging
  as separate physical decisions.
- [ ] Select at most two liquid candidates.
- [ ] Return scope exceeded if evidence requires a third material phase.

Checkpoint: analytic master values and reduced costs pass; a narrow negative
basin is found within its declared Provider-call budget.

## Task 6: implement GREPE Stages 5--8

### Stage 5: simultaneous final NLP

- [ ] Optimize nonlogarithmic phase amounts, normalized phase species amounts,
  and bounded log volumes simultaneously.
- [ ] Enforce global conserved-component balances and per-phase normalization,
  electroneutrality when applicable, and Provider domains.
- [ ] Start exactly from the feasible master state.
- [ ] Assemble exact objective, constraint Jacobian, and Lagrangian Hessian
  blocks from the shared Provider phase evaluator and Equilibrium chain rules.
- [ ] Verify derivatives independently at interior and trace states.

### Stage 6: support and topology refinement

- [ ] Merge only numerical duplicate phases, then rerun the complete NLP and
  certificates.
- [ ] Retire at most one KKT-inactive phase and solve the full reduced problem.
- [ ] Preserve every small positive species amount.

### Stage 7: final state-matched repricing

- [ ] Freeze the final state revision and multiplier checksum.
- [ ] Rerun reactive pricing against those exact final multipliers.
- [ ] Return to phase generation after any independently certified negative
  state.
- [ ] Reject stale or mismatched pricing evidence.

### Stage 8: final certification

- [ ] Independently check balances, charge, normalization, pressure, reaction
  stationarity, transfer stationarity, bound signs, complementarity, local
  minimum evidence, identity, domains, and resource accounting.
- [ ] Accept only evidence sharing the same state revision.

Checkpoint: all manufactured success and fail-closed branches pass through the
one coupled native controller.

## Task 7: prove limiting reductions and manufactured coupled behavior

- [ ] Constrained homogeneous variant reproduces the existing chemical owner.
- [ ] A reaction-free input qualifying under an installed strong-electrolyte
  capability table reproduces HELD2, including its rejection behavior.
- [ ] A reaction-free input matching the reviewed neutral binary fingerprint
  reproduces the existing neutral HELD route, including its rejection behavior.
- [ ] Other nonreactive inputs remain unsupported through the facade.
- [ ] Analytic reactive two-liquid problem recovers expected phases, amounts,
  volumes, multipliers, balances, and reduced costs.
- [ ] Ascani's hypothetical \(A+B\rightleftharpoons C\) PC-SAFT problem
  reproduces the reported one-phase/two-phase topology changes only after its
  exact source parameters are available as an immutable installed Provider
  artifact; until then it is informative design evidence, not an executable
  gate.
- [ ] Multiple-density-root, narrow-pricing-basin, inactive-phase retirement,
  stale-revision, negative-final-repricing, and search-exhaustion cases pass.
- [ ] Native diagnostics and the private Python seam serialize identical
  scientific payloads.

Checkpoint: private neutral reactive LLE is numerically complete but still has
no physical capability claim.

## Task 8: installed Ascani--Senina neutral reactive LLE

- [ ] Provider supplies one immutable source-complete artifact for Ascani's
  four-species PC-SAFT model hypothesis: water, acetic acid, 1-pentanol, and
  pentyl acetate using the source conventions.
- [ ] Equilibrium consumes that installed artifact only through the public SDK.
- [ ] Encode the reaction, conserved-component matrix, \(K_a=43.99\), source
  standard-state identity, \(318.15\ \mathrm K\), and the distinct Ascani
  \(1\ \mathrm{bar}\) and Senina \(101.3\ \mathrm{kPa}\) pressure records with
  immutable provenance.
- [ ] Classify the Senina homogeneous composition used by Ascani to determine
  \(K_a=43.99\) with PC-SAFT activity coefficients as calibration evidence,
  not a validation case.
- [ ] Record the model discrepancy that Senina used less than
  \(2\ \mathrm{wt}\%\) aqueous HCl while Ascani omitted the catalyst from the
  four-species calculation as an approximation.
- [ ] For each selected Senina tie-line, construct a declared overall feed from
  a fixed convex weight applied to the two tabulated endpoint compositions,
  solve it, and compare calculated endpoints with the experimental
  compositions using source- and model-appropriate tolerances.
- [ ] Do not label the result as reproduction of unpublished Ascani calculated
  endpoints; Ascani reports those endpoints only in a figure.
- [ ] Record phase compositions, phase amounts, volumes, reaction/transfer/
  pressure residuals, final pricing, Provider calls, and wall time.
- [ ] Build an immutable Equilibrium artifact and execute the black-box
  campaign in `ePC-SAFT-validation`.

Checkpoint: neutral reactive LLE has one physical installed-artifact evidence
subject. It still does not imply electrolyte reactive LLE or guaranteed
globality.

## Task 9: consider public `equilibrate`

- [ ] Review the stable private problem/result vocabulary and evidence.
- [ ] Define the smallest typed Python facade without starts, phase counts,
  solver controls, backend selectors, or scaling knobs.
- [ ] Preserve `tp_flash` and `chemical_equilibrium` as frozen public
  operations until
  repeated `equilibrate` parity and neutral reactive-LLE campaigns pass.
- [ ] Decide removal or retention of the older symbols in a separate ADR and
  change; do not pre-commit to aliases.

## Later slice: electrolyte reactive LLE

- [ ] Select a source-complete charged reactive-LLE paper and installed
  Provider parameter bundle.
- [ ] Add phase-specific electroneutrality, electrolyte reference transforms,
  and charged trace evidence without changing the neutral GREPE master.
- [ ] Rerun every standalone, neutral reactive, and Perdomo gate.
- [ ] Execute a separate immutable Validation campaign and capability review.

## Completion criteria

The first plan is complete only when:

- current standalone phase and homogeneous chemistry behavior remains frozen;
- HELD2's two certificate gaps are independently closed;
- the shared search mechanism has no duplicated scientific formulation;
- full GREPE Stages 0--8 pass analytic and source-derived manufactured tests;
- the Ascani--Senina installed-Provider campaign is stable and reproducible;
- every finite result retains non-guaranteed globality; and
- no Provider internals, copied controllers, sequential reaction-then-flash,
  positive phase floor, magnitude pruning, or case-specific production branch
  was introduced.
