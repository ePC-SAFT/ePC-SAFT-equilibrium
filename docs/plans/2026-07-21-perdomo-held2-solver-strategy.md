# Perdomo HELD2 Solver-Strategy Implementation Plan

Status: Tasks 0--6 implemented on `main`; Tasks 7--8 remain bounded future work

Authority effect: none

> **For agentic workers:** REQUIRED SUB-SKILL: use
> `superpowers:executing-plans`, `superpowers:test-driven-development`, and
> `superpowers:verification-before-completion`. Stop at every named review
> checkpoint. Do not publish a route, build a candidate artifact, run a sibling
> Validation campaign, or request promotion unless a later bounded assignment
> explicitly authorizes it.

**Goal:** Implement the Perdomo HELD2 numerical strategy with the numerical
method best suited to each stage: deterministic pressure-root enumeration for
density topology, DIRECT-L for low-dimensional bounded exploration, HiGHS for
the Stage-II linear program, and exact-Hessian Ipopt for smooth constrained
local refinement.

**Architecture:** The installed Provider remains the sole EOS and nonlinear
derivative owner. Equilibrium owns the modified-coordinate maps, fixed-pressure
root and lower-envelope construction, finite-search controllers, LP and NLP
formulations, independent physical certificates, diagnostics, and fail-closed
outcomes. Global exploration and local convergence are separate jobs; neither
is a thermodynamic or mathematical globality certificate.

**Target tech stack:** Python 3.13, C++17, pybind11, CMake, Ipopt, pinned NLopt
2.11.0, pinned HiGHS 1.15.1, scikit-build-core, pytest, and one exact installed
Provider artifact. Current `main` fetches the two pinned source archives by
SHA-256 and links them into the private extension.

## Outcome Proof

**Intent:** Implement one source-faithful HELD2 controller whose numerical owner matches each stage and whose certificates remain independent of solver status.
**Current Behavior:** The private integrated controller uses deterministic pressure roots, DIRECT-L Stage I, HiGHS Problem (64), deterministic basin discovery plus exact-Hessian Ipopt Problem (65), and hardened exact-Hessian Ipopt Problem (67). The immutable Khudaida candidate fails closed in Stage II, public electrolyte dispatch is disabled, and electrolyte LLE is unadmitted.
**Expected Outcome:** A separate bounded assignment supplies installed two-liquid Stage-II/III evidence without changing the retained controller contracts.
**Target Output:** One existing `TpFlashResult` or `FlashError` carrying solver, numerical, physical, search, root-completeness, predictive, and globality evidence.
**Owner:** Equilibrium owns HELD2 coordinates, stage solvers, controllers, and certificates; the installed Provider owns EOS values, domains, pressure, packing, and exact tensors.
**Interface:** Existing public `tp_flash(model, T, P, z)` and private solver-neutral HELD2 evaluators; no public stage or optimizer selector.
**Cutover:** Stage-I/II private routing has completed its guarded cutover. Public electrolyte dispatch still requires every remaining checkpoint and exact artifact gate.
**Replaced Path:** The fixed-start Stage-I Ipopt route, analytic Stage-II upper envelope, pre-exploration Stage-II controller, and their legacy baseline fixture are removed. Focused manufactured formulation oracles remain.
**Evidence:** Manufactured root and basin topologies, derivative parity, LP oracle parity, complete search traces, Stage-III certificates, and one immutable installed Khudaida tracer.
**Acceptance Proof:** The target owner map is visible in diagnostics, all hard gates pass without tolerance changes, displaced production paths are unreachable, and finite-search globality remains `not_guaranteed`.
**Stop Criteria:** Stop on missing Provider correction, contradictory pre-change evidence, failed root or certificate accounting, unauthorized runtime scope, or inability to reproduce the frozen regression.
**Avoid:** Do not use SLSQP as a policy substitute, copy EOS equations, return fake penalties, hard-code a chemistry, relax tolerances, or publish finite-search globality.
**Risk:** Solver migration can hide certificate regressions or leave dependency status stale unless each cutover preserves the frozen evaluator, diagnostics, and documentation.

## Implementation Boundaries

**Files To Create:** Add focused internal HELD2 source files only when cohesion requires them; no new public module, route, result family, or local tracker.
**Files To Modify:** `cpp/src/held2.cpp`, `cpp/src/held2.hpp`, `cpp/src/held2_stage_iii.cpp`, `CMakeLists.txt`, `tests/test_held2.py`, and governing documentation after behavior lands.
**Files To Avoid:** Provider implementation, sibling Validation source, archived artifacts, accepted receipts, neutral HELD owners, unrelated public APIs, and user-owned IDE state.
**Source Of Truth:** The Perdomo HELD2 design, Perdomo Problems (63)--(67), the exact Provider contract, frozen regression artifacts, and this numerical strategy.
**Read Path:** Canonical design and authority index, current HELD2 sources and tests, Provider public SDK metadata, and immutable assigned evidence.
**Write Path:** The named Equilibrium HELD2 sources, tests, build metadata, canonical docs, and task-authorized local evidence only.
**Integration Points:** Provider SDK, deterministic pressure-root service, NLopt DIRECT-L, HiGHS LP, Ipopt Stage-II/III TNLPs, Python diagnostics, and later installed-artifact validation.
**Migration Or Cutover:** Land evaluator extraction, root envelope, HiGHS, Stage-I strategy, Stage-II discovery, and Stage-III lifecycle through the ordered checkpoints before any public route changes.
**Replaced Path Handling:** Remove replaced Stage-I Ipopt and analytic upper-envelope production routes after parity; keep only focused manufactured test oracles, without aliases or fallback routing.
**Acceptance Proof Gate:** All Project Truss task use cases, manufactured matrix, frozen pass/fail partitions, dependency diagnostics, and assigned installed tracer must pass under unchanged scientific gates.

## 1. Governing contracts

Read before implementation:

- `AGENTS.md`, `CONTEXT.md`, and `ARCHITECTURE.yaml`;
- `docs/phase-equilibrium.md`;
- `docs/designs/2026-07-21-perdomo-held2.md`;
- Perdomo et al. (2025), Algorithm 1 and Problems (63)--(67); and
- the exact Provider capability and artifact contract assigned to the work.

The design document owns the thermodynamic formulation. This plan owns only
the approved numerical decomposition and implementation order.

The implementation must preserve:

- the Perdomo eliminated-ion and modified-mole coordinates;
- Provider ownership of Helmholtz values, pressure, packing relations,
  admissible domains, exact first derivatives, and exact second derivatives;
- the existing pressure, KKT, complementarity, Step-6, and TPD acceptance
  tolerances unless a separate scientific change is admitted;
- fixed-physical-volume derivatives in Perdomo Eq. (66);
- same-major `UBD` and multiplier use in Stage II;
- composition-**or**-packing candidate distinctness;
- the original-coordinate Stage-II dual-pullback certificate;
- direct total-free-energy Stage III; and
- `globality_certificate="not_guaranteed"` for every finite search.

Do not add:

- a second EOS, packing relation, or density solver;
- an Equilibrium-owned CppAD recording of Provider equations;
- a chemistry-, component-, or Khudaida-specific code path;
- a caller-selectable optimizer or silent solver fallback;
- artificial finite penalties for Provider or root-accounting failures;
- tolerance relaxation to accept a tracer; or
- a claim that DIRECT-L, Sobol sampling, multistart Ipopt, or a completed HELD2
  controller proves global equilibrium.

## 2. Current implementation versus approved target

Current `main` retains a private HELD2 development implementation. The
homogeneous and trial-composition services enumerate pressure roots;
`held2_stage_i_direct.cpp` owns DIRECT-L Stage I;
`held2_stage_ii_upper.cpp` owns the HiGHS Problem-(64) LP; and
`held2_stage_ii_basin.cpp` proposes deterministic physical basin
representatives before `Held2SearchTnlp` refines each representative with
exact-Hessian Ipopt. The pre-exploration Stage-II controller is reachable only
through the manufactured Stage-II test oracle. Stage III retains its separate
`Held2StageIIITnlp`. The installed public-dispatch WIP remains archived at
`archive/held2-pre-strategy-2026-07-21` and has no production authority.

The current owner map is:

| Operation | Target numerical owner | Required interpretation |
| --- | --- | --- |
| Homogeneous reference | deterministic pressure-root enumeration | topology and branch-selection problem |
| Stage I, Problem (63) | NLopt `NLOPT_GN_DIRECT_L` over a pressure-certified reduced TPD envelope | finite search for a certified negative witness |
| Stage II, Problem (64) | HiGHS LP | exact finite linear program with an independent residual audit |
| Stage II, Problem (65), discovery | continuation, physical seeds, Sobol sampling, and stall-triggered DIRECT-L | discover distinct physical basins |
| Stage II, Problem (65), refinement | exact-Hessian Ipopt | converge each retained smooth basin to a local KKT candidate |
| Stage III, Problem (67) | exact-Hessian Ipopt | direct total-free-energy refinement and active-phase lifecycle |
| Scientific acceptance | Equilibrium certificates | independent of optimizer return status |

SLSQP is not the default replacement for Ipopt. It is another local optimizer,
does not improve basin coverage by itself, and would not consume the exact
Provider-derived Lagrangian Hessian used by the current Ipopt path. It may be
used only as a bounded comparator if a later assignment requests an optimizer
A/B study.

## 3. Shared pressure-certified lower-envelope service

Create one internal service used by the homogeneous reference, Stage I, and
Stage-II exploration. For a fixed feasible modified composition `u`, it must:

1. reconstruct the complete explicit-species composition and verify
   normalization and phase electroneutrality;
2. acquire the exact composition-specific Provider molar-volume domain;
3. cover the declared log-volume domain with the canonical deterministic root
   scan, including pressure-residual stationary-point searches;
4. refine every detected ordinary or tangential pressure root;
5. retain the unchanged pressure certificate;
6. classify roots with exact fixed-composition `dP/dV`;
7. evaluate reduced `A + P0 V` at every strict-stable interior root; and
8. return every root, interval-accounting record, branch identity, failure,
   and the selected lower-envelope value.

Strict-stable roots have `dP/dV < 0`. Unstable roots remain topology evidence
but are not lower-envelope candidates. Marginal, boundary, tied, or unresolved
roots must not be silently discarded. A tie makes the homogeneous reference
nonunique. For a derivative-free Stage-I or Stage-II exploration value, retain
all tied strict-stable branches and accept the common lower value only when its
audited objective enclosure resolves the required negative or ordering test.
Incomplete interval accounting remains indeterminate.

The service must report finite numerical completion separately from root
completeness. Pointwise exact derivatives and deterministic subdivision do not
constitute interval proof; without validated interval extensions the retained
status remains `root_completeness="not_proven"`.

Recommended internal responsibilities, subject to a code-surface review:

```text
PressureRootEnumerator
Held2FixedPressureEnvelope
StageITpdEnvelope
StageIILowerEnvelope
```

These are internal responsibilities, not new public APIs or a generic solver
framework. Keep them in the existing HELD2 module unless extraction materially
improves cohesion or testability.

## 4. Stage I: DIRECT-L over the reduced TPD envelope

For each independent modified composition `u`, define

```text
d_hat(u) = min over eligible strict-stable pressure roots V of TPD(V, u).
```

At fixed composition the tangent term is independent of volume, so pressure
stationarity supplies the eligible volume branches. The global explorer sees
only the `C - 2` bounded composition coordinates. For the five-species
Khudaida topology this is a three-dimensional search rather than a joint
composition-and-volume local NLP.

Use deterministic NLopt `NLOPT_GN_DIRECT_L` because it systematically
subdivides a finite bound-constrained search domain and does not require a
derivative of the potentially nonsmooth lower envelope. The existing closed
feasible-simplex forward map supplies the finite box. Exact chart boundaries
may be noninjective; deduplicate evaluations and witnesses in physical
composition and volume coordinates.

The Stage-I terminal logic is asymmetric:

- one finite, feasible, materially distinct, pressure-certified point below
  the existing negative TPD threshold is a valid instability witness;
- no local-optimum or TPD-gradient claim is required for that witness;
- Provider or root-accounting failures encountered elsewhere remain visible
  but do not erase an already certified negative existence witness;
- `no_negative_witness_detected` requires complete accounting of every
  evaluation required by the declared finite search; and
- budget exhaustion or incomplete required evaluations produce
  `indeterminate`, never `stable`.

Do not return a large objective penalty for an invalid Provider evaluation or
an incomplete root scan. Stop the declared search, retain the exact failure,
and apply the asymmetric terminal logic above.

Keep the current deterministic Ipopt multistart strategy as a named regression
oracle during migration. Do not silently substitute the new strategy under an
old resource-profile label. The bounded implementation assignment must freeze
the DIRECT-L evaluation budget and completion contract from evidence; this
plan does not invent new numerical tolerances or resource values.

## 5. Stage II Step 4: HiGHS upper LP

Translate Perdomo Problem (64) and the current cut convention directly into a
sparse HiGHS linear program. Do not rederive or flip cut signs during the
adapter change.

The internal result must retain at least:

```text
solver return and model status
primal and dual feasibility status
UBD and modified multipliers
cut slacks and cut duals
active cut identities
independently recomputed primal and dual residuals
HiGHS version
```

A successful HiGHS call is not the certificate. Recompute the LP objective,
all row activities, bounds, and the required residual norms in Equilibrium.

Retain the present analytic finite-envelope implementation only as a
manufactured test oracle. Production Problem (64) must have one HiGHS owner.

## 6. Stage II Step 5: global basin discovery plus Ipopt

Problem (65) needs accurate local minima, useful cuts, and multiple distinct
phase candidates. A global explorer returning only one best point is
insufficient, while one Ipopt start cannot enumerate basins. Use two layers.

### 6.1 Discovery layer

Build the deterministic start pool in this order:

1. certified minima from the preceding major iteration;
2. active and recently added cut states;
3. retained Stage-I negative witnesses;
4. the homogeneous feed/reference state;
5. source-independent aqueous-rich, organic-rich, and boundary-aware physical
   seed families derived from the feasible coordinates, not component names;
6. dense- and lower-density branches retained by the pressure envelope;
7. deterministic Sobol compositions; and
8. one stall-triggered DIRECT-L exploration of the reduced Problem-(65)
   envelope after the declared no-progress condition.

At a fixed composition the multiplier term is independent of volume, so the
shared pressure-root service may evaluate

```text
L_hat(u; lambda_k) = min over eligible strict-stable pressure roots V
                     of L(V, u; lambda_k).
```

The discovery layer records all valid branch states, clusters them in physical
modified composition and volume or packing fraction, and forwards only
distinct representatives. It does not certify final KKT conditions or claim a
global lower bound.

### 6.2 Local refinement layer

Run the existing smooth Problem-(65) formulation through Ipopt from every
retained basin representative using:

- the current feasible-simplex and packing chart;
- the exact transformed objective gradient;
- the complete exact transformed Lagrangian Hessian, including chart-curvature
  terms;
- unchanged variable bounds and physical tolerances; and
- retained terminal bound and constraint multipliers.

Every terminal then receives the independent original-coordinate certificate,
including Provider and domain validity, pressure closure, primal feasibility,
physical KKT stationarity, multiplier signs, complementarity, exact dual
pullback and reconstruction, and fixed-volume Step-6 derivatives. Ipopt status
is retained separately and does not accept or reject a physical state by
itself.

### 6.3 Cuts and candidates

Maintain separate predicates:

```text
cut_eligible
candidate_phase_eligible
```

A Provider-valid finite state can define a valid affine cut without satisfying
all Perdomo Eq. (66) candidate conditions. Candidate eligibility additionally
requires the same-major dual-value test, fixed-volume component-gradient test,
trace-component rule, and composition-or-packing distinctness.

Repeated starts reaching the same physical basin are robustness evidence, not
additional phases. A completed major with no new cut, basin, or candidate may
trigger the one declared exploration escalation. If escalation also makes no
progress, return a bounded indeterminate/no-progress outcome.

## 7. Stage III: direct total-free-energy Ipopt NLP

Keep the existing general `Held2StageIIITnlp` as the Stage-III numerical owner.
It must continue to solve Perdomo Problem (67), not a residual-only surrogate.

Stage III may call Ipopt more than once to retire a KKT-inactive phase, merge a
certified duplicate, re-solve the active set, or re-enter after Stage-II
feedback. “Ipopt owns Stage III” does not mean exactly one invocation.

A phase may be retired only from phase-amount bound, multiplier sign,
complementarity, reduced-derivative, and remaining-balance evidence. A small
phase amount alone is not sufficient. Retire at most one phase per lifecycle
step and re-solve the reduced active problem.

Final acceptance requires the design-owned material, charge, pressure,
modified-potential, domain, amount, distinctness, KKT, and total-free-energy
checks. A small transformed KKT residual cannot override a failed physical
potential equality.

## 8. Solver-neutral evaluation and dependency boundary

Current `main` provides the shared HELD2 value, gradient, Hessian,
physical-state, and Provider-status evaluation used by Ipopt and the envelope
explorers. New solver work must consume that evaluator rather than create a
parallel objective formula.

The Provider remains the only owner of EOS values and nonlinear tensors.
Equilibrium applies only the exact coordinate transformations

```text
gradient_y = J^T gradient_w

Hessian_y = J^T Hessian_w J
            + sum_k gradient_w[k] Hessian_y(w_k).
```

Current CMake fetches NLopt 2.11.0 and HiGHS 1.15.1 from fixed upstream source
archives with SHA-256 verification. `THIRD_PARTY_NOTICES.md` records their
licenses, and the wheel audit must reject leaked headers, libraries, or extra
extension modules.

## 9. Required status and diagnostic fields

The completed implementation must identify numerical ownership explicitly:

```text
homogeneous_reference_method
stage_i_search_strategy
stage_i_volume_elimination_method
stage_ii_upper_solver
stage_ii_global_explorer
stage_ii_local_solver
stage_iii_solver
derivative_source
provider_fingerprint
equilibrium_commit
globality_certificate
```

Target strategy identifiers are:

```text
homogeneous_reference_method = deterministic_pressure_root_enumeration
stage_i_search_strategy = direct_l_pressure_certified_tpd_envelope
stage_i_volume_elimination_method = strict_stable_pressure_root_lower_envelope
stage_ii_upper_solver = highs_lp
stage_ii_global_explorer = continuation_sobol_direct_l
stage_ii_local_solver = ipopt_exact_hessian
stage_iii_solver = ipopt_exact_hessian
derivative_source = provider_exact
globality_certificate = not_guaranteed
```

Per Stage-II attempt retain the major and start identity, start source family,
internal and physical start and terminal, Ipopt status, Provider status,
objective, pressure residual, chart and physical KKT residuals, terminal duals,
chart conditioning, pulled-back physical multipliers, complementarity, cut and
Step-6 eligibility, and physical basin/candidate identity.

## 10. Implementation tasks

### Task 0: Freeze the exact subject and regression evidence

**Files:** no production edits

**Use Cases:**

- Create visible pre-change evidence before any solver cutover so later changes cannot redefine success.

- [x] Record the Equilibrium commit/tree, Provider artifact hash and source
  identity, model-bundle fingerprint, compiler, Ipopt version, and linear
  solver.
- [x] Replay the current manufactured HELD2 suite and the assigned archived or
  reconstructed Stage-II dual-pullback fixture.
- [x] Freeze every declared start and terminal, including failed attempts.
- [x] Confirm the known pass/fail classification before any solver change. If
  the historical 25-recovered/11-narrow/12-broad partition is not reproducible
  from the assigned artifact, record the evidence gap rather than fabricating
  parity.

**Review checkpoint A:** approve the exact pre-change evidence and immutable regression
artifact before refactoring.

### Task 1: Extract one solver-neutral HELD2 evaluator

**Files:** `cpp/src/held2.cpp`, `cpp/src/held2.hpp`, `tests/test_held2.py`

**Use Cases:**

- Prove evaluator parity and prevent duplicate scientific formulas across solver adapters.

- [x] Add failing parity tests for value, gradient, Hessian, physical lift,
  pressure, and Provider failure propagation.
- [x] Extract the scientific evaluation from `Held2SearchTnlp` without changing
  coordinates, bounds, formulas, or tolerances.
- [x] Make the current Ipopt adapter consume the extracted evaluator.
- [x] Prove start-by-start and derivative parity with the frozen pre-change evidence.

### Task 2: Generalize the deterministic pressure-root envelope

**Files:** `cpp/src/held2.cpp`, `cpp/src/held2.hpp`, `tests/test_held2.py`

**Use Cases:**

- Make density-branch topology and root-completeness evidence visible before Stage-I or Stage-II exploration consumes the envelope.

- [x] Add manufactured tests for one root, stable--unstable--stable roots,
  close roots, a tangential root, a branch switch, a boundary root, an invalid
  interval, a tied stable objective, and deduplication.
- [x] Reuse the canonical fixed-feed root implementation at arbitrary feasible
  modified compositions.
- [x] Return complete interval, root, branch, mechanical, and failure records.
- [x] Keep `root_completeness="not_proven"` unless validated interval evidence
  is later added.

**Review checkpoint B:** approve root topology, failure accounting, and the
absence of a duplicate EOS/density implementation.

### Task 3: Add HiGHS for Stage-II Problem (64)

**Files:** `CMakeLists.txt`, `cpp/src/held2.cpp`, `cpp/src/held2.hpp`,
`tests/test_held2.py`

**Use Cases:**

- Replace the production upper-envelope path only after HiGHS matches the analytic oracle and passes independent LP acceptance checks.

- [x] Verify the installed HiGHS CMake target and record its version.
- [x] Add failing analytic-envelope parity tests covering unique, tied,
  redundant, nearly parallel, infeasible, and unbounded cut systems.
- [x] Implement the sparse LP adapter with primal/dual residual recomputation.
- [x] Retain the analytic implementation only as a test oracle.

**Review checkpoint C:** approve row signs, active cuts, objective parity, and
independent LP certification before production routing changes.

### Task 4: Add Stage-I DIRECT-L as an alternative strategy

**Files:** `CMakeLists.txt`, `cpp/src/held2.cpp`, `cpp/src/held2.hpp`,
`tests/test_held2.py`

**Use Cases:**

- Find a certified negative witness under a declared budget and remove the replaced Ipopt Stage-I route after parity evidence is captured.

- [x] Verify the installed NLopt CMake target, C++ callback/forced-stop API,
  and version.
- [x] Add failing tests for a narrow negative pocket, a nonnegative envelope,
  a pressure-branch switch, physical boundary points, Provider failure, root
  failure, budget exhaustion, and negative-witness asymmetry.
- [x] Implement deterministic `NLOPT_GN_DIRECT_L` over the closed feasible
  composition chart and the shared TPD envelope.
- [x] Remove the replaced fixed-start Ipopt Stage-I strategy after parity evidence was captured.
- [x] Freeze a versioned Stage-I evaluation budget from measured evidence.

**Review checkpoint D:** compare strategies under declared budgets using first
certified negative witness, best certified TPD, Provider calls, failure
accounting, reproducibility, and honest globality wording.

### Task 5: Rerun the corrected Stage-II controller before redesign

**Files:** no strategy edits until evidence is captured

**Use Cases:**

- Distinguish certificate recovery from missing-basin evidence before authorizing any Stage-II discovery change.

- [x] Replay the current Stage-II controller with the accepted dual pullback,
  fixed-volume Step 6, and same-major ordering.
- [x] Determine whether candidate starvation remains after certificate repair.
- [x] Identify missing basins, repeated basins, cut deficiencies, and genuine
  KKT failures separately.

**Review checkpoint E:** authorize only the smallest discovery-layer change
supported by the replay.

### Task 6: Add Stage-II physical basin exploration around Ipopt

**Files:** `cpp/src/held2.cpp`, `cpp/src/held2.hpp`, `tests/test_held2.py`

**Use Cases:**

- Discover and deduplicate physical basins, then require Ipopt and independent certificates before candidate acceptance.

- [x] Add manufactured multiple-basin, same-composition/different-density,
  different-composition/same-density, duplicate-start, stalled-search, and
  Provider-failure tests.
- [x] Add continuation, cut-state, Stage-I witness, physical seed, and Sobol
  start families without chemistry-specific constants.
- [x] Cluster in physical composition and volume/packing coordinates.
- [x] Add one declared stall-triggered DIRECT-L envelope exploration.
- [x] Refine every retained representative through the existing exact-Hessian
  Ipopt Stage-II solve.
- [x] Apply independent physical KKT, cut, and Step-6 candidate predicates to
  every terminal.

**Review checkpoint F:** approve basin coverage, unchanged certificates, and
bounded no-progress behavior before Stage III is exercised on a new pool.

### Task 7: Harden the Stage-III lifecycle

**Files:** `cpp/src/held2_stage_iii.cpp`, `cpp/src/held2.hpp`,
`tests/test_held2.py`

**Use Cases:**

- Retire only KKT-inactive phases, re-solve the active set, and keep failed physical-potential evidence visible despite Ipopt success.

- [x] Add failing tests for KKT-active trace phases, KKT-inactive phase
  retirement, one-at-a-time retirement, active-set re-solve, duplicate merge,
  Stage-II feedback, and physical-potential failure despite Ipopt success.
- [x] Keep the direct total-free-energy objective and exact Hessian.
- [x] Implement only evidence-backed retirement and re-solve semantics.
- [x] Preserve the fail-closed trace-refinement gate until logarithmic trace
  refinement is separately implemented and verified.

### Task 8: Integrate, document, and build an isolated candidate only when authorized

**Files:** `README.md`, `CONTEXT.md`, `ARCHITECTURE.yaml`,
`docs/phase-equilibrium.md`, `docs/designs/2026-07-21-perdomo-held2.md`, and
the assigned receipt/evidence locations

**Use Cases:**

- Complete the governed cutover proof, remove displaced production routing, and stop before Validation, publication, or promotion authority.

- [x] Update current-versus-target documentation after each landed strategy
  component and record the dependency state supported by current CMake.
- [x] Run the compact test suite, derivative checks, manufactured matrix, and
  exact assigned tracer under one immutable installed Provider artifact.
- [x] Retain complete start, branch, cut, candidate, Stage-III lifecycle, and
  failure diagnostics.
- [x] Stop after local candidate evidence. Publication, Validation mutation,
  promotion, and authority transfer require separate authorization.

## 11. Minimum verification matrix

| Test | Expected result | Defect caught |
| --- | --- | --- |
| One pressure root | one unique strict-stable envelope branch | root-service regression |
| Stable--unstable--stable roots | both stable roots retained; middle root unstable | missing or misclassified middle root |
| Easier gas basin, lower dense objective | dense branch selected | selection by numerical ease |
| Tangential pressure root | marginal root detected without sign change | sign-change-only scanning |
| Pressure-branch switch | selected identity changes and kink is recorded | false global smoothness assumption |
| Narrow negative TPD pocket | certified negative witness under declared budget | local-basin dependence |
| Nonnegative TPD envelope | finite no-negative result; no stability proof | globality overclaim |
| Required Provider/root failure | indeterminate no-negative search; no penalty bypass | hidden invalid domain |
| Manufactured Stage-II LP | HiGHS matches analytic optimum and active cuts | LP sign or bound error |
| Two Step-II lower basins | both refined and independently certified | one-best-point discovery |
| Duplicate starts | one physical basin and candidate | start count mistaken for phase count |
| Same composition, different density | distinct basins | composition-only deduplication |
| Different composition, same density | distinct basins | density-only deduplication |
| Frozen dual-pullback terminals | pre-change pass/fail partition unchanged | certificate drift |
| Perturbed physical gradient | KKT certificate fails | flexible multipliers hiding error |
| Fixed-volume Step 6 | only fixed-volume derivative satisfies source test | coordinate-basis recurrence |
| Mixed-major Step 6 | deliberate old/new multiplier mix fails | controller state mixing |
| Tiny phase amount with descent direction | phase retained | numerical cutoff retirement |
| Tiny phase amount with valid bound KKT | phase retired and active set re-solved | failure to retire inactive phase |
| Stage-III modified-potential mismatch | physical acceptance fails despite Ipopt success | transformed-KKT overacceptance |
| Completed finite searches | `globality_certificate="not_guaranteed"` | false global proof |

## 12. Completion boundary

This plan is implemented only when:

- each stage uses the target numerical owner listed in Section 2;
- current and archived evidence remain distinguishable from the new subject;
- the pressure envelope and all solver adapters share one scientific evaluator;
- every failure is retained without artificial penalty substitution;
- Stage-II discovery produces distinct physical basins and Ipopt refines each
  retained representative;
- independent numerical and physical certificates, rather than optimizer
  success, own acceptance;
- all hard scientific tolerances and source-faithful predicates remain
  unchanged;
- no chemistry-specific solver path or duplicate Provider equation exists;
- finite-search, root-completeness, and globality claims remain honest; and
- a separate authority process, not this plan or its implementation, decides
  whether any electrolyte capability is admitted.

## 13. Primary numerical-library references

- [NLopt algorithms: DIRECT and DIRECT-L](https://nlopt.readthedocs.io/en/stable/NLopt_Algorithms/)
- [HiGHS C++ library interface](https://ergo-code.github.io/HiGHS/dev/interfaces/cpp/library/)
- [Ipopt options: exact and limited-memory Hessians](https://coin-or.github.io/Ipopt/OPTIONS.html)
