# Neutral HELD v1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:executing-plans` and `superpowers:test-driven-development`.
> Stop at every named review checkpoint. Do not push, promote, publish, or
> mutate authority.

**Goal:** Replace the unpromoted fixed-two-phase methane/ethane route with one
bounded `tp_flash` operation implementing Pereira HELD Stages I--III in the
existing equilibrium native module.

**Architecture:** The installed provider remains the only EOS kernel and
supplies `Phi=A/(RT n_ref)`, its gradient, and Hessian in component amounts and
volume. Equilibrium owns the Helmholtz-plus-`PV` transforms, Stage I tunneling,
the binary Stage II outer/lower loop, direct Stage III refinement, Ipopt,
certification, and one public result. Internal stage proofs are reviewed but
never promoted separately.

**Tech stack:** Python 3.13, C++17, pybind11, Ipopt, CMake,
scikit-build-core, pytest, and an exact installed provider wheel containing
the neutral-mixture `epcsaft.native_sdk.v1` tail.

## Global constraints

- Read `AGENTS.md`, `CONTEXT.md`, `ARCHITECTURE.yaml`,
  `docs/designs/2026-07-17-neutral-held-v1.md`, and the retained source
  identity named there before editing.
- Keep one native module, one native target, one provider kernel, one
  `ProviderContext`, and one phase/result representation.
- Do not add provider equations, provider source compilation, density closure,
  third derivatives, generic LP/NLP abstractions, route registries, fallback
  solvers, public stage APIs, public solver settings, or compatibility aliases.
- Keep the accepted pure-saturation route unchanged.
- Keep every checkpoint green. Internal Stage I/II code may coexist with the
  current local fixed-route implementation until the final cutover, but it
  must not be publicly exported.
- Freeze this v1 resource profile in native constants and diagnostics:

```text
stage_i_tunneling_starts:       20
stage_ii_major_iterations_max: 100
stage_ii_lower_starts_per_iter: 20
ipopt_iterations_per_solve:    300
composition_bounds:            [1e-8, 1 - 1e-8]
stage_i_ii_molar_volume_bounds_m3_mol: [1e-5, 1e-1]
stage_iii_phase_volume_bounds_m3:     [1e-5, 1e-1]
tpd_negative_threshold:        -1e-8
epsilon_lambda:                0.5
epsilon_b:                     1e-2
epsilon_x:                     1e-3
relative_density_separation:   1e-3
epsilon_mu:                    1e-6
epsilon_g:                     1e-6
solver_constraint_tolerance:   1e-10
local_physical_tolerance:      1e-8
confirmation_tolerance:        1e-7
stage_iii_candidate_radius_mol: 1e-3
```

- A resource cap produces `search_exhausted`; it never produces a stable or
  accepted result.
- Use exact positive path staging and preserve `.idea/` and `.serena/` state.

---

### Task 1: Freeze the Stage I RED contract

**Files:**

- Create: `tests/test_held.py`

**Produces:** focused failing tests for private native Stage I evaluators; no
production edit or public partial stage.

- [ ] Add private-evaluator tests for `g_bar`, pressure stationarity, reduced
  composition gradient, log-volume gradient/Hessian, TPD sign, reference
  selection, deterministic starts, negative confirmation, no-negative-found,
  and indeterminate failure.
- [ ] Bind every real provider state in these tests to the retained May source
  row identity. Keep synthetic values limited to pure algebraic chain-rule
  fixtures that make no model prediction.
- [ ] Run `python -m pytest tests/test_held.py -q` and confirm RED because the
  private evaluators do not exist. Do not commit the failing state.

### Task 2: Implement the common Helmholtz coordinate transform and Stage I

**Files:**

- Create: `cpp/src/held.cpp`
- Create: `cpp/src/held_stage_i.cpp`
- Create: `cpp/src/held.hpp`
- Modify: `cpp/src/module.cpp`
- Modify: `CMakeLists.txt`
- Modify: `tests/test_held.py`

**Consumes:** existing `ProviderContext::evaluate_mixture` only.

- [ ] Add failing private-evaluator tests for `g_bar`, pressure stationarity,
  reduced composition gradient, log-volume gradient/Hessian, and TPD.
- [ ] Implement the one-mole `(x,log V)` to provider `(n_1,n_2,V)` transform
  and exact chain rules from the approved design.
- [ ] Independently compare objective gradients and Hessian-vector products
  with centered directional differences at one retained real May state.
- [ ] Add the homogeneous feed reference search from liquid-like and
  vapor-like volume starts. Retain every attempt; select the lowest locally
  accepted `g_bar` state.
- [ ] Implement the retained exponential tunneling objective with `p*=1e-3`,
  overflow-safe pole handling, and the versioned deterministic 20-start
  profile. A singular composition or provider failure rejects that attempt.
- [ ] Confirm each negative TPD candidate from a materially perturbed start.
  Distinguish negative, no-negative-found, and indeterminate outcomes.
- [ ] Expose only private underscored test bindings for scalar evaluation and
  Stage I diagnostics. Do not export a Python stability operation.
- [ ] Run focused source tests and the unchanged saturation suite.

**Review checkpoint A:** permanent lab review must accept the TPD sign,
normalization, derivative transform, deterministic search semantics, and
absence of a false stability/globality claim before Task 3.

### Task 3: Implement the binary Stage II outer envelope

**Files:**

- Create: `cpp/src/held_stage_ii.cpp`
- Modify: `cpp/src/held.hpp`
- Modify: `cpp/src/module.cpp`
- Modify: `tests/test_held.py`

- [ ] Write failing tests using explicit synthetic cut coefficients for a
  unique vertex, tied vertices, redundant cuts, nearly parallel cuts,
  infeasibility, and an unbounded/missing-endpoint initialization.
- [ ] Implement one binary-only `(v,lambda)` finite-envelope solver by
  enumerating cut intersections, checking every constraint, applying declared
  tie tolerances, and returning active cut identities.
- [ ] Keep the solver private and specific to HELD; add no matrix, simplex, LP,
  or optimization framework.
- [ ] Construct the initial cut set by locally polishing liquid-like and
  vapor-like volume starts at both composition bounds. Require at least one
  valid endpoint cut on each side of the admitted feed.
- [ ] Compare the production envelope result against an independent exhaustive
  Python vertex oracle for every test fixture.
- [ ] Run `python -m pytest tests/test_held.py -k 'outer or cut' -q`.

### Task 4: Implement Stage II lower solves and candidate lifecycle

**Files:**

- Modify: `cpp/src/held_stage_ii.cpp`
- Modify: `cpp/src/held.hpp`
- Modify: `cpp/src/module.cpp`
- Modify: `tests/test_held.py`

- [ ] Write failing tests for one lower solve below the current upper bound,
  one valid but non-tight cut, deterministic start ordering, bound
  monotonicity, candidate deduplication, two distinct candidates, a third
  candidate, and budget exhaustion.
- [ ] Implement the exact Stage II lower objective in `(x,log V)` with
  provider gradient/Hessian transforms and Ipopt settings from the frozen
  profile.
- [ ] Alternate shifted-previous, component-near, and deterministic stratified
  starts. Add every locally accepted point as a valid cut; stop the inner
  search only when `L <= UBD` or the 20-start cap is exhausted.
- [ ] Apply `epsilon_b`, zero-safe `epsilon_lambda`, composition separation,
  and relative-density separation exactly as designed. Merge duplicates by a
  deterministic lowest-objective rule.
- [ ] Return two candidates to Stage III. Return `scope_exceeded` for more than
  two and `search_exhausted` at either declared cap.
- [ ] Retain a compact trace of bounds, active cuts, accepted lower points,
  candidates, and rejection reasons without exporting the entire native
  solver object graph.
- [ ] Run focused Stage II tests plus Stage I and saturation regression tests.

**Review checkpoint B:** permanent lab review must accept the outer/lower bound
semantics, every-cut validity, candidate lifecycle, resource-cap behavior, and
minimality before Task 5.

### Task 5: Refactor direct Stage III and close the feedback loop

**Files:**

- Create: `cpp/src/held_stage_iii.cpp`
- Modify: `cpp/src/held.cpp`
- Modify: `cpp/src/held.hpp`
- Modify: `cpp/src/two_phase_flash.cpp`
- Modify: `cpp/src/module.cpp`
- Modify: `tests/test_held.py`

- [ ] Write failing tests for lever-rule initialization, `+/-1e-3 mol`
  candidate-neighborhood bounds, infeasible balances, exact Stage III
  derivatives, phase collapse, KKT failure, HELD gap failure,
  chemical-potential failure, successful confirmation, and return to Stage II.
- [ ] Move the useful provider phase evaluation, direct total-`g_bar`
  objective, linear material balances, multiplier retention, and exact
  log-volume Hessian logic into internal Stage III ownership.
- [ ] Compute binary phase fractions by the lever rule and center amount bounds
  on Stage II candidates. Never seed a phase amount or composition not derived
  from the candidate pair.
- [ ] Require strict local material, pressure, KKT, inactive-bound,
  distinct-phase, and confirmation checks plus zero-safe `epsilon_mu` and
  `epsilon_g` HELD criteria.
- [ ] On infeasibility or failed HELD criteria, continue Stage II while budget
  remains. Do not return the old local fixed-phase solution as fallback.
- [ ] Reject trace-bound or third-phase candidates as `scope_exceeded`.
- [ ] Run focused Stage III/controller tests and the full source suite.

**Review checkpoint C:** permanent lab review must accept the complete
Stage-I/II/III controller and feedback semantics before public cutover.

### Task 6: Cut over to `tp_flash` and delete the displaced route

**Files:**

- Modify: `cpp/src/module.cpp`
- Modify: `CMakeLists.txt`
- Modify: `src/epcsaft_equilibrium/_api.py`
- Modify: `src/epcsaft_equilibrium/__init__.py`
- Modify: `src/epcsaft_equilibrium/_equilibrium.pyi`
- Modify: `tests/test_held.py`
- Modify: `tests/test_saturation.py`
- Delete: `cpp/src/two_phase_flash.cpp`
- Delete: `cpp/src/two_phase_flash.hpp`
- Delete: `tests/test_two_phase_flash.py`

- [ ] First add the RED public contract tests requiring
  `tp_flash(model,T,P,z) -> TpFlashResult`, existing `PhaseState`, one compact
  `HeldDiagnostics`, existing `FlashError`, the frozen outcome vocabulary, and
  `globality_certificate="not_guaranteed"`. Add input rejection tests for
  units, finiteness, feed, fingerprint, ABI, domain, and forbidden public
  initial guesses/solver options.
- [ ] Bind the complete controller once as private native `_solve_tp_flash` and
  expose only public `tp_flash`.
- [ ] Preserve high-value provider transport, derivative, KKT, confirmation,
  and collapsed-state assertions by moving their behavior checks into
  `tests/test_held.py`; do not bulk-copy the old test file.
- [ ] Delete `two_phase_flash`, `TwoPhaseFlashResult`, duplicate diagnostics,
  native bindings, CMake sources, stubs, exports, and docs in the same commit.
- [ ] Add negative-space checks for retired symbols and prove no compatibility
  alias, second native target, provider source, private provider import,
  density-root path, or public stage API exists.
- [ ] Run the complete source suite and inspect exported Python/native symbols.

### Task 7: Build the isolated artifact and reconcile package documentation

**Files:**

- Modify: `README.md`
- Modify: `CONTEXT.md`
- Modify: `ARCHITECTURE.yaml`
- Modify: `AGENTS.md` only if the active scope sentence requires narrowing
- Create: one candidate receipt under the existing receipt category
- Modify: `docs/designs/2026-07-17-neutral-held-v1.md` status only

- [ ] Update current capability prose to distinguish accepted pure saturation
  from local `neutral-held-v1`; retain `globality=not_guaranteed`, exact
  provider fingerprint/domain, source blob, exclusions, and authority effect
  `none`.
- [ ] Build against the exact retained provider-tail wheel in an isolated
  Python 3.13 environment:

```bash
uv run --isolated --no-project --python 3.13 \
  --with "$PINNED_PROVIDER_WHEEL" --with scikit-build-core --with pybind11 \
  -- uv build --no-build-isolation --wheel .
```

- [ ] Install only the provider and equilibrium wheels in a second isolated
  environment and run the exact source tests from outside the checkout.
- [ ] Inspect wheel RECORD, import origins, dynamic dependencies, native
  symbols, public exports, CMake target count, tracked files, and architecture
  ratchet. Confirm no lab, migration, sibling source, or provider
  implementation dependency.
- [ ] Record exact commits, trees, wheel/header hashes, commands, compiler,
  Ipopt/Python environments, test counts, and known failures.
- [ ] Run cleanup, `git diff --check`, and final status review. Commit only a
  green local candidate and stop for independent review and validation.
- [ ] After the complete local cutover commit, obtain permanent-lab review of
  the whole candidate, including proof that `two_phase_flash` and its displaced
  surface were removed without an alias. Internal checkpoint reviews do not
  substitute for this final post-cutover candidate review.

## Completion boundary

Completion means a strongest truthful local HELD candidate. It does not mean
provider-tail promotion, equilibrium promotion, push, release, publication,
or a guaranteed global phase-equilibrium capability. Those remain migration
and user approval gates.
