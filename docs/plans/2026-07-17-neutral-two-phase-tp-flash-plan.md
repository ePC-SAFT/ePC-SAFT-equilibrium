# Neutral Two-Phase TP Flash Implementation Plan

Status: executed local candidate; validation NON_ADMISSION; closed to further
implementation. Successor design:
`../designs/2026-07-17-neutral-held-v1.md`.

Do not reopen this plan or promote its public route. The successor may reuse
the direct free-energy implementation as internal HELD Stage III and must
remove the displaced `two_phase_flash` public surface in the same complete
candidate.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add one local methane/ethane fixed-two-phase `T,P,z` calculation through direct total-free-energy minimization over an installed provider mixture Hessian.

**Architecture:** Two phase amount/volume blocks minimize `sum(A+PV)` under two linear material balances. The provider owns the Helmholtz value, gradient, and Hessian; equilibrium owns variables, scaling, Ipopt, initialization, confirmation, KKT checks, topology rejection, and typed results.

**Tech Stack:** Python 3.13, C++17, pybind11, Ipopt, CMake, scikit-build-core, pytest, installed `epcsaft.native_sdk.v1`.

## Outcome Proof

**Intent:** Preserve the exact implementation and evidence record for the superseded fixed-two-phase candidate without authorizing its revival.
**Current Behavior:** The `two_phase_flash` route and result family are removed; Validation retains the candidate as `NON_ADMISSION` evidence.
**Expected Outcome:** Readers can reproduce the historical candidate and understand why the successor HELD controller displaced it, while runtime stays on the successor surface.
**Target Output:** Historical `TwoPhaseFlashResult` candidate evidence only; no current public output or compatibility alias.
**Owner:** Equilibrium owns the retained formulation history; the neutral HELD design owns the successor controller and public `tp_flash` route.
**Interface:** Historical `two_phase_flash` and its private native solve seam, both retired from current runtime.
**Cutover:** Completed when neutral HELD introduced `tp_flash` and deleted the fixed-route public and native surfaces.
**Replaced Path:** `two_phase_flash`, `TwoPhaseFlashResult`, duplicate diagnostics, native bindings, sources, tests, stubs, exports, and documentation were removed without redirect or alias.
**Evidence:** The frozen implementation history, installed-artifact campaign, `NON_ADMISSION` result, successor negative-space tests, and canonical designs.
**Acceptance Proof:** Current source exposes no retired symbols, preserves the historical evidence, and leaves the accepted pure-saturation capability unchanged.
**Stop Criteria:** Stop on any proposal to reopen this plan, restore a compatibility path, reinterpret `NON_ADMISSION`, or change runtime from historical documentation work.
**Avoid:** Do not revive fixed-phase public APIs, overwrite evidence, modify Provider ownership, or infer global phase stability.
**Risk:** Treating an executed historical checklist as an active plan could restore a scientifically superseded route or erase the reason for its retirement.

## Implementation Boundaries

**Files To Create:** Historical implementation files and evidence are listed task by task; no file creation is authorized now.
**Files To Modify:** Only canonical status or provenance references may change under a separate documentation assignment.
**Files To Avoid:** Current runtime sources, Provider code, accepted receipts, Validation records, artifacts, and unrelated IDE state.
**Source Of Truth:** The superseded fixed-route design, its immutable campaign evidence, and the neutral HELD successor design.
**Read Path:** Git history, retained design and plan, immutable Validation evidence, and current negative-space tests.
**Write Path:** No active implementation path; documentation-only corrections require their own bounded GitHub leaf.
**Integration Points:** Historical Provider mixture SDK, Ipopt TNLP, Python result conversion, and the successor `tp_flash` cutover.
**Migration Or Cutover:** Neutral HELD owns the completed migration away from this route.
**Replaced Path Handling:** Keep every retired symbol deleted and prove absence rather than adding redirects, aliases, or wrappers.
**Acceptance Proof Gate:** Current exports and native sources must show no fixed-route owner while canonical evidence still explains the `NON_ADMISSION` outcome.

## Global Constraints

- Read `AGENTS.md`, `CONTEXT.md`, `ARCHITECTURE.yaml`, the accepted pure-saturation design/receipt, and the canonical Slice 2 design before editing.
- Do not begin implementation until an installed provider wheel exposes the reviewed mixture value/gradient/Hessian SDK tail and the validation repository has frozen the May 2015 source/tolerance contract.
- Correct stale pending-authority prose before describing the new candidate.
- Admit exactly the approved Gross--Sadowski methane/ethane fingerprint and source-audited domain.
- Keep one native module and target; compile no provider source and link no provider implementation symbol.
- Add no density root, third derivative, association, ion, bubble/dew route, TPD, phase discovery, generic route registry, continuation framework, or fallback solver.
- Return no successful collapsed/single-phase state and no globality claim.
- Use exact provider derivatives and independently check directional objective Hessians.

---

### Task 1: Reconcile authority and freeze the public API

**Files:**
- Modify: `AGENTS.md`
- Modify: `CONTEXT.md`
- Modify: `ARCHITECTURE.yaml`
- Modify: `README.md`
- Modify: `src/epcsaft_equilibrium/_api.py`
- Modify: `src/epcsaft_equilibrium/__init__.py`
- Modify: `src/epcsaft_equilibrium/_equilibrium.pyi`
- Create: `tests/test_two_phase_flash.py`

**Interfaces:**
- Consumes: accepted pure-saturation authority from promotion 0018.
- Produces: `two_phase_flash(...) -> TwoPhaseFlashResult` and one extended `PhaseState`.

**Use Cases:**

- Preserve the original public-contract evidence so the later cutover can prove which surface it retired.

- [ ] **Step 1: Correct stale authority wording**

State that the repository is authoritative for the exact accepted pure-saturation scope. Keep global stability, mixtures, and flash closed until a later receipt.

- [ ] **Step 2: Write failing public contract tests**

Freeze this call:

```python
result = two_phase_flash(
    model,
    230.0 * epcsaft.unit_registry.kelvin,
    2.0 * epcsaft.unit_registry.megapascal,
    (0.5, 0.5),
)
```

Require `TwoPhaseFlashResult`, `FlashDiagnostics`, and `FlashError`. Extend `PhaseState` in place to:

```python
@dataclass(frozen=True)
class PhaseState:
    amount_mol: float
    mole_fractions: tuple[float, ...]
    volume_m3: float
    molar_density_mol_m3: float
    pressure_pa: float
    chemical_potential_over_rt: tuple[float, ...]
```

Update pure-saturation tests to expect `(1.0,)` composition and a one-element chemical-potential tuple. Add no compatibility property for the old scalar field.

- [ ] **Step 3: Add input rejection tests**

Reject bare or wrong-dimensional `T/P`, non-scalar/nonfinite values, composition length other than two, nonpositive entries, sums outside `1 +/- 1e-12`, unsupported model type/fingerprint, and source-domain violations. Do not normalize invalid feeds.

- [ ] **Step 4: Run tests and verify RED**

Run: `python -m pytest tests/test_saturation.py tests/test_two_phase_flash.py -q`

Expected: failure because the new types/operation do not exist and pure chemical potential remains scalar.

- [ ] **Step 5: Add the typed Python shell only**

Implement unit/composition/fingerprint validation and payload conversion. Call a not-yet-existing private `_solve_two_phase_flash` so the valid-path test remains RED while rejection tests pass. Export only the named public types/function.

- [ ] **Step 6: Commit the public contract**

```bash
git add AGENTS.md CONTEXT.md ARCHITECTURE.yaml README.md src/epcsaft_equilibrium tests/test_saturation.py tests/test_two_phase_flash.py
git commit -m "test: freeze local two phase flash contract"
```

### Task 2: Extract one shared provider capsule context

**Files:**
- Create: `cpp/src/provider.hpp`
- Create: `cpp/src/provider.cpp`
- Modify: `cpp/src/saturation.hpp`
- Modify: `cpp/src/saturation.cpp`
- Modify: `cpp/src/module.cpp`
- Modify: `CMakeLists.txt`
- Modify: `tests/test_saturation.py`
- Modify: `tests/test_two_phase_flash.py`

**Interfaces:**
- Consumes: accepted pure SDK prefix plus reviewed mixture SDK tail.
- Produces: one `ProviderContext` that validates/retains the capsule and evaluates pure or mixture phase blocks.

**Use Cases:**

- Demonstrate one Provider boundary and detect a duplicate or private Provider implementation path.

- [ ] **Step 1: Add malformed tail tests**

Construct ctypes capsules with wrong table size, result size, component count, null callback, wrong coordinate count, short buffers, nonfinite tensors, wrong fingerprint, and callback errors. Require the native boundary to reject each with a structured message before Ipopt starts.

- [ ] **Step 2: Run transport tests and verify RED**

Run: `python -m pytest tests/test_two_phase_flash.py -k transport -q`

Expected: failure because the mixture tail is not consumed.

- [ ] **Step 3: Move existing provider ownership without behavior change**

Move `ProviderContext` and capsule lifetime/validation from `saturation.hpp/.cpp` to `provider.hpp/.cpp`. Preserve the pure `evaluate(amount,volume,T)` behavior and all pure tests exactly.

- [ ] **Step 4: Add the mixture evaluation method**

Declare:

```cpp
struct MixturePhaseEvaluation final {
    double value;
    std::vector<double> gradient;
    std::vector<double> hessian;
    double pressure_pa;
};

MixturePhaseEvaluation evaluate_mixture(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3
) const;
```

Allocate buffers inside equilibrium, pass them to the provider callback, and verify exact sizes, finite values, status, error text, and fingerprint on every call.

- [ ] **Step 5: Run pure and mixture transport tests**

Run: `python -m pytest tests/test_saturation.py tests/test_two_phase_flash.py -k transport -q`

Expected: all transport tests pass and pure behavior remains unchanged.

- [ ] **Step 6: Commit the shared boundary**

```bash
git add cpp/src/provider.hpp cpp/src/provider.cpp cpp/src/saturation.hpp cpp/src/saturation.cpp cpp/src/module.cpp CMakeLists.txt tests
git commit -m "refactor: share provider capsule context"
```

### Task 3: Implement the exact direct-minimization TNLP

**Files:**
- Create: `cpp/src/two_phase_flash.hpp`
- Create: `cpp/src/two_phase_flash.cpp`
- Modify: `CMakeLists.txt`
- Modify: `tests/test_two_phase_flash.py`

**Interfaces:**
- Consumes: `ProviderContext::evaluate_mixture`.
- Produces: native objective, gradient, constraints, Jacobian, and exact Lagrangian Hessian in six solver coordinates.

**Use Cases:**

- Retain exact derivative proof for the historical local NLP without converting local convergence into phase-discovery evidence.

- [ ] **Step 1: Add a native evaluation hook test**

The private test hook accepts `(T,P,z,variables)` with variables ordered:

```text
(n_L,0, n_L,1, log(V_L), n_V,0, n_V,1, log(V_V))
```

Require two constraints `n_L,i+n_V,i-z_i`, a constant 2x6 Jacobian, finite objective/gradient/Hessian, and exact symmetry.

- [ ] **Step 2: Add independent directional checks**

At an interior point, compare objective-gradient and gradient-Hessian products with centered differences. The log-volume chain rule must be:

```text
g_logV       = (Phi_V + P/(R*T)) * V
H_n,logV     = Phi_nV * V
H_logV,logV  = Phi_VV * V^2 + (Phi_V + P/(R*T)) * V
```

The balance constraint Hessian is exactly zero.

- [ ] **Step 3: Run derivative tests and verify RED**

Run: `python -m pytest tests/test_two_phase_flash.py -k derivative -q`

Expected: failure because the TNLP is absent.

- [ ] **Step 4: Implement the six-variable TNLP**

Use raw positive component amounts so material balances remain linear and log volumes so volume remains positive. Normalize total feed amount to exactly one mole at the API boundary. Bounds are named constants in the native owner; amount lower bounds are `1e-10 mol`, amount upper bounds are the corresponding feed amount, and volume bounds are declared in physical `m3` before conversion to log space.

Objective is dimensionless:

```cpp
value = liquid.value + vapor.value
      + pressure_pa * (liquid_volume + vapor_volume) / (gas_constant * temperature_k);
```

Assemble exact gradient/Hessian from provider blocks and the chain rules above. Use no finite differences in production.

- [ ] **Step 5: Run derivative and structure tests**

Run: `python -m pytest tests/test_two_phase_flash.py -k 'derivative or structure' -q`

Expected: all exact derivative and linear-balance tests pass.

- [ ] **Step 6: Commit the formulation kernel**

```bash
git add cpp/src/two_phase_flash.hpp cpp/src/two_phase_flash.cpp CMakeLists.txt tests/test_two_phase_flash.py
git commit -m "feat: add direct two phase TP objective"
```

### Task 4: Add bounded attempts, confirmation, and local physical acceptance

**Files:**
- Modify: `cpp/src/two_phase_flash.cpp`
- Modify: `cpp/src/module.cpp`
- Modify: `tests/test_two_phase_flash.py`

**Interfaces:**
- Consumes: exact TNLP and fixed feed.
- Produces: native solve payload with separate solver, numerical, and local-physical evidence.

**Use Cases:**

- Make acceptance evidence visible for local states and reject collapsed or uncertified terminals despite solver success.

- [ ] **Step 1: Add success and adversarial tests**

Require deterministic attempts in both binary enrichment directions, a confirmation attempt from perturbed volumes/phase split, and rejection for collapsed phases, active trace bounds, role reversal after density ordering, failed KKT, callback error, unacceptable Ipopt termination, or confirmation disagreement.

- [ ] **Step 2: Run solve tests and verify RED**

Run: `python -m pytest tests/test_two_phase_flash.py -k 'solve or collapse or confirmation' -q`

Expected: failure because solve orchestration is absent.

- [ ] **Step 3: Implement deterministic seeds without paper identity**

For each enrichment sign, form symmetric binary compositions around `z` with
`delta=min(0.20,0.45*min(z_0,1-z_0))`, phase fractions `0.5/0.5`, liquid seed density `20_000 mol/m3`, and vapor seed volume from `nRT/P`. Record every constant and attempt. Try both enrichment signs; do not select by component name.

- [ ] **Step 4: Implement postsolve checks**

Order phases by density only after solve. Require finite positive amounts/volumes, composition normalization, material residual `<=1e-10`, pressure relative residual per phase `<=1e-8`, interior chemical-potential absolute difference `<=1e-8`, inactive amount/log-volume bounds by `1e-7`, relative phase-density distance `>=1e-3`, and KKT stationarity `<=1e-8`. Require the confirmation solution to agree in phase fractions, compositions, volumes, and objective within named `1e-7` relative/absolute tolerances.

- [ ] **Step 5: Expose one private native solve entry**

Bind `_solve_two_phase_flash(capsule,T,P,z)` and return one fixed-shape payload. Do not expose callbacks, variables, or solver options publicly.

- [ ] **Step 6: Run solve tests**

Run: `python -m pytest tests/test_two_phase_flash.py -q`

Expected: local accepted cases pass; all adversarial cases raise with diagnostics; `globality_certificate` is always false.

- [ ] **Step 7: Commit solve and acceptance**

```bash
git add cpp/src/two_phase_flash.cpp cpp/src/module.cpp tests/test_two_phase_flash.py
git commit -m "feat: solve local neutral two phase flash"
```

### Task 5: Complete typed Python results and preserve saturation

**Files:**
- Modify: `src/epcsaft_equilibrium/_api.py`
- Modify: `src/epcsaft_equilibrium/__init__.py`
- Modify: `src/epcsaft_equilibrium/_equilibrium.pyi`
- Modify: `tests/test_saturation.py`
- Modify: `tests/test_two_phase_flash.py`

**Interfaces:**
- Consumes: fixed native payload.
- Produces: immutable public result and structured `FlashError`.

**Use Cases:**

- Preserve the accepted saturation interface while documenting the fixed-route result that the successor later deleted.

- [ ] **Step 1: Implement exact result conversion**

Define `FlashDiagnostics` fields for attempts, bounds, constraint violation, confirmation deltas, material balance, pressure stationarity, chemical-potential equality, KKT norm, phase distance, exact derivatives, failure reason, and `globality_certificate=False`. Missing or wrong-sized native fields raise `RuntimeError`; Python does not infer acceptance.

- [ ] **Step 2: Update pure `PhaseState` conversion**

Return `(1.0,)` and a one-element chemical-potential tuple for saturation. Update docs/tests in place and add no compatibility property.

- [ ] **Step 3: Run both public suites**

Run: `python -m pytest tests/test_saturation.py tests/test_two_phase_flash.py -q`

Expected: all pure and mixture tests pass.

- [ ] **Step 4: Commit the public result**

```bash
git add src/epcsaft_equilibrium tests
git commit -m "feat: expose local two phase flash result"
```

### Task 6: Installed-artifact source campaign and candidate records

**Files:**
- Create: `docs/designs/2026-07-17-neutral-two-phase-tp-flash.md`
- Modify: `README.md`
- Modify: `CONTEXT.md`
- Modify: `ARCHITECTURE.yaml`
- Create: `reviews/2026-07-17-neutral-two-phase-tp-flash.md` only after independent review

**Interfaces:**
- Consumes: immutable provider/equilibrium wheels and validation-owned May 2015 source/tolerance contract.
- Produces: one local methane/ethane candidate subject for migration review.

**Use Cases:**

- Retain installed-artifact and `NON_ADMISSION` evidence while preventing the superseded candidate from becoming a new migration or promotion path.

- [ ] **Step 1: Document exact capability and equations**

Record ensemble, variables, objective, balances, derivative chain rules, seed contract, tolerances, fingerprint/domain, result schema, and explicit no-globality/no-discovery boundary.

- [ ] **Step 2: Build and test installed wheels in isolation**

Build the equilibrium wheel against the installed provider public header. In a clean environment containing only exact provider/equilibrium wheels and test dependencies, run both package suites. Inspect the extension to prove it links Ipopt/platform runtimes only and contains no provider EOS symbol/source path.

- [ ] **Step 3: Run validation-owned real-data campaign**

Invoke the accepted `may_2015_methane_ethane_flash.py` campaign with exact wheel paths. Expected: all preselected rows and predeclared tolerances pass. A failure stops promotion; do not remove rows or change tolerances.

- [ ] **Step 4: Run final negative-space and cleanup checks**

Prove one target, no sibling/lab/migration runtime dependency, no provider source/private import, no third derivatives, density closure, association, ions, TPD, discovery, route registry, or solver selector. Run `git diff --check`, full tests, tracked-file review, and cleanup hook.

- [ ] **Step 5: Commit candidate documentation and stop**

```bash
git add README.md CONTEXT.md ARCHITECTURE.yaml docs/designs/2026-07-17-neutral-two-phase-tp-flash.md
git commit -m "docs: record neutral two phase flash candidate"
```

Report commit/tree/wheel/test/campaign hashes and remaining exclusions. Do not push, promote, add a release, or begin phase discovery.
