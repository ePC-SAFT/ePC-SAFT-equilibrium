# Paper-Faithful HELD2.0 Rewrite Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` (recommended) or
> `superpowers:executing-plans` to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the charged-mixture HELD2 implementation with one
paper-traceable native Steps 1–10 state machine used identically by the C++
diagnostic and public Python API.

**Architecture:** Each Perdomo step has one focused C++ owner and returns a
typed result to a thin `held2_algorithm` state machine. Installed Provider
access, pressure-root evaluation, progress events, and result serialization
remain shared low-level services. The existing charged controller survives
only as a temporary native comparison subject until the new native route
passes the cutover gates; it is then deleted before Python cutover.

**Tech Stack:** C++17, installed `epcsaft::native_sdk`, HiGHS 1.15.1, Ipopt,
NLopt, pybind11, scikit-build-core, Python 3.13, pytest, Ruff, mypy, optional
ASan/UBSan.

## Global Constraints

- Scientific authority is
  `docs/designs/2026-07-24-held2-paper-algorithm.md` at approved commit
  `9e53ee8f04ba10be5bcd70d13ca54a60214fc3bc`.
- Provider thermodynamics, derivatives, packing fraction, model loading, and
  model lifetime are consumed only through the installed Provider SDK.
- All algorithmic energies and potentials use the one-mole reduced
  \(A/(RT)\), \(G/(RT)\), and \(\mu/(RT)\) basis defined by the specification.
- `globality_certificate` remains `not_guaranteed` and
  `phase_enumeration_certificate` remains `completeness_not_guaranteed`.
- The C++ manufactured checker and installed diagnostic are the first
  acceptance routes. Public Python is exercised only after the equivalent
  native check passes.
- Neutral mixtures retain the accepted neutral HELD implementation behind the
  common flash dispatch.
- No same-major two-candidate quota, candidate-origin gate, synthetic
  Appendix-C replacement, log-volume packing proxy, Python scientific gate,
  Provider parser, EOS duplication, or permanent second controller is allowed.
- Trace enablement may not alter starts, evaluations, transitions, results, or
  serialized diagnostics.
- Step 4's completed LP-solve count is the reported Perdomo `Iter.` value.
  Optimizer iterations, starts, and Provider calls are separate counters.
- Resource exhaustion and incomplete evidence are indeterminate, never
  one-phase stability or successful equilibrium.
- Source-qualified Table 5 inputs are required for final cutover. Molality
  alone may not be converted into an overall feed by inventing butanol,
  water, phase-fraction, or midpoint information.
- The existing task-owned working-tree experiment is preserved until its
  useful primitives are either adopted by a reviewed task or removed by the
  final legacy-deletion task. Unrelated changes are never staged.
- Every task ends with native evidence, focused Python evidence when
  applicable, `git diff --check`, and a narrow commit. No task pushes.
- Execution uses `chemical-engineer` for thermodynamic claims,
  `scientific-coding-and-testing` for numerical evidence,
  `ipopt` and `cmake` at their respective boundaries, and
  `minimize-code-surface` before every checkpoint.
- Tasks 5 and 6 additionally require `cutthroat-code-cleanup`: challenge every
  Stage-II start family, retained field, comparison, allocation, solver call,
  and transition, and remove anything not required by Eqs. (64)–(66),
  certification, or declared diagnostics.

---

## File Map

### New production owners

- `cpp/src/held2_algorithm.hpp`
- `cpp/src/held2_algorithm.cpp`
- `cpp/src/held2_step1.hpp`
- `cpp/src/held2_step3.hpp`
- `cpp/src/held2_step3.cpp`
- `cpp/src/held2_step6.hpp`
- `cpp/src/held2_step6.cpp`
- `cpp/src/held2_step7.hpp`
- `cpp/src/held2_step7.cpp`
- `cpp/src/held2_step8.hpp`
- `cpp/src/held2_step8.cpp`
- `cpp/src/held2_step9.hpp`
- `cpp/src/held2_step9.cpp`
- `cpp/src/held2_step10.hpp`
- `cpp/src/held2_step10.cpp`

### Existing focused owners retained and completed

- `cpp/src/held2.hpp`: common value types and low-level evaluator contracts
  only.
- `cpp/src/held2_step1.cpp`: coordinates, bounds, permutation, lift, and
  reduced Provider derivative transformation.
- `cpp/src/held2_step2.hpp`
- `cpp/src/held2_step2.cpp`
- `cpp/src/held2_step4.hpp`
- `cpp/src/held2_step4.cpp`
- `cpp/src/held2_step5.hpp`
- `cpp/src/held2_step5.cpp`
- `cpp/src/held2_tolerances.hpp`
- `cpp/src/held2_progress.hpp`
- `cpp/src/held2_progress.cpp`
- `cpp/src/flash.hpp`
- `cpp/src/flash.cpp`
- `cpp/src/result_json.hpp`
- `cpp/src/result_json.cpp`
- `cpp/src/diagnostic_main.cpp`
- `cpp/src/manufactured_check_main.cpp`
- `cpp/src/module.cpp`
- `CMakeLists.txt`

### Native test support

- Create `cpp/tests/held2_step_checks.hpp`.
- Create `cpp/tests/held2_step_checks.cpp`.
- Create `cpp/tests/held2_workflow_check.cpp`.

The manufactured-check main only parses `--case` and dispatches to these
checks. Test equations and independent oracles do not live in production
step files.

### Python verification

- Modify `tests/test_held2.py`.
- Modify `tests/test_native_diagnostic.py`.
- Modify `tests/test_perdomo_held2_trace.py`.
- Modify `src/epcsaft_equilibrium/_api.py` only in the final cutover task.
- Modify `src/epcsaft_equilibrium/_equilibrium.pyi` only if an existing private
  diagnostic binding changes; add no public solver controls.

### Legacy files deleted after cutover

- `cpp/src/held2_controller.hpp`
- `cpp/src/held2_controller.cpp`
- `cpp/src/held2_stage_iii.cpp`
- `cpp/src/held2_stage_iii_generic.cpp`
- Any superseded orchestration remaining in `cpp/src/held2.cpp`,
  `cpp/src/flash.cpp`, or `cpp/src/module.cpp`.

`cpp/src/held2.cpp` must either be reduced to genuinely shared low-level
numerics or deleted. It may not remain a second algorithm owner.

---

### Task 1: Establish shared contracts and implement Step 1

**Files:**

- Create: `cpp/src/held2_step1.hpp`
- Modify: `cpp/src/held2.hpp`
- Modify: `cpp/src/held2_step1.cpp`
- Modify: `cpp/src/held2_tolerances.hpp`
- Create: `cpp/tests/held2_step_checks.hpp`
- Create: `cpp/tests/held2_step_checks.cpp`
- Modify: `cpp/src/manufactured_check_main.cpp`
- Modify: `CMakeLists.txt`
- Modify: `tests/test_held2.py`

**Interfaces:**

- Consumes: Provider-order component IDs, charges, overall feed, Provider
  packing-fraction volume bounds.
- Produces:

```cpp
using Held2PackingFractionEvaluator = std::function<double(
    const std::vector<double>& physical_fractions,
    double volume
)>;

using Held2PhysicalEvaluator = std::function<Held2PhysicalPhaseBlock(
    const std::vector<double>& physical_fractions,
    double volume
)>;

using Held2PhysicalVolumeBoundsEvaluator =
    std::function<std::array<double, 2>(
        const std::vector<double>& physical_fractions
    )>;

using Held2VolumeBoundsEvaluator =
    std::function<std::array<double, 2>(
        const std::vector<double>& independent_modified_fractions
    )>;

using Held2PressureRootEvaluator = std::function<Held2PressureEnvelopeResult(
    const std::vector<double>& independent_modified_fractions
)>;

struct Held2StepTiming {
    int invocation_count;
    double wall_seconds;
    double cpu_seconds;
    std::uint64_t provider_evaluations;
    std::uint64_t optimizer_solves;
    std::uint64_t optimizer_iterations;
    std::string terminal_status;
    std::string terminal_reason;
    int next_step;
};

struct Held2Step1Result {
    std::string status;
    std::string reason;
    double temperature_k;
    double pressure_pa;
    std::optional<Held2Coordinates> coordinates;
    std::optional<std::vector<double>> independent_feed;
    std::optional<Held2VolumeBoundsEvaluator> volume_bounds;
    Held2StepTiming timing;
};

struct Held2ResourceProfile {
    int step2_search_budget;
    int step5_start_cap;
    int step7_major_iteration_cap;
};

Held2Step1Result run_held2_step1(
    const std::vector<std::string>& component_ids,
    const std::vector<double>& charges,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& physical_feed,
    const Held2PhysicalVolumeBoundsEvaluator& physical_volume_bounds
);
```

- `Held2Coordinates` stores both Provider-to-paper and paper-to-Provider
  permutations, compact-to-paper indices, \(\alpha_i\), full Step-1 polytope,
  and named lower/upper bounds.
- Step 1 validates finite positive \(T\) and \(P_0\), feed normalization,
  charge balance, coordinate rank, and a nonempty Provider-backed volume
  domain. It lifts every modified composition through the one authoritative
  coordinate map before evaluating Provider volume bounds; clipping or a
  synthetic volume interval is forbidden.

- [ ] **Step 1: Add the failing native Step-1 check**

Add `--case step1` with checks for:

```cpp
require(result.coordinates->eliminated_provider_index == chloride_index);
require(result.coordinates->closure_provider_index == butanol_index);
require(max_abs(charge_residuals) <= 1.0e-9);
require(max_abs(round_trip_feed - feed) <= 1.0e-12);
require(result.coordinates->upper_bounds == corrected_eq_59_60_oracle);
require((*result.volume_bounds)(*result.independent_feed)
        == provider_feed_volume_bounds_oracle);
require(singular_same_charge_topology.reason
        == "unsupported_singular_charge_transformation");
require(invalid_temperature.reason == "invalid_temperature");
require(empty_provider_domain.reason == "empty_physical_volume_domain");
```

Cover interleaved Provider order, LiCl, CaCl2, and a deliberately singular
charge topology.

- [ ] **Step 2: Run the native check and confirm RED**

Run:

```bash
cmake --build build --target epcsaft-equilibrium-manufactured-check -j2
build/epcsaft-equilibrium-manufactured-check --case step1
```

Expected: compilation or dispatch failure because the new Step-1 contract is
absent.

- [ ] **Step 3: Implement the minimum Step-1 owner**

Implement the internal permutation, admissible eliminated-ion selection,
closure selection, Eqs. (23), (28), (30), corrected Eqs. (57)–(61), inverse
lift, complete linear polytope, and reduced derivative pullback. Return typed
indeterminate outcomes for invalid feed, singular coordinates, and empty
domains.

- [ ] **Step 4: Run native GREEN before Python**

Run:

```bash
cmake --build build --target epcsaft-equilibrium-manufactured-check -j2
build/epcsaft-equilibrium-manufactured-check --case step1
```

Expected: exit 0 and one `step1: pass` record.

- [ ] **Step 5: Replace obsolete Python coordinate assertions**

Keep one parameterized test covering round-trip, gauge invariance, asymmetric
charge scaling, and typed singular rejection. Do not retain separate tests
for every helper.

Run:

```bash
.venv/bin/python -m pytest tests/test_held2.py -k "coordinates or step1" -q
```

Expected: pass after the native check.

- [ ] **Step 6: Commit**

```bash
git add CMakeLists.txt cpp/src/held2.hpp cpp/src/held2_step1.hpp \
  cpp/src/held2_step1.cpp cpp/src/held2_tolerances.hpp \
  cpp/src/manufactured_check_main.cpp cpp/tests/held2_step_checks.hpp \
  cpp/tests/held2_step_checks.cpp tests/test_held2.py
git commit -m "feat: implement HELD2 paper step 1"
```

---

### Task 2: Implement Step 2 reference selection and TPD search

**Files:**

- Modify: `cpp/src/held2_step2.hpp`
- Modify: `cpp/src/held2_step2.cpp`
- Modify: `cpp/src/held2_progress.hpp`
- Modify: `cpp/src/held2_progress.cpp`
- Modify: `cpp/tests/held2_step_checks.cpp`
- Modify: `tests/test_held2.py`
- Modify: `tests/test_perdomo_held2_trace.py`

**Interfaces:**

- Consumes: `Held2Step1Result` (including its transformed Provider-backed
  volume-domain evaluator), physical Provider evaluator, declared Step-2
  search budget, observer.
- Produces:

```cpp
enum class Held2Step2Outcome {
    NegativeWitness,
    NoNegativeWitnessDetected,
    Indeterminate,
};

struct Held2Step2Result {
    Held2Step2Outcome outcome;
    std::string reason;
    std::string globality_certificate;
    std::optional<Held2PressureEnvelopeResult> reference_envelope;
    std::optional<Held2StateEvaluation> reference;
    std::optional<Held2StageICandidate> negative_witness;
    std::optional<double> minimum_tpd;
    Held2StepTiming timing;
};

Held2Step2Result run_held2_step2(
    const Held2Step1Result& step1,
    const Held2StateEvaluator& evaluator,
    int search_budget,
    Held2ProgressObserver* observer
);
```

Step 2 uses `step1.volume_bounds`; it may not construct or receive a second
volume-domain policy.

- [ ] **Step 1: Add native RED oracles**

The Step-2 case checks:

```cpp
require(abs(tpd_at_reference) <= 1.0e-12);
require_directional_gradient_matches_centered_difference();
require_hessian_vector_matches_centered_gradient_difference();
require(selected_root_is_unique_lowest_strict_stable());
require(negative_witness.tpd < -1.0e-8);
require(no_negative.outcome == Held2Step2Outcome::NoNegativeWitnessDetected);
require(no_negative.globality_certificate == "not_guaranteed");
require(tied_reference.outcome == Held2Step2Outcome::Indeterminate);
```

The independent TPD oracle uses trial-state \(A\), the complete feed tangent,
and verifies that the printed reference-state Eq. (62) would fail the
curvature check.

- [ ] **Step 2: Run native RED**

```bash
cmake --build build --target epcsaft-equilibrium-manufactured-check -j2
build/epcsaft-equilibrium-manufactured-check --case step2
```

- [ ] **Step 3: Implement Step 2**

Enumerate and certify feed pressure roots, reject boundary truncation and
unresolved ties, form Eq. (55), evaluate corrected Eq. (62), run the declared
deterministic finite search, and stop at the first independently certified
\(d<-10^{-8}\) witness. Keep the finite nonnegative outcome operational, not
certified stability.

- [ ] **Step 4: Run native GREEN and trace identity**

```bash
build/epcsaft-equilibrium-manufactured-check --case step2
build/epcsaft-equilibrium-manufactured-check --case step2 --trace \
  > /tmp/held2-step2-trace.txt
```

The trace run must emit the same terminal result hash as the quiet run.

- [ ] **Step 5: Run focused private Python evidence**

```bash
.venv/bin/python -m pytest \
  tests/test_held2.py \
  tests/test_perdomo_held2_trace.py \
  -k "pressure_envelope or stage_i or observer" -q
```

- [ ] **Step 6: Commit**

```bash
git add cpp/src/held2_step2.hpp cpp/src/held2_step2.cpp \
  cpp/src/held2_progress.hpp cpp/src/held2_progress.cpp \
  cpp/tests/held2_step_checks.cpp tests/test_held2.py \
  tests/test_perdomo_held2_trace.py
git commit -m "feat: implement HELD2 paper step 2"
```

---

### Task 3: Implement Step 3 Appendix-C initialization

**Files:**

- Create: `cpp/src/held2_step3.hpp`
- Create: `cpp/src/held2_step3.cpp`
- Modify: `cpp/tests/held2_step_checks.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

```cpp
struct Held2MPoint {
    std::uint64_t insertion_id;
    std::vector<double> independent_modified_fractions;
    double volume;
    double packing_fraction;
    double reduced_helmholtz;
    std::string origin;
};

struct Held2StageIIIFeedback {
    int source_step;
    std::string reason;
    std::vector<std::uint64_t> candidate_ids;
};

struct Held2PersistentState {
    int major_iteration;
    int upper_solve_count;
    double feed_reduced_gibbs;
    double upper_bound;
    double lower_value;
    std::uint64_t next_start_ordinal;
    std::vector<double> multipliers;
    std::vector<Held2MPoint> M;
    std::optional<Held2StageIIIFeedback> stage_iii_feedback;
};

struct Held2Step3Result {
    std::string status;
    std::string reason;
    std::optional<Held2PersistentState> state;
    Held2StepTiming timing;
};

Held2Step3Result run_held2_step3(
    const Held2Step1Result& step1,
    const Held2Step2Result& step2,
    const Held2PressureRootEvaluator& pressure_roots,
    Held2ProgressObserver* observer
);
```

- [ ] **Step 1: Add the Appendix-C RED case**

Assert exactly `1 + 2 * (C - 2)` points, exact corrected C.1/C.2
compositions, lower/upper bracketing, neutral-limit reduction, deterministic
insertion IDs, pressure closure, and a bounded first Eq. (64) LP. Include a
feed for which `x0_i / 2` violates the finite domain and require
`appendix_c_state_outside_domain` without clipping.

- [ ] **Step 2: Run native RED**

```bash
cmake --build build --target epcsaft-equilibrium-manufactured-check -j2
build/epcsaft-equilibrium-manufactured-check --case step3
```

- [ ] **Step 3: Implement corrected Appendix C**

Implement the \(1/2\) C.2 prefactor and \(z_i\) distinguished index, validate
the complete Step-1 polytope before pressure evaluation, choose only a unique
lowest-objective strict-stable root, and initialize \(k\), `upper_solve_count`,
`next_start_ordinal`, \(UBD^V\), and \(\bar L^V=-\infty\).

- [ ] **Step 4: Run native GREEN**

```bash
build/epcsaft-equilibrium-manufactured-check --case step3
```

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt cpp/src/held2_step3.hpp cpp/src/held2_step3.cpp \
  cpp/tests/held2_step_checks.cpp
git commit -m "feat: implement HELD2 paper step 3"
```

---

### Task 4: Implement Step 4 certified upper LP

**Files:**

- Modify: `cpp/src/held2_step4.hpp`
- Modify: `cpp/src/held2_step4.cpp`
- Modify: `cpp/src/held2_tolerances.hpp`
- Modify: `cpp/tests/held2_step_checks.cpp`

**Interfaces:**

```cpp
struct Held2LpCertificate {
    bool primal_feasible;
    bool dual_feasible;
    double primal_residual_inf;
    double dual_residual_inf;
    double complementarity_inf;
};

struct Held2Step4Result {
    std::string status;
    std::string reason;
    std::optional<double> upper_bound;
    std::optional<std::vector<double>> multipliers;
    std::vector<int> active_cut_ids;
    std::optional<Held2LpCertificate> certificate;
    Held2StepTiming timing;
};

Held2Step4Result run_held2_step4(
    Held2PersistentState& state,
    Held2ProgressObserver* observer
);
```

- [ ] **Step 1: Add native analytic-envelope RED checks**

Compare HiGHS Eq. (64) to an independently enumerated low-dimensional vertex
oracle. Require primal feasibility, dual feasibility, complementarity, active
cuts, and exactly one increment of `upper_solve_count`. Require nonoptimal,
unbounded, or uncertified status to be indeterminate.

- [ ] **Step 2: Run native RED**

```bash
cmake --build build --target epcsaft-equilibrium-manufactured-check -j2
build/epcsaft-equilibrium-manufactured-check --case step4
```

- [ ] **Step 3: Implement the minimum certified LP**

Use all persistent \(\mathcal M\) cuts, the feed Gibbs cap, and unrestricted
mathematical multipliers. Do not add synthetic box bounds. Preserve the full
HiGHS certificate. Increment `upper_solve_count` exactly once, and emit paper
`Iter.`, only after HiGHS returns a certified optimal Eq. (64) solve; failed
or uncertified attempts do not increment it.

- [ ] **Step 4: Run native GREEN**

```bash
build/epcsaft-equilibrium-manufactured-check --case step4
```

- [ ] **Step 5: Commit**

```bash
git add cpp/src/held2_step4.hpp cpp/src/held2_step4.cpp \
  cpp/src/held2_tolerances.hpp cpp/tests/held2_step_checks.cpp
git commit -m "feat: implement HELD2 paper step 4"
```

---

### Task 5: Implement Step 5 full-space lower NLP

**Files:**

- Modify: `cpp/src/held2_step5.hpp`
- Modify: `cpp/src/held2_step5.cpp`
- Modify: `cpp/src/held2.hpp`
- Modify: `cpp/tests/held2_step_checks.cpp`

**Interfaces:**

```cpp
struct Held2LocalCertificate {
    std::uint64_t start_ordinal;
    std::string solver_status;
    bool finite_and_in_domain;
    double pressure_residual;
    double primal_residual_inf;
    double stationarity_residual_inf;
    double dual_sign_violation_inf;
    double complementarity_inf;
    double dual_pullback_residual_inf;
    bool accepted;
};

struct Held2Step5Result {
    std::string status;
    std::string reason;
    std::optional<double> lower_value;
    std::optional<Held2MPoint> terminal;
    std::uint64_t starts_consumed;
    std::vector<Held2LocalCertificate> attempts;
    Held2StepTiming timing;
};

Held2Step5Result run_held2_step5(
    const Held2Step1Result& step1,
    const Held2Step4Result& step4,
    Held2PersistentState& state,
    const Held2StateEvaluator& evaluator,
    const Held2PackingFractionEvaluator& packing_fraction,
    const Held2ResourceProfile& resources,
    Held2ProgressObserver* observer
);
```

- [ ] **Step 1: Add native RED checks**

Check exact objective gradient and Hessian in \(C-2\) modified coordinates
plus volume; original-coordinate KKT, pressure, sign, complementarity, and
dual pullback; best-certified retention; stop on the first best
\(\bar L^{V,k}\le UBD^V\); representation-equivalent insertion; and persistent
start ordinals and epoch state across majors. Failure and resource-exhaustion
results must have `terminal == std::nullopt`; Step 7 may consume a terminal
only after checking the successful Step-5 status.

- [ ] **Step 2: Run native RED**

```bash
cmake --build build --target epcsaft-equilibrium-manufactured-check -j2
build/epcsaft-equilibrium-manufactured-check --case step5
```

- [ ] **Step 3: Implement direct Eq. (65)**

Use one full-space exact-Hessian Ipopt TNLP. Remove the global-basin
controller, same-major candidate quota, and post-stop second-basin search.
The deterministic stream includes declared random interior,
shifted-previous, and near-pure starts and never resets. Solver success alone
does not qualify a terminal.

- [ ] **Step 4: Run native GREEN and bounded replay**

```bash
build/epcsaft-equilibrium-manufactured-check --case step5
build/epcsaft-equilibrium-manufactured-check --case step5 --repeat 2
```

Both repeats must have identical terminals, work counts, and serialized
certificates.

- [ ] **Step 5: Commit**

```bash
git add cpp/src/held2.hpp cpp/src/held2_step5.hpp \
  cpp/src/held2_step5.cpp cpp/tests/held2_step_checks.cpp
git commit -m "feat: implement HELD2 paper step 5"
```

---

### Task 6: Implement Steps 6 and 7 persistent candidate search

**Files:**

- Create: `cpp/src/held2_step6.hpp`
- Create: `cpp/src/held2_step6.cpp`
- Create: `cpp/src/held2_step7.hpp`
- Create: `cpp/src/held2_step7.cpp`
- Modify: `cpp/src/held2_tolerances.hpp`
- Modify: `cpp/tests/held2_step_checks.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

```cpp
struct Held2CandidateDecision {
    std::uint64_t insertion_id;
    bool gap_passed;
    bool derivative_passed;
    bool pairwise_distinct;
    bool retained;
    std::string reason;
};

struct Held2Step6Result {
    std::string status;
    std::string reason;
    std::vector<Held2MPoint> candidates;
    std::vector<Held2CandidateDecision> decisions;
    Held2StepTiming timing;
};

struct Held2Step7Result {
    std::string status;
    std::string reason;
    std::optional<int> next_step;
    Held2StepTiming timing;
};

Held2Step6Result run_held2_step6(
    const Held2Step1Result& step1,
    const Held2Step4Result& step4,
    const Held2PersistentState& state,
    const Held2PackingFractionEvaluator& packing_fraction,
    Held2ProgressObserver* observer
);

Held2Step7Result run_held2_step7(
    Held2PersistentState& state,
    const Held2Step5Result& step5,
    const Held2Step6Result& step6,
    const Held2ResourceProfile& resources,
    Held2ProgressObserver* observer
);
```

- [ ] **Step 1: Add native RED cases**

Require all historical \(\mathcal M\) members to be re-evaluated under current
multipliers, all non-bound independent coordinates to use Eq. (66), Provider
packing fraction to distinguish branches, and a greedy maximal
pairwise-distinct candidate set. Include the nontransitive A-near-B,
B-near-C, A-distinct-C case. Verify `k` increments once, \(\mathcal M\)
persists, and start ordinals persist.

- [ ] **Step 2: Run native RED**

```bash
cmake --build build --target epcsaft-equilibrium-manufactured-check -j2
build/epcsaft-equilibrium-manufactured-check --case step6
build/epcsaft-equilibrium-manufactured-check --case step7
```

- [ ] **Step 3: Implement Steps 6 and 7**

Use \(\epsilon_b=10^{-2}\), \(\epsilon_\lambda=0.5\),
\(\epsilon_\eta=10^{-3}\), and \(\epsilon_x=10^{-3}\) as explicit project
settings with diagnostics. Do not use existing tighter basin tolerances as
paper Eq. (66). Make resource and stagnation exits typed indeterminate.

- [ ] **Step 4: Run native GREEN**

```bash
build/epcsaft-equilibrium-manufactured-check --case step6
build/epcsaft-equilibrium-manufactured-check --case step7
```

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt cpp/src/held2_step6.hpp cpp/src/held2_step6.cpp \
  cpp/src/held2_step7.hpp cpp/src/held2_step7.cpp \
  cpp/src/held2_tolerances.hpp cpp/tests/held2_step_checks.cpp
git commit -m "feat: implement HELD2 paper steps 6 and 7"
```

---

### Task 7: Implement Step 8 total-Gibbs solve and active-set lifecycle

**Files:**

- Create: `cpp/src/held2_step8.hpp`
- Create: `cpp/src/held2_step8.cpp`
- Modify: `cpp/tests/held2_step_checks.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

```cpp
enum class Held2Step8Outcome {
    CertifiedFeasible,
    CertifiedInfeasible,
    Indeterminate,
};

struct Held2Phase {
    std::uint64_t stable_id;
    double phase_fraction;
    std::vector<double> independent_modified_fractions;
    std::vector<double> physical_fractions_provider_order;
    double volume;
    double packing_fraction;
};

struct Held2FeasibilityCertificate {
    std::string solver_status;
    bool feasible;
    bool infeasible;
    bool farkas_certificate_valid;
    double primal_residual_inf;
    double certificate_residual_inf;
};

struct Held2NlpCertificate {
    std::string solver_status;
    double primal_residual_inf;
    double stationarity_residual_inf;
    double dual_sign_violation_inf;
    double complementarity_inf;
    bool accepted;
};

struct Held2LifecycleDecision {
    std::uint64_t stable_id;
    std::string action;
    std::string reason;
    bool reduced_resolve_accepted;
};

struct Held2Step8Result {
    Held2Step8Outcome outcome;
    std::string reason;
    std::optional<double> total_reduced_gibbs;
    std::vector<Held2Phase> active_phases;
    std::optional<Held2FeasibilityCertificate> feasibility;
    std::optional<Held2NlpCertificate> nlp;
    std::vector<Held2LifecycleDecision> lifecycle;
    Held2StepTiming timing;
};

Held2Step8Result run_held2_step8(
    const Held2Step1Result& step1,
    const Held2Step6Result& step6,
    const Held2StateEvaluator& evaluator,
    const Held2PackingFractionEvaluator& packing_fraction,
    Held2ProgressObserver* observer
);
```

- [ ] **Step 1: Add native RED feasibility and lifecycle cases**

Check the exact perspective LP against independent convex-hull cases,
including a validated Farkas certificate. Require Ipopt failure after a
feasible LP to be indeterminate. Check duplicate merge initialization,
one-at-a-time KKT retirement, reduced Problem-(67) re-solves, balance
preservation, and feedback when fewer than two active phases remain.

- [ ] **Step 2: Run native RED**

```bash
cmake --build build --target epcsaft-equilibrium-manufactured-check -j2
build/epcsaft-equilibrium-manufactured-check --case step8
```

- [ ] **Step 3: Implement Step 8**

Implement the exact perspective feasibility LP over the complete Step-1
polytope, then the general \(mp\)-phase exact-Hessian Ipopt Problem (67).
Use stable-ID greedy duplicate decisions and the named physical
composition-plus-packing-fraction identity. Re-solve after every merge or
retirement; never accept arithmetic averaging as final equilibrium.

- [ ] **Step 4: Run native GREEN**

```bash
build/epcsaft-equilibrium-manufactured-check --case step8
```

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt cpp/src/held2_step8.hpp cpp/src/held2_step8.cpp \
  cpp/tests/held2_step_checks.cpp
git commit -m "feat: implement HELD2 paper step 8"
```

---

### Task 8: Implement Step 9 convergence and Stage-II feedback

**Files:**

- Create: `cpp/src/held2_step9.hpp`
- Create: `cpp/src/held2_step9.cpp`
- Modify: `cpp/tests/held2_step_checks.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

```cpp
enum class Held2Step9Outcome {
    Converged,
    PaperConvergenceFailed,
    Indeterminate,
};

struct Held2PotentialComparison {
    std::size_t component_index;
    std::uint64_t left_phase_id;
    std::uint64_t right_phase_id;
    double numerator;
    double denominator;
    double ratio;
    bool passed;
};

struct Held2PhysicalCertificate {
    double modified_balance_inf;
    double ordinary_balance_inf;
    double electroneutrality_inf;
    double pressure_residual_inf;
    double kkt_residual_inf;
    bool accepted;
};

struct Held2Step9Result {
    Held2Step9Outcome outcome;
    std::string reason;
    std::optional<double> free_energy_gap;
    std::vector<Held2PotentialComparison> potential_comparisons;
    std::optional<Held2PhysicalCertificate> physical;
    Held2StepTiming timing;
};

Held2Step9Result run_held2_step9(
    const Held2Step4Result& step4,
    const Held2Step8Result& step8,
    const Held2StateEvaluator& evaluator,
    Held2ProgressObserver* observer
);
```

- [ ] **Step 1: Add native RED checks**

Exercise Eq. (68) pass/fail, Eq. (69) pass/fail, exact zero denominator with
zero numerator, exact zero denominator with nonzero numerator, deterministic
phase ordering, and separation of paper feedback from auxiliary-certificate
indeterminate failure.

- [ ] **Step 2: Run native RED**

```bash
cmake --build build --target epcsaft-equilibrium-manufactured-check -j2
build/epcsaft-equilibrium-manufactured-check --case step9
```

- [ ] **Step 3: Implement Step 9**

Use all \(C-1\) retained modified potentials from Eqs. (28), (32)–(34).
Return to Step 4 only for failed Eqs. (68) or (69). Provider, solver, balance,
charge, pressure, KKT, active-phase, or finite-state certificate failures are
indeterminate.

- [ ] **Step 4: Run native GREEN**

```bash
build/epcsaft-equilibrium-manufactured-check --case step9
```

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt cpp/src/held2_step9.hpp cpp/src/held2_step9.cpp \
  cpp/tests/held2_step_checks.cpp
git commit -m "feat: implement HELD2 paper step 9"
```

---

### Task 9: Implement Step 10 logarithmic trace refinement

**Implemented clarification:** the finite Step-1 floor is a multidimensional
search regularization. Steps 1--9 remain in linear modified-composition and
balance coordinates. Step 10 alone may evaluate a strictly positive charged
trace coordinate below that floor in \(\log_{10}x_i\), while retaining every
other Step-1 polytope constraint. Controller transitions use typed actions;
diagnostic reason strings do not own routing.

**Files:**

- Create: `cpp/src/held2_step10.hpp`
- Create: `cpp/src/held2_step10.cpp`
- Modify: `cpp/tests/held2_step_checks.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

```cpp
struct Held2TraceRefinement {
    std::uint64_t phase_id;
    std::size_t component_index;
    std::uint64_t reference_phase_id;
    double initial_mole_fraction;
    double refined_mole_fraction;
    double final_potential_residual;
    std::string status;
};

struct Held2Step10Result {
    std::string status;
    std::string reason;
    std::vector<Held2Phase> phases;
    std::vector<Held2TraceRefinement> refinements;
    std::optional<Held2PhysicalCertificate> final_certificate;
    Held2StepTiming timing;
};

Held2Step10Result run_held2_step10(
    const Held2Step1Result& step1,
    const Held2Step8Result& step8,
    const Held2Step9Result& step9,
    const Held2StateEvaluator& evaluator,
    Held2ProgressObserver* observer
);
```

- [ ] **Step 1: Add native RED trace cases**

Check no-op without trace components, Eq. (71) use of \(z_j\), log bounds
\([10^{-300},5\times10^{-10}]\), greatest-composition reference selection,
stable-ID tie break, zero-residual solve, absent-reference failure, and
post-refinement balance/certificate rejection.

- [ ] **Step 2: Run native RED**

```bash
cmake --build build --target epcsaft-equilibrium-manufactured-check -j2
build/epcsaft-equilibrium-manufactured-check --case step10
```

- [ ] **Step 3: Implement Step 10**

Solve the bounded one-dimensional chemical-potential residual. Do not
minimize the signed residual toward a bound and do not invent a material
balance correction. Reconstruct Provider-order physical phases and re-run the
complete certificate after all refinements.

- [ ] **Step 4: Run native GREEN**

```bash
build/epcsaft-equilibrium-manufactured-check --case step10
```

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt cpp/src/held2_step10.hpp cpp/src/held2_step10.cpp \
  cpp/tests/held2_step_checks.cpp
git commit -m "feat: implement HELD2 paper step 10"
```

---

### Task 10: Assemble the thin algorithm and native diagnostic

**Files:**

- Create: `cpp/src/held2_algorithm.hpp`
- Create: `cpp/src/held2_algorithm.cpp`
- Create: `cpp/tests/held2_workflow_check.cpp`
- Modify: `cpp/src/flash.hpp`
- Modify: `cpp/src/flash.cpp`
- Modify: `cpp/src/diagnostic_main.cpp`
- Modify: `cpp/src/result_json.hpp`
- Modify: `cpp/src/result_json.cpp`
- Modify: `cpp/src/held2_progress.hpp`
- Modify: `cpp/src/held2_progress.cpp`
- Modify: `cpp/src/manufactured_check_main.cpp`
- Modify: `CMakeLists.txt`
- Modify: `tests/test_native_diagnostic.py`

**Interfaces:**

```cpp
struct Held2Input {
    double temperature_k;
    double pressure_pa;
    std::vector<double> overall_mole_fractions_provider_order;
};

struct Held2ThermodynamicAccess {
    std::vector<std::string> component_ids;
    std::vector<double> charges;
    Held2PhysicalEvaluator evaluate_physical;
    Held2PhysicalVolumeBoundsEvaluator volume_bounds_physical;
    Held2PackingFractionEvaluator packing_fraction;
};

struct Held2AlgorithmResult {
    std::string outcome;
    std::string failure_stage;
    std::string failure_reason;
    std::string globality_certificate;
    std::string phase_enumeration_certificate;
    Held2Step1Result step1;
    std::optional<Held2Step2Result> step2;
    std::optional<Held2Step3Result> step3;
    std::vector<Held2Step4Result> step4_history;
    std::vector<Held2Step5Result> step5_history;
    std::vector<Held2Step6Result> step6_history;
    std::vector<Held2Step7Result> step7_history;
    std::vector<Held2Step8Result> step8_history;
    std::vector<Held2Step9Result> step9_history;
    std::optional<Held2Step10Result> step10;
    std::optional<Held2PersistentState> final_state;
    std::vector<Held2Phase> phases;
    std::vector<Held2StepTiming> step_timings;
    int upper_solve_count;
};

Held2AlgorithmResult run_held2_algorithm(
    const Held2ThermodynamicAccess& thermodynamics,
    const Held2Input& input,
    const Held2ResourceProfile& resources,
    Held2ProgressObserver* observer
);

Held2ThermodynamicAccess make_installed_held2_access(
    const ProviderContext& provider,
    const Held2Input& input
);
```

`held2_algorithm.cpp` contains only the Step 1→10 transitions from the
approved pseudocode. It does not implement equations, solvers, serializers,
or Provider access. The result retains the exact root-selection evidence,
TPD witness/search outcome, every LP and local-NLP certificate, Stage-II
candidate decisions, Step-8 feasibility/KKT/lifecycle evidence, Step-9
potential comparisons, Step-10 refinements, final physical certificate, and
persistent-state counters needed for lossless native JSON and Python parity.
Serializers may omit internal callable evaluators, but may not recompute or
reclassify scientific evidence. Every invocation-dependent field is optional:
absence means the computation was not reached or did not produce that
evidence. Implementations may not substitute NaNs, zero-filled certificates,
empty default states, or fabricated solver records.

- [ ] **Step 1: Add full native workflow RED**

Add a table-driven transition matrix to the native workflow checker. It must
cover:

- Step 2 `NoNegativeWitnessDetected` to the operational one-phase result,
  with both non-guarantee labels;
- every Step 1–6 and Step 8–10 fail-closed/indeterminate terminal;
- successful Step 6 to Step 7, including persistent \(\mathcal M\), start
  epoch, and ordinal state;
- Step 7 resource exhaustion and completed-epoch stagnation;
- Step 8 certified infeasibility, active-set reduction, and feedback to
  Step 4;
- Step 9 paper-convergence feedback to Step 4;
- Step 10 post-refinement certificate rejection; and
- one complete Steps 1–10 accepted multiphase result.

The installed diagnostic must emit real-time step-start/step-end events and
structured per-step wall/CPU/work counts for every invoked step. Each matrix
case asserts the terminal outcome, exact next-step history, retained
evidence, and non-guarantee labels.

- [ ] **Step 2: Run native RED**

```bash
cmake --build build --target \
  epcsaft-equilibrium-manufactured-check \
  epcsaft-equilibrium-diagnostic -j2
build/epcsaft-equilibrium-manufactured-check --case workflow
build/epcsaft-equilibrium-manufactured-check --case workflow-transitions
```

- [ ] **Step 3: Implement the thin controller and serializer**

Have `diagnostic_main.cpp` call `run_held2_algorithm` directly through a
shared Provider-backed `Held2ThermodynamicAccess` factory. Keep
`_equilibrium._solve_tp_flash` and its public Python caller on the legacy
charged engine at this checkpoint. The temporary old/new comparison call is
internal to the native checker and has no public or installed runtime switch.

- [ ] **Step 4: Run installed diagnostic GREEN before Python**

Export a Provider model through `epcsaft.export_native_model`, then run:

```bash
build/epcsaft-equilibrium-diagnostic \
  --model-config /tmp/held2-provider-model.json \
  --temperature 298.15 \
  --pressure 2508.0 \
  --feed 0.083947482323093472,0.083947482323093472,0.83210503535381308 \
  --trace \
  --output /tmp/held2-native-result.json
```

Require deterministic repeated JSON, trace/result identity, Provider
fingerprint, Steps 1–10 timing records, and no Provider re-evaluation during
serialization. Compare the serialized root, TPD, LP, local-certificate,
candidate, lifecycle, convergence, refinement, and final physical-certificate
blocks to the in-memory `Held2AlgorithmResult`; no status-only JSON is
accepted.

- [ ] **Step 5: Run the Python-owned CLI transport smoke test**

```bash
EPCSAFT_EQUILIBRIUM_DIAGNOSTIC="$PWD/build/epcsaft-equilibrium-diagnostic" \
  .venv/bin/python -m pytest \
  tests/test_native_diagnostic.py -k "streams_trace or provider_lifetime" -q
```

This launches the native executable from pytest to check model export,
Provider handle lifetime, trace destination, and JSON parsing. Scientific
native/Python parity is intentionally deferred until both routes use the new
engine in Task 12.

- [ ] **Step 6: Commit**

```bash
git add CMakeLists.txt cpp/src/held2_algorithm.hpp \
  cpp/src/held2_algorithm.cpp cpp/src/flash.hpp cpp/src/flash.cpp \
  cpp/src/diagnostic_main.cpp cpp/src/result_json.hpp \
  cpp/src/result_json.cpp cpp/src/held2_progress.hpp \
  cpp/src/held2_progress.cpp cpp/src/manufactured_check_main.cpp \
  cpp/tests/held2_workflow_check.cpp tests/test_native_diagnostic.py
git commit -m "feat: assemble native HELD2 steps 1 through 10"
```

---

### Task 11: Qualify Table 5 inputs and native performance

**Files:**

- Modify: `tests/test_native_diagnostic.py`
- Create only when source-qualified data exist:
  `tests/data/perdomo_table5_inputs.json`
- Modify: `docs/phase-equilibrium.md`

**Interfaces:**

Each Table 5 record contains this typed schema:

| field | type and required content |
|---|---|
| `molality_mol_per_kg_water` | finite positive number matching one published row |
| `overall_mole_fractions` | four finite nonnegative numbers in declared component order |
| `component_order` | exactly Li+, Cl−, water, 1-butanol Provider IDs |
| `source_classification` | `source_reported` or `uniquely_reconstructed` |
| `source_citation` | article, table/page, and every auxiliary density source |
| `water_molar_mass_g_per_mol` | exact molar mass used in reconstruction |
| `density_inputs` | every value, unit, temperature, and citation used |
| `recomputed_molality_mol_per_kg_water` | independently recomputed value |
| `normalization_residual` | signed residual from one |
| `electroneutrality_residual` | signed charge residual |

No record is created until every field has actual source-qualified content.

- [ ] **Step 1: Resolve feed provenance before execution**

Accept only source-reported or uniquely reconstructed overall feeds for 4.58,
4.95, and 5.74 mol/kg. If butanol amount or stock density remains unknown,
record the gate as blocked and stop this task. Do not create synthetic
acceptance data.

- [ ] **Step 2: Add one parameterized native campaign**

For each qualified feed, require native Steps 1–10 completion, accepted
two-phase output, conservation/certificates, deterministic replay, and
structured per-step timing. Record Perdomo phase endpoints, `Iter.`, and CPU
as comparison fields, not exact ePC-SAFT acceptance values.

- [ ] **Step 3: Profile unexplained work**

Use diagnostic counters to prove:

```text
paper_iter == completed_step4_lp_solves
serializer_provider_evaluations == 0
replayed_optimizer_iterations == 0
step5_start_ordinals_are_unique == true
unchanged_M_epochs_are_bounded == true
```

Any multi-minute case must identify the responsible step, solver calls,
iterations, and Provider evaluations before cutover.

- [ ] **Step 4: Run the native campaign**

```bash
EPCSAFT_EQUILIBRIUM_DIAGNOSTIC="$PWD/build/epcsaft-equilibrium-diagnostic" \
  .venv/bin/python -m pytest \
  tests/test_native_diagnostic.py -k "table5" -q
```

- [ ] **Step 5: Commit only when all three records are qualified**

```bash
git add tests/data/perdomo_table5_inputs.json \
  tests/test_native_diagnostic.py docs/phase-equilibrium.md
git commit -m "test: qualify HELD2 Table 5 native cases"
```

If provenance is blocked, do not make a false completion commit; report the
missing source facts in Issue #67.

---

### Task 12: Cut over Python, delete the legacy controller, and verify wheels

**Files:**

- Modify: `cpp/src/flash.cpp`
- Modify: `cpp/src/flash.hpp`
- Modify: `cpp/src/module.cpp`
- Modify: `cpp/src/result_json.cpp`
- Modify: `src/epcsaft_equilibrium/_api.py`
- Modify: `tests/test_held2.py`
- Modify: `tests/test_native_diagnostic.py`
- Modify: `tests/test_perdomo_held2_trace.py`
- Modify: `ARCHITECTURE.yaml`
- Modify: `CMakeLists.txt`
- Delete: `cpp/src/held2_controller.hpp`
- Delete: `cpp/src/held2_controller.cpp`
- Delete: `cpp/src/held2_stage_iii.cpp`
- Delete: `cpp/src/held2_stage_iii_generic.cpp`
- Delete or reduce: `cpp/src/held2.cpp`

**Interfaces:**

- Public `epcsaft_equilibrium.tp_flash(...)` remains unchanged.
- Python converts the native result to typed public objects but does not
  reclassify phases, impose component-count/charge gates, or recompute
  scientific values.

- [ ] **Step 1: Switch the shared charged dispatch**

Route both the diagnostic and `_equilibrium._solve_tp_flash` through
`run_held2_algorithm`. Remove the temporary native comparison call.

- [ ] **Step 2: Remove Python scientific postprocessing**

Delete the phase-count, neutral/charged, midpoint, or component-specific
scientific gates from `_api.py`. Retain only schema validation, unit
conversion, typed object construction, and fail-closed malformed-payload
handling.

- [ ] **Step 3: Prove native first, then public Python**

Run:

```bash
build/epcsaft-equilibrium-manufactured-check --case workflow
EPCSAFT_EQUILIBRIUM_DIAGNOSTIC="$PWD/build/epcsaft-equilibrium-diagnostic" \
  .venv/bin/python -m pytest tests/test_native_diagnostic.py -q
.venv/bin/python -m pytest \
  tests/test_held2.py tests/test_perdomo_held2_trace.py -q
```

Require field-for-field native/Python scientific payload parity and trace
on/off identity.

- [ ] **Step 4: Delete every legacy scientific owner**

Use `rg` to prove no old controller symbol, Stage-III owner, same-major quota,
origin gate, synthetic Appendix-C point, or legacy result vocabulary remains.
Do not retain compatibility aliases.

- [ ] **Step 5: Run full native, sanitizer, Python, and static verification**

```bash
cmake --build build -j2
.venv/bin/python -m pytest -q
.venv/bin/python -m ruff check .
.venv/bin/python -m mypy src/epcsaft_equilibrium

cmake -S . -B build-sanitize \
  -DEPCSAFT_EQUILIBRIUM_BUILD_DIAGNOSTIC=ON \
  -DEPCSAFT_EQUILIBRIUM_ENABLE_SANITIZERS=ON
cmake --build build-sanitize \
  --target epcsaft-equilibrium-manufactured-check -j2
build-sanitize/epcsaft-equilibrium-manufactured-check --case workflow
```

- [ ] **Step 6: Build and validate the exact wheel**

Build a non-editable Equilibrium wheel, install it beside the exact Provider
wheel in an isolated environment, verify imports resolve only from
`site-packages`, export the Provider model through its public API, and run the
native diagnostic before public `tp_flash`.

- [ ] **Step 7: Run minimality and cleanup audits**

```bash
git diff --check
bash "$HOME/.codex/hooks/codex-cleanup.sh" --repo-root .
```

Review `module.cpp`, `flash.cpp`, and every step file for one owner per
responsibility. Remove obsolete tests and temporary campaign artifacts.

- [ ] **Step 8: Commit**

```bash
git add -- ARCHITECTURE.yaml CMakeLists.txt \
  cpp/src/flash.cpp cpp/src/flash.hpp cpp/src/held2.cpp \
  cpp/src/held2_algorithm.cpp cpp/src/held2_algorithm.hpp \
  cpp/src/held2_controller.cpp cpp/src/held2_controller.hpp \
  cpp/src/held2_stage_iii.cpp cpp/src/held2_stage_iii_generic.cpp \
  cpp/src/module.cpp cpp/src/result_json.cpp \
  src/epcsaft_equilibrium/_api.py \
  tests/test_held2.py tests/test_native_diagnostic.py \
  tests/test_perdomo_held2_trace.py
git commit -m "feat: cut over to paper-faithful HELD2"
```

Do not push, close Issue #67, or claim promotion until separately authorized
remote review and CI complete.

---

## Plan Self-Review

### Specification coverage

- Steps 1–10 each have one production owner and one native acceptance case.
- Eq. (62), Eqs. (59)–(60), Appendix C, and Eq. (71) errata are assigned to
  Tasks 1–3 and 9.
- Persistent \(\mathcal M\), stateful starts, all-\(\mathcal M\) Step 6,
  Step-7 iteration semantics, and stagnation are assigned to Tasks 3–6.
- Certified Step-8 feasibility, duplicate handling, inactive retirement,
  Eq. (68), Eq. (69), and trace refinement are assigned to Tasks 7–9.
- Thin orchestration, live trace, timing, JSON parity, Table 5 provenance,
  public cutover, and legacy deletion are assigned to Tasks 10–12.
- Native-first validation precedes private Python checks in every task and
  public Python cutover in Task 12.

### Test-retention policy

The rewrite replaces obsolete HELD2 assertions rather than accumulating one
pytest function per helper. Durable Python coverage remains concentrated in:

1. one parameterized scientific Step 1–10 contract test;
2. one native/Python payload and trace parity test; and
3. one source-qualified Table 5 campaign.

Fine-grained equations and transition cases live in the native manufactured
checker and are invoked by explicit `--case` names.

### Stop conditions

Stop immediately on:

- an unresolved paper equation or component-index ambiguity;
- an installed Provider contract gap;
- a source-unqualified Table 5 overall feed;
- a failed independent derivative, LP, KKT, balance, pressure, potential, or
  feasibility certificate;
- trace changing a result;
- native/Python scientific payload disagreement;
- a need to preserve two production charged controllers;
- an unexplained long-running step or replayed diagnostic work; or
- any request to push, merge, publish, promote, or close governance without
  separate authorization.
