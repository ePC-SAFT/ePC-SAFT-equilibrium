# Neutral HELD v1 Design

Status: executed v1 candidate frozen at implementation commit
`8318e755d4a8e490822fdf7bb2685d8c5af6436c`; Validation commit
`93ff18541d2fe277a27671e4e6d12b6b009a58ed` is stable `NON_ADMISSION`
evidence; the design-only controller delta in Section 13 is deferred
non-production provenance after the user-authorized Perdomo HELD2 pivot and
authorizes no runtime implementation

Approved in discussion: 2026-07-17

Authority effect: none

Implementation plan: `../plans/2026-07-17-neutral-held-v1-plan.md`

## 1. Purpose

`neutral-held-v1` is the next equilibrium milestone. It implements the three
stages and feedback loop of Pereira's Helmholtz free Energy Lagrangian Dual
(HELD) algorithm for one bounded neutral methane/ethane `T,P,z` capability.
The operation discovers a one- or two-phase result without a user-supplied
phase count or phase guesses.

This is a clean extraction, not a migration of the lab equilibrium framework.
It reuses the installed provider's Helmholtz kernel and the useful direct
free-energy mathematics in the unpromoted fixed-two-phase candidate. It does
not copy provider equations, add a second solver kernel, or claim a guaranteed
global solution from Pereira's finite stochastic search.

## 2. Current boundary and supersession

The accepted runtime capability remains
`promotion-0018-equilibrium-pure-saturation-v1`. The local
`neutral-two-phase-tp-flash-v1` candidate has permanent-lab approval as a local
solver candidate but validation status `NON_ADMISSION`: 12 of 17 May et al.
rows passed the frozen composition contract, four solved rows missed it, and
one row had no package-accepted local state.

The fixed-two-phase candidate therefore receives no promotion. Its direct
`A+PV` formulation becomes HELD Stage III implementation input. The HELD
candidate must delete the unpromoted `two_phase_flash` public operation and
its displaced result surface in the same local checkpoint. No compatibility
alias or parallel fixed-phase route remains.

The neutral-mixture provider SDK tail consumed by that candidate is also
local and unpromoted. It is a required HELD artifact prerequisite, not an
already accepted provider capability. Final cutover requires a separate
ordered provider-tail promotion receipt immediately before any HELD promotion;
provider authority cannot be inferred from a successful consumer build.

HELD improves phase discovery and metastability avoidance. It does not change
the provider model or erase the retained May model/data misses.

## 3. Retained scientific source

The scientific source is the lab Markdown extraction at lab commit
`13ce345b6dcc41d399bb2a4c7b9bedb18f74b45b`, Git blob
`dde7f02d4c93cce86804a8e6b62d37602990ac21`:

```text
docs/papers/md/Equilibrium/
Pereira et al. - 2012 - The HELD algorithm for multicomponent,
multiphase equilibrium calculations with generic equations of.md
```

The slice-owned equations and conventions below are the production
specification. The lab artifact remains evidence and is not a runtime,
build-time, test-time, or sibling-path dependency. Historical lab plans,
algorithm registries, and old HELD implementations are characterization only.

Where the retained paper delegates a detail to earlier work without printing
it, this design names a tested project implementation convention. It does not
invent source attribution or make the omitted paper a hidden prerequisite.

## 4. Provider-to-HELD thermodynamic map

The installed provider owns the dimensionless extensive Helmholtz function

```text
Phi(T,n_1,...,n_C,V) = A(T,n_1,...,n_C,V)/(R T n_ref)
n_ref = 1 mol
```

including ideal mixing and the admitted residual ePC-SAFT contributions. For a
one-mole trial composition, equilibrium sets `n_i = x_i * n_ref`, so the total
volume is also the molar volume. The fixed-pressure phase objective is

```text
g_bar(T,P,n,V) = Phi(T,n,V) + P V/(R T n_ref).
```

Scaling by positive `R T n_ref` does not change a minimizer. Provider
derivatives give

```text
mu_i/(R T) = partial Phi / partial n_i
P_eos       = -R T n_ref * partial Phi / partial V.
```

With `x_C = 1 - sum_(i<C) x_i`, Pereira's independent composition gradient is

```text
(1/(R T)) * partial g / partial x_i = (mu_i - mu_C)/(R T).
```

No fugacity allocation, pressure-density root, or copied EOS expression is
needed.

## 5. Stage I: initialization and stability search

At the normalized feed `z`, equilibrium first minimizes `g_bar(z,V)` from
declared liquid-like and vapor-like starts over molar-volume bounds
`[1e-5,1e-1] m3/mol`. The lowest
locally accepted homogeneous reference supplies `V_z` and the reference
chemical potentials.

The dimensionless tangent-plane distance is

```text
d_bar(x,V;z)
  = g_bar(x,V) - g_bar(z,V_z)
    - sum_(i<C) [(mu_i_ref - mu_C_ref)/(R T)] * (x_i - z_i).
```

Its sign is the sign of Pereira's dimensional tangent-plane distance. Stage I
uses the paper's exponential tunneling form with pole strength `p*=1e-3` and
`10*C = 20` starts for this binary. Starts come from one versioned,
deterministic internal sequence covering feed-near, component-enriched,
liquid-like, vapor-like, and stratified states. The sequence identity and
attempt counts are reported; the public API accepts no seed or initial guess.

For accepted current minimum `(x_star,V_star)` and value `d_star`, the
tunneling objective is the retained source form

```text
f_tunnel(x,V) = [d_bar(x,V) - d_star]
                * exp(p_star / norm(x - x_star)).
```

The implementation must evaluate the pole through an overflow-safe equivalent
and exclude the singular composition neighborhood `norm(x-x_star) = 0`.
Numerical overflow, a singular restart, or provider-domain failure is an
attempt failure, never evidence of stability.

A locally confirmed `d_bar < -1e-8` is sufficient evidence of instability and
initializes Stage II. If all declared attempts complete without such a point,
the result may report one phase with search status
`source_heuristic_complete`. It must not report proven global stability.
Reference or search failure returns `indeterminate`; it is not converted into
a one-phase result.

## 6. Stage II: dual outer/lower loop

For binary composition coordinate `x` and one multiplier `lambda`, the
dimensionless lower problem is

```text
L_bar(x,V,lambda)
  = Phi(T,x,V) + P V/(R T n_ref) + lambda * (z - x).
```

Each evaluated point `(x_m,V_m)` adds the valid outer cut

```text
v <= Phi(T,x_m,V_m) + P V_m/(R T n_ref)
     + lambda * (z - x_m),
v <= g_bar(z,V_z).
```

The initial finite cut set is one declared binary project convention: at each
composition bound, locally polish one liquid-like and one vapor-like volume
start and retain every distinct valid point. At least one valid point at each
opposite composition bound is required. Those endpoint cuts bracket the
admitted feed compositions and make the binary outer envelope bounded;
additional distinct liquid-like or vapor-like endpoint states are retained.
Failure to construct both endpoint cuts returns `indeterminate`; no artificial
multiplier bound is invented.

The admitted outer problem has only two variables, `(v,lambda)`. One internal
binary-envelope solver enumerates feasible cut intersections and selects the
maximum `v`. It is not a reusable general LP package. An independent
test-only exhaustive vertex oracle must agree on value, multiplier, active
cuts, infeasibility, degeneracy, and tie handling.

Each nonconvex lower problem uses Ipopt over `(x,log V)` with molar-volume
bounds `[1e-5,1e-1] m3/mol`. Its deterministic
start classes alternate shifted prior solutions, component-near states, and
stratified states until an accepted lower value no greater than the current
best upper bound is found or the declared search budget is exhausted. Every
accepted local point remains a mathematically valid cut. Search exhaustion is
reported and cannot satisfy the HELD certificate.

Candidate phases satisfy the source upper/lower proximity and multiplier
criteria and are distinct in composition or relative density. Duplicate
candidates are merged deterministically. Two candidates advance to Stage III.
More than two distinct candidates in v1 return `scope_exceeded` rather than
being discarded or forced into a two-phase topology.

The source multiplier criterion is evaluated with the zero-safe project scale
`max(abs(lambda),1)` in the denominator. This preserves a dimensionless
relative test away from zero and prevents a vanishing multiplier from turning
the admission test into an undefined division.

The paper's printed settings are retained where they apply:

```text
composition bounds: 1e-8 to 1 - 1e-8
epsilon_lambda:      0.5
epsilon_b:           1e-2
epsilon_x:           1e-3
epsilon_mu:          1e-6
epsilon_g:           1e-6
```

Pereira's implementation-specific packing-fraction separation is replaced by
the existing relative-density separation `1e-3`. This is a declared project
coordinate convention, not a source claim. The implementation plan must
freeze finite Stage II iteration and per-iteration start budgets before code
is written; exhausting either budget yields `search_exhausted`.

## 7. Stage III: direct total-free-energy refinement

For candidate phases `alpha`, Stage III solves

```text
minimize  sum_alpha [Phi(T,n_alpha,V_alpha)
                     + P V_alpha/(R T n_ref)]
subject to
          sum_alpha n_i_alpha = z_i * n_ref
          n_i_alpha > 0
          V_alpha > 0.
```

Pereira writes the same problem using mass numbers `q_i = M_i n_i`. V1 uses
component amounts directly. This is an exact diagonal change of variables;
the material balances remain linear and no molar-mass owner is introduced
solely for numerical scaling.

For a binary candidate pair, the lever rule supplies initial phase fractions.
Stage III amount bounds are centered on those Stage II candidate amounts with
the source neighborhood `+/- 1e-3 mol`, exactly corresponding to Pereira's
`+/- 1e-3 M_i` mass-number bounds after `q_i=M_i n_i`. If the component
balances are infeasible inside that neighborhood, control returns to Stage II.
Each extensive per-phase volume uses bounds `[1e-5,1e-1] m3`.

The existing fixed-two-phase objective, exact log-volume chain rules, linear
balances, retained Ipopt multipliers, KKT check, distinct-phase check, and
confirmation logic are refactored into this internal stage. Stage III consumes
only the provider value, gradient, and Hessian. Third derivatives are not
required.

If candidate balances are infeasible, the solve fails, or either final HELD
criterion fails, control returns to Stage II within the declared major-
iteration budget:

```text
0 <= UBD - G_star <= epsilon_g
max component/phase relative chemical-potential difference <= epsilon_mu.
```

The chemical-potential comparison uses the symmetric zero-safe denominator
`max(abs(mu_i_alpha),abs(mu_i_beta),1)` in provider `mu/(R T)` units. This
project convention agrees with the printed relative criterion away from zero
and is defined when the chosen Helmholtz reference makes a chemical potential
small.

Final local acceptance also requires material balance, specified-pressure
stationarity, independent KKT stationarity, finite inactive bounds, distinct
positive phase states, and one materially perturbed confirmation solve under
the existing strict package tolerances.

V1's admitted domain contains no source-backed trace-component case. A
candidate at the trace composition bound returns `scope_exceeded`; it is not
silently accepted. Pereira's trace-refinement Step 9 requires a later named
slice and real source-backed evidence before admission.

## 8. Solver coordinates and exact derivatives

Positive phase volume is represented by `y = log V`. This is an exact
one-to-one map of Pereira's volume coordinate, not density closure. The exact
chain rules are

```text
partial f / partial y = V * partial f / partial V
partial^2 f / partial y^2
  = V^2 * partial^2 f / partial V^2 + V * partial f / partial V.
```

Composition reduction and amount-to-phase transformations use explicit linear
Jacobians. Provider Hessians are transformed through those Jacobians and the
log-volume rule. Production uses no finite differences and tapes no Ipopt
iteration history.

CppAD is the provider's derivative mechanism. Using exact Hessians in Ipopt
changes neither Pereira's objectives nor the stationary points; it replaces
the paper implementation's first-derivative local-solver interface with the
already admitted stronger derivative contract.

## 9. Public contract

The only new operation is

```python
tp_flash(
    model: epcsaft.EPCSAFT,
    temperature: pint.Quantity,
    pressure: pint.Quantity,
    overall_mole_fractions: Sequence[float],
) -> TpFlashResult
```

`TpFlashResult` contains an ordered tuple of one or two existing `PhaseState`
objects, phase fractions, total dimensionless free energy, provider
fingerprint, and one compact `HeldDiagnostics`. Diagnostics report stage
outcome, attempt/major-iteration counts, best TPD, final lower and upper
bounds, HELD gap, material/KKT/pressure/chemical-potential checks,
confirmation, deterministic search-profile identity, and
`globality_certificate = "not_guaranteed"`.

One `FlashError` carries the same immutable diagnostics for invalid input,
provider failure, `search_exhausted`, `scope_exceeded`, or `indeterminate`.
There is no public Stage I, Stage II, Stage III, seed, tolerance, solver, or
backend object.

The admitted model fingerprint and rectangular `T,P,z` input domain are
unchanged from the reviewed methane/ethane candidate. Invalid feeds are not
normalized and unsupported states fail before solver dispatch.

## 10. Native and dependency boundary

All HELD code is compiled into the existing
`epcsaft_equilibrium._equilibrium` module and native target. It reuses
`ProviderContext`, `PhaseState`, the provider capsule lifetime owner, Ipopt,
and existing attempt vocabulary.

Internal files follow the scientific stages only when needed to avoid another
giant solver file: one controller/public internal header and focused Stage I,
Stage II, and Stage III implementation owners. They do not expose independent
runtime APIs. No generic NLP hierarchy, route registry, fallback solver,
provider SDK entry, density-root path, continuation system, or new runtime
dependency is allowed.

## 11. Evidence and acceptance

Package evidence is compact and discriminating:

- provider-to-composition/log-volume gradient and Hessian directional checks;
- binary outer-envelope value, active-cut, tie, and infeasibility oracles;
- Stage II candidate addition, deduplication, bound monotonicity, exhaustion,
  and Stage III-return tests;
- direct Stage III material/KKT/certificate and confirmation checks;
- one real May methane/ethane row and the retained collapsed-state rejection;
- one-phase, two-phase, trace-bound, third-candidate, provider-error, and
  unsupported-input outcomes; and
- installed-wheel parity with no lab or sibling source available.

Validation owns the expensive installed-artifact campaign. It uses only the
provider's public pressure-state route to sample the binary Gibbs surface on a
frozen dense composition grid, independently forms the test-only Gibbs value,
and compares HELD's selected phase set and total free energy with the sampled
tangent or chord. It does not reproduce HELD's explicit-volume optimizer or
own a second local-polishing route. This finite scan is an independent sampled
audit, not a continuous globality proof. Validation reruns the complete May
table as honest model/data characterization and retains one observation/model
plot and one representative Gibbs-surface plot with their source rows. Before
implementation, validation must also freeze one traceable methane/ethane
single-phase state for the one-phase terminal path; if no suitable
source-backed state is available, that path remains `indeterminate` and the
candidate cannot be promoted. The four existing composition misses remain
visible and no tolerance is changed.

Promotion requires independent lab review, complete source and installed-
wheel suites, negative-space inspection, validation evidence, a migration
provider-tail promotion receipt, a migration HELD promotion receipt, and
explicit user approval bound to the exact ordered artifacts.

## 12. Exclusions and next boundary

V1 excludes guaranteed deterministic global optimization, arbitrary
multicomponent or multiphase claims, three-phase binary points, trace
refinement, association, ions, electrolyte or reactive equilibrium, bubble or
dew boundary routes, critical continuation, regression, and release work.

After accepted neutral HELD, the next equilibrium capability may extend the
same controller to association-aware provider evaluations. Electrolyte phase
discovery remains a separately formulated capability and must not be labeled
neutral HELD or Pereira HELD2.

## 13. Post-validation binary controller redesign delta

**Deferral:** This design-only delta is retained as historical analysis and is
not the active implementation gate. The user-authorized Perdomo HELD2 pivot
superseded its proposed review/implementation sequence without changing the
frozen neutral runtime, evidence, or authority.

Sections 1--12 and Migration decision D-021 remain the frozen specification
and provenance of the executed v1 candidate. Validation commit
`93ff18541d2fe277a27671e4e6d12b6b009a58ed`, tree
`5aa2bc81941d1e807ba4e579231c4af9b7be15d7`, is accepted as authority-neutral
`NON_ADMISSION` evidence. Its thirteen `third_candidate` exits show that v1
treated loose Stage-II provisional points as confirmed phases before Stage III
could refine them. Its repeated rows 003, 016, and 017 show that a duplicate
below-upper-bound basin stopped the remaining starts without changing the
controller state. This section supersedes only those candidate-lifecycle and
inner-search stopping rules when a future implementation packet is authorized.
It does not reinterpret row 012, which remains a separate Validation-owned
Provider topology question.

### 13.1 Separate cut and candidate-pair ownership

The Stage-II cut pool owns every distinct mathematically valid outer cut
accepted during the current controller solve. Cuts are additive evidence for
the dual envelope: candidate merging, pair rejection, or provisional-candidate
retirement must never remove a cut. A duplicate lower basin may reference an
existing identical cut rather than create a second identity, but its attempt
and rejection evidence remains in the existing trace owner.

The provisional phase-candidate pool is separate controller state. A point may
enter that pool under the existing upper/lower proximity, multiplier,
composition, and relative-density rules without thereby becoming a confirmed
phase. The current binary candidate-pair pool is derived from those provisional
points and owns only canonical candidate-ID pairs plus their refinement state.
Neither pool creates a public type or a second trace/diagnostics owner.
The provisional-candidate hull means the deterministic composition/free-energy
geometry induced by those retained candidates and their feed-bracketing
relations; it is not another outer cut pool or a generic hull solver.

For one current dual/controller state, define progress by the existing
comparison rules: at least one active cut identity, the selected multiplier,
the upper bound, or the provisional-candidate hull changes. Failed-pair memory
is keyed by canonical candidate IDs and retains the exact rejection reason only
for that state. When progress changes the state, feasible pairs are enumerated
again from the updated pool; a rejection is not a permanent blacklist across a
changed envelope or hull.

### 13.2 Deterministic binary pair refinement

Order candidates deterministically by their retained stable identities. For
each distinct pair with compositions `x_a < x_b`, require

```text
x_a <= z <= x_b
beta = (z - x_a)/(x_b - x_a)
phase fractions = (1-beta, beta).
```

Pairs that do not bracket the feed are recorded as nonbracketing and are not
sent to Stage III. Rank feed-bracketing pairs first by increasing lever-rule
provisional free energy

```text
G_pair_0 = (1-beta) * g_bar_a + beta * g_bar_b,
```

then lexicographically by the two stable candidate identities. This ordering
adds no acceptance tolerance or solver setting. Refine each not-yet-rejected
pair in that order through the existing Stage III lever-rule initialization,
candidate neighborhoods, exact Hessian, and unchanged material, pressure,
KKT, inactive-bound, HELD-gap, chemical-potential, distinct-phase, and
confirmation gates.

A Stage-III infeasibility or failed gate rejects only that pair for the current
dual/controller state and records the existing failure evidence. The
controller then considers the next ranked pair and continues the remaining
declared lower starts. More than two provisional points is therefore not by
itself a `scope_exceeded` result.

After refinement evidence exists, phases that collapse under the existing
composition/relative-density distinctness rules are merged in the provisional
pool. A refined point that is equivalent to an already retained provisional
point supersedes the higher-free-energy representative, with stable identity
as the final tie-break. Participating pairs are retired or regenerated as
needed, while every cut remains in the cut pool. If refinement leaves more
than one incompatible certified binary phase pair, or more than two distinct
refined phases that the existing rules cannot merge, the ambiguity is genuine
and remains fail-closed as `scope_exceeded`.

### 13.3 Multistart progress and termination

An accepted lower result below the current upper bound does not stop the
20-start pass when it returns to a duplicate basin and changes no active cut,
selected multiplier, upper bound, or provisional-candidate hull. The attempt
is retained, any already-known cut identity remains valid, and the controller
continues the remaining starts in the unchanged deterministic order.

If one complete 20-start pass changes none of those progress quantities and
produces no newly refinable pair, the controller terminates immediately rather
than repeating the state to the 100-major-iteration cap. It uses the existing
fail-closed `search_exhausted` outcome and the exact failure reason
`no_progress`; no new public outcome or status vocabulary is introduced. If a
start does make progress, the envelope and current candidate-pair pool are
recomputed and the existing major-iteration cap continues to apply.

### 13.4 Frozen contract and future evidence

The redesign changes no thermodynamic equation, derivative transform, Stage
III physics, tolerance, comparison scale, start ordering, start count, major-
iteration cap, provider dependency, binary methane/ethane domain, public
`tp_flash`/`TpFlashResult`/`HeldDiagnostics`/`FlashError` owner, or
`globality_certificate="not_guaranteed"`. The fixed resource profile remains
20 Stage-II lower starts per iteration and 100 major iterations.

A future authorized implementation must first add focused RED evidence for:

- more than two provisional candidates with one valid feed-bracketing pair
  that refines successfully;
- degenerate provisional candidates that merge or retire only after refinement;
- genuine post-refinement ambiguity that remains `scope_exceeded`;
- a failed or nonbracketing pair that leads to alternative starts or pairs
  instead of an identical loop;
- twenty duplicate/nonprogressing starts that terminate promptly with
  `search_exhausted` and exact reason `no_progress`; and
- unchanged derivatives, three-axis diagnostics mapping, globality, retired-
  owner negative space, and exact artifact binding.

No generic N-phase API, retired fixed-two-phase route, Provider implementation,
new resource, or tolerance change is part of this delta.
