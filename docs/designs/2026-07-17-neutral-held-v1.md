# Neutral HELD v1 Design

Status: complete local candidate at implementation commit
`549162a3a9cfd6f02894f8189c624ba1aa2139fb`; awaiting permanent-lab
complete-candidate review and Validation's installed-artifact campaign

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
