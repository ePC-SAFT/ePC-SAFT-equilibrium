# Homogeneous Chemical-Equilibrium Basin Search

**Status:** approved pre-result implementation contract

**Date:** 2026-07-30

**Issue:** [#96](https://github.com/ePC-SAFT/ePC-SAFT-equilibrium/issues/96)

**Result freeze:** This record was committed before executing the new search
against the Hilliard/Böttinger MEA sentinel. Sentinel behavior cannot change
the formulation, start order, budgets, tolerances, or acceptance rules below.

## Scope and authority

This design extends the existing fixed-temperature, fixed-pressure,
single-homogeneous-phase chemical-equilibrium owner with one deterministic
finite basin search. It supersedes only these portions of the earlier
reacting-phase records:

- the one-direct-start rule;
- the Provider `log(V)` description, which the implementation had already
  replaced with inverse `log(packing_fraction)` coordinates; and
- the binary reduced-curvature classification.

The compiler, support classifier, amount chart, Provider boundary, physical
objective, exact derivative ownership, source-reference transformation, and
strict local/non-global claim remain unchanged.

The search returns the lowest observed certified strict local minimum. It does
not establish phase stability, kinetic stability, global thermodynamic
stability, predictive validity, infeasibility, or regression readiness.
Finite-search exhaustion is evidence about the declared attempt sequence only.

MEA supplies one immutable downstream stress subject. No production branch,
start, tolerance, or status depends on a chemistry name, source row, species
name, or application.

## Formulation identity

For physical amounts `n` and volume `V`, Equilibrium minimizes

```text
Phi(n,V) = provider_mechanical_Helmholtz(n,V)
         + P*V/(R*T)
         + dot(g_ref,n)
```

subject to the compiled independent conservation equations. The ionic amount
chart enforces electroneutrality identically. The Provider owns the mechanical
Helmholtz value, state domain, parameter identity, packing geometry, pressure,
and exact derivative tensors. Equilibrium owns the pressure-volume and
reference terms, coordinate pullbacks, constraints, optimization, and
certification.

`Phi` is dimensionless. The compiler's minimum-norm `g_ref` fixes its reported
gauge. A balance-row gauge shift adds the same constant to every feasible state
and cannot change basin ordering. Search receipts record the objective in this
exact gauge.

Every direct, multistart, recovery, and certification evaluation calls the
same compiled objective and constraints. No attempt may substitute an ideal
objective, omit Provider terms, or tape an iterative solve.

## Coordinates, bounds, and reconstruction

Neutral systems use one log amount per species. Ionic systems use log charge
equivalents, cation and anion softmax shares, and neutral log amounts. The
existing analytic amount values, Jacobian, and component Hessians remain the
only coordinate pullback.

Ideal manufactured systems use `log(V)`. Installed Provider systems use
inverse `log(packing_fraction)` geometry supplied by the Provider, including
the exact amount-dependent volume gradient and Hessian. Search code must not
treat that coordinate as `log(V)`.

The max-min initialization LP supplies finite amount upper bounds from the
positive bounding balance. Neutral amount-coordinate lower bounds remain
`log(0.1 * trace_floor)`. Ionic shares retain the Provider-safe
`40 + log(10)` coordinate envelope and the existing trace-aware formulas.
Provider packing bounds are contracted inward with `nextafter`. Every attempt
uses these same bounds.

Each start and terminal is reconstructed in physical `(n,V)` coordinates and
checked against:

- finite positive amounts and volume;
- compiled balances and exact chart charge;
- the declared trace floor;
- artificial coordinate-bound inactivity;
- Provider packing and source-domain limits; and
- the exact physical objective.

Trace-floor contact or artificial-bound contact is
`boundary_unadjudicated`; it cannot become a certified local minimum.

## Local Ipopt attempt

`ReactionTnlp` remains the sole local optimizer. Every primary or recovery
attempt uses:

```text
linear_solver             = mumps
tol                       = 1e-10
acceptable_tol            = 1e-9
acceptable_iter           = 0
constr_viol_tol           = 1e-10
jacobian_approximation    = exact
hessian_approximation     = exact
nlp_scaling_method        = none
bound_relax_factor        = 0
honor_original_bounds     = yes
check_derivatives_for_naninf = yes
```

Only `Solve_Succeeded` satisfies the solver gate. The per-attempt iteration
budget remains 500 unless a private manufactured seam deliberately supplies a
smaller value to test exhaustion. The immutable artifact records the Ipopt
version and linear-solver identity. The implementation target is the current
package environment, Ipopt 3.11.9 with sequential MUMPS 5.6.

The existing strictly interior KKT polish remains part of each attempt.
Balances, reaction affinities, and pressure use the same exact objective and
derivatives as the Ipopt callback.

## Deterministic primary starts

The search constructs at most 25 primary starts in this order:

1. ordinal zero is the deterministic max-min LP state;
2. for each compiled independent reaction row in row order, capped after six
   rows, use extent fractions `0.5` then `0.875`;
3. at each fraction, use the positive direction before the negative direction.

For reaction row `nu`, sign `s`, and fraction `f`, the physical seed is

```text
n_seed = n_maxmin + s*f*lambda_max*nu
```

where `lambda_max` is the largest finite step that keeps every amount strictly
above `trace_floor` and strictly below its max-min LP upper bound. The
implementation rejects a zero, nonfinite, duplicate, or unavailable direction.
It reconstructs the seed through the amount chart and applies the existing
balance retraction when binary64 reconstruction exceeds the balance tolerance.
A failed retraction remains an attempt record and launches no Ipopt solve.

For each accepted physical amount seed, the Provider route recomputes molar
volume bounds and performs a fixed 128-interval scan in log volume. Every
detected sign-changing pressure root is refined by at most 100 bisections with
the stopping test `abs(residual) <= 1e-10 * P`. Distinct roots are ordered by
distance from the log-volume midpoint and interleaved breadth-first across
amount seeds. A failed Provider scan, inverse-packing reconstruction, trace
check, or domain check remains a rejected start record; detected roots omitted
by the fixed primary ceiling are counted as budget-truncated.

The ordinary primary search does not hide a continuation start. If it finds no
certified target-model basin and the caller supplied a typed Provider-model
continuation, the independently receipted recovery route is governed by
[`2026-08-02-provider-model-continuation.md`](2026-08-02-provider-model-continuation.md).
The fixed ordinary basin search and its receipt remain unchanged and mandatory.

## Negative-curvature recovery

The direct/local baseline from PRs #85 through #87 remains unchanged within
each primary attempt:

- use the independently reconstructed exact Lagrangian Hessian;
- try curvature signs `+1` then `-1`;
- test four dyadic displacement lengths per sign;
- retract each candidate onto the exact balance manifold with at most eight
  Newton/KKT corrections and eight dyadic line-search fractions;
- hold the volume coordinate fixed during retraction;
- require exact-objective decrease beyond
  `64*epsilon*max(1,abs(Phi0),abs(Phi_seed))`; and
- launch at most one nonrecursive Ipopt recovery solve per sign.

Every launched recovery is a child attempt with its parent primary ordinal,
sign, seed ordinal, solver status, terminal classification, and objective.
Failed displacement, retraction, domain, and objective screens remain in the
parent attempt's seed accounting.

## Strict local-minimum certificate

The solver, numerical, physical, rank, and second-order gates remain
independent.

The original-coordinate residual limits are:

| Check | Limit |
| --- | ---: |
| balance infinity norm | `1e-9` |
| charge infinity norm | `1e-9` |
| pressure relative residual | `1e-8` |
| reaction-affinity infinity norm | `1e-7` |
| KKT stationarity infinity norm | `1e-7` |
| complementarity infinity norm | `1e-7` |

The certificate independently reconstructs equality multipliers from the
physical state. It forms the same exact KKT matrix used by conditioned
sensitivities, applies four deterministic infinity-norm row/column
equilibration passes, and requires full rank with
`condition_number_inf <= 1e6`. Rank deficiency or worse conditioning is
`second_order_inconclusive`, even when Ipopt converges.

For second-order classification, Equilibrium forms the null space of the
actual transformed equality Jacobian after deterministic row normalization.
The reduced exact Lagrangian Hessian uses diagonal equilibration. Let `d` be a
unit transformed-coordinate direction and
`H_scale = max(1, norm_inf(H_lagrangian))`.

- every normalized Cholesky pivot greater than `1e-10` is
  `certified_local_minimum`;
- a reconstructed direction with
  `dot(d,H_lagrangian*d) < -1e-10*H_scale` is `saddle_observed`; and
- a zero, near-zero, nonfinite, rank-ambiguous, or unclassified pivot is
  `second_order_inconclusive`.

Sensitivity columns and composed observables are emitted only for the selected
certified local minimum. Search failures, saddles, boundary terminals, and
inconclusive terminals suppress them.

## Basin identity and selection

The search retains every materially distinct certified basin. Two certified
terminals are duplicates only when all three checks pass:

```text
max_i abs(n_i - m_i) / max(total(n), total(m), trace_floor) <= 1e-8
abs(log(V_n / V_m))                                      <= 1e-8
abs(Phi_n - Phi_m) <= 256*epsilon*max(1,abs(Phi_n),abs(Phi_m))
```

The earlier ordinal is the canonical representative of a duplicate basin.
Duplicate attempts remain in the receipt and point to that basin ordinal.

The selected basin has the lowest physical objective among certified basin
representatives. Objective ties within the objective tolerance use the lower
basin ordinal. The receipt names the selected value
`lowest_observed_certified_local_value`.

## Budgets and receipt

The fixed primary budget is 25. The actual generated count can be lower when
the reaction dimension is smaller or a direction cannot produce a distinct
strict-interior seed. For an installed Provider, a fixed 128-interval
log-volume scan brackets every detected sign-changing pressure root in the
admitted volume interval. Distinct roots are ordered by distance from the
log-volume midpoint. Starts are then interleaved breadth-first by root rank, so
every constructible composition receives its first root before any composition
receives a second, until the same fixed ceiling is reached; no detected root is
silently chosen as the only density branch. Generated roots omitted only by the
fixed ceiling are reported as budget-truncated; construction-, duplicate-,
infeasible-, and Provider-domain-rejected starts are likewise accounted
separately from evaluated starts. Each evaluated primary attempt can launch at most two recovery
solves. No random seed, adaptive extra start, or result-dependent budget
extension exists.

The receipt reports the clipped nested primary prefixes `1`, `5`, `13`, and
`25`. For each prefix it records the attempted ordinals, certified basin
ordinals, selected basin, and whether selection changed from the preceding
prefix. Prefixes are projections of one fixed full run; they do not trigger
new solves.

Every attempt record includes:

- primary ordinal, kind, parent ordinal, and deterministic start identity;
- start construction, retraction, continuation, and Provider-domain status;
- solver and callback status;
- finite terminal amounts, volume, and objective when available;
- balance, charge, pressure, affinity, KKT, complementarity, trace, rank,
  conditioning, and curvature evidence;
- returned-state central directional checks of the objective gradient,
  constraint Jacobian, and Lagrangian Hessian, including step, scaled error,
  and fixed tolerance;
- terminal classification and duplicate basin reference; and
- recovery seed and launched-solve accounting.

The top-level search status uses only:

- `certified_local_minimum`;
- `saddle_observed`;
- `second_order_inconclusive`;
- `boundary_unadjudicated`;
- `domain_rejected`;
- `no_feasible_start_found`;
- `search_exhausted_no_certified_candidate`; and
- `infeasible_certified`.

`infeasible_certified` is reserved for an independently checked feasibility
certificate. Max-min failure, Ipopt's infeasibility status, all observed
saddles, and budget exhaustion cannot emit it. Until such a certificate owner
exists, no production path emits that status.

When no candidate is certified, the top-level status uses the most specific
observed fail-closed class only when every launched terminal shares it.
Mixed failure classes use `search_exhausted_no_certified_candidate`.

The native and Python result schemas carry the full search receipt on both
success and `ChemicalEquilibriumError`. The existing Regression evaluator v1
ABI remains unchanged: a row is `OK` only when the selected state is a
certified local minimum and requested derivatives are available. Its artifact
and contract fingerprints bind that decision to this search implementation.

## Evidence

Package evidence must cover these claim groups:

1. deterministic start construction, exact-manifold reconstruction, repeat
   ordering, nested prefixes, and budget exhaustion;
2. two certified basins with different objectives, duplicate recognition, and
   lowest-certified selection;
3. saddle recovery success and failure, all-saddle, inconclusive curvature,
   rank/conditioning failure, trace/boundary contact, Provider-domain failure,
   and no feasible start;
4. the existing independent directional checks of objective gradients,
   constraint Jacobians, and Lagrangian Hessians; and
5. unchanged ideal, water self-ionization, observation, and non-MEA public
   behavior.

The frozen Hilliard/Böttinger sentinel runs only after this record is committed.
It either returns a reproducible certified local minimum or remains a blocker
with the complete finite-search receipt.
[Validation issue #18](https://github.com/ePC-SAFT/ePC-SAFT-validation/issues/18)
owns the preregistered MEA matrix over immutable installed artifacts.

### Post-design ownership clarification

The user subsequently confirmed that MEA testing is a downstream application
of this generic capability. The package leaf therefore freezes the immutable
artifact, public receipt, and exact sentinel identity, while
[Validation issue #18](https://github.com/ePC-SAFT/ePC-SAFT-validation/issues/18)
owns both the first installed-artifact sentinel execution and the broader
matrix. This clarification supersedes only the execution owner implied by the
preceding paragraph; it does not change any formulation, start, solver,
tolerance, certification, or failure-accounting rule.

## Out of scope

Reactive phase equilibrium, phase discovery, tangent-plane certification,
global optimization, parameter fitting, application chemistry ownership,
random starts, caller solver controls, alternate derivative backends,
predictive admission, promotion, and release remain outside this design.
