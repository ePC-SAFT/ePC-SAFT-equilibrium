# HELD2.0 for Explicit-Species Electrolyte Phase Equilibrium

Status: publication-facing companion; not a second implementation authority

Date: 2026-07-29

Normative implementation owner:
[`2026-07-24-held2-paper-algorithm.md`](2026-07-24-held2-paper-algorithm.md)

Necessary-condition map:
[`2026-07-28-held2-necessary-condition-map.md`](2026-07-28-held2-necessary-condition-map.md)

Validated working state:
[`2026-07-28-held2-validated-working-baseline.md`](2026-07-28-held2-validated-working-baseline.md)

## Purpose

This document presents the package's revised HELD2.0 method in a form suitable
for a manuscript Methods section or algorithm appendix. It separates the
scientific formulation, algorithm, necessary conditions, and deliberate
revisions from repository layout, historical cutover, and case-specific
performance records.

The longer normative specification remains the sole owner of executable
transitions, named tolerances, resource profiles, diagnostics, and file
responsibilities. If the two documents differ, the normative specification
governs the software and this companion must be corrected.

## Problem and claim

At specified temperature \(T\), pressure \(P_0\), and overall explicit-species
composition \(\mathbf z\), the method searches for one or more homogeneous
electrolyte phases that satisfy:

- normalization and per-phase electroneutrality;
- overall species material balance;
- a common pressure;
- equality of independent electrochemical potentials;
- Provider domain constraints; and
- local first- and second-order numerical acceptance conditions.

The equation of state and all nonlinear thermodynamic derivatives are supplied
by the installed Provider public SDK. Equilibrium owns coordinate
transformations, pressure-root enumeration, finite global exploration, LP and
NLP formulations, and independent evidence.

The phase search is finite. Every successful result therefore retains

```text
globality_certificate = "not_guaranteed"
```

Absence of a discovered negative tangent-plane-distance state is an
operational finite-search result, not a proof of global stability.

## Modified electrolyte coordinates

Let \(C\) denote the number of explicit species and let \(\mathbf q\) contain
their charges. HELD2.0 eliminates one reference species and one charged species
to obtain \(C-2\) independent modified composition coordinates
\(\bar{\mathbf x}^{(EC)}\). The lift

\[
\mathbf x
=
\mathcal L\!\left(\bar{\mathbf x}^{(EC)}\right)
\]

is linear on a complete polytope and satisfies

\[
\mathbf 1^T\mathbf x=1,
\qquad
\mathbf q^T\mathbf x=0.
\]

The coordinate map preserves Provider component order through explicit
permutations. It is not a counterion-pair representation and does not alter
the physical species basis.

Finite lower bounds used during global exploration regularize a bounded
search. They are not physical minimum mole fractions. Ordinary Steps 1--9 use
linear modified coordinates. Step 10 alone may refine a strictly positive
charged trace coordinate in \(\log_{10}x_i\) while retaining every physical
constraint.

For molar volume \(V\), define the dimensionless one-reference-mole Gibbs
quantity

\[
g^{el}(T,P_0,V,\bar{\mathbf x}^{(EC)})
=
\frac{\bar A^{el}(T,V,\bar{\mathbf x}^{(EC)})}{RT}
+
\frac{P_0V}{RT}.
\]

Provider composition-dependent volume bounds are mapped into a normalized
\(\log V\) coordinate whenever a bounded deterministic search is required.

## Reference state and tangent-plane distance

At the feed composition, the method enumerates the complete declared
log-volume pressure envelope, refines every detected pressure root, classifies
mechanical stability, and selects the mechanically stable root with the lowest
\(g^{el}\). Ambiguous root ordering or incomplete required evidence fails
closed.

Let the selected homogeneous reference be
\((\bar{\mathbf x}_0^{(EC)},V_0)\), with value \(g_0^{el}\) and reduced
independent modified potentials

\[
\bar{\boldsymbol\mu}_0
=
\nabla_{\bar{\mathbf x}^{(EC)}}g^{el}
(T,P_0,V_0,\bar{\mathbf x}_0^{(EC)}).
\]

The tangent-plane distance at any feasible trial state is

\[
\operatorname{TPD}(V,\bar{\mathbf x}^{(EC)})
=
g^{el}(T,P_0,V,\bar{\mathbf x}^{(EC)})
-g_0^{el}
-\bar{\boldsymbol\mu}_0^T
\left(
\bar{\mathbf x}^{(EC)}-\bar{\mathbf x}_0^{(EC)}
\right).
\]

A feasible, independently reevaluated state with TPD below the named strict
negative threshold proves that the homogeneous reference is unstable.
Pressure stationarity is not required for this existence witness: the
minimization is joint in composition and volume, and any feasible off-pressure
state below the tangent plane suffices. Pressure refinement is performed only
when a promising witness is needed as a downstream phase seed.

## Algorithm

### Stage I: reference, instability, and initialization

**Step 1: construct the physical search domain.**

Build the exact modified-composition polytope, the Provider-order lift, and
composition-dependent physical volume bounds. Validate \(T\), \(P_0\), feed
normalization, feed charge, component identity, Provider fingerprint, and
every immutable Provider domain restriction.

**Step 2: search for a negative TPD state.**

First enumerate and select the homogeneous feed pressure root. Then search
jointly over the \(C-2\) modified composition coordinates and normalized
\(\log V\). Count actual Provider callbacks, independently reevaluate any
strict-negative candidate, and terminate early only after that
recertification. If the declared search completes without a negative state,
return the operational one-phase result with globality not guaranteed.

**Step 3: initialize the persistent state set.**

Construct the paper's Appendix-C one-dimensional endpoint states, retain the
feed/reference state, and add a Step-2 witness only when it contributes
multivariate information not already represented. Every state receives an
immutable insertion identity. Numerical-copy identity is based on modified
composition and relative molar volume.

### Stage II: finite supporting-plane refinement

Let \(\mathcal M\) be the append-only set of retained physical states.

**Step 4: solve the restricted upper LP.**

Solve

\[
\begin{aligned}
\max_{v,\bar{\boldsymbol\lambda}}\quad &v\\
\text{s.t.}\quad
v\leq&
g^{el}(V^m,\bar{\mathbf x}^{(EC),m})
+
\bar{\boldsymbol\lambda}^T
\left(
\bar{\mathbf x}_0^{(EC)}
-\bar{\mathbf x}^{(EC),m}
\right),
\quad m\in\mathcal M,\\
v\leq&g_0^{el}.
\end{aligned}
\]

Independently audit primal feasibility, dual feasibility, and
complementarity. A solver status without that evidence is not an LP
certificate.

**Step 5: solve the nonconvex lower problem.**

For the current upper-LP multipliers, solve

\[
\min_{V,\bar{\mathbf x}^{(EC)}}
g^{el}(V,\bar{\mathbf x}^{(EC)})
+
\bar{\boldsymbol\lambda}^T
\left(
\bar{\mathbf x}_0^{(EC)}
-\bar{\mathbf x}^{(EC)}
\right)
\]

over the Step-1 domain. Use one persistent deterministic multistart sequence
and exact Provider derivatives. Stop when the best independently acceptable
local value is no greater than the current upper bound. Add the terminal state
to \(\mathcal M\) only when it is not a numerical copy; an equivalent state is
recorded without inventing a new cut.

The completed implementation must certify the terminal in original physical
coordinates, including primal feasibility, KKT stationarity, bound-multiplier
signs, complementarity, and dual reconstruction.

**Step 6: select eligible phase candidates.**

Reevaluate all of \(\mathcal M\) against the current gap, fixed-volume
gradient, and Perdomo distinctness conditions. Candidate distinctness uses two
independent physical axes: modified composition and actual Provider packing
fraction. Packing fraction is not replaced by molar-volume distance.

**Step 7: advance the major iteration.**

If fewer than two eligible distinct phase candidates exist, advance the major
counter, consume the next unused Step-5 start ordinal, and return to Step 4.
Resource exhaustion is indeterminate.

### Stage III: simultaneous phase optimization and final evidence

**Step 8: solve total Gibbs minimization.**

For the selected candidate set, solve

\[
\min_{\boldsymbol\phi,\{\bar{\mathbf x}^m,V^m\}}
\sum_m
\phi^m
g^{el}(V^m,\bar{\mathbf x}^m)
\]

subject to exact modified material balances,
\(\sum_m\phi^m=1\), nonnegative phase amounts, Step-1 domains, and the
candidate composition neighborhoods prescribed by Problem (67).

An exact perspective-form LP first establishes a feasible master-to-NLP start.
Only validated Farkas evidence may classify that LP as infeasible. Ipopt uses
exact derivatives and its terminal state is independently audited in physical
coordinates.

Numerical duplicate phases are merged using physical composition and relative
molar volume. A distinct phase may be retired only from positive
lower-bound-multiplier KKT evidence; retire at most one and solve the complete
reduced problem again. A failed warm solve permits one identical cold retry,
not a changed problem.

**Step 9: apply paper convergence tests.**

Require both the same-major free-energy gap and equality of independent
modified potentials, while preserving separate balance, charge, pressure,
domain, KKT, and phase-identity gates. Passing a paper ratio cannot erase a
failed physical certificate. A failed convergence gate returns to Stage II
with certified phase feedback.

**Step 10: refine charged traces.**

When an otherwise accepted solution contains eligible strictly positive
charged trace fractions, refine only those coordinates in bounded
\(\log_{10}x_i\). Reconstruct the ordinary physical phases and repeat the full
balance, charge, pressure, potential, domain, KKT, and identity audit. Small
accurate values are preserved rather than clipped or retired.

## Controller pseudocode

```text
compile Step-1 coordinates and Provider domain
reference = enumerate_and_select_feed_pressure_root()
witness = joint_tpd_search(reference)

if witness is not independently strict-negative:
    if declared search completed:
        return one_phase(no_negative_witness_detected,
                         globality_not_guaranteed)
    return indeterminate

M = appendix_c_states(feed, witness)
initialize persistent major and start counters

while resource budget remains:
    upper = independently_certified_upper_lp(M)                 # Step 4
    lower = independently_certified_lower_multistart(upper)     # Step 5
    append_if_new(M, lower.state)
    candidates = evaluate_all_members(M, upper)                 # Step 6

    if fewer than two distinct eligible candidates:
        advance_without_replaying_starts()                      # Step 7
        continue

    phases = exact_master_then_total_gibbs_nlp(candidates)      # Step 8
    phases = merge_or_retire_one_then_resolve(phases)
    convergence = paper_and_physical_certificates(upper, phases) # Step 9

    if convergence requests more candidates:
        feed_certified_phase_feedback_to_M()
        continue

    refined = refine_positive_charged_traces(phases)            # Step 10
    if refined passes every final certificate:
        return equilibrium(globality_not_guaranteed)
    return indeterminate

return indeterminate(resource_exhausted)
```

## Independent identities and evidence

Three proximity questions remain separate:

1. persistent-\(\mathcal M\) numerical identity uses modified composition and
   relative molar volume;
2. Step-6 candidate distinctness uses modified composition and Provider
   packing fraction; and
3. Step-8 duplicate-phase merging uses physical composition and relative molar
   volume.

Each scientific predicate owns a named tolerance and retains its raw residual,
scale, threshold, and failure meaning. Solver termination, numerical
acceptance, physical acceptance, finite-search completion, root completeness,
predictive status, and globality are independent evidence axes.

## Deliberate revisions relative to a literal paper transcription

The package adds only the machinery needed to make the paper-style finite
algorithm auditable for the installed Provider:

- complete feed pressure-root enumeration and lowest-objective stable-reference
  selection;
- a joint composition--log-volume Step-2 search rather than a full nested
  pressure envelope at every composition;
- independent negative-witness reevaluation;
- Provider-evaluation resource accounting and deterministic persistent starts;
- immutable candidate identities and nontransitive pairwise-distinctness rules;
- exact master-to-NLP initialization;
- independently audited LP, KKT, balance, pressure, and domain evidence;
- one-at-a-time KKT retirement followed by a complete reduced re-solve;
- evidence invalidation whenever the mathematical problem changes; and
- logarithmic coordinates only for bounded trace refinement.

These additions do not alter the EOS, claim global optimization, or introduce
a sequential chemistry calculation.

## Evidence required before reuse in GREPE

The working HELD2 behavior is frozen by the native Steps 1--10 workflow,
Perdomo numerical regressions, and native/Python diagnostic parity. Before
shared numerical machinery is extracted for GREPE,
two documented evidence gaps must be closed in isolation:

1. the Step-5 original-physical-coordinate KKT and dual-reconstruction audit;
2. validated Farkas evidence for Step-8 LP infeasibility.

Both changes must preserve the frozen numerical results and fail closed when
the additional evidence cannot be established.
