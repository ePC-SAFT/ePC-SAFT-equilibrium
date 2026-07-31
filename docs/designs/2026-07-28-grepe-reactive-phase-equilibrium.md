# GREPE Reactive Phase Equilibrium

Status: normative future implementation contract; no runtime capability or
authority effect

Date: 2026-07-29

Reviewed source: `GREPE_algorithm_notation_refined.md`, SHA-256
`e8c5459890fa086056984373d945d70ac34d4dedb2fc5087b55b39a5ccf2fc4d`

## Purpose and precedence

This document is the Equilibrium-package revision of the reviewed GREPE
manuscript. It preserves GREPE's conserved-component formulation and final
simultaneous reactive multiphase minimization while replacing its generic
phase-discovery policy with the numerical and evidence contracts validated by
the package's current HELD2 workflow.

The nonreactive implementation remains owned by the
[HELD2 algorithm](2026-07-24-held2-paper-algorithm.md),
[publication companion](2026-07-29-held2-publication-algorithm.md),
[condition map](2026-07-28-held2-necessary-condition-map.md), and
[validated baseline](2026-07-28-held2-validated-working-baseline.md).
Reaction-free strong-electrolyte inputs continue to dispatch to the single
`run_held2_algorithm` controller with unchanged Steps 1--10, tolerances,
budgets, Provider contract, and numerical gates.
The condition map remains the sole owner of HELD2's outstanding Step-5 KKT
and Step-8 Farkas evidence gaps.

The implemented reaction compiler, support classification, accessible-face
recompilation, thermodynamic phase block, amount chart, and local Ipopt solve
remain owned by the
[GREPE homogeneous layer](2026-07-27-grepe-homogeneous-chemical-layer.md).
The coupled design composes those owners; it does not copy them.

The public-facade decision is recorded in
[ADR 0001](../adr/0001-one-equilibrate-facade-over-specialized-owners.md).
The first physical evidence subjects are ranked in the
[reactive-LLE benchmark review](2026-07-29-reactive-lle-benchmark-review.md).

## Claim boundary

The first executable slice is fixed-\(T,P\) phase discovery over at most two
liquid phases using one installed Provider model for both phases. The planned
controller returns one liquid only when the declared search completes without
a negative pricing witness and every one-phase evidence check passes. It
returns two liquids when the complete coupled evidence accepts a split and
otherwise fails closed. The caller does not label phases as aqueous or organic
and does not prescribe a phase count.

Every compiler-accessible species is admitted in both liquids. Arbitrarily
small positive amounts remain physical; species are never retired by
magnitude. A phase-specific exact structural zero that cannot be certified
through the installed Provider contract fails closed.

The first physical benchmark is neutral reactive LLE. Electrolyte reactive LLE
is a later capability slice requiring a source-complete charged-species
subject, phase-specific electroneutrality, electrolyte reference-state
evidence, and unchanged neutral results. Neutral evidence cannot imply the
electrolyte claim.

Solids, vapor, surfaces, double layers, kinetics, regression, chemistry
databases, arbitrary phase families, and more than two liquids are out of the
first implementation scope.

Finite deterministic search is not a proof of global phase stability. Until a
separate implementation supplies valid global pricing lower bounds over every
declared phase-family domain, every result retains

```text
globality_certificate = "not_guaranteed"
```

and may claim only local equilibrium plus the declared finite-search result.

## One facade, specialized numerical owners

The future public concept is one typed `equilibrate` operation with distinct
problem variants:

1. a constrained homogeneous reacting problem delegates to the existing
   chemical-equilibrium owner and makes no phase-stability claim;
2. a qualifying nonreactive phase-discovery problem delegates to the existing
   HELD or HELD2 owner with unchanged qualification and rejection behavior;
   and
3. a reactive phase-discovery problem runs this simultaneous GREPE controller,
   even when its accepted result contains one phase.

`tp_flash` and `chemical_equilibrium` remain frozen public operations
during development. The private facade delegates to the same native owners and
must prove parity; it does not copy either implementation. Public `equilibrate`
exposure and any later removal of the existing operations are separate
promotion decisions after physical reactive-LLE evidence.
Using the facade does not admit a nonreactive system that the corresponding
existing route rejects.

A homogeneous chemical result is a Stage-1 seed for reactive phase discovery,
not proof that the reactive system is stable as one phase. GREPE never
speciates once and then treats a phase-only flash as the coupled solution.

## Shared thermodynamic ownership

For phase family \(f\), let the installed Provider supply the complete
mechanical reduced Helmholtz block

\[
\Phi_f(T,V,\mathbf n)
\]

and its required nonlinear derivative tensors in Provider component order.
Let Equilibrium's reaction compiler supply the dimensionless reference vector
\(\mathbf g_{\mathrm{ref}}\), satisfying the compiled reaction constants in
the same Provider Helmholtz basis.

For an extensive phase state, define

\[
J_f(T,P_0,V,\mathbf n)
=
\Phi_f(T,V,\mathbf n)
+\frac{P_0V}{RT}
+\mathbf g_{\mathrm{ref},f}^{T}\mathbf n.
\]

Provider owns \(\Phi_f\), component identity, validity domains, packing
fractions, volume bounds, and nonlinear tensors. Equilibrium owns the linear
reference term, conserved-component and charge coordinates, LPs, NLPs,
phase discovery, and certificates. Equilibrium does not copy Provider
equations or reconstruct Provider parameters.

## Compiled chemistry and conserved totals

The typed chemistry input supplies ordered species feed amounts, source
reactions and equilibrium constants, and a physical conserved-component
matrix. It does not supply optimization bases or reduced coordinates. The
existing compiler audits the supplied conservation semantics, removes
redundant rows and reactions, and derives every internal chart.

On the accessible species face, let \(\mathbf B\) be the compiled full-row-rank
conserved-component matrix, \(\mathbf z\) the charge vector,
\(\boldsymbol\nu\) the independent reaction basis, and
\(\mathbf b=\mathbf B\mathbf n^F\) the feed totals.

The compiler must verify

\[
\operatorname{null}
\begin{bmatrix}
\mathbf B\\
\mathbf z^T
\end{bmatrix}
=
\operatorname{row}(\boldsymbol\nu)
\]

The feed is globally electroneutral. Only independently certified structural
zeros may be removed; unresolved species remain. Compilation, reaction-cycle
checks, standard-state conversion, and accessible-face reconstruction use the
existing homogeneous owner.

## Normalized phase states

Let \(\mathbf s^T\) be a positive conserved scaling row. Reactive GREPE uses
the molar-mass row because total mass is conserved by every admitted reaction.
The strict nonreactive HELD2 reduction instead retains HELD2's existing
one-reference-mole normalization. For a nonzero phase define its scale
\(\theta\), normalized species amounts \(\widehat{\mathbf n}\), and scaled
volume \(v\):

\[
\theta=\mathbf s_f^T\mathbf n,\qquad
\widehat{\mathbf n}=\frac{\mathbf n}{\theta},\qquad
v=\frac{V}{\theta}.
\]

The normalized physical phase domain is

\[
\mathcal Y_f=
\left\{
(\widehat{\mathbf n},v):
\widehat{\mathbf n}\ge0,\;
\mathbf s_f^T\widehat{\mathbf n}=1,\;
\mathbf z_f^T\widehat{\mathbf n}=0,\;
(T,v,\widehat{\mathbf n})\text{ is Provider-valid}
\right\}.
\]

Define the normalized phase value and conserved column

\[
g_f(y)=
\Phi_f(T,v,\widehat{\mathbf n})
+\frac{P_0v}{RT}
+\mathbf g_{\mathrm{ref},f}^{T}\widehat{\mathbf n},
\qquad
\mathbf c_f(y)=\mathbf B_f\widehat{\mathbf n}.
\]

The normalization is a mathematical representation for the phase master and
pricing problem. Changing between two valid positive conserved scaling rows
rescales columns and phase weights but does not alter the underlying extensive
equilibrium. It does not impose a positive physical minimum phase amount.

One formulation-owned nondimensionalization contract is derived from the
physical problem and reused by every GREPE stage. It uses a conserved feed
scale for amounts, normalized Provider-bounded log volume, and thermal-energy
scales for objectives and multipliers. Master constraints, pricing residuals,
the final NLP, derivative tests, and diagnostics all report in that same
contract. Callers cannot provide stage-specific solver scales.

## Semi-infinite master interpretation

The reactive multiphase problem can be written as

\[
\min_{\theta_j,y_j}
\sum_j\theta_j g_{f_j}(y_j)
\]

subject to

\[
\sum_j\theta_j\mathbf c_{f_j}(y_j)=\mathbf b,\qquad
\theta_j\ge0,\qquad
y_j\in\mathcal Y_{f_j}.
\]

For conserved-component multipliers \(\boldsymbol\lambda\), define the phase
reduced cost

\[
r_f(y;\boldsymbol\lambda)
=
g_f(y)-\boldsymbol\lambda^T\mathbf c_f(y)
\]

and pricing value

\[
\rho_f(\boldsymbol\lambda)
=
\min_{y\in\mathcal Y_f}r_f(y;\boldsymbol\lambda).
\]

A feasible state with

\[
r_f(y;\boldsymbol\lambda)<-\epsilon_{\mathrm{price}}
\]

is independently sufficient evidence that the current supporting plane and
candidate phase set are incomplete. A nonnegative best-found value from a
finite search is not a stability proof.

This interpretation explains the upper/lower construction used by HELD2. It
does not replace HELD2's Appendix-C initialization, persistent
\(\mathcal M\), Step-4 upper LP, Step-5 lower solve, Step-6 selection, or
Steps 8--10.

The native manufactured workflow checks the same interpretation on a
two-column problem: the Step-4 value equals the analytic restricted-master
value, active-column reduced costs vanish, and a strict-negative reduced-cost
column lowers the master.

## Exact nonreactive HELD2 reduction

For a nonreactive strong-electrolyte problem, choose the independent conserved
coordinates so that the individual species totals are conserved together
with per-phase electroneutrality. On HELD2's existing modified-mole-fraction
chart, whose reconstructed physical mole fractions sum to one, let

\[
\mathbf x_0
\]

be the homogeneous feed state, let

\[
g_0=g(\mathbf x_0,V_0),
\qquad
\boldsymbol\lambda_x=
\nabla_{\mathbf x}g(\mathbf x_0,V_0),
\]

and absorb the supporting-plane intercept into the normalization multiplier:

\[
\lambda_s
=
g_0-\boldsymbol\lambda_x^T\mathbf x_0.
\]

Then the GREPE reduced cost of any trial state is

\[
\begin{aligned}
r(\mathbf x,V)
&=
g(\mathbf x,V)-\lambda_s-\boldsymbol\lambda_x^T\mathbf x\\
&=
g(\mathbf x,V)-g_0
-\boldsymbol\lambda_x^T(\mathbf x-\mathbf x_0),
\end{aligned}
\]

which is exactly the HELD2 tangent-plane distance implemented by
`evaluate_held2_tpd`.

The equality applies to off-pressure trial volumes. Mechanical stationarity
is required when selecting the homogeneous reference and accepting final
phases, but not to prove that a feasible trial lies strictly below the
reference tangent plane.

The native manufactured workflow owns an executable algebraic check of this
identity. The public Perdomo gate protects the complete nonreactive reduction.

## Pricing search policy and shared mechanism

Reactive phase pricing reuses the successful HELD2 numerical policy:

1. Construct a complete charge- and mass-feasible composition chart for the
   declared phase-family domain.
2. Search jointly over the independent composition coordinates and one
   normalized log-volume coordinate mapped into the Provider's
   composition-dependent physical bounds.
3. Charge the budget in actual Provider callbacks, including bounds queries
   and state evaluations.
4. Stop early only for an independently reevaluated, feasible strict-negative
   reduced-cost state.
5. Keep the raw off-pressure minimum separate; pressure-refine only promising
   states needed for initialization.
6. Return indeterminate on required Provider failure, invalid state,
   incomplete search, or failed recertification.

Complete pressure-root enumeration remains mandatory when a homogeneous feed
reference must be selected. It is a different scientific decision from
searching for one feasible negative reduced-cost state.

HELD2 and GREPE share only a small formulation-neutral deterministic search
mechanism: normalized box exploration, stable candidate ordering, callback and
budget accounting, and termination plumbing. HELD2 owns its modified
composition and TPD; GREPE owns its conserved-coordinate chart and reduced
cost. Each caller independently recertifies its witnesses. Neither controller
calls or translates the other's scientific search.

DIRECT-L is used only when the chart maps the complete assigned finite-search
domain into a true box. A different phase family may require another
deterministic feasible chart, but it may not introduce an infeasible penalty
box and call the result globally complete.

## Revised GREPE workflow

### Stage 0: compile and audit

1. Validate the ordered Provider components, charges, molar masses,
   fingerprint, \(T,P\), and phase families.
2. Convert standard-state inputs to the Provider Helmholtz basis and compile
   reactions, invariants, totals, support, and the accessible face through the
   homogeneous owners.
3. Verify mass, charge, cycles, null-space identity, gauge invariance, and
   Provider domain/derivative contracts.
4. Stop before phase optimization on any unresolved hard contract.

### Stage 1: generate physical seed columns

Use the existing homogeneous reacting-phase solver to generate one or more
locally certified seeds for the one applicable Provider liquid family.
Initialize the homogeneous solves from every distinct mechanically stable
feed-density root and retain distinct locally certified reacting states.
Previous continuation states and source-bound seeds may be added, but they
remain initialization evidence.

In the nonreactive strong-electrolyte reduction, dispatch to the existing
HELD2 Steps 1--10 rather than entering a parallel GREPE seed path.

### Stage 2: establish master feasibility

Solve a physical Phase-I master over the available reactive phase columns.
If artificial balance residuals remain, use phase-family feasibility pricing
to search for another physical column. In search mode, failure to find such a
column is unresolved feasibility, not certified model infeasibility.

### Stage 3: energy master and reactive phase pricing

Solve and independently audit the restricted energy master. Use its conserved
multipliers in joint composition/log-volume pricing, append each independently
recertified strict-negative distinct column under an immutable insertion ID,
and repeat until a column is added or the declared finite search completes
without one.

No finite-search exit is relabeled as guaranteed globality.

### Stage 4: select the candidate phase topology

Keep three decisions separate:

1. persistent column identity uses phase family, true composition, and molar
   volume under a tight numerical-copy tolerance;
2. candidate distinctness uses the formulation-owned independent physical
   axes; and
3. final duplicate merging uses physical composition and molar volume.

Tolerance proximity is not transitive. Visit immutable IDs in stable order
and form a greedy maximal pairwise-distinct set. Do not use union-find or
nearest-composition reconstruction.

For the first slice, select at most two distinct liquid candidates. This is an
output capability bound, not a caller-supplied topology and not a reason to
discard a third discovered state as though the problem were solved. Evidence
for more than two material phases returns a typed scope-exceeded result.

### Stage 5: final simultaneous reactive multiphase Ipopt NLP

Instantiate the semi-infinite master formulation above on the selected finite
phase set and solve all phase amounts, normalized species amounts, and volumes
simultaneously with Ipopt.

Reaction extents and duplicate mass-action constraints are not added.
Stationarity with respect to every admitted invariant-preserving displacement
establishes reaction equilibrium. Common conserved-component multipliers,
phase charge constraints, and volume stationarity establish transfer,
electrochemical, and mechanical equilibrium.

The final NLP must:

- start from the actual feasible master/HiGHS phase amounts and compositions;
- use exact Provider tensors, Equilibrium chain rules, analytic gradients,
  exact Lagrangian Hessians, and the one derived scaling contract;
- keep phase amounts nonlogarithmic so inactive phases can reach zero;
- never clip or retire a phase or species because its value is small;
- retain HELD2's linear ordinary coordinates for its nonreactive reduction and
  reuse the homogeneous positive electroneutral chart for reactive chemistry;
- reserve bounded logarithmic refinement for formulation-owned trace
  coordinates only;
- fail closed on a nonconverged Ipopt result even if a finite iterate exists;
  and
- recompute balances, pressure, charge, KKT conditions, phase identity, and
  Provider-domain evidence independently.

There is no fixed physical \(\theta_{\min}>0\). A candidate phase may be
retired only after a positive lower-bound multiplier supplies KKT inactivity
evidence; retire at most one and solve the complete reduced NLP again. One
cold solve of an identical candidate problem may follow a failed warm solve,
with both attempts counted.

### Stage 6: final support and trace refinement

Refine only retained positive trace coordinates. Do not convert a positive
trace into an exact zero. Merge only numerical duplicate phases and rerun the
complete NLP and physical certificates after every merge. Any support,
phase-set, bound, domain, or solved-state change invalidates downstream
evidence.

### Stage 7: final state-matched repricing

1. Freeze the exact final local state and its evidence key.
2. Extract conserved-component multipliers from that state.
3. Rerun every declared phase-family pricing search with those multipliers.
4. Independently recertify every feasible negative state.
5. If a negative phase exists, append it and return to Stage 3.
6. If the declared finite search completes without a negative state, retain
   `globality_certificate="not_guaranteed"`.

This is a future reactive GREPE requirement. It does not silently add a second
postsolve search to the frozen nonreactive HELD2 controller.

### Stage 8: final certification

Accept only when the final NLP state, support, phase topology, multipliers,
boundary evidence, and postsolve pricing evidence have matching validity keys
and every applicable numerical and physical gate passes.

## Evidence validity and invalidation

Evidence is reusable only when every input that defines its mathematical
problem is unchanged.

### HELD2

HELD2 retains its ordered immutable candidate-ID vector as the exact Step-8
problem key. Any changed candidate problem or uncertified prior result forces
a new solve; Step 9 is always rerun. A controller-wide revision would
needlessly invalidate unchanged Step-8 problems.

### Reactive GREPE

A future GREPE result carries one monotonic `state_revision`. The revision
increments after any change to:

- structural or positive/dormant support;
- active phase IDs or order;
- phase amounts, compositions, or volumes;
- phase-family domain or numerical bounds;
- compiled chemistry or standard-state reference;
- Provider identity or fingerprint; or
- final NLP solution.

Postsolve pricing stores the exact revision plus support and multiplier
checksums. A mismatch invalidates the pricing evidence and forces repricing.
Checksums detect stale evidence; they do not replace physical comparison,
certification, or immutable scientific inputs.

## Evidence gates

### Gate A: frozen standalone behavior

Before HELD2 or chemistry internals are shared, freeze and run:

1. the native manufactured Steps 1--10 workflow;
2. the complete Perdomo numerical matrix;
3. native/Python diagnostic parity;
4. unchanged named tolerances, resource budgets, and
   `globality_certificate="not_guaranteed"`; and
5. the executable HELD2-TPD/GREPE-reduced-cost identity;
6. the complete chemistry compiler, support, reaction-basis, species-order,
   sensitivity, and fail-closed tests; and
7. the public ideal and installed-Provider homogeneous chemical sentinels.

The HELD2 Step-5 original-coordinate KKT audit and validated Step-8 Farkas
evidence must be completed as an isolated, regression-gated prerequisite before
the shared search mechanism is extracted.

### Gate B: coupled manufactured evidence

The private reactive implementation requires:

- reaction-basis, species-order, feed-scale, and gauge invariance;
- a manufactured analytic reactive two-phase equilibrium;
- the source-derived Ascani \(A+B\rightleftharpoons C\) topology case;
- a narrow negative reactive-pricing basin;
- multiple-density-root reference selection;
- exact feasible master-to-NLP initialization;
- inactive-phase KKT retirement followed by a complete reduced re-solve;
- a stale-revision rejection test;
- negative final repricing that returns to phase generation; and
- exact reduction parity for constrained homogeneous and nonreactive
  phase-discovery variants.

Search budgets are expressed in actual Provider callbacks and counted local
solves. Wall time is recorded but not used as a hardware-dependent correctness
assertion. No stage may hide a complete pressure-envelope solve inside every
joint-search trial.

### Gate C: neutral reactive-LLE physical evidence

Before public exposure or a capability claim, reproduce the
Ascani--Senina acetic-acid/1-pentanol/pentyl-acetate/water subject through an
exact installed Provider artifact for Ascani's four-species PC-SAFT model
hypothesis. Record that \(K_a=43.99\) was calibrated by Ascani from one Senina
homogeneous equilibrium composition using PC-SAFT activity coefficients, so
that datum is calibration rather than validation evidence. Also record that
Senina used less than \(2\ \mathrm{wt}\%\) aqueous HCl while Ascani omitted the
catalyst from the four-species calculation as an approximation. Compare
calculated endpoints against the nine tabulated experimental reactive
tie-lines at their source \(T,P\). Construct any flash feed from a declared
convex weight applied to the two tabulated endpoint compositions and retain
that provenance;
do not claim reproduction of unpublished Ascani calculated endpoints. Record
phase compositions, phase amounts, volumes, reaction residuals, transfer
residuals, pressure residuals, final reduced costs, Provider calls, and wall
time in an immutable Validation campaign.

### Gate D: electrolyte extension

Electrolyte reactive LLE begins only after Gate C and requires a separate
source-complete installed-Provider case. It reruns every neutral and
nonreactive regression and additionally proves phase-specific charge,
electrolyte standard-state conversion, charged trace behavior, and charged
coordinate invariance.

## Implementation order

1. Publish the HELD2 paper companion and freeze the standalone gates.
2. Close the two isolated HELD2 evidence gaps.
3. Extract the formulation-neutral deterministic joint-search mechanism and
   prove unchanged HELD2 results and budgets.
4. Add private typed `equilibrate` variants that delegate limiting cases.
5. Implement GREPE Stages 0--4 for the at-most-two-liquid slice.
6. Implement Stages 5--8 with exact derivatives and revision-bound evidence.
7. Complete Gate B, then Gate C in Validation.
8. Consider public `equilibrate` only after stable Gate-C evidence.
9. Treat electrolyte reactive LLE as the separately evidenced Gate-D
   extension.

## Rejected designs

GREPE does not:

- speciate once and treat a subsequent phase-only solve as coupled
  equilibrium;
- copy the HELD2 controller into a generic reactive controller;
- run a complete pressure envelope at every trial composition;
- use a fixed positive phase-mass lower bound;
- prune phases or species by magnitude;
- substitute logarithmic volume distance for Provider packing fraction in
  HELD2 Eq. (66);
- use one tolerance or clustering rule for all phase identities;
- treat Ipopt as a global optimizer;
- treat finite search exhaustion as global stability; or
- depend on Provider internals, a sibling source checkout, or a selectable
  derivative backend.
