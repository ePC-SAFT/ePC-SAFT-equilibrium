# Perdomo HELD2 Strong-Electrolyte Phase Equilibrium

Status: canonical design with private integrated Stage-I/II/III implementation

Authority effect: none

## Status and authority

This document is the package-local scientific and numerical owner for the
Perdomo HELD2 strong-electrolyte formulation. Current `main` contains one
private callback-driven Stage-I/II/III controller. The later
installed public-dispatch and reference-hardening runtime is archived at
`archive/held2-pre-strategy-2026-07-21` as non-production evidence.
Organization doctrine revision 3 and
[the package authority map](../phase-equilibrium.md) govern ownership and the
claim boundary. The canonical numerical decomposition, landed task state, and
guarded execution order are in the
[HELD2 solver-strategy implementation plan](../plans/2026-07-21-perdomo-held2-solver-strategy.md).

The archived subject exposed the controller through the existing public
`tp_flash` operation for qualifying installed Provider capability tables.
Current `main` does not expose that electrolyte dispatch. Neither the archived
public exposure nor this canonical design is capability admission. The only
accepted Equilibrium capability remains the receipt-bound pure-component
saturation slice.

## Scope and nonclaims

The formulation covers fixed-temperature, fixed-pressure, nonreactive
fluid-phase equilibrium for mixtures containing strong, fully dissociated
electrolytes. The archived WIP exercised this scope against installed Provider
contracts. A qualifying topology has at least one molecular species and at
least two charged species and is supported by a complete installed Provider
phase/domain/derivative contract.

The formulation does not cover:

- weak-electrolyte speciation or chemical reaction;
- solid, hydrate, or precipitate equilibrium;
- reactive or coupled phase-chemical equilibrium;
- caller-supplied phase counts, phase guesses, seeds, or solver settings;
- a chemistry-named or case-specific controller branch;
- copied Provider equations or parameter fitting;
- a caller-selectable or silent fallback solver, or numerical production
  derivatives; or
- deterministic global optimization or a global-stability certificate.

Any resumed public implementation must reuse `tp_flash`, `TpFlashResult`,
`HeldDiagnostics`, and `FlashError`. It must also reuse the existing native
module, target, `ProviderContext`, HELD2 result/diagnostic owner, and the
Stage-II/III Ipopt adapters where the solver-strategy plan assigns Ipopt.
Current `main` uses pinned NLopt 2.11.0 for DIRECT-L and pinned HiGHS 1.15.1
for Problem (64). They remain private implementation dependencies, not public
solver choices.

## Scientific sources and notation

Perdomo et al. (2025), Sections 2.2, 3.2, 3.5, and 4, including Algorithm 1
and Problems (63)--(67), is the primary HELD2 source. Pereira et al. (2012),
Sections 3.1--3.3, is the source for the parent three-stage HELD structure.

Perdomo reports an in-house tunnelling method for Stage I, NAG E04UCA for the
nonconvex lower Problem (65), and NAG E04MFA for the upper Problem (64). The
paper fixes the thermodynamic problems and controller logic, not those library
identities. The package strategy deliberately maps the same responsibilities
to DIRECT-L, exact-Hessian Ipopt, and HiGHS respectively. This implementation
choice must not be described as a verbatim reproduction of Perdomo's numerical
software.

The permanent-lab paper Markdown and `docs/latex/equilibrium_formulations.tex`
are retained provenance. The lab M4 Perdomo doctrine is a concise scientific
summary, but its old route and implementation-status fields predate this
package. Migration D-024 through D-026 records implementation chronology;
Validation owns source ledgers and installed-artifact evidence. Neither is the
equation owner.

Let:

- `C` be the number of physical species;
- `z_i` be the integer charge of species `i`;
- `E` be the eliminated charged species;
- `C^(E)` be the species set excluding `E`;
- `q` be the dependent retained molecular species;
- `C^(EC)` be the `C - 2` independent retained species;
- `n_i`, `x_i`, and `V` be physical amounts, fractions, and volume; and
- `n_bar_i`, `x_bar_i`, and `mu_bar_i` be the Perdomo modified quantities.

Quantities divided by `RT` are dimensionless. Phase superscripts are omitted
where the relation applies independently to every phase.

## Modified-mole coordinates

### Eliminated charged species

Perdomo Eq. (9) imposes local electroneutrality by eliminating one charged
species:

```text
n_E = -(1 / z_E) sum_{i != E}(z_i n_i).
```

The implementation chooses `E` from the charged species with largest absolute
charge, as described in Perdomo Stage-I Step 1. This keeps
`|z_i / z_E| <= 1`. Equal-charge or otherwise singular coordinate factors fail
closed rather than silently selecting a different formulation.

### Modified amounts and fractions

Perdomo Eq. (23) defines, for each retained species,

```text
n_bar_i = (1 - z_i / z_E) n_i,
n_bar   = sum_{i in C^(E)} n_bar_i,
x_bar_i = n_bar_i / n_bar.
```

Molecular species have `z_i = 0` and are unchanged. The retained modified
fractions sum to one. One retained molecular species is dependent:

```text
x_bar_q = 1 - sum_{i in C^(EC)} x_bar_i.
```

The optimization therefore has `C - 2` independent modified-composition
coordinates. The physical lift reconstructs the dependent modified fraction,
undoes the modified scaling, reconstructs `n_E`, and verifies:

- finite nonnegative physical amounts and fractions;
- modified and physical normalization;
- ordinary and modified material balance; and
- per-phase electroneutrality.

The back-lift must return the original independent modified state within the
declared numerical contract. A nonpositive eliminated amount, nonpositive
dependent fraction, singular factor, or failed charge/balance reconstruction
is a domain failure.

### Modified potentials

Perdomo Eqs. (28)--(31) give the charge-constrained and scaled modified
potentials:

```text
mu_i^(E)  = mu_tilde_i - (z_i / z_E) mu_tilde_E,
mu_bar_i  = mu_i^(E) / (1 - z_i / z_E).
```

For active phases and retained interior coordinates, equilibrium requires
equality of all `C - 1` retained modified potentials, including the dependent
coordinate recovered through phase-fraction stationarity. It does not require
equality of unconstrained individual-ion chemical potentials. An individual
ionic comparison is meaningful only under an explicit phase-electric or
Galvani-potential convention.

## Provider and derivative contract

### Installed capability boundary

The controller consumes one installed `epcsaft.native_sdk.v1` capsule. Dispatch
requires a compatible strong-electrolyte capability table, finite table sizes,
component identifiers and order, charges, and the exact model fingerprint.
The Provider owns:

- complete ideal-mixing and active residual Helmholtz phase values;
- gradients and Hessians in component-amounts-then-volume order;
- packing-fraction value, gradient, and Hessian;
- temperature and total-ion source-domain metadata; and
- composition-specific molar-volume bounds.

The public feed must be finite, positive, normalized, and identical in length
and order to Provider metadata. Unsupported multicomponent non-electrolyte
tables and incomplete electrolyte tables fail before search.

No Provider call may receive a state outside the source, composition, packing,
or volume domain. Raw Provider bounds remain available for certification and
diagnostics. When strict endpoints are not evaluable, solver bounds are moved
inward only to the next representable value whose round trip is demonstrably
strict; no arbitrary epsilon is introduced.

### Exact coordinate transformations

Let `w = (n_1, ..., n_C, V)` be Provider coordinates and `y` any Equilibrium
coordinate chart. For a scalar phase objective `f(w(y))`, Equilibrium applies
the exact chain rules

```text
gradient_y(f) = J^T gradient_w(f),

Hessian_y(f)  = J^T Hessian_w(f) J
                + sum_k gradient_w(f)_k Hessian_y(w_k),
```

where `J = partial w / partial y`. These are coordinate transformations of the
single Provider tensor owner, not alternate thermodynamic derivatives.

The controller uses related charts at different numerical layers:

- the Stage-I and Stage-II discovery layers use bounded independent modified
  compositions, while the nested pressure-root service searches `log(V)`;
- Stage-II local refinement uses feasible-simplex coordinates plus
  `q = log(eta)`; and
- Stage-III phase fractions, per-phase modified compositions, and per-phase
  `q`, with extensive material bookkeeping.

For `q = log(eta)`, volume is the implicit solution of

```text
log(eta_provider(u, V)) - q = 0.
```

Equilibrium solves the value outside derivative recording and derives exact
first and second sensitivities from the Provider packing gradient and Hessian.
The chart fails closed when the packing root is unbracketed, incomplete, or
singular. Production finite differences are forbidden; centered differences
are test-only derivative evidence.

## Homogeneous reference selection

At the fixed feed and `T,P0`, the homogeneous reference minimizes the reduced
one-phase objective among mechanically stable pressure roots:

```text
g_bar(V) = A(T, n_feed, V) / (RT) + P0 V / (RT).
```

The deterministic search:

1. uses the exact Provider log-volume domain at the fixed feed;
2. evaluates the declared finite scan without silently skipping failures;
3. brackets sign-changing pressure residual roots;
4. searches stationary points of the pressure residual so tangential/double
   roots cannot be classified as absent solely from same-sign endpoints;
5. refines every detected root from Provider phase gradients and Hessians;
6. classifies mechanical stability from the fixed-composition pressure-volume
   derivative, requiring the stable sign `dP/dV < 0`; and
7. selects the stable root with the lowest finite `g_bar`.

The search fails closed if evaluation or refinement is incomplete, a root is
mechanically marginal, a root lies at the Provider domain boundary, or the
lowest stable objectives are tied within the existing selection contract.

Diagnostics retain scan points and intervals, detected roots, stable roots,
stationary points, tangential roots, marginal roots, boundary roots, objective
ties, and evaluation/refinement failures. Current finite detection reports
`root_completeness="not_proven"`; this is independent of solver status,
physical acceptance, and globality.

## Stage I: stability

Stage I transforms the feed into the Perdomo modified coordinates and evaluates
the modified tangent-plane distance in modified composition and volume. The
reference state has zero TPD by construction. Trial states must satisfy the
Provider source, composition, charge, packing, and volume domains.

The archived WIP used a declared deterministic finite Ipopt multistart in
modified composition and volume. That replaced runtime and its baseline
fixture are removed. The implemented Stage-I strategy reduces volume at each
trial composition through the deterministic pressure-root service and defines
the lower TPD envelope over eligible strict-stable pressure roots. It then uses
NLopt `NLOPT_GN_DIRECT_L` over the bounded `C - 2` modified-composition chart.
The lower envelope may be nonsmooth where the selected density branch changes,
so Stage I does not require a global derivative or a local NLP polish to accept
an independently certified negative witness.

The private strategy retains its declared DIRECT-L evaluation budget and
completion accounting. The implementation did not change a physical
tolerance.

Terminal precedence is asymmetric:

- A finite, feasible, materially distinct, independently confirmed point with
  `TPD <` the existing negative threshold establishes homogeneous instability.
  It does not need to be the optimum of every local start. Unrelated incomplete
  starts remain visible as partial search evidence.
- `no_negative_found` requires every declared start needed by that claim to
  complete and certify. It means only that the declared finite search found no
  negative point.
- Without a confirmed negative witness or complete no-negative accounting, the
  result remains indeterminate or resource-limited.

A Provider failure or incomplete pressure-root evaluation is never converted
to a large finite objective penalty. It remains explicit search evidence. Such
a failure prevents a completed no-negative claim, although it does not erase a
separately certified negative existence witness.

`negative_tpd` is an existence claim, not complete phase enumeration or a
globally optimal split.

## Stage II: dual cutting and candidate discovery

### Upper and lower problems

Stage II implements the Perdomo semi-infinite dual strategy. The upper problem
is a finite linear program in the modified multipliers over the complete
retained cut pool. The canonical production owner for this Problem-(64) LP is
HiGHS, with Equilibrium independently recomputing its objective, row
activities, primal residual, and dual residual. The analytic envelope
implementation remains a manufactured test oracle.

Each nonconvex lower Problem-(65) minimizes the reduced phase function in
modified composition and packing/volume at the current multiplier. Its
canonical strategy is global basin discovery followed by accurate local
refinement: continuation states, retained cut states, Stage-I witnesses,
source-independent physical seeds, deterministic Sobol samples, and one
stall-triggered DIRECT-L pressure-envelope exploration propose distinct
physical basins; exact-Hessian Ipopt then refines every retained
representative. SLSQP is not the default replacement because it remains local
and does not use the exact Provider-derived Hessian already available to the
Ipopt path.

Each certified improving lower terminal may add one cut. A certified improving
point is existence evidence: unrelated incomplete starts do not erase it.
However, an incomplete lower search is not a completed lower bound. Complete
accounting remains mandatory before declaring no improving cut or a final
finite-search gap.

The controller retains cut identity, bounds, multipliers, attempt completion,
Provider/domain failures, discovery source, physical basin identity, candidate
status, and bounded no-progress/resource exits. The full cut pool is not the
same concept as the current candidate phase set. Global exploration proposes
basins; it does not certify KKT conditions or prove that the global lower
minimum was found.

### Feasible-simplex and packing charts

Stage-II lower local-refinement solves use one private lower-shifted
feasible-simplex chart. It maps each retained physical basin representative
into an interior cube while preserving all independent lower bounds and the
dependent-fraction upper sum. Arbitrary Ipopt trial points therefore remain
inside the physical modified simplex before a strict Provider evaluation.

The chart supplies exact value, Jacobian, and Hessian chain rules. Terminal
acceptance is recomputed in original `u,q` coordinates. Internal chart KKT is
not a physical certificate.

Finite-barrier cube duals are pulled through the exact chart Jacobian and
deterministically decomposed against the physical independent lower/upper
bounds and modified-simplex inequality. Certification requires:

- correct physical multiplier signs;
- original-coordinate stationarity under the existing threshold;
- primal feasibility;
- physical slack-multiplier complementarity;
- finite Provider/domain evidence; and
- negligible pullback reconstruction residual.

Solver termination remains a separate per-attempt status. A non-success Ipopt
status is neither automatically accepted nor automatically a numerical
failure.

### Source-faithful candidate selection

Perdomo Step 6 runs immediately after every certified improving-cut insertion
and on the no-improvement path. Candidate predicates use the named package
tolerance contract for:

- dual-value proximity;
- fixed-volume composition-gradient agreement with the multiplier over the
  source `I^m` retained component set, using the zero-safe threshold
  `1e-8 + 1e-7*max(abs(gradient), abs(multiplier))`;
- packing/pressure stationarity; and
- phase distinctness.

Candidate gradients for this test are fixed-volume physical gradients. The
`q`-chart gradients belong to the lower TNLP and include a non-negligible chain
term; substituting them changes the source predicate.

Numerical duplicate basins use `1e-7` physical-composition and log-volume
distances. Two candidates are confidently distinct when physical composition
**or** log volume differs by more than `1e-5`. The band between those tests is
unresolved and fails closed. Candidate identity remains independent of the
tighter cut-growth duplicate decision.

Stage III may start only from the ordinary controller-owned candidate vector.
There is no caller-injected witness or case-specific candidate escape hatch.

## Stage III: direct total-free-energy refinement

For `m_p` candidate phases, Perdomo Problem (67) is the direct reduced total
free-energy minimization

```text
min sum_alpha phi_alpha
        [A_bar(T, v_alpha, x_bar_alpha) / (RT)
         + P0 v_alpha / (RT)]
```

subject to

```text
sum_alpha phi_alpha x_bar_i,alpha = x_bar_feed_i
    for every independent modified component i,

sum_alpha phi_alpha = 1,

phi_alpha >= 0,

x_bar_alpha in the modified-composition domain,

eta_alpha in [eta_lower, eta_upper],
```

plus the declared candidate neighborhoods. The omitted dependent modified
balance follows from normalization; the physical lift and electroneutrality
recover ordinary balances.

The implementation is one general exact-Hessian Ipopt `m_p` TNLP. It computes material and
per-phase simplex constraints algebraically before Provider evaluation. Every
phase has one linear dependent-fraction inequality with exact constant
Jacobian and zero Hessian contribution. Its multiplier contributes to the full
Lagrangian-gradient certificate.

Objective, gradient, and Hessian callbacks reject out-of-simplex trial points
as recoverable Ipopt evaluation failures without calling the strict scientific
evaluator. Recoverable trial rejections do not permanently populate a fatal
callback error if Ipopt later reaches a valid terminal. A Provider failure at a
declared-valid state remains separately observable; exhaustion without a valid
terminal fails closed.

Candidate sets that are infeasible, incomplete, nonconverged, or physically
rejected return a Stage-II feedback status under the bounded private
controller policy. Duplicate phases are merged only through the
evidence-backed lifecycle. A small phase amount alone does not authorize
KKT-inactive retirement. Retirement requires a phase-amount bound, correct multiplier sign,
complementarity, a non-descending reduced derivative, remaining-balance
feasibility, one-at-a-time retirement, and an active-set re-solve.

The source requires bound complementarity and logarithmic refinement for
trace-bound components. The controller detects that condition and fails closed
with `complementarity_refinement_required`; it does not yet
implement the final logarithmic trace refinement. Trace components are not
passed to Ipopt near `1e-300` and must not be accepted through an unconditional
interior modified-potential equality.

Step 9 applies Perdomo Eq. (68) to independent controller quantities:

```text
free_energy_gap = same_major_stage_ii_UBD - stage_iii_total_free_energy
```

The installed controller passes the Problem-(64) upper bound from the same
Stage-II major that produced the candidate set. Stage III computes the
Problem-(67) total free energy independently, retains both values and their
provenance, and accepts only when the named free-energy-gap gate passes.
Missing evidence, a default numeric value, or a locally converged Ipopt state
with a failed gap returns to Stage II. The manufactured oracle computes its
upper bound by independent enumeration and includes a perturbed-bound rejection
case. Installed Provider directional-gradient and Hessian-vector tests exercise
the same generic Stage-III derivative owner.

An Ipopt terminal may meet its unscaled local target while its volume
stationarity maps to a relative pressure residual just outside the physical
pressure gate. After the active set is finalized, Stage III therefore applies a
safeguarded Newton correction in log-volume using the exact Provider-derived
pressure residual and derivative. Only phase-activity-certified phases are
corrected; trace or inactive phases remain under their separate fail-closed
logic. The corrected point is accepted only after KKT and complementarity are
recomputed and all unchanged Step-9 gates pass.

## Certification and status axes

An accepted multiphase result requires all applicable existing certificates:

- ordinary and modified material balance;
- physical normalization and per-phase electroneutrality;
- positive finite phase amounts and valid phase identities;
- Provider source, composition, packing, and volume domains;
- imposed pressure in every active phase;
- every retained modified-potential condition, including the dependent
  coordinate from phase-fraction stationarity;
- original-coordinate KKT stationarity and complementarity;
- distinct active phases and no collapsed duplicate accepted as LLE;
- the source Eq. (68) same-major Stage-II/Stage-III free-energy gap; and
- independent confirmation where the controller requires it.

Potential-gap checks use mixed absolute/relative scaling rather than division
by a near-zero modified potential. Solver success never substitutes for these
physical gates.

The existing result owner keeps these independent axes:

| Axis | Question |
| --- | --- |
| Artifact/input | Is the installed SDK and resolved input complete and in scope? |
| Solver | How did the solver or controller terminate? |
| Numerical | Do derivatives, original-coordinate KKT, feasibility, complementarity, and required search accounting pass? |
| Physical | Do balance, charge, pressure, potentials, phase identity, stability, and domains support the result? |
| Search completeness | Did the declared finite search needed for this terminal claim complete? |
| Root completeness | Were pressure roots completely established? Current finite detection is `not_proven`. |
| Predictive | Does an accepted physical result support an external comparison? |
| Globality | Was global phase stability proven? For this controller: `not_guaranteed`. |

`passed`, `failed`, and `not_adjudicated` describe the evidence actually
available on the solver, numerical, and physical axes. Missing evidence must
not be represented by a fake zero, false, or failure. A public `FlashError`
preserves the native terminal diagnostics.

## Archived evidence and future gate

The archived installed public Perdomo Table 3 NaCl/water tracer is accepted as a
homogeneous one-phase package result. Its reference search retains three
detected pressure roots, two mechanically stable roots, the selected
lowest-objective root, complete Stage-I start accounting, and
`root_completeness="not_proven"`. The comparison is explicitly cross-EOS:
Perdomo published SAFT-gamma-Mie outputs while the installed run uses ePC-SAFT.
The resulting source-topology disagreement is evidence, not an electrolyte-LLE
admission or a Perdomo numerical reproduction.

D-026 identified one source-complete, installed, public, genuinely two-liquid
ePC-SAFT case as the next critical gate. The private Stage-I/II strategy now
exists on `main`; that work does not authorize public electrolyte dispatch,
Stage-III completion, Validation, or promotion. A separate bounded assignment
and compatible corrected Provider artifact must rebind the installed gate.
Any resumed public controller must, without a chemistry-specific branch or
retuning:

1. establish Stage-I instability;
2. form eligible Stage-II candidates with truthful search accounting;
3. enter the existing general Stage III;
4. return two physically accepted distinct liquid phases; and
5. retain solver, numerical, physical, predictive, root-completeness, search,
   and globality evidence separately.

Until that gate passes and a distinct acceptance process admits it, the
electrolyte two-liquid capability remains unproven and unpromoted.

## Provenance and deliberately excluded material

The following remain references rather than package doctrine:

- the permanent-lab M4 generalized architecture and algorithm summaries;
- the M4 preservation manifest and admission registry;
- the dated M4 dashboard and issue queue;
- historical Ascani counterion-pair implementations and old `held2` names;
- Migration telemetry, exact WIP commits, artifact hashes, and failed runs;
- Validation case constants, extracted source tables, and campaign results; and
- the broad paper archive.

They are deliberately not copied here because this design owns the durable
formulation and package contract, not historical workflow state. Provider EOS
equations, parameters, and packing formulas remain exclusively Provider-owned.
