# Phase-Equilibrium Documentation Authority

Status: canonical package documentation

Authority effect: none. This document describes ownership and current evidence;
it does not admit, promote, or publish a capability.

## Authority and capability status

Organization doctrine revision 7 defines the ecosystem authority hierarchy.
Within this repository, `AGENTS.md` and `CONTEXT.md` govern package scope,
`ARCHITECTURE.yaml` records the machine-readable architecture, and the design
documents named below own the scientific and numerical contracts of individual
formulations. Plans and receipts record execution and evidence; they do not
replace a formulation owner.

| Formulation | Package design owner | Current capability state |
| --- | --- | --- |
| Pure-component saturation | [Pure-saturation slice](designs/2026-07-17-pure-saturation-slice.md) | Accepted only for the exact methane, ethane, and propane scope in `promotion-0018-equilibrium-pure-saturation-v1` |
| Neutral Pereira HELD | [Neutral HELD v1](designs/2026-07-17-neutral-held-v1.md) | Frozen local candidate; installed campaign retained as `NON_ADMISSION`; controller redesign deferred |
| Strong-electrolyte Perdomo HELD2 | [Paper-faithful Steps 1--10](designs/2026-07-24-held2-paper-algorithm.md), [publication companion](designs/2026-07-29-held2-publication-algorithm.md), with [earlier design provenance](designs/2026-07-21-perdomo-held2.md) | Public development dispatch over one native core with fail-closed installed evidence; no admitted electrolyte LLE capability |
| Homogeneous reacting phase | [Private reacting-phase kernel](designs/2026-07-21-private-reacting-phase-kernel.md), [GREPE homogeneous layer](designs/2026-07-27-grepe-homogeneous-chemical-layer.md), [deterministic basin-search contract](designs/2026-07-30-homogeneous-chemical-basin-search.md), [Provider-model continuation](designs/2026-08-02-provider-model-continuation.md), and [installed-artifact handoff](designs/2026-07-30-chemical-basin-search-validation-handoff.md) | Public typed local value operation with optional conditioned state-input sensitivities and bounded Provider-model continuation; deterministic finite basin search for the lowest observed certified strict local minimum; no predictive, phase-stability, or global admission |
| Superseded fixed two-phase route | [Historical fixed-route design](designs/2026-07-17-neutral-two-phase-tp-flash.md) | Removed without alias; retained only as provenance |
| Ascani counterion-pair electrolyte equilibrium | No current runtime design | Closed future formulation; historical lab evidence only |
| Coupled multiphase chemical equilibrium | [GREPE reactive phase equilibrium](designs/2026-07-28-grepe-reactive-phase-equilibrium.md), [implementation plan](plans/2026-07-29-grepe-reactive-lle-plan.md), and [benchmark review](designs/2026-07-29-reactive-lle-benchmark-review.md) | Normative future implementation contract for an at-most-two-liquid first slice; no public schema or runtime route |

The only accepted capability is `pure-component-saturation-v1`. A public
symbol, a local candidate, an installed campaign, or a converged local solve is
not an authority receipt.

## Public operations

The package exports three equilibrium operations:

- `saturation` owns the accepted, bounded, pure-component local saturation
  boundary.
- `tp_flash` is the sole mixture flash surface. It dispatches the reviewed
  neutral binary fingerprint to neutral HELD and qualifying installed strong-
  electrolyte Provider capability tables to the non-admitted HELD2
  development route.
- `chemical_equilibrium` is the sole fixed-`T,P`, single-homogeneous-phase
  reaction-equilibrium value surface. Typed ideal-gas and installed-Provider
  phase inputs route to the same compiler, Ipopt owner, and certificate model.
  `ChemicalEquilibriumProblem`, `ChemicalEquilibriumConstant`, and
  `ChemicalStandardState` carry physical inputs and provenance;
  `ChemicalEquilibriumResult`, `ChemicalEquilibriumDiagnostics`, and
  `ChemicalEquilibriumError` carry the admitted value or typed failure.
  The declared strict-interior amount floor and admissible packing-fraction
  interval are application-owned formulation/domain bounds; neither exposes
  an optimizer backend, start, iteration limit, or Ipopt option.

`tp_flash` does not accept a phase count, caller seeds, solver settings,
backend selection, or a case-specific mode. Its public result owners are
`TpFlashResult`, `HeldDiagnostics`, `HeldPerformanceDiagnostics`,
`HeldStepTiming`, and `FlashError`. HELD2 performance diagnostics expose
already-measured native work and mechanism evidence without making solver
configuration public. The removed
`two_phase_flash` route and result family have no compatibility alias.

## Formulation owners

### Pure saturation

The pure-saturation design owns its log-density/log-pressure feasibility
problem, exact Provider transformations, local mechanical and phase-separation
checks, retained source anchors, and its explicit lack of phase-discovery or
global-stability evidence.

### Neutral Pereira HELD

The neutral design and plan remain frozen records of the implemented Stage I,
binary Stage II, direct-total-free-energy Stage III, public cutover, and later
`NON_ADMISSION` controller-lifecycle finding. The post-validation redesign is
documented but deferred; it is not current runtime behavior.

### Strong-electrolyte Perdomo HELD2

The Perdomo design owns the eliminated-ion modified-mole formulation, installed
Provider domain and derivative contracts, homogeneous reference selection,
modified-coordinate stability search, complete-cut Stage II, general candidate
set Stage III, and formulation-specific certificates. Perdomo modified moles
must not be replaced by or conflated with Ascani counterion-pair coordinates.
The finite \(10^{-10}\) modified-composition floor is a search regularization,
not a physical equilibrium minimum. Steps 1--9 retain linear modified
composition and material-balance coordinates. Step 10 alone may refine a
strictly positive charged trace component below that floor in
\(\log_{10}x_i\); every upper, closure, eliminated-ion, and Provider-domain
constraint remains active, and molecular lower bounds cannot be bypassed.
The refined phases are reconstructed and recertified in ordinary and modified
linear balances. Provider evaluation and result serialization preserve the
resulting positive trace values rather than clipping them to the search floor.
For finite Provider `total_ion_mole_fraction_max`, the Stage-I cube-to-modified
composition map enforces that immutable source ceiling before requesting a
pressure envelope. In Perdomo coordinates the total physical ion fraction is
the sum of the retained charged modified fractions, so this is an exact
coordinate-domain restriction rather than an EOS approximation or tolerance.
Malformed or infeasible ceilings fail before global exploration; `NaN` retains
the ordinary modified-simplex domain declared by a Provider with no ion cap.
The linked implementation plan assigns deterministic pressure-root
enumeration to feed-reference density topology, DIRECT-L to the joint
modified-composition/log-volume Stage-I search, HiGHS to the Stage-II upper
LP, deterministic capped-multistart exact-Hessian Ipopt to Step 5, and
exact-Hessian Ipopt to Problem (67) in Step 8. One typed
`run_held2_algorithm` state machine owns the closed Steps 1--10 transitions for
both installed and manufactured problems. Stage II has one major loop and one
Eq. (66) decision owner. Each major retains its upper-solve identity, `UBD`,
multipliers, active cuts, local attempts, pressure branches, and certificate
partitions. Replaced HELD2 runtime routes and baseline fixtures are deleted;
focused manufactured numerical oracles call the same transition and step
owners. None of these internal owners is a caller-selectable backend; public
electrolyte dispatch remains non-admitted development behavior.

The current charged Steps 1--10 contract is the
[HELD2.0 Paper Algorithm Specification](designs/2026-07-24-held2-paper-algorithm.md),
with its installed evidence exercised by the Perdomo regression matrix in
`tests/test_perdomo_held2_trace.py`.
The older
[HELD2 Installed Completion](designs/2026-07-22-held2-installed-completion.md)
and its implementation plan are superseded investigation provenance. They
begin from the retained Stage-II-indeterminate artifact and do not reinterpret
it as an equilibrium result. Current source has classified the
artifact's unit-cube callback failure as bounded binary64 contact, preserves
ordered trial/accepted-iterate and upper-LP evidence, and reaches successful
local Ipopt termination for that exact start. The corrected Problem-(65) owner
now follows Provider-certified pressure branches in a composition-only local
solve, reconstructs exact reduced derivatives through the volume Schur
complement, rejects local values above the same-major upper bound, and keeps
dual-cut eligibility separate from Eq. (66) candidate eligibility.
Composition-rich vertex seeds expose both aqueous-rich and organic-rich
basins. Manufactured Step-6 candidates feed the generic Step-8 NLP and pass
the complete Step-9 physical certificate. Step 9 now computes Perdomo
Eq. (68) as the same-major Problem-(64) upper bound minus the independently
solved Problem-(67) total free energy; the result retains both values, the
signed gap, and its provenance. Missing or failed gap evidence returns to
Stage II even when Ipopt and all other physical checks pass. Exact installed
Provider directional-gradient and Hessian-vector tests cover the generic
Stage-III formulation. The private Perdomo Table-5 ePC-SAFT screening exposed a
premature local-search exit that executed only one declared Step-5 attempt per
major. The corrected 24-major/50-attempt profile reaches two same-major Eq.
(66) candidates in major 19, runs the generic Step-8 NLP, applies a bounded
exact-derivative pressure polish on the finalized active phases, and passes the
complete Step-9 physical and Eq. (68) certificates. This is private numerical
evidence under an unadmitted ePC-SAFT parameter hypothesis, not a
source-equivalent reproduction of Perdomo's SAFT-gamma-Mie calculation or a
globality proof.

### HELD2 live progress diagnostic

The installed development controller accepts an internal, nullable progress
observer. It receives already-computed reference roots, DIRECT-L evaluations,
HiGHS bounds, Ipopt iteration metrics, certificate results, and the Step-9
total-free-energy upper bound/objective/gap. The default path supplies no
observer and remains
silent. Observer failures are swallowed and cannot change solver decisions,
structured results, tolerances, budgets, or the
`globality_certificate="not_guaranteed"` label.

The public HELD2 result carries the same completed-work evidence as typed
`HeldPerformanceDiagnostics`: every native Step invocation, Provider and
optimizer counts, Step-5 starts and accepted attempts, dilute-face restart
coordinates, Step-8 retry/candidate totals, and Step-10 trace-refinement
activations and component indices. Exact repeated Step-5 start
ordinals and Step-8 candidate problems are counted explicitly. A benchmark
may classify an aqueous lower-bound edge case only when dilute-face or trace-
refinement evidence is active and the reported coordinate maps to a neutral
aqueous component; point identity or elapsed time alone is not a valid
classifier.

One exact cache is scoped to a single HELD2 flash and shares deterministic
Provider state evaluations across Steps 2, 3, 5, 8, 9, and the ordinary
Step-10 fallback. The distinct trace-domain evaluator is never served from
that cache, and cache hits are excluded from reported Provider-work totals.
Step 8 also reuses a previously certified feasible or Farkas-infeasible
Problem (67) only when its ordered stable candidate IDs, candidate variables,
coordinate transformation, physical feed, phase-coordinate bounds, immutable
Provider context, and requested neighborhood radius are exactly unchanged.
The cache instance is created after Provider access is assembled and owns the
state, volume-bound, and packing-fraction callbacks for its lifetime. Those
callbacks cannot
be rebound through the Step-8 API. The cache is used by one flash and destroyed
with that flash; a separately constructed cache is required for a different
Provider access. Step 1 and the cache carry the same exact Provider parameter
fingerprint, and Step 8 fails closed before evaluation or replay when they do
not match. Provider failures during phase-bound construction, Problem-(67)
evaluation, or final phase reevaluation are returned as counted indeterminate
evidence. An expanded-radius result is retained
only under that expanded radius. The replay
is reported as `cached_problem_67`. Indeterminate and collapsed candidate
solves are path-dependent and are never retained across history. The immediate
previous result may supply a warm start, but is never itself accepted as a
cached adjudication.

### Canonical HELD2 execution-efficiency contract

The performance route is the scientific route; there is no benchmark-only
controller, chemistry switch, point allowlist, expected-answer branch, or
caller-selectable fast mode. The canonical workflow removes work in this
order:

1. Step 2 searches composition and normalized log-volume jointly,
   independently certifies a strict negative witness, and reuses a successful
   targeted pressure refinement instead of rediscovering the same seed.
2. Step 5 consumes one persistent deterministic start stream. Numerical copies
   do not create new cuts, and a return from Stage III does not replay consumed
   starts.
3. The Step-6 discovery gap is the named five-percent working gate. It moves a
   sufficiently developed hull into the simultaneous Stage-III solve without
   weakening any final KKT, balance, charge, pressure, domain, or phase-identity
   certificate.
4. One flash-scoped cache evaluates each exact ordinary thermodynamic state,
   composition-dependent volume bound, and packing-fraction state once.
   Objective, gradient, Hessian,
   pressure-root, certification, and serialization consumers share the complete
   state rather than calling the Provider through parallel value-only paths.
5. Step 8 replays only certified feasible or independently Farkas-certified
   infeasible Problem-(67) evidence under the complete exact key described
   above. An indeterminate solve, collapsed phase set, changed radius, changed
   candidate order, changed coordinate/feed/bounds payload, or changed Provider
   context always requires fresh adjudication.
6. Step 10 uses the ordinary cache only when it uses the ordinary evaluator.
   The logarithmic charged-trace evaluator is a distinct mathematical domain
   and is never served stale ordinary-state evidence.

Elapsed time is validation evidence, not a solver input. The controller never
branches on a clock, and parallel execution of independent feeds is not counted
as a single-call optimization.

Exact installed-artifact Validation at commit
`9c3323abad6813bdc176d8cca0bc0df1640af2e1` retained all 39 Khudaida
Figures 2--7 solutions: 37 completed within the 30-second target, all 39
completed within the approved 35-second ceiling, and the maximum was 31.487
seconds. The representative Figure-3 point-5 call fell from 49.015 to 23.194
seconds while preserving its accepted solution; Step-8 optimizer solves fell
from 46 to 45 and 32 exact certified Problem-(67) decisions were replayed from
the cache. Provider-evaluation totals across those two artifacts are not
compared because the optimized artifact corrected their accounting. The final
EOS and Equilibrium wheel SHA-256 values are respectively
`3073c100e065252d922481f3fda0a173973f1f7da2d86ec380738285df65e227`
and
`c93d2167397dc27626895ab56f05e12eb553aa6546ddef701bce052fe55cb90b`.
Every final physical and KKT gate passed. Native grand AAD was below the
published ePC-SAFT grand AAD at each Figure-2--7 condition; Figures 3--7 remain
descriptive exploratory parameter extrapolations and create no admission.

The smallest source-bound workflow is Perdomo et al. (2025), Table 3,
NaCl-water at 298.15 K, 2508 Pa, and 5.6 mol NaCl per kg water. In the installed
explicit-species ordering `(water, sodium-cation, chloride-anion)`, its frozen
normalized feed is:

```text
(0.8321050353538131, 0.08394748232309347, 0.08394748232309347)
```

Export the installed Provider model configuration, then run its real controller
with live, flushing terminal output using:

```bash
build/epcsaft-equilibrium-diagnostic \
  --model-config MODEL.json --temperature 298.15 --pressure 2508 \
  --feed 0.8321050353538131,0.08394748232309347,0.08394748232309347 \
  --trace
```

Without `--trace`, the executable is quiet and returns the same structured
result. Current ePC-SAFT/Provider evidence selects the lowest of two stable
reference roots. The Figiel Provider caps total ion mole fraction at `0.38`;
the domain-aware Stage-I map therefore keeps every DIRECT-L composition inside
that immutable source domain. The declared 6500-Provider-callback joint
composition/log-volume allowance completes without detecting a TPD below the
material `-1e-8` margin. Stage II and Stage III are consequently skipped
because no negative witness was detected under the completed finite search.
This is not a stability proof, a reproduction of Perdomo's
SAFT-gamma-Mie endpoint, or an admitted electrolyte-LLE result.

### Homogeneous reacting phase

The D-028 design and deterministic basin-search contract own one homogeneous fixed-`T,P` reacting-phase
foundation. It validates ordered species, conservation and independent
reaction ranks, dimensionless source/reference-bound `lnK`, and the exact
Provider identity before solving. It constructs `g_ref` in the Provider
Helmholtz coordinate basis, uses a general positive electroneutral amount
chart, performs max-min initialization, and evaluates one frozen deterministic
finite start sequence against the same true Provider objective.

The public `chemical_equilibrium` operation exposes that owner without solver
settings, chemistry data, aliases, or a second implementation. Its certificate
axes keep artifact/input completeness, Ipopt status, numerical and physical
checks, reduced-Hessian local status, predictive status, finite search, and
globality separate. The manufactured installed-Provider evidence remains
labeled manufactured/nonpredictive. Source-bound input declares exact ordered
components, a salt-free solvent composition, activity convention, standard
molality, source pressure, and activity-scale shifts. Equilibrium calls the
installed public EOS `source_reference_transfer` operation at the actual
system pressure, validates its ordered charge-neutral basis and immutable
fingerprints, and passes only the resulting Provider-basis constants to the
same compiler and `solve_provider_reaction` owner. The source reference
pressure is retained separately as immutable provenance, and the complete
frozen EOS transfer receipt is attached to the public result.
Failures preserve the returned terminal state and expose each scalar numerical
and physical criterion independently, including the first failed criterion,
active bounds, KKT rank/conditioning, and reduced-Hessian evidence. A
`solve_succeeded` optimizer status is never a chemical certificate by itself.
Second-order admission first reconstructs equality multipliers from the actual
chart objective gradient and constraint Jacobian, rather than transferring a
least-squares multiplier from a differently scaled physical system. It then
uses the full chart Lagrangian Hessian. The physical Lagrangian's covariant
congruence remains separate derivative-check evidence. Diagnostics name every
derivative coordinate and constraint row in storage order and declare their
bases; they also expose the equality multipliers and null-space basis needed to
reproduce the projections. Strictly interior states expose the separate
balance/affinity/pressure KKT root-system Jacobian. When any bound is active,
that matrix is explicitly unavailable rather than presented as a complete
active-set KKT system.
The reported raw eigenvalue spectrum, its convergence status, and raw inertia
are separate from the certified inertia produced by the scale-invariant
diagonal congruence and Cholesky test; neither is silently substituted for the
other.
Its source-complete sentinel is `2 H2O <=> H3O+ + OH-` at 298.15 K and 1 bar
using frozen IAPWS R11-07(2019) reaction data and the immutable Held-2008
Provider catalog. This establishes one local, strictly interior, fixed-`T,P`
value result only.

The operation accepts an optional typed sensitivity request. Without it the
response is explicitly `value_only`; the internal solve still uses its exact
Ipopt and Provider derivatives, but no Jacobian columns are returned. With it
the response is `value_plus_jacobian` when exact derivatives are available and
`value_with_unavailable_jacobian` when they fail closed. Available derivatives
are with respect to compiled balance totals (mol), final Provider-basis
reaction constants (dimensionless `ln(K)`), and pressure (Pa) when the
active-set residual Jacobian is full-rank and acceptably conditioned. Amount
rows follow the result species order and volume is the final scalar state. Each
parameter records input and derivative units.

The typed sensitivity receipt reports rank, condition number, active bounds,
amount-chart topology, and the immutable Provider parameter fingerprint. The
result also binds installed Equilibrium and Provider distribution versions and
RECORD fingerprints plus the Provider capsule name, ABI version, table size,
and result-structure sizes. Provider-basis inputs retain exact pressure and
active-parameter sensitivity support. Source-bound sensitivity requests fail
before optimization while the installed source-reference transfer advertises
no derivatives. Redundant source
reaction rows are reduced through the compiler's exact reaction transformation
before entering the KKT column. Provider active-parameter coordinates are
available only when the installed capability table advertises each requested
typed coordinate and one atomic Provider callback supplies the active-model
Helmholtz state derivatives, state-parameter block, pressure and
chemical-potential projections, and packing state
gradient/Hessian. Request order is preserved. Missing blocks, unadvertised
coordinates, active bounds, or unavailable KKT columns fail closed; no missing
column is filled numerically. There is no
chemistry registry, predictive admission, coupled-equilibrium result, or
globality claim.

A source-bound active request requires source-reference transfer derivatives
for the requested coordinates. If that prerequisite is unavailable, the
operation rejects before solving rather than combining an active phase with a
fixed-catalog reference value.

## Shared package contract

The installed `epcsaft` Provider owns resolved thermodynamic input, component
identity and order, charges, model fingerprints, Helmholtz phase values,
packing fractions, validity domains, and nonlinear derivative tensors.
Equilibrium owns:

- fixed-`T,P` equilibrium formulations and coordinate transformations;
- phase stability and discovery controllers;
- deterministic root and finite global-exploration accounting;
- stage-specific LP and NLP construction, including exact Lagrangian
  derivatives where required;
- local numerical and physical certification;
- fail-closed controller outcomes and typed results; and
- package-authored, public-artifact Validation campaigns assigned by Migration.

Equilibrium does not copy Provider equations or parameters, link Provider
implementation symbols, reach through private Provider modules, own density
closure, or create a selectable numerical derivative backend. Exact algebraic
and implicit chain rules assembled from Provider tensors remain Equilibrium
work, not a second EOS derivative owner.

## Evidence and claim axes

Every mixture result keeps the following questions independent:

1. **Artifact and input completeness:** Is the installed Provider table,
   fingerprint, component order, temperature, feed, and model domain suitable?
2. **Solver status:** What did the assigned root search, global explorer, LP
   solver, local NLP solver, or bounded controller terminate with?
3. **Numerical status:** Do original-coordinate feasibility, stationarity,
   complementarity, derivative, and search-accounting checks pass?
4. **Physical status:** Do balance, normalization, charge, pressure, potential,
   stability, phase-identity, distinctness, and source-domain checks support the
   interpretation?
5. **Search completeness:** Did every declared finite start or major iteration
   needed for the terminal claim complete?
6. **Root completeness:** Which homogeneous pressure roots were detected, and
   was completeness established? Archived installed HELD2 evidence reports
   `root_completeness="not_proven"`.
7. **Predictive status:** Is there an accepted physical output that can be
   compared with a source-bound external case?
8. **Globality:** Did the method establish global phase stability? Finite HELD
   and HELD2 searches always report
   `globality_certificate="not_guaranteed"`.

No axis implies a later axis. A finite, independently certified negative TPD
point is an existence certificate for homogeneous instability even when
unrelated search attempts are incomplete. The converse is asymmetric:
`no_negative_found`, no-improvement, and finite-gap claims require complete
declared-search accounting and still do not prove globality.

Missing or unavailable evidence remains `not_adjudicated`; it must not be
converted into a fake numerical value, failure, or success.

## HELD2 numerical tolerance contract

HELD2 uses the named contract in `cpp/src/held2_tolerances.hpp`. A tolerance
belongs to one mathematical quantity and one evidence owner; it is not a
general-purpose accuracy knob. When a material scale exists, the gate is

```text
abs(residual) <= atol + rtol * scale
```

and diagnostics retain the raw residual, scale, both tolerances, the evaluated
threshold, predicate result, category, and failure meaning. Exact-boundary
semantics are defined by each named relation. Optimizer targets are requests
to the solver and never substitute for scientific acceptance.

| Owner | Named quantities | Contract |
|---|---|---|
| Representation | chart contact, composition sum, scaled charge, reconstructed ion | `1e-9` absolute |
| Representation | bound activity | `1e-8` absolute |
| Pressure roots | relative pressure, log-volume width, stationary residual, Provider-boundary and fallback duplicate classification | `1e-8`, `1e-9`, `1e-9`, `1e-8`, `1e-8` respectively |
| Root topology | mechanical margin; stable-objective tie | strict margin `>1e-6`; tie `1e-8 + 1e-9*scale` |
| Stage I | materially negative reduced TPD | `< -1e-8`; `[-1e-8,1e-8]` is inconclusive |
| Stage-II LP | primal, dual, complementarity, active-cut diagnostic | `1e-9 + 1e-8*scale`, `1e-9 + 1e-8*scale`, `1e-8`, `1e-7` |
| Stage-II discovery KKT | primal, dual sign, pullback, stationarity, complementarity | `1e-8`, `1e-9`, `1e-9 + 1e-9*scale`, `1e-6`, `1e-8` |
| Stage-III feasibility LP | Farkas row sign, column dual feasibility, strict contradiction margin | `1e-10`, `1e-10`, `>1e-9 + 1e-9*scale` |
| Step 6 | gap; fixed-volume gradient | `1e-8`; `1e-8 + 1e-7*scale` |
| Generic candidate identity | numerical duplicate; confidently distinct | `<=1e-7`; either physical composition or log-volume distance `>1e-5` |

When Stage III returns to Stage II because an otherwise feasible phase set
fails the non-trace modified-potential check, the failed phases remain as
search feedback and Step 5 must discover a distinct basin before Stage III is
retried. Reusing the same phase set cannot repair chemical-potential equality.
| HELD2 persistent \(\mathcal M\) identity | modified composition and relative molar-volume numerical copy | `<=1e-8` in every modified-composition coordinate and log-volume difference |
| HELD2 Step-6 discovery | Perdomo Eq. (66) upper/lower working gap; fixed-volume derivative agreement | `5e-2`; `0.5*scale + 1e-8` |
| HELD2 Step-6 distinctness | Perdomo Eq. (66) Provider packing fraction or modified composition | either difference `>=1e-3` |
| Stage III | modified/explicit balances, scaled charge, pressure | `1e-8`, `1e-8`, `1e-9`, `1e-8` |
| Stage III | modified potentials; KKT, dual sign, complementarity, free-energy gap | `1e-8 + 1e-7*scale`, `1e-7`, `1e-9`, `1e-8`, `1e-8` |
| Stage III paper convergence | Perdomo Eqs. (68)-(69) free-energy gap; modified-potential ratio | `1e-4`; `1e-2` |
| Phase identity | active, retirement evidence, numerical merge, confidently distinct | `>1e-8`, `>1e-8`, `<=1e-5`, `>1e-4` |
| Three-or-more-phase coalescence | physical composition and log molar-volume distance | both `<=1e-2`, followed by a reduced Stage-III solve and full recertification |
| Step 10 trace | charged physical mole-fraction interval; potential residual | `[1e-300, 5e-10]`; `1e-8` absolute |
| Ipopt | default target; HELD2 Stage-III target; disabled acceptable target; constraint target | `1e-10`; `1e-8`; `1e-9`; `1e-10`; zero bound relaxation |

Candidate distances between `1e-7` and `1e-5`, and phase distances between
`1e-5` and `1e-4`, are unresolved identity bands and fail closed. Marginal
roots, unresolved stable-root ordering, and unavailable evidence likewise do
not become acceptance. Direct invalid user input remains invalid; the
representation allowances apply only to validated transformations and solver
iterates. The Step-10 trace domain relaxes only a charged independent
coordinate's finite-search lower bound; it does not relax any other polytope
constraint. Finite search always retains
`globality_certificate="not_guaranteed"`.

The paper-scale coalescence neighborhood is only considered when Stage III
contains at least three phases. It does not merge the returned states in
place: one member is removed and the smaller problem must independently pass
the same balance, charge, pressure, KKT, and paper-convergence gates. A
two-phase result therefore cannot be collapsed by this recovery rule.

HELD2 keeps three identities separate. Appendix C and persistent
\(\mathcal M\) membership use modified composition plus molar volume. Step 6
alone uses Provider packing fraction plus modified composition for Eq. (66).
Step 8 merges numerical phase copies using physical composition plus molar
volume; log-volume difference is only the dimensionless representation of
relative volume in the two volume-based identities.

## Future formulation boundaries

### Closed Ascani counterion-pair equilibrium

Ascani's independent counterion-pair residual formulation is a separate
scientific family with different coordinates and stationarity conditions. The
lab retains source and historical implementation evidence, but this package
has no current Ascani route, result family, or admission. Future work requires
its own package-local design and bounded capability gate.

### Reactive and coupled equilibrium

The homogeneous value operation does not complete simultaneous
phase-chemical equilibrium. A coupled formulation still requires explicit
phase incidence, global conservation, phase-specific electroneutrality,
source-complete standard-state transformations, phase discovery, and distinct
reaction/transfer/pressure certificates. Staged chemistry followed by a
phase-only solve is initialization evidence, not a coupled equilibrium result.
The package-local
[GREPE reactive phase-equilibrium design](designs/2026-07-28-grepe-reactive-phase-equilibrium.md)
preserves those requirements and revises phase discovery around the validated
HELD2 joint composition/log-volume search, immutable candidate identity,
master-to-NLP initialization, active-set re-solves, and fail-closed evidence
rules. Its first physical subject is the neutral Ascani--Senina esterification
case selected in the
[benchmark review](designs/2026-07-29-reactive-lle-benchmark-review.md);
electrolyte reactive LLE is a separately evidenced later slice. The GREPE file
is the normative future implementation contract, not a runtime or capability
claim. No public coupled reactive-phase schema or runtime route is admitted.

## Historical and scientific provenance

The following sources remain important but are not current package runtime
authority:

- `tannerpolley/ePC-SAFT-lab:docs/latex/equilibrium_formulations.tex` is the
  broad mathematical formulation record. Its old public-route and
  implementation-status text predates this package.
- `tannerpolley/ePC-SAFT-lab:docs/superpowers/milestones/M4-equilibrium/`
  preserves the generalized architecture, Pereira, Perdomo, Ascani, reactive
  doctrine, admission registry, preservation manifest, and dated dashboard.
  The directory explicitly identifies itself as historical archive evidence.
- Pereira et al. (2012) is the primary source for neutral HELD.
- Perdomo et al. (2025) is the primary source for modified-mole HELD2.
- `ePC-SAFT-migration/MIGRATION.md` records D-024 through D-026 sequencing,
  exact checkpoints, artifacts, and review outcomes.
- `ePC-SAFT-validation` owns source ledgers and durable installed-artifact
  campaign evidence. Package-authored evidence cannot accept its own
  promotion.

Historical documents retain their original statements for provenance. Their
old route names, package paths, runtime-status fields, and planning-family
labels must not be used to infer the current public surface or capability
state.
