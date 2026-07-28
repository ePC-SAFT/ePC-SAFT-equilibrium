# Phase-Equilibrium Documentation Authority

Status: canonical package documentation

Authority effect: none. This document describes ownership and current evidence;
it does not admit, promote, or publish a capability.

## Authority and capability status

Organization doctrine revision 4 defines the ecosystem authority hierarchy.
Within this repository, `AGENTS.md` and `CONTEXT.md` govern package scope,
`ARCHITECTURE.yaml` records the machine-readable architecture, and the design
documents named below own the scientific and numerical contracts of individual
formulations. Plans and receipts record execution and evidence; they do not
replace a formulation owner.

| Formulation | Package design owner | Current capability state |
| --- | --- | --- |
| Pure-component saturation | [Pure-saturation slice](designs/2026-07-17-pure-saturation-slice.md) | Accepted only for the exact methane, ethane, and propane scope in `promotion-0018-equilibrium-pure-saturation-v1` |
| Neutral Pereira HELD | [Neutral HELD v1](designs/2026-07-17-neutral-held-v1.md) | Frozen local candidate; installed campaign retained as `NON_ADMISSION`; controller redesign deferred |
| Strong-electrolyte Perdomo HELD2 | [Paper-faithful Steps 1--10](designs/2026-07-24-held2-paper-algorithm.md), with [earlier design provenance](designs/2026-07-21-perdomo-held2.md) | Public development dispatch over one native core with fail-closed installed evidence; no admitted electrolyte LLE capability |
| Superseded fixed two-phase route | [Historical fixed-route design](designs/2026-07-17-neutral-two-phase-tp-flash.md) | Removed without alias; retained only as provenance |
| Ascani counterion-pair electrolyte equilibrium | No current runtime design | Closed future formulation; historical lab evidence only |
| Reactive or coupled phase-chemical equilibrium | No current runtime design | Closed future formulation; no public schema or runtime route |

The only accepted capability is `pure-component-saturation-v1`. A public
symbol, a local candidate, an installed campaign, or a converged local solve is
not an authority receipt.

## Public operations

The package exports two equilibrium operations:

- `saturation` owns the accepted, bounded, pure-component local saturation
  boundary.
- `tp_flash` is the sole mixture flash surface. It dispatches the reviewed
  neutral binary fingerprint to neutral HELD and qualifying installed strong-
  electrolyte Provider capability tables to the non-admitted HELD2
  development route.

`tp_flash` does not accept a phase count, caller seeds, solver settings,
backend selection, or a case-specific mode. Its public result owners are
`TpFlashResult`, `HeldDiagnostics`, and `FlashError`. The removed
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
enumeration to density topology, DIRECT-L to the reduced Stage-I search, HiGHS
to the Stage-II upper LP, deterministic capped-multistart exact-Hessian Ipopt
to Step 5, and exact-Hessian Ipopt to Problem (67) in Step 8. One typed
`run_held2_algorithm` state machine owns the closed Steps 1--10 transitions for
both installed and manufactured problems. Stage II has one major loop and one
Eq. (66) decision owner. Each major retains its upper-solve identity, `UBD`,
multipliers, active cuts, local attempts, pressure branches, and certificate
partitions. Replaced HELD2 runtime routes and baseline fixtures are deleted;
focused manufactured numerical oracles call the same transition and step
owners. None of these internal owners is a caller-selectable backend; public
electrolyte dispatch remains non-admitted development behavior.

The remaining installed-completion contract is specified in
[HELD2 Installed Completion](designs/2026-07-22-held2-installed-completion.md)
and sequenced by the
[HELD2 Installed Completion Plan](plans/2026-07-22-held2-installed-completion.md).
Those documents begin from the retained Stage-II-indeterminate artifact and do
not reinterpret it as an equilibrium result. Current source has classified the
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

The smallest source-bound workflow is Perdomo et al. (2025), Table 3,
NaCl-water at 298.15 K, 2508 Pa, and 5.6 mol NaCl per kg water. In the installed
explicit-species ordering `(water, sodium-cation, chloride-anion)`, its frozen
normalized feed is:

```text
(0.8321050353538131, 0.08394748232309347, 0.08394748232309347)
```

Run its real controller with live, flushing terminal output using:

```bash
pytest -s tests/test_perdomo_held2_trace.py::test_perdomo_table3_nacl_workflow --held2-live
```

Without `--held2-live`, the same test is quiet and asserts the same structured
result. Current ePC-SAFT/Provider evidence selects the lowest of two stable
reference roots. The Figiel Provider caps total ion mole fraction at `0.38`;
the domain-aware Stage-I map therefore changes the first DIRECT-L midpoint from
the inadmissible `0.50` physical ion fraction to `0.1900000001`. All 50 declared
evaluations complete without Provider failure and detect no TPD below the
material `-1e-8` margin; the smallest detected TPD is positive and near zero.
Stage II and Stage III are consequently skipped because no negative witness was
detected under the completed finite search. This is not a stability proof, a
reproduction of Perdomo's SAFT-gamma-Mie endpoint, or an admitted electrolyte-
LLE result.

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
| Stage-II KKT | primal, dual sign, pullback, stationarity, complementarity | `1e-8`, `1e-9`, `1e-9 + 1e-9*scale`, `1e-7`, `1e-8` |
| Step 6 | gap; fixed-volume gradient | `1e-8`; `1e-8 + 1e-7*scale` |
| Candidate identity | numerical duplicate; confidently distinct | `<=1e-7`; either physical composition or log-volume distance `>1e-5` |
| Stage III | modified/explicit balances, scaled charge, pressure | `1e-8`, `1e-8`, `1e-9`, `1e-8` |
| Stage III | modified potentials; KKT, dual sign, complementarity, free-energy gap | `1e-8 + 1e-7*scale`, `1e-7`, `1e-9`, `1e-8`, `1e-8` |
| Phase identity | active, retirement evidence, merge, confidently distinct | `>1e-8`, `>1e-8`, `<=1e-6`, `>1e-4` |
| Step 10 trace | charged physical mole-fraction interval; potential residual | `[1e-300, 5e-10]`; `1e-8` absolute |
| Ipopt | target, disabled acceptable target, constraint target | `1e-10`, `1e-9`, `1e-10`; zero bound relaxation |

Candidate distances between `1e-7` and `1e-5`, and phase distances between
`1e-6` and `1e-4`, are unresolved identity bands and fail closed. Marginal
roots, unresolved stable-root ordering, and unavailable evidence likewise do
not become acceptance. Direct invalid user input remains invalid; the
representation allowances apply only to validated transformations and solver
iterates. The Step-10 trace domain relaxes only a charged independent
coordinate's finite-search lower bound; it does not relax any other polytope
constraint. Finite search always retains
`globality_certificate="not_guaranteed"`.

## Closed future formulations

### Ascani counterion-pair equilibrium

Ascani's independent counterion-pair residual formulation is a separate
scientific family with different coordinates and stationarity conditions. The
lab retains source and historical implementation evidence, but this package
has no current Ascani route, result family, or admission. Future work requires
its own package-local design and bounded capability gate.

### Reactive and coupled equilibrium

Standalone chemical equilibrium and simultaneous phase-chemical equilibrium
require ordered species, phase incidence, stoichiometric and elemental balance
matrices, standard states, source-bound equilibrium constants, and distinct
reaction/transfer certificates. Staged chemistry followed by a phase-only
solve is initialization evidence, not a coupled equilibrium result. No such
public schema or runtime route is currently admitted.

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
