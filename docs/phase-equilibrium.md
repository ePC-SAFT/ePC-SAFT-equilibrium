# Phase-Equilibrium Documentation Authority

Status: canonical package documentation

Authority effect: none. This document describes ownership and current evidence;
it does not admit, promote, or publish a capability.

## Authority and capability status

Organization doctrine revision 3 defines the ecosystem authority hierarchy.
Within this repository, `AGENTS.md` and `CONTEXT.md` govern package scope,
`ARCHITECTURE.yaml` records the machine-readable architecture, and the design
documents named below own the scientific and numerical contracts of individual
formulations. Plans and receipts record execution and evidence; they do not
replace a formulation owner.

| Formulation | Package design owner | Current capability state |
| --- | --- | --- |
| Pure-component saturation | [Pure-saturation slice](designs/2026-07-17-pure-saturation-slice.md) | Accepted only for the exact methane, ethane, and propane scope in `promotion-0018-equilibrium-pure-saturation-v1` |
| Neutral Pereira HELD | [Neutral HELD v1](designs/2026-07-17-neutral-held-v1.md) | Frozen local candidate; installed campaign retained as `NON_ADMISSION`; controller redesign deferred |
| Strong-electrolyte Perdomo HELD2 | [Perdomo HELD2](designs/2026-07-21-perdomo-held2.md) and [solver-strategy plan](plans/2026-07-21-perdomo-held2-solver-strategy.md) | Private integrated Stage-I/II/III development controller with fail-closed installed evidence; no public electrolyte dispatch or admitted electrolyte LLE capability |
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
- `tp_flash` is the sole mixture flash surface. Current `main` dispatches the
  reviewed neutral binary fingerprint to the neutral HELD controller. The
  archived non-production HELD2 subject also dispatched qualifying installed
  strong-electrolyte Provider capability tables, but that experimental runtime
  is not part of current `main`.

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
The linked implementation plan assigns deterministic pressure-root
enumeration to density topology, DIRECT-L to the reduced Stage-I search, HiGHS
to the Stage-II upper LP, basin exploration plus exact-Hessian Ipopt to the
Stage-II lower problem, and exact-Hessian Ipopt to Stage III. The private
integrated controller executes that stage order and retains every bounded
failure as evidence. Replaced HELD2 runtime routes and baseline fixtures are
removed; only focused manufactured formulation oracles remain. None of these
internal owners is a caller-selectable backend or current public electrolyte
behavior.

The remaining installed-completion contract is specified in
[HELD2 Installed Completion](designs/2026-07-22-held2-installed-completion.md)
and sequenced by the
[HELD2 Installed Completion Plan](plans/2026-07-22-held2-installed-completion.md).
Those documents begin from the retained Stage-II-indeterminate artifact and do
not reinterpret it as an equilibrium result. Current source has classified the
artifact's unit-cube callback failure as bounded binary64 contact, preserves
ordered trial/accepted-iterate and upper-LP evidence, and reaches successful
local Ipopt termination for that exact start. It still has no certified
installed Stage-II candidate set or physical Stage-III result.

### HELD2 live progress diagnostic

The private installed controller accepts an internal, nullable progress
observer. It receives already-computed reference roots, DIRECT-L evaluations,
HiGHS bounds, Ipopt iteration metrics, certificate results, stage skips, and
the final controller status. The default path supplies no observer and remains
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
reference roots, then fails closed at the first Stage-I trial composition
because the Provider cannot construct that trial's molar-volume domain. Stage
II and Stage III are reported as skipped with the exact causal reason. This is
a workflow diagnostic, not a reproduction of Perdomo's SAFT-gamma-Mie endpoint
and not an admitted electrolyte-LLE result.

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
| Ipopt | target, disabled acceptable target, constraint target | `1e-10`, `1e-9`, `1e-10`; zero bound relaxation |

Candidate distances between `1e-7` and `1e-5`, and phase distances between
`1e-6` and `1e-4`, are unresolved identity bands and fail closed. Marginal
roots, unresolved stable-root ordering, and unavailable evidence likewise do
not become acceptance. Direct invalid user input remains invalid; the
representation allowances apply only to validated transformations and solver
iterates. Finite search always retains
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
