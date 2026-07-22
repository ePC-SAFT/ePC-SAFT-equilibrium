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
| Strong-electrolyte Perdomo HELD2 | [Perdomo HELD2](designs/2026-07-21-perdomo-held2.md) | Non-production development candidate; public installed dispatch exists, but electrolyte LLE is not admitted |
| Private homogeneous reacting phase | [Private reacting-phase kernel](designs/2026-07-21-private-reacting-phase-kernel.md) | D-028-bound non-production foundation; underscored test seam only; manufactured evidence plus a generic neutral-reference source transform; no admitted source-complete chemistry |
| Superseded fixed two-phase route | [Historical fixed-route design](designs/2026-07-17-neutral-two-phase-tp-flash.md) | Removed without alias; retained only as provenance |
| Ascani counterion-pair electrolyte equilibrium | No current runtime design | Closed future formulation; historical lab evidence only |
| Coupled multiphase chemical equilibrium | No current runtime design | Closed future formulation; no public schema or runtime route |

The only accepted capability is `pure-component-saturation-v1`. A public
symbol, a local candidate, an installed campaign, or a converged local solve is
not an authority receipt.

## Public operations

The package exports two equilibrium operations:

- `saturation` owns the accepted, bounded, pure-component local saturation
  boundary.
- `tp_flash` is the sole mixture flash surface. It dispatches the reviewed
  neutral binary fingerprint to the neutral HELD controller and qualifying
  installed strong-electrolyte Provider capability tables to the Perdomo HELD2
  controller. Unsupported multicomponent non-electrolyte models fail before
  search.

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

### Private homogeneous reacting phase

The D-028 design owns one private homogeneous fixed-`T,P` reacting-phase
foundation. It validates ordered species, conservation and independent
reaction ranks, dimensionless source/reference-bound `lnK`, and the exact
Provider identity before solving. It constructs `g_ref` in the Provider
Helmholtz coordinate basis, uses a general positive electroneutral amount
chart, performs max-min initialization, attempts the true Provider objective
first, and permits only a final `lambda=1` continuation state to pass.

Its certificate axes keep artifact/input completeness, Ipopt status, numerical
and physical checks, reduced-Hessian local status, predictive status, finite
search, and globality separate. The manufactured installed-Provider seam
remains labeled manufactured/nonpredictive. The current private reference seam
consumes the installed Provider's complete charge-neutral basis and transforms
source reaction records with explicit activity scales, exact source/reference
identities, representation residuals, and Provider derivative availability.
The former water-only source-complete profile is unavailable under the current
Provider artifact and must not be inferred from the retained IAPWS source
equation record. There is no public solve route, result type, generalized
chemistry registry, or globality claim.

An authenticated Provider bundle with explicitly unknown applicability may reach
the underscored phase-block evidence seam for diagnostic transport checks. It
cannot enter an Ipopt solve or pass physical validity until Provider supplies a
finite source-temperature domain.

## Shared package contract

The installed `epcsaft` Provider owns resolved thermodynamic input, component
identity and order, charges, model fingerprints, Helmholtz phase values,
packing fractions, validity domains, and nonlinear derivative tensors.
Equilibrium owns:

- fixed-`T,P` equilibrium formulations and coordinate transformations;
- phase stability and discovery controllers;
- Ipopt problem construction and exact Lagrangian derivatives;
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
2. **Solver status:** What did Ipopt or the bounded controller terminate with?
3. **Numerical status:** Do original-coordinate feasibility, stationarity,
   complementarity, derivative, and search-accounting checks pass?
4. **Physical status:** Do balance, normalization, charge, pressure, potential,
   stability, phase-identity, distinctness, and source-domain checks support the
   interpretation?
5. **Search completeness:** Did every declared finite start or major iteration
   needed for the terminal claim complete?
6. **Root completeness:** Which homogeneous pressure roots were detected, and
   was completeness established? Current HELD2 evidence reports
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

## Closed future formulations

### Ascani counterion-pair equilibrium

Ascani's independent counterion-pair residual formulation is a separate
scientific family with different coordinates and stationarity conditions. The
lab retains source and historical implementation evidence, but this package
has no current Ascani route, result family, or admission. Future work requires
its own package-local design and bounded capability gate.

### Reactive and coupled equilibrium

The private homogeneous foundation does not complete simultaneous
phase-chemical equilibrium. A coupled formulation still requires explicit
phase incidence, global conservation, phase-specific electroneutrality,
source-complete standard-state transformations, phase discovery, and distinct
reaction/transfer/pressure certificates. Staged chemistry followed by a
phase-only solve is initialization evidence, not a coupled equilibrium result.
No public reactive schema or runtime route is admitted.

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
