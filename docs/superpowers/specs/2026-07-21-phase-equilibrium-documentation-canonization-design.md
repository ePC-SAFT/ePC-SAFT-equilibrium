# Phase-Equilibrium Documentation Canonization Design

**Status:** approved design for documentation-only implementation

**Authority effect:** none; this design reorganizes and corrects package documentation but does not admit a capability, promote an artifact, or change runtime behavior

## Goal

Make `ePC-SAFT-equilibrium` the clear package-local documentation owner for its
implemented phase-equilibrium formulations while preserving the permanent lab,
Migration, Validation, and primary papers in their correct provenance and
evidence roles.

The result must let a reader answer four questions without consulting a
chronological migration diary:

1. Which document governs each implemented equilibrium formulation?
2. Which capabilities are accepted, non-production candidates, or closed?
3. Which thermodynamic and numerical contracts does the Perdomo HELD2
   implementation follow?
4. Which claims remain unproven, including search completeness, globality, and
   electrolyte two-liquid predictive acceptance?

## Governing authority

The documentation follows Organization doctrine revision 3 in
`../ePC-SAFT-organization/GOVERNANCE.md`:

1. Organization doctrine governs repository ownership and admission rules.
2. Current Equilibrium `AGENTS.md` and `CONTEXT.md` govern package scope.
3. `ARCHITECTURE.yaml` and package-local design documents govern implemented
   architecture and formulation contracts.
4. Approved plans govern bounded implementation work.
5. Receipts and Validation results provide exact evidence and authority state.

Primary papers govern claims about Pereira HELD, Perdomo HELD2, and Ascani
electrolyte equilibrium. They do not govern this package's API, Provider ABI,
diagnostics vocabulary, or capability admission.

The permanent lab M4 documents are broad mathematical and historical
provenance. Their old monorepo routes, activation paths, issue queue, and
runtime-status fields must not be copied into current package authority.

## Chosen structure

### `docs/phase-equilibrium.md`

Create one compact package-level authority map. It owns:

- the documentation hierarchy and one-owner rule;
- the accepted versus candidate versus closed capability table;
- the current public operation map (`saturation` and `tp_flash`);
- links to the pure-saturation, neutral HELD, and Perdomo HELD2 designs;
- the shared Provider/Equilibrium ownership boundary;
- the common solver, numerical, physical, search-completeness,
  root-completeness, predictive, and globality claim boundaries;
- a bounded future-formulation section for Ascani counterion-pair equilibrium
  and reactive/coupled equilibrium, both explicitly closed and noncanonical for
  current runtime behavior; and
- provenance links to the M4 archive, mathematical formulation, Migration,
  Validation, and primary papers.

This file must not duplicate the detailed mathematics of the three implemented
slice designs.

### `docs/designs/2026-07-21-perdomo-held2.md`

Create the missing Equilibrium-owned scientific and numerical design for the
existing non-production Perdomo HELD2 controller. It owns the current package
contract described below.

#### Scope and nonclaims

- Fixed-temperature, fixed-pressure, nonreactive fluid equilibrium for strong,
  fully dissociated electrolyte mixtures supported by an installed Provider
  capability table.
- At least one molecular species and at least two charged species.
- No weak-electrolyte speciation, chemical reactions, solids, hydrates,
  reactive equilibrium, caller-supplied phase count, case-specific route, or
  deterministic global optimization claim.
- The public operation remains `tp_flash`; result and diagnostic owners remain
  `TpFlashResult`, `HeldDiagnostics`, and `FlashError`.
- The implementation remains non-production and unpromoted. The accepted
  capability remains pure-component saturation only.

#### Scientific sources

- Perdomo et al. 2025 is the primary algorithm source.
- Pereira et al. 2012 is the source for the parent three-stage HELD structure.
- The lab `equilibrium_formulations.tex` and M4 algorithm doctrines are
  mathematical provenance, not current package authority.
- Migration D-024 through D-026 and Validation ledgers are implementation and
  evidence provenance, not equation owners.

#### Modified-mole coordinates

For charged species `E`, chosen from the charged species with largest absolute
charge, phase electroneutrality eliminates its amount:

```text
n_E = -(1 / z_E) * sum_{i != E}(z_i n_i)
```

For retained species `i`:

```text
n_bar_i = (1 - z_i / z_E) n_i
x_bar_i = n_bar_i / sum_j(n_bar_j)
```

One retained molecular species is dependent under modified-fraction
normalization. The optimization therefore has `C - 2` independent modified
composition coordinates. Lift and back-lift must preserve ordinary material
balance, normalization, nonnegative physical amounts, and per-phase charge.
Equal-charge singular charts fail closed.

The interphase thermodynamic quantities are the modified potentials

```text
mu_bar_i = (mu_tilde_i - (z_i / z_E) mu_tilde_E) / (1 - z_i / z_E)
```

not unconstrained individual-ion chemical potentials. Individual-ion
comparisons require a separately declared Galvani-potential convention.

#### Provider boundary

Equilibrium consumes only the installed Provider SDK capability table and its
phase, packing, source-domain, component-order, charge, fingerprint, gradient,
and Hessian payloads. It does not copy EOS or packing equations.

The Provider coordinates are component amounts followed by volume. Equilibrium
owns exact algebraic and implicit transformations into modified composition,
`log(V)`, `q = log(eta)`, Stage-II feasible-simplex, and Stage-III extensive
phase coordinates. Production finite differences are forbidden.

Temperature, component order, total-ion limits, packing limits, and
composition-specific volume bounds are Provider-owned input-domain evidence.
Every Provider call must remain inside the declared domain. Raw bounds remain
diagnostic evidence; solver bounds use deterministic representable interiors
when strict endpoint exclusion is required.

#### Homogeneous reference and Stage I

At fixed feed, the reference search scans the exact Provider log-volume domain,
brackets and refines pressure roots using Provider phase value, gradient, and
Hessian data, classifies mechanical stability, and selects the mechanically
stable root with the lowest `A/(RT) + P V/(RT)`.

The search also examines pressure stationary points so tangential roots cannot
be silently classified as absent. Marginal roots, domain-bound roots,
incomplete evaluation/refinement, and objective-tied lowest references fail
closed. Detected-root accounting is reported with
`root_completeness="not_proven"`; it is independent of globality.

Stage I evaluates the modified-composition/volume tangent-plane distance over
the declared deterministic finite multistart. A feasible, materially distinct,
confirmed negative TPD point is an existence certificate for homogeneous
instability even if unrelated starts remain incomplete. Conversely,
`no_negative_found` requires complete declared-start accounting and does not
prove global stability.

#### Stage II

Stage II retains a complete cut pool and alternates:

- the enumerated linear upper problem in the modified multipliers; and
- nonconvex lower problems in modified composition and packing/volume.

Changing-composition lower solves use the private feasible-simplex chart and
`q = log(eta)`. Exact chain rules transform Provider gradients and Hessians.
Physical acceptance is recomputed in original modified-composition and packing
coordinates; internal chart convergence is insufficient.

Cube duals are pulled through the exact chart Jacobian and decomposed against
the original individual-bound and modified-simplex inequality gradients.
Certification requires correct multiplier signs, original-coordinate
stationarity, primal feasibility, complementarity, and negligible
reconstruction residual. Solver termination remains independent.

A certified improving lower point may add a cut even when unrelated starts are
incomplete. Incomplete search must remain visible and cannot be labeled a
completed lower bound. Declaring no improving cut or a final finite-search gap
requires complete accounting.

Step 6 runs immediately after every certified improving-cut insertion and on
the no-improvement path. Candidate eligibility uses fixed-volume composition
gradients, the source `I^m` active-component selection, relative multiplier
scaling, pressure/packing stationarity, and source-faithful eta-or-composition
phase distinctness. The cut-growth duplicate test remains a separate concept.

#### Stage III

For candidate phases `alpha`, phase fractions `phi_alpha`, molar volumes
`v_alpha`, and modified compositions `x_bar_alpha`, the direct reduced
total-free-energy problem is

```text
min sum_alpha phi_alpha * [A_bar(T, v_alpha, x_bar_alpha) + P0 v_alpha]
```

subject to modified material balances, phase-fraction normalization,
nonnegative phase fractions, modified-composition domains, packing/volume
domains, and candidate-neighborhood restrictions.

Stage III is one general `m_p` owner. Its TNLP includes algebraic per-phase
simplex inequalities before Provider evaluation, exact constant Jacobians,
zero simplex-constraint Hessian terms, and simplex multipliers in the final
Lagrangian-gradient certificate. Recoverable invalid Ipopt trials do not
permanently poison terminal diagnostics; a declared-valid Provider failure
remains observable.

Postsolve certification includes ordinary and modified material balance,
per-phase charge, phase fractions, pressure, all retained modified potentials,
original-coordinate KKT and complementarity, packing/source domains, phase
distinctness, duplicate merging, and the Stage-II/Stage-III feedback path.
Trace-bound species use complementarity and logarithmic refinement rather than
unconditional interior potential equality.

No documentation may claim runtime KKT-certified inactive-phase retirement
unless that behavior is present and approved in the current implementation.
A numerical phase-amount cutoff alone is not an equilibrium certificate.

#### Status and claim axes

The design documents the existing result owner and preserves these independent
questions:

- artifact and input completeness;
- solver termination;
- numerical certification;
- physical certification;
- search completeness/resource status;
- reference-root completeness;
- predictive comparison when an accepted physical result exists; and
- globality.

`globality_certificate="not_guaranteed"` is invariant for finite HELD/HELD2
search. No other status implies globality. Missing evidence uses the existing
not-adjudicated vocabulary rather than a fake failure or numeric default.

#### Current evidence boundary

- The public installed Perdomo Table 3 tracer reaches an accepted homogeneous
  one-phase result with three detected pressure roots, two mechanically stable
  roots, complete Stage-I start accounting, and
  `root_completeness="not_proven"`.
- That result is a cross-EOS mathematical extrapolation and source-topology
  disagreement, not an accepted electrolyte LLE or Perdomo numerical
  reproduction.
- D-026 selects one source-complete installed ePC-SAFT two-liquid case as the
  next gate. Until it passes public Stage I, Stage II, Stage III, and every
  physical certificate, electrolyte LLE remains unadmitted.

Exact commits, artifact hashes, case payloads, and campaign results stay in
Migration and Validation rather than becoming timeless formulation constants.

## Existing documents retained without scientific rewrite

- `docs/designs/2026-07-17-pure-saturation-slice.md` remains the design owner
  for the accepted saturation capability.
- `docs/designs/2026-07-17-neutral-held-v1.md` and its plan remain frozen
  neutral-HELD candidate and non-admission provenance.
- `docs/designs/2026-07-17-neutral-two-phase-tp-flash.md` remains explicitly
  superseded history.
- Candidate and promotion receipts remain immutable evidence.

## Governing-file corrections

### `CONTEXT.md`

- Change `governance_doctrine_revision` from 2 to 3.
- Describe the package as implementation owner while keeping accepted
  capability authority receipt-bound.
- Add the Perdomo source/design identity and link the canonical index.
- Record root-completeness independence and the D-026 waiting boundary without
  claiming electrolyte LLE admission.

### `ARCHITECTURE.yaml`

- Change `governance_doctrine_revision` from 2 to 3.
- Clarify implementation-source versus accepted-authority semantics.
- Link the canonical index and Perdomo design.
- Add the current Provider/domain/root/certificate contract at architecture
  level without duplicating the detailed design.
- Record D-026 as the next evidence gate and preserve authority effect `none`.

### `README.md`

- Add the canonical phase-equilibrium index and Perdomo design to the design
  references.
- Explain `root_completeness` separately from globality.
- State plainly that the public electrolyte dispatch is non-production and has
  not yet returned an admitted two-liquid result.

### `AGENTS.md`

No ownership change is required. Add only a short pointer to the canonical
phase-equilibrium index if that improves discoverability without duplicating
policy.

## Material deliberately not copied

- The M4 dated dashboard and issue queue.
- Old `bubble_pressure`, `dew_pressure`, and `single_component_vle` route maps.
- Old monorepo `packages/...` activation-owner paths.
- The preservation manifest's historical Ascani-as-electrolyte-owner invariant.
- Admission-registry family labels as runtime keys.
- Historical private adapters, archived fixed-two-phase routes, or old HELD2
  counterion-pair implementations.
- Validation case constants, Provider parameters, EOS equations, solver traces,
  artifact receipts, or broad paper archives.

## Verification design

Because this is a documentation-only change, verification focuses on authority
and schema consistency rather than new numerical tests:

1. Parse `ARCHITECTURE.yaml` with the repository's installed YAML tooling.
2. Check every new local Markdown link and referenced repository path.
3. Search current governing files for doctrine revision 2, stale old-monorepo
   public routes, false Perdomo-absent claims, or any new accepted electrolyte
   capability claim.
4. Run the existing documentation/static checks and the full source suite if
   repository tooling couples tests to architecture metadata.
5. Review the net diff for one owner per concept, no copied EOS/case constants,
   no runtime/API/build changes, and no broadened authority.
6. Run `git diff --check` and the repository cleanup hook.

## Success criteria

- A reader can find the package's governing phase-equilibrium documentation
  from the root README or `docs/phase-equilibrium.md`.
- Pure saturation, neutral HELD, and Perdomo HELD2 each have one package-local
  detailed design owner.
- All current governing files name doctrine revision 3.
- Lab M4 material is linked as provenance and cannot be mistaken for current
  runtime authority.
- Perdomo and Ascani coordinates remain explicitly non-equivalent.
- Accepted capability remains exactly `pure-component-saturation-v1`.
- HELD and HELD2 continue to report finite-search globality as not guaranteed.
- No source, tests, public symbols, native targets, artifacts, receipts, or
  scientific tolerances change.
