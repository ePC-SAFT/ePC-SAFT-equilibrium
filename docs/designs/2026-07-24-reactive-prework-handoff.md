# Reactive Prework Handoff to One Coupled Ipopt Solve

**Status:** approved design

**Date:** 2026-07-24

**Scope:** private homogeneous chemical-speciation preparation and its future
handoff to reactive multiphase equilibrium

## Decision

Phase-equilibrium preparation and chemical-speciation preparation remain
independent workflows. They produce compatible initialization and evidence
packets for one future simultaneous reactive multiphase Ipopt solve. Neither
preparatory workflow is an accepted coupled equilibrium result.

The current work implements only the source-bound homogeneous-speciation
preparer. It does not modify HELD2, add reactive phase equilibrium, expose a
public reaction API, or claim a predictive weak-electrolyte capability.

## Scientific distinction

The current Perdomo HELD2 formulation treats admitted strong electrolytes as
fully dissociated. The specified totals of those molecular and ionic species
are fixed, and phase equilibrium redistributes the fixed species among phases.

A weak-electrolyte or reactive formulation instead fixes independent elemental
or chemical-moiety totals. Reactions change individual species amounts while
preserving those totals. Strong-electrolyte ions may still participate in
protonation, complexation, ion exchange, ion pairing, precipitation, or other
admitted reactions.

The homogeneous reacting-phase problem is therefore a one-phase restriction of
the future coupled problem, expressed in true species amounts and conserved
chemical material.

## Current HELD2 phase preparer

The future phase-preparation input is the unified HELD2 workflow represented by
Equilibrium branch `codex/issue-65-stage2-handoff` at commit `368ffd3`.

Its relevant sequence is:

1. Select a valid homogeneous pressure root.
2. Run Stage I stability and tunneling searches and require a valid negative
   witness before phase-split work proceeds.
3. Run Stage II Step 4 upper LPs and Step 5 finite multistart lower searches.
4. Admit Step 5 lower states only when the physical/KKT and source stopping
   conditions pass.
5. Admit Step 6 candidates only when the same-major gap, fixed-volume
   gradient, pressure, and distinctness certificates pass.
6. Pass the candidate set and phase-coordinate bounds to Stage III.
7. Solve the active phase split, retire or merge only certified phases, polish
   pressure roots, and recompute KKT, complementarity, and free-energy-gap
   evidence.
8. Accept, return explicit feedback to Stage II, or terminate indeterminate.

Reactive work must preserve this controller logic. The coupled extension does
not fold source chemistry into HELD2's existing strong-electrolyte controller
or treat finite phase search as a global proof.

## Two preparers and one final solve

### Phase-equilibrium preparation

HELD2 supplies:

- selected pressure branches;
- stability witnesses and finite-search status;
- candidate phase compositions;
- candidate phase identities and distinctness;
- phase-coordinate and Provider-domain bounds;
- a phase split or phase-allocation warm start;
- explicit indeterminate and feedback outcomes.

The current strong-electrolyte HELD2 interface requires a species-level feed.
For future reactive initialization, Applications must supply a provisional
species projection consistent with the analytical feed, or Equilibrium may use
the homogeneous speciation preparer's positive state. This projection is an
initialization choice and is recorded as such. It does not replace the global
conserved-total constraints in the coupled problem.

Before the reactive extension exists, this remains strong-electrolyte phase
evidence. In the future coupled workflow it is initialization and topology
evidence, not a certificate for the final chemically minimized Gibbs surface.

### Chemical-speciation preparation

The private homogeneous workflow supplies:

- exact ordered true species and charges;
- the independent conservation matrix and feed totals;
- the independent reaction matrix;
- source equilibrium-constant records and declared standard states;
- an exact transformation to the installed Provider Helmholtz reference basis;
- a positive, exactly electroneutral, reaction-feasible species state;
- a positive volume on a valid Provider pressure branch;
- exact objective, constraint, and KKT derivative blocks;
- balance, charge, pressure, affinity, packing, trace, local-curvature, and
  source-identity evidence;
- explicit solver, numerical, physical, predictive, finite-search, and
  globality statuses.

This state may initialize one candidate phase or be scaled and allocated across
multiple candidates. It is not a sequential substitute for the final coupled
solve.

### Balance-preserving reconciliation

The later coupled implementation maps both preparation packets into one
strictly positive initial state. The reconciliation must:

- preserve global conserved totals exactly;
- preserve per-electrolyte-phase electroneutrality exactly;
- respect phase incidence and Provider component order;
- place each phase on its selected pressure/domain branch;
- keep every evaluated species above the declared trace floor;
- reject incompatible phase and species topologies;
- use max-min feasibility rather than application-fitted speciation vectors.

The reconciliation is initialization logic only. It does not adjudicate the
coupled thermodynamics.

### Final simultaneous reactive multiphase solve

For admitted phase \(p\), Equilibrium evaluates

\[
\mathcal G_p(n_p,V_p)
=
\Phi_p(T,n_p,V_p)
+\frac{P V_p}{R T}
+g_{\mathrm{ref}}^T n_p ,
\]

where the installed Provider owns the total mechanical Helmholtz block
\(\Phi_p\), pressure, packing, applicability, and exact nonlinear
derivatives. Equilibrium owns the chemical reference vector, conservation,
electroneutral coordinates, phase incidence, Ipopt formulation, KKT assembly,
and certificates.

The future coupled problem is

\[
\min_{\{n_p,V_p\}}
\sum_p \mathcal G_p(n_p,V_p)
\]

subject to

\[
\sum_p B n_p=b_{\mathrm{feed}},
\]

per-electrolyte-phase electroneutrality, positive amounts and volumes,
Provider applicability and packing constraints, and explicit phase incidence.

Its KKT conditions impose reaction equilibrium and phase-transfer equilibrium
together. A chemical solve followed by a flash is permitted only as
initialization.

One coupled Ipopt owner solves each proposed phase topology. A postsolve
reactive stability audit may discover another phase and return feedback to
phase preparation. That audit must permit trial phases to re-speciate under
their conserved totals; rerunning today's fixed-species HELD2 test without that
inner chemical freedom would not be a reactive stability certificate. Only the
final coupled solve after the last successful feedback cycle can be accepted.

## Shared numerical kernel

The homogeneous and coupled TNLP wrappers have different variable incidence,
constraint dimensions, and sparsity and therefore remain separate types.

They share:

- the Provider phase-block evaluator;
- the positive electroneutral amount chart and exact first/second chain rules;
- the positive volume coordinate;
- chemical-reference construction;
- conservation and KKT conventions;
- exact objective, constraint, and Lagrangian derivative assembly;
- Ipopt option and termination adjudication;
- callback-error and false-success rejection;
- postsolve physical and numerical certificate logic.

The common Ipopt runner is extracted only when both real consumers exist. This
avoids a speculative generic solver registry or a universal TNLP abstraction.

## Provider neutral-reference transport

The Provider EOS supplies a mechanical Helmholtz free-energy coordinate system.
For nonreactive phase equilibrium, per-species linear reference shifts cancel
because every species total is conserved. For reactions, relative species
reference shifts affect the equilibrium and must be supplied explicitly.

The installed Provider SDK's neutral-reference callback returns:

- a complete basis \(N\) for charge-neutral species combinations in exact
  Provider component order;
- dimensionless reference fugacity contractions for those basis rows;
- the reference composition, temperature, pressure, molality, and convergence
  evidence;
- the Helmholtz basis identity and exact parameter fingerprint.

It does not expose predictive individual-ion standard chemical potentials.
For every charge-balanced reaction row \(\nu_r\), Equilibrium solves

\[
\nu_r=\alpha_r N .
\]

This expresses the required reference correction using observable neutral
combinations, so the arbitrary single-ion electrical gauge cancels.

## Source-standard-state transformation

A source equilibrium constant is usable only when its units, dimensionality,
temperature, pressure convention, species order, activity convention, and
standard state are complete.

Using the retained transformation convention,

\[
\ln K_r^{\mathrm{Provider}}
=
\ln K_r^{\mathrm{source}}
+\nu_r s
+\alpha_r c ,
\]

where \(s\) is the declared source activity-scale shift and \(c\) contains the
Provider neutral-reference contractions.

The transformation must verify:

- exact Provider and source component identity and order;
- exact charges and charge-balanced reactions;
- finite source records bound to the declared state;
- representability of each reaction in the neutral basis;
- rank and conditioning of the neutral-basis representation;
- representation and reconstruction residuals;
- exact Provider basis and parameter fingerprint.

The reaction compiler then constructs a minimum-norm chemical reference vector

\[
\nu g_{\mathrm{ref}}=-\ln K^{\mathrm{Provider}} .
\]

Adding \(B^T q\) to \(g_{\mathrm{ref}}\) is a conservation gauge transformation
and must not alter the constrained equilibrium.

## Why the previous transformation surface was removed

The prior branch contained Provider neutral-reference transport,
source-standard-state transformation results, standalone private bindings, and
many tests. No qualified real source-bound reaction case consumed that stack
end to end. The retained solver evidence supplied synthetic equilibrium
constants already expressed in the Provider basis, while the MEA case remained
blocked by provisional species parameters and unknown applicability.

The same implementation also carried duplicated identity fields, preferred
start contracts, and source-specific test scaffolding. During the approved
code-surface reduction, the dormant layer was deleted instead of retained for
hypothetical future use. The Provider SDK callback remains available.

The transformation returns only as an active path from a complete source
record through compilation and speciation preparation. The old public shape,
standalone selector surface, compatibility behavior, and source-incomplete
fixtures do not return.

## Homogeneous source-bound formulation

For true species amounts \(n\), positive volume \(V\), conservation matrix
\(B\), feed totals \(b\), charge vector \(z\), and independent reaction matrix
\(\nu\), the compiler verifies:

\[
B\nu^T=0
\]

and

\[
z^T\nu^T=0 .
\]

It rejects dependent reactions, inconsistent equilibrium constants, incorrect
component identity or ordering, nonneutral feed, invalid source state,
incomplete reference records, and Provider-domain incompatibility before
Ipopt.

Every Provider evaluation uses the existing positive electroneutral chart:
positive total charge equivalents, cation and anion simplex shares, positive
neutral amounts, and log volume. Exact first and second chain rules are
production requirements.

The existing max-min problem supplies a strictly positive feasible state. The
solver attempts the true Provider objective first. If necessary, adaptive
continuation moves from the realizable ideal/reference objective to the full
residual EOS. Only a final \(\lambda=1\) Ipopt solve can pass.

Postsolve acceptance independently recomputes conservation, exact charge,
pressure, reaction affinities, positivity, Provider domain and packing, KKT,
trace-boundary status, and reduced curvature. Local KKT or continuation success
does not prove globality or predictive agreement.

## First implementation slice

The next local slice:

1. Restores the minimum Provider neutral-reference transport required by a
   source-bound input.
2. Restores source-standard-state transformation as a private function used
   directly by the preparation path.
3. Composes source transformation, reaction compilation, and the existing
   Provider homogeneous solve without adding a second solver.
4. Retains phase-assembly-ready information internally: conserved totals,
   positive species amounts, volume, reference identity, and KKT evidence.
5. Adds one parameterized basis/order/gauge-invariance test and one
   parameterized fail-closed identity/reference/domain test.
6. Adds a real installed-artifact sentinel only when the Provider parameters,
   applicability, source equilibrium constants, and reference transformation
   are complete.

The provisional MEA bundle is not the first accepted sentinel. MEA parameter
qualification remains application and Provider work. The generic preparer must
be ready to consume that bundle later without a chemistry-named native route.

## Deferred work

This design does not authorize:

- reactive HELD/TPD implementation;
- coupled chemical and phase equilibrium;
- reactive bubble equilibrium;
- lithium extraction chemistry or DES-specific species;
- precipitation, stripping, kinetics, or interfacial behavior;
- parameter regression;
- a public reaction API or result family;
- a backend/plugin registry;
- a source-incomplete MEA solve;
- predictive, finite-search-completeness, or globality claims;
- promotion, release, receipt, or authority changes.

## Evidence policy

Persistent package evidence remains compact:

- one source-transform invariance family;
- one source/reference/domain rejection family;
- the existing compiler, chart, derivative, manufactured-solver, and installed
  Provider-block obligations;
- one real source-complete sentinel when available.

Full MEA or lithium sweeps, plots, parameter studies, and application
comparisons belong in application or serialized installed-artifact validation
work, not the package unit-test surface.

## Success criteria

The speciation preparer is ready for later coupled assembly when:

- a complete source record is transformed into the exact installed Provider
  basis without a single-ion or epsilon-composition construction;
- the transformed system reaches the existing homogeneous Ipopt kernel through
  one private path;
- all source, identity, conservation, charge, derivative, domain, KKT, trace,
  and local-curvature gates fail closed;
- the result carries enough internal state for a later balance-preserving
  phase allocation without adding coupled code now;
- manufactured and incomplete application inputs remain explicitly
  nonpredictive;
- no public surface, chemistry-specific native branch, or duplicate solver is
  introduced.
