# Reactive Prework Handoff to One Coupled Ipopt Solve

**Status:** approved architecture; pre-issue sanity review incorporated

**Date:** 2026-07-24

**Scope:** private homogeneous chemical-speciation preparation and its future
handoff to reactive multiphase equilibrium

## Decision

Phase-equilibrium preparation and chemical-speciation preparation remain
independent workflows. They produce compatible initialization and evidence
packets for one future simultaneous reactive multiphase Ipopt solve. Neither
preparatory workflow is an accepted coupled equilibrium result.

The current work first makes the private homogeneous-speciation kernel
trace-ready against the source-derived Belov-Aristova ideal-gas benchmark. A
separate source-bound ePC-SAFT integration slice remains deferred until a real,
source-complete, non-MEA installed Provider consumer exists. Neither slice
modifies HELD2, adds reactive phase equilibrium, exposes a public reaction API,
or claims a predictive weak-electrolyte capability.

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

The current inspection baseline for the future phase-preparation input is the
provisional Equilibrium branch `codex/issue-65-stage2-handoff` at commit
`368ffd3`. It is useful implementation evidence, not a stable dependency or
authority record. The homogeneous chemical-equilibrium work does not depend on that
branch. A later coupled design must bind an exact merged HELD2 artifact before
it consumes any phase-preparation contract.

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
- mass-first independent conservation rows and feed totals;
- supplied and independent reaction matrices, basis maps, and cycle evidence;
- source equilibrium-constant records and declared standard states;
- constants already converted to the installed Provider Helmholtz reference
  basis; source-standard-state transport remains deferred;
- exact homogeneous support classifications and accessible-species maps;
- a positive, exactly electroneutral, reaction-feasible accessible state;
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
+g_{\mathrm{ref},p}^T n_p ,
\]

where the installed Provider owns the total mechanical Helmholtz block
\(\Phi_p\), pressure, packing, applicability, and exact nonlinear
derivatives. Equilibrium owns the chemical reference vectors, conservation,
electroneutral coordinates, phase incidence, Ipopt formulation, KKT assembly,
and certificates. When every phase uses the same Provider model and Helmholtz
basis, the vectors may be one common \(g_{\mathrm{ref}}\). When aqueous and
organic phases use different admitted Provider models, every
\(g_{\mathrm{ref},p}\) must be transformed into one common global chemical
reference convention. Independent per-phase gauges are forbidden because they
would change reaction and transfer equilibria. Only a conservation-space gauge
shift applied consistently across all phases is thermodynamically inert.

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

- exact installed SDK ABI/table/result sizes and a non-null callback;
- exact Provider and source component identity and order;
- exact charges and charge-balanced reactions;
- finite source records bound to the declared state;
- representability of each reaction in the neutral basis;
- rank and conditioning of the neutral-basis representation;
- representation and reconstruction residuals;
- exact returned capacities, basis-row count, and finite callback outputs;
- exact Provider basis and parameter fingerprint.

The per-species activity-scale shift \(s\) is admitted only when the declared
source conversion is representable by such a state-bound linear shift.
Composition-dependent or otherwise nonrepresentable reference changes fail
closed rather than being frozen into constants.

The current installed SDK reports no neutral-reference temperature, pressure,
composition, or parameter derivatives. That is sufficient for the present
fixed-\(T,P\) value transformation. It does not support caloric derivatives,
implicit equilibrium sensitivities, or Regression parameter derivatives.

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

## Homogeneous formulation

For true species amounts \(n\), positive volume \(V\), conservation matrix
\(B\), feed totals \(b\), charge vector \(z\), and supplied reaction matrix
\(\nu^{\mathrm{sup}}\), the compiler first validates mass, balance, charge, and
converted reaction cycles, then records one independent basis \(\nu\) and
verifies:

\[
B\nu^T=0
\]

and

\[
z^T\nu^T=0 .
\]

For the complete closed homogeneous system it also requires

\[
\operatorname{rank}(B)+\operatorname{rank}(\nu)=C ,
\]

where \(C\) is the number of admitted true species. Redundant supplied
reactions are accepted only when their converted constants reconstruct from
the retained basis. The compiler rejects an incomplete
conservation/charge/reaction span, inconsistent cycles, incorrect component
identity or ordering, nonneutral feed, invalid source state, incomplete
reference records, and Provider-domain incompatibility before Ipopt. A future
phase-incidence topology must repeat the rank and span analysis for the species
actually admitted in each phase and the global conservation system; the
homogeneous rank certificate cannot simply be reused.

Every Provider evaluation uses the existing positive electroneutral chart:
positive total charge equivalents, cation and anion simplex shares, positive
neutral amounts, and log volume. Exact first and second chain rules are
production requirements.

The existing deterministic HiGHS max-min LP supplies a strictly positive
feasible state. The solver makes one direct Ipopt solve of the true Provider
objective. It owns no ideal-to-Provider continuation or fallback optimizer.

Postsolve acceptance independently recomputes conservation, exact charge,
pressure, reaction affinities, positivity, Provider domain and packing, KKT,
trace-boundary status, and reduced curvature. Local KKT success does not prove
globality or predictive agreement.

The current logarithmic amount chart still solves only the strictly positive
accessible species. Before solving, a homogeneous support layer exhaustively
enumerates bounded feasible bases with exact binary-rational arithmetic as the
deletion authority. A species is removed only when exhaustive enumeration
finds no positive feasible vertex; an enumeration limit fails compilation.
Reaction combinations are recomputed in the null space of the removed-species
columns, so cancellation through absent intermediates is preserved.

The manufactured ideal restriction supports this exact structural face and
restores zero amounts in original species order. Installed Provider
reduced-component faces remain `BOUNDARY_DIRECTION_UNRESOLVED` because SDK v1
does not expose reduced-topology or physical boundary-direction oracles.
Nonstructural trace and Provider-domain boundaries remain unaccepted.

Belov-Aristova test 1 is mathematically interior even though its smallest
amount is about `1.84e-47 mol`. The implemented generic KKT polish now resolves
that extreme positive trace range after the exact-Hessian Ipopt solve and
independently rebuilds physical affinities, multipliers, stationarity, and
curvature. It does not classify a small positive amount as structurally zero.

## First implementation slice

GitHub issue #69 owns the next local slice:

1. Record the exact authors' `test1.dat` identity, gas-species restriction,
   element matrix, fixed state, and `1 atm` ideal-gas standard state.
2. Convert the source values exactly into the manufactured kernel's
   `1 mol/m3` concentration reference.
3. Retain an independent high-precision element-potential oracle whose expected
   values are not generated by the implementation under test.
4. Add the focused red benchmark showing that current transformed-coordinate
   stationarity can coexist with failed physical reaction affinities.
5. Change only the generic private homogeneous formulation or coordinate
   treatment needed to resolve the extreme positive trace range.
6. Recompute conservation, pressure, all independent affinities,
   physical-coordinate KKT, trace status, and local curvature; do not loosen
   acceptance to accommodate the old pinned result.
7. Prove reaction-basis, species-order, conservation-gauge, and feed-scale
   invariance with the smallest persistent evidence surface.

The implemented generic resolution preserves the exact-Hessian Ipopt
minimization as the sole optimization solve. If that solve is feasible and
strictly above the trace floor but fails the independently recomputed physical
reaction-affinity certificate, a damped Newton pass imposes the square system
of independent balances, reaction affinities, and pressure directly using the
existing exact Jacobian blocks. A polished state is adopted only after physical
affinities pass; Provider-domain inequalities, equality multipliers,
stationarity, complementarity status, and reduced curvature are then
reconstructed independently from that state. Active-bound or active-domain
cases remain unaccepted.

Belov is an ideal-gas algorithm and trace benchmark, not an installed-Provider
or ePC-SAFT validation case. Reaktoro may be used as a disposable independent
cross-check only; it is not a production dependency or oracle of record.

After #69, source-reference transport returns only through a newly shaped leaf
whose exact non-MEA installed Provider artifact and complete source record
exercise the transformation through the existing homogeneous solve. Restoring
that path before a real consumer exists would recreate the dormant surface that
was deliberately removed. MEA qualification remains application and Provider
work and does not define the generic Equilibrium baseline.

## Issue boundary and readiness

GitHub issue #35 owns private homogeneous chemical-equilibrium readiness. Issue
#69 is the executable Belov trace-readiness leaf and is independent of Provider
source-reference work. Historical MEA-specific leaves #46, #47, and #68 do not
define the package baseline; #68 is retired as not planned.

A future source-reference leaf is Ready only when its body names:

- the exact non-editable non-MEA Provider artifact and native SDK contract;
- one complete source record and its standard-state conversion;
- the existing homogeneous solver entry point that consumes the transformed
  values;
- the two compact persistent test families;
- the explicit interior-equilibrium and no-sensitivity limits.

Without that real consumer, no executable restoration leaf should be opened and
implementation must not restore unused contracts or tests.

The worktree currently binds the immutable Provider wheel with SHA-256
`a8b6376193301673429a8d8b648896b6881c6f693ca4d8fc3f6a9d2d16f3b39c`.
That artifact exposes the neutral-reference callback but does not contain the
historical `held-2008-water-self-ionization` catalog bundle. The earlier IAPWS
water self-ionization solve remains useful provenance, but it is not a current
installed-artifact consumer and cannot make the new leaf Ready.

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

- one source-derived Belov trace-readiness benchmark with an independent
  element-potential oracle;
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

- the source-derived Belov homogeneous restriction passes its physical
  affinity, conservation, pressure, trace, KKT, and invariance gates;
- redundant converted reactions, exact structural support, and accessible-face
  reaction combinations retain complete compiler evidence;
- a complete source record is transformed into the exact installed Provider
  basis without a single-ion or epsilon-composition construction;
- the transformed system reaches the existing homogeneous Ipopt kernel through
  one private path;
- all source, identity, conservation, charge, derivative, domain, KKT, trace,
  and local-curvature gates fail closed;
- the existing compiled system and solve result retain the state already
  justified by homogeneous acceptance, without speculative coupled fields;
- manufactured and incomplete application inputs remain explicitly
  nonpredictive;
- no public surface, chemistry-specific native branch, or duplicate solver is
  introduced.
