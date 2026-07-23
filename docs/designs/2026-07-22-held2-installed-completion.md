# HELD2 Installed Completion Specification

## Status and authority

This specification defines the remaining private Equilibrium work after the
integrated HELD2 candidate recorded by issue #22. It does not admit an
electrolyte capability, authorize Validation mutation, alter Provider
ownership, or convert finite search into a global proof.

Perdomo remains authoritative for the strong-electrolyte HELD2 formulation.
The installed Provider remains authoritative for ePC-SAFT values, domains,
pressure, exact first and second derivatives, and volume transformations.
Equilibrium owns modified coordinates, pressure-root enumeration, stage
orchestration, optimizers, independent certificates, and diagnostics.

## Verified starting point

The retained candidate binds Equilibrium implementation commit
`40d01d32dbcda67caa5eafc71b04abe74b80f016`, Equilibrium wheel SHA-256
`606686b766162a4366b2274113a954329ba1496a9b4a27f819b5e2c030979188`,
and Provider wheel SHA-256
`9e4da0d7ba7896bcd2ec096400553d935e0516c61f1bd9f41f2370ab68ab36ea`.

For the frozen Khudaida feed it establishes:

- three detected pressure roots with stable--unstable--stable mechanical
  classes;
- selection of the dense strict-stable reference at
  `3.909560419950766e-5 m3/mol` and reduced objective
  `-6.009742067562037`;
- a pressure-certified negative Stage-I TPD witness with detected value
  `-0.008692863609313706`;
- successful Stage-II reduced-envelope exploration with 19 evaluations, no
  exploration failure, and 37 physical representatives; and
- a fail-closed Stage-II lower-search outcome after the one declared local
  attempt encountered `HELD2 Stage II chart coordinate is outside [0, 1]`.

The retained candidate does not establish a real Provider Stage-II candidate
set, a real Provider Stage-III solve, physical two-phase equilibrium, or
predictive agreement with Khudaida. Its root-completeness and globality claims
remain `not_proven` and `not_guaranteed`, respectively.

## Corrected Step-5 through Step-9 status

The current implementation no longer solves Problem (65) in a fixed volume box
borrowed from the start composition. It minimizes a pressure-root-reduced
objective in the independent modified-composition chart. Each trial
composition obtains its own Provider volume domain, follows a numbered strict
stable pressure branch with a safeguarded exact-derivative continuation, and
falls back to the complete declared pressure-envelope search when branch
continuation cannot be certified. The reduced gradient and Hessian are the
exact pressure-manifold derivatives; the Hessian uses the Schur complement of
the eliminated volume coordinate. The final state is reconstructed at the
polished pressure root before physical KKT and Step-6 certification.

Problem (65) is qualified only when its independently certified local value is
not above the same-major Problem-(64) upper value, apart from the existing
named Step-6 gap allowance. A locally converged point above that upper value is
retained as failed evidence and cannot become a lower bound, cut, or candidate.
The first qualified lower state supplies the next complete-cut iteration.
When a Step-6 state is found, the same major continues through distinct
representatives to look for a second co-minimizer.

Provider-valid pressure-envelope representatives are eligible as affine dual
cuts but not as candidate phases. Candidate eligibility additionally requires
the exact-Hessian Ipopt terminal, pressure closure, original-coordinate KKT,
physical multiplier signs, complementarity, dual reconstruction, same-major
gap, and fixed-volume Eq. (66) gradient tests. Composition-rich simplex-vertex
seeds are included so aqueous-rich and organic-rich basins do not depend on a
fortunate Sobol point. HiGHS solves Problem (64) with a `1e-10` internal
primal/dual target so its cut feasibility is resolved more tightly than the
unchanged `1e-8` Step-6 gap gate; the public tolerance contract is unchanged.

Manufactured end-to-end evidence now passes the Step-6 candidate set directly
to the generic Problem-(67) owner. Step 8 converges with exact Hessians, and
Step 9 accepts two phases only after modified and explicit material balances,
charge, pressure, modified-potential equality, KKT, complementarity, phase
identity, and active-set lifecycle checks pass. Deliberate potential,
infeasible-set, trace-component, duplicate, and inactive-phase cases remain
fail-closed.

A private Perdomo Table-5 LiCl--water--1-butanol screening replay is useful only
as numerical evidence because Perdomo used SAFT-gamma Mie and no
source-equivalent Provider bundle exists. The declared ePC-SAFT hypothesis
combines published solvent and ion records with explicit unfitted screening
assumptions. The corrected search finds and pressure-certifies both an
organic-rich and an aqueous-rich local basin. The prior controller nevertheless
executed only the first Step-5-qualified local solve in each major whenever that
state missed Step 6, despite declaring up to 50 attempts. The corrected
controller continues through distinct representatives until two Step-6 states
are found or the declared attempt cap is exhausted. In the 24-major private
profile, both basins pass Eq. (66) in major 19 with same-major gaps below
`4.1e-9`.

The resulting generic Step-8 solve converges to two active phases. Before Step
9, a bounded exact-derivative corrector polishes only active phase log-volumes
on their already selected strict-stable pressure branches, then recomputes the
NLP KKT and complementarity evidence. It does not alter compositions,
balances, Provider equations, or physical gates. The private screening replay
then passes material and charge closure, pressure, modified-potential equality,
KKT, phase identity, and the independently computed Eq. (68) gap. This is
end-to-end numerical evidence for the implemented controller under an
unadmitted ePC-SAFT parameter hypothesis; it is not reproduction or validation
of Perdomo's SAFT-gamma-Mie result, and its
`globality_certificate` remains `not_guaranteed`.

The subsequent issue-#27 replay classified the retained chart failure as
solver-bound numerical contact, not a physical-simplex excursion. The earliest
rejected trial was an `eval_f` request at `1.0000000000000002`, was not an
accepted iterate, and exceeded the unit bound by one binary64 ULP. Continuing
the same replay exposed only the second and fourth representable values above
`1.0`; the maximum excursion was `8.881784197001252e-16`. Normalizing those
contacts to the exact unit boundary left the mapped physical state unchanged,
removed the callback failure, and let Ipopt terminate successfully. Issue #32
supersedes that ULP-count runtime rule with the canonical named
`chart_contact=1e-9` representation gate. Non-finite or materially larger
excursions remain invalid. The one-attempt controller is still Stage-II
indeterminate because that diagnostic budget does not produce an eligible
candidate set.

The issue-#32 exact-wheel replay records a Stage-II dual-pullback residual of
`1.2304377718370303e-10`, original-coordinate stationarity of
`1.862645149230957e-9`, complementarity of `9.090909078162636e-12`, and
relative pressure residual of `-1.3307506742421539e-9`. Each independently
passes its named gate. The same terminal remains ineligible under Step 6: its
same-major gap is `132.54249236307797` and its fixed-volume gradient residual
is `3853783.6445812103`. The categorized reconstruction allowance therefore
does not convert this local KKT state into a candidate phase.

An offline `0.1x/1x/10x` sensitivity replay of the same immutable raw evidence
certifies `1/3/3` pressure roots and `0/0/0` Stage-II candidates. The tighter
`0.1x` profile changes pressure-root certification because two root residuals
lie between `1e-9` and `1e-8`; nominal and `10x` preserve the detected
stable--unstable--stable topology. No profile changes the accepted candidate
or phase count. This sensitivity result is diagnostic evidence, not a runtime
tuning mechanism.

## Completion invariants

The remaining implementation and evidence must preserve all of the following:

1. No ePC-SAFT equation, density implementation, or model parameter is copied
   into Equilibrium.
2. Every numerical and physical gate uses the categorized contract indexed in
   `docs/phase-equilibrium.md`; no leaf may tune or substitute those values.
   Pressure remains `1e-8` relative, Stage-II stationarity remains `1e-7`,
   complementarity remains `1e-8`, and dual pullback is independently scaled.
3. Fixed-physical-volume Step-6 derivatives and same-major `UBD/lambda` data
   remain authoritative.
4. A chart, Provider, root, solver, resource, certificate, or trace-refinement
   failure remains visible and fail-closed.
5. Ipopt success is neither necessary nor sufficient for physical acceptance.
6. Stage-II candidate distinctness is composition-or-volume/packing
   distinctness in physical coordinates.
7. Stage III retires only KKT-inactive phases, changes at most one active-set
   member per lifecycle step, and re-solves every changed active set.
8. Predictive comparison is prohibited until a physical equilibrium state has
   passed every numerical and thermodynamic certificate.
9. Every finite search reports `globality_certificate=not_guaranteed`.

## Stage-II boundary and diagnostic completion

The first unresolved defect is the exact installed Stage-II trial that leaves
the unit-cube chart. The repair must begin from the retained start, upper bound,
and multipliers and record every local evaluation in order. It must distinguish
an intermediate Ipopt trial from an accepted iterate, quantify the bound
violation, identify whether mapping or solver-bound semantics are responsible,
and retain the last valid physical state.

The implementation must not clamp a scientifically meaningful violation,
relax a bound or tolerance, replace Ipopt, or convert the failure to an
objective penalty. Solver-permitted roundoff may be handled only through a
declared, tested coordinate policy that cannot silently change the physical
state.

Every successful Stage-II upper LP must be recorded before the lower search
begins. Its solver/version, primal and dual feasibility, residuals, cut slacks,
cut duals, active cuts, upper bound, and multipliers must survive any later
lower-search failure or resource exit.

The implemented coordinate policy is `held2_chart_contact_abs_v2`. It applies
only the named `1e-9` representation allowance to a finite solver iterate next
to the closed unit chart; it is not an Ipopt bound relaxation and it cannot
make an invalid direct user composition valid. Ordered local evidence retains
the tolerance audit, raw coordinates, callback owner, trial-versus-accepted
status, maximum violation, mapped physical coordinates, and the last valid
physical state. Problem-(64) evidence is committed to the major history before
Problem (65) starts and marks its lower value unavailable until a certifiable
lower terminal exists.

## Full installed Stage-II candidate contract

After the boundary defect is repaired, the controller must run a declared
full-budget installed Stage-II campaign rather than the one-attempt evidence
profile. Every representative and local terminal, including failures, must be
retained with stable IDs and start-family provenance.

A Stage-II candidate set is acceptable only when at least two physically
distinct states pass:

- Provider and domain validity;
- pressure closure;
- original-coordinate KKT stationarity;
- physical lower, upper, and simplex multiplier signs;
- complementarity and dual-pullback reconstruction;
- same-major lower/upper gap;
- fixed-volume Step-6 component conditions; and
- composition-or-volume/packing distinctness.

Budget exhaustion, an unaccounted Provider region, or inability to produce two
eligible phases remains indeterminate. It is not a one-phase result.

## Installed Stage-III completion contract

The generic Stage-III exact-Hessian NLP must consume the installed Stage-II
candidate set without a manufactured adapter. Its exact derivatives must be
checked directionally through the installed Provider callback.

The final certificate must include:

- modified and explicit material balances;
- phase electroneutrality;
- common pressure;
- equality of every independent non-trace modified potential;
- phase-amount bound KKT, signs, and complementarity;
- phase distinction;
- Provider-domain validity;
- active-set lifecycle and re-solve evidence; and
- a computed source-faithful total-free-energy/HELD gap.

The source Eq. (68) quantity is
`same_major_stage_ii_UBD - stage_iii_total_free_energy`. The generic owner
retains the upper bound, total free energy, signed gap, and provenance. A
default zero is not evidence. The installed derivative seam checks the generic
Problem-(67) gradient and Lagrangian Hessian through the exact Provider
callback. Manufactured tests independently enumerate the upper bound,
demonstrate acceptance at the matching total free energy, reject a perturbed
upper bound while local KKT and physical equalities still pass, and reject
unavailable gap evidence.

The private Perdomo Table-5 ePC-SAFT screening hypothesis reaches two distinct
same-major Eq. (66) candidates in major 19 of the corrected 24-major campaign.
The generic installed Step 8 therefore runs and, after the active-phase
pressure correction described above, passes every Step-9 certificate and the
independently computed Eq. (68) gap. This is numerical controller evidence for
the exact private ePC-SAFT hypothesis, not a source-equivalent SAFT-gamma-Mie
reproduction or a globality proof.

If the exact installed case reaches the declared trace threshold, logarithmic
or complementarity-safe trace refinement requires its own derivative and KKT
evidence. Until that evidence passes, `trace_component_requires_log_refinement`
remains a fail-closed result.

## End-to-end private evidence contract

One immutable Equilibrium wheel and one immutable Provider wheel must execute
the complete private controller at least twice in fresh, non-editable
environments. The retained evidence must include every root, Stage-I
evaluation, Stage-II cut/start/basin/terminal/candidate, Stage-III lifecycle
step, solver status, certificate, failure, and resource count.

Physical acceptance requires at least two retained phases and every certificate
above. Only then may the evidence convert each explicit-ion phase back to the
Khudaida formula-component basis and classify the result as predictive
agreement, predictive disagreement, or model-topology disagreement. That
comparison remains private local evidence and creates no capability authority.

## Required evidence layers

The scientific argument is complete only when all four layers are distinct:

1. **Manufactured formulation evidence:** analytic roots, derivative checks,
   LP oracles, KKT signs, active sets, and deliberate physical failures.
2. **Installed numerical evidence:** exact wheel/import hashes, full resource
   accounting, derivative checks, reproducibility, and certificate residuals.
3. **Physical equilibrium evidence:** balances, charge, pressure, modified
   potentials, active phases, distinction, and free-energy gap.
4. **Predictive evidence:** per-phase basis conversion and comparison with the
   frozen experimental tie line only after physical acceptance.

## Non-goals

This completion work does not authorize public electrolyte dispatch,
Validation repository mutation, capability promotion, release publication,
parameter fitting, tolerance relaxation, a Khudaida-specific solver path, or a
mathematical global-stability claim.
