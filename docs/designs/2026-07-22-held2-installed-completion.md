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

## Completion invariants

The remaining implementation and evidence must preserve all of the following:

1. No ePC-SAFT equation, density implementation, or model parameter is copied
   into Equilibrium.
2. The pressure certificate remains `1e-8` relative, physical Stage-II KKT
   tolerance remains `1e-7`, and complementarity remains `1e-8`.
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

`enumeration_objective_gap` or an equivalent field must be computed from
independent quantities; a default zero is not evidence. A manufactured test
must demonstrate rejection when local KKT and physical equalities pass but the
free-energy gap fails.

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
