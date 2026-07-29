# Private Homogeneous Reacting-Phase Kernel

Status: approved private base design; compiler and structural-boundary sections
superseded by `2026-07-27-grepe-homogeneous-chemical-layer.md`

Authority effect: none

Migration binding: D-028 at `807976e4f6f4b4f5a2c1ff4b1f0c699b5d22ea97`

## Scope

This design owns the first private fixed-temperature, fixed-pressure,
homogeneous reacting-phase foundation in Equilibrium. It adds no public Python
export, selector, result family, receipt, promotion, or predictive claim. The
accepted public capability remains
`pure-component-saturation-v1`; D-026 remains the public HELD2 admission gate.

The first implementation slice contains only:

- a typed reaction/conservation/reference compiler;
- an exact positive electroneutral amount chart;
- a direct-potential phase-block boundary over installed Provider tensors;
- one native Ipopt problem with max-min initialization;
- exact postsolve and sensitivity-ready KKT evidence; and
- analytic and manufactured tests, one installed-Provider manufactured case,
  and one source-complete installed-Provider value case.

Historical chemical-equilibrium code is mathematical evidence only. The
retired public route, selectors, workflows, internal Provider runtime,
`eos_x_gamma` solvent-only reference construction, MEA/Lithium fixtures, and
compatibility surfaces are forbidden inputs to this implementation.

Implementation update (2026-07-27): the private GREPE homogeneous layer now
accepts redundant supplied reaction rows, validates converted reaction cycles,
adds explicit molar-mass conservation, certifies structural support with exact
rational evidence, and recompiles reaction combinations on the accessible
species face. The original strictly positive amount chart and sole Ipopt solve
remain unchanged. Exact structural zeros are supported only for the
manufactured ideal restriction; installed Provider reduced-component faces
remain fail-closed and unresolved.

## Compiled reaction system

Inputs are immutable and ordered: true-species identifiers, integer charges,
positive molar masses, a supplied conservation matrix `B`, a possibly
redundant supplied reaction matrix `nu`, feed amounts, dimensionless
Provider-basis `lnK`, a Provider fingerprint, and complete source/reference
records bound to the fixed temperature and pressure.

Compilation fails before solving unless:

1. installed Provider identifiers, order, charges, and fingerprint agree
   exactly at the Provider solve boundary;
2. all dimensions and scalars are finite and the feed is nonnegative;
3. balances are rank-reduced while retaining molar mass first and charge
   separately;
4. every supplied reaction conserves mass, balances, and charge;
5. converted constants satisfy every redundant reaction cycle before basis
   selection;
6. the independent invariant, charge, and reaction spaces span the declared
   species;
7. the feed is electroneutral; and
8. every equilibrium-constant record is nonempty, dimensionless, finite,
   source-identified, reference-identified, and bound to the fixed `T,P`.

After exact structural-support classification and accessible-face
recompilation, the compiler constructs the minimum-norm Provider-coordinate
standard chemical reference `g_ref` satisfying

```text
nu g_ref = -lnK
```

and verifies the reconstruction residual. `g_ref + B^T c` is a gauge-equivalent
reference and must leave every constrained equilibrium and reaction affinity
unchanged.

## Positive amount coordinates

For a neutral-only system, each amount is `exp(y_i)`. For an ionic system, the
chart requires at least one positive and one negative charge and uses:

```text
Q = exp(q)
n_c = Q alpha_c / z_c
n_a = Q beta_a / abs(z_a)
n_neutral = exp(y_neutral)
```

`alpha` and `beta` are reference-category softmax simplexes. This spans the
strictly positive electroneutral amount manifold without choosing a
case-specific counterion. The chart owns exact amount values, its Jacobian, and
every component Hessian. Production evaluations use only these analytic chain
rules. A finite lower coordinate bound is a numerical trace floor, reported as
a boundary classification and exercised by sensitivity tests; it is not a
chemical absence model.

## Thermodynamic block and objective

The single physical objective is

```text
J(n,V) = Phi(T,n,V) + P V/(R T) + g_ref^T n,
```

subject to `B n = B n_feed`, with charge enforced identically by the chart.
`Phi` is one total mechanical `A/(RT)` phase block in Provider coordinate order
`(n_1, ..., n_C, V)`. A phase block returns exact value, gradient, Hessian,
pressure, volume/packing applicability, and fingerprint evidence. Equilibrium
applies only the exact coordinate chain rule. It does not rebuild pressure,
chemical potentials, association, electrolyte terms, density closure, packing
constants, or EOS derivatives.

Manufactured tests may use the analytic ideal block
`sum_i n_i (log(n_i/V) - 1)`. Installed-Provider tests supply synthetic
`g_ref`/`lnK` and are labeled manufactured and nonpredictive. Provider's current
Helmholtz ideal reference is not relabeled as a reactive standard state.

The source path accepts explicit, state-bound source activity-scale shifts and
uses only charge-neutral contractions from the installed Provider SDK. It
solves each reaction row in the Provider neutral basis, applies the resulting
standard-state offset, rewrites the provenance to the Provider Helmholtz
coordinate basis, and enters the same compiler and solver. The retained
source-complete sentinel is Held-2008/IAPWS R11-07(2019) water
self-ionization at 298.15 K and 1 bar. No absolute single-ion standard is
constructed. The appended SDK tail supplies a certified neutral-reference
pressure derivative, so the implicit pressure column includes the exact
source-reference contribution at the actual system pressure. This does not
supply caloric, composition, or active-parameter reference derivatives.

## Initialization and solve

A deterministic HiGHS max-min LP maximizes the smallest species amount subject
to the compiled balances, charge, and Provider source-domain limit. Its
independently recomputed optimum adjudicates whether a strict-positive start
exists above the declared trace floor. The positive molar-mass balance supplies
finite amount bounds; no per-species optimization or ad hoc epsilon species is
inserted.

The homogeneous TNLP uses the amount chart plus `log(V)`, exact objective and
constraint derivatives, one fixed sparse schema, and exact Lagrangian Hessian.
It makes one direct solve of the true Provider objective. A failed direct solve
remains failed; this private kernel owns no ideal-to-Provider continuation or
fallback optimizer.

## Certificates and claim axes

The private result keeps these axes independent:

- fail-before-solve artifact and input validation;
- raw solver termination and callback state;
- original-coordinate numerical feasibility and KKT stationarity;
- conservation, exact charge, pressure, positivity, Provider domain/packing,
  and reaction-affinity validity;
- trace-floor status;
- reduced-Hessian local-minimum status;
- predictive agreement, normally `not_adjudicated` for this slice;
- finite-search completeness; and
- globality, always `not_guaranteed`.

Ipopt success, finite multistarts, continuation, or a local KKT point never
becomes a global thermodynamic proof. Acceptance is fail-closed: a non-success
status, missing candidate, callback/domain failure, failed independent
recomputation, or certificate failure cannot return a successful kernel result.

The implementation assembles the exact equality-constrained KKT residual and
Jacobian internally. For a certified strict-interior solution it solves that
system for the deterministic parameter order of compiled balance totals,
Provider-basis reaction constants, and pressure, then maps the result through
the exact amount-chart Jacobian and logarithmic volume coordinate. The private
result records the KKT rank, infinity-norm condition number, active variable,
domain-constraint, and trace bounds, chart topology, and immutable Provider
parameter fingerprint. It returns no derivative when the KKT matrix is
singular or unacceptably conditioned, the active set can change, or the primal
result is not certified.

Provider parameters require typed phase-KKT cross derivatives and the
corresponding derivative of the transformed Provider reference vector. The
current tail supplies the pressure/reference column but not an atomic
active-model callback with complete Helmholtz and packing state tensors.
Active Provider parameter families therefore remain explicitly unavailable;
Equilibrium does not synthesize them with finite differences. The public
`chemical_equilibrium` operation optionally exposes the exact conditioned
totals, Provider-basis `ln(K)`, and pressure columns through its typed
sensitivity result.

## Later reuse and explicit deferrals

A later simultaneous phase-chemical formulation may sum the same per-phase
block under global conservation, phase incidence, per-electrolyte-phase charge,
and transfer/reaction/pressure certificates. Sequential speciation then flash
is initialization only. HELD remains the phase-discovery owner.

The required downstream order is: source-complete reacting liquid; reactive
bubble plus exact implicit sensitivities for MEA; simultaneous phase-specific
two-liquid reactions for Lithium; and mixed-observable Regression only after
typed Provider parameter derivatives. None of those subjects is part of this
slice.

Structural support and accessible-face compilation do not authorize coupled
phase equilibrium. The homogeneous result is one local candidate. Phase-family
generation, master pricing, state-matched repricing, simultaneous multiphase
constraints, and rigorous global bounds remain separate work.
