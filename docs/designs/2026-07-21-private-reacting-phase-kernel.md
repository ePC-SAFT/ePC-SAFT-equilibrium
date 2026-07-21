# Private Homogeneous Reacting-Phase Kernel

Status: approved private non-production design

Authority effect: none

Migration binding: D-028 at `807976e4f6f4b4f5a2c1ff4b1f0c699b5d22ea97`

## Scope

This design owns the first private fixed-temperature, fixed-pressure,
homogeneous reacting-phase foundation in Equilibrium. It adds no public Python
export, selector, result family, receipt, promotion, predictive claim, or
source-backed nonideal chemistry. The accepted public capability remains
`pure-component-saturation-v1`; D-026 remains the public HELD2 admission gate.

The first implementation slice contains only:

- a typed reaction/conservation/reference compiler;
- an exact positive electroneutral amount chart;
- a direct-potential phase-block boundary over installed Provider tensors;
- one native Ipopt problem with max-min initialization;
- exact postsolve and sensitivity-ready KKT evidence; and
- analytic, manufactured, and one installed-Provider manufactured test case.

Historical chemical-equilibrium code is mathematical evidence only. The
retired public route, selectors, workflows, internal Provider runtime,
`eos_x_gamma` solvent-only reference construction, MEA/Lithium fixtures, and
compatibility surfaces are forbidden inputs to this implementation.

## Compiled reaction system

Inputs are immutable and ordered: true-species identifiers, Provider species
identifiers, integer charges, conservation matrix `B`, declared independent
rank, independent reaction matrix `nu`, feed amounts, dimensionless `lnK`, and
complete source/reference records bound to the declared temperature and
pressure convention.

Compilation fails before solving unless:

1. Provider identifiers, order, charges, and fingerprint agree exactly;
2. all dimensions and scalars are finite and the feed is nonnegative;
3. `rank(B)` equals the declared rank;
4. `B nu^T = 0` and every reaction conserves charge;
5. `nu` has full row rank;
6. a declared complete closed system satisfies
   `rank(B) + rank(nu) = species_count`;
7. the feed is electroneutral; and
8. every equilibrium-constant record is nonempty, dimensionless, finite,
   source-identified, reference-identified, and bound to the declared `T` and
   optional declared `P`.

The compiler constructs the minimum-norm Provider-coordinate standard chemical
reference `g_ref` satisfying

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

## Initialization, solve, and continuation

A linear max-min Ipopt problem maximizes the smallest species amount subject to
the compiled balances and charge. Its independently recomputed optimum
adjudicates whether a strict-positive start exists above the declared trace
floor. No ad hoc epsilon species is inserted.

The homogeneous TNLP uses the amount chart plus `log(V)`, exact objective and
constraint derivatives, one fixed sparse schema, and exact Lagrangian Hessian.
It attempts the true objective first. A failed direct attempt may use one
adaptive continuation from a realizable analytic ideal/reference block to the
full phase residual. Only a final independently recomputed `lambda=1` attempt
can be accepted. Continuation is robustness evidence, never a nonideal or
global-equilibrium proof.

## Certificates and claim axes

The private result keeps these axes independent:

- artifact and input completeness;
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
Jacobian internally. It exposes no parameter sensitivity until Provider owns
the required typed parameter derivatives and the KKT system is demonstrated
nonsingular and acceptably conditioned.

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
