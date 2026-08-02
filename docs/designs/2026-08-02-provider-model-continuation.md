# Provider-model continuation for homogeneous reaction equilibrium

Status: pre-result implementation contract

This design adds one chemistry-agnostic recovery route to the existing public
fixed-`T,P`, single-phase `chemical_equilibrium` operation. It does not add a
global-equilibrium claim, an ideal-activity fallback, or a second equilibrium
owner. The ordinary deterministic target-model basin search remains mandatory
and runs before continuation.

## Typed endpoints and reaction references

The caller supplies an installed initial `ProviderPhase`. Both endpoint models
must expose the same ordered components and compatible public native-SDK
tensors, temperature domain, and overlapping admitted packing interval. The
receipt retains both immutable model fingerprints.

When the problem owns a source standard state, its one source reaction contract
is independently transformed through each endpoint Provider neutral reference.
When constants are already in Provider coordinates and endpoint fingerprints
differ, the initial endpoint requires its own complete, ordered, dimensionless
Provider-basis constant records. Reusing target-basis constants at the initial
endpoint fails closed. Those records carry the initial Provider fingerprint;
missing, target-bound, or otherwise mismatched endpoint identity fails before
native evaluation.

## Homotopy and corrector

For `lambda` in `[0, 1]`, the residual Helmholtz value, state gradient, state
Hessian, pressure, packing value, packing gradient, packing Hessian, and
compiled reaction reference vector are linearly interpolated between endpoint
evaluations at the same physical `(T, n, V)`. The target Provider inverse-
packing coordinate is only a diffeomorphic volume coordinate; every trial is
evaluated against both endpoint domains.

`lambda=0` must first be an independently certified local equilibrium of the
initial model. The deterministic corrector begins with step `1/4`, doubles an
accepted refined step up to `1/4`, and halves a rejected step. Refinement stops
below `1/256`; at most 32 path attempts are evaluated. These constants are
model-independent parts of this bounded route, not caller controls.

Every accepted continuation step must independently pass the existing balance,
charge, pressure, physical-KKT, reaction-affinity, strict-interior,
Provider-domain, derivative, rank, conditioning, and reduced-Hessian gates.
The merit and acceptance contract therefore cannot terminate on chart
stationarity alone. `lambda=1` is evaluated and certified with the exact target
Provider and target reaction reference; interpolation supplies no endpoint
certificate.

## Search and receipt

The fixed ordinary primary-start ceiling remains 25. Generated, duplicate,
infeasible, evaluated, and Provider-domain-rejected starts are reported
separately. Each continuation corrector and any negative-curvature recovery it
launches remains an attempt record. The endpoint basin is deduplicated against
ordinary target basins, and selection is recomputed using the existing
`lowest_observed_certified_local_value` rule.

Failure is fail-closed. Receipts distinguish initial-model failure, incompatible
endpoint domain or reference contracts, callback or derivative failure,
first-order criterion failure, strict-interior contact, refinement exhaustion,
and attempt-budget exhaustion, and retain the maximum certified `lambda`.
Failures during endpoint preflight retain the initial fingerprint and use a
typed physical-domain, derivative-inconsistency, or unsupported-derivative
blocker rather than a generic native exception.
Continuation through a fold, boundary equilibrium adjudication, phase discovery,
and a globality proof remain outside this contract.
