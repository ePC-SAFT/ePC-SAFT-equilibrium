# Pure-Component Saturation Slice Design

Stage: Phase 6 equilibrium promotion candidate

Approval: The user approved this exact bounded slice in delegation from task
`019f5dba-f699-79c2-94ca-93fdf49d9b4d`. That approval authorizes local
implementation and verification. It does not transfer runtime authority or
accept a promotion receipt.

## Capability

The package will expose one operation:

```python
saturation(model: epcsaft.EPCSAFT, temperature: pint.Quantity) -> SaturationResult
```

It will accept the provider-approved pure methane, ethane, and propane
parameter fingerprints at temperatures inside their Gross and Sadowski (2001)
source domains. A source-domain input does not guarantee a two-phase result.
The route will reject states where it cannot certify a distinct local vapor and
liquid boundary.

The result will contain two immutable phase states, the common pressure,
separate solver, numerical, and physical status fields, residual diagnostics,
the provider fingerprint, and an explicit statement that no globality
certificate exists. Public outputs use named SI fields.

## Ownership and dependency boundary

The installed `epcsaft` wheel owns the fixed-temperature phase function

```text
Phi(T,n,V) = A(T,n,V)/(R T n_ref)
```

and its CppAD gradient, Hessian, and third tensor in `(n_mol, V_m3)`. The
equilibrium extension obtains one model-bound `epcsaft.native_sdk.v1` capsule
through the public `epcsaft.native_sdk(model)` function. Two phase contexts
hold the same capsule table and immutable model context. The extension neither
links provider symbols nor compiles provider source.

The equilibrium extension owns phase identity, logarithmic state variables,
the local coexistence constraints, Ipopt callbacks, initialization, exact
constraint derivatives, diagnostics, certification, and result serialization.
It will compile against the declaration header installed by the provider
wheel.

## Formulation

Each phase uses `n = 1 mol`, `V = 1/rho`, and the coordinates

```text
q_v = ln(rho_v / (1 mol m^-3))
q_l = ln(rho_l / (1 mol m^-3))
r   = ln(P_sat / 1 Pa).
```

The zero-objective Ipopt feasibility problem has three equality constraints:

```text
c_v  = P_v/P_sat - 1
c_l  = P_l/P_sat - 1
c_mu = mu_v/(RT) - mu_l/(RT).
```

The pressure ratios keep low and high saturation pressures on the same scale.
Positive density and pressure follow from the logarithmic coordinates.
Separated vapor and liquid density bounds prevent the optimizer from treating
the identical-state root as a two-phase solution. Postsolve certification also
rejects bound-active or collapsed phase states.

For `q = ln(rho)` and `V = exp(-q)`, the extension transforms the provider
tensors by exact chain rules. It derives `P_q`, `P_qq`, `mu_q`, and `mu_qq`
from the provider Hessian and third tensor, then supplies the exact Jacobian and
lower-triangular Lagrangian Hessian to Ipopt. The production route has no
finite-difference or alternate derivative mode.

## Initialization and numerical acceptance

The route uses a fixed, component-independent grid of vapor and liquid density
seeds. Each pressure seed starts from `R T rho_v`. It stops discovery after the
first solver-converged, physically acceptable local boundary, then reruns two
perturbed seeds around that boundary. Numerical convergence requires the
confirmation solves to agree within the declared pressure and density
tolerances.

The result separates:

- solver convergence: accepted Ipopt termination with finite diagnostics;
- numerical convergence: agreement of perturbed local solves; and
- physical acceptance: pressure and chemical-potential equality, positive and
  distinct phase densities and volumes, bounds, finite values, matching
  fingerprints, and exact derivative evidence.

The public operation raises `SaturationError` with structured diagnostics when
any layer fails. It raises `ValueError` for malformed or out-of-domain inputs.

## Alternatives considered

1. A two-variable Newton or least-squares root solve would use fewer solver
   callbacks. It would omit the required Ipopt integration and exact
   Lagrangian-Hessian contract.
2. Direct total-free-energy minimization fits later specified-`T,P,N` flash and
   phase-discovery work. For an unknown pure saturation pressure it requires an
   artificial overall volume or split constraint.
3. The selected three-variable Ipopt feasibility formulation matches the
   current provider ADR, represents the unknown boundary directly, and keeps
   all nonlinear thermodynamics in the provider.

## Evidence and exclusions

Package tests will check the capsule ABI, sizes, lifetime, fingerprints,
rejection paths, exact directional derivatives, Ipopt callbacks, local
certification, and representative methane, ethane, and propane results. The
scientific comparison uses retained NIST WebBook saturation rows and retained
lab route values as external evidence. Expected values will not come from the
new implementation.

The slice excludes mixture equilibrium, association, electrolytes, reaction,
density-root sensitivities, bubble and dew routes, flash, phase discovery,
TPD or global stability, critical completion, continuation, regression,
provider equations, release machinery, and remote changes.

