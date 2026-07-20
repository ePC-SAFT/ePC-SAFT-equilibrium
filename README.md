# ePC-SAFT Equilibrium

`epcsaft-equilibrium` owns bounded equilibrium formulations over an installed
`epcsaft` provider. The accepted first capability solves a pure-component
saturation boundary for the provider-approved methane, ethane, or propane
model under promotion receipt `promotion-0018-equilibrium-pure-saturation-v1`.
The package also contains a non-authoritative local candidate for one bounded
neutral methane/ethane HELD `T,P,z` calculation.

```python
import epcsaft
from epcsaft_equilibrium import saturation

parameters = epcsaft.ParameterBundle.from_catalog(
    "gross-2001-methane-ethane", version=1
).select(("methane",))
model = epcsaft.EPCSAFT(parameters)

result = saturation(model, 150 * epcsaft.unit_registry.kelvin)
print(result.saturation_pressure_pa)
print(result.vapor.molar_density_mol_m3, result.liquid.molar_density_mol_m3)
```

The operation requires a Pint temperature quantity and accepts these source
domains: methane 97 to 300 K, ethane 90 to 305 K, and propane 85 to 523 K. A
temperature inside a parameter domain does not guarantee a two-phase boundary.
The operation raises `SaturationError` when Ipopt cannot produce a numerically
confirmed state that passes the local physical checks.

The result reports solver convergence, numerical confirmation, and physical
acceptance as separate fields. It checks pressure and chemical-potential
equality, positive distinct phase states, local mechanical stability, bounds,
finite values, and the provider fingerprint. It carries no globality or phase
discovery certificate.

## Local neutral HELD candidate

```python
import epcsaft
from epcsaft_equilibrium import tp_flash

parameters = epcsaft.ParameterBundle.from_catalog(
    "gross-2001-methane-ethane", version=1
).select(("methane", "ethane"))
model = epcsaft.EPCSAFT(parameters)

result = tp_flash(
    model,
    243.61 * epcsaft.unit_registry.kelvin,
    6.691 * epcsaft.unit_registry.megapascal,
    (0.5627, 0.4373),
)
print(result.phases, result.phase_fractions)
print(result.diagnostics.outcome, result.diagnostics.globality_certificate)
```

The controller admits only the reviewed binary fingerprint inside the
rectangular May et al. (2015), Table 5 source domain: 203.22--243.61 K,
2.124--6.885 MPa, and methane feed fraction 0.4661--0.66705. It returns one or
two phases after the declared finite Stage-I/II/III search. `FlashError`
retains invalid-input, provider, exhausted-search, scope, and indeterminate
diagnostics. Every result reports `globality_certificate="not_guaranteed"`;
the accepted pure-saturation authority above is unchanged. HELD diagnostics
classify solver, numerical, and physical evidence independently as `passed`,
`failed`, or `not_adjudicated`; those axes do not change the globality claim.

The candidate fingerprint is
`sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286`.
The retained Pereira source is the permanent-lab Markdown at commit
`13ce345b6dcc41d399bb2a4c7b9bedb18f74b45b`, blob
`dde7f02d4c93cce86804a8e6b62d37602990ac21`; it is provenance, not a build,
test, validation, or runtime dependency.

The earlier fixed-two-phase campaign remains `NON_ADMISSION` under the frozen
May `3*u_c` composition contract: 12 of 17 rows passed, rows 002/009/010/011
were solved model/data misses, and row 012 had no package-accepted local state.
HELD does not erase or reclassify those results. The later installed-artifact
HELD campaign is also stable `NON_ADMISSION`: it returned two one-phase
results, thirteen fail-closed third-candidate results, and three exhausted
searches, so no two-phase composition comparison was evaluable. Permanent-lab
review requires a focused controller-lifecycle redesign before another runtime
candidate. The current wheel remains immutable and unpromoted; the design-only
delta changes no runtime, tolerance, resource, public type, or globality claim.
That neutral redesign is retained as deferred provenance while the separately
authorized Perdomo HELD2 implementation is the active development path.

## Native boundary

The extension calls `epcsaft.native_sdk(model)` and retains one model-bound
`epcsaft.native_sdk.v1` capsule while Ipopt evaluates the phase contexts. The
HELD candidate consumes the reviewed mixture value/gradient/Hessian tail;
the pure route continues to consume the accepted prefix. The extension
compiles against the declaration header installed by the provider wheel. It
does not link provider implementation symbols, compile provider sources, or
import private provider modules.

Source builds require Python 3.13, CMake, a C++17 compiler, pkg-config, Ipopt,
and the non-editable provider wheel installed in the build environment. The
local candidate gate hashes the exact provider wheel before creating an
isolated build environment. Candidate wheels are retained as read-only files
under a commit-bound `artifacts/equilibrium-neutral-held-v1/<commit>/`
directory and are never overwritten by a correction. The exact candidate
commit, declared `SOURCE_DATE_EPOCH`, artifact hashes, review status, and
replayable commands are recorded in
`receipts/promotion/promotion-0002-neutral-held-v1-candidate.yaml`.
One clean build does not support an archive-byte reproducibility claim.
The current corrected subject is implementation commit `8318e755`, retained
under `artifacts/equilibrium-neutral-held-v1/8318e75/`; the rejected
`549162a3` artifact remains immutable provenance.

```text
uv run --isolated --no-project --python 3.13 \
  --with "$PINNED_PROVIDER_WHEEL" --with scikit-build-core --with pybind11 \
  -- uv build --no-build-isolation --wheel .
```

Run the compact package proof and scientific anchors with:

```text
pytest -q
python scripts/validate_saturation.py
```

The design and equations are recorded in
`docs/designs/2026-07-17-pure-saturation-slice.md` and
`docs/designs/2026-07-17-neutral-held-v1.md`. Migration receipt
`promotion-0018-equilibrium-pure-saturation-v1` makes this repository the
production owner of that exact local boundary capability. One local boundary
solve is not a phase-discovery or global-stability proof. The local HELD
candidate has authority effect `none` until separate review, validation,
provider-tail promotion, equilibrium promotion, and explicit user approval.
