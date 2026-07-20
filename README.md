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
the accepted pure-saturation authority above is unchanged.

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
admission gate hashes the exact provider wheel before creating an isolated
build environment; the replayable commands and both wheel hashes are recorded
in the candidate receipt.

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
solve is not a phase-discovery or global-stability proof.
