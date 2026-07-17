# ePC-SAFT Equilibrium

`epcsaft-equilibrium` owns bounded equilibrium formulations over an installed
`epcsaft` provider. The accepted first capability solves a pure-component
saturation boundary for the provider-approved methane, ethane, or propane
model under promotion receipt `promotion-0018-equilibrium-pure-saturation-v1`.
The package also contains a non-authoritative local candidate for one neutral
methane/ethane fixed-two-phase `T,P,z` calculation.

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

## Local methane/ethane flash candidate

```python
import epcsaft
from epcsaft_equilibrium import two_phase_flash

parameters = epcsaft.ParameterBundle.from_catalog(
    "gross-2001-methane-ethane", version=1
).select(("methane", "ethane"))
model = epcsaft.EPCSAFT(parameters)

result = two_phase_flash(
    model,
    243.58 * epcsaft.unit_registry.kelvin,
    3.949 * epcsaft.unit_registry.megapascal,
    (0.48815, 0.51185),
)
print(result.liquid.mole_fractions, result.vapor.mole_fractions)
print(result.diagnostics.kkt_stationarity_max_abs)
```

This route directly minimizes the total two-phase Helmholtz-plus-`PV` energy
with two linear material balances. It admits only the reviewed binary
fingerprint inside the rectangular May et al. (2015), Table 5 source domain:
203.22--243.61 K, 2.124--6.885 MPa, and methane feed fraction
0.4661--0.66705. Being inside that rectangle does not guarantee that a
distinct two-phase local minimum exists. `FlashError` reports rejected local
states with solver, confirmation, material-balance, stationarity, and retained
Ipopt multiplier diagnostics.

The candidate fixes the phase count at two. It performs no phase discovery,
TPD/global-stability analysis, continuation, or fallback solve, and it does
not claim that the returned local result is the globally stable state. The
accepted pure-saturation authority above is unchanged.

## Native boundary

The extension calls `epcsaft.native_sdk(model)` and retains one model-bound
`epcsaft.native_sdk.v1` capsule while Ipopt evaluates the phase contexts. The
flash candidate consumes the reviewed mixture value/gradient/Hessian tail;
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
`docs/designs/2026-07-17-neutral-two-phase-tp-flash.md`. Migration receipt
`promotion-0018-equilibrium-pure-saturation-v1` makes this repository the
production owner of that exact local boundary capability. One local boundary
solve is not a phase-discovery or global-stability proof.
