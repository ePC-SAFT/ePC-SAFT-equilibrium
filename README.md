# ePC-SAFT Equilibrium

`epcsaft-equilibrium` owns bounded equilibrium formulations over an installed
`epcsaft` provider. The first local candidate solves a pure-component
saturation boundary for the provider-approved methane, ethane, or propane
model. It is non-authoritative until a manager accepts its pending promotion
receipt.

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

## Native boundary

The extension calls `epcsaft.native_sdk(model)` and retains one model-bound
`epcsaft.native_sdk.v1` capsule while Ipopt evaluates both phase contexts. It
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
`docs/designs/2026-07-17-pure-saturation-slice.md`. Independent review approved
the candidate for manager review, but no promotion has been accepted, so the
lab retains runtime authority.
