# Pure Saturation Anchor Provenance

`pure_saturation_anchors.csv` is a three-row validation extract, not a fitted
dataset and not an authority book.

## Source identity

- Transitional source repository: `tannerpolley/ePC-SAFT-lab`
- Source commit: `13ce345b6dcc41d399bb2a4c7b9bedb18f74b45b`
- Retained route table:
  `analyses/package_validation/equilibrium_single_component_vle/figures/hydrocarbon_saturation_pressure/results/hydrocarbon_saturation_pressure.csv`
- Route-table SHA-256:
  `d8c6f764a8eafeba908c2056936953de0be8ed003079022488fe6e33f97f51dc`
- Retrieved for this distillation: 2026-07-17

The NIST rows at the same source commit are:

| Component | Repository-relative path | SHA-256 |
| --- | --- | --- |
| methane | `data/reference/pure_component/saturation_properties/methane/saturation_properties.csv` | `a5e16df3bf8ec78483fc340782cddc89ab8b284a9f6dfaecd6cda3ffde579227` |
| ethane | `data/reference/pure_component/saturation_properties/ethane/saturation_properties.csv` | `ed09b8781acfb7025ca505878b884f6353ddd9f3f4bd7aae2e6df88bbe847a67` |
| propane | `data/reference/pure_component/saturation_properties/propane/saturation_properties.csv` | `1cfc73d68061379581cc545d4ac03beb44762aaf162b7572752228cd62481495` |

Each CSV row retains its direct NIST Chemistry WebBook query URL. The source
lab and this repository are both GPL-3.0 licensed. This extract retains three
factual thermophysical observations with citations and three lab route outputs
under that same project license; it does not copy WebBook presentation or
provider equations.

## Transformation record

Before the new implementation was run, one accepted row was selected for each
admitted component: methane 150 K, ethane 240 K, and propane 300 K. Component
labels were lowercased. Temperature, NIST pressure, NIST saturated-liquid mass
density, route pressure, route vapor molar density, route liquid molar density,
and NIST URL were copied without interpolation or recalculation. Catalog names
and fingerprints were attached from the pinned provider artifact. Molar mass
was retained from the route row's mass-density to molar-density conversion and
is used only to compare the new molar-density result with the NIST mass-density
observation.

## Gate rationale

- The `5e-6` relative lab tolerance checks behavioral preservation across a
  rewritten local formulation while allowing solver termination-level changes.
- The `0.5%` NIST pressure and `1%` NIST liquid-density tolerances are explicit
  scientific comparison gates for these three selected observations. They are
  looser than the lab-reproduction gate and cover the retained source-route
  errors at the selected rows; they do not claim universal model accuracy.

The executable owner of these thresholds is `scripts/validate_saturation.py`;
the identical capability-level values are recorded in `ARCHITECTURE.yaml`.
