from __future__ import annotations

import csv
from pathlib import Path

import epcsaft

import epcsaft_equilibrium

REPO_ROOT = Path(__file__).resolve().parents[1]
ANCHORS = REPO_ROOT / "data" / "reference" / "pure_saturation_anchors.csv"
LAB_RELATIVE_TOLERANCE = 5.0e-6
NIST_PRESSURE_RELATIVE_TOLERANCE = 5.0e-3
NIST_LIQUID_DENSITY_RELATIVE_TOLERANCE = 1.0e-2


def relative_error(value: float, reference: float) -> float:
    return abs(value - reference) / abs(reference)


def main() -> int:
    with ANCHORS.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))

    for row in rows:
        component = row["component"]
        model = epcsaft.EPCSAFT(
            epcsaft.ParameterBundle.from_catalog(row["catalog"], version=1).select((component,))
        )
        result = epcsaft_equilibrium.saturation(
            model,
            float(row["T_K"]) * epcsaft.unit_registry.kelvin,
        )
        lab_pressure_error = relative_error(
            result.saturation_pressure_pa,
            float(row["lab_pressure_pa"]),
        )
        lab_vapor_density_error = relative_error(
            result.vapor.molar_density_mol_m3,
            float(row["lab_vapor_density_mol_m3"]),
        )
        lab_liquid_density_error = relative_error(
            result.liquid.molar_density_mol_m3,
            float(row["lab_liquid_density_mol_m3"]),
        )
        nist_pressure_error = relative_error(
            result.saturation_pressure_pa,
            float(row["nist_pressure_pa"]),
        )
        liquid_density_kg_m3 = result.liquid.molar_density_mol_m3 * float(row["molar_mass_kg_mol"])
        nist_liquid_density_error = relative_error(
            liquid_density_kg_m3,
            float(row["nist_liquid_density_kg_m3"]),
        )
        if (
            max(
                lab_pressure_error,
                lab_vapor_density_error,
                lab_liquid_density_error,
            )
            > LAB_RELATIVE_TOLERANCE
        ):
            raise RuntimeError(f"{component} does not reproduce the retained lab boundary")
        if nist_pressure_error > NIST_PRESSURE_RELATIVE_TOLERANCE:
            raise RuntimeError(f"{component} pressure exceeds the retained NIST tolerance")
        if nist_liquid_density_error > NIST_LIQUID_DENSITY_RELATIVE_TOLERANCE:
            raise RuntimeError(f"{component} liquid density exceeds the retained NIST tolerance")
        print(
            f"{component}: P={result.saturation_pressure_pa:.9g} Pa, "
            f"rho_v={result.vapor.molar_density_mol_m3:.9g} mol/m^3, "
            f"rho_l={result.liquid.molar_density_mol_m3:.9g} mol/m^3, "
            f"NIST pressure error={100.0 * nist_pressure_error:.4g}%"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
