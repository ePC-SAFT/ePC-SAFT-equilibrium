from __future__ import annotations

import csv
import hashlib
import json
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent
NATIVE = ROOT / "native"
LITERATURE_SOURCE = ROOT.parent / "al-sahhaf-1997-licl" / "al-sahhaf-1997-table3.csv"
MOLAR_MASS = (6.941, 35.453, 18.01528, 74.1216)


def mass_fractions(mole_fractions: list[float]) -> tuple[float, float, float]:
    masses = [x * molar_mass for x, molar_mass in zip(mole_fractions, MOLAR_MASS, strict=True)]
    total = sum(masses)
    return masses[2] / total, masses[3] / total, (masses[0] + masses[1]) / total


accepted: list[dict[str, object]] = []
outcomes: list[dict[str, object]] = []
for path in sorted(NATIVE.glob("figure3-*.json")):
    result = json.loads(path.read_text())
    outcomes.append(
        {
            "case_id": path.stem,
            "outcome": result["outcome"],
            "failure_stage": result["failure_stage"],
            "failure_reason": result["failure_reason"],
            "upper_solve_count": result["upper_solve_count"],
        }
    )
    if result["outcome"] != "physical_equilibrium_accepted":
        continue
    phases = [
        (phase, mass_fractions(phase["mole_fractions"])) for phase in result["phases"]
    ]
    phases.sort(key=lambda item: item[1][0], reverse=True)
    for label, (phase, fractions) in zip(("aqueous", "organic"), phases, strict=True):
        accepted.append(
            {
                "case_id": path.stem,
                "phase": label,
                "water_mass_fraction": fractions[0],
                "one_butanol_mass_fraction": fractions[1],
                "licl_mass_fraction": fractions[2],
                "phase_fraction": phase["phase_fraction"],
            }
        )

for filename, rows in (
    ("model-accepted-splits.csv", accepted),
    ("case-outcomes.csv", outcomes),
):
    with (ROOT / filename).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

literature_text = LITERATURE_SOURCE.read_text()
(ROOT / "al-sahhaf-1997-table3.csv").write_text(literature_text)
with (ROOT / "al-sahhaf-1997-table3.csv").open(newline="") as stream:
    literature = list(csv.DictReader(stream))

native_campaign = json.loads((ROOT / "native-campaign.json").read_text())
python_campaign = json.loads((ROOT / "python-campaign.json").read_text())
walls = [case["wall_seconds"] for case in native_campaign["cases"]]
python_walls = [case["wall_seconds"] for case in python_campaign["cases"]]
organic_salt = [
    row["licl_mass_fraction"] for row in accepted if row["phase"] == "organic"
]
summary = {
    "parameter_fingerprint": json.loads(next(NATIVE.glob("*.json")).read_text())[
        "parameter_fingerprint"
    ],
    "native_case_count": len(outcomes),
    "accepted_split_count": len(accepted) // 2,
    "one_phase_count": sum(row["outcome"].startswith("one_phase") for row in outcomes),
    "indeterminate_count": sum(row["outcome"] == "indeterminate" for row in outcomes),
    "native_wall_seconds": {
        "campaign": native_campaign["campaign_wall_seconds"],
        "minimum": min(walls),
        "median": statistics.median(walls),
        "mean": statistics.mean(walls),
        "maximum": max(walls),
    },
    "python_native_parity": python_campaign["all_parity"],
    "python_wall_seconds": {
        "campaign": python_campaign["campaign_wall_seconds"],
        "minimum": min(python_walls),
        "median": statistics.median(python_walls),
        "mean": statistics.mean(python_walls),
        "maximum": max(python_walls),
    },
    "model_organic_licl_mass_fraction_range": [min(organic_salt), max(organic_salt)],
    "al_sahhaf_organic_licl_mass_fraction_range": [
        min(float(row["organic_licl_mass_percent"]) / 100.0 for row in literature),
        max(float(row["organic_licl_mass_percent"]) / 100.0 for row in literature),
    ],
    "literature_source_sha256": hashlib.sha256(literature_text.encode()).hexdigest(),
    "comparison_basis": "Al-Sahhaf et al. (1997) Table 3 experimental tie lines only",
    "model_basis": (
        "Figiel et al. (2025) corrected-Born/dielectric model with published "
        "Mohammad et al. (2016) 1-butanol-ion k_ij and l_ij"
    ),
}
(ROOT / "comparison-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
