from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
HEIGHT = math.sqrt(3.0) / 2.0


def xy(water: float, butanol: float, salt: float) -> tuple[float, float]:
    total = water + butanol + salt
    return (butanol + 0.5 * salt) / total, HEIGHT * salt / total


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


literature = rows(ROOT / "al-sahhaf-1997-table3.csv")
model = rows(ROOT / "model-accepted-splits.csv")
model_pairs = {
    case_id: [row for row in model if row["case_id"] == case_id]
    for case_id in dict.fromkeys(row["case_id"] for row in model)
}

plt.rcParams.update({"font.family": "serif", "font.size": 10})
fig, ax = plt.subplots(figsize=(8.2, 7.2), constrained_layout=True)
triangle = [(0.0, 0.0), (1.0, 0.0), (0.5, HEIGHT), (0.0, 0.0)]
ax.plot(*zip(*triangle, strict=True), color="0.12", linewidth=1.2)
for fraction in (0.2, 0.4, 0.6, 0.8):
    grid = "0.86"
    for endpoints in (
        (xy(fraction, 0.0, 1 - fraction), xy(fraction, 1 - fraction, 0.0)),
        (xy(0.0, fraction, 1 - fraction), xy(1 - fraction, fraction, 0.0)),
        (xy(0.0, 1 - fraction, fraction), xy(1 - fraction, 0.0, fraction)),
    ):
        ax.plot(*zip(*endpoints, strict=True), color=grid, linewidth=0.65, zorder=0)

for index, row in enumerate(literature):
    organic = xy(
        float(row["organic_water_mass_percent"]),
        float(row["organic_one_butanol_mass_percent"]),
        float(row["organic_licl_mass_percent"]),
    )
    aqueous = xy(
        float(row["aqueous_water_mass_percent"]),
        float(row["aqueous_one_butanol_mass_percent"]),
        float(row["aqueous_licl_mass_percent"]),
    )
    ax.plot(*zip(aqueous, organic, strict=True), color="0.68", linewidth=0.8, zorder=1)
    ax.scatter(*aqueous, facecolors="white", edgecolors="0.20", s=30, zorder=3,
               label="Al-Sahhaf et al. (1997), Table 3" if index == 0 else None)
    ax.scatter(*organic, facecolors="white", edgecolors="0.20", s=30, zorder=3)

for index, phases in enumerate(model_pairs.values()):
    points = [
        xy(
            float(phase["water_mass_fraction"]),
            float(phase["one_butanol_mass_fraction"]),
            float(phase["licl_mass_fraction"]),
        )
        for phase in phases
    ]
    ax.plot(*zip(*points, strict=True), color="#0072B2", linewidth=1.2, alpha=0.85, zorder=2)
    ax.scatter(*zip(*points, strict=True), color="#0072B2", marker="^", s=38, zorder=4,
               label="Current model, accepted native splits" if index == 0 else None)

detail = ax.inset_axes([0.10, 0.54, 0.31, 0.24])
literature_organic = [
    (
        float(row["organic_water_mass_percent"]) / 100.0,
        float(row["organic_licl_mass_percent"]) / 100.0,
    )
    for row in literature
]
model_organic = [
    (float(row["water_mass_fraction"]), float(row["licl_mass_fraction"]))
    for row in model
    if row["phase"] == "organic"
]
detail.scatter(*zip(*literature_organic, strict=True), facecolors="white", edgecolors="0.20", s=22)
detail.scatter(*zip(*model_organic, strict=True), color="#0072B2", marker="^", s=26)
detail.set(
    title="Organic-phase detail",
    xlabel=r"$w_{\mathrm{H_2O}}$",
    ylabel=r"$w_{\mathrm{LiCl}}$",
)
detail.grid(color="0.88", linewidth=0.6)
detail.spines[["top", "right"]].set_visible(False)

for fraction in (0.2, 0.4, 0.6, 0.8):
    ax.text(fraction, -0.032, f"{fraction:.1f}", ha="center", va="top", color="0.35")
ax.text(-0.02, -0.055, r"Water, $w_{\mathrm{H_2O}}$", ha="left", va="top", weight="bold")
ax.text(1.02, -0.055, r"1-butanol, $w_{\mathrm{BuOH}}$", ha="right", va="top", weight="bold")
ax.text(0.5, HEIGHT + 0.035, r"LiCl, $w_{\mathrm{LiCl}}$", ha="center", weight="bold")
ax.set_title(
    "Water + 1-Butanol + LiCl at 298.15 K\n"
    "Experimental tie lines and current predictive model"
)
ax.set_xlabel("Formula mass-fraction ternary coordinates")
ax.set_ylabel("LiCl formula mass fraction")
ax.legend(loc="upper right", frameon=False)
ax.set_aspect("equal")
ax.set_xlim(-0.06, 1.06)
ax.set_ylim(-0.08, HEIGHT + 0.08)
ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

for suffix in ("svg", "png", "pdf"):
    fig.savefig(ROOT / f"al-sahhaf-1997-current-model-overlay.{suffix}", dpi=240)
plt.close(fig)

summary_path = ROOT / "comparison-summary.json"
summary = json.loads(summary_path.read_text())
summary["plots"] = {
    suffix: {
        "path": str(path),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
    for suffix in ("svg", "png", "pdf")
    if (path := ROOT / f"al-sahhaf-1997-current-model-overlay.{suffix}").exists()
}
summary_path.write_text(json.dumps(summary, indent=2) + "\n")
