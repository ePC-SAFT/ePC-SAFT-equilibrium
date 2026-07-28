from __future__ import annotations

import math


def base_system() -> dict[str, object]:
    temperature_k = 350.0
    pressure_pa = 200_000.0
    return {
        "species_ids": ("A", "B"),
        "charges": (0, 0),
        "provider_fingerprint": "sha256:manufactured",
        "molar_masses_kg_per_mol": (1.0, 1.0),
        "balance_matrix": ((1.0, 1.0),),
        "reaction_matrix": ((-1.0, 1.0),),
        "feed_amounts": (1.0, 0.0),
        "ln_k": (math.log(4.0),),
        "equilibrium_constant_records": (
            {
                "source_id": "manufactured:A-to-B",
                "reference_id": "provider-helmholtz-coordinate-basis",
                "reaction_orientation": "products_positive",
                "conversion_id": "already-provider-basis",
                "dimensionless": True,
                "temperature_k": temperature_k,
                "pressure_pa": pressure_pa,
            },
        ),
        "temperature_k": temperature_k,
        "pressure_pa": pressure_pa,
    }


def bind_records(spec: dict[str, object]) -> None:
    spec["equilibrium_constant_records"] = tuple(
        {
            "source_id": f"manufactured:reaction-{index}",
            "reference_id": "provider-helmholtz-coordinate-basis",
            "reaction_orientation": "products_positive",
            "conversion_id": "already-provider-basis",
            "dimensionless": True,
            "temperature_k": spec["temperature_k"],
            "pressure_pa": spec["pressure_pa"],
        }
        for index in range(len(spec["ln_k"]))  # type: ignore[arg-type]
    )
