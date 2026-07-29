from __future__ import annotations

import math

import epcsaft_equilibrium


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
    balances = spec["balance_matrix"]
    feed = spec["feed_amounts"]
    spec["conserved_totals"] = tuple(
        math.fsum(
            row[index] * feed[index]  # type: ignore[index]
            for index in range(len(feed))  # type: ignore[arg-type]
        )
        for row in balances  # type: ignore[union-attr]
    )
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


def typed_problem(
    spec: dict[str, object],
    *,
    minimum_amount_mol: float = 1.0e-12,
    source_standard_state: epcsaft_equilibrium.ChemicalStandardState | None = None,
) -> epcsaft_equilibrium.ChemicalEquilibriumProblem:
    records = spec["equilibrium_constant_records"]
    ln_k = spec["ln_k"]
    balances = spec["balance_matrix"]
    feed = spec["feed_amounts"]
    return epcsaft_equilibrium.ChemicalEquilibriumProblem(
        species_ids=tuple(spec["species_ids"]),  # type: ignore[arg-type]
        charges=tuple(spec["charges"]),  # type: ignore[arg-type]
        molar_masses_kg_per_mol=tuple(spec["molar_masses_kg_per_mol"]),  # type: ignore[arg-type]
        balance_matrix=tuple(tuple(row) for row in balances),  # type: ignore[union-attr]
        conserved_totals=tuple(
            spec.get(
                "conserved_totals",
                tuple(
                    math.fsum(
                        row[index] * feed[index]  # type: ignore[index]
                        for index in range(len(feed))  # type: ignore[arg-type]
                    )
                    for row in balances  # type: ignore[union-attr]
                ),
            )
        ),
        reaction_matrix=tuple(tuple(row) for row in spec["reaction_matrix"]),  # type: ignore[union-attr]
        feed_amounts_mol=tuple(spec["feed_amounts"]),  # type: ignore[arg-type]
        equilibrium_constants=tuple(
            epcsaft_equilibrium.ChemicalEquilibriumConstant(
                ln_value=ln_k[index],  # type: ignore[index]
                source_id=record["source_id"],
                reference_id=record["reference_id"],
                reaction_orientation=record["reaction_orientation"],
                conversion_id=record["conversion_id"],
                dimensionless=record["dimensionless"],
            )
            for index, record in enumerate(records)  # type: ignore[arg-type]
        ),
        strict_interior_amount_floor_mol=minimum_amount_mol,
        source_standard_state=source_standard_state,
    )
