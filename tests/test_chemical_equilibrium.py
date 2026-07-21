from __future__ import annotations

import copy
import math
from collections.abc import Callable

import pytest

from epcsaft_equilibrium import _equilibrium


def _base_system() -> dict[str, object]:
    temperature_k = 350.0
    pressure_pa = 200_000.0
    return {
        "species_ids": ("A", "B"),
        "charges": (0, 0),
        "provider_component_ids": ("A", "B"),
        "provider_charges": (0, 0),
        "provider_fingerprint": "sha256:manufactured",
        "expected_provider_component_ids": ("A", "B"),
        "expected_provider_charges": (0, 0),
        "expected_provider_fingerprint": "sha256:manufactured",
        "balance_matrix": ((1.0, 1.0),),
        "declared_balance_rank": 1,
        "reaction_matrix": ((-1.0, 1.0),),
        "feed_amounts": (1.0, 0.0),
        "ln_k": (math.log(4.0),),
        "equilibrium_constant_records": (
            {
                "source_id": "manufactured:A-to-B",
                "reference_id": "provider-helmholtz-coordinate-basis",
                "dimensionless": True,
                "temperature_k": temperature_k,
                "pressure_pa": pressure_pa,
                "pressure_binding": "fixed",
            },
        ),
        "temperature_k": temperature_k,
        "pressure_pa": pressure_pa,
        "complete_closed_system": True,
    }


def test_reaction_compiler_reconstructs_minimum_norm_reference() -> None:
    compiled = _equilibrium._chemical_compile_system(_base_system())

    expected = 0.5 * math.log(4.0)
    assert compiled["species_ids"] == ["A", "B"]
    assert compiled["balance_rank"] == 1
    assert compiled["reaction_rank"] == 1
    assert compiled["g_ref"] == pytest.approx((expected, -expected), abs=2.0e-15)
    assert compiled["reference_reconstruction_inf_norm"] <= 2.0e-15
    assert compiled["conservation_reaction_inf_norm"] == 0.0
    assert compiled["charge_reaction_inf_norm"] == 0.0


def _change(path: str, value: object) -> Callable[[dict[str, object]], None]:
    def apply(spec: dict[str, object]) -> None:
        if "." not in path:
            spec[path] = value
            return
        outer, inner = path.split(".", maxsplit=1)
        records = [dict(record) for record in spec[outer]]  # type: ignore[arg-type]
        records[0][inner] = value
        spec[outer] = tuple(records)

    return apply


@pytest.mark.parametrize(
    ("mutate", "message"),
    (
        (_change("provider_component_ids", ("B", "A")), "component order"),
        (_change("provider_charges", (1, -1)), "charges"),
        (_change("provider_fingerprint", "sha256:wrong"), "fingerprint"),
        (_change("declared_balance_rank", 2), "declared balance rank"),
        (_change("balance_matrix", ((1.0, 1.0), (2.0, 2.0))), "balance matrix rank"),
        (_change("reaction_matrix", ((-1.0, 2.0),)), "conserve"),
        (
            _change("reaction_matrix", ((-1.0, 1.0), (-2.0, 2.0))),
            "reaction matrix rank",
        ),
        (_change("complete_closed_system", False), "complete closed system"),
        (_change("equilibrium_constant_records.dimensionless", False), "dimensionless"),
        (_change("equilibrium_constant_records.source_id", ""), "source identity"),
        (_change("equilibrium_constant_records.temperature_k", 351.0), "temperature"),
        (_change("equilibrium_constant_records.pressure_pa", 201_000.0), "pressure"),
    ),
)
def test_reaction_compiler_rejects_inconsistent_contracts(
    mutate: Callable[[dict[str, object]], None], message: str
) -> None:
    spec = copy.deepcopy(_base_system())
    if message == "complete closed system":
        spec["reaction_matrix"] = ()
        spec["ln_k"] = ()
        spec["equilibrium_constant_records"] = ()
        spec["complete_closed_system"] = True
    else:
        mutate(spec)

    with pytest.raises(ValueError, match=message):
        _equilibrium._chemical_compile_system(spec)


def test_reaction_compiler_rejects_non_neutral_feed_and_charge_nonconservation() -> None:
    non_neutral = _base_system()
    non_neutral.update(
        {
            "charges": (1, -1),
            "provider_charges": (1, -1),
            "expected_provider_charges": (1, -1),
            "feed_amounts": (1.0, 0.0),
        }
    )
    with pytest.raises(ValueError, match="electroneutral"):
        _equilibrium._chemical_compile_system(non_neutral)

    nonconserving = copy.deepcopy(non_neutral)
    nonconserving["feed_amounts"] = (1.0, 1.0)
    with pytest.raises(ValueError, match="charge"):
        _equilibrium._chemical_compile_system(nonconserving)
