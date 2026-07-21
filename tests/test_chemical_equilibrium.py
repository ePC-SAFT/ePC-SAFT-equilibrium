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


def _amount_chart(
    charges: tuple[int, ...], coordinates: tuple[float, ...], trace_floor: float = 1.0e-12
) -> dict[str, object]:
    return _equilibrium._chemical_amount_chart(charges, coordinates, trace_floor)


def test_amount_chart_maps_neutral_logs_and_general_ionic_shares_exactly() -> None:
    neutral = _amount_chart((0, 0), (math.log(2.0), math.log(3.0)))
    assert neutral["amounts"] == pytest.approx((2.0, 3.0), abs=2.0e-15)
    assert neutral["charge_residual"] == 0.0

    charges = (0, 2, 1, -1, -2, 0)
    coordinates = (math.log(2.0), math.log(1.0 / 3.0), 0.0, math.log(3.0), math.log(4.0))
    ionic = _amount_chart(charges, coordinates)
    assert ionic["amounts"] == pytest.approx((3.0, 0.25, 1.5, 1.0, 0.5, 4.0), abs=2.0e-15)
    assert ionic["coordinate_count"] == len(charges) - 1
    assert ionic["charge_residual"] == pytest.approx(0.0, abs=2.0e-15)
    assert ionic["trace_status"] == "interior"

    inverse = _equilibrium._chemical_amount_chart_inverse(charges, ionic["amounts"])
    assert inverse == pytest.approx(coordinates, abs=3.0e-15)


def test_amount_chart_is_permutation_and_amount_scale_equivariant() -> None:
    charges = (0, 2, 1, -1, -2, 0)
    amounts = (3.0, 0.25, 1.5, 1.0, 0.5, 4.0)
    permutation = (5, 3, 1, 0, 4, 2)
    permuted_charges = tuple(charges[index] for index in permutation)
    permuted_amounts = tuple(amounts[index] for index in permutation)

    coordinates = _equilibrium._chemical_amount_chart_inverse(permuted_charges, permuted_amounts)
    result = _amount_chart(permuted_charges, tuple(coordinates))
    assert result["amounts"] == pytest.approx(permuted_amounts, abs=3.0e-15)

    scale = 7.0
    scaled = tuple(scale * value for value in permuted_amounts)
    scaled_coordinates = _equilibrium._chemical_amount_chart_inverse(permuted_charges, scaled)
    scaled_result = _amount_chart(permuted_charges, tuple(scaled_coordinates))
    assert scaled_result["amounts"] == pytest.approx(scaled, rel=2.0e-15)


def test_amount_chart_has_exact_directional_first_and_second_derivatives() -> None:
    charges = (0, 2, 1, -1, -2, 0)
    center = (0.2, -0.7, 0.4, -0.3, 0.8)
    direction = (0.3, -0.2, 0.5, 0.1, -0.4)
    step = 2.0e-5

    lower = _amount_chart(
        charges,
        tuple(value - step * delta for value, delta in zip(center, direction, strict=True)),
    )
    result = _amount_chart(charges, center)
    upper = _amount_chart(
        charges,
        tuple(value + step * delta for value, delta in zip(center, direction, strict=True)),
    )

    dimension = len(center)
    jacobian = result["jacobian"]
    hessians = result["amount_hessians"]
    for species in range(len(charges)):
        first_directional = (upper["amounts"][species] - lower["amounts"][species]) / (
            2.0 * step
        )
        exact_first = sum(
            jacobian[species * dimension + column] * direction[column]
            for column in range(dimension)
        )
        assert first_directional == pytest.approx(exact_first, rel=3.0e-10, abs=3.0e-11)

        for row in range(dimension):
            jacobian_directional = (
                upper["jacobian"][species * dimension + row]
                - lower["jacobian"][species * dimension + row]
            ) / (2.0 * step)
            exact_second = sum(
                hessians[
                    species * dimension * dimension + row * dimension + column
                ]
                * direction[column]
                for column in range(dimension)
            )
            assert jacobian_directional == pytest.approx(
                exact_second, rel=2.0e-9, abs=5.0e-11
            )


def test_amount_chart_classifies_trace_floor_without_zeroing_species() -> None:
    trace_floor = 1.0e-10
    result = _amount_chart((0, 0), (math.log(0.5 * trace_floor), 0.0), trace_floor)

    assert result["trace_status"] == "at_or_below_floor"
    assert result["minimum_amount"] == pytest.approx(0.5 * trace_floor, rel=2.0e-15)
    assert all(value > 0.0 for value in result["amounts"])


@pytest.mark.parametrize(
    ("charges", "coordinates", "message"),
    (
        ((1, 0), (0.0,), "both cations and anions"),
        ((1, -1), (0.0, 0.0), "coordinate count"),
        ((1, -1), (math.inf,), "finite"),
    ),
)
def test_amount_chart_rejects_invalid_topology_or_coordinates(
    charges: tuple[int, ...], coordinates: tuple[float, ...], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        _amount_chart(charges, coordinates)
