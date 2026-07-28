from __future__ import annotations

import copy
import math
from collections.abc import Callable

import pytest
from chemical_equilibrium_cases import (
    base_system as _base_system,
)
from chemical_equilibrium_cases import (
    bind_records as _bind_record,
)

from epcsaft_equilibrium import _equilibrium


def test_reaction_compiler_reconstructs_minimum_norm_reference() -> None:
    compiled = _equilibrium._chemical_compile_system(_base_system())

    expected = 0.5 * math.log(4.0)
    assert compiled["species_ids"] == ["A", "B"]
    assert compiled["balance_rank"] == 1
    assert compiled["reaction_rank"] == 1
    assert compiled["g_ref"] == pytest.approx((expected, -expected), abs=2.0e-15)
    assert sum(
        coefficient * reference
        for coefficient, reference in zip(
            _base_system()["reaction_matrix"][0],  # type: ignore[index]
            compiled["g_ref"],
            strict=True,
        )
    ) == pytest.approx(-math.log(4.0), abs=2.0e-15)


def _nearly_dependent_system(
    reaction_scale: float, delta: float, *, reverse_order: bool = False
) -> dict[str, object]:
    reactions: tuple[tuple[float, ...], ...] = (
        (-1.0, 1.0, 0.0, 0.0),
        (
            -reaction_scale,
            reaction_scale,
            reaction_scale * delta,
            -reaction_scale * delta,
        ),
    )
    ln_k: tuple[float, ...] = (2.0, reaction_scale * (2.0 - delta))
    if reverse_order:
        reactions = tuple(reversed(reactions))
        ln_k = tuple(reversed(ln_k))
    spec = {
        **_base_system(),
        "species_ids": ("A", "B", "C", "D"),
        "charges": (0, 0, 0, 0),
        "molar_masses_kg_per_mol": (1.0, 1.0, 1.0, 1.0),
        "balance_matrix": (
            (1.0, 1.0, 1.0, 1.0),
            (1.0, 1.0, 0.0, 0.0),
        ),
        "reaction_matrix": reactions,
        "feed_amounts": (1.0, 1.0, 1.0, 1.0),
        "ln_k": ln_k,
    }
    _bind_record(spec)
    return spec


@pytest.mark.parametrize("inconsistent", (False, True))
def test_redundant_reaction_compiler_and_reaction_constant_cycle(
    inconsistent: bool,
) -> None:
    spec = {
        **_base_system(),
        "species_ids": ("A", "B", "C"),
        "charges": (0, 0, 0),
        "molar_masses_kg_per_mol": (1.0, 1.0, 1.0),
        "balance_matrix": ((1.0, 1.0, 1.0),),
        "reaction_matrix": (
            (-1.0, 1.0, 0.0),
            (0.0, -1.0, 1.0),
            (-1.0, 0.0, 1.0),
        ),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (
            math.log(2.0),
            math.log(3.0),
            math.log(7.0) if inconsistent else math.log(6.0),
        ),
    }
    _bind_record(spec)

    if inconsistent:
        with pytest.raises(ValueError, match="reaction constant cycle"):
            _equilibrium._chemical_compile_system(spec)
        return

    compiled = _equilibrium._chemical_compile_system(spec)

    assert compiled["supplied_reaction_rank"] == 2
    assert compiled["reaction_basis_rows"] == [0, 1]
    for actual, expected in zip(
        compiled["reaction_transform"],
        ((1.0, 0.0), (0.0, 1.0), (1.0, 1.0)),
        strict=True,
    ):
        assert actual == pytest.approx(expected, abs=2.0e-14)
    assert compiled["reaction_matrix"] == [
        [-1.0, 1.0, 0.0],
        [0.0, -1.0, 1.0],
    ]
    assert compiled["ln_k"] == pytest.approx(
        (math.log(2.0), math.log(3.0)), abs=2.0e-14
    )
    for supplied, transform in zip(
        spec["reaction_matrix"], compiled["reaction_transform"], strict=True  # type: ignore[arg-type]
    ):
        reconstructed = tuple(
            sum(
                transform[basis] * compiled["reaction_matrix"][basis][species]
                for basis in range(2)
            )
            for species in range(3)
        )
        assert reconstructed == pytest.approx(supplied, abs=2.0e-14)
    reconstructed_ln_k = tuple(
        sum(transform[basis] * compiled["ln_k"][basis] for basis in range(2))
        for transform in compiled["reaction_transform"]
    )
    assert reconstructed_ln_k == pytest.approx(spec["ln_k"], abs=2.0e-14)


def test_reaction_compiler_is_stable_for_scaled_nearly_dependent_reactions() -> None:
    delta = 1.0e-6
    expected_reference = (1.0, -1.0, 0.5, -0.5)

    for reaction_scale, reverse_order in (
        (1.0, False),
        (1.0, True),
        (1.0e8, False),
        (1.0e8, True),
    ):
        compiled = _equilibrium._chemical_compile_system(
            _nearly_dependent_system(
                reaction_scale, delta, reverse_order=reverse_order
            )
        )

        assert compiled["reaction_rank"] == 2
        assert compiled["g_ref"] == pytest.approx(expected_reference, abs=2.0e-9)
        for row, ln_k in zip(
            compiled["reaction_matrix"], compiled["ln_k"], strict=True
        ):
            assert sum(
                coefficient * reference
                for coefficient, reference in zip(
                    row, compiled["g_ref"], strict=True
                )
            ) == pytest.approx(-ln_k, abs=5.0e-8)


def test_reaction_compiler_rejects_numerically_dependent_reactions() -> None:
    with pytest.raises(ValueError, match="reaction matrix rank"):
        _equilibrium._chemical_compile_system(_nearly_dependent_system(1.0, 1.0e-14))


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
        (_change("molar_masses_kg_per_mol", (1.0, 0.0)), "molar mass"),
        (_change("reaction_matrix", ((-1.0, 2.0),)), "conserve"),
        (_change("equilibrium_constant_records.dimensionless", False), "dimensionless"),
        (
            _change(
                "equilibrium_constant_records.reaction_orientation",
                "reactants_positive",
            ),
            "products_positive",
        ),
        (
            _change(
                "equilibrium_constant_records.conversion_id",
                "unconverted-source-basis",
            ),
            "already-provider-basis",
        ),
    ),
)
def test_reaction_compiler_rejects_inconsistent_contracts(
    mutate: Callable[[dict[str, object]], None], message: str
) -> None:
    spec = copy.deepcopy(_base_system())
    mutate(spec)

    with pytest.raises(ValueError, match=message):
        _equilibrium._chemical_compile_system(spec)


def test_reaction_compiler_rejects_non_neutral_feed_and_charge_nonconservation() -> None:
    non_neutral = _base_system()
    non_neutral.update(
        {
            "charges": (1, -1),
            "feed_amounts": (1.0, 0.0),
        }
    )
    with pytest.raises(ValueError, match="electroneutral"):
        _equilibrium._chemical_compile_system(non_neutral)

    nonconserving = copy.deepcopy(non_neutral)
    nonconserving["feed_amounts"] = (1.0, 1.0)
    with pytest.raises(ValueError, match="charge"):
        _equilibrium._chemical_compile_system(nonconserving)


def _accessible_face_system() -> dict[str, object]:
    spec = {
        **_base_system(),
        "species_ids": ("A", "B", "C", "X"),
        "charges": (0, 0, 0, 0),
        "molar_masses_kg_per_mol": (1.0, 2.0, 1.0, 1.0),
        "balance_matrix": (
            (1.0, 1.0, 1.0, 0.0),
            (0.0, 1.0, 0.0, 1.0),
        ),
        "reaction_matrix": (
            (-1.0, 1.0, 0.0, -1.0),
            (0.0, -1.0, 1.0, 1.0),
        ),
        "feed_amounts": (1.0, 0.0, 0.0, 0.0),
        "ln_k": (math.log(2.0), math.log(3.0)),
    }
    _bind_record(spec)
    return spec


def test_accessible_face_preserves_reaction_combinations_and_certification() -> None:
    spec = _accessible_face_system()

    compiled = _equilibrium._chemical_compile_system(spec)

    assert compiled["original_species_ids"] == ["A", "B", "C", "X"]
    assert compiled["species_ids"] == ["A", "C"]
    assert compiled["retained_species_indices"] == [0, 2]
    assert compiled["removed_species_indices"] == [1, 3]
    assert compiled["supplied_reaction_rank"] == 2
    assert compiled["support_classifications"] == [
        "proved_accessible",
        "proved_structural_zero",
        "proved_accessible",
        "proved_structural_zero",
    ]
    assert len(compiled["accessible_reaction_transform"]) == 1
    assert compiled["accessible_reaction_transform"][0] == pytest.approx(
        (1.0, 1.0), abs=2.0e-14
    )
    assert compiled["reaction_matrix"] == [[-1.0, 1.0]]
    assert compiled["ln_k"] == pytest.approx(
        (math.log(6.0),), abs=2.0e-14
    )
    assert compiled["balance_rank"] + compiled["reaction_rank"] == 2

    result = _equilibrium._chemical_solve_manufactured(spec, {})

    assert result["accepted"] is True
    assert result["amounts"] == pytest.approx(
        (1.0 / 7.0, 0.0, 6.0 / 7.0, 0.0), rel=3.0e-8, abs=0.0
    )
    assert result["amounts"][1] == 0.0
    assert result["amounts"][3] == 0.0
    assert result["retained_species_indices"] == [0, 2]
    assert result["structural_zero_species_indices"] == [1, 3]
    assert result["boundary_status"] == "structural_face"
    assert result["chemical_certification_level"] == "LOCAL_EQUILIBRIUM"
    assert result["support_qualifiers"] == []

    result = _equilibrium._chemical_provider_boundary_guard(spec)

    assert result["accepted"] is False
    assert result["amounts"] == []
    assert result["solver_status"] == "boundary_direction_unresolved"
    assert (
        result["chemical_certification_level"]
        == "BOUNDARY_DIRECTION_UNRESOLVED"
    )
    assert result["boundary_status"] == "boundary_direction_unresolved"
    assert result["structural_zero_species_indices"] == [1, 3]
    assert result["support_qualifiers"] == []
