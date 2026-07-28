from __future__ import annotations

import copy
import math
from collections.abc import Callable

import pytest
from chemical_equilibrium_cases import base_system as _base_system
from chemical_equilibrium_cases import bind_records as _bind_record

from epcsaft_equilibrium import _equilibrium


def _solve(spec: dict[str, object]) -> dict[str, object]:
    return _equilibrium._chemical_solve_manufactured(spec, {})


@pytest.mark.parametrize("inconsistent", (False, True))
def test_redundant_reactions_preserve_equilibrium_or_reject_bad_cycles(
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
            _solve(spec)
        return

    result = _solve(spec)
    assert result["accepted"] is True
    assert result["amounts"] == pytest.approx((1.0 / 9.0, 2.0 / 9.0, 2.0 / 3.0))
    assert result["reaction_affinity_inf_norm"] <= 1.0e-7


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


def test_scaled_nearly_dependent_reactions_are_stable_until_rank_is_lost() -> None:
    expected = (
        2.0 / (1.0 + math.exp(2.0)),
        2.0 * math.exp(2.0) / (1.0 + math.exp(2.0)),
        2.0 / (1.0 + math.e),
        2.0 * math.e / (1.0 + math.e),
    )
    for scale, reverse_order in (
        (1.0, False),
        (1.0, True),
        (1.0e8, False),
        (1.0e8, True),
    ):
        result = _solve(
            _nearly_dependent_system(scale, 1.0e-6, reverse_order=reverse_order)
        )
        assert result["accepted"] is True
        assert result["amounts"] == pytest.approx(expected, rel=2.0e-8)

    with pytest.raises(ValueError, match="reaction matrix rank"):
        _solve(_nearly_dependent_system(1.0, 1.0e-14))


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
def test_reaction_system_rejects_inconsistent_contracts(
    mutate: Callable[[dict[str, object]], None], message: str
) -> None:
    spec = copy.deepcopy(_base_system())
    mutate(spec)
    with pytest.raises(ValueError, match=message):
        _solve(spec)


def test_reaction_system_rejects_charge_inconsistency() -> None:
    spec = _base_system()
    spec.update({"charges": (1, -1), "feed_amounts": (1.0, 0.0)})
    with pytest.raises(ValueError, match="electroneutral"):
        _solve(spec)

    spec["feed_amounts"] = (1.0, 1.0)
    with pytest.raises(ValueError, match="charge"):
        _solve(spec)


def test_accessible_face_preserves_reaction_combinations_and_exact_zeros() -> None:
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

    result = _solve(spec)

    assert result["accepted"] is True
    assert result["amounts"] == pytest.approx(
        (1.0 / 7.0, 0.0, 6.0 / 7.0, 0.0), rel=3.0e-8, abs=0.0
    )
    assert result["structural_zero_species_indices"] == [1, 3]
    assert result["boundary_status"] == "structural_face"
    assert result["chemical_certification_level"] == "LOCAL_EQUILIBRIUM"
