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
from chemical_equilibrium_cases import (
    manufactured_solve as _manufactured_solve,
)

from epcsaft_equilibrium import _equilibrium


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
    assert compiled["independent_reaction_matrix"] == [
        [-1.0, 1.0, 0.0],
        [0.0, -1.0, 1.0],
    ]
    assert compiled["independent_ln_k"] == pytest.approx(
        (math.log(2.0), math.log(3.0)), abs=2.0e-14
    )
    assert compiled["reaction_cycle_inf_norm"] <= 2.0e-14
    assert compiled["reaction_transform_inf_norm"] <= 2.0e-14


def test_reaction_compiler_is_stable_for_scaled_nearly_dependent_reactions() -> None:
    delta = 1.0e-6
    expected_reference = (1.0, -1.0, 0.5, -0.5)
    qr_diagonal_ratios = []

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
        assert compiled["reference_reconstruction_inf_norm"] <= 5.0e-8
        qr_diagonal_ratios.append(compiled["reaction_qr_diagonal_ratio"])

    assert qr_diagonal_ratios == pytest.approx(
        (qr_diagonal_ratios[0],) * 4, rel=2.0e-9
    )
    assert 1.0e-8 < qr_diagonal_ratios[0] < 1.0e-4


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


@pytest.mark.parametrize(
    (
        "balance_matrix",
        "feed_amounts",
        "charges",
        "molar_masses",
        "classifications",
        "expected_average",
    ),
    (
        (
            ((1.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
            (1.0, 0.0, 0.0),
            (0, 0, 0),
            (1.0, 1.0, 1.0),
            ("proved_accessible", "proved_accessible", "proved_structural_zero"),
            (0.5, 0.5, 0.0),
        ),
        (
            ((1.0e8, 1.0e8, 0.0), (0.0, 0.0, 1.0e-8)),
            (1.0, 0.0, 0.0),
            (0, 0, 0),
            (1.0, 1.0, 1.0),
            ("proved_accessible", "proved_accessible", "proved_structural_zero"),
            (0.5, 0.5, 0.0),
        ),
        (
            ((2.0, 1.0, 1.0),),
            (1.0, 0.0, 0.0),
            (0, 1, -1),
            (2.0, 1.0, 1.0),
            ("proved_accessible", "proved_accessible", "proved_accessible"),
            (1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0),
        ),
    ),
)
def test_homogeneous_structural_support_has_exact_primal_and_dual_certificates(
    balance_matrix: tuple[tuple[float, ...], ...],
    feed_amounts: tuple[float, ...],
    charges: tuple[int, ...],
    molar_masses: tuple[float, ...],
    classifications: tuple[str, ...],
    expected_average: tuple[float, ...],
) -> None:
    evidence = _equilibrium._chemical_analyze_homogeneous_support(
        balance_matrix,
        feed_amounts,
        charges,
        molar_masses,
    )

    assert evidence["phase1_status"] == "optimal"
    assert evidence["validation_status"] == "exact_certificates_complete"
    assert tuple(item["classification"] for item in evidence["species"]) \
        == classifications
    assert evidence["equality_inf_norm"] == 0.0

    totals = tuple(
        sum(row[index] * feed_amounts[index] for index in range(len(feed_amounts)))
        for row in balance_matrix
    )
    for item in evidence["species"]:
        if item["primal_validated"]:
            witness = item["witness_amounts"]
            assert all(value >= 0.0 for value in witness)
            assert tuple(
                sum(row[index] * witness[index] for index in range(len(witness)))
                for row in balance_matrix
            ) == pytest.approx(totals, abs=0.0)
            assert sum(
                charges[index] * witness[index] for index in range(len(witness))
            ) == 0.0
        if item["classification"] == "proved_structural_zero":
            assert item["dual_validated"] is True
            assert item["candidate_maximum_mass_fraction"] == pytest.approx(
                0.0, abs=1.0e-12
            )

    average = evidence["witness_average_amounts"]
    assert average == pytest.approx(expected_average, abs=0.0)
    assert tuple(
        sum(row[index] * average[index] for index in range(len(average)))
        for row in balance_matrix
    ) == pytest.approx(totals, abs=0.0)
    assert sum(charges[index] * average[index] for index in range(len(average))) == 0.0
    for index, classification in enumerate(classifications):
        if classification == "proved_accessible":
            assert average[index] > 0.0


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


def test_accessible_face_preserves_reaction_combinations_that_cancel_removed_species() -> None:
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
    assert compiled["accessible_reaction_matrix"] == [[-1.0, 1.0]]
    assert compiled["accessible_ln_k"] == pytest.approx(
        (math.log(6.0),), abs=2.0e-14
    )
    assert compiled["balance_rank"] + compiled["reaction_rank"] == 2
    assert compiled["conservation_reaction_inf_norm"] <= 2.0e-14
    assert compiled["charge_reaction_inf_norm"] == 0.0


def test_manufactured_structural_face_has_exact_zeros_and_local_certification_level() -> None:
    result = _manufactured_solve(_accessible_face_system())

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


def test_provider_boundary_direction_guard_fails_before_any_callback() -> None:
    result = _equilibrium._chemical_provider_boundary_guard(
        _accessible_face_system()
    )

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
