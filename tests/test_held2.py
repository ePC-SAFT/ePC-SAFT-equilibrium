from __future__ import annotations

import math
from collections.abc import Sequence

import pytest

from epcsaft_equilibrium import _equilibrium

FORMULATION_ID = "perdomo-held2.modified-mole.manufactured.v1"
CHARGES = (0.0, 1.0, -1.0)
PHYSICAL_FEED = (0.5, 0.25, 0.25)
CHEMICAL_POTENTIALS = (3.0, -2.0, 4.0)


def _manufactured_helmholtz(composition: float, molar_volume: float) -> float:
    shifted = composition - 0.5
    inner_squared = 0.15**2
    outer_squared = 0.30**2
    composition_energy = 1.0e4 * (
        shifted**6 / 6.0
        - (inner_squared + outer_squared) * shifted**4 / 4.0
        + inner_squared * outer_squared * shifted**2 / 2.0
    )
    volume_energy = 0.5 * 5.0 * (molar_volume - 1.2) ** 2
    return composition_energy + volume_energy


def _manufactured_gradient(composition: float, molar_volume: float) -> tuple[float, float]:
    shifted = composition - 0.5
    inner_squared = 0.15**2
    outer_squared = 0.30**2
    return (
        1.0e4
        * (
            shifted**5
            - (inner_squared + outer_squared) * shifted**3
            + inner_squared * outer_squared * shifted
        ),
        5.0 * (molar_volume - 1.2),
    )


def _objective_oracle(vector: Sequence[float]) -> float:
    u_alpha, u_beta, volume_alpha, volume_beta = vector
    fraction = (u_beta - 0.5) / (u_beta - u_alpha)
    gibbs_alpha = _manufactured_helmholtz(u_alpha, volume_alpha) + volume_alpha
    gibbs_beta = _manufactured_helmholtz(u_beta, volume_beta) + volume_beta
    return fraction * gibbs_alpha + (1.0 - fraction) * gibbs_beta


def test_held2_perdomo_transform_lift_and_bounds_match_the_source_coordinates() -> None:
    result = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        (0.27, 0.73, 0.92, 1.08),
        CHEMICAL_POTENTIALS,
    )

    assert result["formulation_id"] == FORMULATION_ID
    assert result["globality_certificate"] == "not_guaranteed"
    assert result["eliminated_index"] == 2
    assert result["dependent_index"] == 0
    assert result["retained_indices"] == [0, 1]
    assert result["independent_indices"] == [1]
    assert result["modified_factors"] == pytest.approx([1.0, 2.0])
    assert result["modified_feed"] == pytest.approx([0.5, 0.5])
    assert result["independent_lower_bounds"] == pytest.approx([2.0e-10])
    assert result["independent_upper_bounds"] == pytest.approx([1.0])
    assert result["phase_fraction"] == pytest.approx(0.5)
    assert result["physical_phases"][0] == pytest.approx([0.73, 0.135, 0.135])
    assert result["physical_phases"][1] == pytest.approx([0.27, 0.365, 0.365])
    assert result["modified_phases"][0] == pytest.approx([0.73, 0.27])
    assert result["modified_phases"][1] == pytest.approx([0.27, 0.73])
    assert result["phase_charge_residuals"] == pytest.approx([0.0, 0.0], abs=1.0e-15)
    assert result["ordinary_balance"] == pytest.approx(PHYSICAL_FEED, abs=1.0e-15)
    assert result["modified_balance"] == pytest.approx([0.5, 0.5], abs=1.0e-15)
    assert result["transformed_modified_potentials"] == pytest.approx([3.0, 1.0])


def test_held2_adapter_direct_objective_and_gradient_match_oracle() -> None:
    center = (0.27, 0.73, 0.92, 1.08)
    direction = (0.03, -0.04, 0.05, -0.02)
    step = 1.0e-5

    result = _equilibrium._held2_adapter(CHARGES, PHYSICAL_FEED, center, CHEMICAL_POTENTIALS)
    lower_vector = tuple(
        value - step * delta for value, delta in zip(center, direction, strict=True)
    )
    upper_vector = tuple(
        value + step * delta for value, delta in zip(center, direction, strict=True)
    )
    lower = _equilibrium._held2_adapter(CHARGES, PHYSICAL_FEED, lower_vector, CHEMICAL_POTENTIALS)
    upper = _equilibrium._held2_adapter(CHARGES, PHYSICAL_FEED, upper_vector, CHEMICAL_POTENTIALS)

    assert result["objective"] == pytest.approx(_objective_oracle(center), abs=2.0e-14)
    fraction = result["phase_fraction"]
    alpha_gradient = _manufactured_gradient(center[0], center[2])
    beta_gradient = _manufactured_gradient(center[1], center[3])
    assert result["phase_gibbs_gradients"][0] == pytest.approx(
        [alpha_gradient[0], alpha_gradient[1] + 1.0]
    )
    assert result["phase_gibbs_gradients"][1] == pytest.approx(
        [beta_gradient[0], beta_gradient[1] + 1.0]
    )
    assert 0.0 < fraction < 1.0
    numerical_directional = (upper["objective"] - lower["objective"]) / (2.0 * step)
    analytic_directional = sum(
        value * delta for value, delta in zip(result["gradient"], direction, strict=True)
    )
    assert numerical_directional == pytest.approx(analytic_directional, rel=2.0e-9, abs=2.0e-10)


def test_held2_adapter_optimum_matches_oracle_certificate() -> None:
    result = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        (0.2, 0.8, 1.0, 1.0),
        CHEMICAL_POTENTIALS,
    )

    assert result["phase_fraction"] == pytest.approx(0.5, abs=1.0e-12)
    assert result["gradient"] == pytest.approx([0.0, 0.0, 0.0, 0.0], abs=1.0e-12)
    assert result["modified_potential_gap"] < 1.0e-12
    assert result["pressure_stationarity_inf_norm"] < 1.0e-12
    assert result["certificate"]["accepted"] is True
    assert result["certificate"]["independent_evidence"] is True
    assert result["certificate"]["metrics"] == pytest.approx(
        {
            "modified_balance_abs": 0.0,
            "ordinary_balance_inf_norm": 0.0,
            "phase_charge_inf_norm": 0.0,
            "modified_potential_gap": 0.0,
            "pressure_stationarity_inf_norm": 0.0,
            "reduced_kkt_inf_norm": 0.0,
            "enumeration_objective_gap": 0.0,
            "independent_modified_composition_count": 1.0,
        },
        abs=1.0e-10,
    )


@pytest.mark.parametrize(
    ("charges", "feed", "message"),
    [
        ((0.0, 2.0, 2.0, -1.0), (0.4, 0.1, 0.1, 0.4), "singular modified coordinate"),
        (CHARGES, (0.5, 0.3, 0.2), "electroneutral"),
        ((0.0, 0.0, 0.0), (0.3, 0.3, 0.4), "charged species"),
    ],
)
def test_held2_transform_rejects_invalid_electrolyte_coordinates(
    charges: tuple[float, ...],
    feed: tuple[float, ...],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        _equilibrium._held2_adapter(
            charges,
            feed,
            (0.2, 0.8, 1.0, 1.0),
            tuple(float(index) for index in range(len(charges))),
        )


def test_held2_adapter_rejects_nonfinite_or_nonbracketing_stage_iii_state() -> None:
    with pytest.raises(ValueError, match="four finite"):
        _equilibrium._held2_adapter(
            CHARGES,
            PHYSICAL_FEED,
            (0.2, 0.8, math.nan, 1.0),
            CHEMICAL_POTENTIALS,
        )
    with pytest.raises(ValueError, match="straddle"):
        _equilibrium._held2_adapter(
            CHARGES,
            PHYSICAL_FEED,
            (0.2, 0.4, 1.0, 1.0),
            CHEMICAL_POTENTIALS,
        )


def test_held2_modified_potentials_are_galvani_gauge_invariant() -> None:
    baseline = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        (0.2, 0.8, 1.0, 1.0),
        CHEMICAL_POTENTIALS,
    )
    shifted_potentials = tuple(
        potential + charge * 17.0
        for potential, charge in zip(CHEMICAL_POTENTIALS, CHARGES, strict=True)
    )
    shifted = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        (0.2, 0.8, 1.0, 1.0),
        shifted_potentials,
    )

    assert shifted["transformed_modified_potentials"] == pytest.approx(
        baseline["transformed_modified_potentials"], abs=1.0e-12
    )


def test_held2_multivalent_coordinates_select_largest_charge_and_round_trip() -> None:
    charges = (0.0, 1.0, 2.0, -1.0, -2.0)
    feed = (0.6, 0.1, 0.1, 0.1, 0.1)
    result = _equilibrium._held2_adapter(
        charges,
        feed,
        (0.0, 1.0, 2.0, 3.0, 4.0),
    )

    assert result["eliminated_index"] == 4
    assert result["dependent_index"] == 0
    assert result["retained_indices"] == [0, 1, 2, 3]
    assert result["independent_indices"] == [1, 2, 3]
    assert result["modified_factors"] == pytest.approx([1.0, 1.5, 2.0, 0.5])
    assert result["modified_feed"] == pytest.approx([0.6, 0.15, 0.2, 0.05])
    assert result["lifted_feed"] == pytest.approx(feed, abs=1.0e-15)
    assert result["independent_lower_bounds"] == pytest.approx([1.5e-10, 2.0e-10, 5.0e-11])
    assert result["independent_upper_bounds"] == pytest.approx([1.0, 1.0, 1.0 / 3.0])
    assert result["transformed_modified_potentials"] == pytest.approx([0.0, 2.0, 3.0, 2.0])


def test_held2_manufactured_enforces_declared_modified_composition_bounds() -> None:
    with pytest.raises(ValueError, match="declared source bounds"):
        _equilibrium._held2_adapter(
            CHARGES,
            PHYSICAL_FEED,
            (1.0e-12, 0.8, 1.0, 1.0),
            CHEMICAL_POTENTIALS,
        )
    with pytest.raises(ValueError, match="declared source bounds"):
        _equilibrium._held2_adapter(
            CHARGES,
            PHYSICAL_FEED,
            (0.2, 1.0 + 1.0e-12, 1.0, 1.0),
            CHEMICAL_POTENTIALS,
        )
