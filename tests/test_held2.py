from __future__ import annotations

import math
from collections.abc import Sequence

import epcsaft
import pytest

from epcsaft_equilibrium import _equilibrium

FORMULATION_ID = "perdomo-held2.modified-mole.manufactured.v1"
CHARGES = (0.0, 1.0, -1.0)
PHYSICAL_FEED = (0.5, 0.25, 0.25)
CHEMICAL_POTENTIALS = (3.0, -2.0, 4.0)
GAS_CONSTANT_J_PER_MOL_K = 8.31446261815324
KHUDAIDA_COMPONENT_IDS = (
    "water",
    "ethanol",
    "isobutanol",
    "sodium-cation",
    "chloride-anion",
)
KHUDAIDA_AMOUNTS = (0.65905, 0.05095, 0.27015, 0.01985, 0.01985)
KHUDAIDA_CHARGES = (0.0, 0.0, 0.0, 1.0, -1.0)
KHUDAIDA_FINGERPRINT = "sha256:5a59828d86bb29c919513484a26cedaa0f025463aaa7c149ae3d1fbd0eda97ae"


def _figiel_brine_model() -> epcsaft.EPCSAFT:
    parameters = epcsaft.ParameterBundle.from_catalog(
        "figiel-2025-reference-electrolytes", version=1
    ).select(("water", "sodium-cation", "chloride-anion"))
    return epcsaft.EPCSAFT(parameters)


def _khudaida_model() -> epcsaft.EPCSAFT:
    parameters = epcsaft.ParameterBundle.from_catalog(
        "khudaida-2026-figure-2-electrolyte-lle", version=1
    ).select(KHUDAIDA_COMPONENT_IDS)
    return epcsaft.EPCSAFT(parameters)


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


def test_held2_phase_block_affine_transform_has_exact_gradient_and_hessian() -> None:
    charges = (0.0, 1.0, 2.0, -1.0, -2.0)
    center = (0.12, 0.18, 0.08, math.log(0.9))
    direction = (0.03, -0.02, 0.01, -0.04)
    pressure_over_rt = 1.3
    target_pressure_pa = 1.3
    linear = (0.2, -0.1, 0.3, -0.2, 0.4, -0.15)
    hessian = tuple(
        tuple(
            (1.0 if row == column else 0.0) + 0.02 * (min(row, column) + 1) * (max(row, column) + 2)
            for column in range(6)
        )
        for row in range(6)
    )

    def physical_coordinates(values: tuple[float, ...]) -> tuple[float, ...]:
        first, second, third, log_volume = values
        modified = (1.0 - first - second - third, first, second, third)
        physical = [0.0] * 5
        physical[0] = modified[0]
        physical[1] = modified[1] / 1.5
        physical[2] = modified[2] / 2.0
        physical[3] = modified[3] / 0.5
        physical[4] = (physical[1] + 2.0 * physical[2] - physical[3]) / 2.0
        return (*physical, math.exp(log_volume))

    def block(values: tuple[float, ...]) -> tuple[float, tuple[float, ...], tuple[float, ...]]:
        coordinates = physical_coordinates(values)
        gradient = tuple(
            linear[row] + sum(hessian[row][column] * coordinates[column] for column in range(6))
            for row in range(6)
        )
        value = 0.7 + sum(linear[index] * coordinates[index] for index in range(6))
        value += 0.5 * sum(
            coordinates[row] * hessian[row][column] * coordinates[column]
            for row in range(6)
            for column in range(6)
        )
        return value, gradient, tuple(value for row in hessian for value in row)

    def evaluate(values: tuple[float, ...]) -> dict[str, object]:
        value, gradient, flat_hessian = block(values)
        return _equilibrium._held2_adapter(
            charges,
            values[:3],
            values[3],
            pressure_over_rt,
            target_pressure_pa,
            value,
            gradient,
            flat_hessian,
            -gradient[-1],
        )

    step = 2.0e-5
    lower = evaluate(
        tuple(value - step * delta for value, delta in zip(center, direction, strict=True))
    )
    result = evaluate(center)
    upper = evaluate(
        tuple(value + step * delta for value, delta in zip(center, direction, strict=True))
    )

    expected_coordinates = physical_coordinates(center)
    assert result["physical_amounts"] == pytest.approx(expected_coordinates[:5], abs=1.0e-15)
    assert result["volume"] == pytest.approx(expected_coordinates[-1], abs=1.0e-15)
    assert result["objective"] == pytest.approx(
        block(center)[0] + pressure_over_rt * expected_coordinates[-1], abs=1.0e-14
    )
    objective_directional = (upper["objective"] - lower["objective"]) / (2.0 * step)
    assert objective_directional == pytest.approx(
        sum(value * delta for value, delta in zip(result["gradient"], direction, strict=True)),
        rel=2.0e-9,
        abs=2.0e-10,
    )
    for row in range(4):
        gradient_directional = (upper["gradient"][row] - lower["gradient"][row]) / (2.0 * step)
        hessian_directional = sum(
            result["hessian"][4 * row + column] * direction[column] for column in range(4)
        )
        assert gradient_directional == pytest.approx(hessian_directional, rel=3.0e-9, abs=3.0e-10)
    assert tuple(result["hessian"]) == pytest.approx(
        tuple(result["hessian"][4 * column + row] for row in range(4) for column in range(4)),
        abs=2.0e-14,
    )


def test_held2_installed_electrolyte_sdk_phase_block_has_exact_reduced_derivatives() -> None:
    model = _figiel_brine_model()
    capsule = epcsaft.native_sdk(model)
    temperature_k = 298.15
    pressure_pa = 100_000.0
    center = (0.02, math.log(1.0e-3))
    direction = (0.007, -0.03)

    def evaluate(values: tuple[float, float]) -> dict[str, object]:
        return _equilibrium._held2_adapter(
            capsule,
            temperature_k,
            pressure_pa,
            values[0:1],
            values[1],
            model.parameter_fingerprint,
        )

    step = 2.0e-5
    lower = evaluate(
        tuple(value - step * delta for value, delta in zip(center, direction, strict=True))
    )
    result = evaluate(center)
    upper = evaluate(
        tuple(value + step * delta for value, delta in zip(center, direction, strict=True))
    )

    assert result["component_ids"] == ["water", "sodium-cation", "chloride-anion"]
    assert result["charges"] == [0.0, 1.0, -1.0]
    assert result["physical_amounts"] == pytest.approx([0.98, 0.01, 0.01], abs=1.0e-15)
    assert result["parameter_fingerprint"] == model.parameter_fingerprint
    assert result["globality_certificate"] == "not_guaranteed"
    objective_directional = (upper["objective"] - lower["objective"]) / (2.0 * step)
    assert objective_directional == pytest.approx(
        sum(value * delta for value, delta in zip(result["gradient"], direction, strict=True)),
        rel=1.0e-9,
        abs=1.0e-10,
    )
    for row in range(2):
        gradient_directional = (upper["gradient"][row] - lower["gradient"][row]) / (2.0 * step)
        hessian_directional = sum(
            result["hessian"][2 * row + column] * direction[column] for column in range(2)
        )
        assert gradient_directional == pytest.approx(hessian_directional, rel=1.0e-8, abs=1.0e-10)
    packing_directional = (upper["packing_fraction"] - lower["packing_fraction"]) / (2.0 * step)
    assert packing_directional == pytest.approx(
        sum(
            value * delta
            for value, delta in zip(result["packing_gradient"], direction, strict=True)
        ),
        rel=1.0e-9,
        abs=1.0e-11,
    )
    for row in range(2):
        packing_gradient_directional = (
            upper["packing_gradient"][row] - lower["packing_gradient"][row]
        ) / (2.0 * step)
        packing_hessian_directional = sum(
            result["packing_hessian"][2 * row + column] * direction[column] for column in range(2)
        )
        assert packing_gradient_directional == pytest.approx(
            packing_hessian_directional, rel=1.0e-8, abs=1.0e-11
        )
    assert tuple(result["packing_hessian"]) == pytest.approx(
        tuple(
            result["packing_hessian"][2 * column + row] for row in range(2) for column in range(2)
        ),
        abs=2.0e-14,
    )
    assert 1.0e-6 <= result["packing_fraction"] <= 0.74
    assert result["pressure_over_rt"] == pytest.approx(
        pressure_pa / (GAS_CONSTANT_J_PER_MOL_K * temperature_k), rel=2.0e-15
    )
    assert math.isfinite(result["provider_pressure_pa"])


def test_held2_installed_log_packing_coordinate_round_trips_and_has_exact_derivatives() -> None:
    model = _khudaida_model()
    capsule = epcsaft.native_sdk(model)
    temperature_k = 293.15
    pressure_pa = 100_000.0
    states = (
        (
            (0.6462224836985831, 0.04995832720498113, 0.03892729322939648),
            math.log(3.909560419950043e-05),
        ),
        (
            (0.48317652139878337, 0.06575359819988703, 0.012557176718000043),
            math.log(5.2811607028165594e-05),
        ),
        (
            (0.49974816254115306, 0.34140273164275964, 0.10555255266393769),
            -2.0847367866112805,
        ),
    )
    direction = (0.011, -0.008, 0.006, 0.017)
    step = 2.0e-4

    for independent, log_volume in states:
        physical = _equilibrium._held2_adapter(
            capsule,
            temperature_k,
            pressure_pa,
            independent,
            log_volume,
            KHUDAIDA_FINGERPRINT,
        )
        center = (*independent, math.log(physical["packing_fraction"]))

        def evaluate(values: tuple[float, ...]) -> dict[str, object]:
            return _equilibrium._held2_adapter(
                capsule,
                temperature_k,
                pressure_pa,
                values[:-1],
                values[-1],
                KHUDAIDA_FINGERPRINT,
                "log_packing_phase",
            )

        lower = evaluate(
            tuple(value - step * delta for value, delta in zip(center, direction, strict=True))
        )
        result = evaluate(center)
        upper = evaluate(
            tuple(value + step * delta for value, delta in zip(center, direction, strict=True))
        )

        assert result["log_volume"] == pytest.approx(log_volume, abs=2.0e-11)
        assert result["log_packing_residual"] == pytest.approx(0.0, abs=2.0e-13)
        assert math.log(result["packing_fraction"]) == pytest.approx(center[-1], abs=2.0e-13)
        objective_directional = (upper["objective"] - lower["objective"]) / (2.0 * step)
        assert objective_directional == pytest.approx(
            sum(value * delta for value, delta in zip(result["gradient"], direction, strict=True)),
            rel=2.0e-8,
            abs=2.0e-9,
        )
        for row in range(4):
            gradient_directional = (upper["gradient"][row] - lower["gradient"][row]) / (2.0 * step)
            hessian_directional = sum(
                result["hessian"][4 * row + column] * direction[column] for column in range(4)
            )
            assert gradient_directional == pytest.approx(
                hessian_directional, rel=3.0e-8, abs=3.0e-9
            )
        assert tuple(result["hessian"]) == pytest.approx(
            tuple(result["hessian"][4 * column + row] for row in range(4) for column in range(4)),
            abs=3.0e-13,
        )


def test_held2_installed_log_packing_coordinate_completes_former_failed_start_zero() -> None:
    model = _khudaida_model()
    physical_feed = tuple(amount / sum(KHUDAIDA_AMOUNTS) for amount in KHUDAIDA_AMOUNTS)

    result = _equilibrium._held2_adapter(
        epcsaft.native_sdk(model),
        293.15,
        100_000.0,
        physical_feed,
        KHUDAIDA_FINGERPRINT,
        "stage_i_start_0",
    )

    assert result["attempted_start_count"] == 1
    assert result["completed_start_count"] == 1
    assert result["failed_start_count"] == 0
    assert result["failed_start_index"] is None


def test_held2_installed_stage_i_uses_provider_volume_domain_for_khudaida_midpoint() -> None:
    model = _khudaida_model()
    physical_feed = tuple(amount / sum(KHUDAIDA_AMOUNTS) for amount in KHUDAIDA_AMOUNTS)

    result = _equilibrium._held2_adapter(
        epcsaft.native_sdk(model),
        293.15,
        100_000.0,
        physical_feed,
        KHUDAIDA_FINGERPRINT,
        "stage_i",
    )

    assert model.parameter_fingerprint == KHUDAIDA_FINGERPRINT
    assert result["component_ids"] == list(KHUDAIDA_COMPONENT_IDS)
    assert result["charges"] == [0.0, 0.0, 0.0, 1.0, -1.0]
    assert result["parameter_fingerprint"] == KHUDAIDA_FINGERPRINT
    assert result["packing_fraction_bounds"] == pytest.approx([1.0e-6, 0.74])
    assert result["molar_volume_bounds"] == pytest.approx(
        [2.306595376485171e-05, 17.068805785990268], rel=2.0e-15
    )
    assert result["log_molar_volume_bounds"] == pytest.approx(
        [math.log(2.306595376485171e-05), math.log(17.068805785990268)], rel=2.0e-15
    )
    assert result["declared_start_count"] == 10 * len(KHUDAIDA_COMPONENT_IDS)
    assert result["reference_scan_interval_count"] == 50
    assert result["reference_scan_point_count"] == 51
    assert result["reference_evaluation_failure_count"] == 0
    assert result["reference_refinement_failure_count"] == 0
    assert result["reference_root_count"] == 3
    assert result["reference_stable_root_count"] == 2
    assert [root["mechanically_stable"] for root in result["reference_roots"]] == [
        True,
        False,
        True,
    ]
    assert [root["volume"] for root in result["reference_roots"]] == pytest.approx(
        [3.909560419950043e-05, 0.0007107348724168745, 0.02238310519168927],
        rel=2.0e-8,
    )
    assert [root["curvature"] for root in result["reference_roots"]] == pytest.approx(
        [31.411491224899734, -0.5642708323753126, 0.8463122865798167],
        rel=2.0e-8,
    )
    assert result["completed_start_count"] == 11
    assert result["failed_start_count"] == 39
    assert result["completed_start_count"] + result["failed_start_count"] == 50
    assert result["search_completeness"] == "incomplete"
    assert result["failed_start_index"] == 0
    assert result["failed_start_solver_status"] == 3
    assert result["failed_start_solver_converged"] is False
    assert result["failed_start_reason"] == "TPD solve did not return a complete accepted state"
    assert result["failed_start_initial"] == pytest.approx(
        [
            0.49974816254115306,
            0.34140273164275964,
            0.10555255266393769,
            -2.0847367866112805,
        ]
    )
    assert result["outcome"] == "negative_tpd"
    assert result["volume_domain_search_complete"] is True
    assert result["candidate_domain_evaluation_failure_count"] == 0
    assert result["candidate_domain_rejection_count"] == 0
    assert result["reference_volume"] == pytest.approx(3.909560419950043e-05, rel=2.0e-8)
    assert result["minimum_tpd"] == pytest.approx(-0.14499309134029337, abs=2.0e-13)
    assert [candidate["tpd"] for candidate in result["candidates"]] == pytest.approx(
        [-0.00893217913694587], abs=2.0e-13
    )
    assert [candidate["volume"] for candidate in result["candidates"]] == pytest.approx(
        [5.2811607028165594e-05], rel=2.0e-12
    )
    phase_evidence = []
    for candidate in result["candidates"]:
        assert 1.0e-6 <= candidate["packing_fraction"] <= 0.74
        modified = candidate["modified_fractions"]
        phase = _equilibrium._held2_adapter(
            epcsaft.native_sdk(model),
            293.15,
            100_000.0,
            (modified[0], modified[1], modified[3]),
            math.log(candidate["volume"]),
            KHUDAIDA_FINGERPRINT,
        )
        assert sum(phase["physical_amounts"]) == pytest.approx(1.0, abs=2.0e-15)
        assert phase["physical_amounts"][3] == pytest.approx(
            phase["physical_amounts"][4], abs=2.0e-15
        )
        assert (
            max(
                abs(value - feed)
                for value, feed in zip(
                    modified, result["reference_modified_fractions"], strict=True
                )
            )
            > 1.0e-3
        )
        assert candidate["molar_volume_bounds"] == pytest.approx(
            [3.079084484139457e-05, 22.785225182631986], rel=2.0e-14
        )
        assert candidate["molar_volume_bounds"][0] <= candidate["volume"]
        assert candidate["volume"] <= candidate["molar_volume_bounds"][1]
        phase_evidence.append(phase)
    assert abs(phase_evidence[0]["pressure_stationarity_relative"]) <= 1.0e-8
    assert result["candidates"][0]["lower_volume_bound_active"] is False
    assert result["candidates"][0]["upper_volume_bound_active"] is False
    assert result["globality_certificate"] == "not_guaranteed"


def test_held2_manufactured_stage_i_finds_negative_tpd_with_declared_multistart() -> None:
    result = _equilibrium._held2_adapter(CHARGES, PHYSICAL_FEED, "stage_i")

    assert result["profile"] == "perdomo-held2-stage-i-manufactured-v1"
    assert result["outcome"] == "negative_tpd"
    assert result["globality_certificate"] == "not_guaranteed"
    assert result["declared_start_count"] == 10 * len(CHARGES)
    assert result["completed_start_count"] == result["declared_start_count"]
    assert result["failed_start_count"] == 0
    assert result["volume_domain_search_complete"] is True
    assert result["reference_modified_fractions"] == pytest.approx([0.5, 0.5])
    assert result["reference_volume"] == pytest.approx(1.0, abs=2.0e-9)
    assert result["minimum_tpd"] < -1.0e-8
    assert sorted(
        candidate["modified_fractions"][1] for candidate in result["candidates"]
    ) == pytest.approx([0.2, 0.8], abs=2.0e-7)
    assert all(candidate["tpd"] < -1.0e-8 for candidate in result["candidates"])


def test_held2_stage_i_requires_complete_search_before_no_negative_found() -> None:
    result = _equilibrium._held2_adapter(CHARGES, (0.9, 0.05, 0.05), "stage_i")

    assert result["outcome"] == "no_negative_found"
    assert result["completed_start_count"] == result["declared_start_count"] == 30
    assert result["failed_start_count"] == 0
    assert result["volume_domain_search_complete"] is True
    assert result["failed_start_index"] is None
    assert math.isfinite(result["minimum_tpd"])
    assert result["candidates"] == []
    assert result["globality_certificate"] == "not_guaranteed"


def test_held2_manufactured_stage_ii_builds_replayable_candidate_set() -> None:
    result = _equilibrium._held2_adapter(CHARGES, PHYSICAL_FEED, "stage_ii")

    assert result["profile"] == "perdomo-held2-stage-ii-manufactured-v1"
    assert result["outcome"] == "candidate_set"
    assert result["globality_certificate"] == "not_guaranteed"
    assert result["major_iterations"] <= 100
    assert result["lower_starts_per_iteration"] == 30
    assert len(result["bound_history"]) == result["major_iterations"]
    assert all(
        entry["lower_bound"] <= entry["upper_bound"] + 1.0e-10 for entry in result["bound_history"]
    )
    assert result["cut_count"] >= 3
    assert sorted(
        candidate["modified_fractions"][1] for candidate in result["candidates"]
    ) == pytest.approx([0.2, 0.8], abs=2.0e-7)
    assert all(abs(candidate["lower_gap"]) <= 1.0e-8 for candidate in result["candidates"])


def test_held2_stage_ii_certified_improving_cut_advances_with_partial_search() -> None:
    result = _equilibrium._held2_adapter(
        -6.0097420675620405,
        -6.5902644733135745,
        34,
        50,
        "stage_ii_precedence",
    )

    assert result == {
        "decision": "add_improving_cut",
        "search_completeness": "partial",
        "lower_bound_certified": False,
    }


def test_held2_stage_ii_partial_search_without_improving_cut_is_indeterminate() -> None:
    result = _equilibrium._held2_adapter(-6.0, -5.99, 34, 50, "stage_ii_precedence")

    assert result == {
        "decision": "indeterminate",
        "search_completeness": "partial",
        "lower_bound_certified": False,
    }


def test_held2_stage_ii_final_gap_requires_complete_lower_search() -> None:
    result = _equilibrium._held2_adapter(-6.0, -6.0, 50, 50, "stage_ii_precedence")

    assert result == {
        "decision": "evaluate_final_gap",
        "search_completeness": "complete",
        "lower_bound_certified": True,
    }


def test_held2_general_mp_stage_iii_exact_lagrangian_hessian() -> None:
    candidates = ((0.2, 1.0), (0.20000002, 1.0), (0.8, 1.0))
    center = (0.25, 0.2, 1.0, 0.25, 0.21, 1.0, 0.5, 0.795, 1.0)
    multipliers = (0.3, -0.2)
    direction = (0.01, -0.02, 0.03, -0.01, 0.015, -0.02, 0.0, 0.005, 0.01)

    def evaluate(values: tuple[float, ...]) -> dict[str, object]:
        return _equilibrium._held2_adapter(
            CHARGES,
            PHYSICAL_FEED,
            candidates,
            values,
            multipliers,
            "stage_iii_derivatives",
        )

    step = 2.0e-5
    lower = evaluate(
        tuple(value - step * delta for value, delta in zip(center, direction, strict=True))
    )
    result = evaluate(center)
    upper = evaluate(
        tuple(value + step * delta for value, delta in zip(center, direction, strict=True))
    )

    objective_directional = (upper["objective"] - lower["objective"]) / (2.0 * step)
    assert objective_directional == pytest.approx(
        sum(
            value * delta
            for value, delta in zip(result["objective_gradient"], direction, strict=True)
        ),
        rel=2.0e-9,
        abs=2.0e-10,
    )
    size = len(center)
    for row in range(size):
        numerical = (upper["lagrangian_gradient"][row] - lower["lagrangian_gradient"][row]) / (
            2.0 * step
        )
        analytic = sum(
            result["lagrangian_hessian"][size * row + column] * direction[column]
            for column in range(size)
        )
        assert numerical == pytest.approx(analytic, rel=3.0e-8, abs=3.0e-9)
    assert tuple(result["lagrangian_hessian"]) == pytest.approx(
        tuple(
            result["lagrangian_hessian"][size * column + row]
            for row in range(size)
            for column in range(size)
        ),
        abs=2.0e-14,
    )


def test_held2_general_mp_stage_iii_simplex_schema_and_multiplier_mapping() -> None:
    feed = tuple(value / sum(KHUDAIDA_AMOUNTS) for value in KHUDAIDA_AMOUNTS)
    candidates = (
        (0.6588026822616131, 1.108849693428444e-5, 0.34115777594234437, -0.8047587326883016),
        (0.9104972678771838, 3.046088482260374e-5, 0.08850404209802903, -0.7857111519928206),
        (0.9784651953324194, 0.01956780416816369, 2.420187102135475e-5, -0.7872645327604558),
    )
    fraction = 1.0 / len(candidates)
    variables = tuple(value for candidate in candidates for value in (fraction, *candidate))
    equality_multipliers = (0.3, -0.2, 0.4, -0.1)
    zero_simplex = (*equality_multipliers, 0.0, 0.0, 0.0)
    simplex_multipliers = (1.1, -0.7, 0.4)

    baseline = _equilibrium._held2_adapter(
        KHUDAIDA_CHARGES,
        feed,
        candidates,
        variables,
        zero_simplex,
        "stage_iii_general_schema",
    )
    result = _equilibrium._held2_adapter(
        KHUDAIDA_CHARGES,
        feed,
        candidates,
        variables,
        (*equality_multipliers, *simplex_multipliers),
        "stage_iii_general_schema",
    )

    assert result["variable_count"] == 15
    assert result["constraint_count"] == 7
    assert result["jacobian_nonzero_count"] == 30
    assert result["constraint_lower_bounds"] == [0.0] * 4 + [-2.0e19] * 3
    assert result["constraint_upper_bounds"] == pytest.approx(
        [0.0] * 4 + [1.0 - 1.0e-10] * 3,
        abs=0.0,
    )
    assert result["constraints"][-3:] == pytest.approx(
        [sum(candidate[:-1]) for candidate in candidates],
        abs=2.0e-15,
    )

    expected_rows = [0, 0, 0]
    expected_columns = [0, 5, 10]
    expected_values = [1.0, 1.0, 1.0]
    for coordinate in range(3):
        for phase, candidate in enumerate(candidates):
            offset = 5 * phase
            expected_rows.extend((coordinate + 1, coordinate + 1))
            expected_columns.extend((offset, offset + coordinate + 1))
            expected_values.extend((candidate[coordinate], fraction))
    for phase in range(3):
        for coordinate in range(3):
            expected_rows.append(4 + phase)
            expected_columns.append(5 * phase + coordinate + 1)
            expected_values.append(1.0)
    assert result["jacobian_rows"] == expected_rows
    assert result["jacobian_columns"] == expected_columns
    assert result["jacobian_values"] == pytest.approx(expected_values, abs=0.0)

    gradient_difference = tuple(
        value - base
        for value, base in zip(
            result["lagrangian_gradient"],
            baseline["lagrangian_gradient"],
            strict=True,
        )
    )
    expected_difference = [0.0] * len(variables)
    for phase, multiplier in enumerate(simplex_multipliers):
        for coordinate in range(3):
            expected_difference[5 * phase + coordinate + 1] = multiplier
    assert gradient_difference == pytest.approx(expected_difference, abs=2.0e-15)
    assert result["lagrangian_hessian"] == pytest.approx(baseline["lagrangian_hessian"], abs=0.0)


def test_held2_general_mp_stage_iii_rejects_invalid_trial_before_provider() -> None:
    feed = tuple(value / sum(KHUDAIDA_AMOUNTS) for value in KHUDAIDA_AMOUNTS)
    candidates = (
        (0.6588026822616131, 1.108849693428444e-5, 0.34115777594234437, -0.8047587326883016),
        (0.9104972678771838, 3.046088482260374e-5, 0.08850404209802903, -0.7857111519928206),
    )
    invalid_trial = (
        0.5,
        0.6583790364854278,
        0.0010010884969342845,
        0.3407281048600504,
        -0.8053132473085454,
        0.5,
        *candidates[1],
    )

    result = _equilibrium._held2_adapter(
        KHUDAIDA_CHARGES,
        feed,
        candidates,
        invalid_trial,
        (0.0,) * 6,
        "stage_iii_invalid_trial",
    )

    assert result["objective_accepted"] is False
    assert result["scientific_evaluator_call_count"] == 0
    assert result["recoverable_domain_rejection_count"] == 1
    assert result["last_domain_rejection"]
    assert result["fatal_callback_error"] == ""


def test_held2_general_mp_stage_iii_refines_then_merges_duplicate_phases() -> None:
    result = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        ((0.2, 1.0), (0.20000002, 1.0), (0.8, 1.0)),
        "stage_iii",
    )

    assert result["profile"] == "perdomo-held2-stage-iii-manufactured-v1"
    assert result["solver_status"] in {"solve_succeeded", "solved_to_acceptable_level"}
    assert result["numerical_status"] == "converged"
    assert result["physical_status"] == "accepted"
    assert result["globality_certificate"] == "not_guaranteed"
    assert result["feedback"] == "none"
    assert result["input_candidate_count"] == 3
    assert result["retired_duplicate_count"] == 1
    assert len(result["phases"]) == 2
    assert sorted(phase["modified_fractions"][1] for phase in result["phases"]) == pytest.approx(
        [0.2, 0.8], abs=2.0e-7
    )
    assert sorted(phase["phase_fraction"] for phase in result["phases"]) == pytest.approx(
        [0.5, 0.5], abs=2.0e-7
    )
    assert result["modified_balance_inf_norm"] < 1.0e-9
    assert result["ordinary_balance_inf_norm"] < 1.0e-9
    assert result["phase_charge_inf_norm"] < 1.0e-12
    assert result["pressure_stationarity_inf_norm"] < 1.0e-9
    assert result["modified_potential_mixed_gap"] < 1.0e-9
    assert result["certified_modified_potential_count"] == 2
    assert result["trace_component_count"] == 0
    assert result["trace_refinement_status"] == "not_required"
    assert result["minimum_phase_distance"] > 1.0e-3
    assert result["kkt_stationarity_inf_norm"] < 1.0e-7
    assert abs(result["enumeration_objective_gap"]) < 1.0e-9


def test_held2_general_mp_stage_iii_returns_infeasible_set_to_stage_ii() -> None:
    result = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        ((0.1, 1.0), (0.2, 1.0), (0.3, 1.0)),
        "stage_iii",
    )

    assert result["numerical_status"] == "not_adjudicated"
    assert result["physical_status"] == "not_adjudicated"
    assert result["feedback"] == "return_to_stage_ii"
    assert result["failure_reason"] == "candidate_set_does_not_bracket_feed"
    assert result["globality_certificate"] == "not_guaranteed"
