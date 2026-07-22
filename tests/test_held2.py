from __future__ import annotations

import math
from collections.abc import Sequence
from itertools import pairwise

import epcsaft
import pytest

from epcsaft_equilibrium import _equilibrium

FORMULATION_ID = "perdomo-held2.modified-mole.manufactured.v1"
CHARGES = (0.0, 1.0, -1.0)
PHYSICAL_FEED = (0.5, 0.25, 0.25)
CHEMICAL_POTENTIALS = (3.0, -2.0, 4.0)
GAS_CONSTANT_J_PER_MOL_K = 8.31446261815324


def _figiel_brine_model() -> epcsaft.EPCSAFT:
    parameters = epcsaft.ParameterBundle.from_catalog(
        "figiel-2025-reference-electrolytes", version=1
    ).select(("water", "sodium-cation", "chloride-anion"))
    return epcsaft.EPCSAFT(parameters)


def _khudaida_electrolyte_lle_model() -> epcsaft.EPCSAFT:
    parameters = epcsaft.ParameterBundle.from_catalog(
        "khudaida-2026-figure-2-electrolyte-lle", version=1
    ).select(
        (
            "water",
            "ethanol",
            "isobutanol",
            "sodium-cation",
            "chloride-anion",
        )
    )
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
    assert result["pressure_over_rt"] == pytest.approx(
        pressure_pa / (GAS_CONSTANT_J_PER_MOL_K * temperature_k), rel=2.0e-15
    )
    assert math.isfinite(result["provider_pressure_pa"])


def test_held2_manufactured_stage_i_finds_negative_tpd_with_declared_multistart() -> None:
    result = _equilibrium._held2_adapter(CHARGES, PHYSICAL_FEED, "stage_i")

    assert result["profile"] == "perdomo-held2-stage-i-manufactured-v1"
    assert result["outcome"] == "negative_tpd"
    assert result["globality_certificate"] == "not_guaranteed"
    assert result["declared_start_count"] == 10 * len(CHARGES)
    assert result["completed_start_count"] == result["declared_start_count"]
    assert result["failed_start_count"] == 0
    assert result["reference_modified_fractions"] == pytest.approx([0.5, 0.5])
    assert result["reference_volume"] == pytest.approx(1.0, abs=2.0e-10)
    assert result["minimum_tpd"] < -1.0e-8
    assert sorted(
        candidate["modified_fractions"][1] for candidate in result["candidates"]
    ) == pytest.approx([0.2, 0.8], abs=2.0e-7)
    assert all(candidate["tpd"] < -1.0e-8 for candidate in result["candidates"])


def test_held2_stage_i_direct_l_finds_certified_negative_witness_deterministically() -> None:
    first = _equilibrium._held2_stage_i_direct("negative", 80)
    second = _equilibrium._held2_stage_i_direct("negative", 80)

    assert first == second
    assert first["outcome"] == "negative_witness_found"
    assert first["termination_reason"] == "certified_negative_tpd"
    assert first["search_strategy"] == "nlopt_direct_l_pressure_envelope_v1"
    assert first["search_solver"] == "nlopt_gn_direct_l"
    assert first["solver_version"] == "2.11.0"
    assert first["declared_evaluation_budget"] == 80
    assert 0 < first["completed_evaluation_count"] < 80
    assert first["failed_evaluation_count"] == 0
    assert first["globality_certificate"] == "not_guaranteed"
    witness = first["negative_witness"]
    assert witness["tpd"] < -1.0e-8
    assert witness["pressure_certified"] is True
    assert witness["mechanical_class"] == "strict_stable"
    assert witness["root_origin"] == "sign_change"
    assert witness["root_completeness"] == "not_proven"


def test_held2_stage_i_direct_l_no_negative_is_only_finite_search_evidence() -> None:
    result = _equilibrium._held2_stage_i_direct("no_negative", 40)

    assert result["outcome"] == "no_negative_witness_detected"
    assert result["termination_reason"] == "declared_budget_exhausted"
    assert result["completed_evaluation_count"] == 40
    assert result["failed_evaluation_count"] == 0
    assert result["negative_witness"] is None
    assert result["globality_certificate"] == "not_guaranteed"
    assert all(evaluation["tpd"] >= -1.0e-8 for evaluation in result["evaluations"])


def test_held2_stage_i_direct_l_finds_narrow_negative_pocket() -> None:
    result = _equilibrium._held2_stage_i_direct("narrow_negative", 160)

    assert result["outcome"] == "negative_witness_found"
    assert result["negative_witness"]["tpd"] < -1.0e-6
    assert result["negative_witness"]["independent_modified_fractions"][0] == pytest.approx(
        0.73, abs=0.06
    )


@pytest.mark.parametrize(
    ("topology", "failure_reason"),
    [
        ("boundary", "boundary_root"),
        ("provider_failure", "evaluation_failed"),
    ],
)
def test_held2_stage_i_direct_l_fails_closed_on_envelope_failure(
    topology: str,
    failure_reason: str,
) -> None:
    result = _equilibrium._held2_stage_i_direct(topology, 40)

    assert result["outcome"] == "indeterminate"
    assert result["termination_reason"] == "required_envelope_evaluation_failed"
    assert result["failed_evaluation_count"] == 1
    assert result["negative_witness"] is None
    assert result["evaluations"][-1]["failure_reason"] == failure_reason


def test_held2_stage_i_direct_l_retains_pressure_branch_switches() -> None:
    result = _equilibrium._held2_stage_i_direct("branch_switch", 60)
    selected_log_volumes = {
        round(evaluation["selected_root_log_volume"], 6)
        for evaluation in result["evaluations"]
        if evaluation["certified"]
    }

    assert any(log_volume < 0.0 for log_volume in selected_log_volumes)
    assert any(log_volume > 0.0 for log_volume in selected_log_volumes)
    assert all(
        len(evaluation["pressure_envelope"]["roots"]) == 3
        for evaluation in result["evaluations"]
        if evaluation["certified"]
    )


def test_held2_stage_i_direct_l_finds_khudaida_negative_witness() -> None:
    model = _khudaida_electrolyte_lle_model()
    result = _equilibrium._held2_stage_i_direct(
        epcsaft.native_sdk(model),
        293.15,
        100_000.0,
        (
            0.6462224836985831,
            0.04995832720498113,
            0.2648918958670393,
            0.01946364661469824,
            0.01946364661469824,
        ),
        model.parameter_fingerprint,
        50,
    )

    assert result["outcome"] == "negative_witness_found"
    assert result["declared_evaluation_budget"] == 50
    assert result["completed_evaluation_count"] <= 50
    assert result["failed_evaluation_count"] == 0
    assert result["negative_witness"]["tpd"] < -1.0e-8
    assert result["negative_witness"]["pressure_certified"] is True
    assert result["negative_witness"]["mechanical_class"] == "strict_stable"
    assert result["reference_pressure_envelope"]["outcome"] == "selected"
    assert result["reference_pressure_envelope"]["roots"][
        result["reference_pressure_envelope"]["selected_root_index"]
    ]["volume"] == pytest.approx(3.909560419949432e-5, rel=1.0e-8)
    assert result["globality_certificate"] == "not_guaranteed"


def test_held2_integrated_controller_is_fail_closed_at_declared_stage_ii_budget() -> None:
    model = _khudaida_electrolyte_lle_model()
    result = _equilibrium._held2_controller(
        epcsaft.native_sdk(model),
        293.15,
        100_000.0,
        (
            0.6462224836985831,
            0.04995832720498113,
            0.2648918958670393,
            0.01946364661469824,
            0.01946364661469824,
        ),
        model.parameter_fingerprint,
        50,
        1,
        1,
    )

    assert result["controller"] == "perdomo_held2_integrated_private_v1"
    assert result["stage_order"] == (
        "homogeneous_reference",
        "stage_i",
        "stage_ii",
        "stage_iii",
    )
    assert result["stage_i"]["outcome"] == "negative_witness_found"
    assert result["outcome"] == "indeterminate_stage_ii"
    assert result["failure_stage"] == "stage_ii"
    assert result["stage_ii"]["outcome"] in {
        "resource_limit",
        "indeterminate_lower_search",
    }
    assert result["stage_ii"]["local_attempts_truncated"] is True
    assert len(result["stage_ii"]["attempt_trace"]) == 1
    attempt = result["stage_ii"]["attempt_trace"][0]
    assert attempt["solver_converged"] is True
    assert attempt["solver_status"] == "success"
    assert attempt["callback_error"] == ""
    normalized = [
        evaluation
        for evaluation in attempt["evaluation_trace"]
        if evaluation["mapping_status"] == "normalized_boundary_contact"
    ]
    assert normalized
    assert normalized[0]["callback"] == "eval_f"
    assert normalized[0]["accepted_iterate"] is False
    assert normalized[0]["raw_variables"][2] == math.nextafter(1.0, math.inf)
    assert normalized[0]["maximum_bound_violation"] == math.ulp(1.0)
    assert normalized[0]["physical_variables"] == pytest.approx(
        [
            4.553330337150642e-06,
            1.0000123152940647e-10,
            0.9999954464696617,
            -9.391852068856018,
        ]
    )
    assert max(item["maximum_bound_violation"] for item in normalized) == 4.0 * math.ulp(1.0)
    assert any(item["accepted_iterate"] for item in normalized)
    assert attempt["last_valid_physical_variables"] == pytest.approx(
        attempt["internal_terminal"]
    )
    assert len(result["stage_ii"]["bound_history"]) == 1
    upper = result["stage_ii"]["bound_history"][0]
    assert upper["lower_bound_available"] is False
    assert upper["upper_solver"] == "highs_lp"
    assert upper["upper_solver_status"] == "Optimal"
    assert upper["upper_primal_feasible"] is True
    assert upper["upper_dual_feasible"] is True
    assert upper["upper_primal_residual_inf"] == 0.0
    assert upper["upper_dual_residual_inf"] == 0.0
    assert upper["cut_slacks"] == pytest.approx([0.0, 0.0, 0.10764612291119935, 0.0, 0.0])
    assert upper["cut_duals"] == pytest.approx([1.0, 0.0, 0.0, 0.0, 0.0])
    assert upper["active_cut_ids"] == [0, 1, 3, 4]
    assert result["stage_iii"] is None
    assert result["predictive_comparison_status"] == ("not_allowed_before_physical_acceptance")
    assert result["globality_certificate"] == "not_guaranteed"


def test_held2_solver_neutral_search_objective_matches_raw_and_tpd_modes() -> None:
    variables = (0.27, math.log(0.92))
    reference_variables = (0.5, math.log(1.2))

    raw = _equilibrium._held2_adapter(
        CHARGES,
        variables,
        reference_variables,
        False,
        "search_objective",
    )
    reference = _equilibrium._held2_adapter(
        CHARGES,
        reference_variables,
        reference_variables,
        False,
        "search_objective",
    )
    tpd = _equilibrium._held2_adapter(
        CHARGES,
        variables,
        reference_variables,
        True,
        "search_objective",
    )

    expected_tpd = (
        raw["objective"]
        - reference["objective"]
        - sum(
            reference_gradient * (value - reference_value)
            for reference_gradient, value, reference_value in zip(
                reference["gradient"], variables, reference_variables, strict=True
            )
        )
    )
    assert raw["physical_amounts"] == pytest.approx([0.73, 0.135, 0.135])
    assert raw["volume"] == pytest.approx(0.92)
    assert raw["pressure_stationarity_relative"] == pytest.approx(-0.4)
    assert tpd["objective"] == pytest.approx(expected_tpd, abs=2.0e-14)
    assert tpd["gradient"] == pytest.approx(
        [
            value - reference_value
            for value, reference_value in zip(raw["gradient"], reference["gradient"], strict=True)
        ]
    )
    assert tpd["hessian"] == pytest.approx(raw["hessian"])
    assert tpd["modified_fractions"] == pytest.approx(raw["modified_fractions"])
    assert tpd["modified_potentials"] == pytest.approx(raw["modified_potentials"])


def test_held2_solver_neutral_search_objective_fails_closed() -> None:
    with pytest.raises(ValueError, match="coordinate count"):
        _equilibrium._held2_adapter(
            CHARGES,
            (0.27, math.log(0.92), 0.0),
            (0.5, math.log(1.2)),
            True,
            "search_objective",
        )
    with pytest.raises(ValueError, match="phase state has the wrong size"):
        _equilibrium._held2_adapter(
            CHARGES,
            (0.27, math.inf),
            (0.5, math.log(1.2)),
            False,
            "search_objective",
        )


@pytest.mark.parametrize(
    ("raw", "expected", "normalized"),
    [
        (0.25, 0.25, False),
        (math.nextafter(0.0, -math.inf), 0.0, True),
        (math.nextafter(1.0, math.inf), 1.0, True),
        (1.0000000000000009, 1.0, True),
    ],
)
def test_held2_stage_ii_chart_policy_normalizes_only_adjacent_binary64_contact(
    raw: float,
    expected: float,
    normalized: bool,
) -> None:
    evidence = _equilibrium._held2_adapter(raw, "stage_ii_chart_coordinate")

    assert evidence["raw_coordinate"] == raw
    assert evidence["normalized_coordinate"] == expected
    assert evidence["normalized_boundary_contact"] is normalized
    assert evidence["policy"] == "binary64_four_ulp_unit_boundary_v1"


@pytest.mark.parametrize(
    "raw",
    [
        -2.5e-323,
        1.000000000000001,
        -1.0e-12,
        1.0 + 1.0e-12,
        math.inf,
    ],
)
def test_held2_stage_ii_chart_policy_rejects_meaningful_invalid_motion(
    raw: float,
) -> None:
    with pytest.raises(ValueError, match=r"outside|not finite"):
        _equilibrium._held2_adapter(raw, "stage_ii_chart_coordinate")


def test_held2_manufactured_stage_ii_builds_replayable_candidate_set() -> None:
    result = _equilibrium._held2_adapter(CHARGES, PHYSICAL_FEED, "stage_ii")

    assert result["profile"] == "perdomo-held2-stage-ii-manufactured-v1"
    assert result["outcome"] == "candidate_set"
    assert result["globality_certificate"] == "not_guaranteed"
    assert result["search_strategy"] == "continuation_sobol_direct_l_ipopt_v1"
    assert result["global_explorer"] == "continuation_sobol_direct_l"
    assert result["local_solver"] == "ipopt_exact_hessian"
    assert result["exploration_failure_count"] == 0
    assert result["exploration_evaluation_count"] > 0
    assert result["exploration_representative_count"] == len(result["attempt_trace"])
    assert result["distinct_basin_count"] >= 3
    assert result["major_iterations"] <= 100
    assert result["lower_starts_per_iteration"] > 0
    assert len(result["bound_history"]) == result["major_iterations"]
    assert all(entry["lower_bound_available"] for entry in result["bound_history"])
    assert all(
        entry["lower_bound"] <= entry["upper_bound"] + 1.0e-10 for entry in result["bound_history"]
    )
    assert all(entry["upper_solver"] == "highs_lp" for entry in result["bound_history"])
    assert all(entry["upper_solver_version"] == "1.15.1" for entry in result["bound_history"])
    assert all(entry["upper_primal_feasible"] for entry in result["bound_history"])
    assert all(entry["upper_dual_feasible"] for entry in result["bound_history"])
    assert all(entry["upper_primal_residual_inf"] <= 1.0e-8 for entry in result["bound_history"])
    assert all(entry["upper_dual_residual_inf"] <= 1.0e-8 for entry in result["bound_history"])
    assert result["cut_count"] >= 3
    assert sorted(
        candidate["modified_fractions"][1] for candidate in result["candidates"]
    ) == pytest.approx([0.2, 0.8], abs=2.0e-7)
    assert all(abs(candidate["lower_gap"]) <= 1.0e-8 for candidate in result["candidates"])


def test_held2_manufactured_stage_ii_retains_every_lower_attempt() -> None:
    first = _equilibrium._held2_adapter(CHARGES, PHYSICAL_FEED, "stage_ii")
    second = _equilibrium._held2_adapter(CHARGES, PHYSICAL_FEED, "stage_ii")

    assert first["attempt_trace"] == second["attempt_trace"]
    assert len(first["attempt_trace"]) >= first["lower_starts_per_iteration"]

    for expected_attempt_id, attempt in enumerate(first["attempt_trace"]):
        assert attempt["attempt_id"] == expected_attempt_id
        assert attempt["major_iteration"] < first["major_iterations"]
        assert attempt["start_index"] < first["lower_starts_per_iteration"]
        assert attempt["start_source"] in {
            "continuation",
            "cut_state",
            "stage_i_witness",
            "homogeneous_reference",
            "boundary_aware_seed",
            "sobol",
            "direct_l",
        }
        assert len(attempt["internal_start"]) == 2
        assert len(attempt["physical_start_modified_fractions"]) == len(CHARGES) - 1
        assert attempt["physical_start_volume"] > 0.0
        assert attempt["provider_status"] == "manufactured_oracle"
        assert attempt["same_major_upper_bound"] == pytest.approx(
            first["bound_history"][attempt["major_iteration"]]["upper_bound"]
        )
        assert attempt["same_major_multipliers"] == pytest.approx(
            first["bound_history"][attempt["major_iteration"]]["multipliers"]
        )

        if attempt["solver_converged"]:
            assert attempt["solver_status"] in {"success", "acceptable_point"}
            assert len(attempt["internal_terminal"]) == 2
            assert len(attempt["terminal_modified_fractions"]) == len(CHARGES) - 1
            assert attempt["terminal_volume"] > 0.0
            assert len(attempt["lower_bound_multipliers"]) == 2
            assert len(attempt["upper_bound_multipliers"]) == 2
            assert attempt["chart_jacobian_condition"] == pytest.approx(1.0)
            assert attempt["dual_pullback_inf_norm"] == pytest.approx(0.0)
            assert attempt["chart_kkt_inf_norm"] >= 0.0
            assert attempt["physical_kkt_inf_norm"] >= 0.0
            assert attempt["complementarity_inf_norm"] >= 0.0
            assert attempt["pressure_passed"] is True
            assert attempt["dual_signs_valid"] is True
            assert attempt["basin_id"] >= 0
        else:
            assert attempt["cut_eligible"] is False
            assert attempt["step6_eligible"] is False

    classification = first["attempt_classification"]
    assert classification["declared"] == len(first["attempt_trace"])
    assert classification["solver_converged"] == classification["declared"]
    assert classification["solver_failed"] == 0
    assert classification["physical_kkt_passed"] == classification["declared"]
    assert classification["step6_eligible"] > 1
    assert any(
        attempt["cut_eligible"] and not attempt["step6_eligible"]
        for attempt in first["attempt_trace"]
    )
    for major in range(first["major_iterations"]):
        certified_compositions = {
            round(attempt["terminal_modified_fractions"][1], 6)
            for attempt in first["attempt_trace"]
            if attempt["major_iteration"] == major and attempt["physical_kkt_passed"]
        }
        assert len(certified_compositions) == 3


def test_held2_stage_ii_explorer_keeps_same_composition_density_branches_distinct() -> None:
    result = _equilibrium._held2_stage_ii_basin_explorer("same_composition_different_density")

    assert result["outcome"] == "representatives_found"
    assert result["failed_evaluation_count"] == 0
    assert len(result["representatives"]) == 2
    assert {
        tuple(start["independent_modified_fractions"]) for start in result["representatives"]
    } == {(0.5,)}
    assert len({round(start["log_volume"], 8) for start in result["representatives"]}) == 2
    assert all(start["root_completeness"] == "not_proven" for start in result["representatives"])


def test_held2_stage_ii_explorer_keeps_different_compositions_at_same_density() -> None:
    result = _equilibrium._held2_stage_ii_basin_explorer("different_composition_same_density")

    assert result["outcome"] == "representatives_found"
    assert len(result["representatives"]) == 2
    assert sorted(
        start["independent_modified_fractions"][0] for start in result["representatives"]
    ) == pytest.approx([0.25, 0.75])
    assert len({round(start["log_volume"], 8) for start in result["representatives"]}) == 1


def test_held2_stage_ii_explorer_retains_tied_stable_density_branches() -> None:
    result = _equilibrium._held2_stage_ii_basin_explorer("tied_density_branches")

    assert result["outcome"] == "representatives_found"
    assert len(result["representatives"]) == 2
    assert result["evaluations"][0]["pressure_envelope"]["failure_reason"] == "stable_objective_tie"


def test_held2_stage_ii_explorer_deduplicates_physical_starts() -> None:
    result = _equilibrium._held2_stage_ii_basin_explorer("duplicates")

    assert result["outcome"] == "representatives_found"
    assert len(result["representatives"]) == 1
    assert result["duplicate_start_count"] >= 2


def test_held2_stage_ii_explorer_fails_closed_on_provider_or_root_failure() -> None:
    result = _equilibrium._held2_stage_ii_basin_explorer("provider_failure")

    assert result["outcome"] == "indeterminate"
    assert result["termination_reason"] == "required_envelope_evaluation_failed"
    assert result["completed_evaluation_count"] == 0
    assert result["failed_evaluation_count"] == 1
    assert result["evaluations"][0]["failure_reason"] == "evaluation_failed"
    assert result["representatives"] == []
    assert result["globality_certificate"] == "not_guaranteed"


def test_held2_stage_ii_direct_escalation_is_deterministic_and_bounded() -> None:
    first = _equilibrium._held2_stage_ii_basin_explorer("direct", True, 16)
    second = _equilibrium._held2_stage_ii_basin_explorer("direct", True, 16)

    assert first == second
    assert first["outcome"] == "representatives_found"
    assert first["termination_reason"] == "declared_direct_budget_exhausted"
    assert first["direct_escalation_used"] is True
    assert first["direct_solver"] == "nlopt_gn_direct_l"
    assert first["direct_solver_version"] == "2.11.0"
    assert first["declared_direct_budget"] == 16
    assert first["globality_certificate"] == "not_guaranteed"


def test_held2_stage_ii_explorer_uses_exact_installed_provider_envelopes() -> None:
    model = _figiel_brine_model()
    result = _equilibrium._held2_stage_ii_basin_explorer(
        epcsaft.native_sdk(model),
        298.15,
        100_000.0,
        ((0.02,),),
        model.parameter_fingerprint,
        0,
    )

    assert result["outcome"] == "representatives_found"
    assert result["failed_evaluation_count"] == 0
    assert result["completed_evaluation_count"] == 1
    assert result["parameter_fingerprint"] == model.parameter_fingerprint
    assert len(result["representatives"]) == 2
    assert all(
        start["source"] == "external_seed" and start["root_completeness"] == "not_proven"
        for start in result["representatives"]
    )
    assert sorted(start["log_volume"] for start in result["representatives"]) == pytest.approx(
        [-10.929425447212154, -3.78075692619037],
        abs=1.0e-8,
    )


def test_held2_stage_ii_stall_exhaustion_is_indeterminate() -> None:
    result = _equilibrium._held2_adapter(
        CHARGES,
        (0.9, 0.05, 0.05),
        "stage_ii",
    )

    assert result["outcome"] == "indeterminate_finite_search_stalled"
    assert result["direct_escalation_used"] is True
    assert result["globality_certificate"] == "not_guaranteed"


def test_held2_stage_ii_highs_upper_lp_matches_analytic_envelope() -> None:
    result = _equilibrium._held2_stage_ii_upper_lp(
        (1.0, 2.0),
        ((1.0,), (-1.0,)),
    )
    oracle = _equilibrium._held2_stage_ii_upper_lp(
        (1.0, 2.0),
        ((1.0,), (-1.0,)),
        -math.inf,
        "analytic_1d_test_oracle",
    )

    assert result["outcome"] == "optimal"
    assert result["solver"] == "highs_lp"
    assert result["solver_version"] == "1.15.1"
    assert result["upper_bound"] == pytest.approx(1.5)
    assert result["multipliers"] == pytest.approx((0.5,))
    assert result["upper_bound"] == pytest.approx(oracle["upper_bound"])
    assert result["multipliers"] == pytest.approx(oracle["multipliers"])
    assert result["primal_feasible"] is True
    assert result["dual_feasible"] is True
    assert result["primal_residual_inf"] <= 1.0e-8
    assert result["dual_residual_inf"] <= 1.0e-8
    assert result["cut_slacks"] == pytest.approx((0.0, 0.0), abs=1.0e-8)
    assert result["active_cut_ids"] == [0, 1]


def test_held2_stage_ii_highs_upper_lp_retains_tied_and_redundant_cuts() -> None:
    result = _equilibrium._held2_stage_ii_upper_lp(
        (1.0, 1.0, 10.0),
        ((1.0,), (1.0,), (-1.0,)),
    )

    assert result["outcome"] == "optimal"
    assert result["upper_bound"] == pytest.approx(5.5)
    assert result["multipliers"] == pytest.approx((4.5,))
    assert result["cut_slacks"] == pytest.approx((0.0, 0.0, 0.0), abs=1.0e-8)
    assert result["active_cut_ids"] == [0, 1, 2]
    assert len(result["cut_duals"]) == 3


def test_held2_stage_ii_highs_upper_lp_handles_nearly_parallel_cuts() -> None:
    result = _equilibrium._held2_stage_ii_upper_lp(
        (1.0, 1.0),
        ((1.0,), (-1.0 + 1.0e-8,)),
    )

    assert result["outcome"] == "optimal"
    assert result["upper_bound"] == pytest.approx(1.0, abs=1.0e-8)
    assert result["multipliers"] == pytest.approx((0.0,), abs=1.0e-8)
    assert result["primal_feasible"] is True
    assert result["dual_feasible"] is True


@pytest.mark.parametrize(
    ("intercepts", "slopes", "value_lower_bound", "expected"),
    [
        ((1.0,), ((0.0,),), 2.0, "infeasible"),
        ((), (), -math.inf, "unbounded"),
    ],
)
def test_held2_stage_ii_highs_upper_lp_fails_closed(
    intercepts: tuple[float, ...],
    slopes: tuple[tuple[float, ...], ...],
    value_lower_bound: float,
    expected: str,
) -> None:
    result = _equilibrium._held2_stage_ii_upper_lp(
        intercepts,
        slopes,
        value_lower_bound,
    )

    assert result["outcome"] == expected
    assert result["primal_feasible"] is False
    assert result["dual_feasible"] is False
    assert result["upper_bound"] is None
    assert result["multipliers"] == []


def test_held2_pressure_envelope_classifies_three_root_topology() -> None:
    result = _equilibrium._held2_pressure_envelope("three_root", 0.25, 64)

    assert result["outcome"] == "selected"
    assert result["root_completeness"] == "not_proven"
    assert result["selection_scope"] == "lowest_among_detected_strict_stable_roots"
    assert [root["mechanical_class"] for root in result["roots"]] == [
        "strict_stable",
        "unstable",
        "strict_stable",
    ]
    assert [root["log_volume"] for root in result["roots"]] == pytest.approx(
        [-1.0, -0.2, 1.0], abs=1.0e-8
    )
    assert result["selected_root_index"] in {0, 2}
    assert len(result["scan_points"]) > 64
    assert len(result["intervals"]) >= 64
    assert all(interval["status"] for interval in result["intervals"])
    assert result["intervals"][0]["lower_log_volume"] == pytest.approx(-1.5)
    assert result["intervals"][-1]["upper_log_volume"] == pytest.approx(1.5)
    assert all(
        left["upper_log_volume"] == pytest.approx(right["lower_log_volume"])
        for left, right in pairwise(result["intervals"])
    )


@pytest.mark.parametrize(
    ("topology", "outcome", "failure_reason"),
    [
        ("one_root", "selected", ""),
        ("tangential", "indeterminate", "marginal_root"),
        ("boundary", "indeterminate", "boundary_root"),
        ("invalid", "indeterminate", "evaluation_failed"),
        ("tied", "indeterminate", "stable_objective_tie"),
    ],
)
def test_held2_pressure_envelope_fail_closed_topologies(
    topology: str,
    outcome: str,
    failure_reason: str,
) -> None:
    result = _equilibrium._held2_pressure_envelope(topology, 0.25, 64)

    assert result["outcome"] == outcome
    assert result["failure_reason"] == failure_reason
    assert result["root_completeness"] == "not_proven"
    if topology == "invalid":
        assert result["evaluation_failure_count"] > 0
        assert any(not point["valid"] for point in result["scan_points"])
    if topology == "tangential":
        assert result["tangential_root_count"] == 1
        assert result["marginal_root_count"] == 1
    if topology == "boundary":
        assert result["boundary_root_count"] == 1
    if topology == "tied":
        assert result["objective_tie_count"] == 1


def test_held2_pressure_envelope_detects_close_roots_and_deduplicates_nodes() -> None:
    close = _equilibrium._held2_pressure_envelope("close_roots", 0.25, 256)
    duplicate = _equilibrium._held2_pressure_envelope("node_root", 0.25, 64)

    assert close["outcome"] == "indeterminate"
    assert close["failure_reason"] == "stable_objective_tie"
    assert [root["log_volume"] for root in close["roots"]] == pytest.approx(
        [-0.02, 0.0, 0.02], abs=1.0e-8
    )
    assert duplicate["outcome"] == "selected"
    assert len(duplicate["roots"]) == 1
    assert duplicate["deduplicated_root_count"] >= 1


def test_held2_pressure_envelope_records_branch_switch() -> None:
    left = _equilibrium._held2_pressure_envelope("branch_switch", 0.25, 64)
    right = _equilibrium._held2_pressure_envelope("branch_switch", 0.75, 64)

    assert left["outcome"] == "selected"
    assert right["outcome"] == "selected"
    assert left["roots"][left["selected_root_index"]]["log_volume"] < 0.0
    assert right["roots"][right["selected_root_index"]]["log_volume"] > 0.0


def test_held2_pressure_envelope_uses_provider_bounds_and_exact_derivatives() -> None:
    model = _figiel_brine_model()
    result = _equilibrium._held2_pressure_envelope(
        epcsaft.native_sdk(model),
        298.15,
        100_000.0,
        (0.02,),
        model.parameter_fingerprint,
        64,
    )

    assert result["root_completeness"] == "not_proven"
    assert result["globality_certificate"] == "not_guaranteed"
    assert result["molar_volume_bounds"][0] > 0.0
    assert result["molar_volume_bounds"][1] > result["molar_volume_bounds"][0]
    assert result["evaluation_failure_count"] == 0
    assert all(point["valid"] for point in result["scan_points"])
    assert result["outcome"] == "selected"
    assert result["selected_root_index"] == 0
    assert [root["mechanical_class"] for root in result["roots"]] == [
        "strict_stable",
        "unstable",
        "strict_stable",
    ]
    assert result["roots"][0]["objective"] < result["roots"][2]["objective"]
    assert all(abs(root["pressure_residual"]) <= 1.0e-8 for root in result["roots"])
    assert all(
        root["mechanical_class"] in {"strict_stable", "unstable", "marginal"}
        for root in result["roots"]
    )
    for root in result["roots"]:
        if root["mechanical_class"] == "strict_stable":
            assert root["pressure_derivative_log_volume"] < 0.0
            assert root["objective_curvature_log_volume"] > 0.0
        elif root["mechanical_class"] == "unstable":
            assert root["pressure_derivative_log_volume"] > 0.0
            assert root["objective_curvature_log_volume"] < 0.0

    log_volume = result["roots"][0]["log_volume"]
    step = 1.0e-5
    center = _equilibrium._held2_adapter(
        epcsaft.native_sdk(model),
        298.15,
        100_000.0,
        (0.02,),
        log_volume,
        model.parameter_fingerprint,
    )
    lower = _equilibrium._held2_adapter(
        epcsaft.native_sdk(model),
        298.15,
        100_000.0,
        (0.02,),
        log_volume - step,
        model.parameter_fingerprint,
    )
    upper = _equilibrium._held2_adapter(
        epcsaft.native_sdk(model),
        298.15,
        100_000.0,
        (0.02,),
        log_volume + step,
        model.parameter_fingerprint,
    )
    finite_difference = (
        upper["pressure_stationarity_relative"] - lower["pressure_stationarity_relative"]
    ) / (2.0 * step)
    assert center["pressure_stationarity_derivative_log_volume"] == pytest.approx(
        finite_difference,
        rel=2.0e-8,
        abs=2.0e-8,
    )


def test_held2_general_mp_stage_iii_exact_lagrangian_hessian() -> None:
    candidates = ((0.2, 1.0), (0.20000002, 1.0), (0.8, 1.0))
    center = (0.25, 0.2, 0.0, 0.25, 0.21, 0.0, 0.5, 0.795, 0.0)
    multipliers = (0.3, -0.2, 0.0, 0.0, 0.0)
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
    assert result["retired_inactive_count"] == 0
    assert result["stage_iii_solve_count"] == 2
    assert result["active_set_resolve_count"] == 1
    assert [step["action"] for step in result["lifecycle"] if step["action"] != "retain_phase"] == [
        "merge_duplicate",
        "accept_active_set",
    ]
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
    assert result["dual_sign_violation_inf_norm"] < 1.0e-8
    assert result["bound_complementarity_inf_norm"] < 1.0e-8
    assert result["minimum_phase_fraction"] > 1.0e-10
    assert abs(result["enumeration_objective_gap"]) < 1.0e-9


def test_held2_stage_iii_retires_only_kkt_inactive_phases_one_at_a_time() -> None:
    active_trace = _equilibrium._held2_stage_iii_retirement_decision(
        5.0e-11,
        0.0,
        0.0,
        -1.0e-4,
        True,
    )
    inactive = _equilibrium._held2_stage_iii_retirement_decision(
        5.0e-11,
        0.25,
        0.0,
        0.25,
        True,
    )
    balance_required = _equilibrium._held2_stage_iii_retirement_decision(
        5.0e-11,
        0.25,
        0.0,
        0.25,
        False,
    )

    assert active_trace["retire"] is False
    assert active_trace["reason"] == "descent_or_marginal_phase"
    assert inactive["retire"] is True
    assert inactive["reason"] == "kkt_inactive"
    assert balance_required["retire"] is False
    assert balance_required["reason"] == "remaining_balance_infeasible"

    result = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        ((0.2, 1.0), (0.35, 1.0), (0.65, 1.0), (0.8, 1.0)),
        "stage_iii",
    )

    assert result["solver_status"] in {"solve_succeeded", "solved_to_acceptable_level"}
    assert result["numerical_status"] == "converged"
    assert result["physical_status"] == "accepted"
    assert result["retired_inactive_count"] == 2
    assert result["stage_iii_solve_count"] == 3
    assert result["active_set_resolve_count"] == 2
    state_changes = [step for step in result["lifecycle"] if step["action"] != "retain_phase"]
    assert [step["active_candidate_count"] for step in state_changes] == [4, 3, 2]
    assert [step["action"] for step in state_changes] == [
        "retire_kkt_inactive",
        "retire_kkt_inactive",
        "accept_active_set",
    ]
    assert all(
        step["solver_status"] in {"solve_succeeded", "solved_to_acceptable_level"}
        for step in result["lifecycle"]
    )
    assert [step["decision_reason"] for step in state_changes] == [
        "kkt_inactive",
        "kkt_inactive",
        "active_set_certified",
    ]
    assert [
        step["candidate_independent_modified_fractions"][0] for step in state_changes[:-1]
    ] == pytest.approx([0.35, 0.65])
    retained = [step for step in result["lifecycle"] if step["action"] == "retain_phase"]
    assert len(retained) == 4
    assert all(step["decision_reason"] == "phase_amount_active" for step in retained)

    trace = _equilibrium._held2_adapter(
        CHARGES,
        (0.999999998, 1.0e-9, 1.0e-9),
        ((2.0e-10, 1.0), (0.8, 1.0)),
        "stage_iii",
    )
    assert trace["solver_status"] in {"solve_succeeded", "solved_to_acceptable_level"}
    assert trace["numerical_status"] == "converged"
    assert trace["physical_status"] == "not_adjudicated"
    assert trace["feedback"] == "return_to_stage_ii"
    assert trace["failure_reason"] == "trace_component_requires_log_refinement"
    assert trace["trace_refinement_status"] == "complementarity_refinement_required"


def test_held2_stage_iii_ipopt_success_does_not_override_physical_potential_failure() -> None:
    result = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        ((0.1, 1.0), (0.9, 1.0)),
        "stage_iii",
    )

    assert result["solver_status"] in {"solve_succeeded", "solved_to_acceptable_level"}
    assert result["numerical_status"] == "converged"
    assert result["physical_status"] == "rejected"
    assert result["feedback"] == "return_to_stage_ii"
    assert result["failure_reason"] == "stage_iii_physical_certificate_failed"
    assert result["modified_potential_mixed_gap"] > 1.0e-8


def test_held2_general_mp_stage_iii_returns_infeasible_set_to_stage_ii() -> None:
    result = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        ((0.1, 1.0), (0.2, 1.0), (0.3, 1.0)),
        "stage_iii",
    )

    assert result["solver_status"] == "infeasible_problem_detected"
    assert result["numerical_status"] == "not_converged"
    assert result["physical_status"] == "not_adjudicated"
    assert result["feedback"] == "return_to_stage_ii"
    assert result["failure_reason"] == "stage_iii_solver_not_converged"
    assert result["globality_certificate"] == "not_guaranteed"
    assert result["lifecycle"] == [
        {
            "solve_index": 1,
            "active_candidate_count": 3,
            "removed_candidate_index": -1,
            "action": "solve_failed",
            "phase_fraction": 0.0,
            "lower_bound_multiplier": 0.0,
            "reduced_derivative": 0.0,
            "complementarity_inf_norm": 0.0,
            "candidate_independent_modified_fractions": [],
            "candidate_volume": 0.0,
            "solver_status": "infeasible_problem_detected",
            "decision_reason": "stage_iii_solver_not_converged",
        }
    ]
