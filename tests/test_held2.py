from __future__ import annotations

import math

import epcsaft
import pytest

from epcsaft_equilibrium import _equilibrium

CHARGES = (0.0, 1.0, -1.0)
PHYSICAL_FEED = (0.5, 0.25, 0.25)


def _figiel_brine_model() -> epcsaft.EPCSAFT:
    parameters = epcsaft.ParameterBundle.from_catalog(
        "figiel-2025-reference-electrolytes", version=1
    ).select(("water", "sodium-cation", "chloride-anion"))
    return epcsaft.EPCSAFT(parameters)


@pytest.mark.parametrize(
    (
        "charges",
        "physical_feed",
        "expected_eliminated",
        "expected_factors",
        "expected_modified_feed",
    ),
    [
        (
            (0.0, 1.0, -1.0),
            (0.50, 0.25, 0.25),
            2,
            (2.0, 1.0),
            (0.50, 0.50),
        ),
        (
            (0.0, -1.0, 2.0),
            (0.85, 0.10, 0.05),
            2,
            (1.5, 1.0),
            (0.15, 0.85),
        ),
    ],
)
def test_step1_coordinates_round_trip_and_preserve_galvani_gauge(
    charges: tuple[float, ...],
    physical_feed: tuple[float, ...],
    expected_eliminated: int,
    expected_factors: tuple[float, ...],
    expected_modified_feed: tuple[float, ...],
) -> None:
    chemical_potentials = tuple(
        3.0 - 2.0 * index for index in range(len(charges))
    )
    baseline = _equilibrium._held2_coordinates(
        charges,
        physical_feed,
        chemical_potentials,
    )
    shifted = _equilibrium._held2_coordinates(
        charges,
        physical_feed,
        tuple(
            potential + 17.0 * charge
            for potential, charge in zip(
                chemical_potentials, charges, strict=True
            )
        ),
    )

    assert baseline["eliminated_index"] == expected_eliminated
    assert baseline["modified_factors"] == pytest.approx(expected_factors)
    assert baseline["modified_feed"] == pytest.approx(expected_modified_feed)
    assert baseline["lifted_feed"] == pytest.approx(
        physical_feed, abs=1.0e-15
    )
    assert shifted["transformed_modified_potentials"] == pytest.approx(
        baseline["transformed_modified_potentials"], abs=1.0e-12
    )


def test_step1_rejects_singular_charge_transformation() -> None:
    with pytest.raises(ValueError, match="singular"):
        _equilibrium._held2_coordinates(
            (1.0, 1.0, -1.0, -1.0, 0.0),
            (0.10, 0.10, 0.10, 0.10, 0.60),
            (0.0, 0.0, 0.0, 0.0, 0.0),
        )


def test_installed_phase_block_has_exact_reduced_first_and_second_derivatives() -> None:
    model = _figiel_brine_model()
    capsule = epcsaft.native_sdk(model)
    center = (0.02, math.log(1.0e-3))
    direction = (0.007, -0.03)

    def evaluate(values: tuple[float, float]) -> dict[str, object]:
        return _equilibrium._held2_phase_block(
            capsule,
            298.15,
            100_000.0,
            values[0:1],
            values[1],
            model.parameter_fingerprint,
        )

    step = 2.0e-5
    lower = evaluate(
        tuple(
            value - step * delta
            for value, delta in zip(center, direction, strict=True)
        )
    )
    result = evaluate(center)
    upper = evaluate(
        tuple(
            value + step * delta
            for value, delta in zip(center, direction, strict=True)
        )
    )

    numerical = (upper["objective"] - lower["objective"]) / (2.0 * step)
    analytic = sum(
        value * delta
        for value, delta in zip(result["gradient"], direction, strict=True)
    )
    assert numerical == pytest.approx(analytic, rel=1.0e-9, abs=1.0e-10)
    for row in range(2):
        numerical = (
            upper["gradient"][row] - lower["gradient"][row]
        ) / (2.0 * step)
        analytic = sum(
            result["hessian"][2 * row + column] * direction[column]
            for column in range(2)
        )
        assert numerical == pytest.approx(analytic, rel=1.0e-8, abs=1.0e-10)


@pytest.mark.parametrize(
    ("topology", "budget", "outcome", "reason"),
    [
        ("negative", 80, "negative_witness_found", "certified_negative_tpd"),
        (
            "no_negative",
            40,
            "no_negative_witness_detected",
            "declared_budget_exhausted",
        ),
        (
            "boundary",
            40,
            "indeterminate",
            "required_envelope_evaluation_failed",
        ),
        (
            "provider_failure",
            40,
            "indeterminate",
            "required_envelope_evaluation_failed",
        ),
    ],
)
def test_stage_i_reports_only_certified_finite_search_evidence(
    topology: str,
    budget: int,
    outcome: str,
    reason: str,
) -> None:
    result = _equilibrium._held2_stage_i_direct(topology, budget)

    assert result["outcome"] == outcome
    assert result["termination_reason"] == reason
    assert result["globality_certificate"] == "not_guaranteed"
    if outcome == "negative_witness_found":
        assert result["negative_witness"]["tpd"] < -1.0e-8
        assert result["negative_witness"]["pressure_certified"] is True
    elif outcome == "no_negative_witness_detected":
        assert result["completed_evaluation_count"] == budget
        assert result["negative_witness"] is None
    else:
        assert result["failed_evaluation_count"] == 1


def test_highs_problem_64_matches_the_independent_analytic_envelope() -> None:
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
    assert result["upper_bound"] == pytest.approx(1.5)
    assert result["multipliers"] == pytest.approx((0.5,))
    assert result["upper_bound"] == pytest.approx(oracle["upper_bound"])
    assert result["multipliers"] == pytest.approx(oracle["multipliers"])
    assert result["primal_feasible"] is True
    assert result["dual_feasible"] is True
    assert result["active_cut_ids"] == [0, 1]


@pytest.mark.parametrize(
    ("topology", "outcome", "failure"),
    [
        ("one_root", "selected", ""),
        ("three_root", "selected", ""),
        ("tangential", "indeterminate", "marginal_root"),
        ("boundary", "indeterminate", "boundary_root"),
        ("invalid", "indeterminate", "evaluation_failed"),
        ("tied", "indeterminate", "stable_objective_tie"),
    ],
)
def test_pressure_envelope_classifies_or_fails_closed(
    topology: str,
    outcome: str,
    failure: str,
) -> None:
    result = _equilibrium._held2_pressure_envelope(topology, 0.25, 64)

    assert result["outcome"] == outcome
    assert result["failure_reason"] == failure
    assert result["root_completeness"] == "not_proven"
    if topology == "three_root":
        assert [
            root["mechanical_class"] for root in result["roots"]
        ] == ["strict_stable", "unstable", "strict_stable"]
        assert all(
            abs(root["pressure_residual"]) <= 1.0e-8
            for root in result["roots"]
        )


def test_installed_pressure_envelope_uses_provider_bounds_and_stability_sign() -> None:
    model = _figiel_brine_model()
    result = _equilibrium._held2_pressure_envelope(
        epcsaft.native_sdk(model),
        298.15,
        100_000.0,
        (0.02,),
        model.parameter_fingerprint,
        64,
    )

    assert result["outcome"] == "selected"
    assert result["root_completeness"] == "not_proven"
    assert result["molar_volume_bounds"][0] > 0.0
    assert [
        root["mechanical_class"] for root in result["roots"]
    ] == ["strict_stable", "unstable", "strict_stable"]
    for root in result["roots"]:
        expected_sign = (
            -1.0 if root["mechanical_class"] == "strict_stable" else 1.0
        )
        assert expected_sign * root["pressure_derivative_log_volume"] > 0.0
        assert expected_sign * root["objective_curvature_log_volume"] < 0.0


def test_installed_electrolyte_phase_binding_uses_full_component_count() -> None:
    model = _figiel_brine_model()
    capsule = epcsaft.native_sdk(model)
    envelope = _equilibrium._held2_pressure_envelope(
        capsule,
        298.15,
        100_000.0,
        (0.02,),
        model.parameter_fingerprint,
        64,
    )
    volume = math.exp(
        envelope["roots"][envelope["selected_root_index"]]["log_volume"]
    )

    result = _equilibrium.evaluate_electrolyte_phase(
        capsule,
        298.15,
        (0.98, 0.01, 0.01),
        volume,
        model.parameter_fingerprint,
    )

    assert len(result["gradient"]) == 4
    assert len(result["hessian"]) == 16
    assert len(result["chemical_potential_over_rt"]) == 3
    assert math.isfinite(result["pressure_pa"])
    assert result["parameter_fingerprint"] == model.parameter_fingerprint


def test_installed_stage_iii_has_exact_lagrangian_hessian() -> None:
    model = _figiel_brine_model()
    capsule = epcsaft.native_sdk(model)
    envelopes = [
        _equilibrium._held2_pressure_envelope(
            capsule,
            298.15,
            100_000.0,
            (composition,),
            model.parameter_fingerprint,
            64,
        )
        for composition in (0.02, 0.04)
    ]
    log_volumes = tuple(
        envelope["roots"][envelope["selected_root_index"]]["log_volume"]
        for envelope in envelopes
    )
    center = (0.5, 0.02, log_volumes[0], 0.5, 0.04, log_volumes[1])
    direction = (0.0, 0.004, -0.01, 0.0, -0.003, 0.008)
    multipliers = (0.1, -0.2, 0.0, 0.0)

    def evaluate(values: tuple[float, ...]) -> dict[str, object]:
        return _equilibrium._held2_installed_stage_iii_derivatives(
            capsule,
            298.15,
            100_000.0,
            (0.97, 0.015, 0.015),
            2,
            values,
            multipliers,
            model.parameter_fingerprint,
        )

    step = 2.0e-5
    lower = evaluate(
        tuple(
            value - step * delta
            for value, delta in zip(center, direction, strict=True)
        )
    )
    result = evaluate(center)
    upper = evaluate(
        tuple(
            value + step * delta
            for value, delta in zip(center, direction, strict=True)
        )
    )
    size = len(center)
    for row in range(size):
        numerical = (
            upper["lagrangian_gradient"][row]
            - lower["lagrangian_gradient"][row]
        ) / (2.0 * step)
        analytic = sum(
            result["lagrangian_hessian"][size * row + column]
            * direction[column]
            for column in range(size)
        )
        assert numerical == pytest.approx(
            analytic, rel=2.0e-7, abs=2.0e-8
        )
