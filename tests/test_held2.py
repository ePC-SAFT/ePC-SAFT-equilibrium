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


def test_modified_mole_coordinates_preserve_charge_and_galvani_gauge() -> None:
    chemical_potentials = (3.0, -2.0, 4.0)
    baseline = _equilibrium._held2_coordinates(
        CHARGES,
        PHYSICAL_FEED,
        chemical_potentials,
    )
    shifted = _equilibrium._held2_coordinates(
        CHARGES,
        PHYSICAL_FEED,
        tuple(
            potential + 17.0 * charge
            for potential, charge in zip(
                chemical_potentials, CHARGES, strict=True
            )
        ),
    )

    assert baseline["eliminated_index"] == 2
    assert baseline["independent_indices"] == [1]
    assert baseline["modified_factors"] == pytest.approx([1.0, 2.0])
    assert baseline["modified_feed"] == pytest.approx([0.5, 0.5])
    assert baseline["lifted_feed"] == pytest.approx(PHYSICAL_FEED, abs=1.0e-15)
    assert shifted["transformed_modified_potentials"] == pytest.approx(
        baseline["transformed_modified_potentials"], abs=1.0e-12
    )


@pytest.mark.parametrize(
    ("cube", "ceiling", "expected"),
    [
        ((0.0,), 0.38, (2.0e-10,)),
        ((0.5,), 0.38, (0.1900000001,)),
        ((1.0,), 0.38, (0.38,)),
        ((0.5,), math.nan, (0.5,)),
    ],
)
def test_stage_i_domain_map_enforces_only_the_provider_ionic_ceiling(
    cube: tuple[float, ...],
    ceiling: float,
    expected: tuple[float, ...],
) -> None:
    result = _equilibrium._held2_stage_i_domain(CHARGES, cube, ceiling)

    assert result["independent_modified_fractions"] == pytest.approx(expected)
    assert result["physical_total_ion_mole_fraction"] == pytest.approx(
        expected[0]
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


def test_canonical_manufactured_workflow_certifies_all_ten_steps() -> None:
    result = _equilibrium._held2_manufactured_controller(
        CHARGES,
        PHYSICAL_FEED,
    )

    assert result["controller"] == "perdomo_held2_steps_1_to_10_v1"
    assert result["outcome"] == "physical_equilibrium_accepted"
    assert result["next_action"] == "accept_multiphase"
    assert result["completed_step"] == 10
    assert [
        transition["completed_step"] for transition in result["transitions"]
    ] == [1, 3, 7, 10]
    assert result["stage_i"]["outcome"] == "negative_witness_found"

    stage_ii = result["stage_ii"]
    assert stage_ii["outcome"] == "candidate_set"
    assert stage_ii["local_solver"] == "ipopt_exact_hessian"
    assert stage_ii["direct_escalation_used"] is False
    assert sorted(
        candidate["modified_fractions"][1]
        for candidate in stage_ii["candidates"]
    ) == pytest.approx([0.2, 0.8], abs=2.0e-7)
    for context in stage_ii["major_contexts"]:
        attempts = [
            attempt
            for attempt in stage_ii["attempt_trace"]
            if attempt["major_iteration"] == context["major_id"]
        ]
        assert context["lower_attempt_ids"] == [
            attempt["attempt_id"] for attempt in attempts
        ]
        assert context["step6_eligible_attempt_ids"] == [
            attempt["attempt_id"]
            for attempt in attempts
            if attempt["step6_eligible"]
        ]

    stage_iii = result["stage_iii"]
    assert stage_iii["physical_status"] == "accepted"
    assert stage_iii["modified_balance_inf_norm"] < 1.0e-9
    assert stage_iii["ordinary_balance_inf_norm"] < 1.0e-9
    assert stage_iii["pressure_stationarity_inf_norm"] < 1.0e-9
    assert stage_iii["modified_potential_mixed_gap"] < 1.0e-9


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
