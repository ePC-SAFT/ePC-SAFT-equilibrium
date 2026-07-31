from __future__ import annotations

import pytest

from epcsaft_equilibrium import _equilibrium


def _farkas(
    matrix: list[list[float]],
    row_lower: list[float],
    row_upper: list[float],
    column_lower: list[float],
    column_upper: list[float],
    ray: list[float],
) -> dict[str, object]:
    return _equilibrium._held2_audit_farkas_certificate(
        matrix,
        row_lower,
        row_upper,
        column_lower,
        column_upper,
        ray,
    )


@pytest.mark.parametrize(
    ("matrix", "row_lower", "row_upper", "column_lower", "column_upper", "ray"),
    (
        ([[1.0]], [2.0], [float("inf")], [0.0], [1.0], [-1.0]),
        ([[1.0]], [float("-inf")], [-1.0], [0.0], [1.0], [1.0]),
        (
            [[1.0], [1.0]],
            [float("-inf"), 1.0],
            [0.0, float("inf")],
            [0.0],
            [1.0],
            [1.0, -1.0],
        ),
        (
            [[3.0], [2.0]],
            [3.0, float("-inf")],
            [float("inf"), 0.0],
            [0.0],
            [1.0],
            [-1.0 / 3.0, 0.5],
        ),
    ),
)
def test_farkas_audit_certifies_infeasibility_across_bound_and_scaling_conventions(
    matrix: list[list[float]],
    row_lower: list[float],
    row_upper: list[float],
    column_lower: list[float],
    column_upper: list[float],
    ray: list[float],
) -> None:
    audit = _farkas(
        matrix,
        row_lower,
        row_upper,
        column_lower,
        column_upper,
        ray,
    )
    assert audit["accepted"] is True
    assert audit["contradiction_margin"] > audit["contradiction_threshold"]


@pytest.mark.parametrize(
    ("row_lower", "ray", "reason"),
    (
        ([0.5], [-1.0], "contradiction_margin_failed"),
        ([1.0 + 5.0e-10], [-1.0], "contradiction_margin_failed"),
        ([2.0], [1.0], "row_sign_failed"),
        ([2.0], [], "dimension_mismatch"),
        ([2.0], [float("nan")], "nonfinite_evidence"),
    ),
)
def test_farkas_audit_rejects_feasible_marginal_or_malformed_evidence(
    row_lower: list[float],
    ray: list[float],
    reason: str,
) -> None:
    audit = _farkas(
        [[1.0]],
        row_lower,
        [float("inf")],
        [0.0],
        [1.0],
        ray,
    )
    assert audit["accepted"] is False
    assert audit["reason"] == reason


def test_farkas_audit_is_invariant_to_row_permutation() -> None:
    first = _farkas(
        [[1.0], [1.0]],
        [float("-inf"), 1.0],
        [0.0, float("inf")],
        [0.0],
        [1.0],
        [1.0, -1.0],
    )
    second = _farkas(
        [[1.0], [1.0]],
        [1.0, float("-inf")],
        [float("inf"), 0.0],
        [0.0],
        [1.0],
        [-1.0, 1.0],
    )
    assert first["accepted"] is second["accepted"] is True
    assert first["contradiction_margin"] == pytest.approx(
        second["contradiction_margin"]
    )


def test_farkas_audit_is_invariant_to_column_permutation() -> None:
    first = _farkas(
        [[1.0, 2.0]],
        [1.0],
        [float("inf")],
        [float("-inf"), 0.0],
        [0.0, 0.0],
        [-1.0],
    )
    second = _farkas(
        [[2.0, 1.0]],
        [1.0],
        [float("inf")],
        [0.0, float("-inf")],
        [0.0, 0.0],
        [-1.0],
    )
    assert first["accepted"] is second["accepted"] is True
    assert first["contradiction_margin"] == pytest.approx(
        second["contradiction_margin"]
    )


@pytest.mark.parametrize(
    ("solver_infeasible", "has_certificate", "accepted", "expected"),
    (
        (True, True, True, "certified_infeasible"),
        (True, True, False, "indeterminate"),
        (True, False, False, "indeterminate"),
        (False, True, True, "indeterminate"),
    ),
)
def test_farkas_status_and_certificate_must_agree(
    solver_infeasible: bool,
    has_certificate: bool,
    accepted: bool,
    expected: str,
) -> None:
    assert (
        _equilibrium._held2_adjudicate_farkas_status(
            solver_infeasible, has_certificate, accepted
        )
        == expected
    )


@pytest.mark.parametrize(
    ("column_lower", "column_upper", "matrix", "row_lower", "row_upper", "ray"),
    (
        ([0.0], [float("inf")], [[1.0]], [float("-inf")], [-1.0], [1.0]),
        ([float("-inf")], [0.0], [[1.0]], [1.0], [float("inf")], [-1.0]),
        ([0.0], [0.0], [[1.0]], [1.0], [float("inf")], [-1.0]),
        (
            [float("-inf")],
            [float("inf")],
            [[1.0], [1.0]],
            [float("-inf"), 1.0],
            [0.0, float("inf")],
            [1.0, -1.0],
        ),
    ),
)
def test_farkas_audit_covers_one_sided_fixed_and_free_columns(
    column_lower: list[float],
    column_upper: list[float],
    matrix: list[list[float]],
    row_lower: list[float],
    row_upper: list[float],
    ray: list[float],
) -> None:
    audit = _farkas(
        matrix,
        row_lower,
        row_upper,
        column_lower,
        column_upper,
        ray,
    )
    assert audit["accepted"] is True


@pytest.mark.parametrize(
    ("row_lower", "row_upper", "expected_reason"),
    (
        ([float("inf")], [float("inf")], "nonfinite_evidence"),
        ([float("-inf")], [float("-inf")], "nonfinite_evidence"),
    ),
)
def test_farkas_audit_rejects_same_sign_infinite_bound_pairs(
    row_lower: list[float],
    row_upper: list[float],
    expected_reason: str,
) -> None:
    audit = _farkas(
        [[1.0]],
        row_lower,
        row_upper,
        [0.0],
        [1.0],
        [1.0],
    )
    assert audit["accepted"] is False
    assert audit["reason"] == expected_reason


def test_farkas_audit_rejects_overflowed_derived_evidence() -> None:
    audit = _farkas(
        [[1.0e308], [1.0e308]],
        [float("-inf"), float("-inf")],
        [-1.0, -1.0],
        [0.0],
        [1.0],
        [1.0, 1.0],
    )
    assert audit["accepted"] is False
    assert audit["reason"] == "nonfinite_derived_evidence"


def test_farkas_audit_rejects_tiny_wrong_sign_on_unbounded_column_side() -> None:
    audit = _farkas(
        [[1.0e-11]],
        [float("-inf")],
        [-1.0],
        [float("-inf")],
        [-1.0e12],
        [1.0],
    )
    assert audit["accepted"] is False
    assert audit["reason"] == "dual_feasibility_failed"


def test_farkas_audit_rejects_overflowed_contradiction_margin() -> None:
    audit = _farkas(
        [[1.0]],
        [float("-inf")],
        [-1.0e308],
        [1.0e308],
        [float("inf")],
        [1.0],
    )
    assert audit["accepted"] is False
    assert audit["reason"] == "nonfinite_derived_evidence"


def _kkt(
    *,
    variables: list[float],
    lower: list[float],
    upper: list[float],
    gradient: list[float],
    constraints: list[list[float]] | None = None,
    constraint_upper: list[float] | None = None,
    z_lower: list[float] | None = None,
    z_upper: list[float] | None = None,
    multipliers: list[float] | None = None,
    pullback_residual: float = 0.0,
    pullback_scale: float = 1.0,
    pressure_residual: float = 0.0,
    same_major: bool = True,
    step4_binding: bool = True,
    pressure_branch: bool = True,
) -> dict[str, object]:
    return _equilibrium._held2_audit_step5_kkt(
        variables,
        lower,
        upper,
        gradient,
        constraints or [],
        constraint_upper or [],
        z_lower or [0.0] * len(variables),
        z_upper or [0.0] * len(variables),
        multipliers or [],
        pullback_residual,
        pullback_scale,
        pressure_residual,
        same_major,
        step4_binding,
        pressure_branch,
    )


@pytest.mark.parametrize(
    "audit",
    (
        _kkt(
            variables=[0.5, 2.0],
            lower=[0.0, 1.0],
            upper=[1.0, 3.0],
            gradient=[0.0, 0.0],
        ),
        _kkt(
            variables=[0.0],
            lower=[0.0],
            upper=[1.0],
            gradient=[1.0],
            z_lower=[1.0],
        ),
        _kkt(
            variables=[1.0],
            lower=[0.0],
            upper=[1.0],
            gradient=[-1.0],
            z_upper=[1.0],
        ),
        _kkt(
            variables=[0.0],
            lower=[-1.0],
            upper=[1.0],
            gradient=[-1.0],
            constraints=[[1.0]],
            constraint_upper=[0.0],
            multipliers=[1.0],
        ),
    ),
)
def test_step5_kkt_audit_accepts_interior_and_active_bound_terminals(
    audit: dict[str, object],
) -> None:
    assert audit["accepted"] is True
    assert audit["reason"] == "certified"


@pytest.mark.parametrize(
    ("audit", "reason"),
    (
        (
            _kkt(
                variables=[0.0],
                lower=[0.0],
                upper=[1.0],
                gradient=[1.0],
                z_lower=[0.5],
            ),
            "stationarity_failed",
        ),
        (
            _kkt(
                variables=[0.0],
                lower=[-1.0],
                upper=[1.0],
                gradient=[-1.0],
                constraints=[[1.0], [2.0]],
                constraint_upper=[0.0, 0.0],
                multipliers=[0.5, 0.25],
            ),
            "active_jacobian_rank_deficient",
        ),
        (
            _kkt(
                variables=[0.5],
                lower=[0.0],
                upper=[1.0],
                gradient=[0.0],
                pullback_residual=1.0e-5,
            ),
            "pullback_failed",
        ),
        (
            _kkt(
                variables=[0.5],
                lower=[0.0],
                upper=[1.0],
                gradient=[0.0],
                same_major=False,
            ),
            "stale_major_iteration",
        ),
        (
            _kkt(
                variables=[0.5],
                lower=[0.0],
                upper=[1.0],
                gradient=[0.0],
                step4_binding=False,
            ),
            "step4_binding_failed",
        ),
        (
            _kkt(
                variables=[0.5],
                lower=[0.0],
                upper=[1.0],
                gradient=[0.0],
                pressure_branch=False,
            ),
            "pressure_branch_failed",
        ),
        (
            _kkt(
                variables=[0.5],
                lower=[0.0, 0.0],
                upper=[1.0],
                gradient=[0.0],
            ),
            "dimension_mismatch",
        ),
        (
            _kkt(
                variables=[float("nan")],
                lower=[0.0],
                upper=[1.0],
                gradient=[0.0],
            ),
            "nonfinite_evidence",
        ),
    ),
)
def test_step5_kkt_audit_fails_closed_on_corrupt_or_stale_evidence(
    audit: dict[str, object],
    reason: str,
) -> None:
    assert audit["accepted"] is False
    assert audit["reason"] == reason


def test_step5_kkt_audit_is_invariant_to_coordinate_permutation() -> None:
    first = _kkt(
        variables=[0.0, 0.5],
        lower=[0.0, 0.0],
        upper=[1.0, 1.0],
        gradient=[1.0, 0.0],
        z_lower=[1.0, 0.0],
    )
    second = _kkt(
        variables=[0.5, 0.0],
        lower=[0.0, 0.0],
        upper=[1.0, 1.0],
        gradient=[0.0, 1.0],
        z_lower=[0.0, 1.0],
    )
    assert first["accepted"] is second["accepted"] is True
    assert first["stationarity_residual_inf"] == pytest.approx(
        second["stationarity_residual_inf"]
    )
