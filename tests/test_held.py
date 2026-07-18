from __future__ import annotations

import math
import sys
from collections.abc import Callable
from decimal import Decimal, getcontext

import epcsaft
import pytest

from epcsaft_equilibrium import _equilibrium

BINARY_FINGERPRINT = "sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286"
GAS_CONSTANT_J_PER_MOL_K = 8.31446261815324

# May et al. (2015), Table 5, retained source row may2015-ch4-c2h6-001.
MAY_ROW_001_ID = "may2015-ch4-c2h6-001"
MAY_ROW_001_TEMPERATURE_K = 243.58
MAY_ROW_001_PRESSURE_PA = 3_949_000.0
MAY_ROW_001_LIQUID_X_METHANE = 0.3099
MAY_ROW_001_VAPOR_X_METHANE = 0.6664
MAY_ROW_001_FEED_X_METHANE = 0.5 * (MAY_ROW_001_LIQUID_X_METHANE + MAY_ROW_001_VAPOR_X_METHANE)

# The validation plan freezes this liquid-side one-phase inference from retained
# May row may2015-ch4-c2h6-011: 0.5957 - 2 * 0.0165 = 0.5627.
MAY_ROW_011_ID = "may2015-ch4-c2h6-011"
MAY_ROW_011_TEMPERATURE_K = 243.61
MAY_ROW_011_PRESSURE_PA = 6_691_000.0
MAY_ROW_011_LIQUID_SIDE_FEED_X_METHANE = 0.5627

STAGE_I_PROFILE = "held-stage-i-binary-v1"
STAGE_II_PROFILE = "held-stage-ii-binary-v1"
STAGE_III_PROFILE = "held-stage-iii-binary-v1"
TPD_NEGATIVE_THRESHOLD = -1.0e-8
OUTER_NUMERICAL_FACTOR = 256.0

# Provider-local Stage III regression state for retained May row 001. The
# compositions are already frozen by the reviewed fixed-two-phase design; the
# molar volumes and objective are numerical solver evidence, not source data.
MAY_ROW_001_REFINED_CANDIDATES = (
    ("liquid", 0.3025223259589743, 6.614906698837923e-05),
    ("vapor", 0.6703563353120439, 3.58287622219301e-04),
)
MAY_ROW_001_REFINED_G_BAR = 6.46180265452249


def _outer_vertex_oracle(
    cuts: list[tuple[str, float, float]],
) -> tuple[float, float, tuple[str, ...], tuple[float, ...]]:
    """Independent high-precision exhaustive oracle for finite binary fixtures."""
    getcontext().prec = 50
    decimal_cuts = [
        (identity, Decimal(repr(intercept)), Decimal(repr(slope)))
        for identity, intercept, slope in cuts
    ]
    candidates: list[tuple[Decimal, Decimal]] = []
    for left_index, (_, left_intercept, left_slope) in enumerate(decimal_cuts):
        for _, right_intercept, right_slope in decimal_cuts[left_index + 1 :]:
            slope_delta = left_slope - right_slope
            if slope_delta == 0:
                continue
            multiplier = (right_intercept - left_intercept) / slope_delta
            value = left_intercept + left_slope * multiplier
            if all(value <= intercept + slope * multiplier for _, intercept, slope in decimal_cuts):
                candidates.append((value, multiplier))
    assert candidates
    best_value = max(value for value, _ in candidates)
    scale = max(1.0, abs(float(best_value)))
    tolerance = OUTER_NUMERICAL_FACTOR * sys.float_info.epsilon * scale
    tied = sorted(
        float(multiplier)
        for value, multiplier in candidates
        if abs(float(value - best_value)) <= tolerance
    )
    chosen_multiplier = tied[0]
    active = tuple(
        identity
        for identity, intercept, slope in cuts
        if abs(float(best_value) - (intercept + slope * chosen_multiplier)) <= tolerance
    )
    return float(best_value), chosen_multiplier, active, tuple(tied)


OUTER_FINITE_FIXTURES = {
    "unique": [
        ("left", 0.0, 1.0),
        ("right", 2.0, -1.0),
    ],
    "tied": [
        ("flat", 0.0, 0.0),
        ("left", 0.0, 1.0),
        ("right", 1.0, -1.0),
    ],
    "redundant": [
        ("left", 0.0, 1.0),
        ("right", 2.0, -1.0),
        ("redundant", 10.0, 0.0),
    ],
    "nearly_parallel": [
        ("left", 0.0, 1.0),
        ("nearly_parallel_redundant", 1.0e-9, 1.0 + 1.0e-12),
        ("right", 2.0, -1.0),
    ],
}


def _binary_model() -> epcsaft.EPCSAFT:
    parameters = epcsaft.ParameterBundle.from_catalog(
        "gross-2001-methane-ethane", version=1
    ).select(("methane", "ethane"))
    model = epcsaft.EPCSAFT(parameters)
    assert model.parameter_fingerprint == BINARY_FINGERPRINT
    return model


@pytest.mark.parametrize(
    "fixture_name",
    list(OUTER_FINITE_FIXTURES),
)
def test_held_outer_envelope_matches_independent_vertex_oracle(fixture_name: str) -> None:
    cuts = OUTER_FINITE_FIXTURES[fixture_name]
    expected_value, expected_multiplier, expected_active, expected_ties = _outer_vertex_oracle(cuts)

    result = _equilibrium._held_outer_envelope(cuts)

    assert result["status"] == "finite"
    assert result["value"] == pytest.approx(expected_value, rel=0.0, abs=2.0e-12)
    assert result["multiplier"] == pytest.approx(
        expected_multiplier,
        rel=0.0,
        abs=2.0e-12,
    )
    assert tuple(result["active_cut_ids"]) == expected_active
    assert tuple(result["tied_multipliers"]) == pytest.approx(
        expected_ties,
        rel=0.0,
        abs=2.0e-12,
    )
    assert all(
        result["value"]
        <= intercept
        + slope * result["multiplier"]
        + OUTER_NUMERICAL_FACTOR
        * sys.float_info.epsilon
        * max(1.0, abs(result["value"]), abs(intercept))
        for _, intercept, slope in cuts
    )


def test_held_outer_envelope_reports_missing_opposite_endpoint_cut() -> None:
    result = _equilibrium._held_outer_envelope(
        [
            ("lower_endpoint", 0.0, 1.0),
            ("interior", 1.0, 0.25),
        ]
    )

    assert result == {
        "status": "unbounded",
        "failure_reason": "cuts do not bracket the feed with opposite endpoint slopes",
        "active_cut_ids": [],
        "tied_multipliers": [],
    }


def test_held_outer_envelope_reports_infeasible_nonfinite_cut() -> None:
    result = _equilibrium._held_outer_envelope(
        [
            ("lower_endpoint", 0.0, 1.0),
            ("impossible", -math.inf, 0.0),
            ("upper_endpoint", 0.0, -1.0),
        ]
    )

    assert result == {
        "status": "infeasible",
        "failure_reason": "a cut has no finite feasible upper value",
        "active_cut_ids": [],
        "tied_multipliers": [],
    }


def test_held_stage_ii_initial_endpoint_cuts_bracket_may_row_001_feed(
    binary_capsule: object,
) -> None:
    result = _equilibrium._held_stage_ii_initial_cuts(
        binary_capsule,
        MAY_ROW_001_TEMPERATURE_K,
        MAY_ROW_001_PRESSURE_PA,
        MAY_ROW_001_FEED_X_METHANE,
        BINARY_FINGERPRINT,
    )

    assert result["status"] == "initialized"
    assert [attempt["role"] for attempt in result["attempts"]] == [
        "lower_endpoint_liquid",
        "lower_endpoint_vapor",
        "upper_endpoint_liquid",
        "upper_endpoint_vapor",
    ]
    assert all(
        attempt["solver_status"] or attempt["callback_error"] for attempt in result["attempts"]
    )
    cuts = result["cuts"]
    assert any(cut["slope"] > 0.0 for cut in cuts)
    assert any(cut["slope"] < 0.0 for cut in cuts)
    assert {cut["endpoint"] for cut in cuts} == {"lower", "upper"}
    assert all(abs(cut["pressure_stationarity_relative"]) <= 1.0e-8 for cut in cuts)
    assert all(1.0e-5 < cut["volume_m3"] < 1.0e-1 for cut in cuts)
    assert result["outer"]["status"] == "finite"


def test_held_stage_ii_lower_exact_derivatives_at_may_row_001(
    binary_capsule: object,
    unstable_stage_i: dict[str, object],
) -> None:
    reference = unstable_stage_i["reference"]
    multiplier = reference["gradient"][0]
    x_methane = MAY_ROW_001_VAPOR_X_METHANE
    log_volume = math.log(1.0 / 15_000.0)
    direction = (0.03, -0.04)
    step = 3.0e-5

    def evaluate(scale: float) -> dict[str, object]:
        return _equilibrium._held_evaluate_lower(
            binary_capsule,
            MAY_ROW_001_TEMPERATURE_K,
            MAY_ROW_001_PRESSURE_PA,
            MAY_ROW_001_FEED_X_METHANE,
            multiplier,
            x_methane + scale * direction[0],
            log_volume + scale * direction[1],
            BINARY_FINGERPRINT,
        )

    lower = evaluate(-step)
    center = evaluate(0.0)
    upper = evaluate(step)
    state = _equilibrium._held_evaluate_state(
        binary_capsule,
        MAY_ROW_001_TEMPERATURE_K,
        MAY_ROW_001_PRESSURE_PA,
        x_methane,
        log_volume,
        BINARY_FINGERPRINT,
    )
    assert center["objective"] == pytest.approx(
        state["g_bar"] + multiplier * (MAY_ROW_001_FEED_X_METHANE - x_methane),
        rel=2.0e-15,
        abs=2.0e-15,
    )
    assert center["gradient"] == pytest.approx(
        (state["gradient"][0] - multiplier, state["gradient"][1])
    )
    assert center["hessian"] == pytest.approx(state["hessian"])
    objective_directional = (upper["objective"] - lower["objective"]) / (2.0 * step)
    assert objective_directional == pytest.approx(
        sum(value * delta for value, delta in zip(center["gradient"], direction, strict=True)),
        rel=2.0e-8,
        abs=2.0e-10,
    )
    for row in range(2):
        gradient_directional = (upper["gradient"][row] - lower["gradient"][row]) / (2.0 * step)
        hessian_directional = sum(
            center["hessian"][2 * row + column] * direction[column] for column in range(2)
        )
        assert gradient_directional == pytest.approx(
            hessian_directional,
            rel=2.0e-8,
            abs=2.0e-10,
        )


def test_held_stage_ii_lower_search_is_deterministic_and_stops_below_upper(
    binary_capsule: object,
    unstable_stage_i: dict[str, object],
) -> None:
    initialization = _equilibrium._held_stage_ii_initial_cuts(
        binary_capsule,
        MAY_ROW_001_TEMPERATURE_K,
        MAY_ROW_001_PRESSURE_PA,
        MAY_ROW_001_FEED_X_METHANE,
        BINARY_FINGERPRINT,
    )
    result = _equilibrium._held_stage_ii_lower_search(
        binary_capsule,
        MAY_ROW_001_TEMPERATURE_K,
        MAY_ROW_001_PRESSURE_PA,
        MAY_ROW_001_FEED_X_METHANE,
        initialization["outer"]["multiplier"],
        unstable_stage_i["reference"]["g_bar"],
        [(cut["x_methane"], cut["log_volume"]) for cut in initialization["cuts"]],
        BINARY_FINGERPRINT,
    )

    assert result["status"] == "below_upper_found"
    assert 1 <= result["starts_completed"] <= 20
    assert result["accepted_cuts"]
    assert result["accepted_cuts"][-1]["objective"] <= result["upper_bound"]
    assert result["attempts"][-1]["accepted"] is True


def test_held_stage_ii_lower_search_retains_valid_non_tight_cuts_and_cap(
    binary_capsule: object,
) -> None:
    initialization = _equilibrium._held_stage_ii_initial_cuts(
        binary_capsule,
        MAY_ROW_001_TEMPERATURE_K,
        MAY_ROW_001_PRESSURE_PA,
        MAY_ROW_001_FEED_X_METHANE,
        BINARY_FINGERPRINT,
    )
    multiplier = initialization["outer"]["multiplier"]
    result = _equilibrium._held_stage_ii_lower_search(
        binary_capsule,
        MAY_ROW_001_TEMPERATURE_K,
        MAY_ROW_001_PRESSURE_PA,
        MAY_ROW_001_FEED_X_METHANE,
        multiplier,
        -1.0e6,
        [(cut["x_methane"], cut["log_volume"]) for cut in initialization["cuts"]],
        BINARY_FINGERPRINT,
    )

    assert result["status"] == "search_exhausted"
    assert result["starts_completed"] == 20
    assert len(result["planned_starts"]) == 20
    assert [start["class"] for start in result["planned_starts"]] == [
        ("shifted_previous", "component_near", "stratified")[index % 3] for index in range(20)
    ]
    assert result["accepted_cuts"]
    for cut in result["accepted_cuts"]:
        assert cut["objective"] == pytest.approx(
            cut["g_bar"] + multiplier * (MAY_ROW_001_FEED_X_METHANE - cut["x_methane"]),
            rel=2.0e-15,
            abs=2.0e-15,
        )
    objectives = [cut["objective"] for cut in result["accepted_cuts"]]
    assert max(objectives) > min(objectives) + 1.0e-6


def test_held_stage_ii_candidate_lifecycle_deduplicates_by_lowest_objective() -> None:
    result = _equilibrium._held_stage_ii_candidates(
        0.5,
        0.005,
        0.0,
        [
            ("duplicate_higher", 0.005, 0.2, 1.0e-3, 0.0),
            ("duplicate_lower", 0.001, 0.2005, 1.0005e-3, 0.0),
            ("second", 0.0, 0.8, 1.0e-2, 0.0),
        ],
    )

    assert result["status"] == "stage_iii_ready"
    assert result["candidate_ids"] == ["duplicate_lower", "second"]
    assert result["rejections"] == [
        {
            "identity": "duplicate_higher",
            "reason": "duplicate_replaced_by_lower_objective",
        }
    ]


def test_held_stage_ii_candidate_lifecycle_rejects_a_third_distinct_candidate() -> None:
    result = _equilibrium._held_stage_ii_candidates(
        0.5,
        0.005,
        0.0,
        [
            ("first", 0.0, 0.2, 1.0e-3, 0.0),
            ("second", 0.0, 0.5, 1.0e-4, 0.0),
            ("third", 0.0, 0.8, 1.0e-2, 0.0),
        ],
    )

    assert result["status"] == "scope_exceeded"
    assert result["candidate_ids"] == ["first", "second", "third"]


@pytest.mark.parametrize(
    ("major_iterations", "lower_starts", "lower_satisfied"),
    [
        pytest.param(100, 0, True, id="major-iteration-cap"),
        pytest.param(0, 20, False, id="lower-start-cap"),
    ],
)
def test_held_stage_ii_resource_caps_fail_closed(
    major_iterations: int,
    lower_starts: int,
    lower_satisfied: bool,
) -> None:
    assert (
        _equilibrium._held_stage_ii_budget_status(
            major_iterations,
            lower_starts,
            lower_satisfied,
        )
        == "search_exhausted"
    )


def test_held_stage_ii_returns_two_distinct_candidates_with_monotone_bounds(
    unstable_stage_ii: dict[str, object],
) -> None:
    assert unstable_stage_ii["outcome"] == "stage_iii_ready"
    assert unstable_stage_ii["search_status"] == "two_candidates"
    assert unstable_stage_ii["search_profile"] == STAGE_II_PROFILE
    assert unstable_stage_ii["globality_certificate"] == "not_guaranteed"
    assert 1 <= unstable_stage_ii["major_iterations"] <= 100
    candidates = unstable_stage_ii["candidates"]
    assert len(candidates) == 2
    relative_density_difference = abs(
        1.0 / candidates[0]["volume_m3"] - 1.0 / candidates[1]["volume_m3"]
    ) / max(1.0 / candidates[0]["volume_m3"], 1.0 / candidates[1]["volume_m3"])
    assert (
        abs(candidates[0]["x_methane"] - candidates[1]["x_methane"]) >= 1.0e-3
        or relative_density_difference >= 1.0e-3
    )

    trace = unstable_stage_ii["trace"]
    assert len(trace) == unstable_stage_ii["major_iterations"]
    upper_bounds = [entry["upper_bound"] for entry in trace]
    assert upper_bounds == pytest.approx(
        [min(upper_bounds[: index + 1]) for index in range(len(upper_bounds))],
        rel=0.0,
        abs=2.0e-12,
    )
    assert all(entry["active_cut_ids"] for entry in trace)
    assert all(1 <= entry["lower_starts_completed"] <= 20 for entry in trace)
    retained_ids = {cut["identity"] for cut in unstable_stage_ii["cuts"]}
    assert all(cut_id in retained_ids for entry in trace for cut_id in entry["accepted_cut_ids"])
    for cut in unstable_stage_ii["cuts"]:
        assert cut["intercept"] == pytest.approx(cut["g_bar"])
        assert cut["slope"] == pytest.approx(MAY_ROW_001_FEED_X_METHANE - cut["x_methane"])
    assert any(
        rejection["reason"] == "above_upper" for entry in trace for rejection in entry["rejections"]
    )


def test_held_stage_ii_does_not_start_after_stage_i_no_negative(
    binary_capsule: object,
) -> None:
    result = _equilibrium._held_stage_ii(
        binary_capsule,
        MAY_ROW_011_TEMPERATURE_K,
        MAY_ROW_011_PRESSURE_PA,
        MAY_ROW_011_LIQUID_SIDE_FEED_X_METHANE,
        BINARY_FINGERPRINT,
    )

    assert result["outcome"] == "not_required"
    assert result["search_status"] == "stage_i_no_negative"
    assert result["stage_i_outcome"] == "no_negative_found"
    assert result["major_iterations"] == 0
    assert result["endpoint_attempts"] == []
    assert result["cuts"] == []
    assert result["trace"] == []


def test_held_stage_ii_lower_callback_failures_remain_diagnostic(
    binary_capsule: object,
) -> None:
    result = _equilibrium._held_stage_ii_lower_search(
        binary_capsule,
        MAY_ROW_001_TEMPERATURE_K,
        MAY_ROW_001_PRESSURE_PA,
        MAY_ROW_001_FEED_X_METHANE,
        0.0,
        0.0,
        [],
        "sha256:wrong-provider-identity",
    )

    assert result["status"] == "search_exhausted"
    assert result["starts_completed"] == 20
    assert result["accepted_cuts"] == []
    assert all(not attempt["accepted"] for attempt in result["attempts"])
    assert all(
        attempt["solver_status"] or attempt["callback_error"] for attempt in result["attempts"]
    )
    assert all(attempt["callback_error"] for attempt in result["attempts"])


@pytest.fixture(scope="module")
def binary_capsule() -> object:
    return epcsaft.native_sdk(_binary_model())


@pytest.fixture(scope="module")
def unstable_stage_i(binary_capsule: object) -> dict[str, object]:
    return _equilibrium._held_stage_i(
        binary_capsule,
        MAY_ROW_001_TEMPERATURE_K,
        MAY_ROW_001_PRESSURE_PA,
        MAY_ROW_001_FEED_X_METHANE,
        BINARY_FINGERPRINT,
    )


@pytest.fixture(scope="module")
def unstable_stage_ii(binary_capsule: object) -> dict[str, object]:
    return _equilibrium._held_stage_ii(
        binary_capsule,
        MAY_ROW_001_TEMPERATURE_K,
        MAY_ROW_001_PRESSURE_PA,
        MAY_ROW_001_FEED_X_METHANE,
        BINARY_FINGERPRINT,
    )


@pytest.fixture(scope="module")
def one_phase_stage_i(binary_capsule: object) -> dict[str, object]:
    return _equilibrium._held_stage_i(
        binary_capsule,
        MAY_ROW_011_TEMPERATURE_K,
        MAY_ROW_011_PRESSURE_PA,
        MAY_ROW_011_LIQUID_SIDE_FEED_X_METHANE,
        BINARY_FINGERPRINT,
    )


def test_held_stage_iii_lever_rule_candidate_neighborhoods_and_infeasibility() -> None:
    assert MAY_ROW_001_ID == "may2015-ch4-c2h6-001"
    initialized = _equilibrium._held_stage_iii_initialize(
        MAY_ROW_001_FEED_X_METHANE,
        MAY_ROW_001_REFINED_CANDIDATES,
    )

    liquid_fraction = (MAY_ROW_001_REFINED_CANDIDATES[1][1] - MAY_ROW_001_FEED_X_METHANE) / (
        MAY_ROW_001_REFINED_CANDIDATES[1][1] - MAY_ROW_001_REFINED_CANDIDATES[0][1]
    )
    expected_fractions = (liquid_fraction, 1.0 - liquid_fraction)
    expected_amounts = (
        expected_fractions[0] * MAY_ROW_001_REFINED_CANDIDATES[0][1],
        expected_fractions[0] * (1.0 - MAY_ROW_001_REFINED_CANDIDATES[0][1]),
        expected_fractions[1] * MAY_ROW_001_REFINED_CANDIDATES[1][1],
        expected_fractions[1] * (1.0 - MAY_ROW_001_REFINED_CANDIDATES[1][1]),
    )
    assert initialized["status"] == "initialized"
    assert initialized["phase_fractions"] == pytest.approx(expected_fractions)
    assert tuple(
        initialized["initial_variables"][index] for index in (0, 1, 3, 4)
    ) == pytest.approx(expected_amounts)
    assert initialized["initial_variables"][2] == pytest.approx(
        math.log(expected_fractions[0] * MAY_ROW_001_REFINED_CANDIDATES[0][2])
    )
    assert initialized["initial_variables"][5] == pytest.approx(
        math.log(expected_fractions[1] * MAY_ROW_001_REFINED_CANDIDATES[1][2])
    )
    for index, center in zip((0, 1, 3, 4), expected_amounts, strict=True):
        assert initialized["lower_bounds"][index] == pytest.approx(center - 1.0e-3)
        assert initialized["upper_bounds"][index] == pytest.approx(center + 1.0e-3)
    assert tuple(initialized["lower_bounds"][index] for index in (2, 5)) == pytest.approx(
        (math.log(1.0e-5), math.log(1.0e-5))
    )
    assert tuple(initialized["upper_bounds"][index] for index in (2, 5)) == pytest.approx(
        (math.log(1.0e-1), math.log(1.0e-1))
    )

    infeasible = _equilibrium._held_stage_iii_initialize(
        MAY_ROW_001_FEED_X_METHANE,
        (("same_side_1", 0.7, 1.0e-4), ("same_side_2", 0.8, 2.0e-4)),
    )
    assert infeasible["status"] == "return_to_stage_ii"
    assert infeasible["failure_reason"] == "candidate lever rule does not contain the feed"
    assert infeasible["initial_variables"] == []

    trace_bound = _equilibrium._held_stage_iii_initialize(
        MAY_ROW_001_FEED_X_METHANE,
        (("trace", 1.0e-8, 1.0e-4), ("interior", 0.8, 2.0e-4)),
    )
    assert trace_bound["status"] == "scope_exceeded"
    assert trace_bound["failure_reason"] == "trace-bound candidate requires a later slice"


@pytest.mark.parametrize(
    ("override", "status", "reason"),
    (
        (
            {"composition_distance": 5.0e-4, "phase_density_distance": 5.0e-4},
            "return_to_stage_ii",
            "phase_collapse",
        ),
        (
            {"material_balance_max_abs": 2.0e-10},
            "return_to_stage_ii",
            "material_balance_failure",
        ),
        (
            {"pressure_stationarity_max_relative": 2.0e-8},
            "return_to_stage_ii",
            "pressure_failure",
        ),
        ({"kkt_stationarity_max_abs": 2.0e-8}, "return_to_stage_ii", "kkt_failure"),
        ({"inactive_bounds": False}, "return_to_stage_ii", "active_bound_failure"),
        ({"held_gap": 2.0e-6}, "return_to_stage_ii", "held_gap_failure"),
        (
            {"chemical_potential_max_relative": 2.0e-6},
            "return_to_stage_ii",
            "chemical_potential_failure",
        ),
        (
            {"confirmation_succeeded": False},
            "return_to_stage_ii",
            "confirmation_failure",
        ),
        ({}, "accepted", ""),
    ),
)
def test_held_stage_iii_acceptance_is_strict_and_fail_closed(
    override: dict[str, object],
    status: str,
    reason: str,
) -> None:
    evidence = {
        "solver_converged": True,
        "callback_error": "",
        "solver_constraint_violation": 1.0e-12,
        "material_balance_max_abs": 1.0e-12,
        "pressure_stationarity_max_relative": 1.0e-10,
        "kkt_stationarity_max_abs": 1.0e-10,
        "inactive_bounds": True,
        "composition_distance": 0.3,
        "phase_density_distance": 0.5,
        "held_gap": 5.0e-7,
        "chemical_potential_max_relative": 5.0e-7,
        "confirmation_succeeded": True,
        "confirmation_max_difference": 1.0e-9,
    }
    evidence.update(override)

    result = _equilibrium._held_stage_iii_classify(evidence)

    assert result == {"status": status, "failure_reason": reason}


def test_held_stage_iii_zero_safe_chemical_potential_metric() -> None:
    assert _equilibrium._held_stage_iii_mu_difference((0.0, 0.0), (0.0, 0.0)) == 0.0
    assert _equilibrium._held_stage_iii_mu_difference(
        (1.0e-9, -2.0), (-1.0e-9, -2.0)
    ) == pytest.approx(2.0e-9)


def test_held_stage_iii_exact_derivatives_at_may_row_001(
    binary_capsule: object,
) -> None:
    assert MAY_ROW_001_ID == "may2015-ch4-c2h6-001"
    initialized = _equilibrium._held_stage_iii_initialize(
        MAY_ROW_001_FEED_X_METHANE,
        MAY_ROW_001_REFINED_CANDIDATES,
    )
    variables = tuple(initialized["initial_variables"])
    direction = (0.03, -0.02, 0.04, -0.01, 0.02, -0.03)
    step = 2.0e-5

    def evaluate(scale: float) -> dict[str, object]:
        point = tuple(
            value + scale * delta for value, delta in zip(variables, direction, strict=True)
        )
        return _equilibrium._held_evaluate_stage_iii(
            binary_capsule,
            MAY_ROW_001_TEMPERATURE_K,
            MAY_ROW_001_PRESSURE_PA,
            MAY_ROW_001_FEED_X_METHANE,
            point,
            BINARY_FINGERPRINT,
        )

    lower = evaluate(-step)
    center = evaluate(0.0)
    upper = evaluate(step)
    assert center["constraints"] == pytest.approx((0.0, 0.0), abs=2.0e-15)
    assert center["jacobian"] == [
        1.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        1.0,
        0.0,
    ]
    objective_directional = (upper["objective"] - lower["objective"]) / (2.0 * step)
    assert objective_directional == pytest.approx(
        sum(value * delta for value, delta in zip(center["gradient"], direction, strict=True)),
        rel=2.0e-7,
        abs=2.0e-8,
    )
    hessian_lower = center["hessian_lower"]
    hessian = [[0.0] * 6 for _ in range(6)]
    index = 0
    for row in range(6):
        for column in range(row + 1):
            hessian[row][column] = hessian_lower[index]
            hessian[column][row] = hessian_lower[index]
            index += 1
    for row in range(6):
        gradient_directional = (upper["gradient"][row] - lower["gradient"][row]) / (2.0 * step)
        hessian_directional = sum(hessian[row][column] * direction[column] for column in range(6))
        assert gradient_directional == pytest.approx(
            hessian_directional,
            rel=2.0e-6,
            abs=2.0e-7,
        )


def test_held_stage_iii_confirms_retained_may_row_001_refinement(
    binary_capsule: object,
) -> None:
    result = _equilibrium._held_stage_iii(
        binary_capsule,
        MAY_ROW_001_TEMPERATURE_K,
        MAY_ROW_001_PRESSURE_PA,
        MAY_ROW_001_FEED_X_METHANE,
        MAY_ROW_001_REFINED_G_BAR,
        MAY_ROW_001_REFINED_CANDIDATES,
        BINARY_FINGERPRINT,
    )

    assert result["outcome"] == "accepted"
    assert result["search_profile"] == STAGE_III_PROFILE
    assert result["confirmation_solves"] == 1
    assert result["confirmation_succeeded"] is True
    assert result["confirmation_max_difference"] <= 1.0e-7
    assert [attempt["role"] for attempt in result["attempt_log"]] == [
        "refinement",
        "confirmation",
    ]
    assert (
        max(
            abs(
                result["attempt_log"][1]["initial_guess"][index]
                - result["attempt_log"][0]["initial_guess"][index]
            )
            for index in (0, 1, 3, 4)
        )
        >= 1.0e-4
    )
    assert result["material_balance_max_abs"] <= 1.0e-10
    assert result["pressure_stationarity_max_relative"] <= 1.0e-8
    assert result["kkt_stationarity_max_abs"] <= 1.0e-8
    assert result["inactive_bounds"] is True
    assert result["chemical_potential_max_relative"] <= 1.0e-6
    assert 0.0 <= result["held_gap"] <= 1.0e-6
    assert result["total_g_bar"] == pytest.approx(MAY_ROW_001_REFINED_G_BAR, abs=2.0e-10)


def test_held_controller_returns_failed_refinement_to_stage_ii_and_fails_closed(
    binary_capsule: object,
) -> None:
    result = _equilibrium._held(
        binary_capsule,
        MAY_ROW_001_TEMPERATURE_K,
        MAY_ROW_001_PRESSURE_PA,
        MAY_ROW_001_FEED_X_METHANE,
        BINARY_FINGERPRINT,
    )

    assert result["outcome"] == "scope_exceeded"
    assert result["search_status"] == "third_candidate"
    assert result["stage_i_outcome"] == "negative_tpd"
    assert result["stage_iii_attempts"]
    assert result["stage_iii_attempts"][0]["outcome"] == "return_to_stage_ii"
    assert result["major_iterations"] > result["stage_iii_attempts"][0]["major_iteration"]
    trace = result["trace"]
    assert all(entry["active_cut_ids"] for entry in trace)
    assert all(1 <= entry["lower_starts_completed"] <= 20 for entry in trace)
    assert any(entry["accepted_cut_ids"] for entry in trace)
    assert all(
        set(rejection) == {"identity", "reason"}
        for entry in trace
        for rejection in entry["rejections"]
    )
    assert any(
        rejection["reason"] == "above_upper" for entry in trace for rejection in entry["rejections"]
    )
    feedback = [entry for entry in trace if entry["stage_iii_outcome"] == "return_to_stage_ii"]
    assert len(feedback) == 1
    assert feedback[0]["stage_iii_failure_reason"]
    assert result["globality_certificate"] == "not_guaranteed"


def test_held_controller_preserves_stage_i_one_phase_outcome(
    binary_capsule: object,
) -> None:
    result = _equilibrium._held(
        binary_capsule,
        MAY_ROW_011_TEMPERATURE_K,
        MAY_ROW_011_PRESSURE_PA,
        MAY_ROW_011_LIQUID_SIDE_FEED_X_METHANE,
        BINARY_FINGERPRINT,
    )

    assert result["outcome"] == "one_phase"
    assert result["search_status"] == "stage_i_no_negative"
    assert result["stage_i_outcome"] == "no_negative_found"
    assert result["major_iterations"] == 0
    assert result["stage_iii_attempts"] == []


def test_held_one_mole_transform_and_exact_derivatives_at_may_row_001(
    binary_capsule: object,
) -> None:
    assert MAY_ROW_001_ID == "may2015-ch4-c2h6-001"
    x_methane = MAY_ROW_001_FEED_X_METHANE
    log_volume = math.log(1.0 / 15_000.0)
    volume = math.exp(log_volume)
    direction = (0.03, -0.04)
    step = 2.0e-5

    def evaluate(scale: float) -> dict[str, object]:
        return _equilibrium._held_evaluate_state(
            binary_capsule,
            MAY_ROW_001_TEMPERATURE_K,
            MAY_ROW_001_PRESSURE_PA,
            x_methane + scale * direction[0],
            log_volume + scale * direction[1],
            BINARY_FINGERPRINT,
        )

    raw = _equilibrium.evaluate_mixture_phase(
        binary_capsule,
        MAY_ROW_001_TEMPERATURE_K,
        (x_methane, 1.0 - x_methane),
        volume,
        BINARY_FINGERPRINT,
    )
    lower = evaluate(-step)
    center = evaluate(0.0)
    upper = evaluate(step)

    pressure_over_rt = MAY_ROW_001_PRESSURE_PA / (
        GAS_CONSTANT_J_PER_MOL_K * MAY_ROW_001_TEMPERATURE_K
    )
    expected_gradient = (
        raw["gradient"][0] - raw["gradient"][1],
        volume * (raw["gradient"][2] + pressure_over_rt),
    )
    expected_hessian = (
        raw["hessian"][0] - raw["hessian"][1] - raw["hessian"][3] + raw["hessian"][4],
        volume * (raw["hessian"][2] - raw["hessian"][5]),
        volume * (raw["hessian"][6] - raw["hessian"][7]),
        volume * volume * raw["hessian"][8] + volume * (raw["gradient"][2] + pressure_over_rt),
    )
    assert center["amounts_mol"] == pytest.approx((x_methane, 1.0 - x_methane))
    assert center["volume_m3"] == pytest.approx(volume)
    assert center["g_bar"] == pytest.approx(
        raw["helmholtz_over_rt_reference_amount"] + pressure_over_rt * volume
    )
    assert center["gradient"] == pytest.approx(expected_gradient)
    assert center["hessian"] == pytest.approx(expected_hessian)
    assert center["pressure_pa"] == pytest.approx(raw["pressure_pa"])
    assert center["pressure_stationarity_relative"] == pytest.approx(
        (raw["pressure_pa"] - MAY_ROW_001_PRESSURE_PA) / MAY_ROW_001_PRESSURE_PA
    )

    objective_directional = (upper["g_bar"] - lower["g_bar"]) / (2.0 * step)
    assert objective_directional == pytest.approx(
        sum(value * delta for value, delta in zip(center["gradient"], direction, strict=True)),
        rel=2.0e-7,
        abs=2.0e-8,
    )
    for row in range(2):
        gradient_directional = (upper["gradient"][row] - lower["gradient"][row]) / (2.0 * step)
        hessian_directional = sum(
            center["hessian"][2 * row + column] * direction[column] for column in range(2)
        )
        assert gradient_directional == pytest.approx(
            hessian_directional,
            rel=2.0e-6,
            abs=2.0e-7,
        )


def test_held_tpd_and_tunneling_exact_derivatives_at_may_row_001(
    binary_capsule: object,
    unstable_stage_i: dict[str, object],
) -> None:
    assert MAY_ROW_001_ID == "may2015-ch4-c2h6-001"
    reference = unstable_stage_i["reference"]
    minimum = unstable_stage_i["best_state"]
    trial_x_methane = MAY_ROW_001_VAPOR_X_METHANE
    trial_log_volume = math.log(1.0 / 15_000.0)
    direction = (0.03, -0.04)
    step = 3.0e-5

    def evaluate_tpd(scale: float) -> dict[str, object]:
        return _equilibrium._held_evaluate_tpd(
            binary_capsule,
            MAY_ROW_001_TEMPERATURE_K,
            MAY_ROW_001_PRESSURE_PA,
            MAY_ROW_001_FEED_X_METHANE,
            reference["log_volume"],
            trial_x_methane + scale * direction[0],
            trial_log_volume + scale * direction[1],
            BINARY_FINGERPRINT,
        )

    def evaluate_tunneling(scale: float) -> dict[str, object]:
        return _equilibrium._held_evaluate_tunneling(
            binary_capsule,
            MAY_ROW_001_TEMPERATURE_K,
            MAY_ROW_001_PRESSURE_PA,
            MAY_ROW_001_FEED_X_METHANE,
            reference["log_volume"],
            minimum["x_methane"],
            minimum["log_volume"],
            trial_x_methane + scale * direction[0],
            trial_log_volume + scale * direction[1],
            BINARY_FINGERPRINT,
        )

    def check_centered_derivatives(
        evaluate: Callable[[float], dict[str, object]],
        objective_name: str,
    ) -> dict[str, object]:
        lower = evaluate(-step)
        center = evaluate(0.0)
        upper = evaluate(step)

        # Here the directional gradient is O(1e-2) and the Hessian-vector
        # components are O(1e-1). This step balances centered truncation and
        # roundoff; 2e-8 relative plus 2e-10 absolute remains tight against
        # the observed O(1e-11) directional discrepancies.
        objective_directional = (upper[objective_name] - lower[objective_name]) / (2.0 * step)
        exact_objective_directional = sum(
            value * delta for value, delta in zip(center["gradient"], direction, strict=True)
        )
        assert objective_directional == pytest.approx(
            exact_objective_directional,
            rel=2.0e-8,
            abs=2.0e-10,
        )
        for row in range(2):
            gradient_directional = (upper["gradient"][row] - lower["gradient"][row]) / (2.0 * step)
            exact_hessian_directional = sum(
                center["hessian"][2 * row + column] * direction[column] for column in range(2)
            )
            assert gradient_directional == pytest.approx(
                exact_hessian_directional,
                rel=2.0e-8,
                abs=2.0e-10,
            )
        return center

    center_tpd = check_centered_derivatives(evaluate_tpd, "d_bar")
    center_tunneling = check_centered_derivatives(evaluate_tunneling, "objective")
    assert center_tunneling["hessian"][1] == pytest.approx(
        center_tunneling["hessian"][2],
        rel=0.0,
        abs=1.0e-12,
    )
    expected_tunneling = (center_tpd["d_bar"] - unstable_stage_i["best_tpd"]) * math.exp(
        1.0e-3 / abs(trial_x_methane - minimum["x_methane"])
    )
    assert center_tunneling["objective"] == pytest.approx(
        expected_tunneling,
        rel=2.0e-15,
        abs=2.0e-15,
    )


def test_held_reference_selection_and_tpd_sign_at_may_row_001(
    binary_capsule: object,
    unstable_stage_i: dict[str, object],
) -> None:
    reference = unstable_stage_i["reference"]
    reference_attempts = unstable_stage_i["reference_attempts"]
    accepted_references = [attempt for attempt in reference_attempts if attempt["accepted"]]
    assert len(reference_attempts) == 2
    assert accepted_references
    assert reference["g_bar"] == pytest.approx(
        min(attempt["g_bar"] for attempt in accepted_references)
    )
    assert abs(reference["pressure_stationarity_relative"]) <= 1.0e-8

    at_reference = _equilibrium._held_evaluate_tpd(
        binary_capsule,
        MAY_ROW_001_TEMPERATURE_K,
        MAY_ROW_001_PRESSURE_PA,
        MAY_ROW_001_FEED_X_METHANE,
        reference["log_volume"],
        MAY_ROW_001_FEED_X_METHANE,
        reference["log_volume"],
        BINARY_FINGERPRINT,
    )
    at_best = _equilibrium._held_evaluate_tpd(
        binary_capsule,
        MAY_ROW_001_TEMPERATURE_K,
        MAY_ROW_001_PRESSURE_PA,
        MAY_ROW_001_FEED_X_METHANE,
        reference["log_volume"],
        unstable_stage_i["best_state"]["x_methane"],
        unstable_stage_i["best_state"]["log_volume"],
        BINARY_FINGERPRINT,
    )
    assert at_reference["d_bar"] == pytest.approx(0.0, abs=1.0e-12)
    assert at_reference["gradient"] == pytest.approx((0.0, 0.0), abs=1.0e-8)
    assert at_best["d_bar"] == pytest.approx(unstable_stage_i["best_tpd"], abs=1.0e-10)
    assert at_best["d_bar"] < TPD_NEGATIVE_THRESHOLD

    finite_tunnel = _equilibrium._held_evaluate_tunneling(
        binary_capsule,
        MAY_ROW_001_TEMPERATURE_K,
        MAY_ROW_001_PRESSURE_PA,
        MAY_ROW_001_FEED_X_METHANE,
        reference["log_volume"],
        unstable_stage_i["best_state"]["x_methane"],
        unstable_stage_i["best_state"]["log_volume"],
        MAY_ROW_001_FEED_X_METHANE,
        reference["log_volume"],
        BINARY_FINGERPRINT,
    )
    assert math.isfinite(finite_tunnel["objective"])
    numerical_pole_distance = 1.0e-3 / math.log(sys.float_info.max)
    near_pole = _equilibrium._held_evaluate_tunneling(
        binary_capsule,
        MAY_ROW_001_TEMPERATURE_K,
        MAY_ROW_001_PRESSURE_PA,
        MAY_ROW_001_FEED_X_METHANE,
        reference["log_volume"],
        unstable_stage_i["best_state"]["x_methane"],
        unstable_stage_i["best_state"]["log_volume"],
        unstable_stage_i["best_state"]["x_methane"] + 2.0 * numerical_pole_distance,
        unstable_stage_i["best_state"]["log_volume"],
        BINARY_FINGERPRINT,
    )
    assert math.isfinite(near_pole["objective"])
    assert all(math.isfinite(value) for value in near_pole["gradient"])
    assert all(math.isfinite(value) for value in near_pole["hessian"])
    with pytest.raises(ValueError, match="singular composition neighborhood"):
        _equilibrium._held_evaluate_tunneling(
            binary_capsule,
            MAY_ROW_001_TEMPERATURE_K,
            MAY_ROW_001_PRESSURE_PA,
            MAY_ROW_001_FEED_X_METHANE,
            reference["log_volume"],
            unstable_stage_i["best_state"]["x_methane"],
            unstable_stage_i["best_state"]["log_volume"],
            unstable_stage_i["best_state"]["x_methane"] + 0.5 * numerical_pole_distance,
            unstable_stage_i["best_state"]["log_volume"],
            BINARY_FINGERPRINT,
        )


def test_held_stage_i_uses_the_frozen_deterministic_20_start_profile(
    one_phase_stage_i: dict[str, object],
) -> None:
    assert MAY_ROW_011_ID == "may2015-ch4-c2h6-011"
    starts = one_phase_stage_i["planned_starts"]
    assert one_phase_stage_i["search_profile"] == STAGE_I_PROFILE
    assert len(starts) == 20
    assert [start["role"] for start in starts] == [
        "feed_near_lower_liquid",
        "feed_near_lower_vapor",
        "feed_near_upper_liquid",
        "feed_near_upper_vapor",
        "component_1_rich_liquid",
        "component_1_rich_vapor",
        "component_2_rich_liquid",
        "component_2_rich_vapor",
        "stratified_low_liquid",
        "stratified_low_vapor",
        "stratified_high_liquid",
        "stratified_high_vapor",
        "stratified_quarter_liquid",
        "stratified_quarter_vapor",
        "stratified_three_quarter_liquid",
        "stratified_three_quarter_vapor",
        "stratified_inner_low_liquid",
        "stratified_inner_low_vapor",
        "stratified_inner_high_liquid",
        "stratified_inner_high_vapor",
    ]
    starts_by_role = {start["role"]: start for start in starts}
    assert starts_by_role["component_1_rich_liquid"]["x_methane"] == pytest.approx(1.0 - 1.0e-4)
    assert starts_by_role["component_1_rich_vapor"]["x_methane"] == pytest.approx(1.0 - 1.0e-4)
    assert starts_by_role["component_2_rich_liquid"]["x_methane"] == pytest.approx(1.0e-4)
    assert starts_by_role["component_2_rich_vapor"]["x_methane"] == pytest.approx(1.0e-4)
    assert all(1.0e-8 <= start["x_methane"] <= 1.0 - 1.0e-8 for start in starts)
    assert all(1.0e-5 <= math.exp(start["log_volume"]) <= 1.0e-1 for start in starts)


def test_held_stage_i_confirms_negative_tpd_from_a_material_perturbation(
    unstable_stage_i: dict[str, object],
) -> None:
    assert unstable_stage_i["outcome"] == "negative_tpd"
    assert unstable_stage_i["best_tpd"] < TPD_NEGATIVE_THRESHOLD
    assert unstable_stage_i["negative_confirmations"] == 1
    assert unstable_stage_i["confirmation_max_difference"] <= 1.0e-7
    assert unstable_stage_i["globality_certificate"] == "not_guaranteed"
    confirmation_attempts = [
        attempt for attempt in unstable_stage_i["attempt_log"] if attempt["kind"] == "confirmation"
    ]
    assert len(confirmation_attempts) == 1
    assert confirmation_attempts[0]["materially_perturbed"] is True
    assert confirmation_attempts[0]["accepted"] is True
    assert unstable_stage_i["negative_confirmations"] == sum(
        attempt["accepted"] for attempt in confirmation_attempts
    )


def test_held_stage_i_distinguishes_no_negative_found_from_global_stability(
    one_phase_stage_i: dict[str, object],
) -> None:
    assert one_phase_stage_i["outcome"] == "no_negative_found"
    assert one_phase_stage_i["search_status"] == "source_heuristic_complete"
    assert one_phase_stage_i["starts_completed"] == 20
    assert one_phase_stage_i["best_tpd"] >= TPD_NEGATIVE_THRESHOLD
    assert one_phase_stage_i["globality_certificate"] == "not_guaranteed"
    assert one_phase_stage_i["failure_reason"] == ""


def test_held_stage_i_reports_reference_failure_as_indeterminate(
    binary_capsule: object,
) -> None:
    result = _equilibrium._held_stage_i(
        binary_capsule,
        MAY_ROW_001_TEMPERATURE_K,
        MAY_ROW_001_PRESSURE_PA,
        MAY_ROW_001_FEED_X_METHANE,
        "sha256:wrong-provider-identity",
    )

    assert result["outcome"] == "indeterminate"
    assert result["search_status"] == "reference_failed"
    assert result["starts_completed"] == 0
    assert result["globality_certificate"] == "not_guaranteed"
    assert result["failure_reason"]
    assert len(result["reference_attempts"]) == 2
    assert all(not attempt["accepted"] for attempt in result["reference_attempts"])


@pytest.mark.parametrize(
    ("temperature_k", "pressure_pa", "feed_x_methane"),
    [
        pytest.param(math.nan, MAY_ROW_001_PRESSURE_PA, 0.5, id="nan-temperature"),
        pytest.param(MAY_ROW_001_TEMPERATURE_K, math.nan, 0.5, id="nan-pressure"),
        pytest.param(
            MAY_ROW_001_TEMPERATURE_K,
            MAY_ROW_001_PRESSURE_PA,
            math.nan,
            id="nan-feed",
        ),
        pytest.param(0.0, MAY_ROW_001_PRESSURE_PA, 0.5, id="zero-temperature"),
        pytest.param(MAY_ROW_001_TEMPERATURE_K, -1.0, 0.5, id="negative-pressure"),
        pytest.param(
            MAY_ROW_001_TEMPERATURE_K,
            MAY_ROW_001_PRESSURE_PA,
            1.0,
            id="noninterior-feed",
        ),
    ],
)
def test_held_stage_i_rejects_invalid_input_before_planning_starts(
    binary_capsule: object,
    temperature_k: float,
    pressure_pa: float,
    feed_x_methane: float,
) -> None:
    result = _equilibrium._held_stage_i(
        binary_capsule,
        temperature_k,
        pressure_pa,
        feed_x_methane,
        BINARY_FINGERPRINT,
    )

    def assert_finite_diagnostic_numbers(value: object) -> None:
        if isinstance(value, float):
            assert math.isfinite(value)
        elif isinstance(value, dict):
            for item in value.values():
                assert_finite_diagnostic_numbers(item)
        elif isinstance(value, list):
            for item in value:
                assert_finite_diagnostic_numbers(item)

    assert result["outcome"] == "indeterminate"
    assert result["search_status"] == "invalid_input"
    assert result["planned_starts"] == []
    assert result["reference_attempts"] == []
    assert result["attempt_log"] == []
    assert result["starts_completed"] == 0
    assert result["failure_reason"]
    assert_finite_diagnostic_numbers(result)
