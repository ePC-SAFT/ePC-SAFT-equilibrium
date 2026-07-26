from __future__ import annotations

import epcsaft
import pytest

import epcsaft_equilibrium
from epcsaft_equilibrium import _equilibrium

WATER_MOLAR_MASS_G_PER_MOL = 18.0153
PERDOMO_TABLE3_NACL_MOL_PER_KG_WATER = 5.6
_source_amounts = (
    1000.0 / WATER_MOLAR_MASS_G_PER_MOL,
    PERDOMO_TABLE3_NACL_MOL_PER_KG_WATER,
    PERDOMO_TABLE3_NACL_MOL_PER_KG_WATER,
)
PERDOMO_TABLE3_FEED = tuple(value / sum(_source_amounts) for value in _source_amounts)


def _perdomo_table3_model() -> epcsaft.EPCSAFT:
    parameters = epcsaft.ParameterBundle.from_catalog(
        "figiel-2025-reference-electrolytes", version=1
    ).select(("water", "sodium-cation", "chloride-anion"))
    return epcsaft.EPCSAFT(parameters)


def _without_times(value: object) -> object:
    if isinstance(value, dict):
        return {
            key: 0.0 if key in {"wall_seconds", "cpu_seconds"} else _without_times(item)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_without_times(item) for item in value]
    return value


def test_held2_observer_is_quiet_by_default_and_does_not_change_results(capfd) -> None:
    model = _perdomo_table3_model()
    arguments = (
        epcsaft.native_sdk(model),
        298.15,
        2508.0,
        PERDOMO_TABLE3_FEED,
        model.parameter_fingerprint,
    )

    quiet = _equilibrium._solve_tp_flash(*arguments)
    quiet_output = capfd.readouterr()
    traced = _equilibrium._solve_tp_flash(*arguments, trace=True)
    traced_output = capfd.readouterr()

    assert quiet_output.out == ""
    assert quiet_output.err == ""
    assert _without_times(traced) == _without_times(quiet)
    assert "HELD2.0  case=installed-held2-paper-rewrite" in traced_output.out
    assert "STEP 1" in traced_output.out
    assert "STEP 2" in traced_output.out
    assert "REFERENCE PRESSURE ROOTS" in traced_output.out
    assert "reason=declared_search_complete" in traced_output.out


def test_perdomo_table3_nacl_workflow(held2_live: bool) -> None:
    model = _perdomo_table3_model()
    result = _equilibrium._solve_tp_flash(
        epcsaft.native_sdk(model),
        298.15,
        2508.0,
        PERDOMO_TABLE3_FEED,
        model.parameter_fingerprint,
        trace=held2_live,
    )

    expected_feed = tuple(value / sum(_source_amounts) for value in _source_amounts)
    assert PERDOMO_TABLE3_FEED == pytest.approx(expected_feed, abs=1.0e-15)
    assert sum(PERDOMO_TABLE3_FEED) == pytest.approx(1.0, abs=1.0e-15)
    assert PERDOMO_TABLE3_FEED[1] - PERDOMO_TABLE3_FEED[2] == pytest.approx(0.0, abs=1.0e-15)
    assert result["controller"] == "perdomo_held2_paper_steps_1_to_10_v1"
    assert result["outcome"] == "one_phase_no_negative_witness_detected"
    assert result["failure_reason"] is None
    assert result["failure_stage"] is None
    assert len(result["phases"]) == 1
    assert result["phases"][0]["phase_fraction"] == pytest.approx(1.0)
    reference = result["step2"]["reference_pressure_envelope"]
    assert reference["outcome"] == "selected"
    assert [root["mechanical_class"] for root in reference["roots"]] == [
        "strict_stable",
        "unstable",
        "strict_stable",
    ]
    assert reference["selected_root_index"] == 2
    selected = reference["roots"][reference["selected_root_index"]]
    assert selected["pressure_residual"] == pytest.approx(0.0, abs=1.0e-8)
    assert selected["objective"] == min(
        root["objective"]
        for root in reference["roots"]
        if root["mechanical_class"] == "strict_stable"
    )

    assert result["step2"]["reason"] == "declared_search_complete"
    assert result["step2"]["minimum_tpd"] >= -1.0e-8
    assert result["upper_solve_count"] == 0
    assert result["step4_history"] == []
    assert result["step8_history"] == []
    assert result["step10"] is None
    assert result["globality_certificate"] == "not_guaranteed"


def test_public_tp_flash_uses_the_paper_held2_route(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        _equilibrium,
        "evaluate_electrolyte_phase",
        lambda *_args, **_kwargs: pytest.fail(
            "Python must not recompute native phase thermodynamics"
        ),
    )
    result = epcsaft_equilibrium.tp_flash(
        _perdomo_table3_model(),
        298.15 * epcsaft.unit_registry.kelvin,
        2508.0 * epcsaft.unit_registry.pascal,
        PERDOMO_TABLE3_FEED,
    )

    assert result.diagnostics.outcome == "one_phase_no_negative_witness_detected"
    assert result.diagnostics.globality_certificate == "not_guaranteed"
    assert result.overall_mole_fractions == pytest.approx(PERDOMO_TABLE3_FEED)
    assert result.phase_fractions == pytest.approx((1.0,))
    assert result.phases[0].mole_fractions == pytest.approx(PERDOMO_TABLE3_FEED)
