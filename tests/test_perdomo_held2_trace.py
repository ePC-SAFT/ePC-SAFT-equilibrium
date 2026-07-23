from __future__ import annotations

import epcsaft
import pytest

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


def test_held2_terminal_progress_formats_expected_rows() -> None:
    trace = _equilibrium._held2_terminal_progress_fixture()

    assert "HELD2.0  case=perdomo-table3-nacl-water" in trace
    assert "REFERENCE PRESSURE ROOTS" in trace
    assert "STAGE I - DIRECT-L TPD SEARCH" in trace
    assert "eval" in trace
    assert "TPD" in trace
    assert "P rel." in trace
    assert "lower=unavailable" in trace
    assert "STAGE II IPOPT" in trace
    assert "stage=STAGE II STEP 6" in trace
    assert "stage=STAGE III PRESSURE POLISH" in trace
    assert "stage=STAGE III FREE ENERGY" in trace
    assert "objective=-2.000000e+00" in trace
    assert "upper=-2.000000e+00" in trace
    assert "gap= 0.000000e+00" in trace
    assert "FAILURE" in trace
    assert "STAGE II - SKIPPED" in trace
    assert "FINAL" in trace
    assert trace.endswith("\n")


def test_held2_terminal_progress_disables_itself_on_output_failure() -> None:
    assert _equilibrium._held2_terminal_progress_failure_fixture() is True


def test_held2_observer_is_quiet_by_default_and_does_not_change_results(capfd) -> None:
    model = _perdomo_table3_model()
    arguments = (
        epcsaft.native_sdk(model),
        298.15,
        2508.0,
        PERDOMO_TABLE3_FEED,
        model.parameter_fingerprint,
        1,
        1,
        1,
    )

    quiet = _equilibrium._held2_controller(*arguments)
    quiet_output = capfd.readouterr()
    traced = _equilibrium._held2_controller(*arguments, trace=True)
    traced_output = capfd.readouterr()

    assert quiet_output.out == ""
    assert quiet_output.err == ""
    assert traced == quiet
    assert "HELD2.0  case=installed-held2-controller" in traced_output.out
    assert "REFERENCE PRESSURE ROOTS" in traced_output.out
    assert "STAGE I - DIRECT-L TPD SEARCH" in traced_output.out
    assert "x_ion=" in traced_output.out
    assert "limit= 3.800000e-01" in traced_output.out
    assert "STAGE II - SKIPPED" in traced_output.out


def test_perdomo_stage_i_first_trial_respects_provider_ionic_domain() -> None:
    model = _perdomo_table3_model()
    result = _equilibrium._held2_stage_i_direct(
        epcsaft.native_sdk(model),
        298.15,
        2508.0,
        PERDOMO_TABLE3_FEED,
        model.parameter_fingerprint,
        1,
    )

    assert result["total_ion_mole_fraction_max"] == pytest.approx(0.38)
    assert result["completed_evaluation_count"] >= 1
    assert result["failed_evaluation_count"] == 0
    first = result["evaluations"][0]
    assert first["physical_total_ion_mole_fraction"] <= 0.38 + 1.0e-12
    assert "status 3" not in first["failure_reason"]


def test_perdomo_stage_ii_sobol_trials_respect_provider_ionic_domain() -> None:
    model = _perdomo_table3_model()
    coordinate_evidence = _equilibrium._held2_adapter(
        [0.0, 1.0, -1.0],
        PERDOMO_TABLE3_FEED,
        [0.0, 0.0, 0.0],
    )
    feed_independent = [coordinate_evidence["modified_feed"][1]]

    result = _equilibrium._held2_stage_ii_basin_explorer(
        epcsaft.native_sdk(model),
        298.15,
        2508.0,
        [feed_independent],
        model.parameter_fingerprint,
        1,
    )

    assert result["outcome"] == "representatives_found"
    assert result["failed_evaluation_count"] == 0
    assert result["completed_evaluation_count"] == 2
    assert all(
        evaluation["independent_modified_fractions"][0] <= 0.38 + 1.0e-12
        and "status 3" not in evaluation["failure_reason"]
        for evaluation in result["evaluations"]
    )


def test_perdomo_table3_nacl_workflow(held2_live: bool) -> None:
    model = _perdomo_table3_model()
    result = _equilibrium._held2_controller(
        epcsaft.native_sdk(model),
        298.15,
        2508.0,
        PERDOMO_TABLE3_FEED,
        model.parameter_fingerprint,
        50,
        8,
        50,
        trace=held2_live,
    )

    expected_feed = tuple(value / sum(_source_amounts) for value in _source_amounts)
    assert PERDOMO_TABLE3_FEED == pytest.approx(expected_feed, abs=1.0e-15)
    assert sum(PERDOMO_TABLE3_FEED) == pytest.approx(1.0, abs=1.0e-15)
    assert PERDOMO_TABLE3_FEED[1] - PERDOMO_TABLE3_FEED[2] == pytest.approx(
        0.0, abs=1.0e-15
    )
    reference = result["reference_pressure_envelope"]
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

    stage_i = result["stage_i"]
    assert stage_i["outcome"] == "no_negative_witness_detected"
    assert stage_i["termination_reason"] == "declared_budget_exhausted"
    assert stage_i["completed_evaluation_count"] == 50
    assert stage_i["failed_evaluation_count"] == 0
    assert stage_i["minimum_tpd"] >= -1.0e-8
    assert stage_i["total_ion_mole_fraction_max"] == pytest.approx(0.38)
    assert all(
        evaluation["physical_total_ion_mole_fraction"] <= 0.38 + 1.0e-12
        and "status 3" not in evaluation["failure_reason"]
        for evaluation in stage_i["evaluations"]
    )
    assert result["outcome"] == "stage_i_finite_search_without_negative_witness"
    assert result["failure_reason"] == "declared_budget_exhausted"
    assert result["failure_stage"] == "stage_i"
    assert result["stage_ii_skip_reason"] == "stage_i_negative_witness_not_found"
    assert result["stage_iii_skip_reason"] == "stage_i_negative_witness_not_found"
    assert result["stage_ii"] is None
    assert result["stage_iii"] is None
    assert result["globality_certificate"] == "not_guaranteed"
