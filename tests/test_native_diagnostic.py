from __future__ import annotations

import json
import math
import os
import subprocess
from pathlib import Path

import epcsaft
import pytest

from epcsaft_equilibrium import _equilibrium
from tests.test_held import (
    MAY_ROW_001_FEED_X_METHANE,
    MAY_ROW_001_PRESSURE_PA,
    MAY_ROW_001_TEMPERATURE_K,
    MAY_ROW_011_LIQUID_SIDE_FEED_X_METHANE,
    MAY_ROW_011_PRESSURE_PA,
    MAY_ROW_011_TEMPERATURE_K,
    _binary_model,
)
from tests.test_perdomo_held2_trace import (
    PERDOMO_TABLE3_FEED,
    _perdomo_table3_model,
)


def _diagnostic_executable() -> Path:
    configured = os.environ.get("EPCSAFT_EQUILIBRIUM_DIAGNOSTIC")
    executable = (
        Path(configured)
        if configured is not None
        else Path(__file__).parents[1] / "build" / "epcsaft-equilibrium-diagnostic"
    )
    if not executable.is_file():
        pytest.skip(
            "configure with EPCSAFT_EQUILIBRIUM_BUILD_DIAGNOSTIC=ON "
            "or set EPCSAFT_EQUILIBRIUM_DIAGNOSTIC"
        )
    return executable


def _normalize_json(value: object) -> object:
    if isinstance(value, dict):
        return {
            key: 0.0 if key in {"wall_seconds", "cpu_seconds"} else _normalize_json(item)
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [_normalize_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


@pytest.mark.parametrize("route", ("neutral", "neutral_stage_ii", "held2"))
def test_native_diagnostic_matches_python_payload(
    tmp_path: Path,
    route: str,
) -> None:
    if route == "neutral":
        model = _binary_model()
        temperature_k = MAY_ROW_011_TEMPERATURE_K
        pressure_pa = MAY_ROW_011_PRESSURE_PA
        feed = (
            MAY_ROW_011_LIQUID_SIDE_FEED_X_METHANE,
            1.0 - MAY_ROW_011_LIQUID_SIDE_FEED_X_METHANE,
        )
    elif route == "neutral_stage_ii":
        model = _binary_model()
        temperature_k = MAY_ROW_001_TEMPERATURE_K
        pressure_pa = MAY_ROW_001_PRESSURE_PA
        feed = (
            MAY_ROW_001_FEED_X_METHANE,
            1.0 - MAY_ROW_001_FEED_X_METHANE,
        )
    else:
        model = _perdomo_table3_model()
        temperature_k = 298.15
        pressure_pa = 2508.0
        feed = PERDOMO_TABLE3_FEED

    model_config = tmp_path / f"{route}-model.json"
    epcsaft.export_native_model(model, model_config)
    command = (
        str(_diagnostic_executable()),
        "--model-config",
        str(model_config),
        "--temperature",
        repr(temperature_k),
        "--pressure",
        repr(pressure_pa),
        "--feed",
        ",".join(repr(value) for value in feed),
    )
    first = subprocess.run(command, check=True, capture_output=True, text=True)
    second = subprocess.run(command, check=True, capture_output=True, text=True)
    python_payload = _equilibrium._solve_tp_flash(
        epcsaft.native_sdk(model),
        temperature_k,
        pressure_pa,
        feed,
        model.parameter_fingerprint,
    )

    first_payload = json.loads(first.stdout)
    assert _normalize_json(first_payload) == _normalize_json(python_payload)
    assert _normalize_json(json.loads(second.stdout)) == _normalize_json(json.loads(first.stdout))
    assert first.stderr == second.stderr == ""
    if route == "held2":
        assert all(
            isinstance(step8["phase_coalescences"], list)
            for step8 in first_payload["step8_history"]
        )


def test_native_diagnostic_streams_trace_and_writes_json(tmp_path: Path) -> None:
    model = _perdomo_table3_model()
    model_config = tmp_path / "held2-model.json"
    output = tmp_path / "held2-result.json"
    epcsaft.export_native_model(model, model_config)

    completed = subprocess.run(
        (
            str(_diagnostic_executable()),
            "--model-config",
            str(model_config),
            "--temperature",
            "298.15",
            "--pressure",
            "2508.0",
            "--feed",
            ",".join(repr(value) for value in PERDOMO_TABLE3_FEED),
            "--trace",
            "--output",
            str(output),
        ),
        check=True,
        capture_output=True,
        text=True,
    )

    assert completed.stdout == ""
    assert "HELD2.0  case=installed-held2-paper-rewrite" in completed.stderr
    assert "STEP 1" in completed.stderr
    assert "STEP 2" in completed.stderr
    assert "REFERENCE PRESSURE ROOTS" in completed.stderr
    assert "reason=declared_search_complete" in completed.stderr
    assert json.loads(output.read_text())["globality_certificate"] == "not_guaranteed"
