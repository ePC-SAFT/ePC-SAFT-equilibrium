from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
BASELINE_PATH = REPOSITORY_ROOT / "data/reference/perdomo_held2_baseline_v2.json"
CAPTURE_SCRIPT = REPOSITORY_ROOT / "scripts/capture_held2_baseline.py"


def test_held2_baseline_fixture_replays_the_claimed_manufactured_evidence() -> None:
    completed = subprocess.run(
        [sys.executable, str(CAPTURE_SCRIPT), "--verify", str(BASELINE_PATH)],
        cwd=REPOSITORY_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr

    baseline_text = BASELINE_PATH.read_text(encoding="utf-8")
    assert "/home/" not in baseline_text
    baseline = json.loads(baseline_text)
    assert baseline["schema"] == "perdomo-held2-manufactured-baseline-v2"
    assert baseline["source"]["commit"] == "60b9b60e3150e4b295a9c95976449264a8b83907"
    assert baseline["source"]["tree"] == "7dd06a7b35c6394782d31f9b258fca10625396c0"
    assert baseline["artifacts"]["provider_wheel"]["sha256"] == (
        "9e4da0d7ba7896bcd2ec096400553d935e0516c61f1bd9f41f2370ab68ab36ea"
    )
    assert baseline["artifacts"]["equilibrium_wheel"]["sha256"] == (
        "9910205c9f4c90760c3e6c751e98275e93637230dc370170bf032dd679dee04e"
    )
    assert baseline["runtime_evidence"]["stage_i"]["declared_start_count"] == 30
    assert baseline["runtime_evidence"]["stage_i"]["completed_start_count"] == 30
    assert baseline["runtime_evidence"]["stage_i"]["failed_start_count"] == 0
    assert baseline["runtime_evidence"]["stage_ii"]["lower_starts_per_iteration"] == 30
    assert baseline["runtime_evidence"]["stage_ii"]["outcome"] == "candidate_set"
    assert baseline["runtime_evidence"]["stage_iii"]["physical_status"] == "accepted"
    assert baseline["runtime_evidence"]["stage_iii"]["stage_iii_solve_count"] == 2
    assert [
        step["action"]
        for step in baseline["runtime_evidence"]["stage_iii"]["lifecycle"]
        if step["action"] != "retain_phase"
    ] == ["merge_duplicate", "accept_active_set"]
    assert baseline["runtime_evidence"]["globality_certificate"] == "not_guaranteed"
    assert baseline["start_terminal_exposure"]["status"] == "aggregate_only"
    assert baseline["historical_dual_pullback_partition"]["status"] == "evidence_gap"
    assert baseline["historical_dual_pullback_partition"]["reported_counts"] == {
        "broad_failures": 12,
        "narrow_failures": 11,
        "recovered": 25,
    }
