from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
BASELINE_PATH = REPOSITORY_ROOT / "data/reference/perdomo_held2_baseline_v1.json"
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
    assert baseline["schema"] == "perdomo-held2-manufactured-baseline-v1"
    assert baseline["source"]["commit"] == "de28d47dbafa82cfbeabb636fb511ebe149f8c18"
    assert baseline["source"]["tree"] == "8e2e18da45a778eb789647fca38517294c769388"
    assert baseline["artifacts"]["provider_wheel"]["sha256"] == (
        "9e4da0d7ba7896bcd2ec096400553d935e0516c61f1bd9f41f2370ab68ab36ea"
    )
    assert baseline["artifacts"]["equilibrium_wheel"]["sha256"] == (
        "f96c779661abf8ce708d3d013973721620211309b269d1967a6f5f2e88cf4cfd"
    )
    assert baseline["runtime_evidence"]["stage_i"]["declared_start_count"] == 30
    assert baseline["runtime_evidence"]["stage_i"]["completed_start_count"] == 30
    assert baseline["runtime_evidence"]["stage_i"]["failed_start_count"] == 0
    assert baseline["runtime_evidence"]["stage_ii"]["lower_starts_per_iteration"] == 30
    assert baseline["runtime_evidence"]["stage_ii"]["outcome"] == "candidate_set"
    assert baseline["runtime_evidence"]["stage_iii"]["physical_status"] == "accepted"
    assert baseline["runtime_evidence"]["globality_certificate"] == "not_guaranteed"
    assert baseline["start_terminal_exposure"]["status"] == "aggregate_only"
    assert baseline["historical_dual_pullback_partition"]["status"] == "evidence_gap"
    assert baseline["historical_dual_pullback_partition"]["reported_counts"] == {
        "broad_failures": 12,
        "narrow_failures": 11,
        "recovered": 25,
    }
