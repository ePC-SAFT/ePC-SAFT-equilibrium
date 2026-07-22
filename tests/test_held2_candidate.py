from __future__ import annotations

import json
from pathlib import Path

import pytest

CANDIDATE = Path(__file__).parents[1] / "data" / "reference" / "perdomo_held2_candidate_v1.json"


def test_private_held2_candidate_is_exact_and_fail_closed() -> None:
    evidence = json.loads(CANDIDATE.read_text(encoding="utf-8"))

    assert evidence["schema"] == "perdomo-held2-integrated-candidate-v1"
    assert evidence["source"] == {
        "commit": "40d01d32dbcda67caa5eafc71b04abe74b80f016",
        "tree": "2dfec318226194fb556a33c50ce01b501915122b",
    }
    assert evidence["artifacts"]["provider_wheel"]["sha256"] == (
        "9e4da0d7ba7896bcd2ec096400553d935e0516c61f1bd9f41f2370ab68ab36ea"
    )
    assert evidence["artifacts"]["equilibrium_wheel"]["sha256"] == (
        "606686b766162a4366b2274113a954329ba1496a9b4a27f819b5e2c030979188"
    )

    controller = evidence["controller"]
    assert controller["outcome"] == "indeterminate_stage_ii"
    assert controller["failure_stage"] == "stage_ii"
    assert controller["globality_certificate"] == "not_guaranteed"
    assert controller["predictive_comparison_status"] == ("not_allowed_before_physical_acceptance")
    assert controller["stage_iii"] is None

    reference = controller["reference_pressure_envelope"]
    assert reference["outcome"] == "selected"
    assert [root["mechanical_class"] for root in reference["roots"]] == [
        "strict_stable",
        "unstable",
        "strict_stable",
    ]
    selected = reference["roots"][reference["selected_root_index"]]
    assert selected["volume"] == pytest.approx(3.909560419950766e-5)
    assert selected["objective"] == pytest.approx(-6.009742067562037)

    assert controller["stage_i"]["outcome"] == "negative_witness_found"
    assert controller["stage_i"]["minimum_tpd"] < 0.0
    stage_ii = controller["stage_ii"]
    assert stage_ii["outcome"] == "indeterminate_lower_search"
    assert stage_ii["local_attempts_truncated"] is True
    assert stage_ii["attempt_classification"] == {
        "declared": 1,
        "physical_kkt_passed": 0,
        "solver_converged": 0,
        "solver_failed": 1,
        "step6_eligible": 0,
    }
    assert stage_ii["attempt_trace"][0]["provider_status"] == "provider_exact"
    assert stage_ii["attempt_trace"][0]["callback_error"]
