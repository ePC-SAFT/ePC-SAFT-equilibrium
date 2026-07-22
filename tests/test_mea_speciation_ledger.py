from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
LEDGER = REPO_ROOT / "data" / "reference" / "mea_speciation_input_ledger.json"
VALIDATOR = REPO_ROOT / "scripts" / "validate_mea_speciation_ledger.py"


def run_validator(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(VALIDATOR), *args],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )


def test_ledger_is_stoichiometrically_complete_but_not_executable() -> None:
    result = run_validator()

    assert result.returncode == 0, result.stderr
    report = json.loads(result.stdout)
    assert report == {
        "balance_rank": 4,
        "blockers": 7,
        "executable": False,
        "reaction_rank": 5,
        "species_count": 9,
        "status": "blocked",
    }

    executable_result = run_validator("--require-executable")
    assert executable_result.returncode == 2
    assert "ledger is not executable" in executable_result.stderr


@pytest.mark.parametrize(
    ("section", "field"),
    [
        ("species", "identity"),
        ("equilibrium_constant_candidates", "unit"),
        ("equilibrium_constant_candidates", "basis"),
        ("equilibrium_constant_candidates", "source"),
        ("observation_sets", "role"),
    ],
)
def test_validator_rejects_incomplete_records(
    tmp_path: Path,
    section: str,
    field: str,
) -> None:
    document = json.loads(LEDGER.read_text(encoding="utf-8"))
    del document[section][0][field]
    incomplete_ledger = tmp_path / "incomplete.json"
    incomplete_ledger.write_text(json.dumps(document), encoding="utf-8")

    result = run_validator("--ledger", str(incomplete_ledger))

    assert result.returncode == 1
    assert f"{section}[0].{field}" in result.stderr
