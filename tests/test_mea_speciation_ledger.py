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
        "blockers": 6,
        "executable": False,
        "reaction_rank": 5,
        "species_count": 9,
        "status": "blocked",
    }

    executable_result = run_validator("--require-executable")
    assert executable_result.returncode == 2
    assert "ledger is not executable" in executable_result.stderr

    document = json.loads(LEDGER.read_text(encoding="utf-8"))
    provider = document["provider"]
    assert provider["provider_commit"] == "5c4cd54b3596e51331ca9f6c871daec34a72eb4f"
    assert provider["provider_wheel_sha256"] == (
        "3bd36c7bdf946a11eadf407dbfc0262ced73fbd80e8431478a104b835fc55f4a"
    )
    assert provider["application_commit"] == (
        "efdb5bb7372733b3ea16f87823fc6bdbf4d652a8"
    )
    assert provider["bundle_fingerprint"] == (
        "sha256:7b38f85711562a251b2e0e261762b130eaff4366a57e82bb47896787b73222ea"
    )
    assert provider["selected_parameter_fingerprint"] == (
        "sha256:a89481803f4d9f3009feed3c77bc06e815dbe512b592c838eae14505694cc78f"
    )
    assert [record["provider_component_id"] for record in document["species"]] == (
        provider["component_ids"]
    )
    assert [record["charge"] for record in document["species"]] == provider["charges"]
    assert document["feed"]["molar_masses_kg_per_mol"] == {
        "CO2": 0.04401,
        "H2O": 0.01801528,
        "MEA": 0.06108,
    }
    blocker_ids = {record["identity"] for record in document["blockers"]}
    assert "provider-mea-bundle-missing" not in blocker_ids
    assert "feed-molar-mass-authority-missing" not in blocker_ids
    assert "provider-mea-bundle-not-source-complete" in blocker_ids


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
