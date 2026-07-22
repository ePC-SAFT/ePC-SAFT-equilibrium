#!/usr/bin/env python3
"""Capture or verify the private manufactured Perdomo HELD2 baseline."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import platform
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Any

import epcsaft

from epcsaft_equilibrium import _equilibrium

SCHEMA = "perdomo-held2-manufactured-baseline-v1"
CHARGES = (0.0, 1.0, -1.0)
PHYSICAL_FEED = (0.5, 0.25, 0.25)
CHEMICAL_POTENTIALS = (3.0, -2.0, 4.0)


def _canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _command_line(*command: str) -> str:
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    return completed.stdout.strip().splitlines()[0]


def _distribution_origin(name: str) -> dict[str, object]:
    distribution = importlib.metadata.distribution(name)
    direct_url_text = distribution.read_text("direct_url.json") or "{}"
    direct_url = json.loads(direct_url_text)
    editable = bool(direct_url.get("dir_info", {}).get("editable", False))
    location = Path(distribution.locate_file("")).resolve()
    if editable or "site-packages" not in location.parts:
        raise ValueError(f"{name} must be installed non-editably in site-packages")
    return {"editable": False, "kind": "site-packages"}


def _wheel_payload(wheel: Path, member_suffix: str) -> tuple[str, str]:
    with zipfile.ZipFile(wheel) as archive:
        matches = sorted(name for name in archive.namelist() if name.endswith(member_suffix))
        if len(matches) != 1:
            raise ValueError(
                f"{wheel.name} must contain exactly one member ending in {member_suffix!r}"
            )
        member = matches[0]
        return member, _sha256_bytes(archive.read(member))


def _figiel_brine_model() -> epcsaft.EPCSAFT:
    parameters = epcsaft.ParameterBundle.from_catalog(
        "figiel-2025-reference-electrolytes", version=1
    ).select(("water", "sodium-cation", "chloride-anion"))
    return epcsaft.EPCSAFT(parameters)


def _provider_phase_evidence() -> tuple[dict[str, Any], str, dict[str, str]]:
    model = _figiel_brine_model()
    capsule = epcsaft.native_sdk(model)
    result = _equilibrium._held2_adapter(
        capsule,
        298.15,
        100_000.0,
        (0.02,),
        math.log(1.0e-3),
        model.parameter_fingerprint,
    )
    try:
        _equilibrium._held2_adapter(
            capsule,
            1.0,
            100_000.0,
            (0.02,),
            math.log(1.0e-3),
            model.parameter_fingerprint,
        )
    except Exception as error:
        failure = {"type": type(error).__name__, "message": str(error)}
    else:
        raise RuntimeError("the Provider source-domain failure was not propagated")
    return result, model.parameter_fingerprint, failure


def _runtime_evidence() -> dict[str, Any]:
    provider_phase, parameter_fingerprint, provider_failure = _provider_phase_evidence()
    provider_phase = {
        key: value
        for key, value in provider_phase.items()
        if key != "pressure_stationarity_derivative_log_volume"
    }
    direct = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        (0.27, 0.73, 0.92, 1.08),
        CHEMICAL_POTENTIALS,
    )
    stage_i = _equilibrium._held2_adapter(CHARGES, PHYSICAL_FEED, "stage_i")
    stage_ii = _equilibrium._held2_adapter(
        CHARGES, PHYSICAL_FEED, "stage_ii_legacy"
    )
    # Schema v1 was captured before per-attempt Stage-II observability existed.
    # Keep its immutable artifact payload comparable while newer tests freeze
    # the complete deterministic attempt trace independently.
    stage_ii = {
        key: value
        for key, value in stage_ii.items()
        if key
        not in {
            "attempt_classification",
            "attempt_trace",
            "historical_dual_pullback_fixture_status",
        }
    }
    stage_ii["bound_history"] = [
        {
            key: entry[key]
            for key in ("cut_count", "lower_bound", "multiplier", "upper_bound")
        }
        for entry in stage_ii["bound_history"]
    ]
    candidates = ((0.2, 1.0), (0.20000002, 1.0), (0.8, 1.0))
    stage_iii = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        candidates,
        "stage_iii_legacy_baseline",
    )
    stage_iii_failure = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        ((0.1, 1.0), (0.2, 1.0), (0.3, 1.0)),
        "stage_iii_legacy_baseline",
    )
    # Schema v1 predates the active-set lifecycle. Replay its exact private
    # numerical route and omit only fields that did not exist in that schema.
    new_stage_iii_fields = {
        "active_set_resolve_count",
        "bound_complementarity_inf_norm",
        "dual_sign_violation_inf_norm",
        "lifecycle",
        "minimum_phase_fraction",
        "retired_inactive_count",
        "stage_iii_solve_count",
    }
    stage_iii = {
        key: value for key, value in stage_iii.items() if key not in new_stage_iii_fields
    }
    stage_iii_failure = {
        key: value
        for key, value in stage_iii_failure.items()
        if key not in new_stage_iii_fields
    }
    stage_iii_derivatives = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        candidates,
        (0.25, 0.2, 1.0, 0.25, 0.21, 1.0, 0.5, 0.795, 1.0),
        (0.3, -0.2),
        "stage_iii_derivatives",
    )
    return {
        "direct_physical_lift": direct,
        "globality_certificate": "not_guaranteed",
        "parameter_fingerprint": parameter_fingerprint,
        "provider_failure": provider_failure,
        "provider_phase": provider_phase,
        "stage_i": stage_i,
        "stage_ii": stage_ii,
        "stage_iii": stage_iii,
        "stage_iii_derivatives": stage_iii_derivatives,
        "stage_iii_failure": stage_iii_failure,
    }
def _artifact_record(wheel: Path, member_suffix: str, installed_path: Path) -> dict[str, Any]:
    member, member_sha256 = _wheel_payload(wheel, member_suffix)
    installed_sha256 = _sha256_path(installed_path)
    if member_sha256 != installed_sha256:
        raise ValueError(
            f"installed payload {installed_path} does not match {member} from {wheel.name}"
        )
    return {
        "filename": wheel.name,
        "installed_payload_origin": f"site-packages/{member}",
        "member": member,
        "payload_sha256": member_sha256,
        "sha256": _sha256_path(wheel),
        "size_bytes": wheel.stat().st_size,
    }


def _capture(args: argparse.Namespace) -> dict[str, Any]:
    provider_wheel = args.provider_wheel.resolve(strict=True)
    equilibrium_wheel = args.equilibrium_wheel.resolve(strict=True)
    source_root = args.source_root.resolve(strict=True)
    source_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=source_root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    source_tree = subprocess.run(
        ["git", "rev-parse", "HEAD^{tree}"],
        cwd=source_root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if source_commit != args.source_commit or source_tree != args.source_tree:
        raise ValueError("the declared source commit/tree does not match the source root")

    native_matches = sorted(Path(epcsaft.__file__).resolve().parent.glob("*.so"))
    if len(native_matches) != 1:
        raise ValueError("could not resolve the unique installed Provider native module")
    provider_native = native_matches[0]
    equilibrium_native = Path(_equilibrium.__file__).resolve()
    runtime = _runtime_evidence()

    return {
        "schema": SCHEMA,
        "source": {
            "commit": source_commit,
            "tree": source_tree,
            "artifact_build_worktree_clean": args.artifact_build_source_clean,
        },
        "artifacts": {
            "equilibrium_wheel": _artifact_record(
                equilibrium_wheel,
                equilibrium_native.name,
                equilibrium_native,
            ),
            "provider_wheel": _artifact_record(
                provider_wheel,
                provider_native.name,
                provider_native,
            ),
        },
        "environment": {
            "cmake": _command_line("cmake", "--version"),
            "compiler": _command_line("c++", "--version"),
            "equilibrium_distribution_origin": _distribution_origin("epcsaft-equilibrium"),
            "ipopt": _command_line("pkg-config", "--modversion", "ipopt"),
            "ipopt_link": _command_line("pkg-config", "--libs", "ipopt"),
            "platform": platform.platform(),
            "provider_distribution_origin": _distribution_origin("epcsaft"),
            "python": sys.version,
        },
        "capture_command": (
            "python scripts/capture_held2_baseline.py --capture <output> "
            "--source-root <clean-source> --source-commit <commit> --source-tree <tree> "
            "--artifact-build-source-clean --provider-wheel <provider.whl> "
            "--equilibrium-wheel <equilibrium.whl>"
        ),
        "runtime_evidence": runtime,
        "runtime_evidence_sha256": _sha256_bytes(_canonical_bytes(runtime)),
        "start_terminal_exposure": {
            "status": "aggregate_only",
            "available": {
                "stage_i": [
                    "declared_start_count",
                    "completed_start_count",
                    "failed_start_count",
                    "candidates",
                ],
                "stage_ii": [
                    "major_iterations",
                    "lower_starts_per_iteration",
                    "bound_history",
                    "candidates",
                ],
            },
            "unavailable": [
                "per-start initial coordinates",
                "per-start terminal coordinates",
                "per-start solver status",
                "per-start physical KKT certificate",
            ],
            "reason": (
                "The current private manufactured adapter exposes aggregate start accounting, "
                "candidate states, and Stage-II bound history, but not individual start records."
            ),
        },
        "historical_dual_pullback_partition": {
            "status": "evidence_gap",
            "reported_counts": {
                "broad_failures": 12,
                "narrow_failures": 11,
                "recovered": 25,
            },
            "source": "docs/plans/2026-07-21-perdomo-held2-solver-strategy.md#task-0",
            "reason": (
                "No authoritative machine-readable terminal fixture for the historical "
                "25/11/12 partition was assigned to this repository leaf."
            ),
        },
    }


def _verify(path: Path) -> None:
    baseline = json.loads(path.read_text(encoding="utf-8"))
    if baseline.get("schema") != SCHEMA:
        raise ValueError(f"unsupported baseline schema: {baseline.get('schema')!r}")
    runtime = _runtime_evidence()
    actual = _sha256_bytes(_canonical_bytes(runtime))
    expected = baseline.get("runtime_evidence_sha256")
    if actual != expected:
        raise ValueError(f"runtime evidence drift: expected {expected}, got {actual}")
    if runtime != baseline.get("runtime_evidence"):
        raise ValueError("runtime evidence payload differs despite its recorded digest")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--capture", type=Path, metavar="OUTPUT")
    action.add_argument("--verify", type=Path, metavar="FIXTURE")
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--source-commit")
    parser.add_argument("--source-tree")
    parser.add_argument("--artifact-build-source-clean", action="store_true")
    parser.add_argument("--provider-wheel", type=Path)
    parser.add_argument("--equilibrium-wheel", type=Path)
    return parser


def main() -> int:
    parser = _parser()
    args = parser.parse_args()
    if args.verify is not None:
        _verify(args.verify)
        print(f"verified {args.verify}")
        return 0
    required = (
        "source_root",
        "source_commit",
        "source_tree",
        "artifact_build_source_clean",
        "provider_wheel",
        "equilibrium_wheel",
    )
    missing = [name for name in required if getattr(args, name) is None]
    if missing:
        required_flags = ", ".join(f"--{name.replace('_', '-')}" for name in missing)
        parser.error("--capture requires " + required_flags)
    payload = _capture(args)
    args.capture.parent.mkdir(parents=True, exist_ok=True)
    args.capture.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"captured {args.capture}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
