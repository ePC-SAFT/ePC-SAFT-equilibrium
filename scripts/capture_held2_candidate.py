#!/usr/bin/env python3
"""Capture or verify the private installed Perdomo HELD2 candidate evidence."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import subprocess
import zipfile
from pathlib import Path
from typing import Any

import epcsaft

from epcsaft_equilibrium import _equilibrium

SCHEMA = "perdomo-held2-integrated-candidate-v1"
FEED = (
    0.6462224836985831,
    0.04995832720498113,
    0.2648918958670393,
    0.01946364661469824,
    0.01946364661469824,
)


def _sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _installed_distribution(name: str) -> dict[str, Any]:
    distribution = importlib.metadata.distribution(name)
    direct_url = json.loads(distribution.read_text("direct_url.json") or "{}")
    if direct_url.get("dir_info", {}).get("editable", False):
        raise ValueError(f"{name} must be installed non-editably")
    root = Path(distribution.locate_file("")).resolve()
    if "site-packages" not in root.parts:
        raise ValueError(f"{name} is not imported from site-packages")
    return {"version": distribution.version, "location": "site-packages"}


def _wheel_record(wheel: Path, installed_module: Path) -> dict[str, Any]:
    wheel = wheel.resolve(strict=True)
    installed_module = installed_module.resolve(strict=True)
    with zipfile.ZipFile(wheel) as archive:
        matches = [name for name in archive.namelist() if name.endswith(installed_module.name)]
        if len(matches) != 1:
            raise ValueError(f"{wheel.name} does not uniquely contain {installed_module.name}")
        payload = archive.read(matches[0])
    payload_sha256 = hashlib.sha256(payload).hexdigest()
    installed_sha256 = _sha256_path(installed_module)
    if payload_sha256 != installed_sha256:
        raise ValueError(f"installed module {installed_module} does not match {wheel.name}")
    return {
        "filename": wheel.name,
        "member": matches[0],
        "native_payload_sha256": payload_sha256,
        "sha256": _sha256_path(wheel),
        "size_bytes": wheel.stat().st_size,
    }


def _git_value(source_root: Path, revision: str) -> str:
    return subprocess.run(
        ["git", "rev-parse", revision],
        cwd=source_root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _capture(args: argparse.Namespace) -> dict[str, Any]:
    source_root = args.source_root.resolve(strict=True)
    source_commit = _git_value(source_root, "HEAD")
    source_tree = _git_value(source_root, "HEAD^{tree}")
    if source_commit != args.source_commit or source_tree != args.source_tree:
        raise ValueError("declared source commit/tree does not match source root")

    provider_native = next(
        path for path in Path(epcsaft.__file__).resolve().parent.glob("*.so") if path.is_file()
    )
    equilibrium_native = Path(_equilibrium.__file__).resolve()
    model = epcsaft.EPCSAFT(
        epcsaft.ParameterBundle.from_catalog(
            "khudaida-2026-figure-2-electrolyte-lle", version=1
        ).select(
            (
                "water",
                "ethanol",
                "isobutanol",
                "sodium-cation",
                "chloride-anion",
            )
        )
    )
    controller = _equilibrium._held2_controller(
        epcsaft.native_sdk(model),
        293.15,
        100_000.0,
        FEED,
        model.parameter_fingerprint,
        50,
        1,
        1,
    )
    if controller["globality_certificate"] != "not_guaranteed":
        raise ValueError("HELD2 candidate overclaimed globality")
    if controller["predictive_comparison_status"] != ("not_allowed_before_physical_acceptance"):
        raise ValueError("predictive comparison ran before physical acceptance")

    return {
        "schema": SCHEMA,
        "source": {"commit": source_commit, "tree": source_tree},
        "artifacts": {
            "provider_distribution": _installed_distribution("epcsaft"),
            "equilibrium_distribution": _installed_distribution("epcsaft-equilibrium"),
            "provider_wheel": _wheel_record(args.provider_wheel, provider_native),
            "equilibrium_wheel": _wheel_record(args.equilibrium_wheel, equilibrium_native),
        },
        "case": {
            "temperature_k": 293.15,
            "pressure_pa": 100_000.0,
            "physical_feed": FEED,
            "parameter_fingerprint": model.parameter_fingerprint,
        },
        "controller": controller,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--provider-wheel", type=Path, required=True)
    parser.add_argument("--equilibrium-wheel", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--source-tree", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()

    captured = _capture(args)
    serialized = json.dumps(captured, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.verify:
        if args.output.read_text(encoding="utf-8") != serialized:
            raise SystemExit("installed HELD2 candidate evidence differs")
    else:
        args.output.write_text(serialized, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
