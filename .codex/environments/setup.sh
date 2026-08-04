#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

for required_command in git python3 uv cmake ninja pkg-config c++ sha256sum; do
    if ! command -v "$required_command" >/dev/null 2>&1; then
        echo "missing required command: $required_command" >&2
        exit 1
    fi
done

if ! pkg-config --exists ipopt; then
    echo "missing required pkg-config package: ipopt" >&2
    exit 1
fi

git_common_dir="$(git rev-parse --path-format=absolute --git-common-dir)"
canonical_repo_root="$(dirname "$git_common_dir")"
project_root="$(dirname "$canonical_repo_root")"
artifact_tool="${project_root}/ePC-SAFT-governance/tools/artifact_store.py"
if [[ ! -f "$artifact_tool" ]]; then
    echo "missing Governance artifact resolver: $artifact_tool" >&2
    exit 1
fi
eos_wheel="$(python3 "$artifact_tool" resolve --distribution epcsaft)"
if [[ ! -f "$eos_wheel" ]]; then
    echo "missing active hash-addressed EOS wheel: $eos_wheel" >&2
    exit 1
fi

eos_wheel="$(realpath "$eos_wheel")"
if [[ "$eos_wheel" == */dist/* ]]; then
    echo "mutable EOS dist/ wheels are forbidden: $eos_wheel" >&2
    exit 1
fi

eos_wheel_sha256="$(sha256sum "$eos_wheel" | awk '{print $1}')"
eos_artifact_key="$(basename "$(dirname "$eos_wheel")")"
if [[ ! "$eos_artifact_key" =~ ^[0-9a-f]{64}$ ]] \
    || [[ "$eos_wheel_sha256" != "$eos_artifact_key" ]]; then
    echo "EOS wheel is not stored under its full SHA-256: $eos_wheel" >&2
    exit 1
fi

if [[ ! -x .venv/bin/python ]]; then
    uv venv --python 3.13 .venv
fi

uv pip uninstall --python .venv/bin/python \
    epcsaft \
    epcsaft-equilibrium \
    epcsaft-regression

uv pip install --python .venv/bin/python \
    "$eos_wheel" \
    mypy \
    pint \
    pybind11 \
    pytest \
    ruff \
    scikit-build-core
uv pip install --python .venv/bin/python --no-build-isolation --editable .

pybind11_cmake_dir="$(.venv/bin/python -m pybind11 --cmakedir)"
cmake -S . -B .codex-build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DPython_EXECUTABLE="$repo_root/.venv/bin/python" \
    -Dpybind11_DIR="$pybind11_cmake_dir"
cmake --build .codex-build -j2

.venv/bin/python - "$eos_wheel" "$eos_wheel_sha256" <<'PY'
from importlib import metadata
import json
from pathlib import Path
import sys
from urllib.parse import unquote, urlparse

expected_wheel = Path(sys.argv[1]).resolve()
expected_sha256 = sys.argv[2]
eos = metadata.distribution("epcsaft")
direct_url = json.loads(eos.read_text("direct_url.json") or "{}")
if eos.version != "0.2.0.dev0":
    raise SystemExit(f"epcsaft version mismatch: {eos.version}")
if direct_url.get("dir_info", {}).get("editable", False):
    raise SystemExit("editable epcsaft EOS install is forbidden")
installed_url = direct_url.get("url", "")
installed_wheel = Path(unquote(urlparse(installed_url).path)).resolve()
if installed_wheel != expected_wheel:
    raise SystemExit(
        f"installed EOS wheel mismatch: {installed_wheel} != {expected_wheel}"
    )
archive_info = direct_url.get("archive_info", {})
installed_sha256 = archive_info.get("hash", "").removeprefix("sha256=")
if installed_sha256 and installed_sha256 != expected_sha256:
    raise SystemExit(
        f"installed EOS hash mismatch: {installed_sha256} != {expected_sha256}"
    )
PY
