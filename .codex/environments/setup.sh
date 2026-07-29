#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

for required_command in uv cmake ninja pkg-config c++ serena; do
    if ! command -v "$required_command" >/dev/null 2>&1; then
        echo "missing required command: $required_command" >&2
        exit 1
    fi
done

if ! pkg-config --exists ipopt; then
    echo "missing required pkg-config package: ipopt" >&2
    exit 1
fi

skill_root="${HOME}/.agents/skills"
for required_skill in chemical-engineer scientific-coding-and-testing; do
    if [[ ! -f "${skill_root}/${required_skill}/SKILL.md" ]]; then
        echo "missing required Codex skill: ${required_skill}" >&2
        exit 1
    fi
done

git_common_dir="$(git rev-parse --path-format=absolute --git-common-dir)"
canonical_repo_root="$(dirname "$git_common_dir")"
project_root="$(dirname "$canonical_repo_root")"
default_provider_wheel="${project_root}/ePC-SAFT-EoS/dist/epcsaft-0.2.0.dev0-cp313-cp313-linux_x86_64.whl"
provider_wheel="${EPCSAFT_PROVIDER_WHEEL:-$default_provider_wheel}"
if [[ ! -f "$provider_wheel" ]]; then
    echo "missing required Provider wheel: $provider_wheel" >&2
    exit 1
fi

if [[ ! -x .venv/bin/python ]]; then
    uv venv --python 3.13 .venv
fi

uv pip install --python .venv/bin/python \
    "$provider_wheel" \
    mypy \
    pint \
    pybind11 \
    pytest \
    ruff \
    scikit-build-core
uv pip install --python .venv/bin/python --no-build-isolation --editable .

pybind11_cmake_dir="$(.venv/bin/python -m pybind11 --cmakedir)"
cmake -S . -B build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DPython_EXECUTABLE="$repo_root/.venv/bin/python" \
    -Dpybind11_DIR="$pybind11_cmake_dir"
cmake --build build -j2

.venv/bin/python - <<'PY'
from importlib import metadata
import json

provider = metadata.distribution("epcsaft")
direct_url = json.loads(provider.read_text("direct_url.json") or "{}")
if provider.version != "0.2.0.dev0":
    raise SystemExit(f"epcsaft version mismatch: {provider.version}")
if direct_url.get("dir_info", {}).get("editable", False):
    raise SystemExit("editable epcsaft provider is forbidden")
PY
