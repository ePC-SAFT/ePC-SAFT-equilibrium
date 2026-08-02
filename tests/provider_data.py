from __future__ import annotations

import hashlib
import os
import subprocess
from pathlib import Path

import epcsaft

_DATA_COMMIT = "c096285415d4d3198b9d00fc75af48b837dd1305"
_PACKET_FINGERPRINTS = {
    "figiel-2025-reference-electrolytes": (
        "5ea86c44f67ded4c309e17f396b9414cfadcb2b0f00b2ac41f8c48c047af42e5"
    ),
    "gross-2001-methane-ethane": (
        "4b20a82c7ed8f34ec9a60ab547014b78064c3331ce636c49d32c10bfb7298ca7"
    ),
    "gross-2001-propane": (
        "b039bdf6b9622019c46c1cb7342c6c88b64714b8ac64f2dc3d7f413e1b5a720a"
    ),
}


def _packet_fingerprint(packet_root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(item for item in packet_root.rglob("*") if item.is_file()):
        relative = path.relative_to(packet_root).as_posix().encode()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(path.stat().st_size.to_bytes(8, "big"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


def packet_parameters(
    packet_id: str,
    components: tuple[str, ...],
) -> epcsaft.Parameters:
    default_root = Path(__file__).resolve().parents[2] / "ePC-SAFT-data"
    data_root = Path(os.environ.get("EPCSAFT_DATA_ROOT", default_root)).resolve()
    if packet_id not in _PACKET_FINGERPRINTS:
        raise AssertionError(f"public parameter packet is not locked: {packet_id}")
    data_commit = subprocess.run(
        ["git", "-C", str(data_root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if data_commit != _DATA_COMMIT:
        raise AssertionError("public parameter packet checkout commit changed")
    packet = data_root / "packets" / packet_id / "1"
    if _packet_fingerprint(packet) != _PACKET_FINGERPRINTS[packet_id]:
        raise AssertionError(f"public parameter packet bytes changed: {packet_id}")
    bundle = packet / "parameters"
    if not bundle.is_dir():
        raise AssertionError(f"required public parameter packet is unavailable: {packet_id}")
    return epcsaft.Parameters.from_bundle(bundle, components=components)
