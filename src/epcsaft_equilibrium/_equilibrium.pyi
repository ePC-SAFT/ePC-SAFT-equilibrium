from typing import Any

def sdk_info(capsule: object) -> dict[str, object]: ...
def evaluate_phase(
    capsule: object,
    temperature_k: float,
    amount_mol: float,
    volume_m3: float,
    expected_fingerprint: str,
) -> dict[str, object]: ...
def evaluate_nlp(
    capsule: object,
    temperature_k: float,
    expected_fingerprint: str,
    variables: tuple[float, float, float],
    multipliers: tuple[float, float, float],
) -> dict[str, object]: ...
def solve_saturation(
    capsule: object,
    temperature_k: float,
    expected_fingerprint: str,
    liquid_density_upper_mol_m3: float,
) -> dict[str, Any]: ...
def _solve_two_phase_flash(
    capsule: object,
    temperature_k: float,
    pressure_pa: float,
    overall_mole_fractions: tuple[float, float],
) -> dict[str, Any]: ...
