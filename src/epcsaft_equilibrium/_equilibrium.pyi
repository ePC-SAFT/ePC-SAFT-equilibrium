from typing import Any

def sdk_info(capsule: object) -> dict[str, object]: ...
def evaluate_phase(
    capsule: object,
    temperature_k: float,
    amount_mol: float,
    volume_m3: float,
    expected_fingerprint: str,
) -> dict[str, object]: ...
def evaluate_mixture_phase(
    capsule: object,
    temperature_k: float,
    amounts_mol: tuple[float, ...],
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
def _solve_tp_flash(
    capsule: object,
    temperature_k: float,
    pressure_pa: float,
    overall_mole_fractions: tuple[float, ...],
    expected_fingerprint: str,
    *,
    trace: bool = ...,
) -> dict[str, Any]: ...
def _chemical_equilibrium(
    capsule: object | None,
    spec: dict[str, object],
    source_standard_state: dict[str, object] | None,
    packing_fraction_bounds: tuple[float, float] | None,
    trace_floor: float,
    active_parameters: tuple[dict[str, object], ...] | None = ...,
) -> dict[str, Any]: ...
