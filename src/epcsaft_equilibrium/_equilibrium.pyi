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
def _held2_audit_farkas_certificate(
    matrix: list[list[float]],
    row_lower: list[float],
    row_upper: list[float],
    column_lower: list[float],
    column_upper: list[float],
    row_ray: list[float],
) -> dict[str, Any]: ...
def _held2_audit_step5_kkt(
    variables: list[float],
    lower: list[float],
    upper: list[float],
    objective_gradient: list[float],
    constraints: list[list[float]],
    constraint_upper: list[float],
    lower_bound_multipliers: list[float],
    upper_bound_multipliers: list[float],
    constraint_multipliers: list[float],
    pullback_residual: float,
    pullback_scale: float,
    pressure_residual: float,
    same_major_iteration: bool,
    step4_binding_valid: bool,
    pressure_branch_valid: bool,
) -> dict[str, Any]: ...
def _held2_adjudicate_farkas_status(
    solver_infeasible: bool,
    has_certificate: bool,
    certificate_accepted: bool,
) -> str: ...
def _chemical_equilibrium(
    capsule: object | None,
    spec: dict[str, object],
    source_standard_state: dict[str, object] | None,
    packing_fraction_bounds: tuple[float, float] | None,
    trace_floor: float,
    active_parameters: tuple[dict[str, object], ...] | None = ...,
) -> dict[str, Any]: ...

def _chemical_observation_context(
    provider_capsule: object,
    rows: tuple[dict[str, object], ...],
    packing_fraction_bounds: tuple[float, float],
    trace_floor: float,
    parameter_templates: tuple[dict[str, object], ...],
    evaluator_identity: str,
    capability_id: str,
    capability_fingerprint: str,
    provider_artifact_identity: str,
    owner_artifact_identity: str,
    contract_fingerprint: str,
    artifact_identity: str,
) -> object: ...
def _chemical_observation_evaluate(
    context: object,
    parameter_values: tuple[float, ...],
    with_jacobian: bool = ...,
) -> dict[str, Any]: ...
