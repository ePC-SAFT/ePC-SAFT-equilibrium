from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from functools import cache
from importlib import metadata
from types import MappingProxyType
from typing import Any, cast

import epcsaft
from pint import Quantity
from pint.errors import DimensionalityError

from . import _equilibrium


@dataclass(frozen=True)
class PhaseState:
    """One explicit phase state returned by a bounded equilibrium route."""

    amount_mol: float
    mole_fractions: tuple[float, ...]
    volume_m3: float
    molar_density_mol_m3: float
    pressure_pa: float
    chemical_potential_over_rt: tuple[float, ...]


@dataclass(frozen=True)
class SolverAttemptDiagnostics:
    """One deterministic Ipopt search or confirmation attempt."""

    role: str
    initial_guess: tuple[float, ...]
    solver_converged: bool
    solver_status: str
    iterations: int
    constraint_violation: float
    callback_error: str


@dataclass(frozen=True)
class SaturationDiagnostics:
    """Separate numerical and physical evidence for one local boundary."""

    solver_converged: bool
    solver_status: str
    iterations: int
    attempts: int
    attempt_log: tuple[SolverAttemptDiagnostics, ...]
    solver_lower_bounds: tuple[float, float, float]
    solver_upper_bounds: tuple[float, float, float]
    solver_constraint_violation: float
    numerical_converged: bool
    confirmation_solves: int
    confirmation_max_relative_difference: float
    physical_accepted: bool
    pressure_relative_residual: float
    chemical_potential_absolute_residual: float
    phase_density_distance: float
    exact_derivatives: bool
    globality_certificate: bool
    failure_reason: str


@dataclass(frozen=True)
class SaturationResult:
    """Certified local pure-component saturation boundary in named SI fields."""

    temperature_k: float
    saturation_pressure_pa: float
    parameter_fingerprint: str
    vapor: PhaseState
    liquid: PhaseState
    diagnostics: SaturationDiagnostics


class SaturationError(RuntimeError):
    """Raised when no seed produces a numerically confirmed physical boundary."""

    def __init__(self, message: str, diagnostics: Mapping[str, object]) -> None:
        super().__init__(message)
        self.diagnostics = MappingProxyType(dict(diagnostics))


@dataclass(frozen=True)
class HeldStepTiming:
    """One native HELD2 step invocation and its measured work."""

    step: int
    invocation_count: int
    wall_seconds: float
    cpu_seconds: float
    provider_evaluations: int
    optimizer_solves: int
    optimizer_iterations: int
    terminal_status: str
    terminal_reason: str
    next_step: int


@dataclass(frozen=True)
class HeldPerformanceDiagnostics:
    """Measured HELD2 work and mechanism-based recovery evidence."""

    step_timings: tuple[HeldStepTiming, ...]
    provider_evaluations: int
    optimizer_solves: int
    optimizer_iterations: int
    step5_starts_consumed: int
    step5_attempts: int
    step5_accepted_attempts: int
    step5_repeated_start_ordinal_count: int
    dilute_face_restart_attempts: int
    dilute_face_restart_accepts: int
    dilute_face_coordinate_indices: tuple[int, ...]
    dilute_face_component_indices: tuple[int, ...]
    step8_invocations: int
    step8_return_to_stage_ii_count: int
    step8_warm_start_count: int
    step8_cold_fallback_count: int
    step8_provider_state_evaluations: int
    step8_provider_volume_bound_evaluations: int
    step8_provider_packing_evaluations: int
    step8_problem_candidate_count: int
    step8_attempted_candidate_count: int
    step8_repeated_problem_count: int
    step10_invocations: int
    trace_refinement_activated_count: int
    trace_refinement_return_to_stage_ii_count: int
    trace_refinement_component_indices: tuple[int, ...]


@dataclass(frozen=True)
class HeldDiagnostics:
    """Compact evidence from the bounded HELD search and local refinement."""

    outcome: str
    search_status: str
    solver_status: str
    numerical_status: str
    physical_status: str
    attempts: int
    major_iterations: int
    best_tpd: float
    lower_bound: float | None
    upper_bound: float | None
    held_gap: float | None
    material_balance_max_abs: float | None
    pressure_stationarity_max_relative: float | None
    kkt_stationarity_max_abs: float | None
    chemical_potential_max_relative: float | None
    confirmation_succeeded: bool
    confirmation_max_difference: float | None
    search_profiles: tuple[str, ...]
    globality_certificate: str
    failure_reason: str
    performance: HeldPerformanceDiagnostics | None = None


@dataclass(frozen=True)
class TpFlashResult:
    """Certified phase result from the bounded HELD or HELD2 controller."""

    temperature_k: float
    pressure_pa: float
    overall_mole_fractions: tuple[float, ...]
    phases: tuple[PhaseState, ...]
    phase_fractions: tuple[float, ...]
    total_free_energy_over_rt: float
    parameter_fingerprint: str
    diagnostics: HeldDiagnostics


class FlashError(RuntimeError):
    """Raised when the bounded HELD search cannot return an accepted result."""

    def __init__(self, message: str, diagnostics: HeldDiagnostics) -> None:
        super().__init__(message)
        self.diagnostics = diagnostics


@dataclass(frozen=True)
class ChemicalEquilibriumConstant:
    """One dimensionless reaction constant and its immutable provenance."""

    ln_value: float
    source_id: str
    reference_id: str
    reaction_orientation: str
    conversion_id: str
    dimensionless: bool


@dataclass(frozen=True)
class ChemicalStandardState:
    """Source activity convention and p° provenance transformed at each trial P."""

    id: str
    activity_scale_id: str
    log_activity_scale_factors: tuple[float, ...]
    reference_pressure_pa: float
    source_reference_component_ids: tuple[str, ...]
    source_reference_solvent_composition: tuple[float, ...]
    source_reference_activity_convention_id: str
    source_reference_standard_molality_mol_per_kg: float

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "log_activity_scale_factors",
            tuple(self.log_activity_scale_factors),
        )
        object.__setattr__(
            self,
            "source_reference_component_ids",
            tuple(self.source_reference_component_ids),
        )
        object.__setattr__(
            self,
            "source_reference_solvent_composition",
            tuple(self.source_reference_solvent_composition),
        )


@dataclass(frozen=True)
class IdealGasPhase:
    """Source-bound ideal-gas restriction used for nonpredictive benchmarks."""

    model_fingerprint: str
    reference_id: str


@dataclass(frozen=True)
class ProviderPhase:
    """Installed Provider model and application-admitted physical packing domain."""

    model: epcsaft.Mixture
    expected_parameter_fingerprint: str
    admissible_packing_fraction_interval: tuple[float, float]


@dataclass(frozen=True)
class ProviderModelContinuation:
    """A caller-supplied Provider endpoint for deterministic model homotopy."""

    initial_phase: ProviderPhase
    initial_equilibrium_constants: tuple[ChemicalEquilibriumConstant, ...] | None = None
    initial_equilibrium_constants_model_fingerprint: str | None = None


@dataclass(frozen=True)
class ChemicalEquilibriumProblem:
    """One reaction system with an explicit strict-interior admission boundary."""

    species_ids: tuple[str, ...]
    charges: tuple[int, ...]
    molar_masses_kg_per_mol: tuple[float, ...]
    balance_matrix: tuple[tuple[float, ...], ...]
    conserved_totals: tuple[float, ...]
    reaction_matrix: tuple[tuple[float, ...], ...]
    feed_amounts_mol: tuple[float, ...]
    equilibrium_constants: tuple[ChemicalEquilibriumConstant, ...]
    strict_interior_amount_floor_mol: float
    source_standard_state: ChemicalStandardState | None = None


@dataclass(frozen=True)
class ChemicalEquilibriumAttempt:
    """One deterministic primary or recovery terminal in the finite search."""

    ordinal: int
    primary_ordinal: int
    kind: str
    parent_ordinal: int | None
    start_identity: str
    start_construction_status: str
    retraction_status: str
    continuation_status: str
    provider_domain_status: str
    solver_status: str
    callback_error: str
    terminal_status: str
    amounts_mol: tuple[float, ...]
    volume_m3: float | None
    objective: float | None
    balance_inf_norm: float | None
    charge_inf_norm: float | None
    pressure_relative_residual: float | None
    reaction_affinity_inf_norm: float | None
    kkt_stationarity_inf_norm: float | None
    complementarity_inf_norm: float | None
    kkt_dimension: int
    kkt_rank: int
    condition_number_inf: float | None
    local_minimum_status: str
    trace_status: str
    basin_ordinal: int | None
    recovery_seed_count: int
    recovery_solve_count: int


@dataclass(frozen=True)
class ChemicalEquilibriumBasin:
    """One materially distinct certified strict local minimum."""

    ordinal: int
    representative_attempt_ordinal: int
    amounts_mol: tuple[float, ...]
    volume_m3: float
    objective: float


@dataclass(frozen=True)
class ChemicalEquilibriumBudgetPrefix:
    """A projection of the fixed run at one nested primary-start budget."""

    primary_budget: int
    attempted_primary_ordinals: tuple[int, ...]
    basin_ordinals: tuple[int, ...]
    selected_basin_ordinal: int | None
    selection_changed: bool


@dataclass(frozen=True)
class ChemicalEquilibriumSearch:
    """Complete deterministic finite-basin search receipt."""

    status: str
    continuation_status: str
    continuation_blocker: str
    continuation_initial_model_fingerprint: str
    continuation_accepted_lambda: float | None
    continuation_attempt_count: int
    primary_budget: int
    primary_attempt_count: int
    generated_start_count: int
    budget_truncated_start_count: int
    duplicate_start_count: int
    infeasible_start_count: int
    evaluated_start_count: int
    domain_rejected_start_count: int
    construction_rejected_start_count: int
    attempts: tuple[ChemicalEquilibriumAttempt, ...]
    basins: tuple[ChemicalEquilibriumBasin, ...]
    budget_prefixes: tuple[ChemicalEquilibriumBudgetPrefix, ...]
    selected_basin_ordinal: int | None
    selected_objective: float | None
    selection_label: str


@dataclass(frozen=True)
class ChemicalEquilibriumCriterion:
    """One named scalar certification check with its declared threshold."""

    name: str
    status: str
    value: float
    limit: float
    comparison: str


@dataclass(frozen=True)
class ChemicalEquilibriumDiagnostics:
    """Independent local chemical-equilibrium admission axes and residuals."""

    solver_status: str
    callback_error: str
    failure_kind: str
    failure_reason: str
    chemical_certification_level: str
    boundary_status: str
    structural_zero_species_indices: tuple[int, ...]
    numerical_status: str
    physical_status: str
    predictive_status: str
    provider_domain_status: str
    local_minimum_status: str
    negative_curvature_recovery_status: str
    negative_curvature_recovery_attempts: int
    negative_curvature_recovery_selected_sign: int
    trace_status: str
    globality_status: str
    search: ChemicalEquilibriumSearch
    amounts_mol: tuple[float, ...]
    volume_m3: float | None
    balance_inf_norm: float | None
    charge_inf_norm: float | None
    pressure_relative_residual: float | None
    reaction_affinity_inf_norm: float | None
    reaction_affinity_residuals: tuple[float, ...]
    packing_fraction: float | None
    packing_fraction_min: float | None
    packing_fraction_max: float | None
    total_ion_fraction: float | None
    total_ion_fraction_max: float | None
    minimum_amount_mol: float | None
    trace_floor_mol: float | None
    kkt_stationarity_inf_norm: float | None
    physical_stationarity_residuals: tuple[float, ...]
    physical_equality_multipliers: tuple[float, ...]
    physical_constraint_jacobian: tuple[float, ...]
    physical_constraint_shape: tuple[int, int]
    physical_lagrangian_gradient: tuple[float, ...]
    physical_to_chart_jacobian: tuple[float, ...]
    physical_to_chart_shape: tuple[int, int]
    chart_physical_pullback_residual_inf_norm: float | None
    complementarity_inf_norm: float | None
    numerical_criteria: tuple[ChemicalEquilibriumCriterion, ...]
    physical_criteria: tuple[ChemicalEquilibriumCriterion, ...]
    first_failed_numerical_criterion: str | None
    first_failed_physical_criterion: str | None
    kkt_dimension: int
    kkt_rank: int
    condition_number_inf: float | None
    active_lower_bounds: tuple[int, ...]
    active_upper_bounds: tuple[int, ...]
    active_constraint_bounds: tuple[int, ...]
    active_trace_species: tuple[int, ...]
    chart_topology: str
    reduced_hessian_status: str
    reduced_hessian: tuple[float, ...]
    reduced_hessian_nullspace_basis: tuple[float, ...]
    reduced_hessian_nullspace_shape: tuple[int, int]
    reduced_hessian_eigenvalues: tuple[float, ...]
    reduced_hessian_spectrum_status: str
    reduced_hessian_raw_inertia: tuple[int, int, int]
    reduced_hessian_inertia: tuple[int, int, int]
    reduced_hessian_scale: float | None
    reduced_hessian_eigenvalue_tolerance: float | None
    objective_gradient: tuple[float, ...]
    constraint_values: tuple[float, ...]
    constraint_jacobian: tuple[float, ...]
    lagrangian_gradient: tuple[float, ...]
    equality_multipliers: tuple[float, ...]
    chart_stationarity_inf_norm: float | None
    lagrangian_hessian: tuple[float, ...]
    covariant_lagrangian_hessian: tuple[float, ...]
    derivative_check_step: float | None
    objective_gradient_check_relative_error: float | None
    constraint_jacobian_check_relative_error: float | None
    lagrangian_hessian_check_relative_error: float | None
    derivative_check_worst_entry: str
    derivative_check_worst_relative_error: float | None
    derivative_check_worst_analytic_value: float | None
    derivative_check_worst_finite_difference_value: float | None
    derivative_check_worst_step: float | None
    derivative_coordinate_order: tuple[str, ...]
    derivative_objective_basis: str
    derivative_constraint_basis: str
    derivative_constraint_order: tuple[str, ...]
    objective_gradient_shape: tuple[int, ...]
    constraint_jacobian_shape: tuple[int, int]
    lagrangian_hessian_shape: tuple[int, int]
    kkt_root_jacobian: tuple[float, ...]
    kkt_root_shape: tuple[int, int]
    kkt_root_status: str
    kkt_root_jacobian_check_relative_error: float | None
    reference_representation_residual_inf_norm: float | None
    reference_convergence_error: float | None
    reference_derivative_availability: int | None


@dataclass(frozen=True)
class ChemicalEquilibriumActiveParameter:
    """One exact coordinate resolved from the installed Provider descriptor table."""

    family: str
    identity: str
    component_ids: tuple[str, ...]
    value: float
    unit: str


@dataclass(frozen=True)
class ChemicalEquilibriumSensitivityRequest:
    """Request exact operation columns and ordered installed-Provider coordinates."""

    active_parameters: tuple[ChemicalEquilibriumActiveParameter, ...] = ()


@dataclass(frozen=True)
class ChemicalObservationPrimitive:
    """One generic homogeneous-liquid observable evaluated by the installed artifacts."""

    kind: str
    component_id: str


@dataclass(frozen=True)
class ChemicalObservationRow:
    """One source-bound row and immutable homogeneous-liquid state."""

    row_id: str
    state_id: str
    state_schema_id: str
    source_id: str
    transform_id: str
    temperature: Quantity[Any]
    pressure: Quantity[Any]
    problem: ChemicalEquilibriumProblem
    primitive: ChemicalObservationPrimitive


@dataclass(frozen=True)
class ChemicalObservationContext:
    """Process-local C evaluator handle with value and exact-Jacobian transport."""

    _native_capsule: object
    row_ids: tuple[str, ...]
    state_ids: tuple[str, ...]
    state_schema_ids: tuple[str, ...]
    source_ids: tuple[str, ...]
    primitive_ids: tuple[str, ...]
    parameter_ids: tuple[str, ...]
    primitive_units: tuple[str, ...]
    transform_ids: tuple[str, ...]
    reference_ids: tuple[str, ...]
    reference_fingerprints: tuple[str, ...]
    parameter_units: tuple[str, ...]
    evaluator_identity: str
    capability_id: str
    capability_fingerprint: str
    provider_artifact_identity: str
    owner_artifact_identity: str
    contract_fingerprint: str
    artifact_identity: str

    @property
    def native_capsule(self) -> object:
        """Return the installed generic evaluator capsule for C/C++ consumers."""

        return self._native_capsule

    def evaluate(
        self,
        parameter_values: Sequence[float],
        *,
        with_jacobian: bool = True,
    ) -> Mapping[str, object]:
        """Replay the same native callback for diagnostics or installed-artifact tests."""

        values = tuple(float(value) for value in parameter_values)
        return cast(
            Mapping[str, object],
            _equilibrium._chemical_observation_evaluate(
                self._native_capsule,
                values,
                with_jacobian,
            ),
        )


@dataclass(frozen=True)
class ChemicalEquilibriumSensitivityParameter:
    """One ordered scalar input and the units of its returned state derivatives."""

    name: str
    kind: str
    source_index: int | None
    input_unit: str
    amount_derivative_unit: str
    volume_derivative_unit: str
    identity: str | None = None
    component_ids: tuple[str, ...] = ()
    value: float | None = None


@dataclass(frozen=True)
class ChemicalEquilibriumSensitivity:
    """Exact conditioned state-input Jacobian or fail-closed unavailability evidence."""

    status: str
    failure_reason: str
    parameters: tuple[ChemicalEquilibriumSensitivityParameter, ...]
    amount_state_order: tuple[str, ...]
    amount_derivatives: tuple[tuple[float, ...], ...]
    volume_derivatives: tuple[float, ...]
    kkt_dimension: int
    kkt_rank: int
    condition_number_inf: float
    active_lower_bounds: tuple[int, ...]
    active_upper_bounds: tuple[int, ...]
    active_constraint_bounds: tuple[int, ...]
    active_trace_species: tuple[int, ...]
    chart_topology: str
    parameter_fingerprint: str
    provider_parameter_status: str
    provider_parameter_failure_reason: str
    reference_parameter_status: str
    reference_parameter_failure_reason: str


@dataclass(frozen=True)
class ChemicalArtifactIdentity:
    """Installed distributions and Provider ABI consumed by one chemical result."""

    equilibrium_distribution: str
    equilibrium_version: str
    equilibrium_record_sha256: str
    provider_distribution: str | None
    provider_version: str | None
    provider_record_sha256: str | None
    provider_sdk_capsule_name: str | None
    provider_sdk_abi_version: int | None
    provider_sdk_table_size: int | None
    provider_sdk_result_size: int | None
    provider_sdk_mixture_result_size: int | None
    provider_sdk_neutral_reference_result_size: int | None
    provider_sdk_neutral_reference_derivative_result_size: int | None
    provider_sdk_reacting_phase_parameter_result_size: int | None


@dataclass(frozen=True)
class ChemicalEquilibriumResult:
    """Certified local fixed-T,P homogeneous chemical-equilibrium value."""

    temperature_k: float
    pressure_pa: float
    species_ids: tuple[str, ...]
    charges: tuple[int, ...]
    amounts_mol: tuple[float, ...]
    mole_fractions: tuple[float, ...]
    volume_m3: float
    thermodynamic_model: str
    model_fingerprint: str
    provider_component_ids: tuple[str, ...] | None
    provider_parameter_fingerprint: str | None
    equilibrium_constants: tuple[ChemicalEquilibriumConstant, ...]
    source_standard_state: ChemicalStandardState | None
    source_reference_transfer: epcsaft.SourceReferenceTransfer | None
    provider_reference_id: str | None
    standard_offsets: tuple[float, ...] | None
    ln_k_provider_basis: tuple[float, ...] | None
    local_scope: str
    diagnostics: ChemicalEquilibriumDiagnostics
    response_kind: str
    artifact_identity: ChemicalArtifactIdentity
    sensitivity: ChemicalEquilibriumSensitivity | None


class ChemicalEquilibriumError(RuntimeError):
    """Raised when no typed, certified local chemical state can be admitted."""

    def __init__(self, message: str, diagnostics: ChemicalEquilibriumDiagnostics) -> None:
        super().__init__(message)
        self.diagnostics = diagnostics


@dataclass(frozen=True)
class _Scope:
    component: str
    temperature_min_k: float
    temperature_max_k: float
    liquid_density_upper_mol_m3: float


_SCOPES = MappingProxyType(
    {
        "sha256:4a5cde0fad05a150c7fa54bb5ac7db508424f0126cf512596dcc248c284dd9e0": _Scope(
            component="methane",
            temperature_min_k=97.0,
            temperature_max_k=300.0,
            liquid_density_upper_mol_m3=40_000.0,
        ),
        "sha256:6381ef30d4f25c63fee1fd098c7024dd59254f3021d0588ee46e4eecfb31619b": _Scope(
            component="ethane",
            temperature_min_k=90.0,
            temperature_max_k=305.0,
            liquid_density_upper_mol_m3=40_000.0,
        ),
        "sha256:7d893c35288fceec76b7dde4a16fbf2ef95830b7d18827d3205537ce54878140": _Scope(
            component="propane",
            temperature_min_k=85.0,
            temperature_max_k=523.0,
            liquid_density_upper_mol_m3=25_000.0,
        ),
    }
)

_FLASH_FINGERPRINT = "sha256:8347d8daad42af60d61071f0584eb50d8866d98d9636872fd9d173f491ea7947"
_FLASH_TEMPERATURE_DOMAIN_K = (203.22, 243.61)
_FLASH_PRESSURE_DOMAIN_PA = (2_124_000.0, 6_885_000.0)
_FLASH_METHANE_FEED_DOMAIN = (0.4661, 0.66705)


def _temperature_in_kelvin(temperature: Quantity[Any]) -> float:
    if not isinstance(temperature, Quantity):
        raise TypeError("saturation requires a Pint temperature quantity")
    try:
        value = float(temperature.to("kelvin").magnitude)
    except DimensionalityError as error:
        raise ValueError("temperature units must be convertible to kelvin") from error
    except (TypeError, ValueError) as error:
        raise TypeError("saturation requires one scalar Pint temperature quantity") from error
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("temperature must be positive and finite")
    return value


def _triple(payload: object) -> tuple[float, float, float]:
    values = cast(list[float] | tuple[float, ...], payload)
    if len(values) != 3:
        raise RuntimeError("native solver diagnostic triple has the wrong size")
    return (float(values[0]), float(values[1]), float(values[2]))


def _float_tuple(payload: object) -> tuple[float, ...]:
    if isinstance(payload, (int, float)):
        return (float(payload),)
    values = cast(list[float] | tuple[float, ...], payload)
    return tuple(float(value) for value in values)


def _vector(payload: object, size: int, name: str) -> tuple[float, ...]:
    values = _float_tuple(payload)
    if len(values) != size:
        raise RuntimeError(f"native {name} has the wrong size")
    return values


def _attempt(payload: Mapping[str, object]) -> SolverAttemptDiagnostics:
    return SolverAttemptDiagnostics(
        role=str(payload["role"]),
        initial_guess=_triple(payload["initial_guess"]),
        solver_converged=bool(payload["solver_converged"]),
        solver_status=str(payload["solver_status"]),
        iterations=int(cast(int, payload["iterations"])),
        constraint_violation=float(cast(float, payload["constraint_violation"])),
        callback_error=str(payload["callback_error"]),
    )


def _diagnostics(payload: Mapping[str, object]) -> SaturationDiagnostics:
    attempt_payloads = cast(list[Mapping[str, object]], payload["attempt_log"])
    return SaturationDiagnostics(
        solver_converged=bool(payload["solver_converged"]),
        solver_status=str(payload["solver_status"]),
        iterations=int(cast(int, payload["iterations"])),
        attempts=int(cast(int, payload["attempts"])),
        attempt_log=tuple(_attempt(item) for item in attempt_payloads),
        solver_lower_bounds=_triple(payload["solver_lower_bounds"]),
        solver_upper_bounds=_triple(payload["solver_upper_bounds"]),
        solver_constraint_violation=float(cast(float, payload["solver_constraint_violation"])),
        numerical_converged=bool(payload["numerical_converged"]),
        confirmation_solves=int(cast(int, payload["confirmation_solves"])),
        confirmation_max_relative_difference=float(
            cast(float, payload["confirmation_max_relative_difference"])
        ),
        physical_accepted=bool(payload["physical_accepted"]),
        pressure_relative_residual=float(cast(float, payload["pressure_relative_residual"])),
        chemical_potential_absolute_residual=float(
            cast(float, payload["chemical_potential_absolute_residual"])
        ),
        phase_density_distance=float(cast(float, payload["phase_density_distance"])),
        exact_derivatives=bool(payload["exact_derivatives"]),
        globality_certificate=bool(payload["globality_certificate"]),
        failure_reason=str(payload["failure_reason"]),
    )


def _phase(payload: Mapping[str, object]) -> PhaseState:
    mole_fractions = _float_tuple(payload.get("mole_fractions", (1.0,)))
    chemical_potential = _vector(
        payload["chemical_potential_over_rt"],
        len(mole_fractions),
        "phase chemical-potential vector",
    )
    return PhaseState(
        amount_mol=float(cast(float, payload["amount_mol"])),
        mole_fractions=mole_fractions,
        volume_m3=float(cast(float, payload["volume_m3"])),
        molar_density_mol_m3=float(cast(float, payload["molar_density_mol_m3"])),
        pressure_pa=float(cast(float, payload["pressure_pa"])),
        chemical_potential_over_rt=chemical_potential,
    )


def _quantity(quantity: object, units: str, name: str, operation: str) -> float:
    if not isinstance(quantity, Quantity):
        raise TypeError(f"{operation} requires a Pint {name} quantity")
    try:
        value = float(quantity.to(units).magnitude)
    except DimensionalityError as error:
        raise ValueError(f"{name} units must be convertible to {units}") from error
    except (TypeError, ValueError) as error:
        raise TypeError(f"{operation} requires one scalar Pint {name} quantity") from error
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return value


def _tp_flash_feed(
    overall_mole_fractions: Sequence[float],
    component_count: int,
) -> tuple[float, ...]:
    if isinstance(overall_mole_fractions, (str, bytes)):
        raise TypeError("overall mole fractions must be a numeric sequence")
    try:
        values = tuple(float(value) for value in overall_mole_fractions)
    except (TypeError, ValueError) as error:
        raise TypeError("overall mole fractions must be a numeric sequence") from error
    if len(values) != component_count:
        raise ValueError(f"tp_flash requires exactly {component_count} components")
    if not all(math.isfinite(value) and value > 0.0 for value in values):
        raise ValueError("overall mole fractions must be positive and finite")
    if not math.isclose(sum(values), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError("overall mole fractions must sum to one within 1e-12")
    return values


def _optional_float(payload: Mapping[str, object], name: str) -> float | None:
    value = payload[name]
    return None if value is None else float(cast(float, value))


def _held_diagnostics(payload: Mapping[str, object]) -> HeldDiagnostics:
    profiles = cast(Sequence[object], payload["search_profiles"])
    statuses = tuple(
        payload[name] for name in ("solver_status", "numerical_status", "physical_status")
    )
    if not all(
        isinstance(status, str) and status in {"passed", "failed", "not_adjudicated"}
        for status in statuses
    ):
        raise ValueError("native HELD payload has an invalid status vocabulary")
    typed_statuses = cast(tuple[str, str, str], statuses)
    return HeldDiagnostics(
        outcome=str(payload["outcome"]),
        search_status=str(payload["search_status"]),
        solver_status=typed_statuses[0],
        numerical_status=typed_statuses[1],
        physical_status=typed_statuses[2],
        attempts=int(cast(int, payload["attempts"])),
        major_iterations=int(cast(int, payload["major_iterations"])),
        best_tpd=float(cast(float, payload["best_tpd"])),
        lower_bound=_optional_float(payload, "lower_bound"),
        upper_bound=_optional_float(payload, "upper_bound"),
        held_gap=_optional_float(payload, "held_gap"),
        material_balance_max_abs=_optional_float(payload, "material_balance_max_abs"),
        pressure_stationarity_max_relative=_optional_float(
            payload, "pressure_stationarity_max_relative"
        ),
        kkt_stationarity_max_abs=_optional_float(payload, "kkt_stationarity_max_abs"),
        chemical_potential_max_relative=_optional_float(payload, "chemical_potential_max_relative"),
        confirmation_succeeded=bool(payload["confirmation_succeeded"]),
        confirmation_max_difference=_optional_float(payload, "confirmation_max_difference"),
        search_profiles=tuple(str(profile) for profile in profiles),
        globality_certificate=str(payload["globality_certificate"]),
        failure_reason=str(payload["failure_reason"]),
    )


def _failed_held_diagnostics(outcome: str, search_status: str, reason: str) -> HeldDiagnostics:
    return HeldDiagnostics(
        outcome=outcome,
        search_status=search_status,
        solver_status="not_adjudicated",
        numerical_status="not_adjudicated",
        physical_status="not_adjudicated",
        attempts=0,
        major_iterations=0,
        best_tpd=0.0,
        lower_bound=None,
        upper_bound=None,
        held_gap=None,
        material_balance_max_abs=None,
        pressure_stationarity_max_relative=None,
        kkt_stationarity_max_abs=None,
        chemical_potential_max_relative=None,
        confirmation_succeeded=False,
        confirmation_max_difference=None,
        search_profiles=(),
        globality_certificate="not_guaranteed",
        failure_reason=reason,
    )


def _held2_performance(payload: Mapping[str, object]) -> HeldPerformanceDiagnostics:
    timing_payloads = cast(Sequence[object], payload.get("step_timings", ()))
    timings: list[HeldStepTiming] = []
    for raw_timing in timing_payloads:
        if not isinstance(raw_timing, Mapping):
            raise ValueError("native HELD2 step timing must be a mapping")
        timing = cast(Mapping[str, object], raw_timing)
        step = int(cast(int, timing["step"]))
        invocation_count = int(cast(int, timing["invocation_count"]))
        wall_seconds = float(cast(float, timing["wall_seconds"]))
        cpu_seconds = float(cast(float, timing["cpu_seconds"]))
        provider_evaluations = int(cast(int, timing["provider_evaluations"]))
        optimizer_solves = int(cast(int, timing["optimizer_solves"]))
        optimizer_iterations = int(cast(int, timing["optimizer_iterations"]))
        next_step = int(cast(int, timing["next_step"]))
        if (
            step not in range(1, 11)
            or invocation_count < 0
            or not math.isfinite(wall_seconds)
            or wall_seconds < 0.0
            or not math.isfinite(cpu_seconds)
            or cpu_seconds < 0.0
            or min(provider_evaluations, optimizer_solves, optimizer_iterations) < 0
            or next_step not in range(0, 11)
        ):
            raise ValueError("native HELD2 step timing is invalid")
        timings.append(
            HeldStepTiming(
                step=step,
                invocation_count=invocation_count,
                wall_seconds=wall_seconds,
                cpu_seconds=cpu_seconds,
                provider_evaluations=provider_evaluations,
                optimizer_solves=optimizer_solves,
                optimizer_iterations=optimizer_iterations,
                terminal_status=str(timing["terminal_status"]),
                terminal_reason=str(timing["terminal_reason"]),
                next_step=next_step,
            )
        )

    step5_starts_consumed = 0
    step5_attempts = 0
    step5_accepted_attempts = 0
    dilute_face_restart_attempts = 0
    dilute_face_restart_accepts = 0
    dilute_face_coordinate_indices: set[int] = set()
    dilute_face_component_indices: set[int] = set()
    start_ordinals: list[int] = []
    for raw_step5 in cast(Sequence[object], payload.get("step5_history", ())):
        if not isinstance(raw_step5, Mapping):
            raise ValueError("native HELD2 Step 5 history must contain mappings")
        step5 = cast(Mapping[str, object], raw_step5)
        starts_consumed = int(cast(int, step5.get("starts_consumed", 0)))
        if starts_consumed < 0:
            raise ValueError("native HELD2 Step 5 starts consumed is invalid")
        step5_starts_consumed += starts_consumed
        for raw_attempt in cast(Sequence[object], step5.get("attempts", ())):
            if not isinstance(raw_attempt, Mapping):
                raise ValueError("native HELD2 Step 5 attempt must be a mapping")
            attempt = cast(Mapping[str, object], raw_attempt)
            step5_attempts += 1
            step5_accepted_attempts += int(bool(attempt.get("accepted", False)))
            if "start_ordinal" in attempt:
                ordinal = int(cast(int, attempt["start_ordinal"]))
                if ordinal < 0:
                    raise ValueError("native HELD2 Step 5 start ordinal is invalid")
                start_ordinals.append(ordinal)
            raw_restart = attempt.get("dilute_face_restart")
            if not isinstance(raw_restart, Mapping):
                continue
            restart = cast(Mapping[str, object], raw_restart)
            if not bool(restart.get("attempted", False)):
                continue
            dilute_face_restart_attempts += 1
            dilute_face_restart_accepts += int(bool(restart.get("accepted", False)))
            for raw_coordinate in cast(Sequence[object], restart.get("coordinate_indices", ())):
                coordinate = int(cast(int, raw_coordinate))
                if coordinate < 0:
                    raise ValueError("native HELD2 dilute-face coordinate is invalid")
                dilute_face_coordinate_indices.add(coordinate)
            for raw_component in cast(
                Sequence[object], restart.get("provider_component_indices", ())
            ):
                component = int(cast(int, raw_component))
                if component < 0:
                    raise ValueError("native HELD2 dilute-face component is invalid")
                dilute_face_component_indices.add(component)

    return_to_stage_ii_reasons = {
        "problem_67_not_converged",
        "stage_iii_solver_not_converged",
        "stage_iii_active_set_resolve_failed",
        "stage_iii_active_set_balance_failed",
    }
    step8_invocations = 0
    step8_return_to_stage_ii_count = 0
    step8_warm_start_count = 0
    step8_cold_fallback_count = 0
    step8_provider_state_evaluations = 0
    step8_provider_volume_bound_evaluations = 0
    step8_provider_packing_evaluations = 0
    step8_problem_candidate_count = 0
    step8_attempted_candidate_count = 0
    problem_signatures: list[tuple[tuple[int, ...], tuple[float, ...], float]] = []
    for raw_step8 in cast(Sequence[object], payload.get("step8_history", ())):
        if not isinstance(raw_step8, Mapping):
            raise ValueError("native HELD2 Step 8 history must contain mappings")
        step8 = cast(Mapping[str, object], raw_step8)
        step8_invocations += 1
        step8_return_to_stage_ii_count += int(
            str(step8.get("reason", "")) in return_to_stage_ii_reasons
        )
        step8_warm_start_count += int(bool(step8.get("warm_start_used", False)))
        step8_cold_fallback_count += int(bool(step8.get("cold_fallback_used", False)))
        step8_provider_state_evaluations += int(
            cast(int, step8.get("provider_state_evaluations", 0))
        )
        step8_provider_volume_bound_evaluations += int(
            cast(int, step8.get("provider_volume_bound_evaluations", 0))
        )
        step8_provider_packing_evaluations += int(
            cast(int, step8.get("provider_packing_evaluations", 0))
        )
        if (
            min(
                step8_provider_state_evaluations,
                step8_provider_volume_bound_evaluations,
                step8_provider_packing_evaluations,
            )
            < 0
        ):
            raise ValueError("native HELD2 Step 8 Provider work is invalid")
        problem_ids = tuple(
            int(cast(int, value))
            for value in cast(Sequence[object], step8.get("problem_candidate_ids", ()))
        )
        attempted_ids = tuple(
            int(cast(int, value))
            for value in cast(Sequence[object], step8.get("attempted_candidate_ids", ()))
        )
        problem_variables = tuple(
            float(cast(float, value))
            for value in cast(Sequence[object], step8.get("problem_candidate_variables", ()))
        )
        if "neighborhood_radius" not in step8:
            raise ValueError("native HELD2 Step 8 neighborhood radius is missing")
        neighborhood_radius = float(cast(float, step8["neighborhood_radius"]))
        if not math.isfinite(neighborhood_radius) or neighborhood_radius <= 0.0:
            raise ValueError("native HELD2 Step 8 neighborhood radius is invalid")
        step8_problem_candidate_count += len(problem_ids)
        step8_attempted_candidate_count += len(attempted_ids)
        problem_signatures.append((problem_ids, problem_variables, neighborhood_radius))

    step10_timings = tuple(timing for timing in timings if timing.step == 10)
    trace_refinement_component_indices: set[int] = set()
    trace_refinement_activated_count = 0
    raw_step10_history = payload.get("step10_history", ())
    step10_payloads = (
        cast(Sequence[object], raw_step10_history)
        if isinstance(raw_step10_history, Sequence)
        else ()
    )
    if not step10_payloads and isinstance(payload.get("step10"), Mapping):
        step10_payloads = (cast(Mapping[str, object], payload["step10"]),)
    for raw_step10 in step10_payloads:
        if not isinstance(raw_step10, Mapping):
            raise ValueError("native HELD2 Step 10 history must contain mappings")
        step10 = cast(Mapping[str, object], raw_step10)
        refinement_activated = False
        raw_components = step10.get("trace_component_indices")
        if raw_components is None:
            raw_components = tuple(
                cast(Mapping[str, object], refinement)["component_index"]
                for refinement in cast(Sequence[object], step10.get("refinements", ()))
                if isinstance(refinement, Mapping)
            )
        for raw_component in cast(Sequence[object], raw_components):
            component = int(cast(int, raw_component))
            if component < 0:
                raise ValueError("native HELD2 trace component index is invalid")
            trace_refinement_component_indices.add(component)
            refinement_activated = True
        for raw_refinement in cast(Sequence[object], step10.get("refinements", ())):
            if not isinstance(raw_refinement, Mapping):
                raise ValueError("native HELD2 trace refinement must be a mapping")
        trace_refinement_activated_count += int(refinement_activated)

    return HeldPerformanceDiagnostics(
        step_timings=tuple(timings),
        provider_evaluations=sum(timing.provider_evaluations for timing in timings),
        optimizer_solves=sum(timing.optimizer_solves for timing in timings),
        optimizer_iterations=sum(timing.optimizer_iterations for timing in timings),
        step5_starts_consumed=step5_starts_consumed,
        step5_attempts=step5_attempts,
        step5_accepted_attempts=step5_accepted_attempts,
        step5_repeated_start_ordinal_count=(len(start_ordinals) - len(set(start_ordinals))),
        dilute_face_restart_attempts=dilute_face_restart_attempts,
        dilute_face_restart_accepts=dilute_face_restart_accepts,
        dilute_face_coordinate_indices=tuple(sorted(dilute_face_coordinate_indices)),
        dilute_face_component_indices=tuple(sorted(dilute_face_component_indices)),
        step8_invocations=step8_invocations,
        step8_return_to_stage_ii_count=step8_return_to_stage_ii_count,
        step8_warm_start_count=step8_warm_start_count,
        step8_cold_fallback_count=step8_cold_fallback_count,
        step8_provider_state_evaluations=step8_provider_state_evaluations,
        step8_provider_volume_bound_evaluations=(step8_provider_volume_bound_evaluations),
        step8_provider_packing_evaluations=step8_provider_packing_evaluations,
        step8_problem_candidate_count=step8_problem_candidate_count,
        step8_attempted_candidate_count=step8_attempted_candidate_count,
        step8_repeated_problem_count=(len(problem_signatures) - len(set(problem_signatures))),
        step10_invocations=sum(timing.invocation_count for timing in step10_timings),
        trace_refinement_activated_count=trace_refinement_activated_count,
        trace_refinement_return_to_stage_ii_count=sum(
            timing.invocation_count for timing in step10_timings if timing.next_step == 4
        ),
        trace_refinement_component_indices=tuple(sorted(trace_refinement_component_indices)),
    )


def _held2_diagnostics(payload: Mapping[str, object]) -> HeldDiagnostics:
    outcome = str(payload["outcome"])
    accepted = outcome == "physical_equilibrium_accepted"
    one_phase = outcome == "one_phase_no_negative_witness_detected"
    step2 = cast(Mapping[str, object], payload["step2"])
    step9_history = cast(Sequence[Mapping[str, object]], payload["step9_history"])
    step10_value = payload.get("step10")
    step10 = cast(Mapping[str, object], step10_value) if isinstance(step10_value, Mapping) else None
    certificate_value = step10.get("final_certificate") if step10 else None
    certificate = (
        cast(Mapping[str, object], certificate_value)
        if isinstance(certificate_value, Mapping)
        else None
    )
    failure = payload.get("failure_reason")
    gap = step9_history[-1].get("free_energy_gap") if step9_history else None
    status = "passed" if accepted or one_phase else "failed"
    return HeldDiagnostics(
        outcome=outcome,
        search_status=str(payload.get("failure_stage") or outcome),
        solver_status=status,
        numerical_status=status,
        physical_status=status,
        attempts=int(cast(int, payload["upper_solve_count"])),
        major_iterations=int(cast(int, payload["upper_solve_count"])),
        best_tpd=float(cast(float, step2["minimum_tpd"])),
        lower_bound=None,
        upper_bound=None,
        held_gap=None if gap is None else float(cast(float, gap)),
        material_balance_max_abs=(
            None
            if certificate is None
            else max(
                float(cast(float, certificate["modified_balance_inf"])),
                float(cast(float, certificate["ordinary_balance_inf"])),
            )
        ),
        pressure_stationarity_max_relative=(
            None
            if certificate is None
            else float(cast(float, certificate["pressure_residual_inf"]))
        ),
        kkt_stationarity_max_abs=(
            None if certificate is None else float(cast(float, certificate["kkt_residual_inf"]))
        ),
        chemical_potential_max_relative=None,
        confirmation_succeeded=bool(step10 and step10["status"] == "complete"),
        confirmation_max_difference=None,
        search_profiles=("HELD2.0 paper Steps 1-10",),
        globality_certificate=str(payload["globality_certificate"]),
        failure_reason="" if failure is None else str(failure),
        performance=_held2_performance(payload),
    )


def _held2_phase_state(
    component_count: int,
    payload: Mapping[str, object],
) -> PhaseState:
    fractions = _vector(payload["mole_fractions"], component_count, "HELD2 phase composition")
    phase_fraction = float(cast(float, payload["phase_fraction"]))
    molar_volume = float(cast(float, payload["molar_volume_m3_mol"]))
    if not math.isfinite(phase_fraction) or phase_fraction <= 0.0:
        raise ValueError("native HELD2 phase fraction must be positive and finite")
    if not math.isfinite(molar_volume) or molar_volume <= 0.0:
        raise ValueError("native HELD2 phase volume must be positive and finite")
    pressure_pa = float(cast(float, payload["pressure_pa"]))
    if not math.isfinite(pressure_pa):
        raise ValueError("native HELD2 phase pressure must be finite")
    return PhaseState(
        amount_mol=phase_fraction,
        mole_fractions=fractions,
        volume_m3=phase_fraction * molar_volume,
        molar_density_mol_m3=1.0 / molar_volume,
        pressure_pa=pressure_pa,
        chemical_potential_over_rt=_vector(
            payload["chemical_potential_over_rt"],
            component_count,
            "HELD2 chemical-potential vector",
        ),
    )


def _held2_result(
    temperature_k: float,
    pressure_pa: float,
    feed: tuple[float, ...],
    fingerprint: str,
    native: Mapping[str, object],
) -> TpFlashResult:
    diagnostics = _held2_diagnostics(native)
    if diagnostics.globality_certificate != "not_guaranteed":
        raise ValueError("native HELD2 result has an invalid globality certificate")
    if diagnostics.outcome not in {
        "one_phase_no_negative_witness_detected",
        "physical_equilibrium_accepted",
    }:
        reason = diagnostics.failure_reason or "HELD2 did not return a certified equilibrium"
        raise FlashError(reason, diagnostics)
    if str(native["parameter_fingerprint"]) != fingerprint:
        raise ValueError("native HELD2 result has the wrong Provider fingerprint")

    phase_payloads = cast(Sequence[Mapping[str, object]], native["phases"])
    phases = tuple(
        _held2_phase_state(
            len(feed),
            payload,
        )
        for payload in phase_payloads
    )
    phase_fractions = tuple(phase.amount_mol for phase in phases)
    expected_phase_count = 1 if diagnostics.outcome.startswith("one_phase") else 2
    if len(phases) != expected_phase_count or not math.isclose(
        math.fsum(phase_fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-10
    ):
        raise ValueError("native HELD2 phase result is incomplete")
    total_gibbs = float(cast(float, native["total_free_energy_over_rt"]))
    if not math.isfinite(total_gibbs):
        raise ValueError("native HELD2 objective must be finite")
    return TpFlashResult(
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        overall_mole_fractions=feed,
        phases=tuple(phases),
        phase_fractions=phase_fractions,
        total_free_energy_over_rt=total_gibbs,
        parameter_fingerprint=fingerprint,
        diagnostics=diagnostics,
    )


def tp_flash(
    model: epcsaft.Mixture,
    temperature: Quantity[Any],
    pressure: Quantity[Any],
    overall_mole_fractions: Sequence[float],
    *,
    trace: bool = False,
) -> TpFlashResult:
    """Run the bounded HELD or strong-electrolyte HELD2 controller."""

    try:
        if not isinstance(model, epcsaft.Mixture):
            raise TypeError("tp_flash requires an epcsaft.Mixture")
        temperature_k = _quantity(temperature, "kelvin", "temperature", "tp_flash")
        pressure_pa = _quantity(pressure, "pascal", "pressure", "tp_flash")
        component_count = len(model.component_ids)
        feed = _tp_flash_feed(overall_mole_fractions, component_count)
        if model.parameter_fingerprint == _FLASH_FINGERPRINT and not (
            _FLASH_TEMPERATURE_DOMAIN_K[0] <= temperature_k <= _FLASH_TEMPERATURE_DOMAIN_K[1]
        ):
            raise ValueError("temperature is outside the audited May 2015 source domain")
        if model.parameter_fingerprint == _FLASH_FINGERPRINT and not (
            _FLASH_PRESSURE_DOMAIN_PA[0] <= pressure_pa <= _FLASH_PRESSURE_DOMAIN_PA[1]
        ):
            raise ValueError("pressure is outside the audited May 2015 source domain")
        if model.parameter_fingerprint == _FLASH_FINGERPRINT and not (
            _FLASH_METHANE_FEED_DOMAIN[0] <= feed[0] <= _FLASH_METHANE_FEED_DOMAIN[1]
        ):
            raise ValueError("composition is outside the audited May 2015 source domain")
        capsule = epcsaft.native_sdk(model)
    except (TypeError, ValueError) as error:
        diagnostics = _failed_held_diagnostics("invalid_input", "input_rejected", str(error))
        raise FlashError(str(error), diagnostics) from error

    try:
        native = cast(
            Mapping[str, object],
            _equilibrium._solve_tp_flash(
                capsule,
                temperature_k,
                pressure_pa,
                feed,
                model.parameter_fingerprint,
                trace=trace,
            ),
        )
    except (RuntimeError, TypeError, ValueError) as error:
        diagnostics = _failed_held_diagnostics("error", "native_exception", str(error))
        raise FlashError(str(error), diagnostics) from error

    if native.get("controller") == "perdomo_held2_paper_steps_1_to_10_v1":
        try:
            return _held2_result(
                temperature_k,
                pressure_pa,
                feed,
                model.parameter_fingerprint,
                native,
            )
        except FlashError:
            raise
        except (KeyError, TypeError, ValueError) as error:
            diagnostics = _failed_held_diagnostics("error", "payload_error", str(error))
            raise FlashError(
                "native HELD2 payload does not match the typed contract", diagnostics
            ) from error

    try:
        diagnostics = _held_diagnostics(native)
        if diagnostics.globality_certificate != "not_guaranteed":
            raise ValueError("native HELD result has an invalid globality certificate")
        if diagnostics.outcome not in {"one_phase", "accepted"}:
            reason = diagnostics.failure_reason or "HELD search did not return an accepted result"
            raise FlashError(reason, diagnostics)
        overall = _vector(native["overall_mole_fractions"], len(feed), "overall composition")
        phase_payloads = cast(Sequence[Mapping[str, object]], native["phases"])
        phases = tuple(_phase(payload) for payload in phase_payloads)
        phase_fractions = _float_tuple(native["phase_fractions"])
        if len(phases) != len(phase_fractions):
            raise ValueError("native HELD phase count does not match its outcome")
        if not math.isclose(sum(phase_fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-10):
            raise ValueError("native HELD phase fractions do not sum to one")
        if str(native["parameter_fingerprint"]) != model.parameter_fingerprint:
            raise ValueError("native HELD result has the wrong provider fingerprint")
        return TpFlashResult(
            temperature_k=float(cast(float, native["temperature_k"])),
            pressure_pa=float(cast(float, native["pressure_pa"])),
            overall_mole_fractions=overall,
            phases=phases,
            phase_fractions=phase_fractions,
            total_free_energy_over_rt=float(cast(float, native["total_free_energy_over_rt"])),
            parameter_fingerprint=str(native["parameter_fingerprint"]),
            diagnostics=diagnostics,
        )
    except FlashError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        diagnostics = _failed_held_diagnostics("error", "payload_error", str(error))
        raise FlashError(
            "native HELD payload does not match the typed contract", diagnostics
        ) from error


def _empty_chemical_search(status: str = "not_evaluated") -> ChemicalEquilibriumSearch:
    return ChemicalEquilibriumSearch(
        status=status,
        continuation_status="not_used",
        continuation_blocker="",
        continuation_initial_model_fingerprint="",
        continuation_accepted_lambda=None,
        continuation_attempt_count=0,
        primary_budget=0,
        primary_attempt_count=0,
        generated_start_count=0,
        budget_truncated_start_count=0,
        duplicate_start_count=0,
        infeasible_start_count=0,
        evaluated_start_count=0,
        domain_rejected_start_count=0,
        construction_rejected_start_count=0,
        attempts=(),
        basins=(),
        budget_prefixes=(),
        selected_basin_ordinal=None,
        selected_objective=None,
        selection_label="lowest_observed_certified_local_value",
    )


_CHEMICAL_SEARCH_STATUSES = {
    "certified_local_minimum",
    "saddle_observed",
    "second_order_inconclusive",
    "boundary_unadjudicated",
    "domain_rejected",
    "no_feasible_start_found",
    "search_exhausted_no_certified_candidate",
    "infeasible_certified",
    "not_evaluated",
}
_CHEMICAL_TERMINAL_STATUSES = _CHEMICAL_SEARCH_STATUSES | {"solver_failed"}


def _validate_chemical_search(search: ChemicalEquilibriumSearch) -> None:
    if search.status not in _CHEMICAL_SEARCH_STATUSES:
        raise ValueError("native chemical search has an invalid status")
    if search.continuation_status not in {
        "not_used",
        "completed",
        "initial_model_not_certified",
        "endpoint_domain_incompatible",
        "endpoint_derivative_incompatible",
        "no_strictly_feasible_continuation_start",
        "step_refinement_exhausted",
        "attempt_budget_exhausted",
    }:
        raise ValueError("native chemical search has an invalid continuation status")
    if search.selection_label != "lowest_observed_certified_local_value":
        raise ValueError("native chemical search has an invalid selection label")
    if (
        search.continuation_attempt_count < 0
        or (
            search.continuation_accepted_lambda is not None
            and not 0.0 <= search.continuation_accepted_lambda <= 1.0
        )
        or (
            search.continuation_status == "not_used"
            and (
                search.continuation_attempt_count != 0
                or search.continuation_accepted_lambda is not None
                or search.continuation_blocker
                or search.continuation_initial_model_fingerprint
            )
        )
        or (
            search.continuation_status != "not_used"
            and not search.continuation_initial_model_fingerprint
        )
        or (
            search.continuation_status == "completed" and search.continuation_accepted_lambda != 1.0
        )
    ):
        raise ValueError("native chemical continuation evidence is inconsistent")
    if not 0 <= search.primary_attempt_count <= search.primary_budget <= 25:
        raise ValueError("native chemical search has an invalid primary budget")
    if (
        (search.status != "not_evaluated" and search.primary_attempt_count == 0)
        or search.generated_start_count
        != search.primary_attempt_count
        + search.budget_truncated_start_count
        + search.duplicate_start_count
        + search.infeasible_start_count
        or search.budget_truncated_start_count < 0
        or not 0 <= search.evaluated_start_count <= search.primary_attempt_count
        or not 0
        <= search.domain_rejected_start_count
        <= search.primary_attempt_count - search.evaluated_start_count
        or not 0 <= search.construction_rejected_start_count
        or search.evaluated_start_count
        + search.domain_rejected_start_count
        + search.construction_rejected_start_count
        != search.primary_attempt_count
    ):
        raise ValueError("native chemical search start accounting is inconsistent")
    if tuple(attempt.ordinal for attempt in search.attempts) != tuple(range(len(search.attempts))):
        raise ValueError("native chemical search attempt ordinals are inconsistent")
    if len(search.attempts) < search.primary_attempt_count or any(
        attempt.kind != "primary"
        or attempt.primary_ordinal != attempt.ordinal
        or attempt.parent_ordinal is not None
        for attempt in search.attempts[: search.primary_attempt_count]
    ):
        raise ValueError("native chemical search primary lineage is inconsistent")
    continuation_parents = tuple(
        attempt for attempt in search.attempts if attempt.kind == "continuation"
    )
    continuation_preflight = (
        search.continuation_status
        in {"endpoint_domain_incompatible", "endpoint_derivative_incompatible"}
        and search.continuation_attempt_count == 0
        and not search.attempts
    )
    if (search.continuation_status == "not_used" and continuation_parents) or (
        search.continuation_status != "not_used"
        and not continuation_preflight
        and (
            len(continuation_parents) != search.continuation_attempt_count + 1
            or not continuation_parents
            or continuation_parents[0].continuation_status != "initial_model"
        )
    ):
        raise ValueError("native chemical continuation attempt accounting is inconsistent")
    for attempt in search.attempts:
        if attempt.kind not in {
            "primary",
            "recovery",
            "continuation",
            "continuation_recovery",
        }:
            raise ValueError("native chemical search has an invalid attempt kind")
        if (
            attempt.kind not in {"continuation", "continuation_recovery"}
            and not 0 <= attempt.primary_ordinal < search.primary_attempt_count
        ):
            raise ValueError("native chemical search primary ordinal is out of range")
        if attempt.kind == "recovery" and (
            attempt.parent_ordinal != attempt.primary_ordinal
            or attempt.parent_ordinal >= search.primary_attempt_count
        ):
            raise ValueError("native chemical search recovery lineage is inconsistent")
        if attempt.kind == "continuation_recovery" and (
            attempt.parent_ordinal is None
            or attempt.parent_ordinal >= attempt.ordinal
            or search.attempts[attempt.parent_ordinal].kind != "continuation"
        ):
            raise ValueError("native chemical continuation recovery lineage is inconsistent")
        if attempt.start_construction_status not in {
            "accepted",
            "rejected",
            "not_evaluated",
        }:
            raise ValueError("native chemical search has an invalid start status")
        if attempt.retraction_status not in {
            "not_needed",
            "passed",
            "failed",
            "balance_preserved_by_reaction_extent",
        }:
            raise ValueError("native chemical search has an invalid retraction status")
        if attempt.continuation_status not in {
            "not_used",
            "initial_model",
            "step_accepted",
            "step_rejected",
        }:
            raise ValueError("native chemical attempt has an invalid continuation status")
        if attempt.provider_domain_status not in {
            "passed",
            "failed",
            "not_applicable",
            "not_adjudicated",
        }:
            raise ValueError("native chemical attempt has an invalid Provider-domain status")
        if attempt.terminal_status not in _CHEMICAL_TERMINAL_STATUSES:
            raise ValueError("native chemical search has an invalid terminal status")
        if attempt.local_minimum_status not in {
            "passed",
            "saddle_observed",
            "second_order_inconclusive",
            "not_adjudicated",
        }:
            raise ValueError("native chemical attempt has an invalid curvature status")
        if attempt.trace_status not in {
            "interior",
            "at_or_below_floor",
            "not_adjudicated",
        }:
            raise ValueError("native chemical attempt has an invalid trace status")
        optional_values = (
            attempt.volume_m3,
            attempt.objective,
            attempt.balance_inf_norm,
            attempt.charge_inf_norm,
            attempt.pressure_relative_residual,
            attempt.reaction_affinity_inf_norm,
            attempt.kkt_stationarity_inf_norm,
            attempt.complementarity_inf_norm,
            attempt.condition_number_inf,
        )
        if any(value is not None and not math.isfinite(value) for value in optional_values):
            raise ValueError("native chemical attempt contains a nonfinite value")
        if any(not math.isfinite(value) or value < 0.0 for value in attempt.amounts_mol):
            raise ValueError("native chemical attempt contains invalid amounts")
        if attempt.volume_m3 is not None and attempt.volume_m3 <= 0.0:
            raise ValueError("native chemical attempt contains an invalid volume")
        if (
            attempt.kkt_dimension < 0
            or attempt.kkt_rank < 0
            or attempt.kkt_rank > attempt.kkt_dimension
            or attempt.recovery_seed_count < 0
            or attempt.recovery_seed_count < attempt.recovery_solve_count
            or attempt.recovery_solve_count < 0
        ):
            raise ValueError("native chemical attempt counters are inconsistent")
        if attempt.basin_ordinal is not None and not (
            0 <= attempt.basin_ordinal < len(search.basins)
        ):
            raise ValueError("native chemical attempt basin reference is out of range")
    if tuple(basin.ordinal for basin in search.basins) != tuple(range(len(search.basins))):
        raise ValueError("native chemical basin ordinals are inconsistent")
    for basin in search.basins:
        if not 0 <= basin.representative_attempt_ordinal < len(search.attempts):
            raise ValueError("native chemical basin representative is out of range")
        representative = search.attempts[basin.representative_attempt_ordinal]
        if (
            representative.terminal_status != "certified_local_minimum"
            or representative.basin_ordinal != basin.ordinal
            or not basin.amounts_mol
            or any(not math.isfinite(value) or value <= 0.0 for value in basin.amounts_mol)
            or not math.isfinite(basin.volume_m3)
            or basin.volume_m3 <= 0.0
            or not math.isfinite(basin.objective)
        ):
            raise ValueError("native chemical basin evidence is inconsistent")
    if (search.selected_basin_ordinal is None) != (search.selected_objective is None):
        raise ValueError("native chemical search selection is incomplete")
    if search.selected_basin_ordinal is None:
        if search.status == "certified_local_minimum":
            raise ValueError("native chemical search omitted its certified selection")
    else:
        if not 0 <= search.selected_basin_ordinal < len(search.basins):
            raise ValueError("native chemical search selection is out of range")
        selected = search.basins[search.selected_basin_ordinal]
        if (
            search.status != "certified_local_minimum"
            or search.selected_objective != selected.objective
        ):
            raise ValueError("native chemical search selection is inconsistent")
    expected_prefix_budgets = tuple(
        dict.fromkeys(
            min(requested, search.primary_attempt_count)
            for requested in (1, 5, 13, 25)
            if search.primary_attempt_count > 0
        )
    )
    if tuple(prefix.primary_budget for prefix in search.budget_prefixes) != (
        expected_prefix_budgets
    ):
        raise ValueError("native chemical search budget prefixes are incomplete")
    previous_budget = 0
    previous_selection: int | None = None
    for index, prefix in enumerate(search.budget_prefixes):
        if (
            prefix.primary_budget <= previous_budget
            or prefix.primary_budget > search.primary_attempt_count
            or prefix.attempted_primary_ordinals != tuple(range(prefix.primary_budget))
            or len(set(prefix.basin_ordinals)) != len(prefix.basin_ordinals)
            or any(not 0 <= ordinal < len(search.basins) for ordinal in prefix.basin_ordinals)
            or (
                prefix.selected_basin_ordinal is not None
                and prefix.selected_basin_ordinal not in prefix.basin_ordinals
            )
            or prefix.selection_changed
            != (index > 0 and prefix.selected_basin_ordinal != previous_selection)
        ):
            raise ValueError("native chemical search budget prefix is inconsistent")
        previous_budget = prefix.primary_budget
        previous_selection = prefix.selected_basin_ordinal


def _chemical_search(native: Mapping[str, object]) -> ChemicalEquilibriumSearch:
    payload = cast(Mapping[str, object], native["search"])
    attempts = tuple(
        ChemicalEquilibriumAttempt(
            ordinal=int(cast(int, record["ordinal"])),
            primary_ordinal=int(cast(int, record["primary_ordinal"])),
            kind=str(record["kind"]),
            parent_ordinal=(
                int(cast(int, record["parent_ordinal"]))
                if record["parent_ordinal"] is not None
                else None
            ),
            start_identity=str(record["start_identity"]),
            start_construction_status=str(record["start_construction_status"]),
            retraction_status=str(record["retraction_status"]),
            continuation_status=str(record["continuation_status"]),
            provider_domain_status=str(record["provider_domain_status"]),
            solver_status=str(record["solver_status"]),
            callback_error=str(record["callback_error"]),
            terminal_status=str(record["terminal_status"]),
            amounts_mol=tuple(float(value) for value in cast(Sequence[float], record["amounts"])),
            volume_m3=_optional_float(record, "volume_m3"),
            objective=_optional_float(record, "objective"),
            balance_inf_norm=_optional_float(record, "balance_inf_norm"),
            charge_inf_norm=_optional_float(record, "charge_inf_norm"),
            pressure_relative_residual=_optional_float(record, "pressure_relative_residual"),
            reaction_affinity_inf_norm=_optional_float(record, "reaction_affinity_inf_norm"),
            kkt_stationarity_inf_norm=_optional_float(record, "kkt_stationarity_inf_norm"),
            complementarity_inf_norm=_optional_float(record, "complementarity_inf_norm"),
            kkt_dimension=int(cast(int, record["kkt_dimension"])),
            kkt_rank=int(cast(int, record["kkt_rank"])),
            condition_number_inf=_optional_float(record, "condition_number_inf"),
            local_minimum_status=str(record["local_minimum_status"]),
            trace_status=str(record["trace_status"]),
            basin_ordinal=(
                int(cast(int, record["basin_ordinal"]))
                if record["basin_ordinal"] is not None
                else None
            ),
            recovery_seed_count=int(cast(int, record["recovery_seed_count"])),
            recovery_solve_count=int(cast(int, record["recovery_solve_count"])),
        )
        for record in cast(Sequence[Mapping[str, object]], payload["attempts"])
    )
    basins = tuple(
        ChemicalEquilibriumBasin(
            ordinal=int(cast(int, record["ordinal"])),
            representative_attempt_ordinal=int(cast(int, record["representative_attempt_ordinal"])),
            amounts_mol=tuple(float(value) for value in cast(Sequence[float], record["amounts"])),
            volume_m3=float(cast(float, record["volume_m3"])),
            objective=float(cast(float, record["objective"])),
        )
        for record in cast(Sequence[Mapping[str, object]], payload["basins"])
    )
    prefixes = tuple(
        ChemicalEquilibriumBudgetPrefix(
            primary_budget=int(cast(int, record["primary_budget"])),
            attempted_primary_ordinals=tuple(
                int(value) for value in cast(Sequence[int], record["attempted_primary_ordinals"])
            ),
            basin_ordinals=tuple(
                int(value) for value in cast(Sequence[int], record["basin_ordinals"])
            ),
            selected_basin_ordinal=(
                int(cast(int, record["selected_basin_ordinal"]))
                if record["selected_basin_ordinal"] is not None
                else None
            ),
            selection_changed=bool(record["selection_changed"]),
        )
        for record in cast(Sequence[Mapping[str, object]], payload["budget_prefixes"])
    )
    search = ChemicalEquilibriumSearch(
        status=str(payload["status"]),
        continuation_status=str(payload["continuation_status"]),
        continuation_blocker=str(payload["continuation_blocker"]),
        continuation_initial_model_fingerprint=str(
            payload["continuation_initial_model_fingerprint"]
        ),
        continuation_accepted_lambda=_optional_float(payload, "continuation_accepted_lambda"),
        continuation_attempt_count=int(cast(int, payload["continuation_attempt_count"])),
        primary_budget=int(cast(int, payload["primary_budget"])),
        primary_attempt_count=int(cast(int, payload["primary_attempt_count"])),
        generated_start_count=int(cast(int, payload["generated_start_count"])),
        budget_truncated_start_count=int(cast(int, payload["budget_truncated_start_count"])),
        duplicate_start_count=int(cast(int, payload["duplicate_start_count"])),
        infeasible_start_count=int(cast(int, payload["infeasible_start_count"])),
        evaluated_start_count=int(cast(int, payload["evaluated_start_count"])),
        domain_rejected_start_count=int(cast(int, payload["domain_rejected_start_count"])),
        construction_rejected_start_count=int(
            cast(int, payload["construction_rejected_start_count"])
        ),
        attempts=attempts,
        basins=basins,
        budget_prefixes=prefixes,
        selected_basin_ordinal=(
            int(cast(int, payload["selected_basin_ordinal"]))
            if payload["selected_basin_ordinal"] is not None
            else None
        ),
        selected_objective=_optional_float(payload, "selected_objective"),
        selection_label=str(payload["selection_label"]),
    )
    _validate_chemical_search(search)
    return search


def _failed_chemical_diagnostics(
    status: str, failure_reason: str
) -> ChemicalEquilibriumDiagnostics:
    return ChemicalEquilibriumDiagnostics(
        solver_status=status,
        callback_error="",
        failure_kind="input_or_native_error",
        failure_reason=failure_reason,
        chemical_certification_level="NOT_EVALUATED",
        boundary_status="not_adjudicated",
        structural_zero_species_indices=(),
        numerical_status="not_adjudicated",
        physical_status="not_adjudicated",
        predictive_status="not_adjudicated",
        provider_domain_status="not_adjudicated",
        local_minimum_status="not_adjudicated",
        negative_curvature_recovery_status="not_needed",
        negative_curvature_recovery_attempts=0,
        negative_curvature_recovery_selected_sign=0,
        trace_status="not_adjudicated",
        globality_status="not_adjudicated",
        search=_empty_chemical_search(),
        amounts_mol=(),
        volume_m3=None,
        balance_inf_norm=None,
        charge_inf_norm=None,
        pressure_relative_residual=None,
        reaction_affinity_inf_norm=None,
        reaction_affinity_residuals=(),
        packing_fraction=None,
        packing_fraction_min=None,
        packing_fraction_max=None,
        total_ion_fraction=None,
        total_ion_fraction_max=None,
        minimum_amount_mol=None,
        trace_floor_mol=None,
        kkt_stationarity_inf_norm=None,
        physical_stationarity_residuals=(),
        physical_equality_multipliers=(),
        physical_constraint_jacobian=(),
        physical_constraint_shape=(0, 0),
        physical_lagrangian_gradient=(),
        physical_to_chart_jacobian=(),
        physical_to_chart_shape=(0, 0),
        chart_physical_pullback_residual_inf_norm=None,
        complementarity_inf_norm=None,
        numerical_criteria=(),
        physical_criteria=(),
        first_failed_numerical_criterion=None,
        first_failed_physical_criterion=None,
        kkt_dimension=0,
        kkt_rank=0,
        condition_number_inf=None,
        active_lower_bounds=(),
        active_upper_bounds=(),
        active_constraint_bounds=(),
        active_trace_species=(),
        chart_topology="not_adjudicated",
        reduced_hessian_status="not_adjudicated",
        reduced_hessian=(),
        reduced_hessian_nullspace_basis=(),
        reduced_hessian_nullspace_shape=(0, 0),
        reduced_hessian_eigenvalues=(),
        reduced_hessian_spectrum_status="not_adjudicated",
        reduced_hessian_raw_inertia=(0, 0, 0),
        reduced_hessian_inertia=(0, 0, 0),
        reduced_hessian_scale=None,
        reduced_hessian_eigenvalue_tolerance=None,
        objective_gradient=(),
        constraint_values=(),
        constraint_jacobian=(),
        lagrangian_gradient=(),
        equality_multipliers=(),
        chart_stationarity_inf_norm=None,
        lagrangian_hessian=(),
        covariant_lagrangian_hessian=(),
        derivative_check_step=None,
        objective_gradient_check_relative_error=None,
        constraint_jacobian_check_relative_error=None,
        lagrangian_hessian_check_relative_error=None,
        derivative_check_worst_entry="",
        derivative_check_worst_relative_error=None,
        derivative_check_worst_analytic_value=None,
        derivative_check_worst_finite_difference_value=None,
        derivative_check_worst_step=None,
        derivative_coordinate_order=(),
        derivative_objective_basis="dimensionless_fixed_TP_A_plus_PV_plus_reference_over_RT",
        derivative_constraint_basis="ordered_per_row_in_derivative_constraint_order",
        derivative_constraint_order=(),
        objective_gradient_shape=(0,),
        constraint_jacobian_shape=(0, 0),
        lagrangian_hessian_shape=(0, 0),
        kkt_root_jacobian=(),
        kkt_root_shape=(0, 0),
        kkt_root_status="not_evaluated",
        kkt_root_jacobian_check_relative_error=None,
        reference_representation_residual_inf_norm=None,
        reference_convergence_error=None,
        reference_derivative_availability=None,
    )


def _failed_continuation_diagnostics(
    failure_reason: str,
    initial_model_fingerprint: str,
    *,
    failure_kind: str,
    continuation_status: str,
) -> ChemicalEquilibriumDiagnostics:
    search = replace(
        _empty_chemical_search(),
        continuation_status=continuation_status,
        continuation_blocker=failure_reason,
        continuation_initial_model_fingerprint=initial_model_fingerprint,
    )
    return replace(
        _failed_chemical_diagnostics("continuation_preflight_failed", failure_reason),
        failure_kind=failure_kind,
        failure_reason=f"{failure_kind}: {failure_reason}",
        provider_domain_status=(
            "failed" if failure_kind == "physical_domain_failure" else "not_adjudicated"
        ),
        search=search,
    )


def _classified_continuation_failure(
    failure_reason: str, initial_model_fingerprint: str
) -> ChemicalEquilibriumDiagnostics | None:
    """Type only native continuation failures with an unambiguous contract meaning."""
    classifications = (
        (
            "temperature is outside a Provider continuation endpoint domain",
            "physical_domain_failure",
            "endpoint_domain_incompatible",
        ),
        (
            "Provider continuation packing domains do not overlap",
            "physical_domain_failure",
            "endpoint_domain_incompatible",
        ),
        (
            "Provider continuation endpoint derivative bases are incompatible",
            "derivative_inconsistency",
            "endpoint_derivative_incompatible",
        ),
        (
            "Provider continuation inverse-packing dimensions changed",
            "derivative_inconsistency",
            "endpoint_derivative_incompatible",
        ),
        (
            "Provider capsule is missing the reacting-phase SDK contract",
            "unsupported_derivative_capability",
            "endpoint_derivative_incompatible",
        ),
        (
            "Provider capsule is missing the reacting-phase callbacks",
            "unsupported_derivative_capability",
            "endpoint_derivative_incompatible",
        ),
    )
    for native_reason, failure_kind, status in classifications:
        if failure_reason == native_reason:
            return _failed_continuation_diagnostics(
                failure_reason,
                initial_model_fingerprint,
                failure_kind=failure_kind,
                continuation_status=status,
            )
    return None


def _classified_source_reference_failure(
    failure_reason: str,
) -> ChemicalEquilibriumDiagnostics | None:
    classifications = (
        (
            "unsupported-source-reference-transfer",
            "unsupported_source_reference_transfer",
        ),
        ("invalid-source-composition", "invalid_source_reference"),
        ("nonfinite-source-reference-state", "invalid_source_reference"),
        ("component-order-mismatch", "component_order_mismatch"),
        ("source-reference component order", "component_order_mismatch"),
        ("missing-parameter-or-topology", "missing_parameter_or_topology"),
        ("provider-domain-rejection", "provider_domain_rejection"),
        (
            "reference-limit-convergence-failure",
            "reference_limit_convergence_failure",
        ),
        ("derivative-unavailable", "source_reference_derivative_unavailable"),
        (
            "source-reference derivative is unavailable",
            "source_reference_derivative_unavailable",
        ),
    )
    for marker, failure_kind in classifications:
        if marker in failure_reason:
            return replace(
                _failed_chemical_diagnostics(
                    "source_reference_preflight_failed", failure_reason
                ),
                failure_kind=failure_kind,
                provider_domain_status=(
                    "failed"
                    if failure_kind == "provider_domain_rejection"
                    else "not_adjudicated"
                ),
            )
    if "source-reference" in failure_reason or "source " in failure_reason:
        return replace(
            _failed_chemical_diagnostics(
                "source_reference_preflight_failed", failure_reason
            ),
            failure_kind="source_reference_receipt_mismatch",
        )
    return None


def _chemical_diagnostics(
    native: Mapping[str, object],
) -> ChemicalEquilibriumDiagnostics:
    def criteria(name: str) -> tuple[ChemicalEquilibriumCriterion, ...]:
        def passes(record: ChemicalEquilibriumCriterion) -> bool:
            if not math.isfinite(record.value):
                return False
            if record.comparison == "<=":
                return record.value <= record.limit
            if record.comparison == ">":
                return record.value > record.limit
            if record.comparison == ">=":
                return record.value >= record.limit
            raise ValueError(f"native {name} comparison is invalid")

        records = tuple(
            ChemicalEquilibriumCriterion(
                name=str(cast(Mapping[str, object], record)["name"]),
                status=str(cast(Mapping[str, object], record)["status"]),
                value=float(cast(float, cast(Mapping[str, object], record)["value"])),
                limit=float(cast(float, cast(Mapping[str, object], record)["limit"])),
                comparison=str(cast(Mapping[str, object], record)["comparison"]),
            )
            for record in cast(Sequence[object], native[name])
        )
        if any(record.status not in {"passed", "failed"} for record in records):
            raise ValueError(f"native {name} status is invalid")
        if any(
            record.comparison not in {"<=", ">", ">="}
            or not math.isfinite(record.limit)
            or record.limit < 0.0
            or (record.status == "passed") != passes(record)
            for record in records
        ):
            raise ValueError(f"native {name} evidence is inconsistent")
        return records

    numerical_criteria = criteria("numerical_criteria")
    physical_criteria = criteria("physical_criteria")
    first_failed_numerical = next(
        (record.name for record in numerical_criteria if record.status == "failed"), None
    )
    first_failed_physical = next(
        (record.name for record in physical_criteria if record.status == "failed"), None
    )
    reduced_inertia = tuple(
        int(value) for value in cast(Sequence[int], native["reduced_hessian_inertia"])
    )
    if len(reduced_inertia) != 3 or any(value < 0 for value in reduced_inertia):
        raise ValueError("native reduced-Hessian inertia is invalid")
    reduced_raw_inertia = tuple(
        int(value) for value in cast(Sequence[int], native["reduced_hessian_raw_inertia"])
    )
    if len(reduced_raw_inertia) != 3 or any(value < 0 for value in reduced_raw_inertia):
        raise ValueError("native raw reduced-Hessian inertia is invalid")
    failure_kind = str(native["failure_kind"])
    if failure_kind not in {
        "none",
        "genuine_saddle",
        "derivative_inconsistency",
        "rank_deficiency",
        "ill_conditioning",
        "physical_domain_failure",
        "exhausted_multistart_search",
        "unsupported_derivative_capability",
    }:
        raise ValueError("native chemical failure kind is invalid")
    if bool(native["accepted"]) != (failure_kind == "none"):
        raise ValueError("native chemical acceptance and failure kind disagree")
    accepted = bool(native["accepted"])
    certification = str(native["chemical_certification_level"])
    numerical_status = str(native["numerical_status"])
    physical_status = str(native["physical_status"])
    local_minimum_status = str(native["local_minimum_status"])
    search = _chemical_search(native)
    if accepted != (certification == "LOCAL_EQUILIBRIUM") or (
        accepted
        and (
            numerical_status != "passed"
            or physical_status != "passed"
            or local_minimum_status != "passed"
            or search.status != "certified_local_minimum"
            or first_failed_numerical is not None
            or first_failed_physical is not None
        )
    ):
        raise ValueError("native chemical acceptance and certification evidence disagree")
    kkt_dimension = int(cast(int, native["kkt_dimension"]))
    kkt_rank = int(cast(int, native["kkt_rank"]))
    if kkt_dimension < 0 or not 0 <= kkt_rank <= kkt_dimension:
        raise ValueError("native chemical KKT rank evidence is invalid")
    objective_gradient = _float_tuple(native["objective_gradient"])
    reaction_affinity_residuals = _float_tuple(native["reaction_affinity_residuals"])
    physical_stationarity_residuals = _float_tuple(native["physical_stationarity_residuals"])
    physical_equality_multipliers = _float_tuple(native["physical_equality_multipliers"])
    physical_constraint_jacobian = _float_tuple(native["physical_constraint_jacobian"])
    physical_constraint_shape = tuple(
        int(value) for value in cast(Sequence[int], native["physical_constraint_shape"])
    )
    constraint_values = _float_tuple(native["constraint_values"])
    constraint_jacobian = _float_tuple(native["constraint_jacobian"])
    lagrangian_gradient = _float_tuple(native["lagrangian_gradient"])
    equality_multipliers = _float_tuple(native["equality_multipliers"])
    lagrangian_hessian = _float_tuple(native["lagrangian_hessian"])
    covariant_lagrangian_hessian = _float_tuple(native["covariant_lagrangian_hessian"])
    reduced_hessian = _float_tuple(native["reduced_hessian"])
    reduced_hessian_eigenvalues = _float_tuple(native["reduced_hessian_eigenvalues"])
    reduced_hessian_nullspace_basis = _float_tuple(native["reduced_hessian_nullspace_basis"])
    reduced_hessian_nullspace_shape = tuple(
        int(value) for value in cast(Sequence[int], native["reduced_hessian_nullspace_shape"])
    )
    spectrum_status = str(native["reduced_hessian_spectrum_status"])
    derivative_dimension = len(objective_gradient)
    derivative_coordinate_order = tuple(
        str(value) for value in cast(Sequence[str], native["derivative_coordinate_order"])
    )
    objective_gradient_shape = tuple(
        int(value) for value in cast(Sequence[int], native["objective_gradient_shape"])
    )
    constraint_jacobian_shape = tuple(
        int(value) for value in cast(Sequence[int], native["constraint_jacobian_shape"])
    )
    lagrangian_hessian_shape = tuple(
        int(value) for value in cast(Sequence[int], native["lagrangian_hessian_shape"])
    )
    derivative_constraint_order = tuple(
        str(value) for value in cast(Sequence[str], native["derivative_constraint_order"])
    )
    physical_lagrangian_gradient = _float_tuple(native["physical_lagrangian_gradient"])
    physical_to_chart_jacobian = _float_tuple(native["physical_to_chart_jacobian"])
    physical_to_chart_shape = tuple(
        int(value) for value in cast(Sequence[int], native["physical_to_chart_shape"])
    )
    kkt_root_jacobian = _float_tuple(native["kkt_root_jacobian"])
    kkt_root_shape = tuple(int(value) for value in cast(Sequence[int], native["kkt_root_shape"]))
    kkt_root_status = str(native["kkt_root_status"])
    reduced_dimension = math.isqrt(len(reduced_hessian))
    tolerance = _optional_float(native, "reduced_hessian_eigenvalue_tolerance")
    if (
        len(lagrangian_gradient) != derivative_dimension
        or len(constraint_jacobian) != len(constraint_values) * derivative_dimension
        or len(lagrangian_hessian) != derivative_dimension * derivative_dimension
        or len(covariant_lagrangian_hessian) != derivative_dimension * derivative_dimension
        or len(derivative_coordinate_order) != derivative_dimension
        or objective_gradient_shape != (derivative_dimension,)
        or constraint_jacobian_shape != (len(constraint_values), derivative_dimension)
        or lagrangian_hessian_shape != (derivative_dimension, derivative_dimension)
        or len(derivative_constraint_order) != len(constraint_values)
        or len(equality_multipliers) != len(constraint_values)
        or len(reduced_hessian) != reduced_dimension * reduced_dimension
        or len(reduced_hessian_nullspace_shape) != 2
        or reduced_hessian_nullspace_shape != (reduced_dimension, derivative_dimension)
        or len(reduced_hessian_nullspace_basis) != reduced_dimension * derivative_dimension
        or len(kkt_root_shape) != 2
        or len(physical_to_chart_shape) != 2
        or len(physical_lagrangian_gradient) != physical_to_chart_shape[0]
        or len(physical_to_chart_jacobian)
        != physical_to_chart_shape[0] * physical_to_chart_shape[1]
        or len(kkt_root_jacobian) != kkt_root_shape[0] * kkt_root_shape[1]
        or len(physical_constraint_shape) != 2
        or len(physical_equality_multipliers) != physical_constraint_shape[0]
        or len(physical_stationarity_residuals) != physical_constraint_shape[1]
        or len(physical_constraint_jacobian)
        != physical_constraint_shape[0] * physical_constraint_shape[1]
    ):
        raise ValueError("native chemical derivative evidence has inconsistent dimensions")
    if kkt_root_status not in {
        "not_evaluated",
        "interior_no_active_bounds",
        "unavailable_active_bounds",
        "evaluation_failed",
    }:
        raise ValueError("native KKT root-system status is invalid")
    if kkt_root_status == "interior_no_active_bounds":
        if (
            kkt_root_shape != (kkt_dimension, kkt_dimension)
            or len(kkt_root_jacobian) != kkt_dimension * kkt_dimension
        ):
            raise ValueError("native KKT root-system evidence is inconsistent")
    elif kkt_root_shape != (0, 0) or kkt_root_jacobian:
        raise ValueError("native unavailable KKT root-system evidence is not empty")
    if spectrum_status not in {
        "not_adjudicated",
        "not_evaluated",
        "converged",
        "not_converged",
        "nonfinite",
    }:
        raise ValueError("native reduced-Hessian spectrum status is invalid")
    if reduced_dimension and spectrum_status in {"converged", "not_converged"}:
        if len(reduced_hessian_eigenvalues) != reduced_dimension or tolerance is None:
            raise ValueError("native reduced-Hessian spectrum dimensions are invalid")
        classified = (
            sum(value > tolerance for value in reduced_hessian_eigenvalues),
            sum(abs(value) <= tolerance for value in reduced_hessian_eigenvalues),
            sum(value < -tolerance for value in reduced_hessian_eigenvalues),
        )
        if reduced_raw_inertia != classified or sum(reduced_inertia) != reduced_dimension:
            raise ValueError("native reduced-Hessian spectrum and inertia disagree")
    elif reduced_dimension and spectrum_status == "nonfinite":
        if (
            reduced_hessian_eigenvalues
            or reduced_raw_inertia != (0, 0, 0)
            or reduced_inertia != (0, 0, 0)
        ):
            raise ValueError("native nonfinite reduced-Hessian evidence is inconsistent")
    elif reduced_dimension:
        raise ValueError("native reduced-Hessian spectrum was not evaluated")
    elif reduced_raw_inertia != (0, 0, 0) or reduced_inertia != (0, 0, 0):
        raise ValueError("native empty reduced-Hessian inertia is invalid")
    reported_affinity_norm = float(cast(float, native["reaction_affinity_inf_norm"]))
    reported_stationarity_norm = float(cast(float, native["kkt_stationarity_inf_norm"]))
    if (
        reaction_affinity_residuals
        and not math.isclose(
            max(abs(value) for value in reaction_affinity_residuals),
            reported_affinity_norm,
            rel_tol=2.0e-12,
            abs_tol=2.0e-12,
        )
    ) or (
        physical_stationarity_residuals
        and not math.isclose(
            max(abs(value) for value in physical_stationarity_residuals),
            reported_stationarity_norm,
            rel_tol=2.0e-12,
            abs_tol=2.0e-12,
        )
    ):
        raise ValueError("native chemical residual vectors and norms disagree")
    sensitivities = cast(Mapping[str, object], native["sensitivities"])
    return ChemicalEquilibriumDiagnostics(
        solver_status=str(native["solver_status"]),
        callback_error=str(native["callback_error"]),
        failure_kind=failure_kind,
        failure_reason=str(native["failure_reason"]),
        chemical_certification_level=certification,
        boundary_status=str(native["boundary_status"]),
        structural_zero_species_indices=tuple(
            int(index) for index in cast(Sequence[int], native["structural_zero_species_indices"])
        ),
        numerical_status=numerical_status,
        physical_status=physical_status,
        predictive_status="not_adjudicated",
        provider_domain_status=str(native["provider_domain_status"]),
        local_minimum_status=local_minimum_status,
        negative_curvature_recovery_status=str(
            native.get("negative_curvature_recovery_status", "not_needed")
        ),
        negative_curvature_recovery_attempts=int(
            cast(int, native.get("negative_curvature_recovery_attempts", 0))
        ),
        negative_curvature_recovery_selected_sign=int(
            cast(int, native.get("negative_curvature_recovery_selected_sign", 0))
        ),
        trace_status=str(native["trace_status"]),
        globality_status=str(native["globality_certificate"]),
        search=search,
        amounts_mol=_float_tuple(native["amounts"]),
        volume_m3=_optional_float(native, "volume_m3"),
        balance_inf_norm=float(cast(float, native["balance_inf_norm"])),
        charge_inf_norm=float(cast(float, native["charge_inf_norm"])),
        pressure_relative_residual=float(cast(float, native["pressure_relative_residual"])),
        reaction_affinity_inf_norm=reported_affinity_norm,
        reaction_affinity_residuals=reaction_affinity_residuals,
        packing_fraction=float(cast(float, native["packing_fraction"])),
        packing_fraction_min=_optional_float(native, "packing_fraction_min"),
        packing_fraction_max=_optional_float(native, "packing_fraction_max"),
        total_ion_fraction=_optional_float(native, "total_ion_fraction"),
        total_ion_fraction_max=_optional_float(native, "total_ion_fraction_max"),
        minimum_amount_mol=_optional_float(native, "minimum_amount_mol"),
        trace_floor_mol=_optional_float(native, "trace_floor_mol"),
        kkt_stationarity_inf_norm=reported_stationarity_norm,
        physical_stationarity_residuals=physical_stationarity_residuals,
        physical_equality_multipliers=physical_equality_multipliers,
        physical_constraint_jacobian=physical_constraint_jacobian,
        physical_constraint_shape=physical_constraint_shape,
        physical_lagrangian_gradient=physical_lagrangian_gradient,
        physical_to_chart_jacobian=physical_to_chart_jacobian,
        physical_to_chart_shape=physical_to_chart_shape,
        chart_physical_pullback_residual_inf_norm=_optional_float(
            native, "chart_physical_pullback_residual_inf_norm"
        ),
        complementarity_inf_norm=float(cast(float, native["complementarity_inf_norm"])),
        numerical_criteria=numerical_criteria,
        physical_criteria=physical_criteria,
        first_failed_numerical_criterion=first_failed_numerical,
        first_failed_physical_criterion=first_failed_physical,
        kkt_dimension=kkt_dimension,
        kkt_rank=kkt_rank,
        condition_number_inf=_optional_float(native, "condition_number_inf"),
        active_lower_bounds=_index_tuple(native["active_lower_bounds"], "active lower bounds"),
        active_upper_bounds=_index_tuple(native["active_upper_bounds"], "active upper bounds"),
        active_constraint_bounds=_index_tuple(
            native["active_constraint_bounds"], "active constraint bounds"
        ),
        active_trace_species=_index_tuple(
            sensitivities["active_trace_species"], "active trace species"
        ),
        chart_topology=str(sensitivities["chart_topology"]),
        reduced_hessian_status=str(native["reduced_hessian_status"]),
        reduced_hessian=reduced_hessian,
        reduced_hessian_nullspace_basis=reduced_hessian_nullspace_basis,
        reduced_hessian_nullspace_shape=reduced_hessian_nullspace_shape,
        reduced_hessian_eigenvalues=reduced_hessian_eigenvalues,
        reduced_hessian_spectrum_status=spectrum_status,
        reduced_hessian_raw_inertia=reduced_raw_inertia,
        reduced_hessian_inertia=reduced_inertia,
        reduced_hessian_scale=_optional_float(native, "reduced_hessian_scale"),
        reduced_hessian_eigenvalue_tolerance=tolerance,
        objective_gradient=objective_gradient,
        constraint_values=constraint_values,
        constraint_jacobian=constraint_jacobian,
        lagrangian_gradient=lagrangian_gradient,
        equality_multipliers=equality_multipliers,
        chart_stationarity_inf_norm=_optional_float(native, "chart_stationarity_inf_norm"),
        lagrangian_hessian=lagrangian_hessian,
        covariant_lagrangian_hessian=covariant_lagrangian_hessian,
        derivative_check_step=_optional_float(native, "derivative_check_step"),
        objective_gradient_check_relative_error=_optional_float(
            native, "objective_gradient_check_relative_error"
        ),
        constraint_jacobian_check_relative_error=_optional_float(
            native, "constraint_jacobian_check_relative_error"
        ),
        lagrangian_hessian_check_relative_error=_optional_float(
            native, "lagrangian_hessian_check_relative_error"
        ),
        derivative_check_worst_entry=str(native["derivative_check_worst_entry"]),
        derivative_check_worst_relative_error=_optional_float(
            native, "derivative_check_worst_relative_error"
        ),
        derivative_check_worst_analytic_value=_optional_float(
            native, "derivative_check_worst_analytic_value"
        ),
        derivative_check_worst_finite_difference_value=_optional_float(
            native, "derivative_check_worst_finite_difference_value"
        ),
        derivative_check_worst_step=_optional_float(native, "derivative_check_worst_step"),
        derivative_coordinate_order=derivative_coordinate_order,
        derivative_objective_basis=str(native["derivative_objective_basis"]),
        derivative_constraint_basis=str(native["derivative_constraint_basis"]),
        derivative_constraint_order=derivative_constraint_order,
        objective_gradient_shape=objective_gradient_shape,
        constraint_jacobian_shape=constraint_jacobian_shape,
        lagrangian_hessian_shape=lagrangian_hessian_shape,
        kkt_root_jacobian=kkt_root_jacobian,
        kkt_root_shape=kkt_root_shape,
        kkt_root_status=kkt_root_status,
        kkt_root_jacobian_check_relative_error=_optional_float(
            native, "kkt_root_jacobian_check_relative_error"
        ),
        reference_representation_residual_inf_norm=_optional_float(
            native, "reference_representation_residual_inf_norm"
        )
        if "reference_representation_residual_inf_norm" in native
        else None,
        reference_convergence_error=_optional_float(native, "reference_convergence_error")
        if "reference_convergence_error" in native
        else None,
        reference_derivative_availability=int(
            cast(int, native["reference_derivative_availability"])
        )
        if "reference_derivative_availability" in native
        else None,
    )


@cache
def _distribution_identity(distribution_name: str) -> tuple[str, str]:
    distribution = metadata.distribution(distribution_name)
    record = distribution.read_text("RECORD")
    if record is None or not record.strip():
        raise ValueError(f"installed {distribution_name} distribution has no RECORD identity")
    digest = hashlib.sha256(record.encode("utf-8")).hexdigest()
    return distribution.version, f"sha256:{digest}"


def _distribution_file_identity(
    distribution_name: str,
    relative_path: str,
) -> str:
    path = metadata.distribution(distribution_name).locate_file(relative_path)
    payload = path.read_bytes()
    if not payload:
        raise ValueError(f"installed {distribution_name} file identity is unavailable")
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def _chemical_artifact_identity(
    native: Mapping[str, object],
    *,
    uses_provider: bool,
) -> ChemicalArtifactIdentity:
    equilibrium_version, equilibrium_record = _distribution_identity("epcsaft-equilibrium")
    if not uses_provider:
        return ChemicalArtifactIdentity(
            equilibrium_distribution="epcsaft-equilibrium",
            equilibrium_version=equilibrium_version,
            equilibrium_record_sha256=equilibrium_record,
            provider_distribution=None,
            provider_version=None,
            provider_record_sha256=None,
            provider_sdk_capsule_name=None,
            provider_sdk_abi_version=None,
            provider_sdk_table_size=None,
            provider_sdk_result_size=None,
            provider_sdk_mixture_result_size=None,
            provider_sdk_neutral_reference_result_size=None,
            provider_sdk_neutral_reference_derivative_result_size=None,
            provider_sdk_reacting_phase_parameter_result_size=None,
        )
    provider_version, provider_record = _distribution_identity("epcsaft")
    capsule_name = str(native["provider_sdk_capsule_name"])
    abi_version = int(cast(int, native["provider_sdk_abi_version"]))
    sizes = tuple(
        int(cast(int, native[name]))
        for name in (
            "provider_sdk_table_size",
            "provider_sdk_result_size",
            "provider_sdk_mixture_result_size",
            "provider_sdk_neutral_reference_result_size",
            "provider_sdk_neutral_reference_derivative_result_size",
            "provider_sdk_reacting_phase_parameter_result_size",
        )
    )
    if capsule_name != "epcsaft.native_sdk.v1" or abi_version != 1:
        raise ValueError("native result has an incompatible Provider SDK identity")
    if any(size <= 0 for size in sizes[:3]) or any(size < 0 for size in sizes[3:]):
        raise ValueError("native result has invalid Provider SDK structure sizes")
    return ChemicalArtifactIdentity(
        equilibrium_distribution="epcsaft-equilibrium",
        equilibrium_version=equilibrium_version,
        equilibrium_record_sha256=equilibrium_record,
        provider_distribution="epcsaft",
        provider_version=provider_version,
        provider_record_sha256=provider_record,
        provider_sdk_capsule_name=capsule_name,
        provider_sdk_abi_version=abi_version,
        provider_sdk_table_size=sizes[0],
        provider_sdk_result_size=sizes[1],
        provider_sdk_mixture_result_size=sizes[2],
        provider_sdk_neutral_reference_result_size=sizes[3],
        provider_sdk_neutral_reference_derivative_result_size=sizes[4],
        provider_sdk_reacting_phase_parameter_result_size=sizes[5],
    )


def _active_parameter_name(parameter: ChemicalEquilibriumActiveParameter) -> str:
    return (
        f"provider_parameter[{parameter.family};{parameter.identity};"
        f"{','.join(parameter.component_ids)}]"
    )


def _sensitivity_parameter(
    name: str,
    active_parameters: tuple[ChemicalEquilibriumActiveParameter, ...] = (),
) -> ChemicalEquilibriumSensitivityParameter:
    if name == "pressure_pa":
        return ChemicalEquilibriumSensitivityParameter(
            name=name,
            kind="pressure",
            source_index=None,
            input_unit="Pa",
            amount_derivative_unit="mol/Pa",
            volume_derivative_unit="m^3/Pa",
        )
    for prefix, kind, input_unit, amount_unit, volume_unit in (
        (
            "balance_total[",
            "compiled_balance_total",
            "mol",
            "mol/mol",
            "m^3/mol",
        ),
        (
            "ln_k_provider_basis[",
            "provider_basis_ln_k",
            "dimensionless",
            "mol",
            "m^3",
        ),
    ):
        if name.startswith(prefix) and name.endswith("]"):
            index_text = name[len(prefix) : -1]
            if index_text.isdecimal():
                return ChemicalEquilibriumSensitivityParameter(
                    name=name,
                    kind=kind,
                    source_index=int(index_text),
                    input_unit=input_unit,
                    amount_derivative_unit=amount_unit,
                    volume_derivative_unit=volume_unit,
                )
    for parameter in active_parameters:
        if name == _active_parameter_name(parameter):
            return ChemicalEquilibriumSensitivityParameter(
                name=name,
                kind=f"provider_{parameter.identity}_{parameter.family}",
                source_index=None,
                input_unit=parameter.unit,
                amount_derivative_unit=f"mol/{parameter.unit}",
                volume_derivative_unit=f"m^3/{parameter.unit}",
                identity=parameter.identity,
                component_ids=parameter.component_ids,
                value=parameter.value,
            )
    raise ValueError(f"native sensitivity parameter identity is unsupported: {name}")


def _index_tuple(payload: object, name: str) -> tuple[int, ...]:
    values = tuple(int(value) for value in cast(Sequence[int], payload))
    if any(value < 0 for value in values):
        raise ValueError(f"native sensitivity {name} contains a negative index")
    return values


def _chemical_sensitivity(
    native: Mapping[str, object],
    *,
    species_ids: tuple[str, ...],
    model_fingerprint: str,
    active_parameters: tuple[ChemicalEquilibriumActiveParameter, ...],
) -> ChemicalEquilibriumSensitivity:
    payload = cast(Mapping[str, object], native["sensitivities"])
    status = str(payload["status"])
    if status not in {"available", "unavailable"}:
        raise ValueError("native sensitivity status is invalid")
    failure_reason = str(payload["failure_reason"])
    parameter_names = tuple(
        str(value) for value in cast(Sequence[object], payload["parameter_order"])
    )
    parameters = tuple(_sensitivity_parameter(name, active_parameters) for name in parameter_names)
    amount_derivatives = tuple(
        _vector(row, len(species_ids), "chemical sensitivity amount row")
        for row in cast(Sequence[object], payload["amount_derivatives"])
    )
    volume_derivatives = _float_tuple(payload["volume_derivatives"])
    if len(amount_derivatives) != len(parameters) or len(volume_derivatives) != len(parameters):
        raise ValueError("native sensitivity parameter and derivative dimensions disagree")
    if status == "available":
        if failure_reason or not parameters:
            raise ValueError("available native sensitivity has incomplete status evidence")
        balance_parameters = tuple(
            parameter for parameter in parameters if parameter.kind == "compiled_balance_total"
        )
        reaction_parameters = tuple(
            parameter for parameter in parameters if parameter.kind == "provider_basis_ln_k"
        )
        expected_parameters = (
            balance_parameters
            + reaction_parameters
            + (_sensitivity_parameter("pressure_pa"),)
            + tuple(
                _sensitivity_parameter(_active_parameter_name(parameter), active_parameters)
                for parameter in active_parameters
            )
        )
        if (
            parameters != expected_parameters
            or tuple(parameter.source_index for parameter in balance_parameters)
            != tuple(range(len(balance_parameters)))
            or tuple(parameter.source_index for parameter in reaction_parameters)
            != tuple(range(len(reaction_parameters)))
        ):
            raise ValueError("available native sensitivity parameter order is incomplete")
        values = (value for row in amount_derivatives for value in row)
        if not all(math.isfinite(value) for value in (*values, *volume_derivatives)):
            raise ValueError("available native sensitivity contains a non-finite derivative")
    elif parameters or amount_derivatives or volume_derivatives or not failure_reason:
        raise ValueError("unavailable native sensitivity exposes derivative columns")
    kkt_dimension = int(cast(int, payload["kkt_dimension"]))
    kkt_rank = int(cast(int, payload["kkt_rank"]))
    condition_number = float(cast(float, payload["condition_number_inf"]))
    if kkt_dimension < 0 or not 0 <= kkt_rank <= kkt_dimension:
        raise ValueError("native sensitivity KKT rank evidence is invalid")
    if math.isnan(condition_number) or condition_number < 0.0:
        raise ValueError("native sensitivity KKT conditioning evidence is invalid")
    if status == "available" and (
        kkt_dimension == 0 or kkt_rank != kkt_dimension or not math.isfinite(condition_number)
    ):
        raise ValueError("available native sensitivity lacks a conditioned square KKT system")
    if status == "available" and len(parameters) - len(active_parameters) != kkt_dimension:
        raise ValueError("available native sensitivity omits a KKT input column")
    chart_topology = str(payload["chart_topology"])
    parameter_fingerprint = str(payload["parameter_fingerprint"])
    if not chart_topology or parameter_fingerprint != model_fingerprint:
        raise ValueError("native sensitivity topology or parameter fingerprint is invalid")
    active_trace_species = _index_tuple(payload["active_trace_species"], "active trace species")
    if any(index >= len(species_ids) for index in active_trace_species):
        raise ValueError("native sensitivity trace-species index is out of range")
    provider_parameter_status = str(payload["provider_parameter_status"])
    reference_parameter_status = str(payload["reference_parameter_status"])
    if provider_parameter_status not in {"available", "not_applicable", "unavailable"} or (
        reference_parameter_status not in {"available", "not_applicable", "unavailable"}
    ):
        raise ValueError("native sensitivity parameter-support status is invalid")
    provider_parameter_failure = str(payload["provider_parameter_failure_reason"])
    reference_parameter_failure = str(payload["reference_parameter_failure_reason"])
    if (provider_parameter_status == "unavailable") != bool(provider_parameter_failure) or (
        reference_parameter_status == "unavailable"
    ) != bool(reference_parameter_failure):
        raise ValueError("native sensitivity parameter-support evidence is inconsistent")
    active_lower_bounds = _index_tuple(payload["active_lower_bounds"], "active lower bounds")
    active_upper_bounds = _index_tuple(payload["active_upper_bounds"], "active upper bounds")
    active_constraint_bounds = _index_tuple(
        payload["active_constraint_bounds"], "active constraint bounds"
    )
    if status == "available" and (
        active_lower_bounds
        or active_upper_bounds
        or active_constraint_bounds
        or active_trace_species
    ):
        raise ValueError("available native sensitivity has active-set evidence")
    return ChemicalEquilibriumSensitivity(
        status=status,
        failure_reason=failure_reason,
        parameters=parameters,
        amount_state_order=species_ids,
        amount_derivatives=amount_derivatives,
        volume_derivatives=volume_derivatives,
        kkt_dimension=kkt_dimension,
        kkt_rank=kkt_rank,
        condition_number_inf=condition_number,
        active_lower_bounds=active_lower_bounds,
        active_upper_bounds=active_upper_bounds,
        active_constraint_bounds=active_constraint_bounds,
        active_trace_species=active_trace_species,
        chart_topology=chart_topology,
        parameter_fingerprint=parameter_fingerprint,
        provider_parameter_status=provider_parameter_status,
        provider_parameter_failure_reason=provider_parameter_failure,
        reference_parameter_status=reference_parameter_status,
        reference_parameter_failure_reason=reference_parameter_failure,
    )


def _source_reference_transfer(
    phase: ProviderPhase,
    standard_state: ChemicalStandardState,
    *,
    species_ids: tuple[str, ...],
    temperature_k: float,
    pressure_pa: float,
    sensitivities_requested: bool,
    active_parameters_requested: bool,
) -> epcsaft.SourceReferenceTransfer:
    if standard_state.source_reference_component_ids != species_ids:
        raise ValueError(
            "source-reference component order does not match the reaction system"
        )
    if tuple(phase.model.component_ids) != species_ids:
        raise ValueError(
            "source-reference component order does not match the installed Provider"
        )
    if len(standard_state.log_activity_scale_factors) != len(species_ids):
        raise ValueError("source activity-scale factors do not match the species order")
    required_derivatives = (
        ()
        if not sensitivities_requested
        else (("pressure", "parameters") if active_parameters_requested else ("pressure",))
    )
    transfer = epcsaft.source_reference_transfer(
        phase.model,
        epcsaft.SourceReferenceDeclaration(
            reference_state_id=standard_state.id,
            component_ids=standard_state.source_reference_component_ids,
            solvent_mole_fractions=(
                standard_state.source_reference_solvent_composition
            ),
            activity_convention_id=(
                standard_state.source_reference_activity_convention_id
            ),
            standard_molality_mol_per_kg=(
                standard_state.source_reference_standard_molality_mol_per_kg
            ),
            reference_pressure_pa=standard_state.reference_pressure_pa,
            required_derivatives=required_derivatives,
        ),
        T=temperature_k * epcsaft.unit_registry.kelvin,
        P=pressure_pa * epcsaft.unit_registry.pascal,
    )
    if (
        transfer.status != "ok"
        or transfer.domain_status != "admitted"
        or transfer.convergence_status != "converged"
        or transfer.component_ids != species_ids
        or transfer.reference_state_id != standard_state.id
        or transfer.activity_convention_id
        != standard_state.source_reference_activity_convention_id
        or transfer.parameter_fingerprint != phase.expected_parameter_fingerprint
        or transfer.temperature_k != temperature_k
        or transfer.pressure_pa != pressure_pa
        or transfer.source_reference_pressure_pa
        != standard_state.reference_pressure_pa
    ):
        raise ValueError("installed Provider source-reference receipt is incompatible")
    return transfer


def _source_reference_payload(
    transfer: epcsaft.SourceReferenceTransfer,
) -> dict[str, object]:
    return {
        "status": transfer.status,
        "domain_status": transfer.domain_status,
        "convergence_status": transfer.convergence_status,
        "reference_state_id": transfer.reference_state_id,
        "activity_convention_id": transfer.activity_convention_id,
        "component_ids": transfer.component_ids,
        "neutral_basis_row_count": len(transfer.neutral_basis),
        "neutral_basis": tuple(
            value for row in transfer.neutral_basis for value in row
        ),
        "log_fugacity_contractions": transfer.log_fugacity_contractions,
        "activity_scale_log_contractions": (
            transfer.activity_scale_log_contractions
        ),
        "transfer_log_contractions": transfer.transfer_log_contractions,
        "source_composition": transfer.source_composition,
        "derivative_availability": transfer.derivative_availability,
        "pressure_derivatives_per_pa": transfer.pressure_derivatives_per_pa,
        "temperature_k": transfer.temperature_k,
        "pressure_pa": transfer.pressure_pa,
        "source_reference_pressure_pa": transfer.source_reference_pressure_pa,
        "standard_molality_mol_per_kg": (
            transfer.standard_molality_mol_per_kg
        ),
        "solvent_molar_mass_kg_per_mol": (
            transfer.solvent_molar_mass_kg_per_mol
        ),
        "reference_limit_molality_mol_per_kg": (
            transfer.reference_limit_molality_mol_per_kg
        ),
        "reference_convergence_error": transfer.reference_convergence_error,
        "source_temperature_interval_k": transfer.source_temperature_interval_k,
        "source_pressure_interval_pa": transfer.source_pressure_interval_pa,
        "parameter_fingerprint": transfer.parameter_fingerprint,
        "topology_fingerprint": transfer.topology_fingerprint,
        "component_order_fingerprint": transfer.component_order_fingerprint,
        "reference_state_fingerprint": transfer.reference_state_fingerprint,
        "domain_fingerprint": transfer.domain_fingerprint,
        "artifact_fingerprint": transfer.artifact_fingerprint,
        "helmholtz_basis_id": transfer.helmholtz_basis_id,
    }


def chemical_equilibrium(
    phase: IdealGasPhase | ProviderPhase,
    temperature: Quantity[Any],
    pressure: Quantity[Any],
    problem: ChemicalEquilibriumProblem,
    *,
    sensitivity_request: ChemicalEquilibriumSensitivityRequest | None = None,
    continuation: ProviderModelContinuation | None = None,
) -> ChemicalEquilibriumResult:
    """Solve and certify one local fixed-T,P homogeneous reacting phase."""

    try:
        if not isinstance(problem, ChemicalEquilibriumProblem):
            raise TypeError("chemical_equilibrium requires a typed problem")
        if sensitivity_request is not None and not isinstance(
            sensitivity_request, ChemicalEquilibriumSensitivityRequest
        ):
            raise TypeError("sensitivity_request must be a typed sensitivity request")
        if continuation is not None and not isinstance(continuation, ProviderModelContinuation):
            raise TypeError("continuation must be a typed Provider model continuation")
        if continuation is not None and sensitivity_request is not None:
            raise ValueError("Provider model continuation does not support sensitivity requests")
        active_parameters = (
            () if sensitivity_request is None else sensitivity_request.active_parameters
        )
        if any(
            not isinstance(parameter, ChemicalEquilibriumActiveParameter)
            for parameter in active_parameters
        ):
            raise TypeError("active parameters must use the typed request record")
        active_parameter_payload = tuple(
            {
                "family": parameter.family,
                "identity": parameter.identity,
                "component_ids": parameter.component_ids,
                "value": parameter.value,
                "unit": parameter.unit,
            }
            for parameter in active_parameters
        )
        for parameter in active_parameters:
            component_shape_valid = (
                (parameter.identity == "component" and len(parameter.component_ids) == 1)
                or (
                    parameter.identity == "unordered_component_pair"
                    and len(parameter.component_ids) == 2
                    and parameter.component_ids[0] != parameter.component_ids[1]
                )
                or (parameter.identity == "model" and not parameter.component_ids)
            )
            if (
                not parameter.family
                or not parameter.identity
                or not component_shape_valid
                or any(not component_id for component_id in parameter.component_ids)
                or not math.isfinite(parameter.value)
                or not parameter.unit
            ):
                raise ValueError("active Provider parameter request is incomplete")
        temperature_k = _quantity(temperature, "kelvin", "temperature", "chemical_equilibrium")
        pressure_pa = _quantity(pressure, "pascal", "pressure", "chemical_equilibrium")
        standard_state = problem.source_standard_state
        if standard_state is not None:
            if (
                not math.isfinite(standard_state.reference_pressure_pa)
                or standard_state.reference_pressure_pa <= 0.0
            ):
                raise ValueError(
                    "source-standard-state reference pressure must be positive and finite"
                )
            if (
                not standard_state.id
                or not standard_state.activity_scale_id
                or not standard_state.source_reference_activity_convention_id
                or not math.isfinite(
                    standard_state.source_reference_standard_molality_mol_per_kg
                )
                or standard_state.source_reference_standard_molality_mol_per_kg
                <= 0.0
            ):
                raise ValueError("source-standard-state declaration is incomplete")
        if not math.isfinite(problem.strict_interior_amount_floor_mol) or (
            problem.strict_interior_amount_floor_mol <= 0.0
        ):
            raise ValueError("minimum admitted amount must be positive and finite")
        records = tuple(
            {
                "source_id": record.source_id,
                "reference_id": record.reference_id,
                "reaction_orientation": record.reaction_orientation,
                "conversion_id": record.conversion_id,
                "dimensionless": record.dimensionless,
                "temperature_k": temperature_k,
                "pressure_pa": (
                    standard_state.reference_pressure_pa
                    if standard_state is not None
                    else pressure_pa
                ),
            }
            for record in problem.equilibrium_constants
        )
        spec = {
            "species_ids": problem.species_ids,
            "charges": problem.charges,
            "molar_masses_kg_per_mol": problem.molar_masses_kg_per_mol,
            "balance_matrix": problem.balance_matrix,
            "conserved_totals": problem.conserved_totals,
            "reaction_matrix": problem.reaction_matrix,
            "feed_amounts": problem.feed_amounts_mol,
            "ln_k": tuple(record.ln_value for record in problem.equilibrium_constants),
            "equilibrium_constant_records": records,
            "temperature_k": temperature_k,
            "pressure_pa": pressure_pa,
        }
        packing_bounds: tuple[float, float] | None
        source_transfer: epcsaft.SourceReferenceTransfer | None = None
        initial_source_transfer: epcsaft.SourceReferenceTransfer | None = None
        if isinstance(phase, IdealGasPhase):
            if continuation is not None:
                raise ValueError("Provider model continuation requires a ProviderPhase")
            if active_parameters:
                raise ValueError("ideal-gas sensitivity cannot request Provider coordinates")
            if problem.source_standard_state is not None:
                raise ValueError("ideal-gas problems cannot request Provider reference transport")
            if not phase.model_fingerprint or not phase.reference_id:
                raise ValueError("ideal-gas model identity is incomplete")
            if any(
                record.reference_id != phase.reference_id
                for record in problem.equilibrium_constants
            ):
                raise ValueError("reaction constants do not match the ideal-gas reference")
            spec["provider_fingerprint"] = phase.model_fingerprint
            capsule = None
            packing_bounds = None
            thermodynamic_model = "ideal_gas"
            provider_ids = None
            provider_fingerprint = None
            model_fingerprint = phase.model_fingerprint
        elif isinstance(phase, ProviderPhase):
            if not isinstance(phase.model, epcsaft.Mixture):
                raise TypeError("ProviderPhase requires an epcsaft.Mixture")
            if (
                not phase.expected_parameter_fingerprint
                or phase.expected_parameter_fingerprint != phase.model.parameter_fingerprint
            ):
                raise ValueError("installed Provider fingerprint does not match the problem")
            packing_bounds = (
                float(phase.admissible_packing_fraction_interval[0]),
                float(phase.admissible_packing_fraction_interval[1]),
            )
            if (
                not all(math.isfinite(value) for value in packing_bounds)
                or not 0.0 < packing_bounds[0] < packing_bounds[1]
            ):
                raise ValueError("packing-fraction bounds must be finite and increasing")
            model_fingerprint = phase.expected_parameter_fingerprint
            spec["provider_fingerprint"] = model_fingerprint
            capsule = epcsaft.native_sdk(phase.model)
            thermodynamic_model = "installed_provider"
            provider_ids = tuple(phase.model.component_ids)
            provider_fingerprint = model_fingerprint
            if standard_state is not None:
                source_transfer = _source_reference_transfer(
                    phase,
                    standard_state,
                    species_ids=problem.species_ids,
                    temperature_k=temperature_k,
                    pressure_pa=pressure_pa,
                    sensitivities_requested=sensitivity_request is not None,
                    active_parameters_requested=bool(active_parameters),
                )
            if continuation is not None:
                initial_phase = continuation.initial_phase
                if not isinstance(initial_phase, ProviderPhase):
                    raise TypeError("Provider model continuation requires an initial ProviderPhase")
                if not isinstance(initial_phase.model, epcsaft.Mixture):
                    raise TypeError("continuation initial phase requires an epcsaft.Mixture")
                if (
                    not initial_phase.expected_parameter_fingerprint
                    or initial_phase.expected_parameter_fingerprint
                    != initial_phase.model.parameter_fingerprint
                ):
                    raise ValueError("installed continuation Provider fingerprint does not match")
                initial_bounds = tuple(
                    float(value) for value in initial_phase.admissible_packing_fraction_interval
                )
                if (
                    len(initial_bounds) != 2
                    or not all(math.isfinite(value) for value in initial_bounds)
                    or not 0.0 < initial_bounds[0] < initial_bounds[1]
                ):
                    raise ValueError(
                        "continuation packing-fraction bounds must be finite and increasing"
                    )
                initial_constants = continuation.initial_equilibrium_constants
                initial_constants_fingerprint = (
                    continuation.initial_equilibrium_constants_model_fingerprint
                )
                if (initial_constants is None) != (initial_constants_fingerprint is None):
                    raise ValueError(
                        "initial equilibrium constants and their Provider fingerprint "
                        "must be supplied together"
                    )
                if (
                    initial_constants_fingerprint is not None
                    and initial_constants_fingerprint
                    != initial_phase.expected_parameter_fingerprint
                ):
                    raise ValueError(
                        "initial equilibrium constants are not bound to the initial "
                        "Provider fingerprint"
                    )
                if standard_state is not None and initial_constants is not None:
                    raise ValueError(
                        "source-standard continuation transforms one source contract "
                        "at both endpoints"
                    )
                if (
                    standard_state is None
                    and initial_phase.expected_parameter_fingerprint
                    != phase.expected_parameter_fingerprint
                    and initial_constants is None
                ):
                    raise ValueError(
                        "different Provider endpoints require endpoint-bound initial "
                        "equilibrium constants"
                    )
                if initial_constants is not None and (
                    len(initial_constants) != len(problem.reaction_matrix)
                    or any(
                        not isinstance(record, ChemicalEquilibriumConstant)
                        or not math.isfinite(record.ln_value)
                        or not record.source_id
                        or record.reference_id != "provider-helmholtz-coordinate-basis"
                        or record.reaction_orientation != "products_positive"
                        or record.conversion_id != "already-provider-basis"
                        or not record.dimensionless
                        for record in initial_constants
                    )
                ):
                    raise ValueError(
                        "initial equilibrium constants are not a complete "
                        "Provider-basis endpoint contract"
                    )
                if standard_state is not None:
                    initial_source_transfer = _source_reference_transfer(
                        initial_phase,
                        standard_state,
                        species_ids=problem.species_ids,
                        temperature_k=temperature_k,
                        pressure_pa=pressure_pa,
                        sensitivities_requested=False,
                        active_parameters_requested=False,
                    )
        else:
            raise TypeError("chemical_equilibrium requires a typed phase model")
        if standard_state is not None and source_transfer is None:
            raise RuntimeError("source-reference transfer was not evaluated")
        native_standard_state: dict[str, object] | None = (
            None
            if standard_state is None
            else {
                "id": standard_state.id,
                "activity_scale_id": standard_state.activity_scale_id,
                "log_activity_scale_factors": standard_state.log_activity_scale_factors,
                "reference_pressure_pa": standard_state.reference_pressure_pa,
                "source_reference_activity_convention_id": (
                    standard_state.source_reference_activity_convention_id
                ),
                "source_reference_standard_molality_mol_per_kg": (
                    standard_state.source_reference_standard_molality_mol_per_kg
                ),
                "provider_transfer": _source_reference_payload(
                    cast(epcsaft.SourceReferenceTransfer, source_transfer)
                ),
                **(
                    {}
                    if initial_source_transfer is None
                    else {
                        "initial_provider_transfer": _source_reference_payload(
                            initial_source_transfer
                        )
                    }
                ),
            }
        )
        native = cast(
            Mapping[str, object],
            (
                _equilibrium._chemical_equilibrium(
                    capsule,
                    spec,
                    native_standard_state,
                    packing_bounds,
                    problem.strict_interior_amount_floor_mol,
                    None if sensitivity_request is None else active_parameter_payload,
                )
                if continuation is None
                else _equilibrium._chemical_equilibrium_continuation(
                    capsule,
                    epcsaft.native_sdk(continuation.initial_phase.model),
                    spec,
                    continuation.initial_phase.expected_parameter_fingerprint,
                    (
                        None
                        if continuation.initial_equilibrium_constants is None
                        else {
                            "provider_fingerprint": (
                                continuation.initial_equilibrium_constants_model_fingerprint
                            ),
                            "ln_k": tuple(
                                record.ln_value
                                for record in continuation.initial_equilibrium_constants
                            ),
                            "equilibrium_constant_records": tuple(
                                {
                                    "source_id": record.source_id,
                                    "reference_id": record.reference_id,
                                    "reaction_orientation": record.reaction_orientation,
                                    "conversion_id": record.conversion_id,
                                    "dimensionless": record.dimensionless,
                                    "temperature_k": temperature_k,
                                    "pressure_pa": pressure_pa,
                                }
                                for record in continuation.initial_equilibrium_constants
                            ),
                        }
                    ),
                    native_standard_state,
                    cast(tuple[float, float], packing_bounds),
                    continuation.initial_phase.admissible_packing_fraction_interval,
                    problem.strict_interior_amount_floor_mol,
                )
            ),
        )
        diagnostics = _chemical_diagnostics(native)
    except (
        AttributeError,
        IndexError,
        KeyError,
        RuntimeError,
        TypeError,
        ValueError,
    ) as error:
        reason = str(error)
        continuation_fingerprint = (
            getattr(continuation.initial_phase, "expected_parameter_fingerprint", "")
            if continuation is not None
            else ""
        )
        classified = _classified_source_reference_failure(reason) or (
            _classified_continuation_failure(reason, continuation_fingerprint)
            if continuation_fingerprint
            else None
        )
        diagnostics = classified or _failed_chemical_diagnostics(
            "input_or_native_error", reason
        )
        raise ChemicalEquilibriumError(
            reason,
            diagnostics,
        ) from error

    try:
        if not bool(native["accepted"]):
            reason = diagnostics.failure_reason or str(native["callback_error"])
            raise ChemicalEquilibriumError(reason, diagnostics)
        if provider_fingerprint is not None and (
            str(native["parameter_fingerprint"]) != provider_fingerprint
        ):
            raise ValueError("native result has the wrong Provider fingerprint")
        if source_transfer is not None and (
            str(native.get("source_reference_state_fingerprint", ""))
            != source_transfer.reference_state_fingerprint
            or str(native.get("provider_reference_id", ""))
            != source_transfer.helmholtz_basis_id
        ):
            raise ValueError("native result has the wrong source-reference receipt")
        amounts = _vector(native["amounts"], len(problem.species_ids), "chemical amounts")
        total_amount = math.fsum(amounts)
        if total_amount <= 0.0 or not math.isfinite(total_amount):
            raise ValueError("native chemical amounts are invalid")
        volume_m3 = float(cast(float, native["volume_m3"]))
        provider_reference_id = (
            str(native["provider_reference_id"]) if "provider_reference_id" in native else None
        )
        standard_offsets = (
            _float_tuple(native["standard_offsets"]) if "standard_offsets" in native else None
        )
        ln_k_provider_basis = (
            _float_tuple(native["ln_k_provider_basis"]) if "ln_k_provider_basis" in native else None
        )
        artifact_identity = _chemical_artifact_identity(
            native,
            uses_provider=provider_fingerprint is not None,
        )
        sensitivity = (
            None
            if sensitivity_request is None
            else _chemical_sensitivity(
                native,
                species_ids=problem.species_ids,
                model_fingerprint=model_fingerprint,
                active_parameters=active_parameters,
            )
        )
    except ChemicalEquilibriumError:
        raise
    except (KeyError, metadata.PackageNotFoundError, RuntimeError, TypeError, ValueError) as error:
        raise ChemicalEquilibriumError(
            str(error), _failed_chemical_diagnostics("payload_error", str(error))
        ) from error
    return ChemicalEquilibriumResult(
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        species_ids=problem.species_ids,
        charges=problem.charges,
        amounts_mol=amounts,
        mole_fractions=tuple(amount / total_amount for amount in amounts),
        volume_m3=volume_m3,
        thermodynamic_model=thermodynamic_model,
        model_fingerprint=model_fingerprint,
        provider_component_ids=provider_ids,
        provider_parameter_fingerprint=provider_fingerprint,
        equilibrium_constants=problem.equilibrium_constants,
        source_standard_state=standard_state,
        source_reference_transfer=source_transfer,
        provider_reference_id=provider_reference_id,
        standard_offsets=standard_offsets,
        ln_k_provider_basis=ln_k_provider_basis,
        local_scope="fixed_TP_single_homogeneous_phase",
        diagnostics=diagnostics,
        response_kind=(
            "value_only"
            if sensitivity is None
            else (
                "value_plus_jacobian"
                if sensitivity.status == "available"
                else "value_with_unavailable_jacobian"
            )
        ),
        artifact_identity=artifact_identity,
        sensitivity=sensitivity,
    )


def chemical_observation_context(
    phase: ProviderPhase,
    *,
    rows: Sequence[ChemicalObservationRow],
    active_parameters: Sequence[ChemicalEquilibriumActiveParameter],
) -> ChemicalObservationContext:
    """Create one source-bound batch evaluator over immutable liquid states.

    The returned handle owns no chemistry defaults and does not run a Python
    callback. Its native callback re-solves every declared homogeneous state
    and returns ordered values, exact total Jacobian columns, and row status.
    """

    if not isinstance(phase, ProviderPhase):
        raise TypeError("chemical_observation_context requires a ProviderPhase")
    row_records = tuple(rows)
    parameter_records = tuple(active_parameters)
    if not row_records or not parameter_records:
        raise ValueError("chemical-observation context requires rows and parameters")
    if any(not isinstance(item, ChemicalObservationRow) for item in row_records):
        raise TypeError("chemical-observation rows use the typed row record")
    if any(not isinstance(item, ChemicalEquilibriumActiveParameter) for item in parameter_records):
        raise TypeError("active parameters use the typed request record")
    if not isinstance(phase.model, epcsaft.Mixture):
        raise TypeError("ProviderPhase requires an epcsaft.Mixture")
    if (
        not phase.expected_parameter_fingerprint
        or phase.expected_parameter_fingerprint != phase.model.parameter_fingerprint
    ):
        raise ValueError("installed Provider fingerprint does not match the problem")
    if len({item.row_id for item in row_records}) != len(row_records):
        raise ValueError("chemical-observation row identities must be unique")
    packing_bounds = tuple(float(value) for value in phase.admissible_packing_fraction_interval)
    if len(packing_bounds) != 2 or not 0.0 < packing_bounds[0] < packing_bounds[1]:
        raise ValueError("packing-fraction bounds must be finite and increasing")
    native_rows: list[dict[str, object]] = []
    for row in row_records:
        problem = row.problem
        if not isinstance(problem, ChemicalEquilibriumProblem):
            raise TypeError("chemical-observation row requires a typed problem")
        if not isinstance(row.primitive, ChemicalObservationPrimitive):
            raise TypeError("chemical-observation row requires a typed primitive")
        if not row.row_id or not row.state_id or not row.state_schema_id or not row.source_id:
            raise ValueError("chemical-observation row identity is incomplete")
        if row.transform_id not in {"identity", "natural_log"}:
            raise ValueError("unsupported chemical-observation transform")
        temperature_k = _quantity(
            row.temperature,
            "kelvin",
            "temperature",
            "chemical_observation_context",
        )
        pressure_pa = _quantity(
            row.pressure,
            "pascal",
            "pressure",
            "chemical_observation_context",
        )
        if (
            not math.isfinite(problem.strict_interior_amount_floor_mol)
            or problem.strict_interior_amount_floor_mol <= 0.0
        ):
            raise ValueError("minimum admitted amount must be positive and finite")
        constants = tuple(
            {
                "source_id": record.source_id,
                "reference_id": record.reference_id,
                "reaction_orientation": record.reaction_orientation,
                "conversion_id": record.conversion_id,
                "dimensionless": record.dimensionless,
                "temperature_k": temperature_k,
                "pressure_pa": (
                    problem.source_standard_state.reference_pressure_pa
                    if problem.source_standard_state is not None
                    else pressure_pa
                ),
            }
            for record in problem.equilibrium_constants
        )
        problem_record: dict[str, object] = {
            "species_ids": problem.species_ids,
            "charges": problem.charges,
            "provider_fingerprint": phase.expected_parameter_fingerprint,
            "molar_masses_kg_per_mol": problem.molar_masses_kg_per_mol,
            "balance_matrix": problem.balance_matrix,
            "conserved_totals": problem.conserved_totals,
            "reaction_matrix": problem.reaction_matrix,
            "feed_amounts": problem.feed_amounts_mol,
            "ln_k": tuple(record.ln_value for record in problem.equilibrium_constants),
            "equilibrium_constant_records": constants,
            "temperature_k": temperature_k,
            "pressure_pa": pressure_pa,
        }
        standard_state: dict[str, object] | None = None
        if problem.source_standard_state is not None:
            transfer = _source_reference_transfer(
                phase,
                problem.source_standard_state,
                species_ids=problem.species_ids,
                temperature_k=temperature_k,
                pressure_pa=pressure_pa,
                sensitivities_requested=True,
                active_parameters_requested=True,
            )
            standard_state = {
                "id": problem.source_standard_state.id,
                "activity_scale_id": problem.source_standard_state.activity_scale_id,
                "log_activity_scale_factors": (
                    problem.source_standard_state.log_activity_scale_factors
                ),
                "reference_pressure_pa": (problem.source_standard_state.reference_pressure_pa),
                "source_reference_activity_convention_id": (
                    problem.source_standard_state.source_reference_activity_convention_id
                ),
                "source_reference_standard_molality_mol_per_kg": (
                    problem.source_standard_state.source_reference_standard_molality_mol_per_kg
                ),
                "provider_transfer": _source_reference_payload(transfer),
            }
        reference_id = (
            problem.source_standard_state.id
            if problem.source_standard_state is not None
            else "provider-helmholtz-coordinate-basis"
        )
        reference_fingerprint = (
            "sha256:"
            + hashlib.sha256(
                json.dumps(
                    {
                        "reference_id": reference_id,
                        "source_standard_state": standard_state,
                        "equilibrium_constants": constants,
                    },
                    allow_nan=False,
                    separators=(",", ":"),
                    sort_keys=True,
                ).encode("utf-8")
            ).hexdigest()
        )
        native_rows.append(
            {
                "row_id": row.row_id,
                "state_id": row.state_id,
                "state_schema_id": row.state_schema_id,
                "source_id": row.source_id,
                "transform_id": row.transform_id,
                "reference_id": reference_id,
                "reference_fingerprint": reference_fingerprint,
                "problem": problem_record,
                "source_standard_state": standard_state,
                "primitive": {
                    "kind": row.primitive.kind,
                    "component_id": row.primitive.component_id,
                },
            }
        )
    native_parameters = tuple(
        {
            "family": item.family,
            "identity": item.identity,
            "component_ids": item.component_ids,
            "value": item.value,
            "unit": item.unit,
        }
        for item in parameter_records
    )
    provider_version, provider_record = _distribution_identity("epcsaft")
    owner_version, owner_record = _distribution_identity("epcsaft-equilibrium")
    provider_header = _distribution_file_identity(
        "epcsaft",
        "epcsaft/include/epcsaft/native_sdk_v1.h",
    )
    owner_header = _distribution_file_identity(
        "epcsaft-equilibrium",
        "epcsaft_equilibrium/include/epcsaft/regression/evaluator_v1.h",
    )
    provider_artifact_identity = (
        f"epcsaft=={provider_version};RECORD={provider_record};HEADER={provider_header}"
    )
    owner_artifact_identity = (
        f"epcsaft-equilibrium=={owner_version};RECORD={owner_record};HEADER={owner_header}"
    )
    evaluator_identity = "epcsaft-equilibrium.homogeneous-liquid-observation.v1"
    capability_id = "homogeneous-liquid-positive-scalars-v1"
    capability_fingerprint = (
        "sha256:"
        + hashlib.sha256(
            (
                capability_id
                + "\nneutral_component_fugacity_pa:Pa"
                + "\nspecies_mole_fraction:dimensionless"
                + "\nidentity\nnatural_log"
            ).encode("utf-8")
        ).hexdigest()
    )
    contract_fingerprint = (
        "sha256:"
        + hashlib.sha256(
            json.dumps(
                {
                    "rows": native_rows,
                    "parameters": native_parameters,
                },
                allow_nan=False,
                separators=(",", ":"),
                sort_keys=True,
            ).encode("utf-8")
        ).hexdigest()
    )
    artifact_identity = (
        "sha256:"
        + hashlib.sha256(
            (
                evaluator_identity
                + "\n"
                + capability_fingerprint
                + "\n"
                + provider_artifact_identity
                + "\n"
                + owner_artifact_identity
                + "\n"
                + contract_fingerprint
                + "\n"
                + phase.expected_parameter_fingerprint
            ).encode("utf-8")
        ).hexdigest()
    )
    trace_floors = {float(item.problem.strict_interior_amount_floor_mol) for item in row_records}
    if len(trace_floors) != 1:
        raise ValueError("chemical-observation rows require one trace-floor contract")
    capsule = _equilibrium._chemical_observation_context(
        epcsaft.native_sdk(phase.model),
        tuple(native_rows),
        packing_bounds,
        trace_floors.pop(),
        native_parameters,
        evaluator_identity,
        capability_id,
        capability_fingerprint,
        provider_artifact_identity,
        owner_artifact_identity,
        contract_fingerprint,
        artifact_identity,
    )
    return ChemicalObservationContext(
        _native_capsule=capsule,
        row_ids=tuple(item.row_id for item in row_records),
        state_ids=tuple(item.state_id for item in row_records),
        state_schema_ids=tuple(item.state_schema_id for item in row_records),
        source_ids=tuple(item.source_id for item in row_records),
        primitive_ids=tuple(
            f"{item.primitive.kind};{item.primitive.component_id}" for item in row_records
        ),
        parameter_ids=tuple(
            f"{item.family};{item.identity};{','.join(item.component_ids)}"
            for item in parameter_records
        ),
        primitive_units=tuple(
            ("Pa" if item.primitive.kind == "neutral_component_fugacity_pa" else "dimensionless")
            for item in row_records
        ),
        transform_ids=tuple(item.transform_id for item in row_records),
        reference_ids=tuple(
            (
                item.problem.source_standard_state.id
                if item.problem.source_standard_state is not None
                else "provider-helmholtz-coordinate-basis"
            )
            for item in row_records
        ),
        reference_fingerprints=tuple(
            str(record["reference_fingerprint"]) for record in native_rows
        ),
        parameter_units=tuple(item.unit for item in parameter_records),
        evaluator_identity=evaluator_identity,
        capability_id=capability_id,
        capability_fingerprint=capability_fingerprint,
        provider_artifact_identity=provider_artifact_identity,
        owner_artifact_identity=owner_artifact_identity,
        contract_fingerprint=contract_fingerprint,
        artifact_identity=artifact_identity,
    )


def saturation(model: epcsaft.Mixture, temperature: Quantity[Any]) -> SaturationResult:
    """Solve and certify one local pure-component saturation boundary."""

    if not isinstance(model, epcsaft.Mixture):
        raise TypeError("saturation requires an epcsaft.Mixture")
    temperature_k = _temperature_in_kelvin(temperature)
    fingerprint = model.parameter_fingerprint
    scope = _SCOPES.get(fingerprint)
    if scope is None:
        raise ValueError("saturation requires an approved pure-component fingerprint")
    if not scope.temperature_min_k <= temperature_k <= scope.temperature_max_k:
        raise ValueError(
            f"{scope.component} temperature is outside its parameter source domain "
            f"[{scope.temperature_min_k:g}, {scope.temperature_max_k:g}] K"
        )

    capsule = epcsaft.native_sdk(model)
    try:
        native = cast(
            Mapping[str, object],
            _equilibrium.solve_saturation(
                capsule,
                temperature_k,
                fingerprint,
                scope.liquid_density_upper_mol_m3,
            ),
        )
    except (RuntimeError, ValueError) as error:
        failure = {
            "solver_converged": False,
            "solver_status": "native_exception",
            "iterations": 0,
            "attempts": 0,
            "attempt_log": (),
            "solver_lower_bounds": (),
            "solver_upper_bounds": (),
            "solver_constraint_violation": math.inf,
            "numerical_converged": False,
            "confirmation_solves": 0,
            "confirmation_max_relative_difference": math.inf,
            "physical_accepted": False,
            "pressure_relative_residual": math.inf,
            "chemical_potential_absolute_residual": math.inf,
            "phase_density_distance": 0.0,
            "exact_derivatives": False,
            "globality_certificate": False,
            "failure_reason": str(error),
        }
        raise SaturationError(str(error), failure) from error
    diagnostics_payload = cast(Mapping[str, object], native["diagnostics"])
    diagnostics = _diagnostics(diagnostics_payload)
    if not bool(native["accepted"]):
        reason = diagnostics.failure_reason or "local saturation boundary was rejected"
        raise SaturationError(reason, diagnostics_payload)

    return SaturationResult(
        temperature_k=temperature_k,
        saturation_pressure_pa=float(cast(float, native["saturation_pressure_pa"])),
        parameter_fingerprint=fingerprint,
        vapor=_phase(cast(Mapping[str, object], native["vapor"])),
        liquid=_phase(cast(Mapping[str, object], native["liquid"])),
        diagnostics=diagnostics,
    )
