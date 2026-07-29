from __future__ import annotations

import hashlib
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
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
class ChemicalEquilibriumDiagnostics:
    """Independent local chemical-equilibrium admission axes and residuals."""

    solver_status: str
    callback_error: str
    failure_reason: str
    chemical_certification_level: str
    boundary_status: str
    structural_zero_species_indices: tuple[int, ...]
    numerical_status: str
    physical_status: str
    predictive_status: str
    provider_domain_status: str
    local_minimum_status: str
    trace_status: str
    globality_status: str
    balance_inf_norm: float | None
    charge_inf_norm: float | None
    pressure_relative_residual: float | None
    reaction_affinity_inf_norm: float | None
    packing_fraction: float | None
    kkt_stationarity_inf_norm: float | None
    complementarity_inf_norm: float | None
    reference_representation_residual_inf_norm: float | None
    reference_convergence_error: float | None
    reference_derivative_availability: int | None


@dataclass(frozen=True)
class ChemicalEquilibriumActiveParameter:
    """Reserved Provider coordinate request; unavailable without atomic packing tensors."""

    family: str
    identity: str
    component_ids: tuple[str, ...]
    value: float
    unit: str


@dataclass(frozen=True)
class ChemicalEquilibriumSensitivityRequest:
    """Request exact operation columns and optionally reserved Provider coordinates."""

    active_parameters: tuple[ChemicalEquilibriumActiveParameter, ...] = ()


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
        "sha256:905e7a6e22eb1073347575bf833d5aa059d9ccf562e4408cb186d74f580ba36f": _Scope(
            component="methane",
            temperature_min_k=97.0,
            temperature_max_k=300.0,
            liquid_density_upper_mol_m3=40_000.0,
        ),
        "sha256:b81f32e44adb46080dfa91026c6428045e04a219900305767672d0547f9a9fb9": _Scope(
            component="ethane",
            temperature_min_k=90.0,
            temperature_max_k=305.0,
            liquid_density_upper_mol_m3=40_000.0,
        ),
        "sha256:1194db349d0608c89419e70c56ccec9ada2ae0884dd8e64e519e9560e7e8ae42": _Scope(
            component="propane",
            temperature_min_k=85.0,
            temperature_max_k=523.0,
            liquid_density_upper_mol_m3=25_000.0,
        ),
    }
)

_FLASH_FINGERPRINT = "sha256:3a840001adcb8b82f44e48307ad61e566f6a65d9b82d8312299a439dbce09195"
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


def _failed_chemical_diagnostics(
    status: str, failure_reason: str
) -> ChemicalEquilibriumDiagnostics:
    return ChemicalEquilibriumDiagnostics(
        solver_status=status,
        callback_error="",
        failure_reason=failure_reason,
        chemical_certification_level="NOT_EVALUATED",
        boundary_status="not_adjudicated",
        structural_zero_species_indices=(),
        numerical_status="not_adjudicated",
        physical_status="not_adjudicated",
        predictive_status="not_adjudicated",
        provider_domain_status="not_adjudicated",
        local_minimum_status="not_adjudicated",
        trace_status="not_adjudicated",
        globality_status="not_adjudicated",
        balance_inf_norm=None,
        charge_inf_norm=None,
        pressure_relative_residual=None,
        reaction_affinity_inf_norm=None,
        packing_fraction=None,
        kkt_stationarity_inf_norm=None,
        complementarity_inf_norm=None,
        reference_representation_residual_inf_norm=None,
        reference_convergence_error=None,
        reference_derivative_availability=None,
    )


def _chemical_diagnostics(
    native: Mapping[str, object],
) -> ChemicalEquilibriumDiagnostics:
    return ChemicalEquilibriumDiagnostics(
        solver_status=str(native["solver_status"]),
        callback_error=str(native["callback_error"]),
        failure_reason=(
            str(native["callback_error"])
            or f"not certified: {native['chemical_certification_level']}"
        )
        if not bool(native["accepted"])
        else "",
        chemical_certification_level=str(native["chemical_certification_level"]),
        boundary_status=str(native["boundary_status"]),
        structural_zero_species_indices=tuple(
            int(index) for index in cast(Sequence[int], native["structural_zero_species_indices"])
        ),
        numerical_status=str(native["numerical_status"]),
        physical_status=str(native["physical_status"]),
        predictive_status="not_adjudicated",
        provider_domain_status=str(native["provider_domain_status"]),
        local_minimum_status=str(native["local_minimum_status"]),
        trace_status=str(native["trace_status"]),
        globality_status=str(native["globality_certificate"]),
        balance_inf_norm=float(cast(float, native["balance_inf_norm"])),
        charge_inf_norm=float(cast(float, native["charge_inf_norm"])),
        pressure_relative_residual=float(cast(float, native["pressure_relative_residual"])),
        reaction_affinity_inf_norm=float(cast(float, native["reaction_affinity_inf_norm"])),
        packing_fraction=float(cast(float, native["packing_fraction"])),
        kkt_stationarity_inf_norm=float(cast(float, native["kkt_stationarity_inf_norm"])),
        complementarity_inf_norm=float(cast(float, native["complementarity_inf_norm"])),
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
    parameters = tuple(
        _sensitivity_parameter(name, active_parameters) for name in parameter_names
    )
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
                _sensitivity_parameter(
                    _active_parameter_name(parameter), active_parameters
                )
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


def chemical_equilibrium(
    phase: IdealGasPhase | ProviderPhase,
    temperature: Quantity[Any],
    pressure: Quantity[Any],
    problem: ChemicalEquilibriumProblem,
    *,
    sensitivity_request: ChemicalEquilibriumSensitivityRequest | None = None,
) -> ChemicalEquilibriumResult:
    """Solve and certify one local fixed-T,P homogeneous reacting phase."""

    try:
        if not isinstance(problem, ChemicalEquilibriumProblem):
            raise TypeError("chemical_equilibrium requires a typed problem")
        if sensitivity_request is not None and not isinstance(
            sensitivity_request, ChemicalEquilibriumSensitivityRequest
        ):
            raise TypeError("sensitivity_request must be a typed sensitivity request")
        active_parameters = (
            ()
            if sensitivity_request is None
            else sensitivity_request.active_parameters
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
            if (
                not parameter.family
                or not parameter.identity
                or not parameter.component_ids
                or any(not component_id for component_id in parameter.component_ids)
                or not math.isfinite(parameter.value)
                or not parameter.unit
            ):
                raise ValueError("active Provider parameter request is incomplete")
        temperature_k = _quantity(temperature, "kelvin", "temperature", "chemical_equilibrium")
        pressure_pa = _quantity(pressure, "pascal", "pressure", "chemical_equilibrium")
        standard_state = problem.source_standard_state
        if standard_state is not None and (
            not math.isfinite(standard_state.reference_pressure_pa)
            or standard_state.reference_pressure_pa <= 0.0
        ):
            raise ValueError("source-standard-state reference pressure must be positive and finite")
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
        if isinstance(phase, IdealGasPhase):
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
        else:
            raise TypeError("chemical_equilibrium requires a typed phase model")
        native_standard_state: dict[str, object] | None = (
            None
            if standard_state is None
            else {
                "id": standard_state.id,
                "activity_scale_id": standard_state.activity_scale_id,
                "log_activity_scale_factors": standard_state.log_activity_scale_factors,
                "reference_pressure_pa": standard_state.reference_pressure_pa,
            }
        )
        native = cast(
            Mapping[str, object],
            _equilibrium._chemical_equilibrium(
                capsule,
                spec,
                native_standard_state,
                packing_bounds,
                problem.strict_interior_amount_floor_mol,
                None if sensitivity_request is None else active_parameter_payload,
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
        raise ChemicalEquilibriumError(
            str(error),
            _failed_chemical_diagnostics("input_or_native_error", str(error)),
        ) from error

    try:
        if not bool(native["accepted"]):
            reason = str(native["callback_error"]) or (
                f"chemical equilibrium was not certified: "
                f"{diagnostics.chemical_certification_level}"
            )
            raise ChemicalEquilibriumError(reason, diagnostics)
        if provider_fingerprint is not None and (
            str(native["parameter_fingerprint"]) != provider_fingerprint
        ):
            raise ValueError("native result has the wrong Provider fingerprint")
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
