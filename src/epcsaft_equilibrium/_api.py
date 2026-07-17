from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
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
    initial_guess: tuple[float, float, float]
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
class FlashDiagnostics:
    """Separate solver, numerical, and local-physical flash evidence."""

    solver_converged: bool
    solver_status: str
    iterations: int
    attempts: int
    attempt_log: tuple[SolverAttemptDiagnostics, ...]
    solver_lower_bounds: tuple[float, ...]
    solver_upper_bounds: tuple[float, ...]
    solver_constraint_violation: float
    numerical_converged: bool
    confirmation_solves: int
    confirmation_max_difference: float
    physical_accepted: bool
    material_balance_max_abs: float
    pressure_stationarity_max_relative: float
    chemical_potential_max_abs: float
    kkt_stationarity_max_abs: float
    phase_density_distance: float
    exact_derivatives: bool
    globality_certificate: bool
    failure_reason: str


@dataclass(frozen=True)
class TwoPhaseFlashResult:
    """Certified local fixed-two-phase methane/ethane TP result."""

    temperature_k: float
    pressure_pa: float
    overall_mole_fractions: tuple[float, float]
    liquid: PhaseState
    vapor: PhaseState
    liquid_phase_fraction: float
    vapor_phase_fraction: float
    total_free_energy_over_rt: float
    parameter_fingerprint: str
    diagnostics: FlashDiagnostics


class FlashError(RuntimeError):
    """Raised when a local fixed-two-phase calculation is rejected."""

    def __init__(self, message: str, diagnostics: Mapping[str, object]) -> None:
        super().__init__(message)
        self.diagnostics = MappingProxyType(dict(diagnostics))


@dataclass(frozen=True)
class _Scope:
    component: str
    temperature_min_k: float
    temperature_max_k: float
    liquid_density_upper_mol_m3: float


_SCOPES = MappingProxyType(
    {
        "sha256:5f836aa84935df70be2e5cffae51b178a7b797c2cee036e9ff47d8097ca94bbf": _Scope(
            component="methane",
            temperature_min_k=97.0,
            temperature_max_k=300.0,
            liquid_density_upper_mol_m3=40_000.0,
        ),
        "sha256:288fbcaa1304881c16f64c3a784eeed19b75c58cca4558f92a21268e5e91258a": _Scope(
            component="ethane",
            temperature_min_k=90.0,
            temperature_max_k=305.0,
            liquid_density_upper_mol_m3=40_000.0,
        ),
        "sha256:9bfbc8d7789e51609945e61dbdf7a020decc8f9e31b408b0977724c7cb3e1551": _Scope(
            component="propane",
            temperature_min_k=85.0,
            temperature_max_k=523.0,
            liquid_density_upper_mol_m3=25_000.0,
        ),
    }
)

_FLASH_FINGERPRINT = "sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286"
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
    return PhaseState(
        amount_mol=float(cast(float, payload["amount_mol"])),
        mole_fractions=_float_tuple(payload.get("mole_fractions", (1.0,))),
        volume_m3=float(cast(float, payload["volume_m3"])),
        molar_density_mol_m3=float(cast(float, payload["molar_density_mol_m3"])),
        pressure_pa=float(cast(float, payload["pressure_pa"])),
        chemical_potential_over_rt=_float_tuple(payload["chemical_potential_over_rt"]),
    )


def _flash_quantity(quantity: object, units: str, name: str) -> float:
    if not isinstance(quantity, Quantity):
        raise TypeError(f"two_phase_flash requires a Pint {name} quantity")
    try:
        value = float(quantity.to(units).magnitude)
    except DimensionalityError as error:
        raise ValueError(f"{name} units must be convertible to {units}") from error
    except (TypeError, ValueError) as error:
        raise TypeError(f"two_phase_flash requires one scalar Pint {name} quantity") from error
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return value


def _flash_feed(overall_mole_fractions: Sequence[float]) -> tuple[float, float]:
    if isinstance(overall_mole_fractions, (str, bytes)):
        raise TypeError("overall mole fractions must be a numeric sequence")
    try:
        values = tuple(float(value) for value in overall_mole_fractions)
    except (TypeError, ValueError) as error:
        raise TypeError("overall mole fractions must be a numeric sequence") from error
    if len(values) != 2:
        raise ValueError("two_phase_flash requires exactly two components")
    if not all(math.isfinite(value) and value > 0.0 for value in values):
        raise ValueError("overall mole fractions must be positive and finite")
    if not math.isclose(sum(values), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError("overall mole fractions must sum to one within 1e-12")
    return (values[0], values[1])


def two_phase_flash(
    model: epcsaft.EPCSAFT,
    temperature: Quantity[Any],
    pressure: Quantity[Any],
    overall_mole_fractions: Sequence[float],
) -> TwoPhaseFlashResult:
    """Solve one local fixed-two-phase methane/ethane TP calculation."""

    if not isinstance(model, epcsaft.EPCSAFT):
        raise TypeError("two_phase_flash requires an epcsaft.EPCSAFT model")
    temperature_k = _flash_quantity(temperature, "kelvin", "temperature")
    pressure_pa = _flash_quantity(pressure, "pascal", "pressure")
    feed = _flash_feed(overall_mole_fractions)
    if model.parameter_fingerprint != _FLASH_FINGERPRINT:
        raise ValueError("two_phase_flash requires the approved methane/ethane fingerprint")
    if not _FLASH_TEMPERATURE_DOMAIN_K[0] <= temperature_k <= _FLASH_TEMPERATURE_DOMAIN_K[1]:
        raise ValueError("temperature is outside the audited May 2015 source domain")
    if not _FLASH_PRESSURE_DOMAIN_PA[0] <= pressure_pa <= _FLASH_PRESSURE_DOMAIN_PA[1]:
        raise ValueError("pressure is outside the audited May 2015 source domain")
    if not _FLASH_METHANE_FEED_DOMAIN[0] <= feed[0] <= _FLASH_METHANE_FEED_DOMAIN[1]:
        raise ValueError("composition is outside the audited May 2015 source domain")

    capsule = epcsaft.native_sdk(model)
    native = _equilibrium._solve_two_phase_flash(capsule, temperature_k, pressure_pa, feed)
    return cast(TwoPhaseFlashResult, native)


def saturation(model: epcsaft.EPCSAFT, temperature: Quantity[Any]) -> SaturationResult:
    """Solve and certify one local pure-component saturation boundary."""

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
