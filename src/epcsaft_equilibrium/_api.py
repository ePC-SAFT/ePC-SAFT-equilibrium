from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, cast

import epcsaft
from pint import Quantity

from . import _equilibrium


@dataclass(frozen=True)
class PhaseState:
    """One explicit one-mole phase state returned by the local boundary solve."""

    amount_mol: float
    volume_m3: float
    molar_density_mol_m3: float
    pressure_pa: float
    chemical_potential_over_rt: float


@dataclass(frozen=True)
class SaturationDiagnostics:
    """Separate numerical and physical evidence for one local boundary."""

    solver_converged: bool
    solver_status: str
    iterations: int
    attempts: int
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


def _temperature_in_kelvin(temperature: Quantity[Any]) -> float:
    if not isinstance(temperature, Quantity):
        raise TypeError("saturation requires a Pint temperature quantity")
    try:
        value = float(temperature.to("kelvin").magnitude)
    except (TypeError, ValueError) as error:
        raise TypeError("saturation requires one scalar Pint temperature quantity") from error
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("temperature must be positive and finite")
    return value


def _diagnostics(payload: Mapping[str, object]) -> SaturationDiagnostics:
    return SaturationDiagnostics(
        solver_converged=bool(payload["solver_converged"]),
        solver_status=str(payload["solver_status"]),
        iterations=int(cast(int, payload["iterations"])),
        attempts=int(cast(int, payload["attempts"])),
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
        volume_m3=float(cast(float, payload["volume_m3"])),
        molar_density_mol_m3=float(cast(float, payload["molar_density_mol_m3"])),
        pressure_pa=float(cast(float, payload["pressure_pa"])),
        chemical_potential_over_rt=float(cast(float, payload["chemical_potential_over_rt"])),
    )


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
    native = cast(
        Mapping[str, object],
        _equilibrium.solve_saturation(
            capsule,
            temperature_k,
            fingerprint,
            scope.liquid_density_upper_mol_m3,
        ),
    )
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
