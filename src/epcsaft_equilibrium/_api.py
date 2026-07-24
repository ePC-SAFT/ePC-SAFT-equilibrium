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


def _tp_flash_quantity(quantity: object, units: str, name: str) -> float:
    if not isinstance(quantity, Quantity):
        raise TypeError(f"tp_flash requires a Pint {name} quantity")
    try:
        value = float(quantity.to(units).magnitude)
    except DimensionalityError as error:
        raise ValueError(f"{name} units must be convertible to {units}") from error
    except (TypeError, ValueError) as error:
        raise TypeError(f"tp_flash requires one scalar Pint {name} quantity") from error
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
    stage_i = cast(Mapping[str, object], payload["stage_i"])
    stage_ii_value = payload.get("stage_ii")
    stage_iii_value = payload.get("stage_iii")
    stage_ii = (
        cast(Mapping[str, object], stage_ii_value) if isinstance(stage_ii_value, Mapping) else None
    )
    stage_iii = (
        cast(Mapping[str, object], stage_iii_value)
        if isinstance(stage_iii_value, Mapping)
        else None
    )
    statuses = ("not_adjudicated", "not_adjudicated", "not_adjudicated")
    attempts = int(cast(int, stage_i["completed_evaluation_count"]))
    major_iterations = 0
    lower_bound = upper_bound = held_gap = None
    profiles = [str(stage_i["search_strategy"]), str(stage_i["search_solver"])]

    if stage_ii is not None:
        classification = cast(Mapping[str, object], stage_ii["attempt_classification"])
        declared = int(cast(int, classification["declared"]))
        solver_passed = int(cast(int, classification["solver_converged"])) > 0
        numerical_passed = int(cast(int, classification["physical_kkt_passed"])) > 0
        physical_passed = bool(cast(Sequence[object], stage_ii["candidates"]))
        statuses = (
            "passed" if solver_passed else "failed",
            "passed" if numerical_passed else "failed",
            "passed" if physical_passed else "failed",
        )
        attempts += declared
        major_iterations = int(cast(int, stage_ii["major_iterations"]))
        profiles.extend(
            str(stage_ii[name]) for name in ("search_strategy", "global_explorer", "local_solver")
        )
        bounds = cast(Sequence[Mapping[str, object]], stage_ii["bound_history"])
        if bounds:
            latest = bounds[-1]
            upper_bound = float(cast(float, latest["upper_bound"]))
            if bool(latest["lower_bound_available"]):
                lower_bound = float(cast(float, latest["lower_bound"]))
                held_gap = abs(upper_bound - lower_bound)
    else:
        stage_i_passed = (
            int(cast(int, stage_i["failed_evaluation_count"])) == 0
            and str(stage_i["outcome"]) != "indeterminate"
        )
        statuses = (
            "not_adjudicated",
            "passed" if stage_i_passed else "failed",
            "passed" if stage_i["negative_witness"] is not None else "not_adjudicated",
        )

    material_balance = pressure_residual = kkt_residual = potential_gap = None
    if stage_iii is not None:
        statuses = (
            "passed"
            if str(stage_iii["solver_status"]) in {"solve_succeeded", "solved_to_acceptable_level"}
            else "failed",
            "passed" if str(stage_iii["numerical_status"]) == "converged" else "failed",
            "passed" if str(stage_iii["physical_status"]) == "accepted" else "failed",
        )
        material_balance = max(
            float(cast(float, stage_iii["modified_balance_inf_norm"])),
            float(cast(float, stage_iii["ordinary_balance_inf_norm"])),
        )
        pressure_residual = float(cast(float, stage_iii["pressure_stationarity_inf_norm"]))
        kkt_residual = float(cast(float, stage_iii["kkt_stationarity_inf_norm"]))
        potential_gap = float(cast(float, stage_iii["modified_potential_mixed_gap"]))

    minimum_tpd = stage_i["minimum_tpd"]
    failure = payload.get("failure_reason")
    return HeldDiagnostics(
        outcome=str(payload["outcome"]),
        search_status=(
            str(stage_iii["physical_status"])
            if stage_iii is not None
            else str(stage_ii["outcome"])
            if stage_ii is not None
            else str(stage_i["outcome"])
        ),
        solver_status=statuses[0],
        numerical_status=statuses[1],
        physical_status=statuses[2],
        attempts=attempts,
        major_iterations=major_iterations,
        best_tpd=math.nan if minimum_tpd is None else float(cast(float, minimum_tpd)),
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        held_gap=held_gap,
        material_balance_max_abs=material_balance,
        pressure_stationarity_max_relative=pressure_residual,
        kkt_stationarity_max_abs=kkt_residual,
        chemical_potential_max_relative=potential_gap,
        confirmation_succeeded=False,
        confirmation_max_difference=None,
        search_profiles=tuple(profiles),
        globality_certificate=str(payload["globality_certificate"]),
        failure_reason="" if failure is None else str(failure),
    )


def _held2_phase_state(
    capsule: object,
    temperature_k: float,
    component_count: int,
    fingerprint: str,
    payload: Mapping[str, object],
) -> PhaseState:
    fractions = _vector(
        payload["physical_fractions"], component_count, "HELD2 phase composition"
    )
    phase_fraction = float(cast(float, payload["phase_fraction"]))
    volume = float(cast(float, payload["volume"]))
    if not math.isfinite(phase_fraction) or phase_fraction <= 0.0:
        raise ValueError("native HELD2 phase fraction must be positive and finite")
    if not math.isfinite(volume) or volume <= 0.0:
        raise ValueError("native HELD2 phase volume must be positive and finite")
    provider_phase = cast(
        Mapping[str, object],
        _equilibrium.evaluate_electrolyte_phase(
            capsule, temperature_k, fractions, volume, fingerprint
        ),
    )
    if str(provider_phase["parameter_fingerprint"]) != fingerprint:
        raise ValueError("hydrated HELD2 phase has the wrong Provider fingerprint")
    pressure_pa = float(cast(float, provider_phase["pressure_pa"]))
    if not math.isfinite(pressure_pa):
        raise ValueError("hydrated HELD2 phase pressure must be finite")
    return PhaseState(
        amount_mol=phase_fraction,
        mole_fractions=fractions,
        volume_m3=volume,
        molar_density_mol_m3=1.0 / volume,
        pressure_pa=pressure_pa,
        chemical_potential_over_rt=_vector(
            provider_phase["chemical_potential_over_rt"],
            component_count,
            "HELD2 chemical-potential vector",
        )
    )


def _held2_result(
    capsule: object,
    temperature_k: float,
    pressure_pa: float,
    feed: tuple[float, ...],
    fingerprint: str,
    native: Mapping[str, object],
) -> TpFlashResult:
    diagnostics = _held2_diagnostics(native)
    if diagnostics.globality_certificate != "not_guaranteed":
        raise ValueError("native HELD2 result has an invalid globality certificate")
    if diagnostics.outcome != "physical_equilibrium_accepted":
        reason = diagnostics.failure_reason or "HELD2 did not return a certified equilibrium"
        raise FlashError(reason, diagnostics)
    if str(native["parameter_fingerprint"]) != fingerprint:
        raise ValueError("native HELD2 result has the wrong Provider fingerprint")

    stage_iii_value = native.get("stage_iii")
    if not isinstance(stage_iii_value, Mapping):
        raise ValueError("accepted HELD2 result is missing Stage III evidence")
    stage_iii = cast(Mapping[str, object], stage_iii_value)
    if stage_iii.get("physical_status") != "accepted":
        raise ValueError("accepted HELD2 result has a rejected Stage III status")
    phase_payloads = cast(Sequence[Mapping[str, object]], stage_iii["phases"])
    phases = tuple(
        _held2_phase_state(
            capsule,
            temperature_k,
            len(feed),
            fingerprint,
            payload,
        )
        for payload in phase_payloads
    )
    phase_fractions = tuple(phase.amount_mol for phase in phases)
    if len(phases) < 2 or not math.isclose(
        math.fsum(phase_fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-10
    ):
        raise ValueError("native HELD2 phase result is incomplete")
    objective = float(cast(float, stage_iii["objective"]))
    if not math.isfinite(objective):
        raise ValueError("native HELD2 objective must be finite")
    return TpFlashResult(
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        overall_mole_fractions=feed,
        phases=tuple(phases),
        phase_fractions=phase_fractions,
        total_free_energy_over_rt=objective,
        parameter_fingerprint=fingerprint,
        diagnostics=diagnostics,
    )


def tp_flash(
    model: epcsaft.EPCSAFT,
    temperature: Quantity[Any],
    pressure: Quantity[Any],
    overall_mole_fractions: Sequence[float],
    *,
    trace: bool = False,
) -> TpFlashResult:
    """Run the bounded HELD or strong-electrolyte HELD2 controller."""

    try:
        if not isinstance(model, epcsaft.EPCSAFT):
            raise TypeError("tp_flash requires an epcsaft.EPCSAFT model")
        temperature_k = _tp_flash_quantity(temperature, "kelvin", "temperature")
        pressure_pa = _tp_flash_quantity(pressure, "pascal", "pressure")
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

    if native.get("controller") == "perdomo_held2_steps_1_to_10_v1":
        try:
            return _held2_result(
                capsule,
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
        overall = cast(
            tuple[float, float],
            _vector(native["overall_mole_fractions"], 2, "overall composition"),
        )
        phase_payloads = cast(Sequence[Mapping[str, object]], native["phases"])
        phases = tuple(_phase(payload) for payload in phase_payloads)
        phase_fractions = _float_tuple(native["phase_fractions"])
        expected_phase_count = 1 if diagnostics.outcome == "one_phase" else 2
        if len(phases) != expected_phase_count or len(phase_fractions) != expected_phase_count:
            raise ValueError("native HELD phase count does not match its outcome")
        if not math.isclose(sum(phase_fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-10):
            raise ValueError("native HELD phase fractions do not sum to one")
        if str(native["parameter_fingerprint"]) != _FLASH_FINGERPRINT:
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
