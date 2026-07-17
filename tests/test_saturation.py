from __future__ import annotations

import csv
import ctypes
import gc
import math
from pathlib import Path

import epcsaft
import pytest

import epcsaft_equilibrium
from epcsaft_equilibrium import _equilibrium

ANCHORS = Path(__file__).parents[1] / "data" / "reference" / "pure_saturation_anchors.csv"


class _PhaseBlockResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("helmholtz_over_rt_reference_amount", ctypes.c_double),
        ("gradient", ctypes.c_double * 2),
        ("hessian", ctypes.c_double * 4),
        ("third", ctypes.c_double * 8),
        ("pressure_pa", ctypes.c_double),
        ("chemical_potential_over_rt", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class _NativeSdkTable(ctypes.Structure):
    _fields_ = (
        ("abi_version", ctypes.c_uint32),
        ("table_size", ctypes.c_size_t),
        ("result_size", ctypes.c_size_t),
        ("model_context", ctypes.c_void_p),
        ("evaluate_pure_phase", ctypes.c_void_p),
    )


_CAPSULE_NAME_BUFFERS: list[ctypes.Array[ctypes.c_char]] = []
_EVALUATORS: list[object] = []
_SDK_TABLES: list[_NativeSdkTable] = []


def _capsule(table: _NativeSdkTable, name: str = "epcsaft.native_sdk.v1") -> object:
    name_buffer = ctypes.create_string_buffer(name.encode())
    _CAPSULE_NAME_BUFFERS.append(name_buffer)
    _SDK_TABLES.append(table)
    new_capsule = ctypes.pythonapi.PyCapsule_New
    new_capsule.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p)
    new_capsule.restype = ctypes.py_object
    return new_capsule(ctypes.addressof(table), name_buffer, None)


def _model(component: str = "methane") -> epcsaft.EPCSAFT:
    catalog = "gross-2001-propane" if component == "propane" else "gross-2001-methane-ethane"
    parameters = epcsaft.ParameterBundle.from_catalog(catalog, version=1).select((component,))
    return epcsaft.EPCSAFT(parameters)


def test_native_extension_accepts_and_retains_the_public_provider_capsule() -> None:
    model = _model()
    fingerprint = model.parameter_fingerprint
    capsule = epcsaft.native_sdk(model)

    info = _equilibrium.sdk_info(capsule)

    assert info == {
        "capsule_name": "epcsaft.native_sdk.v1",
        "abi_version": 1,
        "table_size": ctypes.sizeof(_NativeSdkTable),
        "result_size": ctypes.sizeof(_PhaseBlockResult),
        "has_model_context": True,
        "has_evaluate_pure_phase": True,
    }

    del model
    gc.collect()
    assert _equilibrium.sdk_info(capsule) == info
    phase = _equilibrium.evaluate_phase(capsule, 150.0, 1.0, 1.0e-3, fingerprint)
    assert phase["parameter_fingerprint"] == fingerprint


def test_native_extension_rejects_malformed_capsule_contracts_and_results() -> None:
    valid_size = ctypes.sizeof(_NativeSdkTable)
    valid_result_size = ctypes.sizeof(_PhaseBlockResult)

    with pytest.raises(ValueError, match="expected capsule"):
        _equilibrium.sdk_info(_capsule(_NativeSdkTable(), "wrong.name"))
    with pytest.raises(ValueError, match="ABI version"):
        _equilibrium.sdk_info(_capsule(_NativeSdkTable(2, valid_size, valid_result_size, 1, 1)))
    with pytest.raises(ValueError, match="table is smaller"):
        _equilibrium.sdk_info(_capsule(_NativeSdkTable(1, valid_size - 1, valid_result_size, 1, 1)))
    with pytest.raises(ValueError, match="result size"):
        _equilibrium.sdk_info(_capsule(_NativeSdkTable(1, valid_size, valid_result_size - 1, 1, 1)))

    fingerprint = _model().parameter_fingerprint
    callback_type = ctypes.CFUNCTYPE(
        ctypes.c_int,
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.POINTER(_PhaseBlockResult),
    )

    def malformed_result(
        _context: int,
        _temperature: float,
        _amount: float,
        _volume: float,
        result: ctypes.POINTER(_PhaseBlockResult),
    ) -> int:
        result.contents.status = 0
        result.contents.struct_size = 0
        result.contents.parameter_fingerprint = fingerprint.encode()
        return 0

    evaluator = callback_type(malformed_result)
    _EVALUATORS.append(evaluator)
    table = _NativeSdkTable(
        1,
        valid_size,
        valid_result_size,
        1,
        ctypes.cast(evaluator, ctypes.c_void_p).value,
    )
    with pytest.raises(ValueError, match="result struct size"):
        _equilibrium.evaluate_phase(_capsule(table), 150.0, 1.0, 1.0e-3, fingerprint)


@pytest.mark.parametrize("component", ("methane", "ethane", "propane"))
def test_native_extension_evaluates_only_the_expected_provider_fingerprint(component: str) -> None:
    model = _model(component)
    capsule = epcsaft.native_sdk(model)

    phase = _equilibrium.evaluate_phase(
        capsule,
        150.0,
        1.0,
        1.0e-3,
        model.parameter_fingerprint,
    )

    assert phase["parameter_fingerprint"] == model.parameter_fingerprint
    assert phase["status"] == 0
    assert phase["amount_mol"] == 1.0
    assert phase["volume_m3"] == 1.0e-3
    assert len(phase["gradient"]) == 2
    assert len(phase["hessian"]) == 4
    assert len(phase["third"]) == 8
    assert all(math.isfinite(value) for value in phase["gradient"])
    assert math.isfinite(phase["pressure_pa"])
    assert math.isfinite(phase["chemical_potential_over_rt"])

    with pytest.raises(ValueError, match="fingerprint"):
        _equilibrium.evaluate_phase(capsule, 150.0, 1.0, 1.0e-3, "sha256:wrong")

    with pytest.raises(ValueError, match="provider phase evaluation failed"):
        _equilibrium.evaluate_phase(
            capsule,
            150.0,
            0.0,
            1.0e-3,
            model.parameter_fingerprint,
        )


@pytest.mark.parametrize(
    ("component", "temperature_k", "variables"),
    (
        ("methane", 150.0, (math.log(1_000.0), math.log(22_000.0), math.log(1.0e6))),
        ("methane", 150.0, (math.log(1.0e-3), math.log(25_000.0), math.log(1.0))),
        ("ethane", 233.15, (math.log(730.0), math.log(20_000.0), math.log(1.0e6))),
        ("ethane", 233.15, (math.log(1.0e-5), math.log(28_000.0), math.log(1.0))),
        ("propane", 300.0, (math.log(500.0), math.log(14_000.0), math.log(1.0e6))),
        ("propane", 300.0, (math.log(1.0e-5), math.log(22_000.0), math.log(1.0))),
    ),
)
def test_saturation_nlp_exact_derivatives_match_independent_directional_differences(
    component: str,
    temperature_k: float,
    variables: tuple[float, float, float],
) -> None:
    model = _model(component)
    capsule = epcsaft.native_sdk(model)
    multipliers = (0.7, -0.4, 1.2)
    direction = (0.3, -0.2, 0.1)
    step = 2.0e-5

    def evaluate(scale: float) -> dict[str, object]:
        point = tuple(
            value + scale * delta for value, delta in zip(variables, direction, strict=True)
        )
        return _equilibrium.evaluate_nlp(
            capsule,
            temperature_k,
            model.parameter_fingerprint,
            point,
            multipliers,
        )

    lower = evaluate(-step)
    center = evaluate(0.0)
    upper = evaluate(step)
    objective_gradient = center["objective_gradient"]
    assert center["objective"] == 0.0
    assert objective_gradient == [0.0, 0.0, 0.0]
    objective_difference = (upper["objective"] - lower["objective"]) / (2.0 * step)
    assert objective_difference == pytest.approx(
        sum(objective_gradient[index] * direction[index] for index in range(3)),
        abs=1.0e-13,
    )

    jacobian = center["jacobian"]
    for row in range(3):
        residual_difference = (upper["constraints"][row] - lower["constraints"][row]) / (2.0 * step)
        jacobian_direction = sum(
            jacobian[row * 3 + column] * direction[column] for column in range(3)
        )
        assert residual_difference == pytest.approx(jacobian_direction, rel=2.0e-7, abs=2.0e-8)

    def lagrangian_gradient(payload: dict[str, object]) -> tuple[float, float, float]:
        payload_jacobian = payload["jacobian"]
        return tuple(
            sum(multipliers[row] * payload_jacobian[row * 3 + column] for row in range(3))
            for column in range(3)
        )

    lower_gradient = lagrangian_gradient(lower)
    upper_gradient = lagrangian_gradient(upper)
    hessian_lower = center["lagrangian_hessian_lower"]
    hessian = (
        (hessian_lower[0], hessian_lower[1], hessian_lower[3]),
        (hessian_lower[1], hessian_lower[2], hessian_lower[4]),
        (hessian_lower[3], hessian_lower[4], hessian_lower[5]),
    )
    for row in range(3):
        gradient_difference = (upper_gradient[row] - lower_gradient[row]) / (2.0 * step)
        hessian_direction = sum(hessian[row][column] * direction[column] for column in range(3))
        assert gradient_difference == pytest.approx(hessian_direction, rel=2.0e-6, abs=2.0e-7)


def test_public_saturation_rejects_noncanonical_or_out_of_scope_inputs() -> None:
    methane = _model("methane")

    with pytest.raises(TypeError, match="Pint temperature quantity"):
        epcsaft_equilibrium.saturation(methane, 150.0)
    with pytest.raises(ValueError, match="convertible to kelvin"):
        epcsaft_equilibrium.saturation(methane, 150.0 * epcsaft.unit_registry.second)
    with pytest.raises(ValueError, match="source domain"):
        epcsaft_equilibrium.saturation(methane, 96.0 * epcsaft.unit_registry.kelvin)

    binary_parameters = epcsaft.ParameterBundle.from_catalog(
        "gross-2001-methane-ethane", version=1
    ).select(("methane", "ethane"))
    binary = epcsaft.EPCSAFT(binary_parameters)
    with pytest.raises(ValueError, match="approved pure-component fingerprint"):
        epcsaft_equilibrium.saturation(binary, 150.0 * epcsaft.unit_registry.kelvin)

    with pytest.raises(epcsaft_equilibrium.SaturationError) as rejected:
        epcsaft_equilibrium.saturation(methane, 250.0 * epcsaft.unit_registry.kelvin)
    assert rejected.value.diagnostics["solver_converged"] is False
    assert rejected.value.diagnostics["globality_certificate"] is False
    assert rejected.value.diagnostics["attempt_log"]


def test_public_saturation_wraps_native_failures_with_structured_diagnostics(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail_native(*_args: object) -> object:
        raise RuntimeError("synthetic native failure")

    monkeypatch.setattr(_equilibrium, "solve_saturation", fail_native)
    with pytest.raises(epcsaft_equilibrium.SaturationError, match="synthetic") as failed:
        epcsaft_equilibrium.saturation(
            _model("methane"),
            150.0 * epcsaft.unit_registry.kelvin,
        )
    assert failed.value.diagnostics["solver_status"] == "native_exception"
    assert failed.value.diagnostics["physical_accepted"] is False


def test_public_ethane_saturation_separates_all_acceptance_layers(
    capfd: pytest.CaptureFixture[str],
) -> None:
    result = epcsaft_equilibrium.saturation(
        _model("ethane"),
        240.0 * epcsaft.unit_registry.kelvin,
    )

    assert result.temperature_k == 240.0
    assert result.parameter_fingerprint == (
        "sha256:288fbcaa1304881c16f64c3a784eeed19b75c58cca4558f92a21268e5e91258a"
    )
    assert result.saturation_pressure_pa == pytest.approx(969_152.1055945412, rel=5.0e-6)
    assert result.vapor.amount_mol == 1.0
    assert result.liquid.amount_mol == 1.0
    assert result.vapor.volume_m3 > result.liquid.volume_m3 > 0.0
    assert result.vapor.molar_density_mol_m3 < result.liquid.molar_density_mol_m3
    assert result.vapor.pressure_pa == pytest.approx(result.saturation_pressure_pa, rel=1.0e-8)
    assert result.liquid.pressure_pa == pytest.approx(result.saturation_pressure_pa, rel=1.0e-8)
    assert result.vapor.chemical_potential_over_rt == pytest.approx(
        result.liquid.chemical_potential_over_rt,
        abs=1.0e-8,
    )
    assert result.diagnostics.solver_converged is True
    assert result.diagnostics.numerical_converged is True
    assert result.diagnostics.physical_accepted is True
    assert result.diagnostics.exact_derivatives is True
    assert result.diagnostics.globality_certificate is False
    assert result.diagnostics.pressure_relative_residual <= 1.0e-8
    assert result.diagnostics.chemical_potential_absolute_residual <= 1.0e-8
    assert result.diagnostics.phase_density_distance > 1.0e-3
    assert len(result.diagnostics.solver_lower_bounds) == 3
    assert len(result.diagnostics.solver_upper_bounds) == 3
    assert len(result.diagnostics.attempt_log) == result.diagnostics.attempts
    assert {attempt.role for attempt in result.diagnostics.attempt_log} == {
        "search",
        "confirmation",
    }
    for attempt in result.diagnostics.attempt_log:
        assert all(
            lower < initial < upper
            for lower, initial, upper in zip(
                result.diagnostics.solver_lower_bounds,
                attempt.initial_guess,
                result.diagnostics.solver_upper_bounds,
                strict=True,
            )
        )
        assert attempt.callback_error == ""
    captured = capfd.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_public_saturation_matches_retained_lab_and_nist_anchors() -> None:
    with ANCHORS.open(encoding="utf-8", newline="") as stream:
        anchors = list(csv.DictReader(stream))

    for anchor in anchors:
        component = anchor["component"]
        model = epcsaft.EPCSAFT(
            epcsaft.ParameterBundle.from_catalog(anchor["catalog"], version=1).select((component,))
        )
        result = epcsaft_equilibrium.saturation(
            model,
            float(anchor["T_K"]) * epcsaft.unit_registry.kelvin,
        )
        assert result.parameter_fingerprint == anchor["parameter_fingerprint"]
        assert result.saturation_pressure_pa == pytest.approx(
            float(anchor["lab_pressure_pa"]), rel=5.0e-6
        )
        assert result.vapor.molar_density_mol_m3 == pytest.approx(
            float(anchor["lab_vapor_density_mol_m3"]), rel=5.0e-6
        )
        assert result.liquid.molar_density_mol_m3 == pytest.approx(
            float(anchor["lab_liquid_density_mol_m3"]), rel=5.0e-6
        )
        assert result.saturation_pressure_pa == pytest.approx(
            float(anchor["nist_pressure_pa"]), rel=5.0e-3
        )
        liquid_density_kg_m3 = result.liquid.molar_density_mol_m3 * float(
            anchor["molar_mass_kg_mol"]
        )
        assert liquid_density_kg_m3 == pytest.approx(
            float(anchor["nist_liquid_density_kg_m3"]), rel=1.0e-2
        )
