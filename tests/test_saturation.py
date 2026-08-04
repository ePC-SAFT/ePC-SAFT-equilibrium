from __future__ import annotations

import csv
import ctypes
import math
from pathlib import Path

import epcsaft
import pytest
from parameter_dictionaries import PURE_SATURATION_PARAMETERS

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
        ("parameterized_result_size", ctypes.c_size_t),
        ("evaluate_pure_phase_parameters", ctypes.c_void_p),
        ("component_count", ctypes.c_size_t),
        ("mixture_result_size", ctypes.c_size_t),
        ("evaluate_mixture_phase", ctypes.c_void_p),
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


def _model(component: str = "methane") -> epcsaft.Mixture:
    parameters = epcsaft.Parameters.from_dictionary(PURE_SATURATION_PARAMETERS[component])
    return epcsaft.Mixture(parameters)


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
    assert failed.value.diagnostics["exact_derivatives"] is False


def test_public_ethane_saturation_separates_all_acceptance_layers(
    capfd: pytest.CaptureFixture[str],
) -> None:
    result = epcsaft_equilibrium.saturation(
        _model("ethane"),
        240.0 * epcsaft.unit_registry.kelvin,
    )

    assert result.temperature_k == 240.0
    assert result.parameter_fingerprint == (
        "sha256:73aea4044ad3a49a8045861ba88c8e5966ad7b9db99e28f83df90c1d3d456223"
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
    assert len(result.vapor.chemical_potential_over_rt) == 1
    assert result.diagnostics.solver_converged is True
    assert math.isfinite(result.diagnostics.solver_constraint_violation)
    assert result.diagnostics.solver_constraint_violation <= 1.0e-10
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
    assert sorted(anchor["component"] for anchor in anchors) == [
        "ethane",
        "methane",
        "propane",
    ]

    for anchor in anchors:
        component = anchor["component"]
        model = _model(component)
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
