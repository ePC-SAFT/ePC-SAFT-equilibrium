from __future__ import annotations

import ctypes
import math

import epcsaft
import pytest

import epcsaft_equilibrium
from epcsaft_equilibrium import _equilibrium

BINARY_FINGERPRINT = "sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286"
SOURCE_TEMPERATURE_K = 243.58
SOURCE_PRESSURE_PA = 3_949_000.0
SOURCE_LIQUID_X_METHANE = 0.3099
SOURCE_VAPOR_Y_METHANE = 0.6664
SOURCE_FEED = (
    0.5 * (SOURCE_LIQUID_X_METHANE + SOURCE_VAPOR_Y_METHANE),
    1.0 - 0.5 * (SOURCE_LIQUID_X_METHANE + SOURCE_VAPOR_Y_METHANE),
)
COLLAPSED_SOURCE_TEMPERATURE_K = 243.60
COLLAPSED_SOURCE_PRESSURE_PA = 6_885_000.0
COLLAPSED_SOURCE_FEED = (
    0.5 * (0.6218 + 0.7123),
    1.0 - 0.5 * (0.6218 + 0.7123),
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


class _ExtendedNativeSdkTable(ctypes.Structure):
    _fields_ = (*_NativeSdkTable._fields_, ("future_tail", ctypes.c_uint64))


_CAPSULE_NAME_BUFFERS: list[ctypes.Array[ctypes.c_char]] = []
_SDK_TABLES: list[ctypes.Structure] = []


def _capsule(table: ctypes.Structure) -> object:
    name_buffer = ctypes.create_string_buffer(b"epcsaft.native_sdk.v1")
    _CAPSULE_NAME_BUFFERS.append(name_buffer)
    _SDK_TABLES.append(table)
    new_capsule = ctypes.pythonapi.PyCapsule_New
    new_capsule.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p)
    new_capsule.restype = ctypes.py_object
    return new_capsule(ctypes.addressof(table), name_buffer, None)


def _malformed_mixture_capsule(mixture_result_size: int, evaluator: int | None) -> object:
    provider_info = _equilibrium.sdk_info(epcsaft.native_sdk(_binary_model()))
    table = _NativeSdkTable(
        1,
        ctypes.sizeof(_NativeSdkTable),
        int(provider_info["result_size"]),
        1,
        1,
        0,
        None,
        2,
        mixture_result_size,
        evaluator,
    )
    return _capsule(table)


def _binary_model() -> epcsaft.EPCSAFT:
    parameters = epcsaft.ParameterBundle.from_catalog(
        "gross-2001-methane-ethane", version=1
    ).select(("methane", "ethane"))
    model = epcsaft.EPCSAFT(parameters)
    assert model.parameter_fingerprint == BINARY_FINGERPRINT
    return model


def test_native_mixture_transport_uses_reviewed_sdk_tail() -> None:
    model = _binary_model()
    capsule = epcsaft.native_sdk(model)

    info = _equilibrium.sdk_info(capsule)
    assert info["component_count"] == 2
    assert info["has_evaluate_mixture_phase"] is True
    phase = _equilibrium.evaluate_mixture_phase(
        capsule,
        SOURCE_TEMPERATURE_K,
        (0.4, 0.6),
        1.25e-4,
        BINARY_FINGERPRINT,
    )

    assert phase["parameter_fingerprint"] == BINARY_FINGERPRINT
    assert len(phase["gradient"]) == 3
    assert len(phase["hessian"]) == 9
    assert all(math.isfinite(value) for value in phase["gradient"])
    assert all(math.isfinite(value) for value in phase["hessian"])
    assert math.isfinite(phase["helmholtz_over_rt_reference_amount"])
    assert math.isfinite(phase["pressure_pa"])
    with pytest.raises(ValueError, match="component count"):
        _equilibrium.evaluate_mixture_phase(
            capsule, SOURCE_TEMPERATURE_K, (1.0,), 1.25e-4, BINARY_FINGERPRINT
        )
    with pytest.raises(ValueError, match="provider mixture phase evaluation failed"):
        _equilibrium.evaluate_mixture_phase(
            capsule, SOURCE_TEMPERATURE_K, (0.0, 1.0), 1.25e-4, BINARY_FINGERPRINT
        )
    with pytest.raises(ValueError, match="fingerprint"):
        _equilibrium.evaluate_mixture_phase(
            capsule, SOURCE_TEMPERATURE_K, (0.4, 0.6), 1.25e-4, "sha256:wrong"
        )


def test_native_mixture_transport_accepts_fixed_prefix_and_extended_v1_tables() -> None:
    model = _binary_model()
    provider_capsule = epcsaft.native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    provider_pointer = get_pointer(provider_capsule, b"epcsaft.native_sdk.v1")
    provider_table = ctypes.cast(provider_pointer, ctypes.POINTER(_NativeSdkTable)).contents
    prefix_size = _NativeSdkTable.evaluate_mixture_phase.offset + ctypes.sizeof(ctypes.c_void_p)

    values = [getattr(provider_table, name) for name, _ in _NativeSdkTable._fields_]
    prefix_table = _NativeSdkTable(*values)
    prefix_table.table_size = prefix_size
    extended_table = _ExtendedNativeSdkTable(*values, 0xA5A5A5A5A5A5A5A5)
    extended_table.table_size = ctypes.sizeof(_ExtendedNativeSdkTable)

    for table in (prefix_table, extended_table):
        capsule = _capsule(table)
        info = _equilibrium.sdk_info(capsule)
        assert info["mixture_prefix_size"] == prefix_size
        assert info["has_evaluate_mixture_phase"] is True
        phase = _equilibrium.evaluate_mixture_phase(
            capsule,
            SOURCE_TEMPERATURE_K,
            (0.4, 0.6),
            1.25e-4,
            BINARY_FINGERPRINT,
        )
        assert phase["parameter_fingerprint"] == BINARY_FINGERPRINT


@pytest.mark.parametrize(
    ("result_size_offset", "evaluator", "message"),
    ((-1, 1, "mixture result size"), (0, None, "mixture phase evaluator")),
)
def test_native_flash_rejects_malformed_mixture_tail_before_solving(
    result_size_offset: int,
    evaluator: int | None,
    message: str,
) -> None:
    provider_info = _equilibrium.sdk_info(epcsaft.native_sdk(_binary_model()))
    capsule = _malformed_mixture_capsule(
        int(provider_info["mixture_result_size"]) + result_size_offset,
        evaluator,
    )

    with pytest.raises(ValueError, match=message):
        _equilibrium._solve_two_phase_flash(
            capsule,
            SOURCE_TEMPERATURE_K,
            SOURCE_PRESSURE_PA,
            SOURCE_FEED,
        )


def test_two_phase_flash_objective_exact_derivatives_and_linear_structure() -> None:
    model = _binary_model()
    capsule = epcsaft.native_sdk(model)
    variables = (
        0.5 * SOURCE_LIQUID_X_METHANE,
        0.5 * (1.0 - SOURCE_LIQUID_X_METHANE),
        math.log(0.5 / 18_000.0),
        0.5 * SOURCE_VAPOR_Y_METHANE,
        0.5 * (1.0 - SOURCE_VAPOR_Y_METHANE),
        math.log(0.5 * 8.31446261815324 * SOURCE_TEMPERATURE_K / SOURCE_PRESSURE_PA),
    )
    direction = (0.03, -0.02, 0.04, -0.01, 0.02, -0.03)
    step = 2.0e-5

    def evaluate(scale: float) -> dict[str, object]:
        point = tuple(
            value + scale * delta for value, delta in zip(variables, direction, strict=True)
        )
        return _equilibrium.evaluate_two_phase_flash_nlp(
            capsule,
            SOURCE_TEMPERATURE_K,
            SOURCE_PRESSURE_PA,
            SOURCE_FEED,
            point,
            BINARY_FINGERPRINT,
        )

    lower = evaluate(-step)
    center = evaluate(0.0)
    upper = evaluate(step)
    assert center["constraints"] == pytest.approx((0.0, 0.0), abs=1.0e-15)
    assert center["jacobian"] == [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0]
    objective_difference = (upper["objective"] - lower["objective"]) / (2.0 * step)
    assert objective_difference == pytest.approx(
        sum(value * delta for value, delta in zip(center["gradient"], direction, strict=True)),
        rel=2.0e-7,
        abs=2.0e-8,
    )

    hessian_lower = center["hessian_lower"]
    hessian = [[0.0] * 6 for _ in range(6)]
    index = 0
    for row in range(6):
        for column in range(row + 1):
            hessian[row][column] = hessian_lower[index]
            hessian[column][row] = hessian_lower[index]
            index += 1
    for row in range(6):
        gradient_difference = (upper["gradient"][row] - lower["gradient"][row]) / (2.0 * step)
        hessian_direction = sum(hessian[row][column] * direction[column] for column in range(6))
        assert gradient_difference == pytest.approx(hessian_direction, rel=2.0e-6, abs=2.0e-7)


def test_native_two_phase_flash_source_anchor_solve_and_local_certificate() -> None:
    model = _binary_model()
    native = _equilibrium._solve_two_phase_flash(
        epcsaft.native_sdk(model),
        SOURCE_TEMPERATURE_K,
        SOURCE_PRESSURE_PA,
        SOURCE_FEED,
    )

    assert native["accepted"] is True
    assert native["parameter_fingerprint"] == BINARY_FINGERPRINT
    diagnostics = native["diagnostics"]
    assert diagnostics["solver_converged"] is True
    assert diagnostics["solver_constraint_violation"] <= 1.0e-10
    assert diagnostics["numerical_converged"] is True
    assert diagnostics["confirmation_solves"] == 1
    assert diagnostics["physical_accepted"] is True
    assert diagnostics["material_balance_max_abs"] <= 1.0e-10
    assert diagnostics["pressure_stationarity_max_relative"] <= 1.0e-8
    assert diagnostics["chemical_potential_max_abs"] <= 1.0e-8
    assert diagnostics["kkt_stationarity_max_abs"] <= 1.0e-8
    assert diagnostics["phase_density_distance"] >= 1.0e-3
    assert diagnostics["exact_derivatives"] is True
    assert diagnostics["globality_certificate"] is False
    assert len(diagnostics["equality_multipliers"]) == 2
    assert len(diagnostics["lower_bound_multipliers"]) == 6
    assert len(diagnostics["upper_bound_multipliers"]) == 6
    assert native["liquid"]["molar_density_mol_m3"] > native["vapor"]["molar_density_mol_m3"]
    assert abs(native["liquid"]["mole_fractions"][0] - SOURCE_LIQUID_X_METHANE) <= 0.0111
    assert abs(native["vapor"]["mole_fractions"][0] - SOURCE_VAPOR_Y_METHANE) <= 0.0114


def test_native_two_phase_flash_rejects_collapsed_may_row_012() -> None:
    native = _equilibrium._solve_two_phase_flash(
        epcsaft.native_sdk(_binary_model()),
        COLLAPSED_SOURCE_TEMPERATURE_K,
        COLLAPSED_SOURCE_PRESSURE_PA,
        COLLAPSED_SOURCE_FEED,
    )

    assert native["accepted"] is False
    diagnostics = native["diagnostics"]
    assert diagnostics["solver_converged"] is True
    assert diagnostics["physical_accepted"] is False
    assert 0.0 <= diagnostics["phase_density_distance"] < 1.0e-3
    assert diagnostics["failure_reason"] == "Ipopt solution failed local physical acceptance"


def test_public_two_phase_flash_returns_typed_source_anchor_result() -> None:
    model = _binary_model()
    native = _equilibrium._solve_two_phase_flash(
        epcsaft.native_sdk(model),
        SOURCE_TEMPERATURE_K,
        SOURCE_PRESSURE_PA,
        SOURCE_FEED,
    )
    result = epcsaft_equilibrium.two_phase_flash(
        model,
        SOURCE_TEMPERATURE_K * epcsaft.unit_registry.kelvin,
        SOURCE_PRESSURE_PA * epcsaft.unit_registry.pascal,
        SOURCE_FEED,
    )

    assert isinstance(result, epcsaft_equilibrium.TwoPhaseFlashResult)
    assert result.overall_mole_fractions == pytest.approx(SOURCE_FEED)
    assert result.liquid.mole_fractions == pytest.approx(native["liquid"]["mole_fractions"])
    assert result.vapor.mole_fractions == pytest.approx(native["vapor"]["mole_fractions"])
    assert result.liquid_phase_fraction == pytest.approx(native["liquid_phase_fraction"])
    assert result.vapor_phase_fraction == pytest.approx(native["vapor_phase_fraction"])
    assert result.total_free_energy_over_rt == pytest.approx(native["total_free_energy_over_rt"])
    assert result.liquid_phase_fraction + result.vapor_phase_fraction == pytest.approx(1.0)
    assert result.diagnostics.material_balance_max_abs <= 1.0e-10
    assert result.diagnostics.kkt_stationarity_max_abs <= 1.0e-8
    assert result.diagnostics.globality_certificate is False


def test_public_two_phase_flash_wraps_native_rejection(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        _equilibrium,
        "_solve_two_phase_flash",
        lambda *_args: {
            "accepted": False,
            "diagnostics": {
                "solver_converged": False,
                "solver_status": "synthetic_rejection",
                "iterations": 0,
                "attempts": 0,
                "attempt_log": [],
                "solver_lower_bounds": [0.0] * 6,
                "solver_upper_bounds": [1.0] * 6,
                "solver_constraint_violation": math.inf,
                "numerical_converged": False,
                "confirmation_solves": 0,
                "confirmation_max_difference": math.inf,
                "physical_accepted": False,
                "material_balance_max_abs": math.inf,
                "pressure_stationarity_max_relative": math.inf,
                "chemical_potential_max_abs": math.inf,
                "kkt_stationarity_max_abs": math.inf,
                "phase_density_distance": 0.0,
                "equality_multipliers": [0.0] * 2,
                "lower_bound_multipliers": [0.0] * 6,
                "upper_bound_multipliers": [0.0] * 6,
                "exact_derivatives": True,
                "globality_certificate": False,
                "failure_reason": "synthetic rejection",
            },
        },
    )
    with pytest.raises(epcsaft_equilibrium.FlashError, match="synthetic rejection") as rejected:
        epcsaft_equilibrium.two_phase_flash(
            _binary_model(),
            SOURCE_TEMPERATURE_K * epcsaft.unit_registry.kelvin,
            SOURCE_PRESSURE_PA * epcsaft.unit_registry.pascal,
            SOURCE_FEED,
        )
    assert rejected.value.diagnostics["solver_status"] == "synthetic_rejection"


def test_public_two_phase_flash_wraps_provider_sdk_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def unavailable_sdk(_model: epcsaft.EPCSAFT) -> object:
        raise RuntimeError("provider SDK unavailable")

    monkeypatch.setattr(epcsaft, "native_sdk", unavailable_sdk)
    with pytest.raises(
        epcsaft_equilibrium.FlashError, match="provider SDK unavailable"
    ) as rejected:
        epcsaft_equilibrium.two_phase_flash(
            _binary_model(),
            SOURCE_TEMPERATURE_K * epcsaft.unit_registry.kelvin,
            SOURCE_PRESSURE_PA * epcsaft.unit_registry.pascal,
            SOURCE_FEED,
        )
    assert rejected.value.diagnostics["solver_status"] == "native_exception"
    assert rejected.value.diagnostics["attempts"] == 0


@pytest.mark.parametrize(
    ("temperature", "pressure", "feed", "error", "message"),
    (
        (243.58, 3_949_000.0, SOURCE_FEED, TypeError, "Pint temperature"),
        (
            243.58 * epcsaft.unit_registry.kelvin,
            3_949_000.0,
            SOURCE_FEED,
            TypeError,
            "Pint pressure",
        ),
        (
            243.58 * epcsaft.unit_registry.kelvin,
            3_949_000.0 * epcsaft.unit_registry.pascal,
            (0.5,),
            ValueError,
            "two components",
        ),
        (
            243.58 * epcsaft.unit_registry.kelvin,
            3_949_000.0 * epcsaft.unit_registry.pascal,
            (0.0, 1.0),
            ValueError,
            "positive",
        ),
        (
            243.58 * epcsaft.unit_registry.kelvin,
            3_949_000.0 * epcsaft.unit_registry.pascal,
            (0.4, 0.5),
            ValueError,
            "sum to one",
        ),
    ),
)
def test_public_two_phase_flash_rejects_invalid_inputs(
    temperature: object,
    pressure: object,
    feed: tuple[float, ...],
    error: type[Exception],
    message: str,
) -> None:
    with pytest.raises(error, match=message):
        epcsaft_equilibrium.two_phase_flash(_binary_model(), temperature, pressure, feed)


def test_public_two_phase_flash_rejects_wrong_fingerprint_and_source_domain() -> None:
    pure = epcsaft.EPCSAFT(
        epcsaft.ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1).select(
            ("methane",)
        )
    )
    with pytest.raises(ValueError, match="approved methane/ethane fingerprint"):
        epcsaft_equilibrium.two_phase_flash(
            pure,
            SOURCE_TEMPERATURE_K * epcsaft.unit_registry.kelvin,
            SOURCE_PRESSURE_PA * epcsaft.unit_registry.pascal,
            SOURCE_FEED,
        )
    with pytest.raises(ValueError, match="source domain"):
        epcsaft_equilibrium.two_phase_flash(
            _binary_model(),
            300.1 * epcsaft.unit_registry.kelvin,
            SOURCE_PRESSURE_PA * epcsaft.unit_registry.pascal,
            SOURCE_FEED,
        )
