from __future__ import annotations

import ctypes
import math
from collections.abc import Sequence

import epcsaft
import pytest

import epcsaft_equilibrium
from epcsaft_equilibrium import _equilibrium

FORMULATION_ID = "perdomo-held2.modified-mole.manufactured.v1"
CHARGES = (0.0, 1.0, -1.0)
PHYSICAL_FEED = (0.5, 0.25, 0.25)
CHEMICAL_POTENTIALS = (3.0, -2.0, 4.0)
GAS_CONSTANT_J_PER_MOL_K = 8.31446261815324
KHUDAIDA_COMPONENT_IDS = (
    "water",
    "ethanol",
    "isobutanol",
    "sodium-cation",
    "chloride-anion",
)
KHUDAIDA_AMOUNTS = (0.65905, 0.05095, 0.27015, 0.01985, 0.01985)
KHUDAIDA_CHARGES = (0.0, 0.0, 0.0, 1.0, -1.0)
KHUDAIDA_FINGERPRINT = "sha256:5a59828d86bb29c919513484a26cedaa0f025463aaa7c149ae3d1fbd0eda97ae"
KHUDAIDA_INDEPENDENT_LOWER = (1.0e-10, 1.0e-10, 2.0e-10)
KHUDAIDA_INDEPENDENT_UPPER = (1.0, 1.0, 1.0)
KHUDAIDA_COMPOSITION_SUM_UPPER = 1.0 - 1.0e-10
TABLE3_TEMPERATURE_K = 298.15
TABLE3_PRESSURE_PA = 2_508.0
TABLE3_FEED = (
    0.8321050353538130581,
    0.0839474823230934710,
    0.0839474823230934710,
)


class _Held2MixtureResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("coordinate_count", ctypes.c_size_t),
        ("gradient_capacity", ctypes.c_size_t),
        ("hessian_capacity", ctypes.c_size_t),
        ("gradient", ctypes.POINTER(ctypes.c_double)),
        ("hessian", ctypes.POINTER(ctypes.c_double)),
        ("helmholtz_over_rt_reference_amount", ctypes.c_double),
        ("pressure_pa", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class _Held2NativeSdk(ctypes.Structure):
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
        ("evaluate_mixture_phase_kij", ctypes.c_void_p),
        ("component_ids", ctypes.POINTER(ctypes.c_char_p)),
        ("component_charges", ctypes.POINTER(ctypes.c_int32)),
        ("evaluate_electrolyte_phase", ctypes.c_void_p),
        ("evaluate_molar_volume_bounds", ctypes.c_void_p),
        ("evaluate_packing_fraction", ctypes.c_void_p),
        ("source_temperature_min_k", ctypes.c_double),
        ("source_temperature_max_k", ctypes.c_double),
        ("total_ion_mole_fraction_max", ctypes.c_double),
    )


_ELECTROLYTE_CALLBACK = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.c_double,
    ctypes.POINTER(_Held2MixtureResult),
)
_VOLUME_BOUNDS_CALLBACK = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
)
_PACKING_CALLBACK = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
)


def _native_sdk_table(capsule: object) -> _Held2NativeSdk:
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    pointer = get_pointer(capsule, b"epcsaft.native_sdk.v1")
    assert pointer
    return _Held2NativeSdk.from_buffer_copy(
        ctypes.string_at(pointer, ctypes.sizeof(_Held2NativeSdk))
    )


def _sdk_capsule(
    table: _Held2NativeSdk,
    *owners: object,
) -> tuple[object, list[object]]:
    name = ctypes.create_string_buffer(b"epcsaft.native_sdk.v1")
    new_capsule = ctypes.pythonapi.PyCapsule_New
    new_capsule.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p)
    new_capsule.restype = ctypes.py_object
    capsule = new_capsule(ctypes.addressof(table), name, None)
    return capsule, [table, name, capsule, *owners]


def _domain_spy_capsule(model: epcsaft.EPCSAFT) -> tuple[object, dict[str, int], list[object]]:
    provider_capsule = epcsaft.native_sdk(model)
    table = _native_sdk_table(provider_capsule)
    charges = tuple(table.component_charges[index] for index in range(table.component_count))
    cap = table.total_ion_mole_fraction_max
    calls = {
        "electrolyte": 0,
        "volume_bounds": 0,
        "packing": 0,
        "outside": 0,
        "outside_volume": 0,
        "at_or_below_volume_lower": 0,
        "at_or_above_volume_upper": 0,
    }

    def record_domain(amounts: ctypes.POINTER(ctypes.c_double), count: int) -> None:
        values = tuple(amounts[index] for index in range(count))
        total = sum(values)
        charged = sum(value for value, charge in zip(values, charges, strict=True) if charge)
        if math.isfinite(cap) and charged / total > cap:
            calls["outside"] += 1

    original_electrolyte = _ELECTROLYTE_CALLBACK(table.evaluate_electrolyte_phase)
    original_volume_bounds = _VOLUME_BOUNDS_CALLBACK(table.evaluate_molar_volume_bounds)
    original_packing = _PACKING_CALLBACK(table.evaluate_packing_fraction)

    def record_volume_domain(
        context: int,
        fingerprint: bytes,
        temperature: float,
        amounts: ctypes.POINTER(ctypes.c_double),
        count: int,
        volume: float,
    ) -> None:
        bounds = (ctypes.c_double * 2)()
        status = original_volume_bounds(
            context,
            fingerprint,
            temperature,
            amounts,
            count,
            1.0e-6,
            0.74,
            bounds,
            2,
        )
        if status == 0:
            below = volume <= bounds[0]
            above = volume >= bounds[1]
            calls["at_or_below_volume_lower"] += int(below)
            calls["at_or_above_volume_upper"] += int(above)
            calls["outside_volume"] += int(below or above)

    @_ELECTROLYTE_CALLBACK
    def electrolyte(
        context: int, temperature: float, amounts: object, count: int, volume: float, result: object
    ) -> int:
        calls["electrolyte"] += 1
        record_domain(amounts, count)
        record_volume_domain(
            context,
            model.parameter_fingerprint.encode(),
            temperature,
            amounts,
            count,
            volume,
        )
        return original_electrolyte(context, temperature, amounts, count, volume, result)

    @_VOLUME_BOUNDS_CALLBACK
    def volume_bounds(
        context: int,
        fingerprint: bytes,
        temperature: float,
        amounts: object,
        count: int,
        eta_min: float,
        eta_max: float,
        bounds: object,
        bound_count: int,
    ) -> int:
        calls["volume_bounds"] += 1
        record_domain(amounts, count)
        return original_volume_bounds(
            context, fingerprint, temperature, amounts, count, eta_min, eta_max, bounds, bound_count
        )

    @_PACKING_CALLBACK
    def packing(
        context: int,
        fingerprint: bytes,
        temperature: float,
        amounts: object,
        count: int,
        volume: float,
        value: object,
        gradient: object,
        gradient_count: int,
        hessian: object,
        hessian_count: int,
    ) -> int:
        calls["packing"] += 1
        record_domain(amounts, count)
        record_volume_domain(context, fingerprint, temperature, amounts, count, volume)
        return original_packing(
            context,
            fingerprint,
            temperature,
            amounts,
            count,
            volume,
            value,
            gradient,
            gradient_count,
            hessian,
            hessian_count,
        )

    table.evaluate_electrolyte_phase = ctypes.cast(electrolyte, ctypes.c_void_p).value
    table.evaluate_molar_volume_bounds = ctypes.cast(volume_bounds, ctypes.c_void_p).value
    table.evaluate_packing_fraction = ctypes.cast(packing, ctypes.c_void_p).value
    name = ctypes.create_string_buffer(b"epcsaft.native_sdk.v1")
    new_capsule = ctypes.pythonapi.PyCapsule_New
    new_capsule.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p)
    new_capsule.restype = ctypes.py_object
    capsule = new_capsule(ctypes.addressof(table), name, None)
    owners: list[object] = [
        model,
        provider_capsule,
        table,
        name,
        original_electrolyte,
        original_volume_bounds,
        original_packing,
        electrolyte,
        volume_bounds,
        packing,
    ]
    return capsule, calls, owners


def _held2_step6(
    feed: tuple[float, float, float],
    upper_bound: float,
    multiplier: tuple[float, float, float],
    cuts: tuple[
        tuple[
            tuple[float, float, float],
            float,
            float,
            float,
            tuple[float, float, float, float],
            tuple[float, float, float],
        ],
        ...,
    ],
) -> dict[str, object]:
    return _equilibrium._held2_adapter(
        KHUDAIDA_CHARGES,
        feed,
        upper_bound,
        multiplier,
        cuts,
        "stage_ii_step6",
    )


def _held2_step6_cut(
    feed: tuple[float, float, float],
    multiplier: tuple[float, float, float],
    independent: tuple[float, float, float],
    packing_fraction: float,
    *,
    q_gradient: tuple[float, float, float, float] | None = None,
    fixed_volume_gradient: tuple[float, float, float] | None = None,
    volume: float = 1.0e-4,
) -> tuple[
    tuple[float, float, float],
    float,
    float,
    float,
    tuple[float, float, float, float],
    tuple[float, float, float],
]:
    objective = -sum(
        value * (feed[index] - independent[index]) for index, value in enumerate(multiplier)
    )
    return (
        independent,
        packing_fraction,
        volume,
        objective,
        q_gradient or (*multiplier, 0.0),
        fixed_volume_gradient or multiplier,
    )


def _figiel_brine_model() -> epcsaft.EPCSAFT:
    parameters = epcsaft.ParameterBundle.from_catalog(
        "figiel-2025-reference-electrolytes", version=1
    ).select(("water", "sodium-cation", "chloride-anion"))
    return epcsaft.EPCSAFT(parameters)


def test_public_tp_flash_dispatches_table3_electrolyte_to_held2() -> None:
    model = _figiel_brine_model()

    result = epcsaft_equilibrium.tp_flash(
        model,
        TABLE3_TEMPERATURE_K * epcsaft.unit_registry.kelvin,
        TABLE3_PRESSURE_PA * epcsaft.unit_registry.pascal,
        TABLE3_FEED,
    )

    assert isinstance(result, epcsaft_equilibrium.TpFlashResult)
    assert result.parameter_fingerprint == model.parameter_fingerprint
    assert result.overall_mole_fractions == pytest.approx(TABLE3_FEED)
    assert result.phase_fractions == pytest.approx((1.0,))
    assert len(result.phases) == 1
    phase = result.phases[0]
    assert phase.amount_mol == pytest.approx(1.0)
    assert phase.mole_fractions == pytest.approx(TABLE3_FEED)
    assert len(phase.chemical_potential_over_rt) == len(TABLE3_FEED)
    assert phase.volume_m3 == pytest.approx(0.9849669199245724, rel=2.0e-12)
    assert phase.pressure_pa == pytest.approx(TABLE3_PRESSURE_PA, rel=1.0e-8)
    assert result.diagnostics.outcome == "one_phase"
    assert result.diagnostics.search_status == "complete_no_negative_found"
    assert result.diagnostics.attempts == 30
    assert result.diagnostics.major_iterations == 0
    assert result.diagnostics.best_tpd == pytest.approx(-1.6139519381498581e-12, abs=2.0e-14)
    assert result.diagnostics.solver_status == "passed"
    assert result.diagnostics.numerical_status == "passed"
    assert result.diagnostics.physical_status == "passed"
    assert result.diagnostics.globality_certificate == "not_guaranteed"


def test_public_tp_flash_rejects_electrolyte_component_mismatch_before_native(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = _figiel_brine_model()
    monkeypatch.setattr(
        epcsaft,
        "native_sdk",
        lambda _model: pytest.fail("component mismatch reached native dispatch"),
    )

    with pytest.raises(epcsaft_equilibrium.FlashError) as failed:
        epcsaft_equilibrium.tp_flash(
            model,
            TABLE3_TEMPERATURE_K * epcsaft.unit_registry.kelvin,
            TABLE3_PRESSURE_PA * epcsaft.unit_registry.pascal,
            TABLE3_FEED[:-1],
        )

    assert failed.value.diagnostics.outcome == "invalid_input"
    assert failed.value.diagnostics.search_status == "input_rejected"
    assert failed.value.diagnostics.solver_status == "not_adjudicated"
    assert failed.value.diagnostics.numerical_status == "not_adjudicated"
    assert failed.value.diagnostics.physical_status == "not_adjudicated"


def test_public_tp_flash_rejects_source_domain_before_provider_phase_calls(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = _figiel_brine_model()
    table = _native_sdk_table(epcsaft.native_sdk(model))
    capsule, calls, owners = _domain_spy_capsule(model)
    assert owners
    monkeypatch.setattr(epcsaft, "native_sdk", lambda _model: capsule)

    with pytest.raises(epcsaft_equilibrium.FlashError) as failed:
        epcsaft_equilibrium.tp_flash(
            model,
            math.nextafter(table.source_temperature_min_k, -math.inf)
            * epcsaft.unit_registry.kelvin,
            TABLE3_PRESSURE_PA * epcsaft.unit_registry.pascal,
            TABLE3_FEED,
        )

    assert failed.value.diagnostics.outcome == "error"
    assert failed.value.diagnostics.search_status == "native_exception"
    assert sum(calls.values()) == 0


def test_public_tp_flash_rejects_provider_component_order_mismatch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = _figiel_brine_model()
    provider_capsule = epcsaft.native_sdk(model)
    table = _native_sdk_table(provider_capsule)
    swapped_ids = (ctypes.c_char_p * 3)(
        b"water",
        b"chloride-anion",
        b"sodium-cation",
    )
    table.component_ids = swapped_ids
    capsule, owners = _sdk_capsule(table, provider_capsule, swapped_ids, model)
    assert owners
    monkeypatch.setattr(epcsaft, "native_sdk", lambda _model: capsule)

    with pytest.raises(epcsaft_equilibrium.FlashError) as failed:
        epcsaft_equilibrium.tp_flash(
            model,
            TABLE3_TEMPERATURE_K * epcsaft.unit_registry.kelvin,
            TABLE3_PRESSURE_PA * epcsaft.unit_registry.pascal,
            TABLE3_FEED,
        )

    assert failed.value.diagnostics.outcome == "error"
    assert failed.value.diagnostics.search_status == "native_exception"
    assert "component order" in str(failed.value)


def test_public_tp_flash_rejects_multicomponent_provider_without_electrolyte_contract(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = _figiel_brine_model()
    provider_capsule = epcsaft.native_sdk(model)
    table = _native_sdk_table(provider_capsule)
    table.table_size = _Held2NativeSdk.evaluate_mixture_phase.offset + ctypes.sizeof(
        ctypes.c_void_p
    )
    capsule, owners = _sdk_capsule(table, provider_capsule, model)
    assert owners
    monkeypatch.setattr(epcsaft, "native_sdk", lambda _model: capsule)

    with pytest.raises(epcsaft_equilibrium.FlashError) as failed:
        epcsaft_equilibrium.tp_flash(
            model,
            TABLE3_TEMPERATURE_K * epcsaft.unit_registry.kelvin,
            TABLE3_PRESSURE_PA * epcsaft.unit_registry.pascal,
            TABLE3_FEED,
        )

    assert failed.value.diagnostics.outcome == "error"
    assert failed.value.diagnostics.search_status == "native_exception"
    assert "electrolyte SDK tail" in str(failed.value)


def test_public_tp_flash_maps_general_mp_accepted_and_error_payloads(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = _figiel_brine_model()
    fractions = (0.2, 0.3, 0.5)
    phases = tuple(
        {
            "amount_mol": fraction,
            "mole_fractions": composition,
            "volume_m3": fraction * molar_volume,
            "molar_density_mol_m3": 1.0 / molar_volume,
            "pressure_pa": TABLE3_PRESSURE_PA,
            "chemical_potential_over_rt": (1.0, 2.0, 3.0),
        }
        for fraction, composition, molar_volume in zip(
            fractions,
            ((0.9, 0.05, 0.05), (0.8, 0.1, 0.1), (0.7, 0.15, 0.15)),
            (0.1, 0.01, 0.001),
            strict=True,
        )
    )
    payload = {
        "outcome": "accepted",
        "search_status": "stage_iii_accepted",
        "failure_reason": "",
        "attempts": 30,
        "major_iterations": 4,
        "best_tpd": -1.0e-3,
        "lower_bound": -25.0,
        "upper_bound": None,
        "held_gap": None,
        "material_balance_max_abs": 1.0e-12,
        "pressure_stationarity_max_relative": 1.0e-10,
        "kkt_stationarity_max_abs": 1.0e-10,
        "chemical_potential_max_relative": 1.0e-10,
        "confirmation_succeeded": False,
        "confirmation_max_difference": None,
        "search_profiles": (
            "perdomo-held2-stage-i-installed-v1",
            "perdomo-held2-stage-ii-installed-v1",
            "perdomo-held2-stage-iii-installed-v1",
        ),
        "solver_status": "passed",
        "numerical_status": "passed",
        "physical_status": "passed",
        "globality_certificate": "not_guaranteed",
        "temperature_k": TABLE3_TEMPERATURE_K,
        "pressure_pa": TABLE3_PRESSURE_PA,
        "overall_mole_fractions": TABLE3_FEED,
        "parameter_fingerprint": model.parameter_fingerprint,
        "phases": phases,
        "phase_fractions": fractions,
        "total_free_energy_over_rt": -25.0,
    }
    monkeypatch.setattr(_equilibrium, "_solve_tp_flash", lambda *_args: payload)

    result = epcsaft_equilibrium.tp_flash(
        model,
        TABLE3_TEMPERATURE_K * epcsaft.unit_registry.kelvin,
        TABLE3_PRESSURE_PA * epcsaft.unit_registry.pascal,
        TABLE3_FEED,
    )

    assert result.phase_fractions == pytest.approx(fractions)
    assert len(result.phases) == 3
    assert all(len(phase.mole_fractions) == 3 for phase in result.phases)
    assert result.diagnostics.outcome == "accepted"
    assert result.diagnostics.globality_certificate == "not_guaranteed"

    payload |= {
        "outcome": "search_exhausted",
        "search_status": "resource_limit",
        "solver_status": "failed",
        "numerical_status": "not_adjudicated",
        "physical_status": "not_adjudicated",
        "failure_reason": "HELD2 Stage II resource limit",
    }
    with pytest.raises(epcsaft_equilibrium.FlashError) as failed:
        epcsaft_equilibrium.tp_flash(
            model,
            TABLE3_TEMPERATURE_K * epcsaft.unit_registry.kelvin,
            TABLE3_PRESSURE_PA * epcsaft.unit_registry.pascal,
            TABLE3_FEED,
        )
    assert failed.value.diagnostics.outcome == "search_exhausted"
    assert failed.value.diagnostics.search_status == "resource_limit"
    assert failed.value.diagnostics.solver_status == "failed"
    assert failed.value.diagnostics.numerical_status == "not_adjudicated"
    assert failed.value.diagnostics.physical_status == "not_adjudicated"


def test_public_tp_flash_does_not_call_private_python_held2_adapter() -> None:
    import inspect

    source = inspect.getsource(epcsaft_equilibrium.tp_flash)
    assert "_held2_adapter" not in source


def _khudaida_model() -> epcsaft.EPCSAFT:
    parameters = epcsaft.ParameterBundle.from_catalog(
        "khudaida-2026-figure-2-electrolyte-lle", version=1
    ).select(KHUDAIDA_COMPONENT_IDS)
    return epcsaft.EPCSAFT(parameters)


def _held2_simplex_chart(values: tuple[float, float, float], stage: str) -> dict[str, object]:
    return _equilibrium._held2_adapter(
        KHUDAIDA_INDEPENDENT_LOWER,
        KHUDAIDA_INDEPENDENT_UPPER,
        KHUDAIDA_COMPOSITION_SUM_UPPER,
        values,
        (),
        (),
        (),
        (0.0, 1.0),
        stage,
    )


def _held2_chart_bound_duals(
    lower: tuple[float, float, float],
    upper: tuple[float, float, float],
    composition_sum_upper: float,
    variables: tuple[float, float, float, float],
    physical_bound_contribution: tuple[float, float, float, float],
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    inverse = _equilibrium._held2_adapter(
        lower,
        upper,
        composition_sum_upper,
        variables[:-1],
        (),
        (),
        (),
        (0.0, 1.0),
        "stage_ii_simplex_inverse",
    )
    jacobian = tuple(inverse["jacobian"])
    dimension = len(lower)
    chart_contribution = (
        *(
            sum(
                jacobian[row * dimension + column] * physical_bound_contribution[row]
                for row in range(dimension)
            )
            for column in range(dimension)
        ),
        physical_bound_contribution[-1],
    )
    return (
        tuple(max(0.0, -value) for value in chart_contribution),
        tuple(max(0.0, value) for value in chart_contribution),
    )


def _held2_dual_pullback_kkt(
    variables: tuple[float, float, float, float],
    physical_gradient: tuple[float, float, float, float],
    master_multiplier: tuple[float, float, float],
    chart_lower_multipliers: tuple[float, ...],
    chart_upper_multipliers: tuple[float, ...],
    *,
    lower: tuple[float, float, float] = KHUDAIDA_INDEPENDENT_LOWER,
    upper: tuple[float, float, float] = KHUDAIDA_INDEPENDENT_UPPER,
    composition_sum_upper: float = KHUDAIDA_COMPOSITION_SUM_UPPER,
) -> dict[str, object]:
    return _equilibrium._held2_adapter(
        lower,
        upper,
        composition_sum_upper,
        variables,
        physical_gradient,
        (),
        master_multiplier,
        (math.log(1.0e-6), math.log(0.74)),
        "stage_ii_physical_kkt",
        chart_lower_multipliers,
        chart_upper_multipliers,
    )


def _manufactured_helmholtz(composition: float, molar_volume: float) -> float:
    shifted = composition - 0.5
    inner_squared = 0.15**2
    outer_squared = 0.30**2
    composition_energy = 1.0e4 * (
        shifted**6 / 6.0
        - (inner_squared + outer_squared) * shifted**4 / 4.0
        + inner_squared * outer_squared * shifted**2 / 2.0
    )
    volume_energy = 0.5 * 5.0 * (molar_volume - 1.2) ** 2
    return composition_energy + volume_energy


def _manufactured_gradient(composition: float, molar_volume: float) -> tuple[float, float]:
    shifted = composition - 0.5
    inner_squared = 0.15**2
    outer_squared = 0.30**2
    return (
        1.0e4
        * (
            shifted**5
            - (inner_squared + outer_squared) * shifted**3
            + inner_squared * outer_squared * shifted
        ),
        5.0 * (molar_volume - 1.2),
    )


def _objective_oracle(vector: Sequence[float]) -> float:
    u_alpha, u_beta, volume_alpha, volume_beta = vector
    fraction = (u_beta - 0.5) / (u_beta - u_alpha)
    gibbs_alpha = _manufactured_helmholtz(u_alpha, volume_alpha) + volume_alpha
    gibbs_beta = _manufactured_helmholtz(u_beta, volume_beta) + volume_beta
    return fraction * gibbs_alpha + (1.0 - fraction) * gibbs_beta


def test_held2_perdomo_transform_lift_and_bounds_match_the_source_coordinates() -> None:
    result = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        (0.27, 0.73, 0.92, 1.08),
        CHEMICAL_POTENTIALS,
    )

    assert result["formulation_id"] == FORMULATION_ID
    assert result["globality_certificate"] == "not_guaranteed"
    assert result["eliminated_index"] == 2
    assert result["dependent_index"] == 0
    assert result["retained_indices"] == [0, 1]
    assert result["independent_indices"] == [1]
    assert result["modified_factors"] == pytest.approx([1.0, 2.0])
    assert result["modified_feed"] == pytest.approx([0.5, 0.5])
    assert result["independent_lower_bounds"] == pytest.approx([2.0e-10])
    assert result["independent_upper_bounds"] == pytest.approx([1.0])
    assert result["phase_fraction"] == pytest.approx(0.5)
    assert result["physical_phases"][0] == pytest.approx([0.73, 0.135, 0.135])
    assert result["physical_phases"][1] == pytest.approx([0.27, 0.365, 0.365])
    assert result["modified_phases"][0] == pytest.approx([0.73, 0.27])
    assert result["modified_phases"][1] == pytest.approx([0.27, 0.73])
    assert result["phase_charge_residuals"] == pytest.approx([0.0, 0.0], abs=1.0e-15)
    assert result["ordinary_balance"] == pytest.approx(PHYSICAL_FEED, abs=1.0e-15)
    assert result["modified_balance"] == pytest.approx([0.5, 0.5], abs=1.0e-15)
    assert result["transformed_modified_potentials"] == pytest.approx([3.0, 1.0])


def test_held2_adapter_direct_objective_and_gradient_match_oracle() -> None:
    center = (0.27, 0.73, 0.92, 1.08)
    direction = (0.03, -0.04, 0.05, -0.02)
    step = 1.0e-5

    result = _equilibrium._held2_adapter(CHARGES, PHYSICAL_FEED, center, CHEMICAL_POTENTIALS)
    lower_vector = tuple(
        value - step * delta for value, delta in zip(center, direction, strict=True)
    )
    upper_vector = tuple(
        value + step * delta for value, delta in zip(center, direction, strict=True)
    )
    lower = _equilibrium._held2_adapter(CHARGES, PHYSICAL_FEED, lower_vector, CHEMICAL_POTENTIALS)
    upper = _equilibrium._held2_adapter(CHARGES, PHYSICAL_FEED, upper_vector, CHEMICAL_POTENTIALS)

    assert result["objective"] == pytest.approx(_objective_oracle(center), abs=2.0e-14)
    fraction = result["phase_fraction"]
    alpha_gradient = _manufactured_gradient(center[0], center[2])
    beta_gradient = _manufactured_gradient(center[1], center[3])
    assert result["phase_gibbs_gradients"][0] == pytest.approx(
        [alpha_gradient[0], alpha_gradient[1] + 1.0]
    )
    assert result["phase_gibbs_gradients"][1] == pytest.approx(
        [beta_gradient[0], beta_gradient[1] + 1.0]
    )
    assert 0.0 < fraction < 1.0
    numerical_directional = (upper["objective"] - lower["objective"]) / (2.0 * step)
    analytic_directional = sum(
        value * delta for value, delta in zip(result["gradient"], direction, strict=True)
    )
    assert numerical_directional == pytest.approx(analytic_directional, rel=2.0e-9, abs=2.0e-10)


def test_held2_adapter_optimum_matches_oracle_certificate() -> None:
    result = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        (0.2, 0.8, 1.0, 1.0),
        CHEMICAL_POTENTIALS,
    )

    assert result["phase_fraction"] == pytest.approx(0.5, abs=1.0e-12)
    assert result["gradient"] == pytest.approx([0.0, 0.0, 0.0, 0.0], abs=1.0e-12)
    assert result["modified_potential_gap"] < 1.0e-12
    assert result["pressure_stationarity_inf_norm"] < 1.0e-12
    assert result["certificate"]["accepted"] is True
    assert result["certificate"]["independent_evidence"] is True
    assert result["certificate"]["metrics"] == pytest.approx(
        {
            "modified_balance_abs": 0.0,
            "ordinary_balance_inf_norm": 0.0,
            "phase_charge_inf_norm": 0.0,
            "modified_potential_gap": 0.0,
            "pressure_stationarity_inf_norm": 0.0,
            "reduced_kkt_inf_norm": 0.0,
            "enumeration_objective_gap": 0.0,
            "independent_modified_composition_count": 1.0,
        },
        abs=1.0e-10,
    )


@pytest.mark.parametrize(
    ("charges", "feed", "message"),
    [
        ((0.0, 2.0, 2.0, -1.0), (0.4, 0.1, 0.1, 0.4), "singular modified coordinate"),
        (CHARGES, (0.5, 0.3, 0.2), "electroneutral"),
        ((0.0, 0.0, 0.0), (0.3, 0.3, 0.4), "charged species"),
    ],
)
def test_held2_transform_rejects_invalid_electrolyte_coordinates(
    charges: tuple[float, ...],
    feed: tuple[float, ...],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        _equilibrium._held2_adapter(
            charges,
            feed,
            (0.2, 0.8, 1.0, 1.0),
            tuple(float(index) for index in range(len(charges))),
        )


def test_held2_adapter_rejects_nonfinite_or_nonbracketing_stage_iii_state() -> None:
    with pytest.raises(ValueError, match="four finite"):
        _equilibrium._held2_adapter(
            CHARGES,
            PHYSICAL_FEED,
            (0.2, 0.8, math.nan, 1.0),
            CHEMICAL_POTENTIALS,
        )
    with pytest.raises(ValueError, match="straddle"):
        _equilibrium._held2_adapter(
            CHARGES,
            PHYSICAL_FEED,
            (0.2, 0.4, 1.0, 1.0),
            CHEMICAL_POTENTIALS,
        )


def test_held2_modified_potentials_are_galvani_gauge_invariant() -> None:
    baseline = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        (0.2, 0.8, 1.0, 1.0),
        CHEMICAL_POTENTIALS,
    )
    shifted_potentials = tuple(
        potential + charge * 17.0
        for potential, charge in zip(CHEMICAL_POTENTIALS, CHARGES, strict=True)
    )
    shifted = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        (0.2, 0.8, 1.0, 1.0),
        shifted_potentials,
    )

    assert shifted["transformed_modified_potentials"] == pytest.approx(
        baseline["transformed_modified_potentials"], abs=1.0e-12
    )


def test_held2_multivalent_coordinates_select_largest_charge_and_round_trip() -> None:
    charges = (0.0, 1.0, 2.0, -1.0, -2.0)
    feed = (0.6, 0.1, 0.1, 0.1, 0.1)
    result = _equilibrium._held2_adapter(
        charges,
        feed,
        (0.0, 1.0, 2.0, 3.0, 4.0),
    )

    assert result["eliminated_index"] == 4
    assert result["dependent_index"] == 0
    assert result["retained_indices"] == [0, 1, 2, 3]
    assert result["independent_indices"] == [1, 2, 3]
    assert result["modified_factors"] == pytest.approx([1.0, 1.5, 2.0, 0.5])
    assert result["modified_feed"] == pytest.approx([0.6, 0.15, 0.2, 0.05])
    assert result["lifted_feed"] == pytest.approx(feed, abs=1.0e-15)
    assert result["independent_lower_bounds"] == pytest.approx([1.5e-10, 2.0e-10, 5.0e-11])
    assert result["independent_upper_bounds"] == pytest.approx([1.0, 1.0, 1.0 / 3.0])
    assert result["transformed_modified_potentials"] == pytest.approx([0.0, 2.0, 3.0, 2.0])


def test_held2_manufactured_enforces_declared_modified_composition_bounds() -> None:
    with pytest.raises(ValueError, match="declared source bounds"):
        _equilibrium._held2_adapter(
            CHARGES,
            PHYSICAL_FEED,
            (1.0e-12, 0.8, 1.0, 1.0),
            CHEMICAL_POTENTIALS,
        )
    with pytest.raises(ValueError, match="declared source bounds"):
        _equilibrium._held2_adapter(
            CHARGES,
            PHYSICAL_FEED,
            (0.2, 1.0 + 1.0e-12, 1.0, 1.0),
            CHEMICAL_POTENTIALS,
        )


def test_held2_phase_block_affine_transform_has_exact_gradient_and_hessian() -> None:
    charges = (0.0, 1.0, 2.0, -1.0, -2.0)
    center = (0.12, 0.18, 0.08, math.log(0.9))
    direction = (0.03, -0.02, 0.01, -0.04)
    pressure_over_rt = 1.3
    target_pressure_pa = 1.3
    linear = (0.2, -0.1, 0.3, -0.2, 0.4, -0.15)
    hessian = tuple(
        tuple(
            (1.0 if row == column else 0.0) + 0.02 * (min(row, column) + 1) * (max(row, column) + 2)
            for column in range(6)
        )
        for row in range(6)
    )

    def physical_coordinates(values: tuple[float, ...]) -> tuple[float, ...]:
        first, second, third, log_volume = values
        modified = (1.0 - first - second - third, first, second, third)
        physical = [0.0] * 5
        physical[0] = modified[0]
        physical[1] = modified[1] / 1.5
        physical[2] = modified[2] / 2.0
        physical[3] = modified[3] / 0.5
        physical[4] = (physical[1] + 2.0 * physical[2] - physical[3]) / 2.0
        return (*physical, math.exp(log_volume))

    def block(values: tuple[float, ...]) -> tuple[float, tuple[float, ...], tuple[float, ...]]:
        coordinates = physical_coordinates(values)
        gradient = tuple(
            linear[row] + sum(hessian[row][column] * coordinates[column] for column in range(6))
            for row in range(6)
        )
        value = 0.7 + sum(linear[index] * coordinates[index] for index in range(6))
        value += 0.5 * sum(
            coordinates[row] * hessian[row][column] * coordinates[column]
            for row in range(6)
            for column in range(6)
        )
        return value, gradient, tuple(value for row in hessian for value in row)

    def evaluate(values: tuple[float, ...]) -> dict[str, object]:
        value, gradient, flat_hessian = block(values)
        return _equilibrium._held2_adapter(
            charges,
            values[:3],
            values[3],
            pressure_over_rt,
            target_pressure_pa,
            value,
            gradient,
            flat_hessian,
            -gradient[-1],
        )

    step = 2.0e-5
    lower = evaluate(
        tuple(value - step * delta for value, delta in zip(center, direction, strict=True))
    )
    result = evaluate(center)
    upper = evaluate(
        tuple(value + step * delta for value, delta in zip(center, direction, strict=True))
    )

    expected_coordinates = physical_coordinates(center)
    assert result["physical_amounts"] == pytest.approx(expected_coordinates[:5], abs=1.0e-15)
    assert result["volume"] == pytest.approx(expected_coordinates[-1], abs=1.0e-15)
    assert result["objective"] == pytest.approx(
        block(center)[0] + pressure_over_rt * expected_coordinates[-1], abs=1.0e-14
    )
    objective_directional = (upper["objective"] - lower["objective"]) / (2.0 * step)
    assert objective_directional == pytest.approx(
        sum(value * delta for value, delta in zip(result["gradient"], direction, strict=True)),
        rel=2.0e-9,
        abs=2.0e-10,
    )
    for row in range(4):
        gradient_directional = (upper["gradient"][row] - lower["gradient"][row]) / (2.0 * step)
        hessian_directional = sum(
            result["hessian"][4 * row + column] * direction[column] for column in range(4)
        )
        assert gradient_directional == pytest.approx(hessian_directional, rel=3.0e-9, abs=3.0e-10)
    assert tuple(result["hessian"]) == pytest.approx(
        tuple(result["hessian"][4 * column + row] for row in range(4) for column in range(4)),
        abs=2.0e-14,
    )


def test_held2_installed_electrolyte_sdk_phase_block_has_exact_reduced_derivatives() -> None:
    model = _figiel_brine_model()
    capsule = epcsaft.native_sdk(model)
    temperature_k = 298.15
    pressure_pa = 100_000.0
    center = (0.02, math.log(1.0e-3))
    direction = (0.007, -0.03)

    def evaluate(values: tuple[float, float]) -> dict[str, object]:
        return _equilibrium._held2_adapter(
            capsule,
            temperature_k,
            pressure_pa,
            values[0:1],
            values[1],
            model.parameter_fingerprint,
        )

    step = 2.0e-5
    lower = evaluate(
        tuple(value - step * delta for value, delta in zip(center, direction, strict=True))
    )
    result = evaluate(center)
    upper = evaluate(
        tuple(value + step * delta for value, delta in zip(center, direction, strict=True))
    )

    assert result["component_ids"] == ["water", "sodium-cation", "chloride-anion"]
    assert result["charges"] == [0.0, 1.0, -1.0]
    assert result["physical_amounts"] == pytest.approx([0.98, 0.01, 0.01], abs=1.0e-15)
    assert result["parameter_fingerprint"] == model.parameter_fingerprint
    assert result["globality_certificate"] == "not_guaranteed"
    objective_directional = (upper["objective"] - lower["objective"]) / (2.0 * step)
    assert objective_directional == pytest.approx(
        sum(value * delta for value, delta in zip(result["gradient"], direction, strict=True)),
        rel=1.0e-9,
        abs=1.0e-10,
    )
    for row in range(2):
        gradient_directional = (upper["gradient"][row] - lower["gradient"][row]) / (2.0 * step)
        hessian_directional = sum(
            result["hessian"][2 * row + column] * direction[column] for column in range(2)
        )
        assert gradient_directional == pytest.approx(hessian_directional, rel=1.0e-8, abs=1.0e-10)
    packing_directional = (upper["packing_fraction"] - lower["packing_fraction"]) / (2.0 * step)
    assert packing_directional == pytest.approx(
        sum(
            value * delta
            for value, delta in zip(result["packing_gradient"], direction, strict=True)
        ),
        rel=1.0e-9,
        abs=1.0e-11,
    )
    for row in range(2):
        packing_gradient_directional = (
            upper["packing_gradient"][row] - lower["packing_gradient"][row]
        ) / (2.0 * step)
        packing_hessian_directional = sum(
            result["packing_hessian"][2 * row + column] * direction[column] for column in range(2)
        )
        assert packing_gradient_directional == pytest.approx(
            packing_hessian_directional, rel=1.0e-8, abs=1.0e-11
        )
    assert tuple(result["packing_hessian"]) == pytest.approx(
        tuple(
            result["packing_hessian"][2 * column + row] for row in range(2) for column in range(2)
        ),
        abs=2.0e-14,
    )
    assert 1.0e-6 <= result["packing_fraction"] <= 0.74
    assert result["pressure_over_rt"] == pytest.approx(
        pressure_pa / (GAS_CONSTANT_J_PER_MOL_K * temperature_k), rel=2.0e-15
    )
    assert math.isfinite(result["provider_pressure_pa"])


def test_held2_installed_log_packing_coordinate_round_trips_and_has_exact_derivatives() -> None:
    model = _khudaida_model()
    capsule = epcsaft.native_sdk(model)
    temperature_k = 293.15
    pressure_pa = 100_000.0
    states = (
        (
            (0.6462224836985831, 0.04995832720498113, 0.03892729322939648),
            math.log(3.909560419950043e-05),
        ),
        (
            (0.48317652139878337, 0.06575359819988703, 0.012557176718000043),
            math.log(5.2811607028165594e-05),
        ),
        (
            (0.49974816254115306, 0.34140273164275964, 0.10555255266393769),
            -2.0847367866112805,
        ),
    )
    direction = (0.011, -0.008, 0.006, 0.017)
    step = 2.0e-4

    for independent, log_volume in states:
        physical = _equilibrium._held2_adapter(
            capsule,
            temperature_k,
            pressure_pa,
            independent,
            log_volume,
            KHUDAIDA_FINGERPRINT,
        )
        center = (*independent, math.log(physical["packing_fraction"]))

        def evaluate(values: tuple[float, ...]) -> dict[str, object]:
            return _equilibrium._held2_adapter(
                capsule,
                temperature_k,
                pressure_pa,
                values[:-1],
                values[-1],
                KHUDAIDA_FINGERPRINT,
                "log_packing_phase",
            )

        lower = evaluate(
            tuple(value - step * delta for value, delta in zip(center, direction, strict=True))
        )
        result = evaluate(center)
        upper = evaluate(
            tuple(value + step * delta for value, delta in zip(center, direction, strict=True))
        )

        assert result["log_volume"] == pytest.approx(log_volume, abs=2.0e-11)
        assert result["log_packing_residual"] == pytest.approx(0.0, abs=2.0e-13)
        assert math.log(result["packing_fraction"]) == pytest.approx(center[-1], abs=2.0e-13)
        objective_directional = (upper["objective"] - lower["objective"]) / (2.0 * step)
        assert objective_directional == pytest.approx(
            sum(value * delta for value, delta in zip(result["gradient"], direction, strict=True)),
            rel=2.0e-8,
            abs=2.0e-9,
        )
        for row in range(4):
            gradient_directional = (upper["gradient"][row] - lower["gradient"][row]) / (2.0 * step)
            hessian_directional = sum(
                result["hessian"][4 * row + column] * direction[column] for column in range(4)
            )
            assert gradient_directional == pytest.approx(
                hessian_directional, rel=3.0e-8, abs=3.0e-9
            )
        assert tuple(result["hessian"]) == pytest.approx(
            tuple(result["hessian"][4 * column + row] for row in range(4) for column in range(4)),
            abs=3.0e-13,
        )


def test_held2_installed_stage_i_log_packing_discriminator_remains_fail_closed() -> None:
    model = _khudaida_model()
    physical_feed = tuple(amount / sum(KHUDAIDA_AMOUNTS) for amount in KHUDAIDA_AMOUNTS)

    result = _equilibrium._held2_adapter(
        epcsaft.native_sdk(model),
        293.15,
        100_000.0,
        physical_feed,
        KHUDAIDA_FINGERPRINT,
        "stage_i_start_0",
    )

    assert result["outcome"] == "indeterminate"
    assert result["attempted_start_count"] == 1
    assert result["completed_start_count"] == 0
    assert result["failed_start_count"] == 1
    assert result["failed_start_index"] == 0
    assert result["failed_start_solver_status"] == 3
    assert result["failed_start_solver_converged"] is False
    assert result["failed_start_reason"] == "TPD solve did not return a complete accepted state"
    assert result["search_completeness"] == "incomplete"


def test_held2_provider_source_domain_metadata_is_prefix_bound() -> None:
    figiel = _figiel_brine_model()
    khudaida = _khudaida_model()

    figiel_info = _equilibrium.sdk_info(epcsaft.native_sdk(figiel))
    khudaida_info = _equilibrium.sdk_info(epcsaft.native_sdk(khudaida))

    assert figiel_info["source_domain_prefix_size"] == ctypes.sizeof(_Held2NativeSdk)
    assert figiel_info["table_size"] >= figiel_info["source_domain_prefix_size"]
    assert figiel_info["source_temperature_min_k"] == figiel_info["source_temperature_max_k"]
    assert math.isfinite(figiel_info["total_ion_mole_fraction_max"])
    assert khudaida_info["table_size"] >= khudaida_info["source_domain_prefix_size"]
    assert khudaida_info["source_temperature_min_k"] == khudaida_info["source_temperature_max_k"]
    assert khudaida_info["total_ion_mole_fraction_max"] is None


def test_held2_stage_i_rejects_source_temperature_and_feed_before_provider() -> None:
    model = _figiel_brine_model()
    table = _native_sdk_table(epcsaft.native_sdk(model))
    cap = table.total_ion_mole_fraction_max
    interior_ion_fraction = 0.5 * cap
    interior_feed = (
        1.0 - interior_ion_fraction,
        0.5 * interior_ion_fraction,
        0.5 * interior_ion_fraction,
    )
    capsule, calls, owners = _domain_spy_capsule(model)
    assert owners

    with pytest.raises(ValueError, match="temperature is outside the Provider source domain"):
        _equilibrium._held2_adapter(
            capsule,
            math.nextafter(table.source_temperature_min_k, -math.inf),
            2_508.0,
            interior_feed,
            model.parameter_fingerprint,
            "stage_i",
        )
    assert sum(calls.values()) == 0

    outside = math.nextafter(cap, math.inf)
    outside_feed = (1.0 - outside, 0.5 * outside, 0.5 * outside)
    with pytest.raises(ValueError, match="ion mole fraction exceeds the Provider source domain"):
        _equilibrium._held2_adapter(
            capsule,
            table.source_temperature_min_k,
            2_508.0,
            outside_feed,
            model.parameter_fingerprint,
            "stage_i",
        )
    assert sum(calls.values()) == 0


@pytest.mark.parametrize("direction", (0.0, -math.inf))
def test_held2_stage_i_keeps_ion_boundary_trials_inside_provider_volume_domain(
    direction: float,
) -> None:
    model = _figiel_brine_model()
    table = _native_sdk_table(epcsaft.native_sdk(model))
    ion_fraction = (
        table.total_ion_mole_fraction_max
        if direction == 0.0
        else math.nextafter(table.total_ion_mole_fraction_max, direction)
    )
    feed = (1.0 - ion_fraction, 0.5 * ion_fraction, 0.5 * ion_fraction)
    capsule, calls, owners = _domain_spy_capsule(model)
    assert owners

    result = _equilibrium._held2_adapter(
        capsule,
        table.source_temperature_min_k,
        2_508.0,
        feed,
        model.parameter_fingerprint,
        "stage_i_start_0",
    )

    assert result["attempted_start_count"] == 1
    assert calls["electrolyte"] > 0
    assert calls["volume_bounds"] > 0
    assert calls["packing"] > 0
    assert calls["outside"] == 0
    assert calls["outside_volume"] == 0
    assert calls["at_or_below_volume_lower"] == 0
    assert calls["at_or_above_volume_upper"] == 0


def test_held2_stage_i_maps_all_declared_starts_deterministically_inside_source_domain() -> None:
    planned_profiles: list[list[list[float]]] = []
    for _ in range(2):
        model = _figiel_brine_model()
        table = _native_sdk_table(epcsaft.native_sdk(model))
        ion_fraction = 0.5 * table.total_ion_mole_fraction_max
        feed = (1.0 - ion_fraction, 0.5 * ion_fraction, 0.5 * ion_fraction)
        capsule, calls, owners = _domain_spy_capsule(model)
        assert owners

        result = _equilibrium._held2_adapter(
            capsule,
            table.source_temperature_min_k,
            2_508.0,
            feed,
            model.parameter_fingerprint,
            "stage_i",
        )

        assert result["declared_start_count"] == 30
        assert result["attempted_start_count"] == 30
        assert calls["outside"] == 0
        assert calls["outside_volume"] == 0
        assert calls["at_or_below_volume_lower"] == 0
        assert calls["at_or_above_volume_upper"] == 0
        assert result["completed_start_count"] == 30
        assert result["failed_start_count"] == 0
        assert len(result["planned_starts"]) == 30
        assert all(
            start[0] <= table.total_ion_mole_fraction_max for start in result["planned_starts"]
        )
        planned_profiles.append(result["planned_starts"])

    assert planned_profiles[0] == planned_profiles[1]


def test_held2_installed_stage_i_uses_provider_volume_domain_for_khudaida_midpoint() -> None:
    model = _khudaida_model()
    physical_feed = tuple(amount / sum(KHUDAIDA_AMOUNTS) for amount in KHUDAIDA_AMOUNTS)

    result = _equilibrium._held2_adapter(
        epcsaft.native_sdk(model),
        293.15,
        100_000.0,
        physical_feed,
        KHUDAIDA_FINGERPRINT,
        "stage_i",
    )

    assert model.parameter_fingerprint == KHUDAIDA_FINGERPRINT
    assert result["component_ids"] == list(KHUDAIDA_COMPONENT_IDS)
    assert result["charges"] == [0.0, 0.0, 0.0, 1.0, -1.0]
    assert result["parameter_fingerprint"] == KHUDAIDA_FINGERPRINT
    assert result["packing_fraction_bounds"] == pytest.approx([1.0e-6, 0.74])
    assert result["molar_volume_bounds"] == pytest.approx(
        [2.306595376485171e-05, 17.068805785990268], rel=2.0e-15
    )
    assert result["log_molar_volume_bounds"] == pytest.approx(
        [math.log(2.306595376485171e-05), math.log(17.068805785990268)], rel=2.0e-15
    )
    assert result["declared_start_count"] == 10 * len(KHUDAIDA_COMPONENT_IDS)
    assert result["reference_scan_interval_count"] == 50
    assert result["reference_scan_point_count"] == 51
    assert result["reference_evaluation_failure_count"] == 0
    assert result["reference_refinement_failure_count"] == 0
    assert result["reference_root_count"] == 3
    assert result["reference_stable_root_count"] == 2
    assert [root["mechanically_stable"] for root in result["reference_roots"]] == [
        True,
        False,
        True,
    ]
    assert [root["volume"] for root in result["reference_roots"]] == pytest.approx(
        [3.909560419950043e-05, 0.0007107348724168745, 0.02238310519168927],
        rel=2.0e-8,
    )
    assert [root["curvature"] for root in result["reference_roots"]] == pytest.approx(
        [31.411491224899734, -0.5642708323753126, 0.8463122865798167],
        rel=2.0e-8,
    )
    # The source-domain guard removes the two former Provider-domain escapes without
    # changing the retained negative witness or the Khudaida scientific outcome.
    assert result["completed_start_count"] == 50
    assert result["failed_start_count"] == 0
    assert result["completed_start_count"] + result["failed_start_count"] == 50
    assert result["search_completeness"] == "complete"
    assert result["failed_start_index"] is None
    assert result["failed_start_solver_status"] is None
    assert result["failed_start_solver_converged"] is None
    assert result["failed_start_reason"] is None
    assert result["outcome"] == "negative_tpd"
    assert result["volume_domain_search_complete"] is True
    assert result["candidate_domain_evaluation_failure_count"] == 0
    assert result["candidate_domain_rejection_count"] == 2
    assert result["reference_volume"] == pytest.approx(3.909560419950043e-05, rel=2.0e-8)
    assert result["minimum_tpd"] == pytest.approx(-0.14499309134029337, abs=2.0e-13)
    assert [candidate["tpd"] for candidate in result["candidates"]] == pytest.approx(
        [-0.00893217913694587], abs=2.0e-13
    )
    assert [candidate["volume"] for candidate in result["candidates"]] == pytest.approx(
        [5.2811607028165594e-05], rel=2.0e-12
    )
    phase_evidence = []
    for candidate in result["candidates"]:
        assert 1.0e-6 <= candidate["packing_fraction"] <= 0.74
        modified = candidate["modified_fractions"]
        phase = _equilibrium._held2_adapter(
            epcsaft.native_sdk(model),
            293.15,
            100_000.0,
            (modified[0], modified[1], modified[3]),
            math.log(candidate["volume"]),
            KHUDAIDA_FINGERPRINT,
        )
        assert sum(phase["physical_amounts"]) == pytest.approx(1.0, abs=2.0e-15)
        assert phase["physical_amounts"][3] == pytest.approx(
            phase["physical_amounts"][4], abs=2.0e-15
        )
        assert (
            max(
                abs(value - feed)
                for value, feed in zip(
                    modified, result["reference_modified_fractions"], strict=True
                )
            )
            > 1.0e-3
        )
        assert candidate["molar_volume_bounds"] == pytest.approx(
            [3.079084484139457e-05, 22.785225182631986], rel=2.0e-14
        )
        assert candidate["molar_volume_bounds"][0] <= candidate["volume"]
        assert candidate["volume"] <= candidate["molar_volume_bounds"][1]
        phase_evidence.append(phase)
    assert abs(phase_evidence[0]["pressure_stationarity_relative"]) <= 1.0e-8
    assert result["candidates"][0]["lower_volume_bound_active"] is False
    assert result["candidates"][0]["upper_volume_bound_active"] is False
    assert result["globality_certificate"] == "not_guaranteed"


def test_held2_manufactured_stage_i_finds_negative_tpd_with_declared_multistart() -> None:
    result = _equilibrium._held2_adapter(CHARGES, PHYSICAL_FEED, "stage_i")

    assert result["profile"] == "perdomo-held2-stage-i-manufactured-v1"
    assert result["outcome"] == "negative_tpd"
    assert result["globality_certificate"] == "not_guaranteed"
    assert result["declared_start_count"] == 10 * len(CHARGES)
    assert result["completed_start_count"] == result["declared_start_count"]
    assert result["failed_start_count"] == 0
    assert result["volume_domain_search_complete"] is True
    assert result["reference_modified_fractions"] == pytest.approx([0.5, 0.5])
    assert result["reference_volume"] == pytest.approx(1.0, abs=2.0e-9)
    assert result["minimum_tpd"] < -1.0e-8
    assert sorted(
        candidate["modified_fractions"][1] for candidate in result["candidates"]
    ) == pytest.approx([0.2, 0.8], abs=2.0e-7)
    assert all(candidate["tpd"] < -1.0e-8 for candidate in result["candidates"])


def test_held2_stage_i_requires_complete_search_before_no_negative_found() -> None:
    result = _equilibrium._held2_adapter(CHARGES, (0.9, 0.05, 0.05), "stage_i")

    assert result["outcome"] == "no_negative_found"
    assert result["completed_start_count"] == result["declared_start_count"] == 30
    assert result["failed_start_count"] == 0
    assert result["volume_domain_search_complete"] is True
    assert result["failed_start_index"] is None
    assert math.isfinite(result["minimum_tpd"])
    assert result["candidates"] == []
    assert result["globality_certificate"] == "not_guaranteed"


def test_held2_manufactured_stage_ii_builds_replayable_candidate_set() -> None:
    result = _equilibrium._held2_adapter(CHARGES, PHYSICAL_FEED, "stage_ii")

    assert result["profile"] == "perdomo-held2-stage-ii-manufactured-v1"
    assert result["outcome"] == "candidate_set"
    assert result["globality_certificate"] == "not_guaranteed"
    assert result["major_iterations"] <= 100
    assert result["lower_starts_per_iteration"] == 30
    assert len(result["bound_history"]) == result["major_iterations"]
    assert all(
        entry["lower_bound"] <= entry["upper_bound"] + 1.0e-10 for entry in result["bound_history"]
    )
    assert result["cut_count"] >= 3
    assert sorted(
        candidate["modified_fractions"][1] for candidate in result["candidates"]
    ) == pytest.approx([0.2, 0.8], abs=2.0e-7)
    assert all(abs(candidate["lower_gap"]) <= 1.0e-8 for candidate in result["candidates"])


def test_held2_stage_ii_certified_improving_cut_advances_with_partial_search() -> None:
    result = _equilibrium._held2_adapter(
        -6.0097420675620405,
        -6.5902644733135745,
        34,
        50,
        "stage_ii_precedence",
    )

    assert result == {
        "decision": "add_improving_cut",
        "search_completeness": "partial",
        "lower_bound_certified": False,
    }


def test_held2_stage_ii_partial_search_without_improving_cut_is_indeterminate() -> None:
    result = _equilibrium._held2_adapter(-6.0, -5.99, 34, 50, "stage_ii_precedence")

    assert result == {
        "decision": "indeterminate",
        "search_completeness": "partial",
        "lower_bound_certified": False,
    }


def test_held2_stage_ii_final_gap_requires_complete_lower_search() -> None:
    result = _equilibrium._held2_adapter(-6.0, -6.0, 50, 50, "stage_ii_precedence")

    assert result == {
        "decision": "evaluate_final_gap",
        "search_completeness": "complete",
        "lower_bound_certified": True,
    }


def test_held2_stage_ii_simplex_chart_round_trips_profile_sized_points_and_preserves_domain() -> (
    None
):
    profile_sized_chart_points = tuple(
        (
            (17 * index % 97 + 0.5) / 98.0,
            (31 * index % 89 + 0.5) / 90.0,
            (43 * index % 83 + 0.5) / 84.0,
        )
        for index in range(50)
    )
    adversarial_chart_points = (
        (0.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (1.0, 0.0, 0.0),
        (math.nextafter(1.0, 0.0), math.nextafter(1.0, 0.0), 1.0),
    )

    for chart in (*profile_sized_chart_points, *adversarial_chart_points):
        forward = _held2_simplex_chart(chart, "stage_ii_simplex_forward")
        physical = tuple(forward["physical"])
        assert all(
            value >= lower
            for value, lower in zip(physical, KHUDAIDA_INDEPENDENT_LOWER, strict=True)
        )
        assert sum(physical) <= KHUDAIDA_COMPOSITION_SUM_UPPER
        if not forward["singular"]:
            inverse = _held2_simplex_chart(physical, "stage_ii_simplex_inverse")
            assert inverse["chart"] == pytest.approx(chart, abs=2.0e-14)
            assert inverse["physical"] == pytest.approx(physical, abs=2.0e-15)

    assert _held2_simplex_chart((1.0, 0.0, 0.0), "stage_ii_simplex_forward")["singular"]

    active_lower = _held2_simplex_chart((0.0, 0.4, 0.3), "stage_ii_simplex_forward")
    active_lower_jacobian = active_lower["jacobian"]
    assert active_lower["physical"][0] == KHUDAIDA_INDEPENDENT_LOWER[0]
    assert active_lower["physical"][1] - KHUDAIDA_INDEPENDENT_LOWER[1] == pytest.approx(
        active_lower_jacobian[4] * 0.4,
        abs=2.0e-16,
    )
    active_simplex = _held2_simplex_chart((0.2, 0.4, 1.0), "stage_ii_simplex_forward")
    active_simplex_jacobian = active_simplex["jacobian"]
    assert KHUDAIDA_COMPOSITION_SUM_UPPER - sum(active_simplex["physical"]) == pytest.approx(
        active_simplex_jacobian[-1] * (1.0 - 1.0),
        abs=2.0e-16,
    )


@pytest.mark.parametrize(
    "independent,log_packing",
    (
        ((0.6462224836985831, 0.04995832720498113, 0.03892729322939648), -0.828757),
        ((0.48317652138888467, 0.06575359820030886, 0.012557176717123258), -0.840619),
        ((0.4406535569650898, 0.08882925535454135, 0.0004486805424638745), -0.8436444167759523),
    ),
)
def test_held2_stage_ii_simplex_chart_exact_provider_gradient_and_hvp(
    independent: tuple[float, float, float], log_packing: float
) -> None:
    model = _khudaida_model()
    capsule = epcsaft.native_sdk(model)
    chart = tuple(_held2_simplex_chart(independent, "stage_ii_simplex_inverse")["chart"])
    center = (*chart, log_packing)
    direction = (0.013, -0.017, 0.011, -0.019)
    master_multiplier = (1.3, -0.7, 2.1)

    def evaluate(values: tuple[float, ...]) -> dict[str, object]:
        mapped = _held2_simplex_chart(
            tuple(values[:-1]),
            "stage_ii_simplex_forward",
        )
        physical = _equilibrium._held2_adapter(
            capsule,
            293.15,
            100_000.0,
            tuple(mapped["physical"]),
            values[-1],
            KHUDAIDA_FINGERPRINT,
            "log_packing_phase",
        )
        lower_gradient = tuple(
            value - master_multiplier[index] if index < 3 else value
            for index, value in enumerate(physical["gradient"])
        )
        transformed = _equilibrium._held2_adapter(
            KHUDAIDA_INDEPENDENT_LOWER,
            KHUDAIDA_INDEPENDENT_UPPER,
            KHUDAIDA_COMPOSITION_SUM_UPPER,
            values,
            lower_gradient,
            tuple(physical["hessian"]),
            (),
            (0.0, 1.0),
            "stage_ii_simplex_chain",
        )
        return {
            "physical": mapped["physical"],
            "objective": physical["objective"]
            - sum(
                value * mapped["physical"][index] for index, value in enumerate(master_multiplier)
            ),
            "packing_fraction": physical["packing_fraction"],
            "gradient": transformed["gradient"],
            "hessian": transformed["hessian"],
        }

    step = 3.0e-6
    lower = evaluate(
        tuple(value - step * delta for value, delta in zip(center, direction, strict=True))
    )
    result = evaluate(center)
    upper = evaluate(
        tuple(value + step * delta for value, delta in zip(center, direction, strict=True))
    )

    assert result["physical"] == pytest.approx(independent, abs=3.0e-15)
    assert result["packing_fraction"] == pytest.approx(math.exp(log_packing), rel=2.0e-12)
    objective_directional = (upper["objective"] - lower["objective"]) / (2.0 * step)
    assert objective_directional == pytest.approx(
        sum(value * delta for value, delta in zip(result["gradient"], direction, strict=True)),
        rel=3.0e-8,
        abs=3.0e-9,
    )
    size = len(center)
    for row in range(size):
        numerical = (upper["gradient"][row] - lower["gradient"][row]) / (2.0 * step)
        analytic = sum(
            result["hessian"][size * row + column] * direction[column] for column in range(size)
        )
        assert numerical == pytest.approx(analytic, rel=8.0e-8, abs=5.0e-9)
    assert tuple(result["hessian"]) == pytest.approx(
        tuple(
            result["hessian"][size * column + row] for row in range(size) for column in range(size)
        ),
        abs=2.0e-14,
    )


def test_held2_stage_ii_simplex_chart_fails_closed_for_nonredundant_upper_bound() -> None:
    with pytest.raises(ValueError, match="requires redundant independent upper bounds"):
        _equilibrium._held2_adapter(
            (1.5e-10, 2.0e-10, 0.5e-10),
            (1.0, 1.0, 1.0 / 3.0),
            1.0 - 1.0e-10,
            (0.2, 0.3, 0.4),
            (),
            (),
            (),
            (0.0, 1.0),
            "stage_ii_simplex_forward",
        )


@pytest.mark.parametrize(
    "variables,bound_contribution",
    (
        ((0.2, 0.3, 0.1, -0.8), (0.0, 0.0, 0.0, 0.0)),
        ((5.1e-9, 0.3, 0.1, -0.8), (-1.0, 0.0, 0.0, 0.0)),
        ((0.2, 0.3, 0.4999999949, -0.8), (0.75, 0.75, 0.75, 0.0)),
    ),
)
def test_held2_stage_ii_dual_pullback_certifies_interior_lower_and_simplex_faces(
    variables: tuple[float, float, float, float],
    bound_contribution: tuple[float, float, float, float],
) -> None:
    master_multiplier = (1.0, 2.0, 3.0)
    physical_gradient = (
        *(master_multiplier[index] - bound_contribution[index] for index in range(3)),
        -bound_contribution[-1],
    )
    chart_lower, chart_upper = _held2_chart_bound_duals(
        KHUDAIDA_INDEPENDENT_LOWER,
        KHUDAIDA_INDEPENDENT_UPPER,
        KHUDAIDA_COMPOSITION_SUM_UPPER,
        variables,
        bound_contribution,
    )

    result = _held2_dual_pullback_kkt(
        variables,
        physical_gradient,
        master_multiplier,
        chart_lower,
        chart_upper,
    )

    assert result["dual_signs_valid"] is True
    assert result["stationarity_inf_norm"] == pytest.approx(0.0, abs=2.0e-14)
    assert result["complementarity"] <= 1.0e-8
    assert result["reconstruction_inf_norm"] <= 1.0e-12


@pytest.mark.parametrize("permutation", ((0, 1, 2), (2, 0, 1), (1, 2, 0)))
def test_held2_stage_ii_dual_pullback_obeys_exact_chain_rule_under_permutation(
    permutation: tuple[int, int, int],
) -> None:
    base_lower = (1.0e-10, 2.0e-10, 3.0e-10)
    base_variables = (2.1e-9, 2.2e-9, 2.3e-9, -0.8)
    base_master = (1.25, -0.5, 2.75)
    base_contribution = (-0.6, -0.2, -0.4, 0.0)
    lower = tuple(base_lower[index] for index in permutation)
    variables = (*(base_variables[index] for index in permutation), base_variables[-1])
    master = tuple(base_master[index] for index in permutation)
    contribution = (*(base_contribution[index] for index in permutation), 0.0)
    physical_gradient = (*(master[index] - contribution[index] for index in range(3)), 0.0)
    chart_lower, chart_upper = _held2_chart_bound_duals(
        lower,
        KHUDAIDA_INDEPENDENT_UPPER,
        KHUDAIDA_COMPOSITION_SUM_UPPER,
        variables,
        contribution,
    )

    result = _held2_dual_pullback_kkt(
        variables,
        physical_gradient,
        master,
        chart_lower,
        chart_upper,
        lower=lower,
    )

    assert result["stationarity_inf_norm"] == pytest.approx(0.0, abs=2.0e-14)
    assert result["reconstruction_inf_norm"] <= 1.0e-12


def test_held2_stage_ii_dual_pullback_rejects_wrong_sign_and_excessive_complementarity() -> None:
    variables = (2.01e-8, 0.3, 0.1, -0.8)
    master = (1.0, 2.0, 3.0)
    contribution = (-1.0, 0.0, 0.0, 0.0)
    gradient = (2.0, 2.0, 3.0, 0.0)
    chart_lower, chart_upper = _held2_chart_bound_duals(
        KHUDAIDA_INDEPENDENT_LOWER,
        KHUDAIDA_INDEPENDENT_UPPER,
        KHUDAIDA_COMPOSITION_SUM_UPPER,
        variables,
        contribution,
    )
    excessive = _held2_dual_pullback_kkt(
        variables,
        gradient,
        master,
        chart_lower,
        chart_upper,
    )
    assert excessive["complementarity"] > 1.0e-8

    wrong_sign_lower = list(chart_lower)
    nonzero = next(index for index, value in enumerate(wrong_sign_lower) if value > 0.0)
    wrong_sign_lower[nonzero] *= -1.0
    wrong_sign = _held2_dual_pullback_kkt(
        variables,
        gradient,
        master,
        tuple(wrong_sign_lower),
        chart_upper,
    )
    assert wrong_sign["dual_signs_valid"] is False
    assert math.isinf(wrong_sign["stationarity_inf_norm"])


def test_held2_stage_ii_dual_pullback_rejects_reconstruction_or_stationarity_failure() -> None:
    variables = (0.2, 0.3, 0.1, -0.8)
    master = (1.0, 2.0, 3.0)
    zero_duals = (0.0, 0.0, 0.0, 0.0)
    stationarity_failure = _held2_dual_pullback_kkt(
        variables,
        (1.0 + 1.01e-7, 2.0, 3.0, 0.0),
        master,
        zero_duals,
        zero_duals,
    )
    assert stationarity_failure["stationarity_inf_norm"] > 1.0e-7

    reconstruction_failure = _held2_dual_pullback_kkt(
        variables,
        (1.0, 2.0, 3.0, 0.0),
        master,
        zero_duals,
        (1.0e16, 0.0, 1.0e16, 0.0),
    )
    assert reconstruction_failure["reconstruction_inf_norm"] > 1.0e-12


def test_held2_stage_ii_dual_pullback_rejects_nonfinite_or_out_of_domain_evidence() -> None:
    zero_duals = (0.0, 0.0, 0.0, 0.0)
    with pytest.raises(ValueError, match="must be finite"):
        _held2_dual_pullback_kkt(
            (0.2, 0.3, 0.1, -0.8),
            (math.nan, 2.0, 3.0, 0.0),
            (1.0, 2.0, 3.0),
            zero_duals,
            zero_duals,
        )

    with pytest.raises(ValueError, match=r"outside (its physical domain|the feasible simplex)"):
        _held2_dual_pullback_kkt(
            (0.6, 0.3, 0.2, -0.8),
            (1.0, 2.0, 3.0, 0.0),
            (1.0, 2.0, 3.0),
            zero_duals,
            zero_duals,
        )


def test_held2_stage_ii_step6_selector_admits_multiple_distinct_cuts() -> None:
    feed = (0.4, 0.1, 0.04)
    multiplier = (2.0, -3.0, 4.0)
    cuts = tuple(
        _held2_step6_cut(feed, multiplier, independent, packing)
        for independent, packing in (
            ((0.2, 0.1, 0.05), 0.2),
            ((0.4, 0.2, 0.10), 0.4),
            ((0.6, 0.05, 0.02), 0.6),
            (feed, 0.5),
        )
    )

    result = _held2_step6(
        feed,
        0.0,
        multiplier,
        cuts,
    )

    assert result["candidate_count"] == 4
    for actual, cut in zip(result["independent_modified_fractions"], cuts, strict=True):
        assert actual == pytest.approx(cut[0])


def test_held2_stage_ii_step6_uses_fixed_volume_gradient() -> None:
    feed = (0.4, 0.1, 0.04)
    multiplier = (2.0, -3.0, 4.0)
    cut = _held2_step6_cut(
        feed,
        multiplier,
        (0.3, 0.2, 0.05),
        0.3,
        q_gradient=(12.0, -13.0, 14.0, 0.0),
    )

    result = _held2_step6(feed, 0.0, multiplier, (cut,))

    assert result["candidate_count"] == 1


def test_held2_stage_ii_step6_applies_relative_scaling_and_excludes_lower_bound() -> None:
    feed = (0.4, 0.1, 0.04)
    multiplier = (1.0e8, 2.0, 3.0)
    cut = _held2_step6_cut(
        feed,
        multiplier,
        (0.3, 1.0e-10, 0.05),
        0.3,
        fixed_volume_gradient=(1.0e8 + 0.5, 100.0, 3.0 + 2.0e-8),
    )

    result = _held2_step6(feed, 0.0, multiplier, (cut,))

    assert result["candidate_count"] == 1


def test_held2_stage_ii_step6_clusters_only_eta_and_composition_duplicates() -> None:
    feed = (0.4, 0.1, 0.04)
    multiplier = (2.0, -3.0, 4.0)
    cuts = tuple(
        _held2_step6_cut(feed, multiplier, independent, packing)
        for independent, packing in (
            ((0.3, 0.1, 0.05), 0.3),
            ((0.3006, 0.1006, 0.0506), 0.3004),
            ((0.31, 0.1002, 0.0494), 0.3004),
        )
    )

    result = _held2_step6(feed, 0.0, multiplier, cuts)

    assert result["candidate_count"] == 2
    for actual, cut in zip(
        result["independent_modified_fractions"], (cuts[0], cuts[2]), strict=True
    ):
        assert actual == pytest.approx(cut[0])


def test_held2_stage_ii_step6_keeps_exact_c1_c29_replay_ineligible() -> None:
    feed = (0.6462224836985831, 0.04995832720498113, 0.03892729322939648)
    multiplier = (1.357695519654953, -0.196864482419187, -141.8379142414196)
    model = _khudaida_model()
    capsule = epcsaft.native_sdk(model)
    retained = (
        (
            (0.91049726787718377, 3.0460884822603739e-05, 0.088504042098029034),
            -0.7857111519928206,
            1.747197390420499e-05,
        ),
        (
            (0.4406535569650898, 0.08882925535454135, 0.0004486805424638745),
            -0.8436444167759523,
            5.613949558296684e-05,
        ),
    )
    cuts = []
    for independent, q, volume in retained:
        fixed_volume = _equilibrium._held2_adapter(
            capsule,
            293.15,
            100_000.0,
            independent,
            math.log(volume),
            KHUDAIDA_FINGERPRINT,
        )
        q_chart = _equilibrium._held2_adapter(
            capsule,
            293.15,
            100_000.0,
            independent,
            q,
            KHUDAIDA_FINGERPRINT,
            "log_packing_phase",
        )
        cuts.append(
            (
                independent,
                q_chart["packing_fraction"],
                q_chart["volume"],
                q_chart["objective"],
                tuple(q_chart["gradient"]),
                tuple(fixed_volume["gradient"][:3]),
            )
        )

    result = _held2_step6(feed, -6.114583161912342, multiplier, tuple(cuts))

    assert result["candidate_count"] == 0


def test_held2_general_mp_stage_iii_exact_lagrangian_hessian() -> None:
    candidates = ((0.2, 1.0), (0.20000002, 1.0), (0.8, 1.0))
    center = (0.25, 0.2, 1.0, 0.25, 0.21, 1.0, 0.5, 0.795, 1.0)
    multipliers = (0.3, -0.2)
    direction = (0.01, -0.02, 0.03, -0.01, 0.015, -0.02, 0.0, 0.005, 0.01)

    def evaluate(values: tuple[float, ...]) -> dict[str, object]:
        return _equilibrium._held2_adapter(
            CHARGES,
            PHYSICAL_FEED,
            candidates,
            values,
            multipliers,
            "stage_iii_derivatives",
        )

    step = 2.0e-5
    lower = evaluate(
        tuple(value - step * delta for value, delta in zip(center, direction, strict=True))
    )
    result = evaluate(center)
    upper = evaluate(
        tuple(value + step * delta for value, delta in zip(center, direction, strict=True))
    )

    objective_directional = (upper["objective"] - lower["objective"]) / (2.0 * step)
    assert objective_directional == pytest.approx(
        sum(
            value * delta
            for value, delta in zip(result["objective_gradient"], direction, strict=True)
        ),
        rel=2.0e-9,
        abs=2.0e-10,
    )
    size = len(center)
    for row in range(size):
        numerical = (upper["lagrangian_gradient"][row] - lower["lagrangian_gradient"][row]) / (
            2.0 * step
        )
        analytic = sum(
            result["lagrangian_hessian"][size * row + column] * direction[column]
            for column in range(size)
        )
        assert numerical == pytest.approx(analytic, rel=3.0e-8, abs=3.0e-9)
    assert tuple(result["lagrangian_hessian"]) == pytest.approx(
        tuple(
            result["lagrangian_hessian"][size * column + row]
            for row in range(size)
            for column in range(size)
        ),
        abs=2.0e-14,
    )


def test_held2_general_mp_stage_iii_simplex_schema_and_multiplier_mapping() -> None:
    feed = tuple(value / sum(KHUDAIDA_AMOUNTS) for value in KHUDAIDA_AMOUNTS)
    candidates = (
        (0.6588026822616131, 1.108849693428444e-5, 0.34115777594234437, -0.8047587326883016),
        (0.9104972678771838, 3.046088482260374e-5, 0.08850404209802903, -0.7857111519928206),
        (0.9784651953324194, 0.01956780416816369, 2.420187102135475e-5, -0.7872645327604558),
    )
    fraction = 1.0 / len(candidates)
    variables = tuple(value for candidate in candidates for value in (fraction, *candidate))
    equality_multipliers = (0.3, -0.2, 0.4, -0.1)
    zero_simplex = (*equality_multipliers, 0.0, 0.0, 0.0)
    simplex_multipliers = (1.1, -0.7, 0.4)

    baseline = _equilibrium._held2_adapter(
        KHUDAIDA_CHARGES,
        feed,
        candidates,
        variables,
        zero_simplex,
        "stage_iii_general_schema",
    )
    result = _equilibrium._held2_adapter(
        KHUDAIDA_CHARGES,
        feed,
        candidates,
        variables,
        (*equality_multipliers, *simplex_multipliers),
        "stage_iii_general_schema",
    )

    assert result["variable_count"] == 15
    assert result["constraint_count"] == 7
    assert result["jacobian_nonzero_count"] == 30
    assert result["constraint_lower_bounds"] == [0.0] * 4 + [-2.0e19] * 3
    assert result["constraint_upper_bounds"] == pytest.approx(
        [0.0] * 4 + [1.0 - 1.0e-10] * 3,
        abs=0.0,
    )
    assert result["constraints"][-3:] == pytest.approx(
        [sum(candidate[:-1]) for candidate in candidates],
        abs=2.0e-15,
    )

    expected_rows = [0, 0, 0]
    expected_columns = [0, 5, 10]
    expected_values = [1.0, 1.0, 1.0]
    for coordinate in range(3):
        for phase, candidate in enumerate(candidates):
            offset = 5 * phase
            expected_rows.extend((coordinate + 1, coordinate + 1))
            expected_columns.extend((offset, offset + coordinate + 1))
            expected_values.extend((candidate[coordinate], fraction))
    for phase in range(3):
        for coordinate in range(3):
            expected_rows.append(4 + phase)
            expected_columns.append(5 * phase + coordinate + 1)
            expected_values.append(1.0)
    assert result["jacobian_rows"] == expected_rows
    assert result["jacobian_columns"] == expected_columns
    assert result["jacobian_values"] == pytest.approx(expected_values, abs=0.0)

    gradient_difference = tuple(
        value - base
        for value, base in zip(
            result["lagrangian_gradient"],
            baseline["lagrangian_gradient"],
            strict=True,
        )
    )
    expected_difference = [0.0] * len(variables)
    for phase, multiplier in enumerate(simplex_multipliers):
        for coordinate in range(3):
            expected_difference[5 * phase + coordinate + 1] = multiplier
    assert gradient_difference == pytest.approx(expected_difference, abs=2.0e-15)
    assert result["lagrangian_hessian"] == pytest.approx(baseline["lagrangian_hessian"], abs=0.0)


def test_held2_general_mp_stage_iii_rejects_invalid_trial_before_provider() -> None:
    feed = tuple(value / sum(KHUDAIDA_AMOUNTS) for value in KHUDAIDA_AMOUNTS)
    candidates = (
        (0.6588026822616131, 1.108849693428444e-5, 0.34115777594234437, -0.8047587326883016),
        (0.9104972678771838, 3.046088482260374e-5, 0.08850404209802903, -0.7857111519928206),
    )
    invalid_trial = (
        0.5,
        0.6583790364854278,
        0.0010010884969342845,
        0.3407281048600504,
        -0.8053132473085454,
        0.5,
        *candidates[1],
    )

    result = _equilibrium._held2_adapter(
        KHUDAIDA_CHARGES,
        feed,
        candidates,
        invalid_trial,
        (0.0,) * 6,
        "stage_iii_invalid_trial",
    )

    assert result["objective_accepted"] is False
    assert result["scientific_evaluator_call_count"] == 0
    assert result["recoverable_domain_rejection_count"] == 1
    assert result["last_domain_rejection"]
    assert result["fatal_callback_error"] == ""


def test_held2_general_mp_stage_iii_refines_then_merges_duplicate_phases() -> None:
    result = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        ((0.2, 1.0), (0.20000002, 1.0), (0.8, 1.0)),
        "stage_iii",
    )

    assert result["profile"] == "perdomo-held2-stage-iii-manufactured-v1"
    assert result["solver_status"] in {"solve_succeeded", "solved_to_acceptable_level"}
    assert result["numerical_status"] == "converged"
    assert result["physical_status"] == "accepted"
    assert result["globality_certificate"] == "not_guaranteed"
    assert result["feedback"] == "none"
    assert result["input_candidate_count"] == 3
    assert result["retired_duplicate_count"] == 1
    assert len(result["phases"]) == 2
    assert sorted(phase["modified_fractions"][1] for phase in result["phases"]) == pytest.approx(
        [0.2, 0.8], abs=2.0e-7
    )
    assert sorted(phase["phase_fraction"] for phase in result["phases"]) == pytest.approx(
        [0.5, 0.5], abs=2.0e-7
    )
    assert result["modified_balance_inf_norm"] < 1.0e-9
    assert result["ordinary_balance_inf_norm"] < 1.0e-9
    assert result["phase_charge_inf_norm"] < 1.0e-12
    assert result["pressure_stationarity_inf_norm"] < 1.0e-9
    assert result["modified_potential_mixed_gap"] < 1.0e-9
    assert result["certified_modified_potential_count"] == 2
    assert result["trace_component_count"] == 0
    assert result["trace_refinement_status"] == "not_required"
    assert result["minimum_phase_distance"] > 1.0e-3
    assert result["kkt_stationarity_inf_norm"] < 1.0e-7
    assert abs(result["enumeration_objective_gap"]) < 1.0e-9


def test_held2_general_mp_stage_iii_returns_infeasible_set_to_stage_ii() -> None:
    result = _equilibrium._held2_adapter(
        CHARGES,
        PHYSICAL_FEED,
        ((0.1, 1.0), (0.2, 1.0), (0.3, 1.0)),
        "stage_iii",
    )

    assert result["numerical_status"] == "not_adjudicated"
    assert result["physical_status"] == "not_adjudicated"
    assert result["feedback"] == "return_to_stage_ii"
    assert result["failure_reason"] == "candidate_set_does_not_bracket_feed"
    assert result["globality_certificate"] == "not_guaranteed"
