from __future__ import annotations

import copy
import ctypes
import json
import math
from collections.abc import Callable
from pathlib import Path

import epcsaft
import pytest

from epcsaft_equilibrium import _equilibrium


class _StandardReferenceResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("temperature_k", ctypes.c_double),
        ("pressure_pa", ctypes.c_double),
        ("formula_unit_infinite_dilution_log_fugacity", ctypes.c_double),
        ("pure_solvent_log_fugacity_coefficient", ctypes.c_double),
        ("solvent_molar_mass_kg_per_mol", ctypes.c_double),
        ("reference_molality_mol_per_kg", ctypes.c_double),
        ("reference_convergence_error", ctypes.c_double),
        ("pure_solvent_molar_volume_m3_per_mol", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("helmholtz_basis_id", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class _StandardReferenceSdk(ctypes.Structure):
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
        ("ion_solvation_born_result_size", ctypes.c_size_t),
        ("evaluate_ion_solvation_born", ctypes.c_void_p),
        ("helmholtz_basis_id", ctypes.c_char_p),
        ("reference_amount_mol", ctypes.c_double),
        ("reference_number_density_mol_per_m3", ctypes.c_double),
        ("source_pressure_min_pa", ctypes.c_double),
        ("source_pressure_max_pa", ctypes.c_double),
        ("electrolyte_standard_reference_result_size", ctypes.c_size_t),
        ("evaluate_electrolyte_standard_reference", ctypes.c_void_p),
    )


_STANDARD_REFERENCE_CALLBACK = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(_StandardReferenceResult),
)


def _held_water_ionization_model() -> epcsaft.EPCSAFT:
    parameters = epcsaft.ParameterBundle.from_catalog(
        "held-2008-water-self-ionization", version=1
    ).select(("water", "hydronium-cation", "hydroxide-anion"))
    return epcsaft.EPCSAFT(parameters)


def _water_ionization_reference_record() -> dict[str, object]:
    path = (
        Path(__file__).parents[1]
        / "data/reference/water_self_ionization_iapws_r11_24.yaml"
    )
    return json.loads(path.read_text())


def _manufactured_water_standard_reference() -> dict[str, object]:
    return {
        "formula_unit_log_fugacity": -1.75,
        "pure_solvent_log_fugacity": -0.125,
        "solvent_molar_mass_kg_per_mol": 0.01801528,
        "reference_molality_mol_per_kg": 1.0e-12,
        "convergence_error": 2.0e-6,
        "pure_solvent_molar_volume_m3_per_mol": 1.8e-5,
        "basis_id": "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3",
        "parameter_fingerprint": "sha256:manufactured-water-reference",
        "component_ids": ["water", "hydronium-cation", "hydroxide-anion"],
        "charges": [0, 1, -1],
        "temperature_k": 298.15,
        "pressure_pa": 100_000.0,
    }


def _water_self_ionization_spec(
    model: epcsaft.EPCSAFT, *, feed_scale: float = 1.0
) -> dict[str, object]:
    return {
        "species_ids": ("water", "hydronium-cation", "hydroxide-anion"),
        "charges": (0, 1, -1),
        "provider_component_ids": (
            "water",
            "hydronium-cation",
            "hydroxide-anion",
        ),
        "provider_charges": (0, 1, -1),
        "provider_fingerprint": model.parameter_fingerprint,
        "expected_provider_component_ids": (
            "water",
            "hydronium-cation",
            "hydroxide-anion",
        ),
        "expected_provider_charges": (0, 1, -1),
        "expected_provider_fingerprint": model.parameter_fingerprint,
        "balance_matrix": ((2.0, 3.0, 1.0), (1.0, 1.0, 1.0)),
        "declared_balance_rank": 2,
        "reaction_matrix": ((-2.0, 1.0, 1.0),),
        "feed_amounts": (feed_scale, 0.0, 0.0),
        "temperature_k": 298.15,
        "pressure_pa": 100_000.0,
        "complete_closed_system": True,
        "water_self_ionization_reference": _water_ionization_reference_record(),
    }


def _iapws_p_kw(record: dict[str, object], temperature_k: float, density: float) -> float:
    correlation = record["correlation"]
    assert isinstance(correlation, dict)
    alpha = correlation["alpha"]
    beta = correlation["beta"]
    ideal = correlation["ideal_gas_p_kw_coefficients"]
    assert isinstance(alpha, list)
    assert isinstance(beta, list)
    assert isinstance(ideal, list)
    z_value = density * math.exp(
        alpha[0]
        + alpha[1] / temperature_k
        + alpha[2] * density ** (2.0 / 3.0) / temperature_k**2
    )
    p_kw_ideal = (
        ideal[0]
        + ideal[1] / temperature_k
        + ideal[2] / temperature_k**2
        + ideal[3] / temperature_k**3
    )
    coordination = correlation["ion_coordination_number"]
    molar_mass = correlation["water_molar_mass_g_per_mol"]
    return (
        -2.0
        * coordination
        * (
            math.log10(1.0 + z_value)
            - z_value / (z_value + 1.0) * density
            * (beta[0] + beta[1] / temperature_k + beta[2] * density)
        )
        + p_kw_ideal
        + 2.0 * math.log10(molar_mass / 1000.0)
    )


def _copied_sdk_capsule(
    model: epcsaft.EPCSAFT,
) -> tuple[object, _StandardReferenceSdk, tuple[object, ...]]:
    provider_capsule = epcsaft.native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    source = get_pointer(provider_capsule, b"epcsaft.native_sdk.v1")
    assert source
    table = _StandardReferenceSdk()
    ctypes.memmove(ctypes.addressof(table), source, ctypes.sizeof(table))
    name = ctypes.create_string_buffer(b"epcsaft.native_sdk.v1")
    new_capsule = ctypes.pythonapi.PyCapsule_New
    new_capsule.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p)
    new_capsule.restype = ctypes.py_object
    capsule = new_capsule(ctypes.addressof(table), name, None)
    return capsule, table, (provider_capsule, name)


def _base_system() -> dict[str, object]:
    temperature_k = 350.0
    pressure_pa = 200_000.0
    return {
        "species_ids": ("A", "B"),
        "charges": (0, 0),
        "provider_component_ids": ("A", "B"),
        "provider_charges": (0, 0),
        "provider_fingerprint": "sha256:manufactured",
        "expected_provider_component_ids": ("A", "B"),
        "expected_provider_charges": (0, 0),
        "expected_provider_fingerprint": "sha256:manufactured",
        "balance_matrix": ((1.0, 1.0),),
        "declared_balance_rank": 1,
        "reaction_matrix": ((-1.0, 1.0),),
        "feed_amounts": (1.0, 0.0),
        "ln_k": (math.log(4.0),),
        "equilibrium_constant_records": (
            {
                "source_id": "manufactured:A-to-B",
                "reference_id": "provider-helmholtz-coordinate-basis",
                "dimensionless": True,
                "temperature_k": temperature_k,
                "pressure_pa": pressure_pa,
                "pressure_binding": "fixed",
            },
        ),
        "temperature_k": temperature_k,
        "pressure_pa": pressure_pa,
        "complete_closed_system": True,
    }


def test_reaction_compiler_reconstructs_minimum_norm_reference() -> None:
    compiled = _equilibrium._chemical_compile_system(_base_system())

    expected = 0.5 * math.log(4.0)
    assert compiled["species_ids"] == ["A", "B"]
    assert compiled["balance_rank"] == 1
    assert compiled["reaction_rank"] == 1
    assert compiled["g_ref"] == pytest.approx((expected, -expected), abs=2.0e-15)
    assert compiled["reference_reconstruction_inf_norm"] <= 2.0e-15
    assert compiled["conservation_reaction_inf_norm"] == 0.0
    assert compiled["charge_reaction_inf_norm"] == 0.0


def _change(path: str, value: object) -> Callable[[dict[str, object]], None]:
    def apply(spec: dict[str, object]) -> None:
        if "." not in path:
            spec[path] = value
            return
        outer, inner = path.split(".", maxsplit=1)
        records = [dict(record) for record in spec[outer]]  # type: ignore[arg-type]
        records[0][inner] = value
        spec[outer] = tuple(records)

    return apply


@pytest.mark.parametrize(
    ("mutate", "message"),
    (
        (_change("provider_component_ids", ("B", "A")), "component order"),
        (_change("provider_charges", (1, -1)), "charges"),
        (_change("provider_fingerprint", "sha256:wrong"), "fingerprint"),
        (_change("declared_balance_rank", 2), "declared balance rank"),
        (_change("balance_matrix", ((1.0, 1.0), (2.0, 2.0))), "balance matrix rank"),
        (_change("reaction_matrix", ((-1.0, 2.0),)), "conserve"),
        (
            _change("reaction_matrix", ((-1.0, 1.0), (-2.0, 2.0))),
            "reaction matrix rank",
        ),
        (_change("complete_closed_system", False), "complete closed system"),
        (_change("equilibrium_constant_records.dimensionless", False), "dimensionless"),
        (_change("equilibrium_constant_records.source_id", ""), "source identity"),
        (_change("equilibrium_constant_records.temperature_k", 351.0), "temperature"),
        (_change("equilibrium_constant_records.pressure_pa", 201_000.0), "pressure"),
    ),
)
def test_reaction_compiler_rejects_inconsistent_contracts(
    mutate: Callable[[dict[str, object]], None], message: str
) -> None:
    spec = copy.deepcopy(_base_system())
    if message == "complete closed system":
        spec["reaction_matrix"] = ()
        spec["ln_k"] = ()
        spec["equilibrium_constant_records"] = ()
        spec["complete_closed_system"] = True
    else:
        mutate(spec)

    with pytest.raises(ValueError, match=message):
        _equilibrium._chemical_compile_system(spec)


def test_reaction_compiler_rejects_non_neutral_feed_and_charge_nonconservation() -> None:
    non_neutral = _base_system()
    non_neutral.update(
        {
            "charges": (1, -1),
            "provider_charges": (1, -1),
            "expected_provider_charges": (1, -1),
            "feed_amounts": (1.0, 0.0),
        }
    )
    with pytest.raises(ValueError, match="electroneutral"):
        _equilibrium._chemical_compile_system(non_neutral)

    nonconserving = copy.deepcopy(non_neutral)
    nonconserving["feed_amounts"] = (1.0, 1.0)
    with pytest.raises(ValueError, match="charge"):
        _equilibrium._chemical_compile_system(nonconserving)


def test_reaction_compiler_rejects_nonfinite_temperature_only_pressure_metadata() -> None:
    spec = _base_system()
    record = dict(spec["equilibrium_constant_records"][0])  # type: ignore[index]
    record["pressure_binding"] = "temperature_only"
    record["pressure_pa"] = math.nan
    spec["equilibrium_constant_records"] = (record,)

    with pytest.raises(ValueError, match="pressure must be finite"):
        _equilibrium._chemical_compile_system(spec)


def _amount_chart(
    charges: tuple[int, ...], coordinates: tuple[float, ...], trace_floor: float = 1.0e-12
) -> dict[str, object]:
    return _equilibrium._chemical_amount_chart(charges, coordinates, trace_floor)


def test_amount_chart_maps_neutral_logs_and_general_ionic_shares_exactly() -> None:
    neutral = _amount_chart((0, 0), (math.log(2.0), math.log(3.0)))
    assert neutral["amounts"] == pytest.approx((2.0, 3.0), abs=2.0e-15)
    assert neutral["charge_residual"] == 0.0

    charges = (0, 2, 1, -1, -2, 0)
    coordinates = (math.log(2.0), math.log(1.0 / 3.0), 0.0, math.log(3.0), math.log(4.0))
    ionic = _amount_chart(charges, coordinates)
    assert ionic["amounts"] == pytest.approx((3.0, 0.25, 1.5, 1.0, 0.5, 4.0), abs=2.0e-15)
    assert ionic["coordinate_count"] == len(charges) - 1
    assert ionic["charge_residual"] == pytest.approx(0.0, abs=2.0e-15)
    assert ionic["trace_status"] == "interior"

    inverse = _equilibrium._chemical_amount_chart_inverse(charges, ionic["amounts"])
    assert inverse == pytest.approx(coordinates, abs=3.0e-15)


def test_amount_chart_is_permutation_and_amount_scale_equivariant() -> None:
    charges = (0, 2, 1, -1, -2, 0)
    amounts = (3.0, 0.25, 1.5, 1.0, 0.5, 4.0)
    permutation = (5, 3, 1, 0, 4, 2)
    permuted_charges = tuple(charges[index] for index in permutation)
    permuted_amounts = tuple(amounts[index] for index in permutation)

    coordinates = _equilibrium._chemical_amount_chart_inverse(permuted_charges, permuted_amounts)
    result = _amount_chart(permuted_charges, tuple(coordinates))
    assert result["amounts"] == pytest.approx(permuted_amounts, abs=3.0e-15)

    scale = 7.0
    scaled = tuple(scale * value for value in permuted_amounts)
    scaled_coordinates = _equilibrium._chemical_amount_chart_inverse(permuted_charges, scaled)
    scaled_result = _amount_chart(permuted_charges, tuple(scaled_coordinates))
    assert scaled_result["amounts"] == pytest.approx(scaled, rel=2.0e-15)


def test_amount_chart_has_exact_directional_first_and_second_derivatives() -> None:
    charges = (0, 2, 1, -1, -2, 0)
    center = (0.2, -0.7, 0.4, -0.3, 0.8)
    direction = (0.3, -0.2, 0.5, 0.1, -0.4)
    step = 2.0e-5

    lower = _amount_chart(
        charges,
        tuple(value - step * delta for value, delta in zip(center, direction, strict=True)),
    )
    result = _amount_chart(charges, center)
    upper = _amount_chart(
        charges,
        tuple(value + step * delta for value, delta in zip(center, direction, strict=True)),
    )

    dimension = len(center)
    jacobian = result["jacobian"]
    hessians = result["amount_hessians"]
    for species in range(len(charges)):
        first_directional = (upper["amounts"][species] - lower["amounts"][species]) / (
            2.0 * step
        )
        exact_first = sum(
            jacobian[species * dimension + column] * direction[column]
            for column in range(dimension)
        )
        assert first_directional == pytest.approx(exact_first, rel=3.0e-10, abs=3.0e-11)

        for row in range(dimension):
            jacobian_directional = (
                upper["jacobian"][species * dimension + row]
                - lower["jacobian"][species * dimension + row]
            ) / (2.0 * step)
            exact_second = sum(
                hessians[
                    species * dimension * dimension + row * dimension + column
                ]
                * direction[column]
                for column in range(dimension)
            )
            assert jacobian_directional == pytest.approx(
                exact_second, rel=2.0e-9, abs=5.0e-11
            )


def test_amount_chart_classifies_trace_floor_without_zeroing_species() -> None:
    trace_floor = 1.0e-10
    result = _amount_chart((0, 0), (math.log(0.5 * trace_floor), 0.0), trace_floor)

    assert result["trace_status"] == "at_or_below_floor"
    assert result["minimum_amount"] == pytest.approx(0.5 * trace_floor, rel=2.0e-15)
    assert all(value > 0.0 for value in result["amounts"])


@pytest.mark.parametrize(
    ("charges", "coordinates", "message"),
    (
        ((1, 0), (0.0,), "both cations and anions"),
        ((1, -1), (0.0, 0.0), "coordinate count"),
        ((1, -1), (math.inf,), "finite"),
        ((1, 1, -1), (0.0, -1000.0), "strictly positive representable"),
    ),
)
def test_amount_chart_rejects_invalid_topology_or_coordinates(
    charges: tuple[int, ...], coordinates: tuple[float, ...], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        _amount_chart(charges, coordinates)


def _bind_record(spec: dict[str, object]) -> None:
    records = []
    for index in range(len(spec["ln_k"])):  # type: ignore[arg-type]
        records.append(
            {
                "source_id": f"manufactured:reaction-{index}",
                "reference_id": "provider-helmholtz-coordinate-basis",
                "dimensionless": True,
                "temperature_k": spec["temperature_k"],
                "pressure_pa": spec["pressure_pa"],
                "pressure_binding": "fixed",
            }
        )
    spec["equilibrium_constant_records"] = tuple(records)


def _manufactured_solve(
    spec: dict[str, object], options: dict[str, object] | None = None
) -> dict[str, object]:
    return _equilibrium._chemical_solve_manufactured(spec, options or {})


@pytest.mark.parametrize(
    ("spec", "expected_amounts", "expected_volume"),
    (
        (_base_system(), (0.2, 0.8), None),
        (
            {
                **_base_system(),
                "balance_matrix": ((1.0, 2.0),),
                "reaction_matrix": ((-2.0, 1.0),),
                "feed_amounts": (2.0, 0.0),
                "ln_k": (math.log(0.75),),
                "temperature_k": 300.0,
                "pressure_pa": 8.31446261815324 * 300.0,
            },
            (1.0, 0.5),
            1.5,
        ),
    ),
)
def test_manufactured_ideal_reactions_match_independent_analytic_states(
    spec: dict[str, object],
    expected_amounts: tuple[float, ...],
    expected_volume: float | None,
) -> None:
    _bind_record(spec)
    result = _manufactured_solve(spec)

    pressure_over_rt = spec["pressure_pa"] / (8.31446261815324 * spec["temperature_k"])
    volume = expected_volume or sum(expected_amounts) / pressure_over_rt
    assert result["accepted"] is True
    assert result["profile"] == "manufactured_ideal_nonpredictive"
    assert result["amounts"] == pytest.approx(expected_amounts, rel=2.0e-8, abs=2.0e-10)
    assert result["volume_m3"] == pytest.approx(volume, rel=2.0e-8)
    assert result["solver_status"] == "solve_succeeded"
    assert result["artifact_input_status"] == "passed"
    assert result["numerical_status"] == "passed"
    assert result["physical_status"] == "passed"
    assert result["local_minimum_status"] == "passed"
    assert result["trace_status"] == "interior"
    assert result["predictive_status"] == "not_adjudicated"
    assert result["globality_certificate"] == "not_guaranteed"
    assert result["final_lambda"] == 1.0
    assert result["balance_inf_norm"] <= 1.0e-9
    assert result["charge_inf_norm"] <= 1.0e-12
    assert result["pressure_relative_residual"] <= 1.0e-8
    assert result["reaction_affinity_inf_norm"] <= 1.0e-7
    assert result["kkt_stationarity_inf_norm"] <= 1.0e-7
    assert result["complementarity_inf_norm"] <= 1.0e-7
    assert result["kkt_scope"] == "equality_kkt_on_strict_interior"
    assert len(result["kkt_residual"]) + len(result["kkt_jacobian"]) > 0


def test_manufactured_charged_reaction_uses_exact_electroneutral_chart() -> None:
    temperature_k = 300.0
    spec = {
        **_base_system(),
        "species_ids": ("A", "C+", "D-"),
        "charges": (0, 1, -1),
        "provider_component_ids": ("A", "C+", "D-"),
        "provider_charges": (0, 1, -1),
        "expected_provider_component_ids": ("A", "C+", "D-"),
        "expected_provider_charges": (0, 1, -1),
        "balance_matrix": ((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
        "declared_balance_rank": 2,
        "reaction_matrix": ((-1.0, 1.0, 1.0),),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (math.log(1.0 / 3.0),),
        "temperature_k": temperature_k,
        "pressure_pa": 8.31446261815324 * temperature_k,
    }
    _bind_record(spec)

    result = _manufactured_solve(spec)

    assert result["accepted"] is True
    assert result["amounts"] == pytest.approx((0.5, 0.5, 0.5), rel=3.0e-8)
    assert result["volume_m3"] == pytest.approx(1.5, rel=3.0e-8)
    assert result["charge_inf_norm"] <= 2.0e-15
    assert result["active_balance_constraint_count"] == 1


def test_manufactured_equilibrium_is_gauge_scale_and_reaction_basis_invariant() -> None:
    base = _base_system()
    base["feed_amounts"] = (1.0, 0.0)
    _bind_record(base)
    plain = _manufactured_solve(base)
    gauged = _manufactured_solve(base, {"gauge_coefficients": (3.25,)})
    assert gauged["amounts"] == pytest.approx(plain["amounts"], rel=2.0e-10)

    scaled = copy.deepcopy(base)
    scaled["feed_amounts"] = (7.0, 0.0)
    scaled_result = _manufactured_solve(scaled)
    assert scaled_result["amounts"] == pytest.approx(
        tuple(7.0 * value for value in plain["amounts"]), rel=3.0e-8
    )
    assert scaled_result["volume_m3"] == pytest.approx(7.0 * plain["volume_m3"], rel=3.0e-8)

    three = {
        **_base_system(),
        "species_ids": ("A", "B", "C"),
        "charges": (0, 0, 0),
        "provider_component_ids": ("A", "B", "C"),
        "provider_charges": (0, 0, 0),
        "expected_provider_component_ids": ("A", "B", "C"),
        "expected_provider_charges": (0, 0, 0),
        "balance_matrix": ((1.0, 1.0, 1.0),),
        "reaction_matrix": ((-1.0, 1.0, 0.0), (0.0, -1.0, 1.0)),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (math.log(2.0), math.log(3.0)),
        "complete_closed_system": True,
    }
    _bind_record(three)
    first_basis = _manufactured_solve(three)
    changed_basis = copy.deepcopy(three)
    changed_basis["reaction_matrix"] = ((-1.0, 0.0, 1.0), (0.0, -1.0, 1.0))
    changed_basis["ln_k"] = (math.log(6.0), math.log(3.0))
    _bind_record(changed_basis)
    second_basis = _manufactured_solve(changed_basis)
    assert first_basis["amounts"] == pytest.approx(second_basis["amounts"], rel=3.0e-9)


def test_max_min_initialization_fails_closed_without_strict_positive_state() -> None:
    result = _equilibrium._chemical_max_min_initialization(
        ((1.0, 0.0), (0.0, 1.0)),
        (1.0, 0.0),
        (0, 0),
        1.0e-10,
    )

    assert result["solver_status"] == "solve_succeeded"
    assert result["strict_positive_feasible"] is False
    assert result["max_min_amount"] <= 1.0e-10
    assert result["reason"] == "no_strict_positive_state_above_trace_floor"

    unbounded = _equilibrium._chemical_max_min_initialization(
        ((0.0, 0.0),),
        (0.0, 0.0),
        (0, 0),
        1.0e-10,
    )
    assert unbounded["strict_positive_feasible"] is False


def test_max_min_initialization_has_no_artificial_amount_cap() -> None:
    result = _equilibrium._chemical_max_min_initialization(
        ((1000.0, 1.0, 0.0), (1.0, 0.0, 1.0)),
        (1.0, 0.0, 0.0),
        (0, 0, 0),
        0.1,
    )

    assert result["strict_positive_feasible"] is True
    assert result["max_min_amount"] == pytest.approx(0.5, rel=2.0e-8)
    assert result["amounts"][1] > 100.0


def test_manufactured_reaction_solve_has_no_feed_scaled_amount_cap() -> None:
    target = (0.5, 500.0, 0.5)
    volume = sum(target)
    spec = {
        **_base_system(),
        "species_ids": ("A", "B", "C"),
        "charges": (0, 0, 0),
        "provider_component_ids": ("A", "B", "C"),
        "provider_charges": (0, 0, 0),
        "expected_provider_component_ids": ("A", "B", "C"),
        "expected_provider_charges": (0, 0, 0),
        "balance_matrix": ((1000.0, 1.0, 0.0), (1.0, 0.0, 1.0)),
        "declared_balance_rank": 2,
        "reaction_matrix": ((-1.0, 1000.0, 1.0),),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (1000.0 * math.log(target[1] / volume),),
        "temperature_k": 300.0,
        "pressure_pa": 8.31446261815324 * 300.0,
    }
    _bind_record(spec)

    result = _manufactured_solve(spec, {"trace_floor": 0.1})

    assert result["accepted"] is True
    assert result["amounts"] == pytest.approx(target, rel=2.0e-7, abs=2.0e-8)


def test_manufactured_nlp_has_exact_directional_gradient_jacobian_and_hessian() -> None:
    spec = _base_system()
    _bind_record(spec)
    pressure_over_rt = spec["pressure_pa"] / (8.31446261815324 * spec["temperature_k"])
    center = (math.log(0.3), math.log(0.7), math.log(1.0 / pressure_over_rt))
    direction = (0.2, -0.4, 0.3)
    multipliers = (0.37,)
    step = 2.0e-5

    def evaluate(variables: tuple[float, ...]) -> dict[str, object]:
        return _equilibrium._chemical_evaluate_manufactured_nlp(
            spec, variables, multipliers
        )

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
            result["objective_gradient"][index] * direction[index]
            for index in range(len(direction))
        ),
        rel=2.0e-9,
        abs=2.0e-10,
    )
    constraint_directional = (upper["constraints"][0] - lower["constraints"][0]) / (
        2.0 * step
    )
    assert constraint_directional == pytest.approx(
        sum(
            result["constraint_jacobian"][index] * direction[index]
            for index in range(len(direction))
        ),
        rel=2.0e-10,
        abs=2.0e-11,
    )
    dimension = len(direction)
    for row in range(dimension):
        gradient_directional = (
            upper["lagrangian_gradient"][row] - lower["lagrangian_gradient"][row]
        ) / (2.0 * step)
        hessian_directional = sum(
            result["lagrangian_hessian"][row * dimension + column] * direction[column]
            for column in range(dimension)
        )
        assert gradient_directional == pytest.approx(
            hessian_directional, rel=3.0e-9, abs=3.0e-10
        )


def test_manufactured_solver_rejects_indeterminate_and_false_success_terminals() -> None:
    spec = _base_system()
    _bind_record(spec)
    indeterminate = _manufactured_solve(spec, {"test_max_iterations": 0})
    assert indeterminate["solver_status"] == "maximum_iterations_exceeded"
    assert indeterminate["accepted"] is False
    assert indeterminate["numerical_status"] == "failed"
    assert indeterminate["final_lambda"] is None

    trace = copy.deepcopy(spec)
    trace["ln_k"] = (math.log(1.0e-10),)
    _bind_record(trace)
    false_success = _manufactured_solve(trace, {"trace_floor": 1.0e-8})
    assert false_success["accepted"] is False
    assert false_success["physical_status"] == "failed"
    assert false_success["trace_status"] == "at_or_below_floor"
    assert false_success["globality_certificate"] == "not_guaranteed"


def test_manufactured_charged_solution_is_species_order_invariant() -> None:
    temperature_k = 300.0
    original = {
        **_base_system(),
        "species_ids": ("A", "C+", "D-"),
        "charges": (0, 1, -1),
        "provider_component_ids": ("A", "C+", "D-"),
        "provider_charges": (0, 1, -1),
        "expected_provider_component_ids": ("A", "C+", "D-"),
        "expected_provider_charges": (0, 1, -1),
        "balance_matrix": ((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
        "declared_balance_rank": 2,
        "reaction_matrix": ((-1.0, 1.0, 1.0),),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (math.log(1.0 / 3.0),),
        "temperature_k": temperature_k,
        "pressure_pa": 8.31446261815324 * temperature_k,
    }
    _bind_record(original)
    baseline = _manufactured_solve(original)
    permutation = (2, 0, 1)
    permuted = copy.deepcopy(original)
    for field in (
        "species_ids",
        "charges",
        "provider_component_ids",
        "provider_charges",
        "expected_provider_component_ids",
        "expected_provider_charges",
        "feed_amounts",
    ):
        values = original[field]
        permuted[field] = tuple(values[index] for index in permutation)
    permuted["balance_matrix"] = tuple(
        tuple(row[index] for index in permutation) for row in original["balance_matrix"]
    )
    permuted["reaction_matrix"] = tuple(
        tuple(row[index] for index in permutation) for row in original["reaction_matrix"]
    )
    _bind_record(permuted)
    reordered = _manufactured_solve(permuted)
    inverse = tuple(permutation.index(index) for index in range(len(permutation)))
    assert tuple(reordered["amounts"][index] for index in inverse) == pytest.approx(
        baseline["amounts"], rel=4.0e-8
    )


def test_installed_provider_standard_reference_tail_is_consumed_once() -> None:
    model = _held_water_ionization_model()
    result = _equilibrium._chemical_evaluate_provider_standard_reference(
        epcsaft.native_sdk(model),
        298.15,
        100_000.0,
        model.parameter_fingerprint,
    )
    public = model.reference_state(
        temperature=298.15 * epcsaft.unit_registry.kelvin,
        pressure=100_000.0 * epcsaft.unit_registry.pascal,
        solvent_mole_fractions=(1.0, 0.0, 0.0),
        phase="liquid",
    )
    pure_water = model.evaluate(
        temperature=298.15 * epcsaft.unit_registry.kelvin,
        pressure=100_000.0 * epcsaft.unit_registry.pascal,
        mole_fractions=(1.0, 0.0, 0.0),
        phase="liquid",
    )

    assert result["formula_unit_log_fugacity"] == pytest.approx(
        math.fsum(public.infinite_dilution_log_fugacity), abs=2.0e-12
    )
    assert result["pure_solvent_log_fugacity"] == pytest.approx(
        public._pure_water_log_fugacity, abs=2.0e-12
    )
    assert result["solvent_molar_mass_kg_per_mol"] == pytest.approx(0.01801528)
    assert result["reference_molality_mol_per_kg"] == 1.0e-12
    assert 0.0 <= result["convergence_error"] <= 5.0e-5
    assert result["pure_solvent_molar_volume_m3_per_mol"] == pytest.approx(
        1.0 / pure_water.molar_density.to("mole / meter**3").magnitude,
        rel=2.0e-12,
    )
    assert result["basis_id"] == (
        "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3"
    )
    assert result["parameter_fingerprint"] == model.parameter_fingerprint
    assert result["component_ids"] == ["water", "hydronium-cation", "hydroxide-anion"]
    assert result["charges"] == [0, 1, -1]


def test_retained_iapws_record_reproduces_source_equations() -> None:
    record = _water_ionization_reference_record()
    state = record["state"]
    values = record["values"]
    assert isinstance(state, dict)
    assert isinstance(values, dict)

    assert _iapws_p_kw(record, 300.0, 1.0) == pytest.approx(13.906672, abs=5.0e-7)
    retained_p_kw = _iapws_p_kw(
        record,
        state["temperature_k"],
        state["density_kg_per_m3"] / 1000.0,
    )
    assert retained_p_kw == pytest.approx(values["p_kw"], abs=5.0e-13)
    assert -math.log(10.0) * retained_p_kw == pytest.approx(
        values["ln_kw"], abs=1.0e-13
    )


def test_iapws_mixed_standard_state_transform_uses_thermodynamic_standard_molality() -> None:
    record = _water_ionization_reference_record()
    reference = _manufactured_water_standard_reference()

    transformed = _equilibrium._chemical_transform_water_self_ionization_standard_state(
        record, reference
    )

    standard_state = record["standard_state"]
    values = record["values"]
    assert isinstance(standard_state, dict)
    assert isinstance(values, dict)
    expected_delta = (
        2.0
        * math.log(
            reference["solvent_molar_mass_kg_per_mol"]
            * standard_state["standard_molality_mol_per_kg"]
        )
        + reference["formula_unit_log_fugacity"]
        - 2.0 * reference["pure_solvent_log_fugacity"]
    )
    wrong_terminal_delta = (
        2.0
        * math.log(
            reference["solvent_molar_mass_kg_per_mol"]
            * reference["reference_molality_mol_per_kg"]
        )
        + reference["formula_unit_log_fugacity"]
        - 2.0 * reference["pure_solvent_log_fugacity"]
    )
    assert transformed["delta_standard_offset"] == pytest.approx(
        expected_delta, abs=2.0e-15
    )
    assert transformed["delta_standard_offset"] != pytest.approx(wrong_terminal_delta)
    assert transformed["ln_k_provider_basis"] == pytest.approx(
        values["ln_kw"] + expected_delta, abs=5.0e-15
    )
    assert transformed["transformation_residual"] <= 2.0e-15
    assert transformed["provider_basis_id"] == reference["basis_id"]
    assert transformed["source_id"] == "iapws-r11-24"


def test_installed_provider_reference_transforms_iapws_ln_kw_in_exact_basis() -> None:
    model = _held_water_ionization_model()
    provider_reference = _equilibrium._chemical_evaluate_provider_standard_reference(
        epcsaft.native_sdk(model),
        298.15,
        100_000.0,
        model.parameter_fingerprint,
    )
    record = _water_ionization_reference_record()

    transformed = _equilibrium._chemical_transform_water_self_ionization_standard_state(
        record, provider_reference
    )

    standard_state = record["standard_state"]
    values = record["values"]
    assert isinstance(standard_state, dict)
    assert isinstance(values, dict)
    expected_delta = (
        2.0
        * math.log(
            provider_reference["solvent_molar_mass_kg_per_mol"]
            * standard_state["standard_molality_mol_per_kg"]
        )
        + provider_reference["formula_unit_log_fugacity"]
        - 2.0 * provider_reference["pure_solvent_log_fugacity"]
    )
    assert transformed["delta_standard_offset"] == pytest.approx(
        expected_delta, abs=2.0e-15
    )
    assert transformed["ln_k_provider_basis"] == pytest.approx(
        values["ln_kw"] + expected_delta, abs=5.0e-15
    )
    assert transformed["provider_basis_id"] == provider_reference["basis_id"]
    assert transformed["parameter_fingerprint"] == model.parameter_fingerprint


@pytest.mark.parametrize(
    ("section", "field", "value", "message"),
    (
        ("standard_state", "solvent_molar_mass_unit", "g/mol", "molar-mass unit"),
        (
            "standard_state",
            "transformation_sign",
            "lnK_provider=lnKw_mixed-delta_standard_offset",
            "transformation sign",
        ),
        ("standard_state", "standard_molality_mol_per_kg", 1.0e-12, "standard molality"),
        ("provider_binding", "basis_id", "wrong-basis", "basis identity"),
        ("provider_binding", "reaction_stoichiometry", (-1.0, 1.0, 1.0), "1:1"),
        ("state", "phase", "vapor", "state identity"),
        ("state", "pressure_binding", "temperature_only", "pressure binding"),
        ("source", "id", "unidentified", "source identity"),
        ("values", "ln_kw", 32.2231929374538, "pKw and lnKw"),
    ),
)
def test_iapws_mixed_standard_state_transform_rejects_record_mutations(
    section: str, field: str, value: object, message: str
) -> None:
    record = _water_ionization_reference_record()
    nested = record[section]
    assert isinstance(nested, dict)
    nested[field] = value

    with pytest.raises(ValueError, match=message):
        _equilibrium._chemical_transform_water_self_ionization_standard_state(
            record, _manufactured_water_standard_reference()
        )


@pytest.mark.parametrize(
    ("field", "value", "message"),
    (
        ("basis_id", "wrong-basis", "basis identity"),
        ("component_ids", ("hydronium-cation", "water", "hydroxide-anion"), "order"),
        ("charges", (0, -1, 1), "charges"),
        ("temperature_k", 299.0, "temperature binding"),
        ("pressure_pa", 101_325.0, "pressure binding"),
        ("solvent_molar_mass_kg_per_mol", 18.01528, "molar mass"),
        ("reference_molality_mol_per_kg", 0.0, "terminal molality"),
        ("convergence_error", 1.0e-3, "convergence"),
        ("pure_solvent_molar_volume_m3_per_mol", 0.0, "molar volume"),
        ("parameter_fingerprint", "", "fingerprint"),
    ),
)
def test_iapws_mixed_standard_state_transform_rejects_provider_reference_mutations(
    field: str, value: object, message: str
) -> None:
    reference = _manufactured_water_standard_reference()
    reference[field] = value

    with pytest.raises(ValueError, match=message):
        _equilibrium._chemical_transform_water_self_ionization_standard_state(
            _water_ionization_reference_record(), reference
        )


def test_source_complete_water_self_ionization_solves_and_recomputes_iapws_activity() -> None:
    model = _held_water_ionization_model()
    capsule = epcsaft.native_sdk(model)
    spec = _water_self_ionization_spec(model)

    result = _equilibrium._chemical_solve_provider_water_self_ionization(
        capsule,
        spec,
        {"packing_fraction_bounds": (1.0e-6, 0.74), "trace_floor": 1.0e-12},
    )

    assert result["profile"] == "installed_provider_source_complete_consistency"
    assert result["accepted"] is True
    assert result["solver_status"] == "solve_succeeded"
    assert result["artifact_input_status"] == "passed"
    assert result["numerical_status"] == "passed"
    assert result["physical_status"] == "passed"
    assert result["provider_domain_status"] == "passed"
    assert result["local_minimum_status"] == "passed"
    assert result["trace_status"] == "interior"
    assert result["predictive_status"] == "not_adjudicated"
    assert result["globality_certificate"] == "not_guaranteed"
    assert result["final_lambda"] == 1.0
    assert result["balance_inf_norm"] <= 1.0e-9
    assert result["charge_inf_norm"] <= 1.0e-9
    assert result["pressure_relative_residual"] <= 1.0e-8
    assert result["reaction_affinity_inf_norm"] <= 1.0e-7
    assert result["kkt_stationarity_inf_norm"] <= 1.0e-7
    assert result["complementarity_inf_norm"] <= 1.0e-7
    assert result["reference_reconstruction_inf_norm"] <= 1.0e-12
    assert result["standard_state_transformation_residual"] <= 2.0e-15

    provider_reference = _equilibrium._chemical_evaluate_provider_standard_reference(
        capsule, 298.15, 100_000.0, model.parameter_fingerprint
    )
    phase = _equilibrium._chemical_evaluate_provider_block(
        capsule,
        298.15,
        result["amounts"],
        result["volume_m3"],
        model.parameter_fingerprint,
    )
    total_amount = math.fsum(result["amounts"])
    mole_fractions = tuple(amount / total_amount for amount in result["amounts"])
    log_fugacity = tuple(
        phase["gradient"][index]
        - math.log(mole_fractions[index] * 100_000.0 / (8.31446261815324 * 298.15))
        for index in range(3)
    )
    water_amount, hydronium_amount, hydroxide_amount = result["amounts"]
    water_mass_kg = water_amount * provider_reference["solvent_molar_mass_kg_per_mol"]
    log_hydronium_molality = math.log(hydronium_amount / water_mass_kg)
    log_hydroxide_molality = math.log(hydroxide_amount / water_mass_kg)
    ion_activity_sum = (
        log_hydronium_molality
        + log_hydroxide_molality
        + log_fugacity[1]
        + log_fugacity[2]
        - provider_reference["formula_unit_log_fugacity"]
    )
    log_water_activity = (
        math.log(mole_fractions[0])
        + log_fugacity[0]
        - provider_reference["pure_solvent_log_fugacity"]
    )
    reference_record = spec["water_self_ionization_reference"]
    assert isinstance(reference_record, dict)
    values = reference_record["values"]
    assert isinstance(values, dict)
    independently_recomputed = ion_activity_sum - 2.0 * log_water_activity
    assert independently_recomputed == pytest.approx(values["ln_kw"], abs=1.0e-8)


@pytest.mark.parametrize("trace_floor", (1.0e-18, 1.0e-12))
def test_source_complete_water_self_ionization_is_trace_floor_stable(
    trace_floor: float,
) -> None:
    model = _held_water_ionization_model()
    result = _equilibrium._chemical_solve_provider_water_self_ionization(
        epcsaft.native_sdk(model),
        _water_self_ionization_spec(model),
        {"packing_fraction_bounds": (1.0e-6, 0.74), "trace_floor": trace_floor},
    )

    assert result["accepted"] is True
    assert result["amounts"][1] == pytest.approx(1.8139519e-9, rel=2.0e-7)
    assert result["reaction_affinity_inf_norm"] <= 1.0e-7


def test_source_complete_water_self_ionization_is_scale_and_basis_invariant() -> None:
    model = _held_water_ionization_model()
    capsule = epcsaft.native_sdk(model)
    options = {"packing_fraction_bounds": (1.0e-6, 0.74), "trace_floor": 1.0e-12}
    baseline = _equilibrium._chemical_solve_provider_water_self_ionization(
        capsule, _water_self_ionization_spec(model), options
    )

    scaled_feed = _equilibrium._chemical_solve_provider_water_self_ionization(
        capsule, _water_self_ionization_spec(model, feed_scale=10.0), options
    )
    scaled_basis_spec = _water_self_ionization_spec(model)
    scaled_basis_spec["reaction_matrix"] = ((-4.0, 2.0, 2.0),)
    scaled_basis = _equilibrium._chemical_solve_provider_water_self_ionization(
        capsule, scaled_basis_spec, options
    )

    assert baseline["accepted"] is True
    assert scaled_feed["accepted"] is True
    assert scaled_basis["accepted"] is True
    assert scaled_basis["amounts"] == pytest.approx(baseline["amounts"], rel=2.0e-8)
    assert scaled_feed["amounts"] == pytest.approx(
        tuple(10.0 * amount for amount in baseline["amounts"]), rel=2.0e-8
    )
    assert scaled_feed["volume_m3"] == pytest.approx(
        10.0 * baseline["volume_m3"], rel=2.0e-8
    )
    assert scaled_basis["ln_k_provider_basis"] == pytest.approx(
        2.0 * baseline["ln_k_provider_basis"]
    )


def test_source_complete_water_self_ionization_rejects_stale_or_incomplete_inputs() -> None:
    model = _held_water_ionization_model()
    capsule = epcsaft.native_sdk(model)
    options = {"packing_fraction_bounds": (1.0e-6, 0.74), "trace_floor": 1.0e-12}

    stale = _water_self_ionization_spec(model)
    stale["provider_fingerprint"] = "sha256:stale"
    stale["expected_provider_fingerprint"] = "sha256:stale"
    with pytest.raises(ValueError, match="fingerprint"):
        _equilibrium._chemical_solve_provider_water_self_ionization(
            capsule, stale, options
        )

    wrong_basis = _water_self_ionization_spec(model)
    reference = wrong_basis["water_self_ionization_reference"]
    assert isinstance(reference, dict)
    provider_binding = reference["provider_binding"]
    assert isinstance(provider_binding, dict)
    provider_binding["basis_id"] = "stale-basis"
    with pytest.raises(ValueError, match="basis identity"):
        _equilibrium._chemical_solve_provider_water_self_ionization(
            capsule, wrong_basis, options
        )

    missing_source = _water_self_ionization_spec(model)
    reference = missing_source["water_self_ionization_reference"]
    assert isinstance(reference, dict)
    source = reference["source"]
    assert isinstance(source, dict)
    del source["id"]
    with pytest.raises((KeyError, ValueError)):
        _equilibrium._chemical_solve_provider_water_self_ionization(
            capsule, missing_source, options
        )

    nonneutral = _water_self_ionization_spec(model)
    nonneutral["feed_amounts"] = (1.0, 1.0e-8, 0.0)
    with pytest.raises(ValueError, match="pure-water neutral feed"):
        _equilibrium._chemical_solve_provider_water_self_ionization(
            capsule, nonneutral, options
        )

    nonconverged = _water_self_ionization_spec(model)
    reference = nonconverged["water_self_ionization_reference"]
    assert isinstance(reference, dict)
    provider_binding = reference["provider_binding"]
    assert isinstance(provider_binding, dict)
    provider_binding["maximum_reference_convergence_error"] = 1.0e-8
    with pytest.raises(ValueError, match="convergence"):
        _equilibrium._chemical_solve_provider_water_self_ionization(
            capsule, nonconverged, options
        )

    with pytest.raises(ValueError, match="outside packing bounds"):
        _equilibrium._chemical_solve_provider_water_self_ionization(
            capsule,
            _water_self_ionization_spec(model),
            {"packing_fraction_bounds": (0.5, 0.74)},
        )


def test_source_complete_water_self_ionization_rejects_indeterminate_and_false_success() -> None:
    model = _held_water_ionization_model()
    capsule = epcsaft.native_sdk(model)
    spec = _water_self_ionization_spec(model)

    indeterminate = _equilibrium._chemical_solve_provider_water_self_ionization(
        capsule,
        spec,
        {"packing_fraction_bounds": (1.0e-6, 0.74), "test_max_iterations": 0},
    )
    assert indeterminate["accepted"] is False
    assert indeterminate["solver_status"] == "maximum_iterations_exceeded"
    assert indeterminate["final_lambda"] is None

    false_success = _equilibrium._chemical_solve_provider_water_self_ionization(
        capsule,
        spec,
        {"packing_fraction_bounds": (1.0e-6, 0.74), "trace_floor": 1.0e-8},
    )
    assert false_success["solver_status"] == "solve_succeeded"
    assert false_success["accepted"] is False
    assert false_success["physical_status"] == "failed"
    assert false_success["trace_status"] == "at_or_below_floor"


def test_provider_standard_reference_tail_rejects_abi_identity_and_domain_mutations() -> None:
    model = _held_water_ionization_model()

    capsule, table, owners = _copied_sdk_capsule(model)
    assert owners
    table.table_size = _StandardReferenceSdk.source_pressure_max_pa.offset + ctypes.sizeof(
        ctypes.c_double
    )
    with pytest.raises(ValueError, match="standard-reference ABI"):
        _equilibrium._chemical_evaluate_provider_standard_reference(
            capsule, 298.15, 100_000.0, model.parameter_fingerprint
        )

    capsule, table, owners = _copied_sdk_capsule(model)
    assert owners
    table.evaluate_electrolyte_standard_reference = None
    with pytest.raises(ValueError, match="standard-reference ABI"):
        _equilibrium._chemical_evaluate_provider_standard_reference(
            capsule, 298.15, 100_000.0, model.parameter_fingerprint
        )

    capsule, table, owners = _copied_sdk_capsule(model)
    assert owners
    table.electrolyte_standard_reference_result_size -= 1
    with pytest.raises(ValueError, match="standard-reference ABI"):
        _equilibrium._chemical_evaluate_provider_standard_reference(
            capsule, 298.15, 100_000.0, model.parameter_fingerprint
        )

    capsule, table, owners = _copied_sdk_capsule(model)
    wrong_basis = ctypes.create_string_buffer(b"wrong-basis")
    table.helmholtz_basis_id = ctypes.cast(wrong_basis, ctypes.c_char_p)
    with pytest.raises(ValueError, match="standard-reference identity"):
        _equilibrium._chemical_evaluate_provider_standard_reference(
            capsule, 298.15, 100_000.0, model.parameter_fingerprint
        )
    assert owners and wrong_basis

    with pytest.raises(ValueError, match="standard-reference identity"):
        _equilibrium._chemical_evaluate_provider_standard_reference(
            epcsaft.native_sdk(model), 298.15, 100_000.0, "sha256:wrong"
        )

    with pytest.raises(ValueError, match="standard-reference source domain"):
        _equilibrium._chemical_evaluate_provider_standard_reference(
            epcsaft.native_sdk(model), 299.0, 100_000.0, model.parameter_fingerprint
        )


def test_provider_standard_reference_tail_distinguishes_provider_failure() -> None:
    model = _held_water_ionization_model()
    capsule, table, owners = _copied_sdk_capsule(model)

    @_STANDARD_REFERENCE_CALLBACK
    def fail(
        _context: int,
        _temperature: float,
        _pressure: float,
        result: ctypes.POINTER(_StandardReferenceResult),
    ) -> int:
        result.contents.status = 4
        result.contents.error = b"forced reference failure"
        return 4

    table.evaluate_electrolyte_standard_reference = ctypes.cast(
        fail, ctypes.c_void_p
    ).value
    with pytest.raises(ValueError, match="Provider standard-reference evaluation failed"):
        _equilibrium._chemical_evaluate_provider_standard_reference(
            capsule, 298.15, 100_000.0, model.parameter_fingerprint
        )
    assert owners and fail


def _figiel_provider_model() -> epcsaft.EPCSAFT:
    parameters = epcsaft.ParameterBundle.from_catalog(
        "figiel-2025-reference-electrolytes", version=1
    ).select(("water", "sodium-cation", "chloride-anion"))
    return epcsaft.EPCSAFT(parameters)


def test_installed_provider_manufactured_reaction_consumes_exact_phase_and_domain_blocks() -> None:
    model = _figiel_provider_model()
    capsule = epcsaft.native_sdk(model)
    temperature_k = 298.15
    target_amounts = (0.82, 0.08, 0.08)
    target_volume = 1.0e-3
    target = _equilibrium._chemical_evaluate_provider_block(
        capsule,
        temperature_k,
        target_amounts,
        target_volume,
        model.parameter_fingerprint,
    )
    ln_k = -target["gradient"][0] + target["gradient"][1] + target["gradient"][2]
    spec = {
        "species_ids": ("water", "sodium-cation", "chloride-anion"),
        "charges": (0, 1, -1),
        "provider_component_ids": ("water", "sodium-cation", "chloride-anion"),
        "provider_charges": (0, 1, -1),
        "provider_fingerprint": model.parameter_fingerprint,
        "expected_provider_component_ids": ("water", "sodium-cation", "chloride-anion"),
        "expected_provider_charges": (0, 1, -1),
        "expected_provider_fingerprint": model.parameter_fingerprint,
        "balance_matrix": ((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
        "declared_balance_rank": 2,
        "reaction_matrix": ((-1.0, 1.0, 1.0),),
        "feed_amounts": (0.8, 0.1, 0.1),
        "ln_k": (ln_k,),
        "temperature_k": temperature_k,
        "pressure_pa": target["pressure_pa"],
        "complete_closed_system": True,
    }
    _bind_record(spec)

    result = _equilibrium._chemical_solve_provider_manufactured(
        capsule,
        spec,
        {"packing_fraction_bounds": (1.0e-6, 0.74)},
    )

    assert result["profile"] == "installed_provider_manufactured_nonpredictive"
    assert result["accepted"] is True
    final = _equilibrium._chemical_evaluate_provider_block(
        capsule,
        temperature_k,
        result["amounts"],
        result["volume_m3"],
        model.parameter_fingerprint,
    )
    assert 2.0 * result["amounts"][0] + result["amounts"][1] + result["amounts"][2] \
        == pytest.approx(1.8, abs=2.0e-9)
    assert -final["gradient"][0] + final["gradient"][1] + final["gradient"][2] \
        == pytest.approx(ln_k, abs=2.0e-7)
    assert final["pressure_pa"] == pytest.approx(target["pressure_pa"], rel=1.0e-8)
    assert result["parameter_fingerprint"] == model.parameter_fingerprint
    assert result["provider_domain_status"] == "passed"
    assert result["packing_fraction_bounds"] == pytest.approx((1.0e-6, 0.74))
    assert 1.0e-6 < result["packing_fraction"] < 0.74
    assert result["predictive_status"] == "manufactured_nonpredictive"
    assert result["globality_certificate"] == "not_guaranteed"


def test_installed_provider_phase_block_derivatives_are_exact_and_contract_bound() -> None:
    model = _figiel_provider_model()
    capsule = epcsaft.native_sdk(model)
    center = (0.82, 0.08, 0.08, 1.0e-3)
    direction = (0.03, -0.01, -0.01, 2.0e-5)
    step = 2.0e-5

    def evaluate(values: tuple[float, ...]) -> dict[str, object]:
        return _equilibrium._chemical_evaluate_provider_block(
            capsule,
            298.15,
            values[:-1],
            values[-1],
            model.parameter_fingerprint,
        )

    lower = evaluate(
        tuple(value - step * delta for value, delta in zip(center, direction, strict=True))
    )
    result = evaluate(center)
    upper = evaluate(
        tuple(value + step * delta for value, delta in zip(center, direction, strict=True))
    )
    value_directional = (upper["value"] - lower["value"]) / (2.0 * step)
    assert value_directional == pytest.approx(
        sum(value * delta for value, delta in zip(result["gradient"], direction, strict=True)),
        rel=2.0e-8,
        abs=2.0e-9,
    )
    for row in range(len(center)):
        gradient_directional = (upper["gradient"][row] - lower["gradient"][row]) / (
            2.0 * step
        )
        hessian_directional = sum(
            result["hessian"][row * len(center) + column] * direction[column]
            for column in range(len(center))
        )
        assert gradient_directional == pytest.approx(
            hessian_directional, rel=3.0e-8, abs=3.0e-8
        )
    assert result["component_ids"] == ["water", "sodium-cation", "chloride-anion"]
    assert result["charges"] == [0, 1, -1]


def test_provider_manufactured_reaction_rejects_capsule_identity_and_source_domain() -> None:
    model = _figiel_provider_model()
    capsule = epcsaft.native_sdk(model)
    spec = {
        "species_ids": ("water", "sodium-cation", "chloride-anion"),
        "charges": (0, 1, -1),
        "provider_component_ids": ("water", "sodium-cation", "chloride-anion"),
        "provider_charges": (0, 1, -1),
        "provider_fingerprint": model.parameter_fingerprint,
        "expected_provider_component_ids": ("chloride-anion", "sodium-cation", "water"),
        "expected_provider_charges": (-1, 1, 0),
        "expected_provider_fingerprint": model.parameter_fingerprint,
        "balance_matrix": ((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
        "declared_balance_rank": 2,
        "reaction_matrix": ((-1.0, 1.0, 1.0),),
        "feed_amounts": (0.8, 0.1, 0.1),
        "ln_k": (0.0,),
        "temperature_k": 298.15,
        "pressure_pa": 100_000.0,
        "complete_closed_system": True,
    }
    _bind_record(spec)
    with pytest.raises(ValueError, match="component order"):
        _equilibrium._chemical_solve_provider_manufactured(capsule, spec, {})

    fingerprint = copy.deepcopy(spec)
    fingerprint["expected_provider_component_ids"] = fingerprint["provider_component_ids"]
    fingerprint["expected_provider_charges"] = fingerprint["provider_charges"]
    fingerprint["provider_fingerprint"] = "sha256:wrong"
    fingerprint["expected_provider_fingerprint"] = "sha256:wrong"
    with pytest.raises(ValueError, match="fingerprint"):
        _equilibrium._chemical_solve_provider_manufactured(
            capsule, fingerprint, {"packing_fraction_bounds": (1.0e-6, 0.74)}
        )

    outside = copy.deepcopy(spec)
    outside["expected_provider_component_ids"] = outside["provider_component_ids"]
    outside["expected_provider_charges"] = outside["provider_charges"]
    outside["feed_amounts"] = (0.02, 0.49, 0.49)
    with pytest.raises(ValueError, match="source domain"):
        _equilibrium._chemical_solve_provider_manufactured(capsule, outside, {})
