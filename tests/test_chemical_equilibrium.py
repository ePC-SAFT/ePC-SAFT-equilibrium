from __future__ import annotations

import copy
import ctypes
import json
import math
import subprocess
import sys
import textwrap
from pathlib import Path

import epcsaft
import pytest
from chemical_equilibrium_cases import (
    base_system as _base_system,
)
from chemical_equilibrium_cases import (
    bind_records as _bind_record,
)
from chemical_equilibrium_cases import typed_problem as _typed_problem

import epcsaft_equilibrium
from epcsaft_equilibrium import _equilibrium

# gvbelov/Heterogeneous-Equilibrium commit c74b87545a3418262e8e38a7c7a2e31e1b12966a,
# test1.dat blob 5490caab3f7a3f8cdea0d1cb883cce4323902657, SHA-256 below.
# The element potentials and expected amounts are the independent 80-decimal
# calculation recorded in docs/research/2026-07-24-belov-aristova-chemical-equilibrium.md.
_BELOV_ELEMENT_POTENTIALS = (
    -10.0776517753450153959143507476192439107332024990741,
    44.8480568890243540815059574241188924984510845101041,
)
_BELOV_SOURCE_GIBBS = (
    -7.072468897041624,
    -28.765987200068952,
    -26.51300119486301,
    21.44955454258845,
    21.26761422175228,
    -34.35125913287487,
    -55.36542486137457,
    -17.10224518423043,
)
_HELD_WATER_IONIZATION_FINGERPRINT = (
    "sha256:6af6c7aec1106723cf6fa536391b7ba08f2e70ef9fe7064bcb6bb61db18644e8"
)
_FINGERPRINT_CAPACITY = 72
_PROVIDER_ERROR_CAPACITY = 160


class _NeutralReferenceResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("component_count", ctypes.c_size_t),
        ("neutral_basis_row_count", ctypes.c_size_t),
        ("neutral_basis_capacity", ctypes.c_size_t),
        ("contraction_capacity", ctypes.c_size_t),
        ("reference_composition_capacity", ctypes.c_size_t),
        ("neutral_basis", ctypes.POINTER(ctypes.c_double)),
        ("log_fugacity_contractions", ctypes.POINTER(ctypes.c_double)),
        ("reference_composition", ctypes.POINTER(ctypes.c_double)),
        ("derivative_availability", ctypes.c_uint32),
        ("temperature_k", ctypes.c_double),
        ("pressure_pa", ctypes.c_double),
        ("solvent_molar_mass_kg_per_mol", ctypes.c_double),
        ("reference_amount_mol", ctypes.c_double),
        ("reference_number_density_mol_per_m3", ctypes.c_double),
        ("reference_molality_mol_per_kg", ctypes.c_double),
        ("reference_convergence_error", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * _FINGERPRINT_CAPACITY),
        ("helmholtz_basis_id", ctypes.c_char * _FINGERPRINT_CAPACITY),
        ("error", ctypes.c_char * _PROVIDER_ERROR_CAPACITY),
    )


class _NeutralReferenceDerivativeResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("value", _NeutralReferenceResult),
        ("active_parameter_count", ctypes.c_size_t),
        ("pressure_derivative_capacity", ctypes.c_size_t),
        ("parameter_derivative_capacity", ctypes.c_size_t),
        ("pressure_derivatives_per_pa", ctypes.POINTER(ctypes.c_double)),
        ("parameter_derivatives", ctypes.POINTER(ctypes.c_double)),
        ("derivative_availability", ctypes.c_uint32),
        ("source_pressure_min_pa", ctypes.c_double),
        ("source_pressure_max_pa", ctypes.c_double),
        ("maximum_root_residual_pa", ctypes.c_double),
        ("minimum_pressure_density_derivative_pa_m3_per_mol", ctypes.c_double),
        ("maximum_density_condition_number", ctypes.c_double),
        ("reference_derivative_convergence_error", ctypes.c_double),
        ("maximum_relative_root_bracket_width", ctypes.c_double),
        ("maximum_relative_root_density_step", ctypes.c_double),
        ("stable_root_count", ctypes.c_size_t),
        ("selected_stable_root_index", ctypes.c_size_t),
        ("reference_branch", ctypes.c_char * 16),
        ("parameter_fingerprint", ctypes.c_char * _FINGERPRINT_CAPACITY),
        ("topology_fingerprint", ctypes.c_char * _FINGERPRINT_CAPACITY),
        ("helmholtz_basis_id", ctypes.c_char * _FINGERPRINT_CAPACITY),
        ("error", ctypes.c_char * _PROVIDER_ERROR_CAPACITY),
    )


class _ChemicalSdkTable(ctypes.Structure):
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
        ("component_ids", ctypes.c_void_p),
        ("component_charges", ctypes.c_void_p),
        ("evaluate_electrolyte_phase", ctypes.c_void_p),
        ("evaluate_molar_volume_bounds", ctypes.c_void_p),
        ("evaluate_packing_fraction", ctypes.c_void_p),
        ("source_temperature_min_k", ctypes.c_double),
        ("source_temperature_max_k", ctypes.c_double),
        ("total_ion_mole_fraction_max", ctypes.c_double),
        ("ion_solvation_born_result_size", ctypes.c_size_t),
        ("evaluate_ion_solvation_born", ctypes.c_void_p),
        ("aqueous_miac_kij_result_size", ctypes.c_size_t),
        ("evaluate_aqueous_miac_kij", ctypes.c_void_p),
        ("neutral_reference_basis_row_count", ctypes.c_size_t),
        ("neutral_reference_result_size", ctypes.c_size_t),
        ("evaluate_neutral_reference", ctypes.c_void_p),
        ("aqueous_miac_solvation_factor_result_size", ctypes.c_size_t),
        ("evaluate_aqueous_miac_solvation_factor", ctypes.c_void_p),
        ("evaluate_aqueous_miac_kij_batch", ctypes.c_void_p),
        ("evaluate_aqueous_miac_solvation_factor_batch", ctypes.c_void_p),
        ("evaluation_budget_size", ctypes.c_size_t),
        ("evaluate_aqueous_miac_kij_batch_bounded", ctypes.c_void_p),
        ("evaluate_aqueous_miac_solvation_factor_batch_bounded", ctypes.c_void_p),
        ("capability_count", ctypes.c_size_t),
        ("capabilities", ctypes.c_void_p),
        ("evaluate_mixture_phase_lij", ctypes.c_void_p),
        ("evaluate_pure_phase_parameter", ctypes.c_void_p),
        ("ion_fraction_suppression_result_size", ctypes.c_size_t),
        ("evaluate_ion_fraction_suppression", ctypes.c_void_p),
        ("ion_solvation_kij_result_size", ctypes.c_size_t),
        ("evaluate_ion_solvation_kij", ctypes.c_void_p),
        ("ion_solvation_ionic_permittivity_result_size", ctypes.c_size_t),
        ("evaluate_ion_solvation_ionic_permittivity", ctypes.c_void_p),
        ("ion_solvation_solvent_permittivity_result_size", ctypes.c_size_t),
        ("evaluate_ion_solvation_solvent_permittivity", ctypes.c_void_p),
        ("associating_parameterized_result_size", ctypes.c_size_t),
        ("evaluate_diameter_basis_associating_pure_phase_parameters", ctypes.c_void_p),
        ("electrolyte_phase_value_result_size", ctypes.c_size_t),
        ("evaluate_electrolyte_phase_value", ctypes.c_void_p),
        ("neutral_reference_derivative_result_size", ctypes.c_size_t),
        ("evaluate_neutral_reference_derivatives", ctypes.c_void_p),
        ("reacting_phase_parameter_result_size", ctypes.c_size_t),
        ("evaluate_reacting_phase_parameters", ctypes.c_void_p),
    )


_NeutralReferenceCallback = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(_NeutralReferenceResult),
)
_NeutralReferenceDerivativeCallback = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_void_p,
    ctypes.c_size_t,
    ctypes.POINTER(_NeutralReferenceDerivativeResult),
)


def _water_ionization_source() -> dict[str, object]:
    path = Path(__file__).parents[1] / "data/reference/water_self_ionization_iapws_r11_07_2019.json"
    return json.loads(path.read_text())


def _iapws_p_kw(record: dict[str, object], temperature_k: float, density: float) -> float:
    correlation = record["correlation"]
    standard_state = record["standard_state"]
    assert isinstance(correlation, dict)
    assert isinstance(standard_state, dict)
    alpha = correlation["alpha"]
    beta = correlation["beta"]
    ideal = correlation["ideal_gas_p_kw_coefficients"]
    assert isinstance(alpha, list)
    assert isinstance(beta, list)
    assert isinstance(ideal, list)
    q_value = density * math.exp(
        alpha[0] + alpha[1] / temperature_k + alpha[2] * density ** (2.0 / 3.0) / temperature_k**2
    )
    p_kw_ideal = (
        ideal[0]
        + ideal[1] / temperature_k
        + ideal[2] / temperature_k**2
        + ideal[3] / temperature_k**3
    )
    return (
        -2.0
        * correlation["ion_coordination_number"]
        * (
            math.log10(1.0 + q_value)
            - q_value
            / (q_value + 1.0)
            * density
            * (beta[0] + beta[1] / temperature_k + beta[2] * density)
        )
        + p_kw_ideal
        + 2.0
        * math.log10(
            standard_state["standard_molality_mol_per_kg"]
            * standard_state["solvent_molar_mass_kg_per_mol"]
        )
    )


def _belov_aristova_gas_system() -> tuple[dict[str, object], tuple[float, ...]]:
    temperature_k = 2000.0
    pressure_pa = 100_000.0
    standard_pressure_pa = 101_325.0
    reactions = (
        (-2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (-3.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, -2.0, 1.0, 0.0, 0.0, 0.0),
        (-1.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0),
        (-2.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0),
        (-1.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 1.0),
    )
    concentration_shift = math.log(8.31446261815324 * temperature_k / standard_pressure_pa)
    manufactured_gibbs = tuple(value + concentration_shift for value in _BELOV_SOURCE_GIBBS)
    ln_k = tuple(
        -sum(row[index] * manufactured_gibbs[index] for index in range(8)) for row in reactions
    )
    spec: dict[str, object] = {
        "species_ids": ("O", "O2", "O3", "C", "C2", "CO", "CO2", "C2O"),
        "charges": (0,) * 8,
        "provider_fingerprint": (
            "sha256:2f24904e5db64e96aa2ee96c53bd80863c9518b6c2ddd914e266181fe7459f84"
        ),
        "molar_masses_kg_per_mol": (
            0.016,
            0.032,
            0.048,
            0.012,
            0.024,
            0.028,
            0.044,
            0.040,
        ),
        "balance_matrix": (
            (0.0, 0.0, 0.0, 1.0, 2.0, 1.0, 1.0, 2.0),
            (1.0, 2.0, 3.0, 0.0, 0.0, 1.0, 2.0, 1.0),
        ),
        "reaction_matrix": reactions,
        "feed_amounts": (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0),
        "ln_k": ln_k,
        "temperature_k": temperature_k,
        "pressure_pa": pressure_pa,
    }
    _bind_record(spec)
    expected = (
        5.969938280702114e-17,
        5.249059580600665e-27,
        1.838054910648896e-47,
        1.7496236329525488e-5,
        4.996074090543112e-1,
        9.992323142113352e-1,
        4.453885591412746e-11,
        7.676856995870245e-4,
    )
    return spec, expected


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
        first_directional = (upper["amounts"][species] - lower["amounts"][species]) / (2.0 * step)
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
                hessians[species * dimension * dimension + row * dimension + column]
                * direction[column]
                for column in range(dimension)
            )
            assert jacobian_directional == pytest.approx(exact_second, rel=2.0e-9, abs=5.0e-11)


def test_amount_chart_classifies_trace_floor_without_zeroing_species() -> None:
    trace_floor = 1.0e-10
    result = _amount_chart((0, 0), (math.log(0.5 * trace_floor), 0.0), trace_floor)

    assert result["trace_status"] == "at_or_below_floor"
    assert result["minimum_amount"] == pytest.approx(0.5 * trace_floor, rel=2.0e-15)
    assert all(value > 0.0 for value in result["amounts"])


def test_amount_chart_rejects_one_sided_ionic_topology() -> None:
    with pytest.raises(ValueError, match="both cations and anions"):
        _amount_chart((1, 0), (0.0,))


def _manufactured_solve(
    spec: dict[str, object], options: dict[str, object] | None = None
) -> epcsaft_equilibrium.ChemicalEquilibriumResult:
    options = options or {}
    return epcsaft_equilibrium.chemical_equilibrium(
        epcsaft_equilibrium.IdealGasPhase(
            model_fingerprint=spec["provider_fingerprint"],
            reference_id="provider-helmholtz-coordinate-basis",
        ),
        spec["temperature_k"] * epcsaft.unit_registry.kelvin,
        spec["pressure_pa"] * epcsaft.unit_registry.pascal,
        _typed_problem(
            spec,
            minimum_amount_mol=options.get("trace_floor", 1.0e-12),
        ),
    )


@pytest.mark.parametrize(
    ("spec", "expected_amounts", "expected_volume"),
    (
        (_base_system(), (0.2, 0.8), None),
        (
            {
                **_base_system(),
                "balance_matrix": ((1.0, 2.0),),
                "reaction_matrix": ((-2.0, 1.0),),
                "molar_masses_kg_per_mol": (1.0, 2.0),
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
    diagnostics = result.diagnostics
    assert result.thermodynamic_model == "ideal_gas"
    assert result.amounts_mol == pytest.approx(expected_amounts, rel=2.0e-8, abs=2.0e-10)
    assert result.volume_m3 == pytest.approx(volume, rel=2.0e-8)
    assert diagnostics.solver_status == "solve_succeeded"
    assert diagnostics.numerical_status == "passed"
    assert diagnostics.physical_status == "passed"
    assert diagnostics.local_minimum_status == "passed"
    assert diagnostics.trace_status == "interior"
    assert diagnostics.chemical_certification_level == "LOCAL_EQUILIBRIUM"
    assert diagnostics.boundary_status == "strict_interior"
    assert diagnostics.globality_status == "not_guaranteed"
    assert diagnostics.balance_inf_norm <= 1.0e-9
    assert diagnostics.charge_inf_norm <= 1.0e-12
    assert diagnostics.pressure_relative_residual <= 1.0e-8
    assert diagnostics.reaction_affinity_inf_norm <= 1.0e-7
    assert diagnostics.kkt_stationarity_inf_norm <= 1.0e-7
    assert diagnostics.complementarity_inf_norm <= 1.0e-7


def test_manufactured_reaction_reports_exact_conditioned_implicit_sensitivities() -> None:
    spec = _base_system()
    _bind_record(spec)
    native = _equilibrium._chemical_equilibrium(None, spec, None, None, 1.0e-12)

    sensitivities = native["sensitivities"]
    assert sensitivities["status"] == "available"
    assert sensitivities["parameter_order"] == (
        "balance_total[0]",
        "ln_k_provider_basis[0]",
        "pressure_pa",
    )
    assert sensitivities["parameter_fingerprint"] == "sha256:manufactured"
    assert sensitivities["chart_topology"] == "neutral_log_amounts[2]+log_volume"
    assert sensitivities["kkt_rank"] == sensitivities["kkt_dimension"] == 3
    assert sensitivities["condition_number_inf"] < 1.0e6
    assert sensitivities["active_lower_bounds"] == ()
    assert sensitivities["active_upper_bounds"] == ()
    assert sensitivities["active_constraint_bounds"] == ()

    amount_derivatives = sensitivities["amount_derivatives"]
    assert amount_derivatives[0] == pytest.approx((0.2, 0.8), abs=2.0e-9)
    assert amount_derivatives[1] == pytest.approx((-0.16, 0.16), abs=2.0e-9)
    assert amount_derivatives[2] == pytest.approx((0.0, 0.0), abs=2.0e-12)

    expected_volume = 8.31446261815324 * 350.0 / 200_000.0
    assert sensitivities["volume_derivatives"] == pytest.approx(
        (
            expected_volume,
            0.0,
            -expected_volume / 200_000.0,
        ),
        rel=2.0e-8,
        abs=2.0e-13,
    )

    steps = (1.0e-5, 1.0e-5, 5.0)
    for parameter, step in enumerate(steps):
        perturbed: list[dict[str, object]] = []
        for direction in (-1.0, 1.0):
            trial = copy.deepcopy(spec)
            if parameter == 0:
                trial["feed_amounts"] = (
                    1.0 + direction * step,
                    0.0,
                )
            elif parameter == 1:
                trial["ln_k"] = (
                    spec["ln_k"][0] + direction * step,  # type: ignore[index]
                )
            else:
                trial["pressure_pa"] = spec["pressure_pa"] + direction * step
            _bind_record(trial)
            perturbed.append(
                _equilibrium._chemical_equilibrium(
                    None, trial, None, None, 1.0e-12
                )
            )
        finite_difference_amounts = tuple(
            (
                perturbed[1]["amounts"][species]
                - perturbed[0]["amounts"][species]
            )
            / (2.0 * step)
            for species in range(2)
        )
        finite_difference_volume = (
            perturbed[1]["volume_m3"] - perturbed[0]["volume_m3"]
        ) / (2.0 * step)
        assert amount_derivatives[parameter] == pytest.approx(
            finite_difference_amounts,
            rel=2.0e-6,
            abs=2.0e-10,
        )
        assert sensitivities["volume_derivatives"][parameter] == pytest.approx(
            finite_difference_volume,
            rel=2.0e-6,
            abs=2.0e-12,
        )


def test_manufactured_sensitivity_is_invariant_to_reaction_scaling() -> None:
    reaction_scale = 1.0e-11
    base = _base_system()
    _bind_record(base)
    scaled = copy.deepcopy(base)
    scaled["reaction_matrix"] = ((-reaction_scale, reaction_scale),)
    scaled["ln_k"] = (reaction_scale * math.log(4.0),)
    _bind_record(scaled)

    base_result = _equilibrium._chemical_equilibrium(
        None, base, None, None, 1.0e-12
    )["sensitivities"]
    scaled_result = _equilibrium._chemical_equilibrium(
        None, scaled, None, None, 1.0e-12
    )["sensitivities"]

    assert base_result["status"] == scaled_result["status"] == "available"
    assert scaled_result["condition_number_inf"] == pytest.approx(
        base_result["condition_number_inf"], rel=2.0e-10
    )
    assert scaled_result["amount_derivatives"][0] == pytest.approx(
        base_result["amount_derivatives"][0], abs=2.0e-9
    )
    assert tuple(
        reaction_scale * value
        for value in scaled_result["amount_derivatives"][1]
    ) == pytest.approx(base_result["amount_derivatives"][1], abs=2.0e-9)
    assert scaled_result["amount_derivatives"][2] == pytest.approx(
        base_result["amount_derivatives"][2], abs=2.0e-12
    )


def test_manufactured_sensitivity_tracks_exact_species_permutation() -> None:
    base = _base_system()
    _bind_record(base)
    permutation = (1, 0)
    permuted = copy.deepcopy(base)
    for field in (
        "species_ids",
        "charges",
        "molar_masses_kg_per_mol",
        "feed_amounts",
    ):
        values = base[field]
        permuted[field] = tuple(values[index] for index in permutation)
    permuted["balance_matrix"] = tuple(
        tuple(row[index] for index in permutation)
        for row in base["balance_matrix"]  # type: ignore[union-attr]
    )
    permuted["reaction_matrix"] = tuple(
        tuple(row[index] for index in permutation)
        for row in base["reaction_matrix"]  # type: ignore[union-attr]
    )
    _bind_record(permuted)

    base_result = _equilibrium._chemical_equilibrium(
        None, base, None, None, 1.0e-12
    )["sensitivities"]
    permuted_result = _equilibrium._chemical_equilibrium(
        None, permuted, None, None, 1.0e-12
    )["sensitivities"]

    assert permuted_result["parameter_order"] == base_result["parameter_order"]
    assert permuted_result["condition_number_inf"] == pytest.approx(
        base_result["condition_number_inf"], rel=2.0e-10
    )
    for parameter in range(3):
        assert permuted_result["amount_derivatives"][parameter] == pytest.approx(
            tuple(
                base_result["amount_derivatives"][parameter][index]
                for index in permutation
            ),
            abs=2.0e-9,
        )


def test_manufactured_reaction_with_ill_conditioned_kkt_has_no_sensitivity() -> None:
    balance_separation = 8.0e-6
    spec = _base_system()
    spec["species_ids"] = ("A", "B", "C")
    spec["charges"] = (0, 0, 0)
    spec["molar_masses_kg_per_mol"] = (1.0, 1.0, 1.0)
    spec["balance_matrix"] = (
        (1.0, 1.0, 1.0 + balance_separation),
    )
    spec["reaction_matrix"] = ((-1.0, 1.0, 0.0),)
    spec["feed_amounts"] = (0.5, 0.5, 1.0)
    spec["ln_k"] = (0.0,)
    spec["temperature_k"] = 300.0
    spec["pressure_pa"] = 8.31446261815324 * 300.0
    _bind_record(spec)

    native = _equilibrium._chemical_equilibrium(None, spec, None, None, 1.0e-12)

    sensitivities = native["sensitivities"]
    assert sensitivities["status"] == "unavailable"
    assert sensitivities["failure_reason"] == "ill_conditioned_kkt_jacobian"
    assert sensitivities["kkt_rank"] == sensitivities["kkt_dimension"] == 4
    assert sensitivities["condition_number_inf"] > 1.0e6
    assert sensitivities["parameter_order"] == ()
    assert sensitivities["amount_derivatives"] == ()
    assert sensitivities["volume_derivatives"] == ()


@pytest.mark.parametrize("trace_floor", (1.0e-50, 1.0e-55))
def test_belov_aristova_gas_restriction_resolves_extreme_positive_traces(
    trace_floor: float,
) -> None:
    spec, expected = _belov_aristova_gas_system()

    result = _manufactured_solve(spec, {"trace_floor": trace_floor})

    assert result.diagnostics.trace_status == "interior"
    amounts = result.amounts_mol
    assert tuple(map(math.log, amounts)) == pytest.approx(
        tuple(map(math.log, expected)), abs=2.0e-6
    )
    balances = spec["balance_matrix"]
    assert tuple(
        sum(row[index] * amounts[index] for index in range(8))
        for row in balances  # type: ignore[union-attr]
    ) == pytest.approx((2.0, 1.0), abs=1.0e-9)
    assert (8.31446261815324 * 2000.0 * sum(amounts) / result.volume_m3) == pytest.approx(
        100_000.0, rel=1.0e-8
    )
    concentration_shift = math.log(8.31446261815324 * 2000.0 / 101_325.0)
    potentials = tuple(
        _BELOV_SOURCE_GIBBS[index]
        + concentration_shift
        + math.log(amounts[index] / result.volume_m3)
        for index in range(8)
    )
    reactions = spec["reaction_matrix"]
    assert tuple(
        sum(row[index] * potentials[index] for index in range(8))
        for row in reactions  # type: ignore[union-attr]
    ) == pytest.approx((0.0,) * 6, abs=1.0e-7)
    assert tuple(
        potentials[index]
        + balances[0][index] * _BELOV_ELEMENT_POTENTIALS[0]  # type: ignore[index]
        + balances[1][index] * _BELOV_ELEMENT_POTENTIALS[1]  # type: ignore[index]
        for index in range(8)
    ) == pytest.approx((0.0,) * 8, abs=1.0e-7)
    diagnostics = result.diagnostics
    assert diagnostics.balance_inf_norm <= 1.0e-9
    assert diagnostics.pressure_relative_residual <= 1.0e-8
    assert diagnostics.reaction_affinity_inf_norm <= 1.0e-7
    assert diagnostics.kkt_stationarity_inf_norm <= 1.0e-7
    assert diagnostics.local_minimum_status == "passed"
    assert diagnostics.chemical_certification_level == "LOCAL_EQUILIBRIUM"
    assert diagnostics.boundary_status == "strict_interior"


@pytest.mark.parametrize(
    "variant", ("reaction_basis", "species_order", "conservation_gauge", "feed_scale")
)
def test_belov_aristova_trace_solution_is_coordinate_invariant(variant: str) -> None:
    spec, expected = _belov_aristova_gas_system()
    options: dict[str, object] = {"trace_floor": 1.0e-50}

    if variant == "reaction_basis":
        reactions = spec["reaction_matrix"]
        ln_k = spec["ln_k"]
        spec["reaction_matrix"] = (
            tuple(
                reactions[0][index] + reactions[1][index]  # type: ignore[index]
                for index in range(8)
            ),
            *reactions[1:],  # type: ignore[index]
        )
        spec["ln_k"] = (ln_k[0] + ln_k[1], *ln_k[1:])  # type: ignore[index]
        _bind_record(spec)
    elif variant == "species_order":
        permutation = (7, 5, 3, 1, 6, 4, 2, 0)
        for field in (
            "species_ids",
            "charges",
            "molar_masses_kg_per_mol",
            "feed_amounts",
        ):
            values = spec[field]
            spec[field] = tuple(values[index] for index in permutation)  # type: ignore[index]
        spec["balance_matrix"] = tuple(
            tuple(row[index] for index in permutation)
            for row in spec["balance_matrix"]  # type: ignore[union-attr]
        )
        spec["reaction_matrix"] = tuple(
            tuple(row[index] for index in permutation)
            for row in spec["reaction_matrix"]  # type: ignore[union-attr]
        )
        expected = tuple(expected[index] for index in permutation)
    elif variant == "conservation_gauge":
        balances = spec["balance_matrix"]
        spec["balance_matrix"] = (
            tuple(
                balances[0][index] + balances[1][index]  # type: ignore[index]
                for index in range(8)
            ),
            tuple(
                2.0 * balances[0][index] - balances[1][index]  # type: ignore[index]
                for index in range(8)
            ),
        )
        _bind_record(spec)
    else:
        scale = 3.0
        spec["feed_amounts"] = tuple(
            scale * value
            for value in spec["feed_amounts"]  # type: ignore[union-attr]
        )
        expected = tuple(scale * value for value in expected)
        _bind_record(spec)

    result = _manufactured_solve(spec, options)

    assert tuple(map(math.log, result.amounts_mol)) == pytest.approx(
        tuple(map(math.log, expected)), abs=3.0e-6
    )
    assert result.diagnostics.reaction_affinity_inf_norm <= 1.0e-7


def test_manufactured_charged_reaction_uses_exact_electroneutral_chart() -> None:
    temperature_k = 300.0
    spec = {
        **_base_system(),
        "species_ids": ("A", "C+", "D-"),
        "charges": (0, 1, -1),
        "molar_masses_kg_per_mol": (2.0, 1.0, 1.0),
        "balance_matrix": ((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
        "reaction_matrix": ((-1.0, 1.0, 1.0),),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (math.log(1.0 / 3.0),),
        "temperature_k": temperature_k,
        "pressure_pa": 8.31446261815324 * temperature_k,
    }
    _bind_record(spec)

    result = _manufactured_solve(spec)

    assert result.amounts_mol == pytest.approx((0.5, 0.5, 0.5), rel=3.0e-8)
    assert result.volume_m3 == pytest.approx(1.5, rel=3.0e-8)
    assert result.diagnostics.charge_inf_norm <= 2.0e-15


def test_manufactured_ultra_trace_charged_share_stays_differentiable() -> None:
    """Both charged softmax references and explicit shares stay interior."""

    temperature_k = 300.0
    trace_floor = 1.0e-18
    records = tuple(
        {
            "source_id": f"manufactured:reaction-{index}",
            "reference_id": "provider-helmholtz-coordinate-basis",
            "reaction_orientation": "products_positive",
            "conversion_id": "already-provider-basis",
            "dimensionless": True,
            "temperature_k": temperature_k,
            "pressure_pa": 8.31446261815324 * temperature_k,
        }
        for index in range(2)
    )

    def solve(log_ratio_constants: tuple[float, float], floor: float) -> dict[str, object]:
        spec = {
            **_base_system(),
            "species_ids": ("C1+", "C2+", "D1-", "D2-"),
            "charges": (1, 1, -1, -1),
            "molar_masses_kg_per_mol": (1.0, 1.0, 1.0, 1.0),
            "balance_matrix": ((1.0, 1.0, 1.0, 1.0), (1.0, 1.0, -1.0, -1.0)),
            "reaction_matrix": ((1.0, -1.0, 0.0, 0.0), (0.0, 0.0, 1.0, -1.0)),
            "feed_amounts": (1.0, 0.0, 1.0, 0.0),
            "ln_k": log_ratio_constants,
            "temperature_k": temperature_k,
            "pressure_pa": 8.31446261815324 * temperature_k,
            "equilibrium_constant_records": records,
        }
        spec["conserved_totals"] = (2.0, 0.0)
        return _equilibrium._chemical_equilibrium(None, spec, None, None, floor)

    explicit_trace = solve((-40.0, 0.0), trace_floor)
    assert explicit_trace["amounts"][0] > trace_floor
    assert explicit_trace["sensitivities"]["status"] == "available"
    assert explicit_trace["sensitivities"]["active_lower_bounds"] == ()
    assert explicit_trace["sensitivities"]["active_upper_bounds"] == ()

    reference_trace = solve((40.0, 40.0), trace_floor)
    assert reference_trace["amounts"][1] > trace_floor
    assert reference_trace["amounts"][3] > trace_floor
    assert reference_trace["sensitivities"]["status"] == "available"
    assert reference_trace["sensitivities"]["active_lower_bounds"] == ()
    assert reference_trace["sensitivities"]["active_upper_bounds"] == ()

    rejected_trace = solve((40.0, 40.0), 1.0e-12)
    assert rejected_trace["sensitivities"]["status"] == "unavailable"
    assert rejected_trace["sensitivities"]["failure_reason"] == (
        "active_set_change_not_differentiable"
    )
    assert rejected_trace["sensitivities"]["active_trace_species"] == (1, 3)


def test_manufactured_equilibrium_is_gauge_scale_and_reaction_basis_invariant() -> None:
    base = _base_system()
    base["feed_amounts"] = (1.0, 0.0)
    _bind_record(base)
    plain = _manufactured_solve(base)
    gauge = copy.deepcopy(base)
    gauge["balance_matrix"] = ((3.25, 3.25),)
    _bind_record(gauge)
    gauged = _manufactured_solve(gauge)
    assert gauged.amounts_mol == pytest.approx(plain.amounts_mol, rel=2.0e-10)

    scaled = copy.deepcopy(base)
    scaled["feed_amounts"] = (7.0, 0.0)
    _bind_record(scaled)
    scaled_result = _manufactured_solve(scaled)
    assert scaled_result.amounts_mol == pytest.approx(
        tuple(7.0 * value for value in plain.amounts_mol), rel=3.0e-8
    )
    assert scaled_result.volume_m3 == pytest.approx(7.0 * plain.volume_m3, rel=3.0e-8)

    three = {
        **_base_system(),
        "species_ids": ("A", "B", "C"),
        "charges": (0, 0, 0),
        "molar_masses_kg_per_mol": (1.0, 1.0, 1.0),
        "balance_matrix": ((1.0, 1.0, 1.0),),
        "reaction_matrix": ((-1.0, 1.0, 0.0), (0.0, -1.0, 1.0)),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (math.log(2.0), math.log(3.0)),
    }
    _bind_record(three)
    first_basis = _manufactured_solve(three)
    changed_basis = copy.deepcopy(three)
    changed_basis["reaction_matrix"] = ((-1.0, 0.0, 1.0), (0.0, -1.0, 1.0))
    changed_basis["ln_k"] = (math.log(6.0), math.log(3.0))
    _bind_record(changed_basis)
    second_basis = _manufactured_solve(changed_basis)
    assert first_basis.amounts_mol == pytest.approx(second_basis.amounts_mol, rel=3.0e-9)


def test_manufactured_reaction_solve_has_no_feed_scaled_amount_cap() -> None:
    target = (0.5, 500.0, 0.5)
    volume = sum(target)
    spec = {
        **_base_system(),
        "species_ids": ("A", "B", "C"),
        "charges": (0, 0, 0),
        "molar_masses_kg_per_mol": (1001.0, 1.0, 1.0),
        "balance_matrix": ((1000.0, 1.0, 0.0), (1.0, 0.0, 1.0)),
        "reaction_matrix": ((-1.0, 1000.0, 1.0),),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (1000.0 * math.log(target[1] / volume),),
        "temperature_k": 300.0,
        "pressure_pa": 8.31446261815324 * 300.0,
    }
    _bind_record(spec)

    result = _manufactured_solve(spec, {"trace_floor": 0.1})

    assert result.amounts_mol == pytest.approx(target, rel=2.0e-7, abs=2.0e-8)


def test_manufactured_nlp_has_exact_directional_gradient_jacobian_and_hessian() -> None:
    spec = _base_system()
    _bind_record(spec)
    pressure_over_rt = spec["pressure_pa"] / (8.31446261815324 * spec["temperature_k"])
    center = (math.log(0.3), math.log(0.7), math.log(1.0 / pressure_over_rt))
    direction = (0.2, -0.4, 0.3)
    multipliers = (0.37,)
    step = 2.0e-5

    def evaluate(variables: tuple[float, ...]) -> dict[str, object]:
        return _equilibrium._chemical_evaluate_manufactured_nlp(spec, variables, multipliers)

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
    constraint_directional = (upper["constraints"][0] - lower["constraints"][0]) / (2.0 * step)
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
        assert gradient_directional == pytest.approx(hessian_directional, rel=3.0e-9, abs=3.0e-10)


def test_manufactured_solver_rejects_trace_false_success() -> None:
    spec = _base_system()
    trace = copy.deepcopy(spec)
    trace["ln_k"] = (math.log(1.0e-10),)
    _bind_record(trace)
    with pytest.raises(epcsaft_equilibrium.ChemicalEquilibriumError) as failed:
        _manufactured_solve(trace, {"trace_floor": 1.0e-8})
    diagnostics = failed.value.diagnostics
    assert diagnostics.physical_status == "failed"
    assert diagnostics.trace_status == "at_or_below_floor"
    assert diagnostics.chemical_certification_level == "FEASIBLE_ONLY"
    assert diagnostics.globality_status == "not_guaranteed"
    native = _equilibrium._chemical_equilibrium(
        None, trace, None, None, 1.0e-8
    )
    sensitivities = native["sensitivities"]
    assert sensitivities["status"] == "unavailable"
    assert sensitivities["failure_reason"] == "active_set_change_not_differentiable"
    assert sensitivities["active_trace_species"] == (1,)
    assert sensitivities["parameter_order"] == ()


def test_manufactured_charged_solution_is_species_order_invariant() -> None:
    temperature_k = 300.0
    original = {
        **_base_system(),
        "species_ids": ("A", "C+", "D-"),
        "charges": (0, 1, -1),
        "molar_masses_kg_per_mol": (2.0, 1.0, 1.0),
        "balance_matrix": ((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
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
        "molar_masses_kg_per_mol",
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
    assert tuple(reordered.amounts_mol[index] for index in inverse) == pytest.approx(
        baseline.amounts_mol, rel=4.0e-8
    )


def _figiel_provider_model(
    components: tuple[str, ...] = ("water", "sodium-cation", "chloride-anion"),
) -> epcsaft.Mixture:
    parameters = epcsaft.Parameters.from_catalog(
        "figiel-2025-reference-electrolytes",
        components=components,
        version=1,
    )
    return epcsaft.Mixture(parameters)


@pytest.mark.parametrize("with_sensitivity", (False, True))
def test_unreachable_provider_state_fails_closed_in_subprocess(
    with_sensitivity: bool,
) -> None:
    script = textwrap.dedent(
        f"""
        import json

        import epcsaft
        import epcsaft_equilibrium

        component_ids = ("water", "sodium-cation", "chloride-anion")
        parameters = epcsaft.Parameters.from_catalog(
            "figiel-2025-reference-electrolytes",
            components=component_ids,
            version=1,
        )
        model = epcsaft.Mixture(parameters)
        phase = epcsaft_equilibrium.ProviderPhase(
            model=model,
            expected_parameter_fingerprint=model.parameter_fingerprint,
            admissible_packing_fraction_interval=(1.0e-6, 0.74),
        )
        problem = epcsaft_equilibrium.ChemicalEquilibriumProblem(
            species_ids=component_ids,
            charges=(0, 1, -1),
            molar_masses_kg_per_mol=(2.0, 1.0, 1.0),
            balance_matrix=((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
            conserved_totals=(1.8, 0.0),
            reaction_matrix=((-1.0, 1.0, 1.0),),
            feed_amounts_mol=(0.8, 0.1, 0.1),
            equilibrium_constants=(
                epcsaft_equilibrium.ChemicalEquilibriumConstant(
                    ln_value=0.0,
                    source_id="issue-81-minimal-reproducer",
                    reference_id="provider-helmholtz-coordinate-basis",
                    reaction_orientation="products_positive",
                    conversion_id="already-provider-basis",
                    dimensionless=True,
                ),
            ),
            strict_interior_amount_floor_mol=1.0e-12,
        )
        sensitivity = (
            epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest()
            if {with_sensitivity!r}
            else None
        )
        try:
            epcsaft_equilibrium.chemical_equilibrium(
                phase,
                298.15 * epcsaft.unit_registry.kelvin,
                100000.0 * epcsaft.unit_registry.pascal,
                problem,
                sensitivity_request=sensitivity,
            )
        except epcsaft_equilibrium.ChemicalEquilibriumError as error:
            print(json.dumps({{
                "kind": type(error).__name__,
                "solver_status": error.diagnostics.solver_status,
                "failure_reason": error.diagnostics.failure_reason,
                "callback_error": error.diagnostics.callback_error,
            }}, sort_keys=True))
            raise SystemExit(0)
        raise SystemExit(3)
        """
    )
    completed = subprocess.run(
        (sys.executable, "-X", "faulthandler", "-c", script),
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    payload = json.loads(completed.stdout)
    assert payload["kind"] == "ChemicalEquilibriumError"
    assert payload["solver_status"]
    assert "ion mole fraction exceeds the parameter source domain" in payload["callback_error"]


def _held_water_ionization_problem() -> tuple[
    epcsaft.Mixture, dict[str, object], dict[str, object]
]:
    source = _water_ionization_source()
    state = source["state"]
    standard_state = source["standard_state"]
    values = source["values"]
    assert isinstance(state, dict)
    assert isinstance(standard_state, dict)
    assert isinstance(values, dict)
    components = ("water", "hydronium-cation", "hydroxide-anion")
    parameters = epcsaft.Parameters.from_catalog(
        "held-2008-water-self-ionization",
        components=components,
        version=1,
    )
    model = epcsaft.Mixture(parameters)
    assert model.parameter_fingerprint == _HELD_WATER_IONIZATION_FINGERPRINT
    conversion = math.log(
        standard_state["standard_molality_mol_per_kg"]
        * standard_state["solvent_molar_mass_kg_per_mol"]
    )
    spec: dict[str, object] = {
        "species_ids": components,
        "charges": (0, 1, -1),
        "provider_fingerprint": model.parameter_fingerprint,
        "molar_masses_kg_per_mol": (0.01801528, 0.01902322, 0.01700734),
        "balance_matrix": ((2.0, 3.0, 1.0), (1.0, 1.0, 1.0)),
        "reaction_matrix": ((-2.0, 1.0, 1.0),),
        "feed_amounts": (1.0, 0.0, 0.0),
        "ln_k": (values["ln_kw"],),
        "equilibrium_constant_records": (
            {
                "source_id": source["source"]["id"],
                "reference_id": standard_state["id"],
                "reaction_orientation": "products_positive",
                "conversion_id": "source-standard-state-to-provider-neutral-reference",
                "dimensionless": True,
                "temperature_k": state["temperature_k"],
                "pressure_pa": state["pressure_pa"],
            },
        ),
        "temperature_k": state["temperature_k"],
        "pressure_pa": state["pressure_pa"],
    }
    return (
        model,
        spec,
        {
            "id": standard_state["id"],
            "activity_scale_id": standard_state["activity_scale_id"],
            "log_activity_scale_factors": (0.0, conversion, conversion),
            "reference_pressure_pa": state["pressure_pa"],
        },
    )


def _provider_solve(
    model: epcsaft.Mixture,
    spec: dict[str, object],
    source_standard_state: dict[str, object] | None = None,
    sensitivity_request: epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest
    | None = None,
) -> epcsaft_equilibrium.ChemicalEquilibriumResult:
    standard_state = (
        None
        if source_standard_state is None
        else epcsaft_equilibrium.ChemicalStandardState(
            id=source_standard_state["id"],
            activity_scale_id=source_standard_state["activity_scale_id"],
            log_activity_scale_factors=tuple(source_standard_state["log_activity_scale_factors"]),
            reference_pressure_pa=source_standard_state["reference_pressure_pa"],
        )
    )
    return epcsaft_equilibrium.chemical_equilibrium(
        epcsaft_equilibrium.ProviderPhase(
            model=model,
            expected_parameter_fingerprint=spec["provider_fingerprint"],
            admissible_packing_fraction_interval=(1.0e-6, 0.74),
        ),
        spec["temperature_k"] * epcsaft.unit_registry.kelvin,
        spec["pressure_pa"] * epcsaft.unit_registry.pascal,
        _typed_problem(spec, source_standard_state=standard_state),
        sensitivity_request=sensitivity_request,
    )


def _linear_reference_capsule(
    model: epcsaft.Mixture,
    pressure_derivatives_per_pa: tuple[float, ...],
) -> tuple[object, tuple[object, ...]]:
    """Test-only broad linear reference layered over one installed Provider point."""
    source_pressure_pa = 100_000.0
    provider_capsule = epcsaft.native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    provider_pointer = get_pointer(provider_capsule, b"epcsaft.native_sdk.v1")
    provider_table = ctypes.cast(
        provider_pointer, ctypes.POINTER(_ChemicalSdkTable)
    ).contents
    assert provider_table.table_size == ctypes.sizeof(_ChemicalSdkTable)
    assert len(pressure_derivatives_per_pa) == provider_table.neutral_reference_basis_row_count

    table = _ChemicalSdkTable.from_buffer_copy(
        ctypes.string_at(provider_pointer, provider_table.table_size)
    )
    provider_value = _NeutralReferenceCallback(table.evaluate_neutral_reference)
    provider_derivatives = _NeutralReferenceDerivativeCallback(
        table.evaluate_neutral_reference_derivatives
    )

    def apply_linear_reference(
        pressure_pa: float,
        value: _NeutralReferenceResult,
    ) -> None:
        for row, derivative in enumerate(pressure_derivatives_per_pa):
            value.log_fugacity_contractions[row] += derivative * (
                pressure_pa - source_pressure_pa
            )
        value.pressure_pa = pressure_pa

    def evaluate_value(
        context: int,
        fingerprint: bytes,
        temperature_k: float,
        pressure_pa: float,
        result: ctypes.POINTER(_NeutralReferenceResult),
    ) -> int:
        status = provider_value(
            context,
            fingerprint,
            temperature_k,
            source_pressure_pa,
            result,
        )
        if status == 0:
            apply_linear_reference(pressure_pa, result.contents)
        return status

    def evaluate_derivatives(
        context: int,
        fingerprint: bytes,
        topology_fingerprint: bytes,
        temperature_k: float,
        pressure_pa: float,
        active_parameters: int,
        active_parameter_count: int,
        result: ctypes.POINTER(_NeutralReferenceDerivativeResult),
    ) -> int:
        status = provider_derivatives(
            context,
            fingerprint,
            topology_fingerprint,
            temperature_k,
            source_pressure_pa,
            active_parameters,
            active_parameter_count,
            result,
        )
        if status == 0:
            native = result.contents
            apply_linear_reference(pressure_pa, native.value)
            for row, derivative in enumerate(pressure_derivatives_per_pa):
                native.pressure_derivatives_per_pa[row] = derivative
            native.derivative_availability = 2
            native.source_pressure_min_pa = 50_000.0
            native.source_pressure_max_pa = 150_000.0
        return status

    value_callback = _NeutralReferenceCallback(evaluate_value)
    derivative_callback = _NeutralReferenceDerivativeCallback(evaluate_derivatives)
    table.evaluate_neutral_reference = ctypes.cast(
        value_callback, ctypes.c_void_p
    ).value
    table.evaluate_neutral_reference_derivatives = ctypes.cast(
        derivative_callback, ctypes.c_void_p
    ).value
    name = ctypes.create_string_buffer(b"epcsaft.native_sdk.v1")
    new_capsule = ctypes.pythonapi.PyCapsule_New
    new_capsule.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p)
    new_capsule.restype = ctypes.py_object
    capsule = new_capsule(ctypes.addressof(table), name, None)
    return capsule, (
        provider_capsule,
        provider_value,
        provider_derivatives,
        value_callback,
        derivative_callback,
        table,
        name,
    )


def _provider_native(
    model: epcsaft.Mixture,
    spec: dict[str, object],
    source_standard_state: dict[str, object] | None = None,
) -> dict[str, object]:
    native_spec = copy.deepcopy(spec)
    if "conserved_totals" not in native_spec:
        balances = native_spec["balance_matrix"]
        feed = native_spec["feed_amounts"]
        native_spec["conserved_totals"] = tuple(
            math.fsum(
                row[index] * feed[index]  # type: ignore[index]
                for index in range(len(feed))  # type: ignore[arg-type]
            )
            for row in balances  # type: ignore[union-attr]
        )
    return _equilibrium._chemical_equilibrium(
        epcsaft.native_sdk(model),
        native_spec,
        source_standard_state,
        (1.0e-6, 0.74),
        1.0e-12,
    )


def test_iapws_2019_source_record_reproduces_release_equations() -> None:
    source = _water_ionization_source()
    state = source["state"]
    values = source["values"]
    assert isinstance(state, dict)
    assert isinstance(values, dict)
    assert _iapws_p_kw(source, 300.0, 1.0) == pytest.approx(13.906565, abs=5.0e-7)
    p_kw = _iapws_p_kw(
        source,
        state["temperature_k"],
        state["density_kg_per_m3"] / 1000.0,
    )
    assert p_kw == pytest.approx(values["p_kw"], abs=5.0e-13)
    assert -math.log(10.0) * p_kw == pytest.approx(values["ln_kw"], abs=1.0e-13)


def test_held_water_self_ionization_consumes_source_reference_and_provider() -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    capsule = epcsaft.native_sdk(model)

    result = _provider_solve(
        model,
        spec,
        source_standard_state,
        epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(),
    )

    assert result.thermodynamic_model == "installed_provider"
    amounts = result.amounts_mol
    assert 2.0 * amounts[0] + 3.0 * amounts[1] + amounts[2] == pytest.approx(2.0, abs=2.0e-9)
    assert sum(amounts) == pytest.approx(1.0, abs=2.0e-9)
    assert amounts[1] == pytest.approx(amounts[2], abs=2.0e-12)
    final = _equilibrium._chemical_evaluate_provider_block(
        capsule,
        spec["temperature_k"],
        amounts,
        result.volume_m3,
        model.parameter_fingerprint,
    )
    provider_affinity = -2.0 * final["gradient"][0] + final["gradient"][1] + final["gradient"][2]
    assert provider_affinity - result.standard_offsets[0] == pytest.approx(
        spec["ln_k"][0], abs=2.0e-7
    )
    assert result.source_standard_state.id == ("iapws-molality-ions-mole-fraction-water")
    assert result.provider_reference_id == (
        "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3"
    )
    diagnostics = result.diagnostics
    assert diagnostics.reference_derivative_availability == 2
    assert diagnostics.reference_convergence_error <= 5.0e-5
    assert diagnostics.reference_representation_residual_inf_norm <= 1.0e-12
    assert (
        result.source_standard_state.activity_scale_id == source_standard_state["activity_scale_id"]
    )
    assert result.standard_offsets == pytest.approx((-21.200377331250067,), abs=2.0e-12)
    assert result.ln_k_provider_basis == pytest.approx((-53.42391974920632,), abs=2.0e-12)
    assert diagnostics.chemical_certification_level == "LOCAL_EQUILIBRIUM"
    assert diagnostics.boundary_status == "strict_interior"
    assert diagnostics.trace_status == "interior"
    assert diagnostics.provider_domain_status == "passed"
    assert diagnostics.local_minimum_status == "passed"
    assert diagnostics.reaction_affinity_inf_norm <= 1.0e-12
    assert diagnostics.kkt_stationarity_inf_norm <= 1.0e-12
    assert diagnostics.globality_status == "not_guaranteed"
    sensitivity = result.sensitivity
    assert sensitivity is not None
    assert result.response_kind == "value_plus_jacobian"
    assert sensitivity.status == "available"
    assert sensitivity.failure_reason == ""
    assert sensitivity.reference_parameter_status == "available"
    assert tuple(parameter.name for parameter in sensitivity.parameters) == (
        "balance_total[0]",
        "ln_k_provider_basis[0]",
        "pressure_pa",
    )
    pressure_index = tuple(
        parameter.name for parameter in sensitivity.parameters
    ).index("pressure_pa")
    assert sensitivity.amount_derivatives[pressure_index] == pytest.approx(
        (
            6.04641378378937e-18,
            -3.0241437366821673e-18,
            -3.0241437366821673e-18,
        ),
        rel=2.0e-11,
        abs=2.0e-24,
    )
    assert sensitivity.volume_derivatives[pressure_index] == pytest.approx(
        -2.4677833433192133e-7,
        rel=2.0e-11,
        abs=2.0e-15,
    )
    assert result.source_standard_state.reference_pressure_pa == 100_000.0
    artifact = result.artifact_identity
    assert artifact.provider_distribution == "epcsaft"
    assert artifact.provider_record_sha256.startswith("sha256:")
    assert artifact.provider_sdk_capsule_name == "epcsaft.native_sdk.v1"
    assert artifact.provider_sdk_abi_version == 1
    assert artifact.provider_sdk_table_size > 0
    assert artifact.provider_sdk_mixture_result_size > 0
    assert artifact.provider_sdk_neutral_reference_derivative_result_size > 0
    assert artifact.provider_sdk_reacting_phase_parameter_result_size > 0


def test_source_reference_pressure_is_distinct_from_trial_pressure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    source_pressure_pa = 100_000.0
    trial_pressure_pa = 100_010.0
    pressure_step_pa = 2.0
    reference_derivatives = (1.0e-4, -2.0e-4)
    expected_reaction_derivative = (
        -2.0 * reference_derivatives[0] + reference_derivatives[1]
    )
    source_standard_state["id"] = "manufactured-linear-reference-standard"
    source_standard_state["reference_pressure_pa"] = source_pressure_pa
    source_record = spec["equilibrium_constant_records"][0]
    spec["equilibrium_constant_records"] = (
        {
            **source_record,
            "source_id": "manufactured:linear-reference-pressure",
            "reference_id": source_standard_state["id"],
            "pressure_pa": source_pressure_pa,
        },
    )
    spec["ln_k"] = (spec["ln_k"][0] + 50.0,)
    spec["pressure_pa"] = trial_pressure_pa
    capsule, keepalive = _linear_reference_capsule(model, reference_derivatives)
    monkeypatch.setattr(
        epcsaft,
        "native_sdk",
        lambda _model: (keepalive, capsule)[1],
    )

    result = _provider_solve(
        model,
        spec,
        source_standard_state,
        epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(),
    )
    perturbed: list[epcsaft_equilibrium.ChemicalEquilibriumResult] = []
    for direction in (-1.0, 1.0):
        trial = copy.deepcopy(spec)
        trial["pressure_pa"] = trial_pressure_pa + direction * pressure_step_pa
        perturbed.append(_provider_solve(model, trial, source_standard_state))

    assert result.pressure_pa == trial_pressure_pa
    assert result.source_standard_state is not None
    assert result.source_standard_state.reference_pressure_pa == source_pressure_pa
    assert result.sensitivity is not None
    pressure_index = tuple(
        parameter.name for parameter in result.sensitivity.parameters
    ).index("pressure_pa")
    finite_difference_amounts = tuple(
        (
            perturbed[1].amounts_mol[species]
            - perturbed[0].amounts_mol[species]
        )
        / (2.0 * pressure_step_pa)
        for species in range(len(result.amounts_mol))
    )
    assert result.sensitivity.amount_derivatives[pressure_index] == pytest.approx(
        finite_difference_amounts,
        rel=2.0e-5,
        abs=2.0e-10,
    )
    finite_difference_ln_k = (
        perturbed[1].ln_k_provider_basis[0]
        - perturbed[0].ln_k_provider_basis[0]
    ) / (2.0 * pressure_step_pa)
    assert finite_difference_ln_k == pytest.approx(
        expected_reaction_derivative,
        rel=1.0e-11,
        abs=1.0e-14,
    )


def test_source_pressure_derivatives_follow_compiled_reaction_basis() -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    reaction = spec["reaction_matrix"][0]
    source_ln_k = spec["ln_k"][0]
    source_record = spec["equilibrium_constant_records"][0]
    spec["reaction_matrix"] = (
        reaction,
        tuple(2.0 * coefficient for coefficient in reaction),
    )
    spec["ln_k"] = (source_ln_k, 2.0 * source_ln_k)
    spec["equilibrium_constant_records"] = (
        source_record,
        {
            **source_record,
            "source_id": f"{source_record['source_id']}:scaled-by-two",
        },
    )

    result = _provider_solve(
        model,
        spec,
        source_standard_state,
        epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(),
    )

    assert result.sensitivity is not None
    assert result.sensitivity.status == "available"
    assert tuple(parameter.name for parameter in result.sensitivity.parameters) == (
        "balance_total[0]",
        "ln_k_provider_basis[0]",
        "pressure_pa",
    )


def _held_active_parameter(
    family: str = "segment_diameter",
    value: float = 2.2740,
    unit: str = "angstrom",
) -> epcsaft_equilibrium.ChemicalEquilibriumActiveParameter:
    return epcsaft_equilibrium.ChemicalEquilibriumActiveParameter(
        family=family,
        identity="component",
        component_ids=("hydronium-cation",),
        value=value,
        unit=unit,
    )


def test_installed_provider_active_request_is_atomic_exact_and_ordered() -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    value_only = _provider_solve(model, spec, source_standard_state)
    empty_request = _provider_solve(
        model,
        spec,
        source_standard_state,
        epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(),
    )
    assert empty_request.amounts_mol == pytest.approx(
        value_only.amounts_mol, rel=1.0e-9, abs=2.0e-22
    )
    assert empty_request.volume_m3 == pytest.approx(
        value_only.volume_m3, rel=2.0e-13
    )

    parameter = _held_active_parameter()
    request = epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(
        active_parameters=(parameter,)
    )
    result = _provider_solve(model, spec, source_standard_state, request)
    assert result.amounts_mol == pytest.approx(
        value_only.amounts_mol, rel=1.0e-9, abs=2.0e-22
    )
    assert result.volume_m3 == pytest.approx(value_only.volume_m3, rel=2.0e-13)
    assert result.sensitivity is not None
    assert result.sensitivity.status == "available"
    assert result.sensitivity.provider_parameter_status == "available"
    assert result.sensitivity.reference_parameter_status == "available"
    assert tuple(item.name for item in result.sensitivity.parameters) == (
        "balance_total[0]",
        "ln_k_provider_basis[0]",
        "pressure_pa",
        "provider_parameter[segment_diameter;component;hydronium-cation]",
    )

    step = 2.0e-5
    shifted = tuple(
        _provider_solve(
            model,
            spec,
            source_standard_state,
            epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(
                active_parameters=(
                    _held_active_parameter(value=parameter.value + direction * step),
                )
            ),
        )
        for direction in (-1.0, 1.0)
    )
    exact = result.sensitivity.amount_derivatives[-1]
    for ion in (1, 2):
        numerical = (
            shifted[1].amounts_mol[ion] - shifted[0].amounts_mol[ion]
        ) / (2.0 * step)
        assert exact[ion] == pytest.approx(numerical, rel=3.0e-7, abs=1.0e-18)

    dispersion = _held_active_parameter(
        family="dispersion_energy_over_k",
        value=1616.4939,
        unit="kelvin",
    )
    forward = _provider_solve(
        model,
        spec,
        source_standard_state,
        epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(
            active_parameters=(parameter, dispersion)
        ),
    )
    reversed_result = _provider_solve(
        model,
        spec,
        source_standard_state,
        epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(
            active_parameters=(dispersion, parameter)
        ),
    )
    assert forward.sensitivity is not None
    assert reversed_result.sensitivity is not None
    for forward_row, reversed_row in zip(
        forward.sensitivity.amount_derivatives[-2:],
        reversed(reversed_result.sensitivity.amount_derivatives[-2:]),
        strict=True,
    ):
        assert forward_row == pytest.approx(
            reversed_row, rel=2.0e-13, abs=1.0e-25
        )
    assert forward.sensitivity.volume_derivatives[-2:] == pytest.approx(
        tuple(reversed(reversed_result.sensitivity.volume_derivatives[-2:])),
        rel=2.0e-13,
        abs=1.0e-28,
    )


def test_installed_generic_observation_handle_batches_rows_and_exact_columns() -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    standard_state = epcsaft_equilibrium.ChemicalStandardState(
        id=source_standard_state["id"],
        activity_scale_id=source_standard_state["activity_scale_id"],
        log_activity_scale_factors=tuple(
            source_standard_state["log_activity_scale_factors"]
        ),
        reference_pressure_pa=source_standard_state["reference_pressure_pa"],
    )
    problem = _typed_problem(spec, source_standard_state=standard_state)
    phase = epcsaft_equilibrium.ProviderPhase(
        model=model,
        expected_parameter_fingerprint=model.parameter_fingerprint,
        admissible_packing_fraction_interval=(1.0e-6, 0.74),
    )
    state_schema = "fixed_TP_homogeneous_liquid_v1"
    rows = (
        epcsaft_equilibrium.ChemicalObservationRow(
            row_id="held-water-fugacity",
            state_id="held-water-ionization",
            state_schema_id=state_schema,
            source_id="Held2008",
            transform_id="natural_log",
            temperature=spec["temperature_k"] * epcsaft.unit_registry.kelvin,
            pressure=spec["pressure_pa"] * epcsaft.unit_registry.pascal,
            problem=problem,
            primitive=epcsaft_equilibrium.ChemicalObservationPrimitive(
                kind="neutral_component_fugacity_pa",
                component_id="water",
            ),
        ),
        epcsaft_equilibrium.ChemicalObservationRow(
            row_id="held-hydronium-mole-fraction",
            state_id="held-water-ionization",
            state_schema_id=state_schema,
            source_id="Held2008",
            transform_id="natural_log",
            temperature=spec["temperature_k"] * epcsaft.unit_registry.kelvin,
            pressure=spec["pressure_pa"] * epcsaft.unit_registry.pascal,
            problem=problem,
            primitive=epcsaft_equilibrium.ChemicalObservationPrimitive(
                kind="species_mole_fraction",
                component_id="hydronium-cation",
            ),
        ),
    )
    parameter = _held_active_parameter()
    context = epcsaft_equilibrium.chemical_observation_context(
        phase,
        rows=rows,
        active_parameters=(parameter,),
    )

    value_only = context.evaluate((parameter.value,), with_jacobian=False)
    exact = context.evaluate((parameter.value,), with_jacobian=True)

    assert context.row_ids == tuple(row.row_id for row in rows)
    assert value_only["status"] == 0
    assert value_only["values"] == pytest.approx(exact["values"], rel=2.0e-13)
    assert value_only["jacobian"] == []
    assert exact["status"] == 0
    assert len(exact["values"]) == 2
    assert all(math.isfinite(value) and value > 0.0 for value in exact["values"])
    assert len(exact["jacobian"]) == 2
    assert all(math.isfinite(value) for value in exact["jacobian"])
    assert [row["status"] for row in exact["row_results"]] == [0, 0]
    assert [row["solver_status"] for row in exact["row_results"]] == [
        "solve_succeeded",
        "solve_succeeded",
    ]
    assert [row["numerical_status"] for row in exact["row_results"]] == [
        "passed",
        "passed",
    ]
    assert [row["physical_status"] for row in exact["row_results"]] == [
        "passed",
        "passed",
    ]
    assert [row["derivative_status"] for row in exact["row_results"]] == [
        "available",
        "available",
    ]
    assert exact["parameter_ids"] == [
        "segment_diameter;component;hydronium-cation"
    ]
    assert exact["value_only_avoids_derivative_work"] is False
    assert exact["provider_artifact_identity"].startswith("epcsaft==")
    assert exact["owner_artifact_identity"].startswith("epcsaft-equilibrium==")
    assert ";HEADER=sha256:" in exact["provider_artifact_identity"]
    assert ";HEADER=sha256:" in exact["owner_artifact_identity"]
    assert exact["contract_fingerprint"].startswith("sha256:")
    assert exact["capability_fingerprint"].startswith("sha256:")
    assert exact["expected_provider_topology_fingerprint"].startswith("sha256:")
    assert exact["transform_ids"] == ["natural_log", "natural_log"]
    assert all(
        fingerprint.startswith("sha256:")
        for fingerprint in exact["reference_fingerprints"]
    )
    assert exact["artifact_identity"].startswith("sha256:")

    step = 2.0e-3
    shifted = tuple(
        context.evaluate((parameter.value + direction * step,), with_jacobian=False)
        for direction in (-1.0, 1.0)
    )
    finite_difference_log = tuple(
        (
            math.log(shifted[1]["values"][index])
            - math.log(shifted[0]["values"][index])
        )
        / (2.0 * step)
        for index in range(2)
    )
    exact_log = tuple(
        exact["jacobian"][index] / exact["values"][index] for index in range(2)
    )
    assert exact_log == pytest.approx(
        finite_difference_log,
        rel=5.0e-4,
        abs=1.0e-12,
    )

    with pytest.raises(
        ValueError,
        match="descriptor is missing or incompatible",
    ):
        epcsaft_equilibrium.chemical_observation_context(
            phase,
            rows=rows,
            active_parameters=(_held_active_parameter(unit="meter"),),
        )


def test_installed_provider_basis_active_request_has_zero_lnk_cross_block() -> None:
    model, source_spec, source_standard_state = _held_water_ionization_problem()
    transformed = _provider_solve(
        model,
        source_spec,
        source_standard_state,
        epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(),
    )
    spec = copy.deepcopy(source_spec)
    spec["ln_k"] = transformed.ln_k_provider_basis
    spec["equilibrium_constant_records"] = tuple(
        {
            **record,
            "reference_id": "provider-helmholtz-coordinate-basis",
            "conversion_id": "already-provider-basis",
        }
        for record in spec["equilibrium_constant_records"]
    )
    result = _provider_solve(
        model,
        spec,
        None,
        epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(
            active_parameters=(_held_active_parameter(),)
        ),
    )
    assert result.sensitivity is not None
    assert result.sensitivity.status == "available"
    assert result.sensitivity.provider_parameter_status == "available"
    assert result.sensitivity.reference_parameter_status == "not_applicable"


def test_installed_provider_active_request_rejects_unadvertised_unit() -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    with pytest.raises(
        epcsaft_equilibrium.ChemicalEquilibriumError,
        match="descriptor is missing or incompatible",
    ):
        _provider_solve(
            model,
            spec,
            source_standard_state,
            epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(
                active_parameters=(_held_active_parameter(unit="meter"),)
            ),
        )


def test_active_provider_substatus_tracks_kkt_unavailability() -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    balances = spec["balance_matrix"]
    feed = spec["feed_amounts"]
    spec["conserved_totals"] = tuple(
        math.fsum(
            row[index] * feed[index]  # type: ignore[index]
            for index in range(len(feed))  # type: ignore[arg-type]
        )
        for row in balances  # type: ignore[union-attr]
    )
    native = _equilibrium._chemical_equilibrium(
        epcsaft.native_sdk(model),
        spec,
        source_standard_state,
        (1.0e-6, 0.74),
        1.0e-8,
        (
            {
                "family": "segment_diameter",
                "identity": "component",
                "component_ids": ("hydronium-cation",),
                "value": 2.2740,
                "unit": "angstrom",
            },
        ),
    )
    sensitivities = native["sensitivities"]
    assert sensitivities["status"] == "unavailable"
    assert sensitivities["provider_parameter_status"] == "unavailable"
    assert sensitivities["provider_parameter_failure_reason"] == sensitivities["failure_reason"]


@pytest.mark.parametrize("variant", ("species_order", "reaction_orientation"))
def test_held_source_transform_is_coordinate_invariant(variant: str) -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    if variant == "species_order":
        permutation = (2, 0, 1)
        components = tuple(spec["species_ids"][index] for index in permutation)
        parameters = epcsaft.Parameters.from_catalog(
            "held-2008-water-self-ionization",
            components=components,
            version=1,
        )
        model = epcsaft.Mixture(parameters)
        for field in (
            "species_ids",
            "charges",
            "molar_masses_kg_per_mol",
            "feed_amounts",
        ):
            values = spec[field]
            spec[field] = tuple(values[index] for index in permutation)
        spec["balance_matrix"] = tuple(
            tuple(row[index] for index in permutation) for row in spec["balance_matrix"]
        )
        spec["reaction_matrix"] = tuple(
            tuple(row[index] for index in permutation) for row in spec["reaction_matrix"]
        )
        spec["provider_fingerprint"] = model.parameter_fingerprint
        scales = source_standard_state["log_activity_scale_factors"]
        source_standard_state["log_activity_scale_factors"] = tuple(
            scales[index] for index in permutation
        )
        water_index, cation_index, anion_index = (1, 2, 0)
    else:
        spec["reaction_matrix"] = tuple(
            tuple(-value for value in row) for row in spec["reaction_matrix"]
        )
        spec["ln_k"] = tuple(-value for value in spec["ln_k"])
        water_index, cation_index, anion_index = (0, 1, 2)

    result = _provider_solve(model, spec, source_standard_state)

    assert result.amounts_mol[water_index] == pytest.approx(0.9999999999956604, rel=2.0e-10)
    assert result.amounts_mol[cation_index] == pytest.approx(
        result.amounts_mol[anion_index], abs=2.0e-12
    )


@pytest.mark.parametrize(
    "invalid_field",
    ("conversion", "fingerprint", "temperature", "standard_state", "activity_scale"),
)
def test_held_source_transform_fails_closed_on_identity_mismatch(
    invalid_field: str,
) -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    if invalid_field == "conversion":
        record = spec["equilibrium_constant_records"][0]
        spec["equilibrium_constant_records"] = (
            {**record, "conversion_id": "already-provider-basis"},
        )
        match = "provenance"
    elif invalid_field == "fingerprint":
        spec["provider_fingerprint"] = "sha256:wrong"
        match = "fingerprint"
    elif invalid_field == "temperature":
        spec["temperature_k"] = 300.0
        match = "source domain"
    elif invalid_field == "standard_state":
        source_standard_state["id"] = "wrong-standard-state"
        match = "inconsistent"
    else:
        source_standard_state["activity_scale_id"] = ""
        match = "incomplete"

    with pytest.raises(epcsaft_equilibrium.ChemicalEquilibriumError, match=match):
        _provider_solve(model, spec, source_standard_state)


def test_provider_structural_face_fails_before_reduced_topology_evaluation() -> None:
    components = (
        "water",
        "ethanol",
        "isobutanol",
        "sodium-cation",
        "chloride-anion",
    )
    parameters = epcsaft.Parameters.from_catalog(
        "khudaida-2026-figure-2-electrolyte-lle",
        components=components,
        version=1,
    )
    model = epcsaft.Mixture(parameters)
    spec = {
        "species_ids": components,
        "charges": (0, 0, 0, 1, -1),
        "molar_masses_kg_per_mol": (1.0, 1.0, 1.0, 1.0, 1.0),
        "provider_fingerprint": model.parameter_fingerprint,
        "balance_matrix": (
            (1.0, 1.0, 1.0, 1.0, 1.0),
            (0.0, 0.0, 0.0, 1.0, 1.0),
        ),
        "reaction_matrix": (
            (-1.0, 1.0, 0.0, 0.0, 0.0),
            (-1.0, 0.0, 1.0, 0.0, 0.0),
        ),
        "feed_amounts": (1.0, 0.0, 0.0, 0.0, 0.0),
        "ln_k": (0.0, 0.0),
        "temperature_k": 293.15,
        "pressure_pa": 100_000.0,
    }
    _bind_record(spec)

    with pytest.raises(epcsaft_equilibrium.ChemicalEquilibriumError) as failed:
        _provider_solve(model, spec)

    assert failed.value.diagnostics.chemical_certification_level == "BOUNDARY_DIRECTION_UNRESOLVED"
    assert failed.value.diagnostics.structural_zero_species_indices == (3, 4)


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
        "molar_masses_kg_per_mol": (2.0, 1.0, 1.0),
        "provider_fingerprint": model.parameter_fingerprint,
        "balance_matrix": ((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
        "reaction_matrix": ((-1.0, 1.0, 1.0),),
        "feed_amounts": (0.8, 0.1, 0.1),
        "ln_k": (ln_k,),
        "temperature_k": temperature_k,
        "pressure_pa": target["pressure_pa"],
    }
    _bind_record(spec)

    result = _provider_solve(model, spec)

    assert result.thermodynamic_model == "installed_provider"
    final = _equilibrium._chemical_evaluate_provider_block(
        capsule,
        temperature_k,
        result.amounts_mol,
        result.volume_m3,
        model.parameter_fingerprint,
    )
    assert 2.0 * result.amounts_mol[0] + result.amounts_mol[1] + result.amounts_mol[
        2
    ] == pytest.approx(1.8, abs=2.0e-9)
    assert -final["gradient"][0] + final["gradient"][1] + final["gradient"][2] == pytest.approx(
        ln_k, abs=2.0e-7
    )
    assert final["pressure_pa"] == pytest.approx(target["pressure_pa"], rel=1.0e-8)
    assert result.provider_parameter_fingerprint == model.parameter_fingerprint
    assert result.diagnostics.provider_domain_status == "passed"
    assert result.diagnostics.chemical_certification_level == "LOCAL_EQUILIBRIUM"
    assert result.diagnostics.boundary_status == "strict_interior"
    assert 1.0e-6 < result.diagnostics.packing_fraction < 0.74
    assert result.diagnostics.globality_status == "not_guaranteed"


def test_installed_provider_pressure_sensitivity_matches_independent_resolves() -> None:
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
    spec = {
        "species_ids": ("water", "sodium-cation", "chloride-anion"),
        "charges": (0, 1, -1),
        "molar_masses_kg_per_mol": (2.0, 1.0, 1.0),
        "provider_fingerprint": model.parameter_fingerprint,
        "balance_matrix": ((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
        "reaction_matrix": ((-1.0, 1.0, 1.0),),
        "feed_amounts": (0.8, 0.1, 0.1),
        "ln_k": (
            -target["gradient"][0] + target["gradient"][1] + target["gradient"][2],
        ),
        "temperature_k": temperature_k,
        "pressure_pa": target["pressure_pa"],
    }
    _bind_record(spec)
    native = _provider_native(model, spec)
    sensitivities = native["sensitivities"]
    public = _provider_solve(
        model,
        spec,
        sensitivity_request=epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(),
    )
    public_sensitivity = public.sensitivity
    assert public_sensitivity is not None
    assert public_sensitivity.status == "available"
    assert tuple(parameter.name for parameter in public_sensitivity.parameters) == (
        sensitivities["parameter_order"]
    )
    assert public_sensitivity.amount_derivatives == tuple(
        tuple(row) for row in sensitivities["amount_derivatives"]
    )
    assert public_sensitivity.volume_derivatives == pytest.approx(
        sensitivities["volume_derivatives"]
    )
    assert sensitivities["status"] == "available"
    assert sensitivities["parameter_fingerprint"] == model.parameter_fingerprint
    assert sensitivities["provider_parameter_status"] == "not_applicable"
    assert sensitivities["provider_parameter_failure_reason"] == ""
    pressure_index = sensitivities["parameter_order"].index("pressure_pa")

    step_pa = 5.0
    perturbed: list[dict[str, object]] = []
    for direction in (-1.0, 1.0):
        trial = copy.deepcopy(spec)
        trial["pressure_pa"] = spec["pressure_pa"] + direction * step_pa
        _bind_record(trial)
        perturbed.append(_provider_native(model, trial))
    finite_difference_amounts = tuple(
        (
            perturbed[1]["amounts"][species]
            - perturbed[0]["amounts"][species]
        )
        / (2.0 * step_pa)
        for species in range(3)
    )
    finite_difference_volume = (
        perturbed[1]["volume_m3"] - perturbed[0]["volume_m3"]
    ) / (2.0 * step_pa)
    assert sensitivities["amount_derivatives"][pressure_index] == pytest.approx(
        finite_difference_amounts,
        rel=2.0e-5,
        abs=2.0e-12,
    )
    assert sensitivities["volume_derivatives"][pressure_index] == pytest.approx(
        finite_difference_volume,
        rel=2.0e-5,
        abs=2.0e-13,
    )


def test_provider_manufactured_reaction_rejects_capsule_identity_and_source_domain() -> None:
    model = _figiel_provider_model()
    spec = {
        "species_ids": ("water", "sodium-cation", "chloride-anion"),
        "charges": (0, 1, -1),
        "molar_masses_kg_per_mol": (2.0, 1.0, 1.0),
        "provider_fingerprint": model.parameter_fingerprint,
        "balance_matrix": ((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
        "reaction_matrix": ((-1.0, 1.0, 1.0),),
        "feed_amounts": (0.8, 0.1, 0.1),
        "ln_k": (0.0,),
        "temperature_k": 298.15,
        "pressure_pa": 100_000.0,
    }
    _bind_record(spec)
    wrong_order = copy.deepcopy(spec)
    wrong_order["species_ids"] = ("chloride-anion", "sodium-cation", "water")
    wrong_order["charges"] = (-1, 1, 0)
    with pytest.raises(epcsaft_equilibrium.ChemicalEquilibriumError, match="component order"):
        _provider_solve(model, wrong_order)

    fingerprint = copy.deepcopy(spec)
    fingerprint["provider_fingerprint"] = "sha256:wrong"
    with pytest.raises(epcsaft_equilibrium.ChemicalEquilibriumError, match="fingerprint") as failed:
        _provider_solve(model, fingerprint)
    assert "fingerprint" in failed.value.diagnostics.failure_reason
    assert failed.value.diagnostics.numerical_status == "not_adjudicated"
    assert failed.value.diagnostics.balance_inf_norm is None

    outside = copy.deepcopy(spec)
    outside["feed_amounts"] = (0.02, 0.49, 0.49)
    _bind_record(outside)
    with pytest.raises(epcsaft_equilibrium.ChemicalEquilibriumError, match="source domain"):
        _provider_solve(model, outside)
