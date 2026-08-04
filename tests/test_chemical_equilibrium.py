from __future__ import annotations

import copy
import functools
import hashlib
import json
import math
import os
import subprocess
import sys
import textwrap
import tomllib
from dataclasses import replace
from pathlib import Path

import epcsaft
import numpy as np
import pytest
from chemical_equilibrium_cases import (
    base_system as _base_system,
)
from chemical_equilibrium_cases import (
    bind_records as _bind_record,
)
from chemical_equilibrium_cases import typed_problem as _typed_problem
from parameter_dictionaries import FIGIEL_REFERENCE_ELECTROLYTE_PARAMETERS

import epcsaft_equilibrium
from epcsaft_equilibrium import _api as equilibrium_api
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
    "sha256:8499a0cedeb7e8e34f7d70fbbe4c03180aea0018acae932e831736ce293e2aca"
)
def _packet_fingerprint(packet_root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(item for item in packet_root.rglob("*") if item.is_file()):
        relative = path.relative_to(packet_root).as_posix().encode()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(path.stat().st_size.to_bytes(8, "big"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


@functools.cache
def _held_parameter_packet() -> Path:
    with Path(__file__).with_name("data-lock.toml").open("rb") as stream:
        lock = tomllib.load(stream)["held_cameretti_sadowski_2008"]
    expected_lock = {
        "data_repository": "ePC-SAFT/ePC-SAFT-data",
        "data_commit": "c096285415d4d3198b9d00fc75af48b837dd1305",
        "packet_path": "packets/held-cameretti-sadowski-2008/1",
        "packet_id": "held-cameretti-sadowski-2008",
        "packet_version": "1",
        "packet_fingerprint": ("6850a420751cc323f2700450aef8b900af95aa388b6acb003b1d8eaad04e8dbb"),
        "materialization": "byte-for-byte",
    }
    if lock != expected_lock:
        raise AssertionError("HELD parameter packet lock changed")
    default_root = Path(__file__).resolve().parents[2] / "ePC-SAFT-data"
    data_root = Path(os.environ.get("EPCSAFT_DATA_ROOT", default_root)).resolve()
    data_commit = subprocess.run(
        ["git", "-C", str(data_root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if data_commit != lock["data_commit"]:
        raise AssertionError("HELD parameter packet checkout commit changed")
    packet = data_root / lock["packet_path"]
    if _packet_fingerprint(packet) != lock["packet_fingerprint"]:
        raise AssertionError("HELD parameter packet bytes changed")
    return packet


def _held_parameters(components: tuple[str, ...]) -> epcsaft.Parameters:
    return epcsaft.Parameters.from_bundle(
        _held_parameter_packet() / "parameters",
        components=components,
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


def test_ionic_chart_pullback_spans_complete_physical_reaction_tangent() -> None:
    charges = np.asarray((0.0, 1.0, 1.0, -1.0, -1.0))
    chart = _amount_chart(tuple(charges.astype(int)), (0.2, -0.7, 0.4, -0.3))
    jacobian = np.asarray(chart["jacobian"], dtype=float).reshape(5, 4)
    balance = np.ones((1, 5))
    reactions = np.asarray(
        (
            (-2.0, 1.0, 0.0, 1.0, 0.0),
            (-2.0, 0.0, 1.0, 0.0, 1.0),
            (0.0, 1.0, -1.0, -1.0, 1.0),
        )
    )

    _, _, chart_null_rows = np.linalg.svd(balance @ jacobian)
    chart_feasible = jacobian @ chart_null_rows[1:].T
    chart_basis, _ = np.linalg.qr(chart_feasible)
    reaction_basis, _ = np.linalg.qr(reactions.T)
    principal_cosines = np.linalg.svd(
        chart_basis[:, :3].T @ reaction_basis[:, :3], compute_uv=False
    )

    assert np.linalg.matrix_rank(jacobian) == 4
    assert np.linalg.matrix_rank(chart_feasible) == 3
    assert np.linalg.matrix_rank(reactions) == 3
    assert np.max(np.abs(charges @ jacobian)) <= 2.0e-16
    assert np.max(np.abs(balance @ reactions.T)) <= 2.0e-16
    assert np.max(np.abs(charges @ reactions.T)) <= 2.0e-16
    assert principal_cosines == pytest.approx((1.0, 1.0, 1.0), abs=2.0e-15)

    chemical_potentials = 1.7 * balance[0] - 0.4 * charges
    assert reactions @ chemical_potentials == pytest.approx((0.0, 0.0, 0.0))
    chart_gradient = jacobian.T @ chemical_potentials
    chart_balance_gradient = jacobian.T @ balance[0]
    assert chart_gradient - 1.7 * chart_balance_gradient == pytest.approx(
        (0.0, 0.0, 0.0, 0.0), abs=2.0e-15
    )

    trace_chart = _amount_chart(
        tuple(charges.astype(int)),
        (0.0, -28.0, -27.0, math.log(2.0)),
        trace_floor=1.0e-18,
    )
    trace_jacobian = np.asarray(trace_chart["jacobian"], dtype=float).reshape(5, 4)
    _, _, trace_null_rows = np.linalg.svd(balance @ trace_jacobian)
    trace_pullback = trace_jacobian @ trace_null_rows[1:].T
    trace_singular_values = np.linalg.svd(trace_pullback, compute_uv=False)
    assert trace_singular_values[0] / trace_singular_values[-1] > 1.0e10


def test_incomplete_tangent_can_hide_a_failed_physical_affinity() -> None:
    balance = np.ones(5)
    charges = np.asarray((0.0, 1.0, 1.0, -1.0, -1.0))
    reactions = np.asarray(
        (
            (-2.0, 1.0, 0.0, 1.0, 0.0),
            (-2.0, 0.0, 1.0, 0.0, 1.0),
            (0.0, 1.0, -1.0, -1.0, 1.0),
        )
    )
    incomplete = np.vstack((balance, charges, reactions[:2]))
    _, _, right = np.linalg.svd(incomplete)
    chemical_potentials = right[-1]

    assert incomplete @ chemical_potentials == pytest.approx((0.0, 0.0, 0.0, 0.0), abs=8.0e-16)
    assert abs(reactions[2] @ chemical_potentials) > 0.5


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
    assert max(map(abs, diagnostics.reaction_affinity_residuals)) == pytest.approx(
        diagnostics.reaction_affinity_inf_norm
    )
    assert diagnostics.kkt_stationarity_inf_norm <= 1.0e-7
    assert max(map(abs, diagnostics.physical_stationarity_residuals)) == pytest.approx(
        diagnostics.kkt_stationarity_inf_norm
    )
    physical_rows, physical_columns = diagnostics.physical_constraint_shape
    derivative_dimension = len(diagnostics.derivative_coordinate_order)
    assert physical_columns == len(diagnostics.physical_stationarity_residuals)
    assert len(diagnostics.physical_equality_multipliers) == physical_rows
    assert len(diagnostics.physical_constraint_jacobian) == (physical_rows * physical_columns)
    physical_chart_rows, physical_chart_columns = diagnostics.physical_to_chart_shape
    assert physical_chart_rows == len(diagnostics.physical_lagrangian_gradient)
    assert physical_chart_columns == derivative_dimension
    assert len(diagnostics.physical_to_chart_jacobian) == (
        physical_chart_rows * physical_chart_columns
    )
    assert diagnostics.chart_physical_pullback_residual_inf_norm is not None
    assert diagnostics.chart_physical_pullback_residual_inf_norm <= 1.0e-10
    assert diagnostics.complementarity_inf_norm <= 1.0e-7
    assert diagnostics.objective_gradient_shape == (derivative_dimension,)
    assert diagnostics.constraint_jacobian_shape == (
        len(diagnostics.constraint_values),
        derivative_dimension,
    )
    assert diagnostics.lagrangian_hessian_shape == (
        derivative_dimension,
        derivative_dimension,
    )
    assert diagnostics.derivative_objective_basis == (
        "dimensionless_fixed_TP_A_plus_PV_plus_reference_over_RT"
    )
    assert diagnostics.derivative_constraint_basis == (
        "ordered_per_row_in_derivative_constraint_order"
    )
    assert len(diagnostics.derivative_constraint_order) == len(diagnostics.constraint_values)
    assert diagnostics.derivative_coordinate_order[-1] == "log_volume_m3"
    assert diagnostics.chart_stationarity_inf_norm <= 1.0e-7
    assert diagnostics.derivative_check_step is not None
    assert diagnostics.derivative_check_step > 0.0
    assert diagnostics.objective_gradient_check_relative_error is not None
    assert diagnostics.objective_gradient_check_relative_error <= 1.0e-5
    assert diagnostics.constraint_jacobian_check_relative_error is not None
    assert diagnostics.constraint_jacobian_check_relative_error <= 1.0e-5
    assert diagnostics.lagrangian_hessian_check_relative_error is not None
    assert diagnostics.lagrangian_hessian_check_relative_error <= 2.0e-4
    assert diagnostics.derivative_check_worst_entry
    assert diagnostics.derivative_check_worst_relative_error is not None
    assert diagnostics.derivative_check_worst_analytic_value is not None
    assert diagnostics.derivative_check_worst_finite_difference_value is not None
    assert diagnostics.derivative_check_worst_step is not None
    assert diagnostics.kkt_root_status == "interior_no_active_bounds"
    assert diagnostics.kkt_root_shape == (
        diagnostics.kkt_dimension,
        diagnostics.kkt_dimension,
    )
    assert len(diagnostics.kkt_root_jacobian) == diagnostics.kkt_dimension**2
    assert diagnostics.kkt_root_jacobian_check_relative_error is not None
    assert diagnostics.kkt_root_jacobian_check_relative_error <= 2.0e-4
    assert diagnostics.reduced_hessian_spectrum_status == "converged"
    assert sum(diagnostics.reduced_hessian_raw_inertia) == len(
        diagnostics.reduced_hessian_eigenvalues
    )
    reduced_dimension, full_dimension = diagnostics.reduced_hessian_nullspace_shape
    assert full_dimension == derivative_dimension
    assert reduced_dimension**2 == len(diagnostics.reduced_hessian)
    assert len(diagnostics.reduced_hessian_nullspace_basis) == (reduced_dimension * full_dimension)
    for row in range(reduced_dimension):
        for column in range(reduced_dimension):
            projected = math.fsum(
                diagnostics.reduced_hessian_nullspace_basis[row * full_dimension + left]
                * diagnostics.lagrangian_hessian[left * full_dimension + right]
                * diagnostics.reduced_hessian_nullspace_basis[column * full_dimension + right]
                for left in range(full_dimension)
                for right in range(full_dimension)
            )
            assert projected == pytest.approx(
                diagnostics.reduced_hessian[row * reduced_dimension + column],
                abs=2.0e-10,
            )
    assert diagnostics.minimum_amount_mol is not None
    assert diagnostics.trace_floor_mol == 1.0e-12
    trace_criterion = next(
        criterion
        for criterion in diagnostics.physical_criteria
        if criterion.name == "minimum_amount_mol"
    )
    assert trace_criterion.status == "passed"
    assert trace_criterion.comparison == ">"


def test_manufactured_nonconvex_saddle_recovers_certified_lower_minimum() -> None:
    """A deterministic tangent displacement escapes a stationary saddle only."""

    spec = {**_base_system(), "ln_k": (0.0,)}
    _bind_record(spec)
    native = _equilibrium._chemical_solve_manufactured_nonconvex(
        spec,
        trace_floor=1.0e-12,
        max_iterations=200,
    )

    assert native["accepted"] is True
    assert native["solver_status"] == "solve_succeeded"
    assert native["numerical_status"] == "passed"
    assert native["physical_status"] == "passed"
    assert native["local_minimum_status"] == "passed"
    assert native["negative_curvature_recovery_status"] == "recovered"
    assert native["negative_curvature_recovery_attempts"] == 2
    assert native["negative_curvature_recovery_selected_sign"] in (-1, 1)
    attempts = native["search"]["attempts"]
    assert attempts[0]["kind"] == "primary"
    assert attempts[0]["terminal_status"] == "saddle_observed"
    recovery_attempts = tuple(attempt for attempt in attempts if attempt["kind"] == "recovery")
    assert tuple(attempt["parent_ordinal"] for attempt in recovery_attempts) == (
        0,
        0,
    )
    assert tuple(attempt["start_identity"] for attempt in recovery_attempts) == (
        "negative_curvature;sign=1;seed=0",
        "negative_curvature;sign=-1;seed=0",
    )
    assert all(
        attempt["terminal_status"] == "certified_local_minimum" for attempt in recovery_attempts
    )
    assert native["amounts"][0] != pytest.approx(native["amounts"][1], abs=1.0e-5)
    pressure_over_rt = spec["pressure_pa"] / (8.31446261815324 * spec["temperature_k"])

    def objective(amounts: tuple[float, ...], volume: float) -> float:
        value = math.fsum(amount * (math.log(amount / volume) - 1.0) for amount in amounts)
        difference = amounts[0] - amounts[1]
        value += -2.3 * difference * difference + 2.0 * difference**4
        return value + pressure_over_rt * volume

    saddle_volume = 1.0 / pressure_over_rt
    saddle_coordinates = (
        math.log(0.5),
        math.log(0.5),
        math.log(saddle_volume),
    )
    lower = (
        math.log(0.1e-12),
        math.log(0.1e-12),
        math.log(saddle_volume) - 30.0,
    )
    upper = (0.0, 0.0, math.log(saddle_volume) + 30.0)
    tangent = (1.0 / math.sqrt(2.0), -1.0 / math.sqrt(2.0), 0.0)
    coarse_seed = _equilibrium._chemical_manufactured_recovery_displacement(
        saddle_coordinates, lower, upper, tangent, 1, 0
    )
    backtracked_seed = _equilibrium._chemical_manufactured_recovery_displacement(
        saddle_coordinates, lower, upper, tangent, 1, 3
    )
    assert objective(
        tuple(math.exp(value) for value in coarse_seed[:2]), saddle_volume
    ) > objective((0.5, 0.5), saddle_volume)
    assert objective(
        tuple(math.exp(value) for value in backtracked_seed[:2]), saddle_volume
    ) < objective((0.5, 0.5), saddle_volume)
    assert objective(tuple(native["amounts"]), native["volume_m3"]) < objective(
        (0.5, 0.5), saddle_volume
    )


def test_finite_but_inconsistent_callback_hessian_fails_spanning_audit() -> None:
    spec = {**_base_system(), "ln_k": (0.0,)}
    _bind_record(spec)

    native = _equilibrium._chemical_solve_manufactured_inconsistent_derivatives(
        spec,
        trace_floor=1.0e-12,
        max_iterations=200,
    )

    assert native["accepted"] is False
    assert native["failure_kind"] == "derivative_inconsistency"
    assert native["numerical_status"] == "failed"
    assert "spanning" in native["callback_error"]
    assert native["derivative_check_worst_entry"].startswith("lagrangian_hessian[")
    assert native["derivative_check_worst_analytic_value"] is not None
    assert native["derivative_check_worst_finite_difference_value"] is not None
    hessian_check = next(
        criterion
        for criterion in native["numerical_criteria"]
        if criterion["name"] == "lagrangian_hessian_spanning_relative_error"
    )
    assert hessian_check["status"] == "failed"


def test_manufactured_nonconvex_search_retains_and_selects_distinct_basins() -> None:
    spec = {**_base_system(), "ln_k": (math.log(1.1),)}
    _bind_record(spec)

    native = _equilibrium._chemical_solve_manufactured_nonconvex(
        spec,
        trace_floor=1.0e-12,
        max_iterations=200,
    )

    search = native["search"]
    assert search["status"] == "certified_local_minimum"
    assert search["continuation_status"] == "not_used"
    assert search["primary_budget"] == 25
    assert search["primary_attempt_count"] == 5
    assert search["generated_start_count"] == 5
    assert search["evaluated_start_count"] == 5
    assert search["duplicate_start_count"] == 0
    assert search["infeasible_start_count"] == 0
    assert search["domain_rejected_start_count"] == 0
    assert search["construction_rejected_start_count"] == 0
    assert tuple(prefix["primary_budget"] for prefix in search["budget_prefixes"]) == (
        1,
        5,
    )
    assert len(search["attempts"]) >= search["primary_attempt_count"]
    assert len(search["basins"]) == 2
    certified_attempts = tuple(
        attempt
        for attempt in search["attempts"]
        if attempt["terminal_status"] == "certified_local_minimum"
    )
    assert len(certified_attempts) > len(search["basins"])
    assert len({attempt["basin_ordinal"] for attempt in certified_attempts}) == 2
    objectives = tuple(basin["objective"] for basin in search["basins"])
    assert objectives[0] != pytest.approx(objectives[1], abs=1.0e-10)
    selected = search["selected_basin_ordinal"]
    assert search["selected_objective"] == pytest.approx(min(objectives))
    assert search["basins"][selected]["objective"] == pytest.approx(min(objectives))
    assert all(attempt["terminal_status"] for attempt in search["attempts"])
    repeated = _equilibrium._chemical_solve_manufactured_nonconvex(
        spec,
        trace_floor=1.0e-12,
        max_iterations=200,
    )
    assert repeated["search"] == search


def test_manufactured_search_distinguishes_no_feasible_start_from_infeasibility() -> None:
    spec = _base_system()
    _bind_record(spec)

    native = _equilibrium._chemical_equilibrium(
        None,
        spec,
        None,
        None,
        0.6,
    )

    search = native["search"]
    assert native["accepted"] is False
    assert search["status"] == "no_feasible_start_found"
    assert search["status"] != "infeasible_certified"
    assert search["primary_attempt_count"] == 1
    assert search["attempts"][0]["terminal_status"] == "boundary_unadjudicated"
    assert search["basins"] == []


@pytest.mark.parametrize(
    ("hessian", "status", "inertia"),
    (
        ((2.0, 0.0, 0.0, 1.0), "certified_local_minimum", (2, 0, 0)),
        ((1.0, 0.0, 0.0, 0.0), "second_order_inconclusive", (1, 1, 0)),
        ((1.0, 0.0, 0.0, -0.25), "saddle_observed", (1, 0, 1)),
    ),
)
def test_reduced_hessian_reports_spectrum_and_inertia(
    hessian: tuple[float, ...],
    status: str,
    inertia: tuple[int, int, int],
) -> None:
    evidence = _equilibrium._chemical_analyze_manufactured_reduced_hessian(hessian)

    assert evidence["status"] == status
    assert tuple(evidence["inertia"]) == inertia
    assert tuple(evidence["raw_inertia"]) == inertia
    assert evidence["spectrum_status"] == "converged"
    assert tuple(evidence["nullspace_shape"]) == (2, 2)
    assert tuple(evidence["nullspace_basis"]) == pytest.approx((1.0, 0.0, 0.0, 1.0), abs=2.0e-13)
    assert tuple(evidence["eigenvalues"]) == pytest.approx(tuple(sorted(hessian[::3])), abs=2.0e-13)
    assert evidence["eigenvalue_tolerance"] == pytest.approx(
        1.0e-10 * max(1.0, *(abs(value) for value in hessian[::3]))
    )
    if status == "saddle_observed":
        assert evidence["curvature"] < 0.0
        assert len(evidence["negative_direction"]) == 2


def test_reduced_hessian_separates_raw_spectrum_from_scaled_certification() -> None:
    evidence = _equilibrium._chemical_analyze_manufactured_reduced_hessian((1.0e20, 0.0, 0.0, 1.0))

    assert evidence["status"] == "certified_local_minimum"
    assert tuple(evidence["inertia"]) == (2, 0, 0)
    assert tuple(evidence["raw_inertia"]) == (1, 1, 0)
    assert evidence["spectrum_status"] == "converged"


def test_reduced_hessian_nullspace_is_invariant_to_constraint_row_scaling() -> None:
    hessian = (2.0, 0.25, 0.0, 0.25, 1.5, 0.0, 0.0, 0.0, 0.75)
    baseline = _equilibrium._chemical_analyze_manufactured_reduced_hessian(
        hessian, (1.0, -2.0, 0.5), 1
    )
    scaled = _equilibrium._chemical_analyze_manufactured_reduced_hessian(
        hessian, (1.0e-18, -2.0e-18, 0.5e-18), 1
    )

    assert tuple(baseline["nullspace_shape"]) == (2, 3)
    assert tuple(scaled["nullspace_shape"]) == (2, 3)
    assert scaled["status"] == baseline["status"]
    assert tuple(scaled["inertia"]) == tuple(baseline["inertia"])
    assert tuple(scaled["eigenvalues"]) == pytest.approx(
        tuple(baseline["eigenvalues"]), rel=2.0e-12, abs=2.0e-12
    )


def test_inactive_inequality_must_not_remove_negative_curvature_direction() -> None:
    hessian = (1.0, 0.0, 0.0, -1.0)
    jacobian = (1.0, 0.0, 0.0, 1.0)

    equality_tangent = _equilibrium._chemical_analyze_manufactured_reduced_hessian(
        hessian, jacobian[:2], 1
    )
    incorrectly_active_tangent = _equilibrium._chemical_analyze_manufactured_reduced_hessian(
        hessian, jacobian, 2
    )

    assert equality_tangent["status"] == "saddle_observed"
    assert tuple(equality_tangent["nullspace_shape"]) == (1, 2)
    assert incorrectly_active_tangent["status"] == "certified_local_minimum"
    assert tuple(incorrectly_active_tangent["nullspace_shape"]) == (0, 2)


def test_recovery_displacement_uses_mixed_sign_bound_room() -> None:
    variables = (0.0, 0.0)
    lower = (-1.0, -0.02)
    upper = (1.0, 1.0)
    direction = (1.0, -1.0)

    positive = _equilibrium._chemical_manufactured_recovery_displacement(
        variables, lower, upper, direction, 1
    )
    negative = _equilibrium._chemical_manufactured_recovery_displacement(
        variables, lower, upper, direction, -1
    )

    assert positive == pytest.approx((0.005, -0.005))
    assert negative == pytest.approx((-0.25, 0.25))
    assert all(lower[i] < positive[i] < upper[i] for i in range(2))
    assert all(lower[i] < negative[i] < upper[i] for i in range(2))


def test_balance_retraction_removes_raw_log_seed_balance_error() -> None:
    spec = _base_system()
    _bind_record(spec)
    pressure_over_rt = spec["pressure_pa"] / (8.31446261815324 * spec["temperature_k"])
    saddle_volume = 1.0 / pressure_over_rt
    variables = (math.log(0.5), math.log(0.5), math.log(saddle_volume))
    lower = (math.log(0.1e-12), math.log(0.1e-12), math.log(saddle_volume) - 30.0)
    upper = (0.0, 0.0, math.log(saddle_volume) + 30.0)
    raw = _equilibrium._chemical_manufactured_recovery_displacement(
        variables,
        lower,
        upper,
        (1.0 / math.sqrt(2.0), -1.0 / math.sqrt(2.0), 0.0),
        1,
        0,
    )
    raw_balance_error = abs(math.exp(raw[0]) + math.exp(raw[1]) - 1.0)
    assert raw_balance_error > 1.0e-4

    retracted = _equilibrium._chemical_retract_manufactured_balance(
        spec, raw, lower, upper, 1.0e-12
    )
    assert retracted
    assert abs(math.exp(retracted[0]) + math.exp(retracted[1]) - 1.0) <= 1.0e-9
    assert retracted[2] == variables[2]


def test_failed_balance_retraction_makes_no_recovery_attempt() -> None:
    spec = {**_base_system(), "ln_k": (0.0,)}
    _bind_record(spec)
    native = _equilibrium._chemical_solve_manufactured_nonconvex(
        spec,
        trace_floor=0.49,
        max_iterations=200,
    )

    assert native["accepted"] is False
    assert native["negative_curvature_recovery_status"] == "unresolved"
    assert native["negative_curvature_recovery_attempts"] == 0
    assert native["search"]["status"] == "search_exhausted_no_certified_candidate"
    primary = native["search"]["attempts"][0]
    assert primary["terminal_status"] == "saddle_observed"
    assert primary["recovery_seed_count"] == 8
    assert primary["recovery_solve_count"] == 0


def test_manufactured_search_reports_exhausted_solver_budget() -> None:
    spec = _base_system()
    _bind_record(spec)

    native = _equilibrium._chemical_solve_manufactured_nonconvex(
        spec,
        trace_floor=1.0e-12,
        max_iterations=0,
        quadratic_strength=0.0,
    )

    assert native["accepted"] is False
    assert native["search"]["status"] == "search_exhausted_no_certified_candidate"
    assert native["search"]["basins"] == []
    terminal_statuses = {attempt["terminal_status"] for attempt in native["search"]["attempts"]}
    assert "solver_failed" in terminal_statuses
    assert "certified_local_minimum" not in terminal_statuses


def test_manufactured_ideal_case_does_not_attempt_negative_curvature_recovery() -> None:
    spec = _base_system()
    _bind_record(spec)
    native = _equilibrium._chemical_equilibrium(
        None,
        spec,
        None,
        None,
        1.0e-12,
        None,
    )

    assert native["accepted"] is True
    assert native["negative_curvature_recovery_status"] == "not_needed"
    assert native["negative_curvature_recovery_attempts"] == 0
    assert native["negative_curvature_recovery_selected_sign"] == 0


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
            perturbed.append(_equilibrium._chemical_equilibrium(None, trial, None, None, 1.0e-12))
        finite_difference_amounts = tuple(
            (perturbed[1]["amounts"][species] - perturbed[0]["amounts"][species]) / (2.0 * step)
            for species in range(2)
        )
        finite_difference_volume = (perturbed[1]["volume_m3"] - perturbed[0]["volume_m3"]) / (
            2.0 * step
        )
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

    base_result = _equilibrium._chemical_equilibrium(None, base, None, None, 1.0e-12)[
        "sensitivities"
    ]
    scaled_result = _equilibrium._chemical_equilibrium(None, scaled, None, None, 1.0e-12)[
        "sensitivities"
    ]

    assert base_result["status"] == scaled_result["status"] == "available"
    assert scaled_result["condition_number_inf"] == pytest.approx(
        base_result["condition_number_inf"], rel=2.0e-10
    )
    assert scaled_result["amount_derivatives"][0] == pytest.approx(
        base_result["amount_derivatives"][0], abs=2.0e-9
    )
    assert tuple(
        reaction_scale * value for value in scaled_result["amount_derivatives"][1]
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

    base_result = _equilibrium._chemical_equilibrium(None, base, None, None, 1.0e-12)[
        "sensitivities"
    ]
    permuted_result = _equilibrium._chemical_equilibrium(None, permuted, None, None, 1.0e-12)[
        "sensitivities"
    ]

    assert permuted_result["parameter_order"] == base_result["parameter_order"]
    assert permuted_result["condition_number_inf"] == pytest.approx(
        base_result["condition_number_inf"], rel=2.0e-10
    )
    for parameter in range(3):
        assert permuted_result["amount_derivatives"][parameter] == pytest.approx(
            tuple(base_result["amount_derivatives"][parameter][index] for index in permutation),
            abs=2.0e-9,
        )


def test_manufactured_reaction_with_ill_conditioned_kkt_has_no_sensitivity() -> None:
    balance_separation = 3.0e-6
    spec = _base_system()
    spec["species_ids"] = ("A", "B", "C")
    spec["charges"] = (0, 0, 0)
    spec["molar_masses_kg_per_mol"] = (1.0, 1.0, 1.0)
    spec["balance_matrix"] = ((1.0, 1.0, 1.0 + balance_separation),)
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
    assert native["accepted"] is False
    assert native["local_minimum_status"] == "second_order_inconclusive"
    assert native["reduced_hessian_status"] == "second_order_inconclusive"
    assert native["search"]["status"] == "second_order_inconclusive"
    assert native["failure_kind"] == "ill_conditioning"


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


def test_manufactured_inverse_log_packing_has_exact_directional_pullback() -> None:
    spec = _base_system()
    _bind_record(spec)
    center = (math.log(0.3), math.log(0.7), 0.15)
    direction = (0.2, -0.4, 0.3)
    multipliers = (0.37,)
    step = 2.0e-5

    def evaluate(variables: tuple[float, ...]) -> dict[str, object]:
        return _equilibrium._chemical_evaluate_manufactured_inverse_log_packing_nlp(
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
        rel=3.0e-9,
        abs=3.0e-10,
    )
    constraint_directional = (upper["constraints"][0] - lower["constraints"][0]) / (2.0 * step)
    assert constraint_directional == pytest.approx(
        sum(
            result["constraint_jacobian"][index] * direction[index]
            for index in range(len(direction))
        ),
        rel=3.0e-10,
        abs=3.0e-11,
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
        assert gradient_directional == pytest.approx(hessian_directional, rel=5.0e-9, abs=5.0e-10)
    hessian = result["lagrangian_hessian"]
    for row in range(dimension):
        for column in range(dimension):
            scale = max(1.0, abs(hessian[row * dimension + column]))
            assert (
                abs(hessian[row * dimension + column] - hessian[column * dimension + row])
                <= 2.0e-13 * scale
            )
    assert result["volume_m3"] == pytest.approx(math.exp(-center[-1]), rel=2.0e-15)
    kkt_rhs = result["kkt_backtransform_rhs"]
    kkt_solution = result["kkt_backtransform_solution"]
    kkt_dimension = len(kkt_solution)
    assert kkt_dimension == dimension + 1
    for row in range(dimension):
        reconstructed = (
            sum(
                hessian[row * dimension + column] * kkt_solution[column]
                for column in range(dimension)
            )
            + result["constraint_jacobian"][row] * kkt_solution[-1]
        )
        assert reconstructed == pytest.approx(kkt_rhs[row], rel=2.0e-12, abs=2.0e-13)
    constraint_reconstructed = sum(
        result["constraint_jacobian"][column] * kkt_solution[column] for column in range(dimension)
    )
    assert constraint_reconstructed == pytest.approx(kkt_rhs[-1], rel=2.0e-12, abs=2.0e-13)
    zero_rhs = _equilibrium._chemical_evaluate_manufactured_inverse_log_packing_nlp(
        spec,
        center,
        multipliers,
        zero_kkt_rhs=True,
    )
    assert zero_rhs["kkt_backtransform_rhs"] == [0.0] * kkt_dimension
    assert zero_rhs["kkt_backtransform_solution"] == [0.0] * kkt_dimension


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
    assert diagnostics.search.status == "boundary_unadjudicated"
    assert diagnostics.search.selected_basin_ordinal is None
    assert all(
        attempt.terminal_status == "boundary_unadjudicated"
        for attempt in diagnostics.search.attempts
    )
    native = _equilibrium._chemical_equilibrium(None, trace, None, None, 1.0e-8)
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
    if components != tuple(FIGIEL_REFERENCE_ELECTROLYTE_PARAMETERS["components"]):
        raise ValueError("the compact Figiel test dictionary has one canonical component order")
    parameters = epcsaft.Parameters.from_dictionary(FIGIEL_REFERENCE_ELECTROLYTE_PARAMETERS)
    return epcsaft.Mixture(parameters)


@pytest.mark.parametrize("polar", (False, True), ids=("nonpolar", "full-polar"))
def test_public_provider_phase_block_matches_independent_finite_differences(
    polar: bool,
) -> None:
    mapping: dict[str, object] = {
        "schema": "epcsaft.parameters",
        "schema_version": 1,
        "components": ("solvent", "cation", "anion"),
        "parameters": {
            "mw": (0.018, None, None),
            "m": (1.2, 1.0, 1.0),
            "s": (2.8, 2.8, 2.8),
            "e": (350.0, 230.0, 170.0),
            "z": (0, 1, -1),
            "d_born": (None, 3.4, 4.1),
            "f_solv": (1.5, 1.0, 1.0),
            "epsilon_r": (78.0, None, None),
            "k_ij": (
                (0.0, -0.3, -0.3),
                (-0.3, 0.0, 0.8),
                (-0.3, 0.8, 0.0),
            ),
            "sites": (
                {
                    "component_id": "solvent",
                    "site_id": "a",
                    "site_class": "donor",
                    "multiplicity": 1,
                },
                {
                    "component_id": "solvent",
                    "site_id": "b",
                    "site_class": "acceptor",
                    "multiplicity": 1,
                },
            ),
            "association": (
                {
                    "component_id_a": "solvent",
                    "site_id_a": "a",
                    "component_id_b": "solvent",
                    "site_id_b": "b",
                    "association_energy_over_k": 2425.7,
                    "association_volume": 0.04509,
                },
            ),
        },
        "options": {
            "epsilon_r_ion": 8.0,
            "a_dh": 7.01,
            "permittivity_model": "ion-fraction-suppression",
        },
        "validity": {
            "kind": "reported-conditions",
            "temperature_min_k": 298.0,
            "temperature_max_k": 350.0,
            "pressure_min_pa": 10_000.0,
            "pressure_max_pa": 200_000.0,
            "ion_mole_fraction_max": 0.38,
        },
    }
    if polar:
        parameters = mapping["parameters"]
        options = mapping["options"]
        assert isinstance(parameters, dict)
        assert isinstance(options, dict)
        parameters["dipole_moment"] = (1.85, None, None)
        parameters["quadrupole_moment"] = (4.1, None, None)
        options["polar_formulation"] = "gross-vrabec-point-multipole"
    model = epcsaft.Mixture.from_dictionary(mapping)
    capsule = epcsaft.native_sdk(model)
    temperature_k = 313.15
    point = (0.98, 0.01, 0.01, 1.0e-3)

    def evaluate(values: tuple[float, float, float, float]) -> dict[str, object]:
        return _equilibrium._chemical_evaluate_provider_block(
            capsule,
            temperature_k,
            values[:3],
            values[3],
            model.parameter_fingerprint,
        )

    center = evaluate(point)
    gradient = tuple(center["gradient"])
    hessian = tuple(center["hessian"])
    assert center["component_ids"] == ["solvent", "cation", "anion"]
    assert center["density_transformation"] == "rho_i_mol_per_m3=amount_i_mol/volume_m3"
    assert len(gradient) == 4
    assert len(hessian) == 16
    directions = (
        ((1.0, 0.0, 0.0, 0.0), 2.0e-6),
        # The public electrolyte phase contract admits only electroneutral
        # states, so the ionic derivative is checked on its complete physical
        # tangent rather than through inadmissible single-ion perturbations.
        ((0.0, 1.0, 1.0, 0.0), 2.0e-7),
        ((0.0, 0.0, 0.0, 1.0), 2.0e-8),
    )
    for direction, step in directions:
        lower = tuple(value - step * delta for value, delta in zip(point, direction, strict=True))
        upper = tuple(value + step * delta for value, delta in zip(point, direction, strict=True))
        below = evaluate(lower)
        above = evaluate(upper)
        finite_difference_gradient = (float(above["value"]) - float(below["value"])) / (2.0 * step)
        analytic_directional_gradient = math.fsum(
            value * delta for value, delta in zip(gradient, direction, strict=True)
        )
        assert analytic_directional_gradient == pytest.approx(
            finite_difference_gradient, rel=3.0e-6, abs=3.0e-7
        )
        for row in range(4):
            finite_difference_hessian = (above["gradient"][row] - below["gradient"][row]) / (
                2.0 * step
            )
            analytic_directional_hessian = math.fsum(
                hessian[row * 4 + column] * direction[column] for column in range(4)
            )
            assert analytic_directional_hessian == pytest.approx(
                finite_difference_hessian, rel=3.0e-5, abs=3.0e-5
            )

def test_public_nonpolar_to_polar_recovery_input_preserves_certified_target() -> None:
    def model(polar: bool) -> epcsaft.Mixture:
        parameters: dict[str, object] = {
            "mw": (0.018, None, None),
            "m": (1.2, 1.0, 1.0),
            "s": (2.8, 2.8, 2.8),
            "e": (350.0, 230.0, 170.0),
            "z": (0, 1, -1),
            "d_born": (None, 3.4, 4.1),
            "f_solv": (1.5, 1.0, 1.0),
            "epsilon_r": (78.0, None, None),
            "k_ij": (
                (0.0, -0.3, -0.3),
                (-0.3, 0.0, 0.8),
                (-0.3, 0.8, 0.0),
            ),
            "sites": (
                {
                    "component_id": "solvent",
                    "site_id": "a",
                    "site_class": "donor",
                    "multiplicity": 1,
                },
                {
                    "component_id": "solvent",
                    "site_id": "b",
                    "site_class": "acceptor",
                    "multiplicity": 1,
                },
            ),
            "association": (
                {
                    "component_id_a": "solvent",
                    "site_id_a": "a",
                    "component_id_b": "solvent",
                    "site_id_b": "b",
                    "association_energy_over_k": 2425.7,
                    "association_volume": 0.04509,
                },
            ),
        }
        options: dict[str, object] = {
            "epsilon_r_ion": 8.0,
            "a_dh": 7.01,
            "permittivity_model": "ion-fraction-suppression",
        }
        if polar:
            parameters["dipole_moment"] = (0.2, None, None)
            parameters["quadrupole_moment"] = (0.2, None, None)
            options["polar_formulation"] = "gross-vrabec-point-multipole"
        return epcsaft.Mixture.from_dictionary(
            {
                "schema": "epcsaft.parameters",
                "schema_version": 1,
                "components": ("solvent", "cation", "anion"),
                "parameters": parameters,
                "options": options,
                "validity": {
                    "kind": "reported-conditions",
                    "temperature_min_k": 298.0,
                    "temperature_max_k": 350.0,
                    "pressure_min_pa": 1.0,
                    "pressure_max_pa": 100_000_000.0,
                },
            }
        )

    initial_model = model(False)
    target_model = model(True)
    initial_block = _equilibrium._chemical_evaluate_provider_block(
        epcsaft.native_sdk(initial_model),
        313.15,
        (0.98, 0.01, 0.01),
        0.025,
        initial_model.parameter_fingerprint,
    )
    initial_gradient = tuple(float(value) for value in initial_block["gradient"])
    initial_ln_k = -initial_gradient[0] + initial_gradient[1] + initial_gradient[2]
    target_block = _equilibrium._chemical_evaluate_provider_block(
        epcsaft.native_sdk(target_model),
        313.15,
        (0.98, 0.01, 0.01),
        0.025,
        target_model.parameter_fingerprint,
    )
    target_gradient = tuple(float(value) for value in target_block["gradient"])
    target_ln_k = -target_gradient[0] + target_gradient[1] + target_gradient[2]
    path_pressure_pa = float(initial_block["pressure_pa"])
    initial_phase = epcsaft_equilibrium.ProviderPhase(
        initial_model, initial_model.parameter_fingerprint, (1.0e-6, 0.74)
    )
    target_phase = epcsaft_equilibrium.ProviderPhase(
        target_model, target_model.parameter_fingerprint, (1.0e-6, 0.74)
    )
    problem = epcsaft_equilibrium.ChemicalEquilibriumProblem(
        species_ids=("solvent", "cation", "anion"),
        charges=(0, 1, -1),
        molar_masses_kg_per_mol=(0.018, 0.009, 0.009),
        balance_matrix=((2.0, 1.0, 1.0), (0.0, 1.0, -1.0)),
        conserved_totals=(1.98, 0.0),
        reaction_matrix=((-1.0, 1.0, 1.0),),
        feed_amounts_mol=(0.98, 0.01, 0.01),
        equilibrium_constants=(
            epcsaft_equilibrium.ChemicalEquilibriumConstant(
                ln_value=target_ln_k,
                source_id="synthetic-model-family",
                reference_id="provider-helmholtz-coordinate-basis",
                reaction_orientation="products_positive",
                conversion_id="already-provider-basis",
                dimensionless=True,
            ),
        ),
        strict_interior_amount_floor_mol=1.0e-30,
    )

    with pytest.raises(
        epcsaft_equilibrium.ChemicalEquilibriumError,
        match="endpoint-bound initial equilibrium constants",
    ) as missing_endpoint_constants:
        epcsaft_equilibrium.chemical_equilibrium(
            target_phase,
            313.15 * epcsaft.unit_registry.kelvin,
            path_pressure_pa * epcsaft.unit_registry.pascal,
            problem,
            continuation=epcsaft_equilibrium.ProviderModelContinuation(initial_phase),
        )
    assert missing_endpoint_constants.value.diagnostics.failure_kind == "input_or_native_error"

    initial_constants = (
        epcsaft_equilibrium.ChemicalEquilibriumConstant(
            ln_value=initial_ln_k,
            source_id="synthetic-model-family-initial-endpoint",
            reference_id="provider-helmholtz-coordinate-basis",
            reaction_orientation="products_positive",
            conversion_id="already-provider-basis",
            dimensionless=True,
        ),
    )
    with pytest.raises(
        epcsaft_equilibrium.ChemicalEquilibriumError,
        match="not bound to the initial Provider fingerprint",
    ) as wrong_endpoint_fingerprint:
        epcsaft_equilibrium.chemical_equilibrium(
            target_phase,
            313.15 * epcsaft.unit_registry.kelvin,
            path_pressure_pa * epcsaft.unit_registry.pascal,
            problem,
            continuation=epcsaft_equilibrium.ProviderModelContinuation(
                initial_phase, initial_constants, target_model.parameter_fingerprint
            ),
        )
    assert wrong_endpoint_fingerprint.value.diagnostics.failure_kind == "input_or_native_error"

    incompatible_initial_phase = epcsaft_equilibrium.ProviderPhase(
        initial_model,
        initial_model.parameter_fingerprint,
        (0.80, 0.90),
    )
    with pytest.raises(epcsaft_equilibrium.ChemicalEquilibriumError) as preflight:
        epcsaft_equilibrium.chemical_equilibrium(
            target_phase,
            313.15 * epcsaft.unit_registry.kelvin,
            path_pressure_pa * epcsaft.unit_registry.pascal,
            problem,
            continuation=epcsaft_equilibrium.ProviderModelContinuation(
                incompatible_initial_phase,
                initial_constants,
                initial_model.parameter_fingerprint,
            ),
        )
    assert preflight.value.diagnostics.failure_kind == "physical_domain_failure"
    assert preflight.value.diagnostics.search.continuation_status == "endpoint_domain_incompatible"
    assert preflight.value.diagnostics.search.continuation_blocker

    solved = epcsaft_equilibrium.chemical_equilibrium(
        target_phase,
        313.15 * epcsaft.unit_registry.kelvin,
        path_pressure_pa * epcsaft.unit_registry.pascal,
        problem,
        continuation=epcsaft_equilibrium.ProviderModelContinuation(
            initial_phase,
            initial_constants,
            initial_model.parameter_fingerprint,
        ),
    )

    assert solved.model_fingerprint == target_model.parameter_fingerprint
    assert solved.diagnostics.chemical_certification_level == "LOCAL_EQUILIBRIUM"
    assert solved.diagnostics.search.continuation_status == "not_used"
    assert not solved.diagnostics.search.continuation_initial_model_fingerprint
    assert solved.diagnostics.search.selected_objective == min(
        basin.objective for basin in solved.diagnostics.search.basins
    )
    assert solved.diagnostics.search.continuation_attempt_count == 0
    assert not any(
        attempt.kind == "continuation" for attempt in solved.diagnostics.search.attempts
    )

    record = {
        "source_id": "synthetic-model-family",
        "reference_id": "provider-helmholtz-coordinate-basis",
        "reaction_orientation": "products_positive",
        "conversion_id": "already-provider-basis",
        "dimensionless": True,
        "temperature_k": 313.15,
        "pressure_pa": path_pressure_pa,
    }
    native_spec = {
        "species_ids": problem.species_ids,
        "charges": problem.charges,
        "molar_masses_kg_per_mol": problem.molar_masses_kg_per_mol,
        "provider_fingerprint": target_model.parameter_fingerprint,
        "balance_matrix": problem.balance_matrix,
        "conserved_totals": problem.conserved_totals,
        "reaction_matrix": problem.reaction_matrix,
        "feed_amounts": problem.feed_amounts_mol,
        "ln_k": (target_ln_k,),
        "equilibrium_constant_records": (record,),
        "temperature_k": 313.15,
        "pressure_pa": path_pressure_pa,
    }
    forced = _equilibrium._chemical_equilibrium_continuation(
        epcsaft.native_sdk(target_model),
        epcsaft.native_sdk(initial_model),
        native_spec,
        initial_model.parameter_fingerprint,
        {
            "provider_fingerprint": initial_model.parameter_fingerprint,
            "ln_k": (initial_ln_k,),
            "equilibrium_constant_records": (
                {**record, "source_id": "synthetic-model-family-initial-endpoint"},
            ),
        },
        None,
        (1.0e-6, 0.74),
        (1.0e-6, 0.74),
        problem.strict_interior_amount_floor_mol,
        True,
    )
    forced_diagnostics = equilibrium_api._chemical_diagnostics(forced)
    assert forced_diagnostics.search.continuation_status == "completed"
    assert forced_diagnostics.search.continuation_accepted_lambda == 1.0
    assert forced_diagnostics.search.continuation_attempt_count > 0
    assert any(
        attempt.kind == "continuation"
        and attempt.start_identity == "model_continuation:lambda=1.000000"
        for attempt in forced_diagnostics.search.attempts
    )


def test_public_provider_certifies_multineutral_reactive_minimum() -> None:
    components = ("solvent", "cosolvent", "cation", "anion")
    model = epcsaft.Mixture.from_dictionary(
        {
            "schema": "epcsaft.parameters",
            "schema_version": 1,
            "components": components,
            "parameters": {
                "mw": (0.018, 0.030, None, None),
                "m": (1.2, 1.1, 1.0, 1.0),
                "s": (2.8, 3.0, 2.8, 2.8),
                "e": (350.0, 280.0, 230.0, 170.0),
                "z": (0, 0, 1, -1),
                "d_born": (None, None, 3.4, 4.1),
                "f_solv": (1.5, 1.2, 1.0, 1.0),
                "epsilon_r": (78.0, 40.0, None, None),
                "k_ij": tuple(tuple(0.0 for _ in components) for _ in components),
                "sites": (
                    {
                        "component_id": "solvent",
                        "site_id": "a",
                        "site_class": "donor",
                        "multiplicity": 1,
                    },
                    {
                        "component_id": "solvent",
                        "site_id": "b",
                        "site_class": "acceptor",
                        "multiplicity": 1,
                    },
                ),
                "association": (
                    {
                        "component_id_a": "solvent",
                        "site_id_a": "a",
                        "component_id_b": "solvent",
                        "site_id_b": "b",
                        "association_energy_over_k": 2425.7,
                        "association_volume": 0.04509,
                    },
                ),
            },
            "options": {
                "epsilon_r_ion": 8.0,
                "a_dh": 7.01,
                "permittivity_model": "ion-fraction-suppression",
            },
            "validity": {
                "kind": "reported-conditions",
                "temperature_min_k": 298.0,
                "temperature_max_k": 350.0,
                "pressure_min_pa": 1.0,
                "pressure_max_pa": 100_000_000.0,
            },
        }
    )
    phase = epcsaft_equilibrium.ProviderPhase(
        model, model.parameter_fingerprint, (1.0e-6, 0.74)
    )
    target_amounts = (0.8, 0.2, 0.1, 0.1)
    target_volume = 0.025
    target = _equilibrium._chemical_evaluate_provider_block(
        epcsaft.native_sdk(model),
        313.15,
        target_amounts,
        target_volume,
        model.parameter_fingerprint,
    )
    target_gradient = tuple(float(value) for value in target["gradient"])
    target_ln_k = -target_gradient[0] + target_gradient[2] + target_gradient[3]
    problem = epcsaft_equilibrium.ChemicalEquilibriumProblem(
        species_ids=components,
        charges=(0, 0, 1, -1),
        molar_masses_kg_per_mol=(0.018, 0.030, 0.009, 0.009),
        balance_matrix=(
            (1.0, 0.0, 1.0, 0.0),
            (1.0, 0.0, 0.0, 1.0),
            (0.0, 1.0, 0.0, 0.0),
        ),
        conserved_totals=(0.9, 0.9, 0.2),
        reaction_matrix=((-1.0, 0.0, 1.0, 1.0),),
        feed_amounts_mol=(0.8, 0.2, 0.1, 0.1),
        equilibrium_constants=(
            epcsaft_equilibrium.ChemicalEquilibriumConstant(
                ln_value=target_ln_k,
                source_id="synthetic-multineutral-source",
                reference_id="provider-helmholtz-coordinate-basis",
                reaction_orientation="products_positive",
                conversion_id="already-provider-basis",
                dimensionless=True,
            ),
        ),
        strict_interior_amount_floor_mol=1.0e-30,
    )

    result = epcsaft_equilibrium.chemical_equilibrium(
        phase,
        313.15 * epcsaft.unit_registry.kelvin,
        float(target["pressure_pa"]) * epcsaft.unit_registry.pascal,
        problem,
    )

    assert result.diagnostics.chemical_certification_level == "LOCAL_EQUILIBRIUM"
    assert result.diagnostics.provider_domain_status == "passed"
    assert result.diagnostics.local_minimum_status == "passed"
    assert min(result.amounts_mol) > problem.strict_interior_amount_floor_mol
    assert result.amounts_mol[1] == pytest.approx(0.2, abs=2.0e-10)
    assert result.amounts_mol[2] == pytest.approx(result.amounts_mol[3], abs=2.0e-10)


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
        parameters = epcsaft.Parameters.from_dictionary(
            {FIGIEL_REFERENCE_ELECTROLYTE_PARAMETERS!r}
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
    parameters = _held_parameters(components)
    model = epcsaft.Mixture(parameters)
    assert model.parameter_fingerprint == _HELD_WATER_IONIZATION_FINGERPRINT
    # The installed EOS evaluates the declared pure-solvent molecular mass.
    conversion = math.log(standard_state["standard_molality_mol_per_kg"] * 0.01801528)
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
            "source_reference_component_ids": components,
            "source_reference_solvent_composition": (1.0, 0.0, 0.0),
            "source_reference_activity_convention_id": (
                "molality-infinite-dilution-v1"
            ),
            "source_reference_standard_molality_mol_per_kg": 1.0,
        },
    )


def _provider_solve(
    model: epcsaft.Mixture,
    spec: dict[str, object],
    source_standard_state: dict[str, object] | None = None,
    sensitivity_request: epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest | None = None,
) -> epcsaft_equilibrium.ChemicalEquilibriumResult:
    standard_state = (
        None
        if source_standard_state is None
        else epcsaft_equilibrium.ChemicalStandardState(
            id=source_standard_state["id"],
            activity_scale_id=source_standard_state["activity_scale_id"],
            log_activity_scale_factors=tuple(source_standard_state["log_activity_scale_factors"]),
            reference_pressure_pa=source_standard_state["reference_pressure_pa"],
            source_reference_component_ids=tuple(
                source_standard_state["source_reference_component_ids"]
            ),
            source_reference_solvent_composition=tuple(
                source_standard_state["source_reference_solvent_composition"]
            ),
            source_reference_activity_convention_id=source_standard_state[
                "source_reference_activity_convention_id"
            ],
            source_reference_standard_molality_mol_per_kg=(
                source_standard_state[
                    "source_reference_standard_molality_mol_per_kg"
                ]
            ),
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


def _provider_basis_spec(
    model: epcsaft.Mixture,
    source_spec: dict[str, object],
    source_standard_state: dict[str, object],
) -> tuple[dict[str, object], epcsaft_equilibrium.ChemicalEquilibriumResult]:
    transformed = _provider_solve(model, source_spec, source_standard_state)
    spec = copy.deepcopy(source_spec)
    spec["ln_k"] = transformed.ln_k_provider_basis
    spec["equilibrium_constant_records"] = tuple(
        {
            **record,
            "reference_id": "provider-helmholtz-coordinate-basis",
            "conversion_id": "already-provider-basis",
        }
        for record in spec["equilibrium_constant_records"]  # type: ignore[union-attr]
    )
    return spec, transformed


def _provider_native(
    model: epcsaft.Mixture,
    spec: dict[str, object],
    source_standard_state: dict[str, object] | None = None,
    packing_fraction_bounds: tuple[float, float] = (1.0e-6, 0.74),
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
        packing_fraction_bounds,
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
    assert diagnostics.reference_derivative_availability == 0
    assert diagnostics.reference_convergence_error <= 5.0e-5
    assert diagnostics.reference_representation_residual_inf_norm <= 1.0e-12
    assert (
        result.source_standard_state.activity_scale_id == source_standard_state["activity_scale_id"]
    )
    assert result.standard_offsets == pytest.approx((-21.20037599919826,), abs=2.0e-12)
    assert result.ln_k_provider_basis == pytest.approx((-53.42391841715451,), abs=2.0e-12)
    assert diagnostics.chemical_certification_level == "LOCAL_EQUILIBRIUM"
    assert diagnostics.boundary_status == "strict_interior"
    assert diagnostics.trace_status == "interior"
    assert diagnostics.provider_domain_status == "passed"
    assert diagnostics.local_minimum_status == "passed"
    assert diagnostics.reaction_affinity_inf_norm <= 2.0e-12
    assert diagnostics.kkt_stationarity_inf_norm <= 2.0e-12
    assert diagnostics.globality_status == "not_guaranteed"
    assert result.sensitivity is None
    assert result.response_kind == "value_only"
    transfer = result.source_reference_transfer
    assert transfer is not None
    assert transfer.status == "ok"
    assert transfer.domain_status == "admitted-domain"
    assert transfer.convergence_status == "converged-limit"
    assert transfer.component_ids == tuple(model.component_ids)
    assert transfer.pressure_pa == spec["pressure_pa"]
    assert transfer.source_reference_pressure_pa == 100_000.0
    assert transfer.parameter_fingerprint == model.parameter_fingerprint
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
) -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    source_standard_state["reference_pressure_pa"] = 100_001.0

    value = _provider_solve(model, spec, source_standard_state)

    assert value.pressure_pa == 100_000.0
    assert value.source_reference_transfer is not None
    assert value.source_reference_transfer.pressure_pa == 100_000.0
    assert value.source_standard_state is not None
    assert value.source_standard_state.reference_pressure_pa == 100_001.0
    with pytest.raises(
        epcsaft_equilibrium.ChemicalEquilibriumError,
        match="derivative-unavailable",
    ) as captured:
        _provider_solve(
            model,
            spec,
            source_standard_state,
            epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(),
        )
    assert (
        captured.value.diagnostics.failure_kind
        == "source_reference_derivative_unavailable"
    )


def test_source_transfer_failure_prevents_native_optimization(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    native_called = False

    def reject_transfer(*args: object, **kwargs: object) -> object:
        del args, kwargs
        raise RuntimeError("unsupported-source-reference-transfer")

    def native_sentinel(*args: object, **kwargs: object) -> object:
        del args, kwargs
        nonlocal native_called
        native_called = True
        raise AssertionError("native optimization must not run")

    monkeypatch.setattr(epcsaft, "source_reference_transfer", reject_transfer)
    monkeypatch.setattr(_equilibrium, "_chemical_equilibrium", native_sentinel)

    with pytest.raises(
        epcsaft_equilibrium.ChemicalEquilibriumError,
        match="unsupported-source-reference-transfer",
    ) as captured:
        _provider_solve(model, spec, source_standard_state)
    assert native_called is False
    assert (
        captured.value.diagnostics.failure_kind
        == "unsupported_source_reference_transfer"
    )


def test_source_continuation_evaluates_both_provider_endpoints(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    original = epcsaft.source_reference_transfer
    calls: list[str] = []

    def recording_transfer(*args: object, **kwargs: object) -> object:
        calls.append(args[0].parameter_fingerprint)  # type: ignore[union-attr]
        return original(*args, **kwargs)

    monkeypatch.setattr(epcsaft, "source_reference_transfer", recording_transfer)
    standard_state = epcsaft_equilibrium.ChemicalStandardState(
        id=source_standard_state["id"],
        activity_scale_id=source_standard_state["activity_scale_id"],
        log_activity_scale_factors=tuple(
            source_standard_state["log_activity_scale_factors"]
        ),
        reference_pressure_pa=source_standard_state["reference_pressure_pa"],
        source_reference_component_ids=tuple(
            source_standard_state["source_reference_component_ids"]
        ),
        source_reference_solvent_composition=tuple(
            source_standard_state["source_reference_solvent_composition"]
        ),
        source_reference_activity_convention_id=source_standard_state[
            "source_reference_activity_convention_id"
        ],
        source_reference_standard_molality_mol_per_kg=source_standard_state[
            "source_reference_standard_molality_mol_per_kg"
        ],
    )
    phase = epcsaft_equilibrium.ProviderPhase(
        model=model,
        expected_parameter_fingerprint=model.parameter_fingerprint,
        admissible_packing_fraction_interval=(1.0e-6, 0.74),
    )

    result = epcsaft_equilibrium.chemical_equilibrium(
        phase,
        spec["temperature_k"] * epcsaft.unit_registry.kelvin,
        spec["pressure_pa"] * epcsaft.unit_registry.pascal,
        _typed_problem(spec, source_standard_state=standard_state),
        continuation=epcsaft_equilibrium.ProviderModelContinuation(phase),
    )

    assert calls == [model.parameter_fingerprint, model.parameter_fingerprint]
    assert result.source_reference_transfer is not None
    assert result.diagnostics.chemical_certification_level == "LOCAL_EQUILIBRIUM"


def test_source_transfer_rejects_nonfinite_contractions(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model, spec, source_standard_state = _held_water_ionization_problem()
    original = epcsaft.source_reference_transfer

    def malformed_transfer(*args: object, **kwargs: object) -> object:
        transfer = original(*args, **kwargs)
        return replace(
            transfer,
            activity_scale_log_contractions=(
                math.nan,
                *transfer.activity_scale_log_contractions[1:],
            ),
        )

    monkeypatch.setattr(epcsaft, "source_reference_transfer", malformed_transfer)

    with pytest.raises(
        epcsaft_equilibrium.ChemicalEquilibriumError,
        match="source activity-scale contractions must be finite",
    ) as captured:
        _provider_solve(model, spec, source_standard_state)
    assert (
        captured.value.diagnostics.failure_kind
        == "source_reference_receipt_mismatch"
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
    model, source_spec, source_standard_state = _held_water_ionization_problem()
    spec, value_only = _provider_basis_spec(
        model, source_spec, source_standard_state
    )
    with pytest.raises(
        epcsaft_equilibrium.ChemicalEquilibriumError,
        match="derivative-unavailable",
    ):
        _provider_solve(
            model,
            source_spec,
            source_standard_state,
            epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(),
        )
    empty_request = _provider_solve(
        model,
        spec,
        None,
        epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(),
    )
    assert empty_request.amounts_mol == pytest.approx(
        value_only.amounts_mol, rel=1.0e-9, abs=2.0e-22
    )
    assert empty_request.volume_m3 == pytest.approx(value_only.volume_m3, rel=2.0e-13)

    parameter = _held_active_parameter()
    request = epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(
        active_parameters=(parameter,)
    )
    result = _provider_solve(model, spec, None, request)
    assert result.amounts_mol == pytest.approx(value_only.amounts_mol, rel=1.0e-9, abs=2.0e-22)
    assert result.volume_m3 == pytest.approx(value_only.volume_m3, rel=2.0e-13)
    assert result.sensitivity is not None
    assert result.sensitivity.status == "available"
    assert result.sensitivity.provider_parameter_status == "available"
    assert result.sensitivity.reference_parameter_status == "not_applicable"
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
            None,
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
        numerical = (shifted[1].amounts_mol[ion] - shifted[0].amounts_mol[ion]) / (2.0 * step)
        assert exact[ion] == pytest.approx(numerical, rel=3.0e-7, abs=1.0e-18)

    dispersion = _held_active_parameter(
        family="dispersion_energy_over_k",
        value=1616.4939,
        unit="kelvin",
    )
    forward = _provider_solve(
        model,
        spec,
        None,
        epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(
            active_parameters=(parameter, dispersion)
        ),
    )
    reversed_result = _provider_solve(
        model,
        spec,
        None,
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
        assert forward_row == pytest.approx(reversed_row, rel=2.0e-13, abs=1.0e-25)
    assert forward.sensitivity.volume_derivatives[-2:] == pytest.approx(
        tuple(reversed(reversed_result.sensitivity.volume_derivatives[-2:])),
        rel=2.0e-13,
        abs=1.0e-28,
    )


def test_installed_generic_observation_handle_batches_rows_and_exact_columns() -> None:
    model, source_spec, source_standard_state = _held_water_ionization_problem()
    spec, _ = _provider_basis_spec(model, source_spec, source_standard_state)
    problem = _typed_problem(spec)
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
    assert value_only["values"] == pytest.approx(
        exact["values"], rel=2.0e-13
    ), json.dumps(exact, sort_keys=True)
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
    assert exact["parameter_ids"] == ["segment_diameter;component;hydronium-cation"]
    assert exact["value_only_avoids_derivative_work"] is False
    assert exact["provider_artifact_identity"].startswith("epcsaft==")
    assert exact["owner_artifact_identity"].startswith("epcsaft-equilibrium==")
    assert ";HEADER=sha256:" in exact["provider_artifact_identity"]
    assert ";HEADER=sha256:" in exact["owner_artifact_identity"]
    assert exact["contract_fingerprint"].startswith("sha256:")
    assert exact["capability_fingerprint"].startswith("sha256:")
    assert exact["expected_provider_topology_fingerprint"].startswith("sha256:")
    assert exact["transform_ids"] == ["natural_log", "natural_log"]
    assert all(fingerprint.startswith("sha256:") for fingerprint in exact["reference_fingerprints"])
    assert exact["artifact_identity"].startswith("sha256:")

    step = 2.0e-3
    shifted = tuple(
        context.evaluate((parameter.value + direction * step,), with_jacobian=False)
        for direction in (-1.0, 1.0)
    )
    finite_difference_log = tuple(
        (math.log(shifted[1]["values"][index]) - math.log(shifted[0]["values"][index]))
        / (2.0 * step)
        for index in range(2)
    )
    exact_log = tuple(exact["jacobian"][index] / exact["values"][index] for index in range(2))
    assert exact_log == pytest.approx(
        finite_difference_log,
        rel=2.0e-3,
        abs=3.0e-12,
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
    spec, _ = _provider_basis_spec(
        model, source_spec, source_standard_state
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
    model, source_spec, source_standard_state = _held_water_ionization_problem()
    spec, _ = _provider_basis_spec(model, source_spec, source_standard_state)
    with pytest.raises(
        epcsaft_equilibrium.ChemicalEquilibriumError,
        match="descriptor is missing or incompatible",
    ):
        _provider_solve(
            model,
            spec,
            None,
            epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(
                active_parameters=(_held_active_parameter(unit="meter"),)
            ),
        )


def test_active_provider_substatus_tracks_kkt_unavailability() -> None:
    model, source_spec, source_standard_state = _held_water_ionization_problem()
    spec, _ = _provider_basis_spec(model, source_spec, source_standard_state)
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
        None,
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
        parameters = _held_parameters(components)
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
        reference_ids = source_standard_state["source_reference_component_ids"]
        source_standard_state["source_reference_component_ids"] = tuple(
            reference_ids[index] for index in permutation
        )
        solvent_composition = source_standard_state[
            "source_reference_solvent_composition"
        ]
        source_standard_state["source_reference_solvent_composition"] = tuple(
            solvent_composition[index] for index in permutation
        )
        water_index, cation_index, anion_index = (1, 2, 0)
    else:
        spec["reaction_matrix"] = tuple(
            tuple(-value for value in row) for row in spec["reaction_matrix"]
        )
        spec["ln_k"] = tuple(-value for value in spec["ln_k"])
        water_index, cation_index, anion_index = (0, 1, 2)

    result = _provider_solve(model, spec, source_standard_state)

    assert result.amounts_mol[water_index] == pytest.approx(0.9999999963727326, rel=2.0e-10)
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
        match = "provider-domain-rejection"
    elif invalid_field == "standard_state":
        source_standard_state["id"] = "wrong-standard-state"
        match = "inconsistent"
    else:
        source_standard_state["activity_scale_id"] = ""
        match = "incomplete"

    with pytest.raises(epcsaft_equilibrium.ChemicalEquilibriumError, match=match):
        _provider_solve(model, spec, source_standard_state)


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

    phase = epcsaft_equilibrium.ProviderPhase(
        model=model,
        expected_parameter_fingerprint=model.parameter_fingerprint,
        admissible_packing_fraction_interval=(1.0e-6, 0.74),
    )
    direct_with_recovery_input = epcsaft_equilibrium.chemical_equilibrium(
        phase,
        temperature_k * epcsaft.unit_registry.kelvin,
        target["pressure_pa"] * epcsaft.unit_registry.pascal,
        _typed_problem(spec),
        continuation=epcsaft_equilibrium.ProviderModelContinuation(phase),
    )

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
    assert direct_with_recovery_input.diagnostics.search.continuation_status == "not_used"
    assert direct_with_recovery_input.amounts_mol == pytest.approx(result.amounts_mol)

    rejected = _provider_native(
        model,
        spec,
        packing_fraction_bounds=(0.73, 0.74),
    )
    assert rejected["accepted"] is False
    assert rejected["search"]["status"] == "domain_rejected"
    assert rejected["search"]["primary_attempt_count"] == 1
    assert rejected["search"]["attempts"][0]["start_construction_status"] == "rejected"
    assert rejected["search"]["attempts"][0]["provider_domain_status"] == "failed"
    assert rejected["search"]["attempts"][0]["terminal_status"] == "domain_rejected"


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
        "ln_k": (-target["gradient"][0] + target["gradient"][1] + target["gradient"][2],),
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
    assert (
        tuple(parameter.name for parameter in public_sensitivity.parameters)
        == (sensitivities["parameter_order"])
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
        (perturbed[1]["amounts"][species] - perturbed[0]["amounts"][species]) / (2.0 * step_pa)
        for species in range(3)
    )
    finite_difference_volume = (perturbed[1]["volume_m3"] - perturbed[0]["volume_m3"]) / (
        2.0 * step_pa
    )
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
