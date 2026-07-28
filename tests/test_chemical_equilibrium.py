from __future__ import annotations

import copy
import math

import epcsaft
import pytest
from chemical_equilibrium_cases import (
    base_system as _base_system,
)
from chemical_equilibrium_cases import (
    bind_records as _bind_record,
)

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
    concentration_shift = math.log(
        8.31446261815324 * temperature_k / standard_pressure_pa
    )
    manufactured_gibbs = tuple(
        value + concentration_shift for value in _BELOV_SOURCE_GIBBS
    )
    ln_k = tuple(
        -sum(row[index] * manufactured_gibbs[index] for index in range(8))
        for row in reactions
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

    inverse = _equilibrium._chemical_amount_chart_inverse(charges, ionic["amounts"])
    assert inverse == pytest.approx(coordinates, abs=3.0e-15)


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


def test_amount_chart_rejects_one_sided_ionic_topology() -> None:
    with pytest.raises(ValueError, match="both cations and anions"):
        _amount_chart((1, 0), (0.0,))


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
    assert result["accepted"] is True
    assert result["profile"] == "manufactured_ideal_nonpredictive"
    assert result["amounts"] == pytest.approx(expected_amounts, rel=2.0e-8, abs=2.0e-10)
    assert result["volume_m3"] == pytest.approx(volume, rel=2.0e-8)
    assert result["solver_status"] == "solve_succeeded"
    assert result["numerical_status"] == "passed"
    assert result["physical_status"] == "passed"
    assert result["local_minimum_status"] == "passed"
    assert result["trace_status"] == "interior"
    assert result["chemical_certification_level"] == "LOCAL_EQUILIBRIUM"
    assert result["boundary_status"] == "strict_interior"
    assert result["support_qualifiers"] == []
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


@pytest.mark.parametrize("trace_floor", (1.0e-50, 1.0e-55))
def test_belov_aristova_gas_restriction_resolves_extreme_positive_traces(
    trace_floor: float,
) -> None:
    spec, expected = _belov_aristova_gas_system()

    result = _manufactured_solve(spec, {"trace_floor": trace_floor})

    assert result["accepted"] is True
    assert result["trace_status"] == "interior"
    amounts = tuple(result["amounts"])
    assert tuple(map(math.log, amounts)) == pytest.approx(
        tuple(map(math.log, expected)), abs=2.0e-6
    )
    balances = spec["balance_matrix"]
    assert tuple(
        sum(row[index] * amounts[index] for index in range(8))
        for row in balances  # type: ignore[union-attr]
    ) == pytest.approx((2.0, 1.0), abs=1.0e-9)
    assert (
        8.31446261815324 * 2000.0 * sum(amounts) / result["volume_m3"]
    ) == pytest.approx(100_000.0, rel=1.0e-8)
    concentration_shift = math.log(8.31446261815324 * 2000.0 / 101_325.0)
    potentials = tuple(
        _BELOV_SOURCE_GIBBS[index] + concentration_shift
        + math.log(amounts[index] / result["volume_m3"])
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
    assert result["balance_inf_norm"] <= 1.0e-9
    assert result["pressure_relative_residual"] <= 1.0e-8
    assert result["reaction_affinity_inf_norm"] <= 1.0e-7
    assert result["kkt_stationarity_inf_norm"] <= 1.0e-7
    assert result["local_minimum_status"] == "passed"
    assert result["chemical_certification_level"] == "LOCAL_EQUILIBRIUM"
    assert result["boundary_status"] == "strict_interior"
    assert result["support_qualifiers"] == []


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
        options["gauge_coefficients"] = (3.25, -1.75)
    else:
        scale = 3.0
        spec["feed_amounts"] = tuple(
            scale * value for value in spec["feed_amounts"]  # type: ignore[union-attr]
        )
        expected = tuple(scale * value for value in expected)

    result = _manufactured_solve(spec, options)

    assert result["accepted"] is True
    assert tuple(map(math.log, result["amounts"])) == pytest.approx(
        tuple(map(math.log, expected)), abs=3.0e-6
    )
    assert result["reaction_affinity_inf_norm"] <= 1.0e-7


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

    assert result["accepted"] is True
    assert result["amounts"] == pytest.approx((0.5, 0.5, 0.5), rel=3.0e-8)
    assert result["volume_m3"] == pytest.approx(1.5, rel=3.0e-8)
    assert result["charge_inf_norm"] <= 2.0e-15


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
    assert indeterminate["chemical_certification_level"] == "FEASIBLE_ONLY"
    assert indeterminate["final_lambda"] is None

    trace = copy.deepcopy(spec)
    trace["ln_k"] = (math.log(1.0e-10),)
    _bind_record(trace)
    false_success = _manufactured_solve(trace, {"trace_floor": 1.0e-8})
    assert false_success["accepted"] is False
    assert false_success["physical_status"] == "failed"
    assert false_success["trace_status"] == "at_or_below_floor"
    assert false_success["chemical_certification_level"] == "FEASIBLE_ONLY"
    assert false_success["globality_certificate"] == "not_guaranteed"


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
    assert tuple(reordered["amounts"][index] for index in inverse) == pytest.approx(
        baseline["amounts"], rel=4.0e-8
    )


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
    assert result["chemical_certification_level"] == "LOCAL_EQUILIBRIUM"
    assert result["boundary_status"] == "strict_interior"
    assert result["support_qualifiers"] == []
    assert result["packing_fraction_bounds"] == pytest.approx((1.0e-6, 0.74))
    assert 1.0e-6 < result["packing_fraction"] < 0.74
    assert result["predictive_status"] == "manufactured_nonpredictive"
    assert result["globality_certificate"] == "not_guaranteed"


def test_provider_manufactured_reaction_rejects_capsule_identity_and_source_domain() -> None:
    model = _figiel_provider_model()
    capsule = epcsaft.native_sdk(model)
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
    with pytest.raises(ValueError, match="component order"):
        _equilibrium._chemical_solve_provider_manufactured(capsule, wrong_order, {})

    fingerprint = copy.deepcopy(spec)
    fingerprint["provider_fingerprint"] = "sha256:wrong"
    with pytest.raises(ValueError, match="fingerprint"):
        _equilibrium._chemical_solve_provider_manufactured(
            capsule, fingerprint, {"packing_fraction_bounds": (1.0e-6, 0.74)}
        )

    outside = copy.deepcopy(spec)
    outside["feed_amounts"] = (0.02, 0.49, 0.49)
    with pytest.raises(ValueError, match="source domain"):
        _equilibrium._chemical_solve_provider_manufactured(capsule, outside, {})
