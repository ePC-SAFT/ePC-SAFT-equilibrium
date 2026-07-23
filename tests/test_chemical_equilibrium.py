from __future__ import annotations

import copy
import math
from collections.abc import Callable

import epcsaft
import pytest

from epcsaft_equilibrium import _equilibrium


def _base_system() -> dict[str, object]:
    temperature_k = 350.0
    pressure_pa = 200_000.0
    return {
        "species_ids": ("A", "B"),
        "charges": (0, 0),
        "provider_fingerprint": "sha256:manufactured",
        "balance_matrix": ((1.0, 1.0),),
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
            },
        ),
        "temperature_k": temperature_k,
        "pressure_pa": pressure_pa,
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


def _bind_record(spec: dict[str, object]) -> None:
    spec["equilibrium_constant_records"] = tuple(
        {
            "source_id": f"manufactured:reaction-{index}",
            "reference_id": "provider-helmholtz-coordinate-basis",
            "dimensionless": True,
            "temperature_k": spec["temperature_k"],
            "pressure_pa": spec["pressure_pa"],
        }
        for index in range(len(spec["ln_k"]))  # type: ignore[arg-type]
    )


def _nearly_dependent_system(
    reaction_scale: float, delta: float, *, reverse_order: bool = False
) -> dict[str, object]:
    reactions: tuple[tuple[float, ...], ...] = (
        (-1.0, 1.0, 0.0, 0.0),
        (
            -reaction_scale,
            reaction_scale,
            reaction_scale * delta,
            -reaction_scale * delta,
        ),
    )
    ln_k: tuple[float, ...] = (2.0, reaction_scale * (2.0 - delta))
    if reverse_order:
        reactions = tuple(reversed(reactions))
        ln_k = tuple(reversed(ln_k))
    spec = {
        **_base_system(),
        "species_ids": ("A", "B", "C", "D"),
        "charges": (0, 0, 0, 0),
        "balance_matrix": (
            (1.0, 1.0, 1.0, 1.0),
            (1.0, 1.0, 0.0, 0.0),
        ),
        "reaction_matrix": reactions,
        "feed_amounts": (1.0, 1.0, 1.0, 1.0),
        "ln_k": ln_k,
    }
    _bind_record(spec)
    return spec


def test_reaction_compiler_is_stable_for_scaled_nearly_dependent_reactions() -> None:
    delta = 1.0e-6
    expected_reference = (1.0, -1.0, 0.5, -0.5)
    qr_diagonal_ratios = []

    for reaction_scale, reverse_order in (
        (1.0, False),
        (1.0, True),
        (1.0e8, False),
        (1.0e8, True),
    ):
        compiled = _equilibrium._chemical_compile_system(
            _nearly_dependent_system(
                reaction_scale, delta, reverse_order=reverse_order
            )
        )

        assert compiled["reaction_rank"] == 2
        assert compiled["g_ref"] == pytest.approx(expected_reference, abs=2.0e-9)
        assert compiled["reference_reconstruction_inf_norm"] <= 5.0e-8
        qr_diagonal_ratios.append(compiled["reaction_qr_diagonal_ratio"])

    assert qr_diagonal_ratios == pytest.approx(
        (qr_diagonal_ratios[0],) * 4, rel=2.0e-9
    )
    assert 1.0e-8 < qr_diagonal_ratios[0] < 1.0e-4


def test_reaction_compiler_rejects_numerically_dependent_reactions() -> None:
    with pytest.raises(ValueError, match="reaction matrix rank"):
        _equilibrium._chemical_compile_system(_nearly_dependent_system(1.0, 1.0e-14))


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
        (_change("balance_matrix", ((1.0, 1.0), (2.0, 2.0))), "balance matrix rank"),
        (_change("reaction_matrix", ((-1.0, 2.0),)), "conserve"),
        (_change("equilibrium_constant_records.dimensionless", False), "dimensionless"),
    ),
)
def test_reaction_compiler_rejects_inconsistent_contracts(
    mutate: Callable[[dict[str, object]], None], message: str
) -> None:
    spec = copy.deepcopy(_base_system())
    mutate(spec)

    with pytest.raises(ValueError, match=message):
        _equilibrium._chemical_compile_system(spec)


def test_reaction_compiler_rejects_non_neutral_feed_and_charge_nonconservation() -> None:
    non_neutral = _base_system()
    non_neutral.update(
        {
            "charges": (1, -1),
            "feed_amounts": (1.0, 0.0),
        }
    )
    with pytest.raises(ValueError, match="electroneutral"):
        _equilibrium._chemical_compile_system(non_neutral)

    nonconserving = copy.deepcopy(non_neutral)
    nonconserving["feed_amounts"] = (1.0, 1.0)
    with pytest.raises(ValueError, match="charge"):
        _equilibrium._chemical_compile_system(nonconserving)


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
