from __future__ import annotations

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
