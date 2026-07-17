from __future__ import annotations

import epcsaft
import pytest

import epcsaft_equilibrium

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
