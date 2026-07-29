from __future__ import annotations

import epcsaft
import pytest

import epcsaft_equilibrium

COMPONENT_IDS = (
    "water",
    "ethanol",
    "isobutanol",
    "sodium-cation",
    "chloride-anion",
)
PARAMETER_FINGERPRINT = (
    "sha256:b43fac77754d9d5cca8b3db2cbe709892a786d97b756084b167ce126ab4c3007"
)
FIGURE2_TIE_LINE5_FEED = (
    0.7005034356224062,
    0.03286847384962303,
    0.21275572005079824,
    0.026936185238586256,
    0.026936185238586256,
)
EXPECTED_PHASE_FRACTIONS = (
    0.5928509426897598,
    0.4071490573102402,
)
EXPECTED_MOLE_FRACTIONS = (
    (
        0.9065542174872768,
        0.0018825054461472383,
        0.0011236248861048459,
        0.04521982609023558,
        0.04521982609023558,
    ),
    (
        0.40047228467841944,
        0.07798723379349802,
        0.5209138377452583,
        0.00031332189141210004,
        0.00031332189141210004,
    ),
)
EXPECTED_PHASE_VOLUMES_M3 = (
    1.0397037817069201e-5,
    2.4232849287077738e-5,
)


def test_khudaida_figure2_tie_line5_public_numerical_gate() -> None:
    parameters = epcsaft.Parameters.from_catalog(
        "khudaida-2026-figure-2-electrolyte-lle",
        components=COMPONENT_IDS,
        version=1,
    )
    assert parameters.fingerprint == PARAMETER_FINGERPRINT
    model = epcsaft.Mixture(parameters)

    result = epcsaft_equilibrium.tp_flash(
        model,
        293.15 * epcsaft.unit_registry.kelvin,
        100000.0 * epcsaft.unit_registry.pascal,
        FIGURE2_TIE_LINE5_FEED,
    )

    assert result.parameter_fingerprint == PARAMETER_FINGERPRINT
    assert result.overall_mole_fractions == pytest.approx(
        FIGURE2_TIE_LINE5_FEED, abs=1.0e-12
    )
    assert result.diagnostics.outcome == "physical_equilibrium_accepted"
    assert result.diagnostics.solver_status == "passed"
    assert result.diagnostics.numerical_status == "passed"
    assert result.diagnostics.physical_status == "passed"
    assert result.diagnostics.globality_certificate == "not_guaranteed"
    assert result.diagnostics.attempts <= 16
    assert result.diagnostics.material_balance_max_abs is not None
    assert result.diagnostics.material_balance_max_abs <= 1.0e-10
    assert result.diagnostics.pressure_stationarity_max_relative is not None
    assert result.diagnostics.pressure_stationarity_max_relative <= 1.0e-8
    assert result.diagnostics.kkt_stationarity_max_abs is not None
    assert result.diagnostics.kkt_stationarity_max_abs <= 1.0e-8
    assert result.total_free_energy_over_rt == pytest.approx(
        -8.156553792022185, abs=1.0e-10, rel=1.0e-10
    )
    assert result.phase_fractions == pytest.approx(
        EXPECTED_PHASE_FRACTIONS, abs=1.0e-7, rel=1.0e-7
    )
    assert tuple(phase.volume_m3 for phase in result.phases) == pytest.approx(
        EXPECTED_PHASE_VOLUMES_M3, abs=1.0e-10, rel=1.0e-7
    )
    assert len(result.phases) == 2
    for phase, expected in zip(
        result.phases, EXPECTED_MOLE_FRACTIONS, strict=True
    ):
        assert phase.mole_fractions == pytest.approx(
            expected, abs=1.0e-7, rel=1.0e-7
        )
