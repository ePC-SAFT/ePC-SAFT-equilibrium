from __future__ import annotations

from dataclasses import dataclass

import epcsaft
import pytest
from parameter_dictionaries import KHUDAIDA_HELD2_PARAMETERS

import epcsaft_equilibrium

SOLUTION_ATOL = 1.0e-7
SOLUTION_RTOL = 1.0e-7
CONSERVATION_ATOL = 1.0e-10
RESIDUAL_ATOL = 1.0e-8
MAX_PROVIDER_EVALUATIONS = 8_000
MAX_OPTIMIZER_SOLVES = 250
MAX_OPTIMIZER_ITERATIONS = 4_500
KHUDAIDA_PARAMETER_FINGERPRINT = (
    "sha256:cd6ca3c2f91a299a44cf0f57b69d481c4049544a726f301ea1533a134400b961"
)


@dataclass(frozen=True)
class _KhudaidaCase:
    figure: int
    point: int
    temperature_k: float
    feed: tuple[float, ...]
    free_energy_over_rt: float
    phase_fractions: tuple[float, float]
    phase_volumes_m3: tuple[float, float]
    phase_mole_fractions: tuple[tuple[float, ...], tuple[float, ...]]


# Two fastest retained native cases covering the lower- and higher-salt
# 293.15 K families. These checks guarantee numerical continuity only.
KHUDAIDA_CASES = (
    _KhudaidaCase(
        figure=2,
        point=7,
        temperature_k=293.15,
        feed=(
            0.646925133827633,
            0.07045345803743744,
            0.22564721508764596,
            0.028487096523641797,
            0.028487096523641797,
        ),
        free_energy_over_rt=-8.663320942103068,
        phase_fractions=(0.5047970252838168, 0.49520297471618324),
        phase_volumes_m3=(8.767003154558872e-06, 2.8256828355602866e-05),
        phase_mole_fractions=(
            (
                0.8849069190200202,
                0.0022309085496053107,
                0.0007347408661007027,
                0.0560637157821369,
                0.0560637157821369,
            ),
            (
                0.4043326952307936,
                0.13999774956447253,
                0.4549171381960011,
                0.000376208504366339,
                0.000376208504366339,
            ),
        ),
    ),
    _KhudaidaCase(
        figure=5,
        point=5,
        temperature_k=293.15,
        feed=(
            0.6761039014329819,
            0.006153787871534955,
            0.2108900424131834,
            0.053426134141149835,
            0.053426134141149835,
        ),
        free_energy_over_rt=-15.63759219087301,
        phase_fractions=(0.6742281650132003, 0.32577183498679974),
        phase_volumes_m3=(1.1426008799368972e-05, 2.167669749037526e-05),
        phase_mole_fractions=(
            (
                0.8414381724180008,
                0.00011240761666243223,
                0.00046428990040603413,
                0.07899256503246538,
                0.07899256503246538,
            ),
            (
                0.33392262555767094,
                0.01865722827347986,
                0.6463941399173783,
                0.0005130031257354991,
                0.0005130031257354991,
            ),
        ),
    ),
)


@pytest.fixture(scope="module")
def khudaida_model() -> epcsaft.Mixture:
    model = epcsaft.Mixture(epcsaft.Parameters.from_dictionary(KHUDAIDA_HELD2_PARAMETERS))
    assert model.parameter_fingerprint == KHUDAIDA_PARAMETER_FINGERPRINT
    return model


def _aqueous_first(
    result: epcsaft_equilibrium.TpFlashResult,
) -> tuple[tuple[float, ...], tuple[epcsaft_equilibrium.PhaseState, ...]]:
    ordered = sorted(
        zip(result.phase_fractions, result.phases, strict=True),
        key=lambda pair: pair[1].mole_fractions[0],
        reverse=True,
    )
    return tuple(pair[0] for pair in ordered), tuple(pair[1] for pair in ordered)


@pytest.mark.parametrize(
    "case",
    KHUDAIDA_CASES,
    ids=lambda case: f"figure-{case.figure}-point-{case.point}",
)
def test_representative_khudaida_cases_retain_native_solution_and_convergence(
    khudaida_model: epcsaft.Mixture,
    case: _KhudaidaCase,
) -> None:
    result = epcsaft_equilibrium.tp_flash(
        khudaida_model,
        case.temperature_k * epcsaft.unit_registry.kelvin,
        100000.0 * epcsaft.unit_registry.pascal,
        case.feed,
    )
    assert result.diagnostics.outcome == "physical_equilibrium_accepted"
    assert result.diagnostics.solver_status == "passed"
    assert result.diagnostics.numerical_status == "passed"
    assert result.diagnostics.physical_status == "passed"
    assert result.diagnostics.globality_certificate == "not_guaranteed"
    assert result.parameter_fingerprint == KHUDAIDA_PARAMETER_FINGERPRINT
    assert len(result.phases) == 2

    phase_fractions, phases = _aqueous_first(result)
    assert result.overall_mole_fractions == pytest.approx(
        case.feed, abs=SOLUTION_ATOL, rel=SOLUTION_RTOL
    )
    assert result.total_free_energy_over_rt == pytest.approx(
        case.free_energy_over_rt, abs=SOLUTION_ATOL, rel=SOLUTION_RTOL
    )
    assert phase_fractions == pytest.approx(
        case.phase_fractions, abs=SOLUTION_ATOL, rel=SOLUTION_RTOL
    )
    assert tuple(phase.amount_mol for phase in phases) == pytest.approx(
        case.phase_fractions, abs=SOLUTION_ATOL, rel=SOLUTION_RTOL
    )
    assert tuple(phase.volume_m3 for phase in phases) == pytest.approx(
        case.phase_volumes_m3, abs=SOLUTION_ATOL, rel=SOLUTION_RTOL
    )
    for phase, expected in zip(phases, case.phase_mole_fractions, strict=True):
        assert phase.mole_fractions == pytest.approx(expected, abs=SOLUTION_ATOL, rel=SOLUTION_RTOL)

    reconstructed_feed = tuple(
        sum(
            fraction * phase.mole_fractions[index]
            for fraction, phase in zip(phase_fractions, phases, strict=True)
        )
        for index in range(len(case.feed))
    )
    assert (
        max(
            abs(actual - expected)
            for actual, expected in zip(reconstructed_feed, case.feed, strict=True)
        )
        <= CONSERVATION_ATOL
    )
    assert (
        max(abs(phase.mole_fractions[3] - phase.mole_fractions[4]) for phase in phases)
        <= CONSERVATION_ATOL
    )
    assert max(abs(phase.pressure_pa - 100000.0) / 100000.0 for phase in phases) <= RESIDUAL_ATOL
    assert result.diagnostics.kkt_stationarity_max_abs <= RESIDUAL_ATOL
    performance = result.diagnostics.performance
    assert performance is not None
    assert performance.provider_evaluations <= MAX_PROVIDER_EVALUATIONS
    assert performance.optimizer_solves <= MAX_OPTIMIZER_SOLVES
    assert performance.optimizer_iterations <= MAX_OPTIMIZER_ITERATIONS
