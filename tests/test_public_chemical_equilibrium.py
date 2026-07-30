from __future__ import annotations

import inspect
import math

import epcsaft
import pytest

import epcsaft_equilibrium


def _ideal_problem() -> epcsaft_equilibrium.ChemicalEquilibriumProblem:
    return epcsaft_equilibrium.ChemicalEquilibriumProblem(
        species_ids=("A", "B"),
        charges=(0, 0),
        molar_masses_kg_per_mol=(1.0, 1.0),
        balance_matrix=((1.0, 1.0),),
        conserved_totals=(1.0,),
        reaction_matrix=((-1.0, 1.0),),
        feed_amounts_mol=(1.0, 0.0),
        equilibrium_constants=(
            epcsaft_equilibrium.ChemicalEquilibriumConstant(
                ln_value=math.log(4.0),
                source_id="analytic:A-to-B",
                reference_id="provider-helmholtz-coordinate-basis",
                reaction_orientation="products_positive",
                conversion_id="already-provider-basis",
                dimensionless=True,
            ),
        ),
        strict_interior_amount_floor_mol=1.0e-12,
    )


def _ideal_phase() -> epcsaft_equilibrium.IdealGasPhase:
    return epcsaft_equilibrium.IdealGasPhase(
        model_fingerprint="sha256:analytic-a-to-b",
        reference_id="provider-helmholtz-coordinate-basis",
    )


def test_public_chemical_equilibrium_solves_one_typed_ideal_problem() -> None:
    temperature = 350.0 * epcsaft.unit_registry.kelvin
    pressure = 200_000.0 * epcsaft.unit_registry.pascal
    problem = _ideal_problem()

    result = epcsaft_equilibrium.chemical_equilibrium(
        _ideal_phase(),
        temperature,
        pressure,
        problem,
    )
    with_jacobian = epcsaft_equilibrium.chemical_equilibrium(
        _ideal_phase(),
        temperature,
        pressure,
        problem,
        sensitivity_request=epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(),
    )

    assert result.amounts_mol == pytest.approx((0.2, 0.8), rel=2.0e-8)
    assert result.mole_fractions == pytest.approx((0.2, 0.8), rel=2.0e-8)
    assert result.diagnostics.chemical_certification_level == "LOCAL_EQUILIBRIUM"
    assert result.diagnostics.solver_status == "solve_succeeded"
    assert result.diagnostics.numerical_status == "passed"
    assert result.diagnostics.physical_status == "passed"
    assert result.diagnostics.predictive_status == "not_adjudicated"
    assert result.diagnostics.globality_status == "not_guaranteed"
    assert result.diagnostics.search.status == "certified_local_minimum"
    assert result.diagnostics.search.continuation_status == "not_used"
    assert result.diagnostics.search.primary_attempt_count == 5
    assert result.diagnostics.search.selected_objective is not None
    assert result.diagnostics.search.selected_basin_ordinal == 0
    assert len(result.diagnostics.search.attempts) == 5
    assert len(result.diagnostics.search.basins) == 1
    assert result.local_scope == "fixed_TP_single_homogeneous_phase"
    assert result.provider_parameter_fingerprint is None
    assert result.response_kind == "value_only"
    assert result.sensitivity is None
    assert result.artifact_identity.equilibrium_distribution == "epcsaft-equilibrium"
    assert result.artifact_identity.equilibrium_record_sha256.startswith("sha256:")
    assert result.artifact_identity.provider_distribution is None

    assert with_jacobian.amounts_mol == pytest.approx(result.amounts_mol)
    assert with_jacobian.response_kind == "value_plus_jacobian"
    sensitivity = with_jacobian.sensitivity
    assert sensitivity is not None
    assert sensitivity.status == "available"
    assert tuple(parameter.name for parameter in sensitivity.parameters) == (
        "balance_total[0]",
        "ln_k_provider_basis[0]",
        "pressure_pa",
    )
    assert tuple(parameter.input_unit for parameter in sensitivity.parameters) == (
        "mol",
        "dimensionless",
        "Pa",
    )
    assert sensitivity.amount_state_order == ("A", "B")
    assert sensitivity.amount_derivatives[0] == pytest.approx((0.2, 0.8), abs=2.0e-9)
    assert sensitivity.amount_derivatives[1] == pytest.approx((-0.16, 0.16), abs=2.0e-9)
    assert sensitivity.amount_derivatives[2] == pytest.approx((0.0, 0.0), abs=2.0e-12)
    expected_volume_per_mol = 8.31446261815324 * 350.0 / 200_000.0
    assert sensitivity.volume_derivatives == pytest.approx(
        (
            expected_volume_per_mol,
            0.0,
            -expected_volume_per_mol / 200_000.0,
        ),
        rel=2.0e-8,
        abs=2.0e-13,
    )


def test_public_chemical_equilibrium_rejects_incomplete_jacobian_payload(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    native_solve = epcsaft_equilibrium._equilibrium._chemical_equilibrium

    def incomplete_jacobian(*args: object) -> dict[str, object]:
        native = native_solve(*args)
        sensitivities = dict(native["sensitivities"])
        sensitivities["parameter_order"] = sensitivities["parameter_order"][:-1]
        native["sensitivities"] = sensitivities
        return native

    monkeypatch.setattr(
        epcsaft_equilibrium._equilibrium,
        "_chemical_equilibrium",
        incomplete_jacobian,
    )
    with pytest.raises(
        epcsaft_equilibrium.ChemicalEquilibriumError,
        match="parameter and derivative dimensions disagree",
    ):
        epcsaft_equilibrium.chemical_equilibrium(
            _ideal_phase(),
            350.0 * epcsaft.unit_registry.kelvin,
            200_000.0 * epcsaft.unit_registry.pascal,
            _ideal_problem(),
            sensitivity_request=epcsaft_equilibrium.ChemicalEquilibriumSensitivityRequest(),
        )


def test_public_chemical_equilibrium_exports_typed_value_and_jacobian_contract() -> None:
    expected = {
        "ChemicalArtifactIdentity",
        "ChemicalEquilibriumActiveParameter",
        "ChemicalEquilibriumConstant",
        "ChemicalEquilibriumDiagnostics",
        "ChemicalEquilibriumError",
        "ChemicalEquilibriumProblem",
        "ChemicalEquilibriumResult",
        "ChemicalEquilibriumSensitivity",
        "ChemicalEquilibriumSensitivityParameter",
        "ChemicalEquilibriumSensitivityRequest",
        "ChemicalEquilibriumAttempt",
        "ChemicalEquilibriumBasin",
        "ChemicalEquilibriumBudgetPrefix",
        "ChemicalEquilibriumSearch",
        "ChemicalStandardState",
        "IdealGasPhase",
        "ProviderPhase",
        "chemical_equilibrium",
    }
    assert expected <= set(epcsaft_equilibrium.__all__)
    assert tuple(inspect.signature(epcsaft_equilibrium.chemical_equilibrium).parameters) == (
        "phase",
        "temperature",
        "pressure",
        "problem",
        "sensitivity_request",
    )
    for retired in (
        "_chemical_solve_manufactured",
        "_chemical_solve_provider_manufactured",
        "_chemical_solve_provider_source",
    ):
        assert not hasattr(epcsaft_equilibrium._equilibrium, retired)
