from __future__ import annotations

import inspect
import math

import epcsaft
import pytest

import epcsaft_equilibrium


def test_public_chemical_equilibrium_solves_one_typed_ideal_problem() -> None:
    temperature = 350.0 * epcsaft.unit_registry.kelvin
    pressure = 200_000.0 * epcsaft.unit_registry.pascal
    problem = epcsaft_equilibrium.ChemicalEquilibriumProblem(
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

    result = epcsaft_equilibrium.chemical_equilibrium(
        epcsaft_equilibrium.IdealGasPhase(
            model_fingerprint="sha256:analytic-a-to-b",
            reference_id="provider-helmholtz-coordinate-basis",
        ),
        temperature,
        pressure,
        problem,
    )

    assert result.amounts_mol == pytest.approx((0.2, 0.8), rel=2.0e-8)
    assert result.mole_fractions == pytest.approx((0.2, 0.8), rel=2.0e-8)
    assert result.diagnostics.chemical_certification_level == "LOCAL_EQUILIBRIUM"
    assert result.diagnostics.solver_status == "solve_succeeded"
    assert result.diagnostics.numerical_status == "passed"
    assert result.diagnostics.physical_status == "passed"
    assert result.diagnostics.predictive_status == "not_adjudicated"
    assert result.diagnostics.globality_status == "not_guaranteed"
    assert result.local_scope == "fixed_TP_single_homogeneous_phase"
    assert result.provider_parameter_fingerprint is None


def test_public_chemical_equilibrium_exports_only_typed_value_contract() -> None:
    expected = {
        "ChemicalEquilibriumConstant",
        "ChemicalEquilibriumDiagnostics",
        "ChemicalEquilibriumError",
        "ChemicalEquilibriumProblem",
        "ChemicalEquilibriumResult",
        "ChemicalStandardState",
        "IdealGasPhase",
        "ProviderPhase",
        "chemical_equilibrium",
    }
    assert expected <= set(epcsaft_equilibrium.__all__)
    assert tuple(
        inspect.signature(epcsaft_equilibrium.chemical_equilibrium).parameters
    ) == ("phase", "temperature", "pressure", "problem")
    for retired in (
        "_chemical_solve_manufactured",
        "_chemical_solve_provider_manufactured",
        "_chemical_solve_provider_source",
    ):
        assert not hasattr(epcsaft_equilibrium._equilibrium, retired)
