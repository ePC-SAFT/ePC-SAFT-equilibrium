from __future__ import annotations

import copy
import inspect
import math
from dataclasses import replace

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
    assert result.diagnostics.kkt_root_status == "interior_no_active_bounds"
    assert result.diagnostics.first_failed_numerical_criterion is None
    assert tuple(
        criterion.name for criterion in result.diagnostics.numerical_criteria
    )[-2:] == ("chart_stationarity_inf_norm", "derivative_evidence_finite")
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


def test_chart_multiplier_reconstruction_is_invariant_to_balance_row_scale() -> None:
    baseline_problem = epcsaft_equilibrium.ChemicalEquilibriumProblem(
        species_ids=("A", "B", "C"),
        charges=(0, 0, 0),
        molar_masses_kg_per_mol=(1.0, 1.0, 2.0),
        balance_matrix=((1.0, 0.0, 1.0),),
        conserved_totals=(2.0,),
        reaction_matrix=((-1.0, -1.0, 1.0),),
        feed_amounts_mol=(2.0, 1.0, 0.0),
        equilibrium_constants=(
            epcsaft_equilibrium.ChemicalEquilibriumConstant(
                ln_value=0.0,
                source_id="analytic:a-plus-b-to-c",
                reference_id="provider-helmholtz-coordinate-basis",
                reaction_orientation="products_positive",
                conversion_id="already-provider-basis",
                dimensionless=True,
            ),
        ),
        strict_interior_amount_floor_mol=1.0e-12,
    )
    scaled_problem = replace(
        baseline_problem,
        balance_matrix=((1.0e-18, 0.0, 1.0e-18),),
        conserved_totals=(2.0e-18,),
    )
    baseline = epcsaft_equilibrium.chemical_equilibrium(
        _ideal_phase(),
        350.0 * epcsaft.unit_registry.kelvin,
        200_000.0 * epcsaft.unit_registry.pascal,
        baseline_problem,
    )
    scaled = epcsaft_equilibrium.chemical_equilibrium(
        _ideal_phase(),
        350.0 * epcsaft.unit_registry.kelvin,
        200_000.0 * epcsaft.unit_registry.pascal,
        scaled_problem,
    )

    assert scaled.amounts_mol == pytest.approx(baseline.amounts_mol, abs=2.0e-11)
    assert scaled.diagnostics.local_minimum_status == "passed"
    assert scaled.diagnostics.chart_stationarity_inf_norm is not None
    assert scaled.diagnostics.chart_stationarity_inf_norm <= 1.0e-7
    assert scaled.diagnostics.kkt_stationarity_inf_norm is not None
    assert scaled.diagnostics.kkt_stationarity_inf_norm <= 1.0e-7


def test_reaction_compilation_is_invariant_to_row_scale() -> None:
    baseline_problem = _ideal_problem()
    reaction = baseline_problem.equilibrium_constants[0]
    scale = 1.0e-12
    scaled_problem = replace(
        baseline_problem,
        reaction_matrix=((-scale, scale),),
        equilibrium_constants=(
            replace(reaction, ln_value=scale * reaction.ln_value),
        ),
    )

    baseline = epcsaft_equilibrium.chemical_equilibrium(
        _ideal_phase(),
        350.0 * epcsaft.unit_registry.kelvin,
        200_000.0 * epcsaft.unit_registry.pascal,
        baseline_problem,
    )
    scaled = epcsaft_equilibrium.chemical_equilibrium(
        _ideal_phase(),
        350.0 * epcsaft.unit_registry.kelvin,
        200_000.0 * epcsaft.unit_registry.pascal,
        scaled_problem,
    )

    assert scaled.amounts_mol == pytest.approx(baseline.amounts_mol, abs=2.0e-11)
    assert scaled.diagnostics.local_minimum_status == "passed"


def test_scaled_balance_row_rejects_inconsistent_declared_total() -> None:
    problem = epcsaft_equilibrium.ChemicalEquilibriumProblem(
        species_ids=("A", "B", "C"),
        charges=(0, 0, 0),
        molar_masses_kg_per_mol=(1.0, 1.0, 2.0),
        balance_matrix=((1.0e-18, 0.0, 1.0e-18),),
        conserved_totals=(1.0e-13,),
        reaction_matrix=((-1.0, -1.0, 1.0),),
        feed_amounts_mol=(2.0, 1.0, 0.0),
        equilibrium_constants=(
            epcsaft_equilibrium.ChemicalEquilibriumConstant(
                ln_value=0.0,
                source_id="analytic:a-plus-b-to-c",
                reference_id="provider-helmholtz-coordinate-basis",
                reaction_orientation="products_positive",
                conversion_id="already-provider-basis",
                dimensionless=True,
            ),
        ),
        strict_interior_amount_floor_mol=1.0e-12,
    )

    with pytest.raises(
        epcsaft_equilibrium.ChemicalEquilibriumError,
        match="conserved totals do not match",
    ):
        epcsaft_equilibrium.chemical_equilibrium(
            _ideal_phase(),
            350.0 * epcsaft.unit_registry.kelvin,
            200_000.0 * epcsaft.unit_registry.pascal,
            problem,
        )


def test_scaled_dependent_reaction_rejects_inconsistent_constant_cycle() -> None:
    reaction = _ideal_problem().equilibrium_constants[0]
    scale = 1.0e-12
    problem = replace(
        _ideal_problem(),
        reaction_matrix=((-1.0, 1.0), (-scale, scale)),
        equilibrium_constants=(
            reaction,
            replace(
                reaction,
                ln_value=scale * reaction.ln_value + 1.0e-12,
                source_id="analytic:inconsistent-scaled-cycle",
            ),
        ),
    )

    with pytest.raises(
        epcsaft_equilibrium.ChemicalEquilibriumError,
        match="reaction constant cycle is inconsistent",
    ):
        epcsaft_equilibrium.chemical_equilibrium(
            _ideal_phase(),
            350.0 * epcsaft.unit_registry.kelvin,
            200_000.0 * epcsaft.unit_registry.pascal,
            problem,
        )


@pytest.mark.parametrize(
    ("charges", "feed", "reaction", "message"),
    (
        (
            (0, 0),
            (1.0, 0.0),
            (-1.0e-12, 2.0e-12),
            "reaction stoichiometry does not conserve",
        ),
        (
            (1, -1),
            (1.0, 1.0),
            (1.0e-12, -1.0e-12),
            "reaction stoichiometry does not conserve charge",
        ),
    ),
)
def test_scaled_invalid_reaction_contract_is_rejected_at_compilation(
    charges: tuple[int, int],
    feed: tuple[float, float],
    reaction: tuple[float, float],
    message: str,
) -> None:
    problem = replace(
        _ideal_problem(),
        charges=charges,
        feed_amounts_mol=feed,
        conserved_totals=(sum(feed),),
        reaction_matrix=(reaction,),
        equilibrium_constants=(
            replace(_ideal_problem().equilibrium_constants[0], ln_value=0.0),
        ),
    )

    with pytest.raises(epcsaft_equilibrium.ChemicalEquilibriumError, match=message):
        epcsaft_equilibrium.chemical_equilibrium(
            _ideal_phase(),
            350.0 * epcsaft.unit_registry.kelvin,
            200_000.0 * epcsaft.unit_registry.pascal,
            problem,
        )


def test_structural_face_reduction_is_invariant_to_reaction_row_scale() -> None:
    def problem(scale: float) -> epcsaft_equilibrium.ChemicalEquilibriumProblem:
        constants = _ideal_problem().equilibrium_constants[0]
        return epcsaft_equilibrium.ChemicalEquilibriumProblem(
            species_ids=("A", "B", "C", "D"),
            charges=(0, 0, 0, 0),
            molar_masses_kg_per_mol=(1.0, 1.0, 1.0, 1.0),
            balance_matrix=((0.0, 0.0, 1.0, 1.0),),
            conserved_totals=(0.0,),
            reaction_matrix=(
                (-scale, scale, -scale, scale),
                (0.0, 0.0, -scale, scale),
            ),
            feed_amounts_mol=(1.0, 0.0, 0.0, 0.0),
            equilibrium_constants=(
                replace(constants, ln_value=scale * math.log(4.0)),
                replace(
                    constants,
                    ln_value=0.0,
                    source_id="analytic:c-to-d",
                ),
            ),
            strict_interior_amount_floor_mol=1.0e-12,
        )

    baseline = epcsaft_equilibrium.chemical_equilibrium(
        _ideal_phase(),
        350.0 * epcsaft.unit_registry.kelvin,
        200_000.0 * epcsaft.unit_registry.pascal,
        problem(1.0),
    )
    scaled = epcsaft_equilibrium.chemical_equilibrium(
        _ideal_phase(),
        350.0 * epcsaft.unit_registry.kelvin,
        200_000.0 * epcsaft.unit_registry.pascal,
        problem(1.0e-12),
    )

    assert baseline.amounts_mol == pytest.approx((0.2, 0.8, 0.0, 0.0), abs=2.0e-10)
    assert scaled.amounts_mol == pytest.approx(baseline.amounts_mol, abs=2.0e-10)
    assert scaled.diagnostics.boundary_status == "structural_face"
    assert scaled.diagnostics.chemical_certification_level == "LOCAL_EQUILIBRIUM"


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


@pytest.mark.parametrize(
    ("failure_kind", "failure_reason"),
    (
        (
            "exhausted_multistart_search",
            "exhausted_multistart_search: no certified local minimum was observed",
        ),
        (
            "derivative_inconsistency",
            "derivative inconsistency: Lagrangian Hessian check failed",
        ),
    ),
)
def test_public_failure_exposes_typed_criterion_and_search_blocker(
    monkeypatch: pytest.MonkeyPatch,
    failure_kind: str,
    failure_reason: str,
) -> None:
    native_solve = epcsaft_equilibrium._equilibrium._chemical_equilibrium

    def failed_certificate(*args: object) -> dict[str, object]:
        native = native_solve(*args)
        native["accepted"] = False
        native["chemical_certification_level"] = "FEASIBLE_ONLY"
        native["failure_kind"] = failure_kind
        native["failure_reason"] = failure_reason
        native["numerical_criteria"][2]["value"] = 56.0
        native["numerical_criteria"][2]["status"] = "failed"
        native["numerical_status"] = "failed"
        native["physical_criteria"][2]["value"] = 2.0e-7
        native["physical_criteria"][2]["status"] = "failed"
        native["physical_status"] = "failed"
        return native

    monkeypatch.setattr(
        epcsaft_equilibrium._equilibrium,
        "_chemical_equilibrium",
        failed_certificate,
    )
    with pytest.raises(epcsaft_equilibrium.ChemicalEquilibriumError) as failed:
        epcsaft_equilibrium.chemical_equilibrium(
            _ideal_phase(),
            350.0 * epcsaft.unit_registry.kelvin,
            200_000.0 * epcsaft.unit_registry.pascal,
            _ideal_problem(),
        )

    diagnostics = failed.value.diagnostics
    assert diagnostics.failure_kind == failure_kind
    assert diagnostics.failure_reason == failure_reason
    assert diagnostics.first_failed_numerical_criterion == "kkt_stationarity_inf_norm"
    assert diagnostics.first_failed_physical_criterion == "pressure_relative_residual"
    assert diagnostics.numerical_criteria[2].value == 56.0
    assert diagnostics.numerical_criteria[2].limit == 1.0e-7
    assert diagnostics.physical_criteria[2].value == 2.0e-7
    assert diagnostics.physical_criteria[2].limit == 1.0e-8


def test_public_nonfinite_reduced_hessian_preserves_derivative_blocker(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    native_solve = epcsaft_equilibrium._equilibrium._chemical_equilibrium

    def nonfinite_hessian(*args: object) -> dict[str, object]:
        native = native_solve(*args)
        native["accepted"] = False
        native["chemical_certification_level"] = "FEASIBLE_ONLY"
        native["failure_kind"] = "derivative_inconsistency"
        native["failure_reason"] = (
            "derivative inconsistency: reduced Lagrangian Hessian spectrum is nonfinite"
        )
        native["callback_error"] = native["failure_reason"]
        native["numerical_status"] = "failed"
        native["numerical_criteria"][4]["value"] = 0.0
        native["numerical_criteria"][4]["status"] = "failed"
        native["reduced_hessian_spectrum_status"] = "nonfinite"
        native["reduced_hessian_eigenvalues"] = ()
        native["reduced_hessian_raw_inertia"] = (0, 0, 0)
        native["reduced_hessian_inertia"] = (0, 0, 0)
        return native

    monkeypatch.setattr(
        epcsaft_equilibrium._equilibrium,
        "_chemical_equilibrium",
        nonfinite_hessian,
    )
    with pytest.raises(epcsaft_equilibrium.ChemicalEquilibriumError) as failed:
        epcsaft_equilibrium.chemical_equilibrium(
            _ideal_phase(),
            350.0 * epcsaft.unit_registry.kelvin,
            200_000.0 * epcsaft.unit_registry.pascal,
            _ideal_problem(),
        )

    assert failed.value.diagnostics.failure_kind == "derivative_inconsistency"
    assert failed.value.diagnostics.reduced_hessian_spectrum_status == "nonfinite"
    assert (
        failed.value.diagnostics.first_failed_numerical_criterion
        == "derivative_evidence_finite"
    )


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        ("unknown_status", "invalid status"),
        ("nonfinite_attempt", "nonfinite value"),
        ("inconsistent_prefix", "budget prefix"),
    ),
)
def test_public_chemical_equilibrium_rejects_malformed_search_receipt(
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
    message: str,
) -> None:
    native_solve = epcsaft_equilibrium._equilibrium._chemical_equilibrium

    def malformed_search(*args: object) -> dict[str, object]:
        native = native_solve(*args)
        search = copy.deepcopy(native["search"])
        if mutation == "unknown_status":
            search["status"] = "unknown"
        elif mutation == "nonfinite_attempt":
            search["attempts"][0]["objective"] = math.nan
        else:
            search["budget_prefixes"][-1]["attempted_primary_ordinals"] = (0,)
        native["search"] = search
        return native

    monkeypatch.setattr(
        epcsaft_equilibrium._equilibrium,
        "_chemical_equilibrium",
        malformed_search,
    )
    with pytest.raises(
        epcsaft_equilibrium.ChemicalEquilibriumError,
        match=message,
    ):
        epcsaft_equilibrium.chemical_equilibrium(
            _ideal_phase(),
            350.0 * epcsaft.unit_registry.kelvin,
            200_000.0 * epcsaft.unit_registry.pascal,
            _ideal_problem(),
        )


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        ("spectrum_status", "spectrum status"),
        ("raw_inertia", "spectrum and inertia"),
        ("nullspace_shape", "inconsistent dimensions"),
    ),
)
def test_public_chemical_equilibrium_rejects_malformed_curvature_receipt(
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
    message: str,
) -> None:
    native_solve = epcsaft_equilibrium._equilibrium._chemical_equilibrium

    def malformed_curvature(*args: object) -> dict[str, object]:
        native = native_solve(*args)
        if mutation == "spectrum_status":
            native["reduced_hessian_spectrum_status"] = "unknown"
        elif mutation == "raw_inertia":
            native["reduced_hessian_raw_inertia"] = (0, 0, 0)
        else:
            native["reduced_hessian_nullspace_shape"] = (1, 1)
        return native

    monkeypatch.setattr(
        epcsaft_equilibrium._equilibrium,
        "_chemical_equilibrium",
        malformed_curvature,
    )
    with pytest.raises(epcsaft_equilibrium.ChemicalEquilibriumError, match=message):
        epcsaft_equilibrium.chemical_equilibrium(
            _ideal_phase(),
            350.0 * epcsaft.unit_registry.kelvin,
            200_000.0 * epcsaft.unit_registry.pascal,
            _ideal_problem(),
        )


def test_public_chemical_equilibrium_exports_typed_value_and_jacobian_contract() -> None:
    expected = {
        "ChemicalArtifactIdentity",
        "ChemicalEquilibriumActiveParameter",
        "ChemicalEquilibriumConstant",
        "ChemicalEquilibriumCriterion",
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
