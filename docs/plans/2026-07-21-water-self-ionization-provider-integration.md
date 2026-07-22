# Water Self-Ionization Provider Integration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` or `superpowers:executing-plans` to
> implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for
> tracking.

**Goal:** Prove one source-complete, nonideal, homogeneous water
self-ionization solve at `298.15 K` and `100000 Pa` through immutable installed
Provider and Equilibrium wheels.

**Architecture:** Provider adds one Held-2008 Debye--Huckel model bundle and an
append-only native reference/basis tail. Equilibrium transforms the IAPWS
mixed-standard-state `ln Kw` into that exact Helmholtz basis and reuses the
existing private reaction compiler, electroneutral chart, Ipopt TNLP, and
certificates.

**Tech Stack:** Python 3.13, C++20, pybind11, CMake/scikit-build-core, CppAD,
Ipopt, Pint, pytest, Ruff, mypy, immutable wheels.

## Global Constraints

- Work from the clean Equilibrium branch `codex/chemical-equilibrium` and a new
  isolated Provider branch `codex/water-self-ionization-reference`; do not edit
  either ordinary checkout when it contains user-owned IDE state.
- Provider is the only EOS, association, Debye--Huckel, density, packing,
  reference-sequence, and nonlinear derivative owner.
- Equilibrium is the only reaction, conservation, standard-state
  transformation, Ipopt, KKT, and certificate owner.
- Cross-package compilation and tests use an exact installed Provider wheel and
  its public native SDK; never compile Provider source into Equilibrium.
- The Held model is `ionic-held-2008-debye-huckel`; Born and SSM+DS are absent,
  not zeroed by a caller option and not borrowed from another source bundle.
- The only admitted state in this slice is `298.15 K` and `100000 Pa`.
- Keep all new Equilibrium routes private and underscored. Add no public solve,
  result family, selector, backend registry, compatibility shim, receipt,
  promotion, merge, or release.
- D-026 remains the public HELD2 admission gate. Predictive agreement remains
  `not_adjudicated`; globality remains `not_guaranteed`.
- Apply TDD to every runtime change. Production code uses exact CppAD or exact
  analytic linear transformations; directional finite differences exist only
  in tests.
- End every repository checkpoint with its focused tests, static checks,
  `git diff --check`, and `bash "$HOME/.codex/hooks/codex-cleanup.sh" --repo-root .`.

---

## Outcome Proof

**Intent:** Establish the first real reacting-liquid input without reviving the
retired chemical API or inventing a Provider reference state.

**Current Behavior:** Equilibrium solves analytic and installed-Provider
manufactured reactions, while Provider has no reactive basis metadata, no
Held-2008 hydronium/hydroxide model, and no gauge-invariant mixed-standard-state
reference tail.

**Expected Outcome:** An installed Provider artifact supplies the exact
Held-2008 phase tensors and neutral reference quantities; Equilibrium transforms
IAPWS `ln Kw`, solves the conserved charged system, and independently closes
the activity product and KKT certificates.

**Target Output:** Two stable working-branch commits and immutable wheel hashes,
with one private Equilibrium test result labeled
`installed_provider_source_complete_consistency`.

**Owner:** Provider owns model/reference evaluation; Equilibrium owns reaction
transformation/solve; Applications own later chemistry selection and datasets.

**Interface:** Append-only `epcsaft.native_sdk.v1` table fields plus
`ProviderContext::evaluate_standard_reference`, consumed only through an
installed Provider capsule.

**Cutover:** None; this appends a private candidate and leaves every accepted
public route unchanged.

**Replaced Path:** The installed-Provider manufactured reaction remains a
test-only diagnostic; it is not renamed or treated as the source-backed case.
The retired lab `eos_x_gamma` construction remains unused.

**Evidence:** Provider package tests, Equilibrium package tests, exact
directional derivative checks, IAPWS activity recomputation, wheel symbol/file
inspection, and isolated cross-wheel execution.

**Acceptance Proof:** The final true-`lambda=1` result has balance, exact charge,
pressure, packing/domain, reaction-affinity, KKT, complementarity, and reduced-
Hessian certificates within the existing private kernel limits, while the
independently recomputed mixed-standard activity product matches retained
`ln Kw` within `1e-8`.

**Stop Criteria:** Stop before solving on any source, formula, order, charge,
fingerprint, basis, standard-state, domain, reference-convergence, or wheel-
isolation mismatch; reject any indeterminate Ipopt or false-success path.

**Avoid:** Provider reaction ownership, copied EOS equations, finite-difference
production derivatives, source mixing, epsilon species, public claims, and MEA
or Lithium fixtures.

**Risk:** The Held-2008 model topology and mixed-standard sign convention are
the primary scientific risks; mutation tests and independent activity
recomputation own their detection.

## Implementation Boundaries

**Files To Create:** One Provider bundle directory, one Equilibrium IAPWS data
record, Provider and Equilibrium focused tests, and this slice's canonical
documentation updates.

**Files To Modify:** Provider record resolution/native model/reference/SDK
owners and Equilibrium Provider wrapper/chemical compiler/bindings/tests.

**Files To Avoid:** Provider equations outside its repository, public
Equilibrium Python exports, Regression, Migration, Validation, retired Windows
worktree runtime, and source-incomplete MEA/Lithium data.

**Source Of Truth:** Held et al. 2008 DOI `10.1016/j.fluid.2008.06.010`, IAPWS
R11-24, Provider immutable parameter records, and the approved design at
`docs/designs/2026-07-21-water-self-ionization-provider-integration.md`.

**Read Path:** Source records -> Provider `ParameterBundle`/`ParameterSet` ->
installed capsule -> reference tail and phase tensors -> Equilibrium
standard-state transformer -> existing compiler/TNLP -> certificates.

**Write Path:** TDD commits first to the isolated Provider branch, build and
hash its wheel, install that wheel into Equilibrium's isolated build/test
environment, then commit Equilibrium integration and evidence.

**Integration Points:** `epcsaft.native_sdk.v1`,
`epcsaft_equilibrium::ProviderContext`, `compile_reaction_system`, and
`solve_provider_manufactured_reaction`'s existing physical block.

**Migration Or Cutover:** None; Migration remains read-only and no authority
record changes.

**Replaced Path Handling:** Preserve the manufactured private seam with its
`manufactured_nonpredictive` label; add a distinct source-complete private seam
and prove no retired selector or `eos_x_gamma` code is reachable.

**Acceptance Proof Gate:** Fresh Provider and Equilibrium wheel builds, hashes,
full package tests, static checks, artifact-content inspection, and an isolated
environment containing only those wheels and runtime dependencies.

### Task 1: Add the source-bound Held-2008 Provider model

**Repository:** `ePC-SAFT/ePC-SAFT`

**Use Cases:**

- Exact water/hydronium/hydroxide records resolve in source order with formulas,
  charges, temperature domain, and fingerprint evidence.
- Missing formulation identity, ion parameter, water association record, or
  explicit ion--ion noninteraction fails before native construction.
- The new formulation appends one model; it does not redirect or duplicate any
  Figiel or Ascani path.

**Files:**

- Create: `src/epcsaft/data/bundles/held-2008-water-self-ionization/1/bundle.toml`
- Create: `src/epcsaft/data/bundles/held-2008-water-self-ionization/1/components.csv`
- Create: `src/epcsaft/data/bundles/held-2008-water-self-ionization/1/single.csv`
- Create: `src/epcsaft/data/bundles/held-2008-water-self-ionization/1/pair.csv`
- Create: `src/epcsaft/data/bundles/held-2008-water-self-ionization/1/sites.csv`
- Create: `src/epcsaft/data/bundles/held-2008-water-self-ionization/1/association.csv`
- Create: `src/epcsaft/data/bundles/held-2008-water-self-ionization/1/correlations.toml`
- Create: `src/epcsaft/data/bundles/held-2008-water-self-ionization/1/model.csv`
- Modify: `src/epcsaft/records.py`
- Modify: `src/epcsaft/_resolve.py`
- Modify: `src/epcsaft/_native/model.hpp`
- Modify: `src/epcsaft/_native/eos.cpp`
- Modify: `src/epcsaft/_native/helmholtz.hpp`
- Modify: `src/epcsaft/_native/electrolyte.hpp`
- Modify: `src/epcsaft/_native/born.hpp`
- Test: `tests/test_parameters.py`
- Test: `tests/test_eos.py`
- Test: `tests/test_debye_huckel.py`
- Test: `tests/test_native_sdk.py`

**Interfaces:**

- Produces: resolved topology literal `ionic-held-2008-debye-huckel` and native
  constant `ELECTROLYTE_MODEL_HELD_2008_DEBYE_HUCKEL`.
- Preserves: existing topology inference and fingerprints for every installed
  Figiel, Khudaida, Ascani, and Gross bundle.

- [ ] **Step 1: Write failing source and mutation tests.** Add a catalog test
  that selects exactly

  ```python
  ids = ("water", "hydronium-cation", "hydroxide-anion")
  parameters = ParameterBundle.from_catalog(
      "held-2008-water-self-ionization", version=1
  ).select(ids)
  assert parameters.component_ids == ids
  assert tuple(record.formula for record in parameters.components) == (
      "H2O", "H3O+", "OH-"
  )
  assert EPCSAFT(parameters)._model.topology == "ionic-held-2008-debye-huckel"
  ```

  Parameterize mutations removing the formulation record, `H+`/`OH-` Table 2
  values, a 2B water site/pair, or the explicit ion--ion `k_ij = 1` record and
  require `ParameterError` with the exact missing family.

- [ ] **Step 2: Run the focused tests and verify RED.** Run
  `python -m pytest tests/test_parameters.py tests/test_eos.py -q`. Expected:
  collection succeeds and the new catalog/topology tests fail because the
  bundle and topology do not exist.

- [ ] **Step 3: Add the immutable source records.** Encode Held Table 1 water
  values `m=1.204659`, `sigma=2.792700 A`, exponential terms
  `(10.1100,-0.01775)` and `(-1.41700,-0.01146)`,
  `epsilon/k=353.9449 K`, 2B association `2425.6714 K` and `0.0450989`;
  Table 2 hydrated-ion values `H+=2.2740 A,1616.4939 K` and
  `OH-=1.6401 A,2444.7555 K`; charges `(0,+1,-1)`; and exact source locators.
  Use the explicit true species names and formulas from the approved design.

- [ ] **Step 4: Add one typed formulation identity.** Extend the model-record
  parser with the closed value
  `held-2008-debye-huckel` and branch `_resolve.py` only when that immutable
  record is present. The branch requires neutral-solvent permittivity,
  association records, all pair records, opposite charges, no Born diameter,
  no solvation factor, no dielectric-suppression record, and no ionic-region
  permittivity.

- [ ] **Step 5: Implement the native formulation branch.** Add
  `ELECTROLYTE_MODEL_HELD_2008_DEBYE_HUCKEL = 3`; use the paper's constant
  solvent dielectric in `electrolyte_contribution`; retain the existing ion
  diameter rule; and return exact zero from `born_contribution` only for this
  immutable model. Keep the existing Figiel and Ascani branches byte-for-byte
  behaviorally unchanged.

- [ ] **Step 6: Prove contribution ownership and exact derivatives.** Add tests
  asserting finite hard-chain, dispersion, association, and Debye--Huckel
  contributions, exact zero Born, nonzero charge response, and directional
  agreement of the native phase gradient/Hessian at a positive electroneutral
  state. Run
  `python -m pytest tests/test_debye_huckel.py tests/test_native_sdk.py -q`;
  expected: PASS.

- [ ] **Step 7: Run Provider package checks and commit.** Run the complete
  Provider pytest suite, Ruff, strict mypy, `git diff --check`, and cleanup.
  Commit only the listed files as
  `feat: add Held water-ionization provider model`.

### Task 2: Append machine-readable Helmholtz and pressure-domain metadata

**Repository:** `ePC-SAFT/ePC-SAFT`

**Use Cases:**

- Native consumers can fail closed unless the exact mechanical basis is
  `A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3`.
- Source pressure bounds propagate from immutable records; undeclared bounds
  are explicit `NaN`, not hidden infinities.
- Old SDK prefix callbacks retain offsets, sizes, and behavior.

**Files:**

- Modify: `src/epcsaft/_resolve.py`
- Modify: `src/epcsaft/eos.py`
- Modify: `src/epcsaft/_native/module.cpp`
- Modify: `src/epcsaft/_native/native_sdk.cpp`
- Modify: `src/epcsaft/include/epcsaft/native_sdk_v1.h`
- Test: `tests/test_parameters.py`
- Test: `tests/test_native_sdk.py`

**Interfaces:**

- Produces appended table fields
  `helmholtz_basis_id`, `reference_amount_mol`,
  `reference_number_density_mol_per_m3`, `source_pressure_min_pa`, and
  `source_pressure_max_pa`.

- [ ] **Step 1: Write failing ABI and domain tests.** Extend the ctypes table
  definition only after the existing active-Born tail and assert the stable
  basis string, both reference scalars equal `1.0`, and exact/`NaN` pressure
  semantics for Held, Khudaida, and Gross models. Add a mutation with reversed
  pressure bounds and require catalog rejection.

- [ ] **Step 2: Run RED.** Run
  `python -m pytest tests/test_parameters.py tests/test_native_sdk.py -q`.
  Expected: new fields are absent and the ctypes size assertion fails.

- [ ] **Step 3: Resolve pressure bounds without a second domain owner.** Change
  `_domain_limits` to return `(T_min, T_max, P_min, P_max, ion_x_max)` from the
  selected records' existing `ValidityDomain` objects. Thread those two values
  through `_ResolvedModel`, `EPCSAFT.native_sdk`, pybind, and
  `NativeSdkEndpoint`; convert only missing values to quiet `NaN`.

- [ ] **Step 4: Append the native fields.** Define
  `EPCSAFT_NATIVE_SDK_V1_HELMHOLTZ_BASIS_ID` in the public header and append the
  five fields after the current tail. Initialize them from constants and the
  resolved model; never insert them into the fixed prefix.

- [ ] **Step 5: Run tests and commit.** Run the focused tests, full Provider
  suite, Ruff, mypy, diff check, and cleanup. Commit as
  `feat: expose provider Helmholtz reference metadata`.

### Task 3: Add the gauge-invariant electrolyte standard-reference callback

**Repository:** `ePC-SAFT/ePC-SAFT`

**Use Cases:**

- A fixed-`T,P` one-solvent/one-cation/one-anion model returns a converged
  formula-unit infinite-dilution log-fugacity sum and pure-solvent reference.
- Permuted, multi-solvent, multiple-cation/anion, non-1:1, wrong-temperature,
  wrong-pressure, or nonconverged paths fail closed.
- The public reference-state path and native callback share one sequence owner;
  no reference iteration is copied.

**Files:**

- Modify: `src/epcsaft/_native/electrolyte_reference.hpp`
- Modify: `src/epcsaft/_native/eos.hpp`
- Modify: `src/epcsaft/_native/eos.cpp`
- Modify: `src/epcsaft/_native/module.cpp`
- Modify: `src/epcsaft/_native/native_sdk.cpp`
- Modify: `src/epcsaft/include/epcsaft/native_sdk_v1.h`
- Modify: `src/epcsaft/reference.py`
- Test: `tests/test_reference.py`
- Test: `tests/test_native_sdk.py`

**Interfaces:**

- Produces `epcsaft_electrolyte_standard_reference_result_v1` with
  `formula_unit_infinite_dilution_log_fugacity`,
  `pure_solvent_log_fugacity_coefficient`, `solvent_molar_mass_kg_per_mol`,
  `reference_molality_mol_per_kg`, `reference_convergence_error`, state,
  fingerprint, basis ID, and error/status.
- Produces appended callback `evaluate_electrolyte_standard_reference`.

- [ ] **Step 1: Write failing reference tests.** For the Held model, call both
  the public reference state and the proposed SDK callback at `298.15 K` and
  `100000 Pa`; require the native formula-unit sum to equal the sum of public
  ion infinite-dilution log fugacities within `2e-12`, the pure-water value to
  match a direct public pure-water evaluation within `2e-12`, and convergence
  error at most `5e-5`.

- [ ] **Step 2: Write failure tests.** Assert `ABI_MISMATCH` for the wrong result
  size, `UNSUPPORTED_MODEL` for neutral/multi-ion models, `DOMAIN_ERROR` outside
  source `T`, and no partial success payload after any failure.

- [ ] **Step 3: Run RED.** Run
  `python -m pytest tests/test_reference.py tests/test_native_sdk.py -q`.
  Expected: the result struct and callback do not exist.

- [ ] **Step 4: Extend the sole reference-sequence result.** Add the pure-
  solvent log fugacity coefficient to the native sequence result already used
  by `reference.py`, calculated through the existing pressure-state evaluator.
  Preserve the dilution schedule `(1e-6,1e-8,1e-10,1e-12) mol kg^-1` and
  terminal `5e-5` convergence rule.

- [ ] **Step 5: Implement the append-only callback.** Validate exact component
  topology/order/charges, source domain, positive state, exact result size, and
  finite returned scalars. Return only the cation+anion sum, never label either
  individual ion as a predictive standard chemical potential.

- [ ] **Step 6: Run Provider verification and commit.** Run focused/full tests,
  Ruff, mypy, diff check, cleanup, and commit as
  `feat: expose electrolyte standard-reference tail`.

### Task 4: Build and freeze the Provider checkpoint artifact

**Repository:** `ePC-SAFT/ePC-SAFT`

**Use Cases:**

- Equilibrium receives one immutable Provider wheel and public header from the
  exact stable Provider commit.
- The artifact contains no Equilibrium, Regression, Migration, Validation, or
  lab source and does not displace an earlier wheel.

**Files:**

- Modify: `README.md`
- Modify: `CONTEXT.md`
- Modify: `ARCHITECTURE.yaml`
- Modify: `docs/science/full-provider-eos.md`

**Interfaces:**

- Produces an immutable wheel path and SHA-256 consumed by Tasks 5--8.

- [ ] **Step 1: Document exact candidate scope.** Record the Held topology,
  source, basis tail, reference callback, pressure metadata, all exclusions,
  and `authority_effect: none`. Keep reaction equilibrium explicitly
  Equilibrium-owned.

- [ ] **Step 2: Run architecture/schema checks.** Parse `ARCHITECTURE.yaml`, run
  the full Provider tests and static checks, inspect `git diff --check`, and run
  cleanup. Commit as `docs: record water-ionization provider candidate`.

- [ ] **Step 3: Build without overwriting artifacts.** Build one wheel from the
  exact commit using the repository's documented isolated build command. Place
  it under a new commit-keyed project artifact directory, compute SHA-256, and
  inspect the wheel file list plus exported native header.

- [ ] **Step 4: Push the stable Provider branch.** Push
  `codex/water-self-ionization-reference`, then record exact commit, tree,
  wheel path, wheel hash, SDK table size, and full test count for the
  Equilibrium checkpoint.

### Task 5: Consume and validate the Provider reference tail in Equilibrium

**Repository:** `ePC-SAFT/ePC-SAFT-equilibrium`

**Use Cases:**

- Equilibrium rejects capsules lacking the appended basis/reference contract,
  wrong basis identity, wrong state/domain, wrong order/charges, or mismatched
  fingerprint.
- The wrapper returns the exact neutral reference quantities without rebuilding
  a Provider density or fugacity equation.
- Existing Provider consumers retain their behavior; the source-backed seam is
  additive and private.

**Files:**

- Modify: `cpp/src/provider.hpp`
- Modify: `cpp/src/provider.cpp`
- Modify: `cpp/src/chemical_equilibrium_bindings.cpp`
- Test: `tests/test_chemical_equilibrium.py`

**Interfaces:**

- Produces:

  ```cpp
  struct StandardReferenceEvaluation {
      double formula_unit_log_fugacity;
      double pure_solvent_log_fugacity;
      double solvent_molar_mass_kg_per_mol;
      double reference_molality_mol_per_kg;
      double convergence_error;
      std::string basis_id;
      std::string parameter_fingerprint;
  };

  StandardReferenceEvaluation ProviderContext::evaluate_standard_reference(
      double temperature_k,
      double pressure_pa
  ) const;
  ```

- [ ] **Step 1: Install the exact Provider wheel in the isolated Equilibrium
  environment.** Verify `epcsaft.__file__` resolves inside that environment and
  record the installed distribution version and artifact hash.

- [ ] **Step 2: Write failing capsule tests.** Add one success expectation and
  mutations for old table size, null callback, wrong result size, wrong basis,
  wrong fingerprint, and out-of-domain state. Expected failure messages must
  distinguish ABI, identity, source-domain, and Provider-evaluation failures.

- [ ] **Step 3: Run RED.** Run the focused chemical test. Expected: the wrapper
  method and private binding are absent.

- [ ] **Step 4: Implement the minimal wrapper.** Require the table through the
  new callback tail, exact basis ID and reference scalars, call the Provider
  function once, check return/result statuses, and validate every finite scalar
  plus fingerprint. Add only an underscored test binding.

- [ ] **Step 5: Run focused/static checks and commit.** Run the chemical tests,
  Ruff, strict mypy, diff check, cleanup, and commit as
  `feat: consume provider standard-reference tail`.

### Task 6: Transform the IAPWS mixed standard state exactly

**Repository:** `ePC-SAFT/ePC-SAFT-equilibrium`

**Use Cases:**

- The retained IAPWS input is dimensionless, source-identified, fixed to
  `298.15 K`/`100000 Pa`, and independent of Provider-predicted water density.
- The transformation sign, molar-mass units, and basis identity are protected
  by independent algebra and mutation tests.
- No generic standard-state registry or temperature correlation is added.

**Files:**

- Create: `data/reference/water_self_ionization_iapws_r11_24.yaml`
- Modify: `cpp/src/chemical_equilibrium.hpp`
- Modify: `cpp/src/chemical_equilibrium.cpp`
- Modify: `cpp/src/chemical_equilibrium_bindings.cpp`
- Test: `tests/test_chemical_equilibrium.py`

**Interfaces:**

- Produces:

  ```cpp
  struct MixedStandardStateResult {
      double delta_standard_offset;
      double ln_k_provider_basis;
  };

  MixedStandardStateResult transform_water_self_ionization_standard_state(
      double ln_kw_mixed_standard,
      const StandardReferenceEvaluation& reference
  );
  ```

- [ ] **Step 1: Retain the exact IAPWS record.** Store IAPWS R11-24 parameters,
  `T=298.15 K`, `P=100000 Pa`, independently sourced ordinary-water density
  `997.047039 kg m^-3`, `pKw=13.99435488204`, and
  `lnKw=-32.2231929374538`, with the molal-ion/mole-fraction-water convention
  and equation locator. Record the density source separately from Provider.

- [ ] **Step 2: Write RED algebra tests.** Given a manufactured Provider
  reference, assert

  ```python
  delta = (
      2.0 * log(
          reference["solvent_molar_mass_kg_per_mol"]
          * reference["reference_molality_mol_per_kg"]
      )
      + reference["formula_unit_log_fugacity"]
      - 2.0 * reference["pure_solvent_log_fugacity"]
  )
  assert transformed["ln_k_provider_basis"] == approx(ln_kw + delta)
  ```

  Mutate `kg/mol` to `g/mol`, reverse the sign, change basis ID, state, source
  identity, or pressure binding and require pre-solve rejection.

- [ ] **Step 3: Run RED.** Run the focused tests. Expected: transformation and
  IAPWS loader bindings are absent.

- [ ] **Step 4: Implement the exact linear transformation.** Use only `log`,
  additions, and scalar validation. Require the Provider basis ID and a 1:1
  reference result. Extend `ReactionSystemInput`/compiled evidence with basis
  identity and transformation residual; do not change the existing SVD
  construction of `g_ref`.

- [ ] **Step 5: Run tests and commit.** Run focused tests, YAML parsing, Ruff,
  mypy, diff check, cleanup, and commit as
  `feat: transform IAPWS water-ionization standard state`.

### Task 7: Solve and certify the source-complete water reaction

**Repository:** `ePC-SAFT/ePC-SAFT-equilibrium`

**Use Cases:**

- Pure-water feed with zero ionic feed compiles under `B=[[2,3,1],[1,1,1]]`,
  `nu=[[-2,1,1]]`, and exact charge `(0,+1,-1)`.
- The existing TNLP reaches a final true `lambda=1` point and reports every
  certificate axis independently.
- Trace-floor, indeterminate solver, callback, domain, reference, pressure,
  packing, affinity, KKT, complementarity, or reduced-Hessian failures cannot
  return success.

**Files:**

- Modify: `cpp/src/chemical_equilibrium.hpp`
- Modify: `cpp/src/chemical_equilibrium_solver.cpp`
- Modify: `cpp/src/chemical_equilibrium_bindings.cpp`
- Test: `tests/test_chemical_equilibrium.py`

**Interfaces:**

- Produces private binding
  `_chemical_solve_provider_water_self_ionization(capsule, spec, options)` and
  profile label `installed_provider_source_complete_consistency`.
- Reuses `solve_provider_manufactured_reaction`'s phase evaluator/TNLP; it does
  not create a second `NlpProblem` owner.

- [ ] **Step 1: Write the failing end-to-end test.** Build the Held Provider
  model from the installed catalog, obtain the capsule/reference, transform the
  retained `lnKw`, compile the exact reaction system, solve, and assert accepted
  final `lambda=1`, exact charge, all existing strict certificate limits,
  `predictive_status == "not_adjudicated"`, and
  `globality_certificate == "not_guaranteed"`.

- [ ] **Step 2: Write independent activity recomputation.** From the final
  Provider fixed-pressure fugacity coefficients and amounts, compute molality
  and water mole-fraction activities outside the production transformation;
  require
  `abs(log(a_hydronium) + log(a_hydroxide) - 2*log(a_water) - lnKw) <= 1e-8`.

- [ ] **Step 3: Add invariance and negative tests.** Cover gauge, species
  permutation, feed amount scaling, reaction-basis scaling, trace floors from
  `1e-18` through `1e-12 mol`, stale fingerprint/basis, missing IAPWS field,
  nonneutral feed, reference nonconvergence, Ipopt indeterminate status, and a
  manufactured false-success callback.

- [ ] **Step 4: Run RED.** Run the focused end-to-end tests. Expected: the
  source-complete binding/profile does not exist.

- [ ] **Step 5: Reuse the existing solver owner.** Factor only the shared
  Provider TNLP call needed to accept precompiled `g_ref`; retain max-min,
  direct-first/adaptive continuation, final recomputation, KKT, and local-
  minimum code paths unchanged. Set the new profile label only after all input
  and physical certificates pass.

- [ ] **Step 6: Run tests and commit.** Run all chemical tests plus the full
  Equilibrium suite, Ruff, strict mypy, diff check, cleanup, and commit as
  `feat: solve source-complete water self-ionization`.

### Task 8: Freeze cross-wheel proof and private status

**Repository:** `ePC-SAFT/ePC-SAFT-equilibrium`

**Use Cases:**

- A fresh environment with only immutable wheels reproduces the result and
  proves no sibling Provider source or lab runtime was imported.
- Documentation reports exact artifact/test evidence without promotion or a
  predictive/global claim.
- The manufactured seam remains visibly nonpredictive; no duplicate public
  route or displaced-path ambiguity remains.

**Files:**

- Modify: `README.md`
- Modify: `CONTEXT.md`
- Modify: `ARCHITECTURE.yaml`
- Modify: `docs/phase-equilibrium.md`
- Modify: `docs/designs/2026-07-21-water-self-ionization-provider-integration.md`

**Interfaces:**

- Produces final Provider/Equilibrium commit, tree, wheel SHA-256, test counts,
  SDK basis ID, and unresolved-science evidence for handoff.

- [ ] **Step 1: Update private candidate documentation.** Record the exact
  source-complete input status, consistency-only result, artifact IDs, explicit
  exclusions, next MEA reference/species blockers, and `authority_effect: none`.

- [ ] **Step 2: Build the Equilibrium wheel.** Build from the stable commit in a
  clean environment containing the frozen Provider wheel, hash it, inspect its
  file/symbol inventory, and prove it contains no Provider implementation
  source or copied EOS symbol.

- [ ] **Step 3: Run isolated installed-artifact proof.** Create a fresh virtual
  environment, install only the two wheels and declared runtime dependencies,
  print both module paths and distribution versions, run the source-complete
  test, and scan imported module paths for lab/sibling-source roots. Expected:
  one accepted consistency result and zero forbidden paths.

- [ ] **Step 4: Run final verification.** Run the full Equilibrium suite, Ruff,
  strict mypy, YAML parsing, staged diff review, `git diff --check`, and cleanup.
  Retain exact command outputs and test counts.

- [ ] **Step 5: Commit and push the stable Equilibrium checkpoint.** Commit as
  `docs: record water-ionization checkpoint evidence`, push
  `codex/chemical-equilibrium`, and report exact commits, trees, artifacts,
  hashes, tests, claim axes, and remaining MEA blockers. Do not merge, promote,
  edit Migration, or write Validation.
