# GREPE Homogeneous Chemical Layer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` (recommended) or
> `superpowers:executing-plans` to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Deepen the existing private fixed-\(T,P\) chemical-equilibrium module
with GREPE reaction compilation, validated homogeneous structural support,
accessible-face recompilation, and explicit local or unresolved certification.

**Architecture:** The private compiler accepts unreduced reaction rows and
source records, selects and records one independent basis, validates reaction
cycles before reduction, and constructs the complete invariant/charge
decomposition. A new private support implementation uses HiGHS only to propose
LP solutions and exact rational arithmetic to validate deletion or
accessibility certificates. The existing amount chart and
`solve_provider_reaction` remain the only optimization path.

**Tech Stack:** C++17, pybind11, HiGHS 1.15.1, Ipopt, Boost header-only exact
integer/rational arithmetic, Python 3.13, pytest.

## Global Constraints

- Keep the chemical capability private and underscored; add no public Python
  export or result family.
- Consume only the exact non-editable installed Provider through
  `epcsaft/native_sdk_v1.h`; do not compile Provider source or copy EOS logic.
- Preserve `solve_provider_reaction` as the sole Provider-backed homogeneous
  chemical solver.
- Convert source reaction constants before cycle checking; this slice accepts
  only records already converted to the declared Provider convention and fails
  closed on other conventions.
- Delete a species only after an independently validated exact zero
  certificate. Unresolved support remains retained.
- Do not represent exact structural zeros with artificial epsilon amounts.
- Emit no phase-pricing, `SEARCH_STABLE`, `CERTIFIED_EPS_GLOBAL`, coupled
  phase-equilibrium, sensitivity, release, promotion, or authority claim.
- Keep persistent scientific evidence to the three compact parameterized
  families approved in the design.

---

## File Structure

- Modify `cpp/src/chemical_equilibrium.hpp`: private compiler, support, and
  result evidence types.
- Modify `cpp/src/chemical_equilibrium.cpp`: reaction/invariant reduction,
  cycle validation, accessible-face compilation, and existing amount chart.
- Create `cpp/src/chemical_equilibrium_support.cpp`: homogeneous HiGHS LP
  candidates and exact rational certificate validation.
- Modify `cpp/src/chemical_equilibrium_solver.cpp`: reduced-face solve,
  full-species result expansion, and GREPE-compatible local evidence.
- Modify `cpp/src/chemical_equilibrium_bindings.cpp`: parse and serialize the
  deeper private contract.
- Modify `CMakeLists.txt`: compile the support implementation in the existing
  native module.
- Modify `tests/test_chemical_equilibrium.py`: compact compiler, support, and
  end-to-end evidence.
- Modify `CONTEXT.md`: record the implemented private GREPE homogeneous layer
  and its fail-closed limits.

### Task 1: Compile Redundant Reaction Sets Before Solving

**Files:**

- Modify: `tests/test_chemical_equilibrium.py`
- Modify: `cpp/src/chemical_equilibrium.hpp`
- Modify: `cpp/src/chemical_equilibrium.cpp`
- Modify: `cpp/src/chemical_equilibrium_bindings.cpp`

**Interfaces:**

- Consumes: private Python specification with `species_ids`, `charges`,
  `molar_masses_kg_per_mol`, `balance_matrix`, supplied `reaction_matrix`,
  supplied `ln_k`, source records, feed, \(T\), \(P\), and Provider
  fingerprint.
- Produces: `CompiledReactionSystem` with supplied and independent matrices,
  `reaction_transform`, selected basis rows, converted constants, mass-first
  independent balances, and cycle/reconstruction residuals.

- [ ] **Step 1: Extend the test specification with explicit molar masses and
  converted-record provenance**

  Update `_base_system()` and `_bind_record()` so every record contains:

  ```python
  {
      "source_id": "manufactured:reaction-0",
      "reference_id": "provider-helmholtz-coordinate-basis",
      "reaction_orientation": "products_positive",
      "conversion_id": "already-provider-basis",
      "dimensionless": True,
      "temperature_k": temperature_k,
      "pressure_pa": pressure_pa,
  }
  ```

  Add `molar_masses_kg_per_mol` to every changed species topology. Select
  positive values satisfying every supplied reaction exactly within the
  compiler's scaled floating-point audit.

- [ ] **Step 2: Write the failing redundant-reaction and cycle tests**

  Add one parameterized compiler test using:

  ```python
  reaction_matrix = (
      (-1.0, 1.0, 0.0),
      (0.0, -1.0, 1.0),
      (-1.0, 0.0, 1.0),
  )
  ln_k = (math.log(2.0), math.log(3.0), math.log(6.0))
  ```

  Require rank two, a three-by-two transform, reconstruction of all three
  supplied rows/constants, and a small cycle residual. Parameterize one
  inconsistent third constant and require rejection containing
  `"reaction constant cycle"`.

- [ ] **Step 3: Run the focused test and verify RED**

  Build the editable native module from this worktree, then run:

  ```bash
  .venv/bin/python -m pytest -q \
    tests/test_chemical_equilibrium.py \
    -k "redundant_reaction or reaction_constant_cycle"
  ```

  Expected: failure because the compiler still rejects dependent reaction
  rows and the new record fields are not represented.

- [ ] **Step 4: Add the unreduced compiler evidence types**

  Extend `EquilibriumConstantRecord` with:

  ```cpp
  std::string reaction_orientation;
  std::string conversion_id;
  ```

  Extend `ReactionSystemInput` with:

  ```cpp
  std::vector<double> molar_masses_kg_per_mol;
  ```

  Extend `CompiledReactionSystem` with:

  ```cpp
  std::vector<double> molar_masses_kg_per_mol;
  DenseMatrix supplied_balance_matrix;
  DenseMatrix supplied_reaction_matrix;
  std::vector<double> supplied_ln_k;
  DenseMatrix reaction_transform;
  std::vector<std::size_t> reaction_basis_rows;
  double reaction_cycle_inf_norm = 0.0;
  double reaction_transform_inf_norm = 0.0;
  ```

- [ ] **Step 5: Implement basis selection and pre-reduction cycle validation**

  In `chemical_equilibrium.cpp`:

  1. validate positive finite molar masses and mass conservation on every
     supplied reaction;
  2. require every record to be dimensionless,
     `products_positive`, and already bound to the Provider basis;
  3. select numerically independent supplied rows with scaled pivoting;
  4. solve each supplied row as a combination of retained rows to form
     `reaction_transform`;
  5. use retained source constants as the independent constants;
  6. compute `reaction_transform * independent_ln_k - supplied_ln_k`; and
  7. reject an inconsistent cycle before constructing `g_ref`.

  Rank-reduce the supplied balance rows while treating charge as a separately
  seeded row and forcing molar mass to be the first retained balance row.
  Require the final invariant/charge rank plus reaction rank to equal the
  species count.

- [ ] **Step 6: Serialize compiler evidence through the private binding**

  Parse the new input fields and return:

  ```python
  {
      "supplied_reaction_rank": ...,
      "reaction_basis_rows": ...,
      "reaction_transform": ...,
      "independent_reaction_matrix": ...,
      "independent_ln_k": ...,
      "reaction_cycle_inf_norm": ...,
      "reaction_transform_inf_norm": ...,
      "balance_matrix": ...,
      "molar_masses_kg_per_mol": ...,
  }
  ```

  Keep `_chemical_compile_system` underscored.

- [ ] **Step 7: Rebuild and verify GREEN**

  Run the same focused pytest selection. Expected: all selected tests pass and
  inconsistent converted cycles fail before reference reconstruction.

- [ ] **Step 8: Review the Task 1 diff**

  Run:

  ```bash
  git diff --check
  git diff -- cpp/src/chemical_equilibrium.hpp \
    cpp/src/chemical_equilibrium.cpp \
    cpp/src/chemical_equilibrium_bindings.cpp \
    tests/test_chemical_equilibrium.py
  ```

  Confirm that the old already-independent path has been replaced rather than
  retained as a second compiler.

### Task 2: Classify Homogeneous Structural Support With Exact Certificates

**Files:**

- Modify: `tests/test_chemical_equilibrium.py`
- Modify: `cpp/src/chemical_equilibrium.hpp`
- Create: `cpp/src/chemical_equilibrium_support.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**

- Consumes:

  ```cpp
  HomogeneousSupportAnalysis analyze_homogeneous_support(
      const DenseMatrix& balance_matrix,
      const std::vector<double>& balance_totals,
      const std::vector<int>& charges,
      const std::vector<double>& molar_masses_kg_per_mol
  );
  ```

- Produces one `SpeciesSupportEvidence` per original species plus Phase-I
  status and a witness-average state.

- [ ] **Step 1: Write the failing structural-support tests**

  Add one parameterized family containing:

  - `A <=> B` with an independently conserved absent `C`, requiring `C` to be
    `proved_structural_zero`;
  - neutral `A` and a jointly accessible `C+`/`D-` pair, requiring both ions
    to have validated positive witnesses; and
  - a scaled equivalent of one case, requiring the same exact support
    classification and witness after normalization.

  Assert exact balance and charge residuals for every validated witness and for
  the witness-average state.

- [ ] **Step 2: Run the support tests and verify RED**

  ```bash
  .venv/bin/python -m pytest -q \
    tests/test_chemical_equilibrium.py -k "structural_support"
  ```

  Expected: failure because no support evidence exists.

- [ ] **Step 3: Add private support evidence types**

  Add:

  ```cpp
  struct SpeciesSupportEvidence {
      std::string classification;
      double candidate_maximum_mass_fraction = 0.0;
      bool primal_validated = false;
      bool dual_validated = false;
      std::vector<double> witness_amounts;
      std::vector<double> dual_multipliers;
  };

  struct HomogeneousSupportAnalysis {
      std::string phase1_status;
      std::string validation_status;
      std::vector<SpeciesSupportEvidence> species;
      std::vector<double> witness_average_amounts;
      double equality_inf_norm = 0.0;
  };
  ```

- [ ] **Step 4: Implement HiGHS candidate LPs**

  In `chemical_equilibrium_support.cpp`, build the independent equality system
  from balances and charge. Run:

  1. a zero-objective feasibility LP; and
  2. one bounded maximum-mass-fraction LP per species.

  Use one deterministic HiGHS thread and disable console output. Record model
  status, primal candidate, row duals, column duals, and objective.

- [ ] **Step 5: Implement exact rational validation**

  Use:

  ```cpp
  using ExactInteger = boost::multiprecision::cpp_int;
  using ExactRational = boost::rational<ExactInteger>;
  ```

  Convert every finite `double` to its exact binary rational value. Reconstruct
  candidate basic primal and dual solutions by exact Gaussian elimination on
  independent active columns. Validate:

  - exact equality feasibility;
  - exact nonnegativity;
  - a strictly positive target objective for accessibility;
  - exact dual inequalities for the reachability upper bound; and
  - exact zero dual objective before declaring structural zero.

  If either certificate cannot be reconstructed or validated, classify the
  species `unresolved`.

- [ ] **Step 6: Construct and validate the witness average**

  Average all exact accessible witnesses in floating output coordinates.
  Independently recompute balance and charge residuals. Require every
  `proved_accessible` coordinate to be positive. Keep unresolved coordinates
  retained but do not claim the witness average proves them accessible.

- [ ] **Step 7: Compile and bind the support evidence**

  Add `chemical_equilibrium_support.cpp` to the existing `_equilibrium` target.
  Add one underscored binding:

  ```python
  _chemical_analyze_homogeneous_support(
      balance_matrix,
      feed_amounts,
      charges,
      molar_masses_kg_per_mol,
  )
  ```

  The binding exists only for the package's compact white-box tests.

- [ ] **Step 8: Rebuild and verify GREEN**

  Run the structural-support selection and require all selected tests to pass.

### Task 3: Recompile the Accessible Chemical Face

**Files:**

- Modify: `tests/test_chemical_equilibrium.py`
- Modify: `cpp/src/chemical_equilibrium.hpp`
- Modify: `cpp/src/chemical_equilibrium.cpp`
- Modify: `cpp/src/chemical_equilibrium_bindings.cpp`

**Interfaces:**

- Consumes: independent full-species reaction/invariant system and validated
  support classifications.
- Produces: accessible species map, accessible reaction transform, reduced
  invariant/reaction system, and exact-zero restoration map.

- [ ] **Step 1: Write the failing cancellation test**

  Use species `(A, B, C, X)` with masses `(1, 2, 1, 1)`, feed `A=1`, and:

  ```python
  reaction_matrix = (
      (-1.0, 1.0, 0.0, -1.0),  # A + X <=> B
      (0.0, -1.0, 1.0, 1.0),   # B <=> C + X
  )
  ```

  Require `B` and `X` to be proved structural zeros while the accessible
  combination `A <=> C` survives with the summed equilibrium constant.
  Assert the null-space identity on the reduced system.

- [ ] **Step 2: Run the cancellation test and verify RED**

  ```bash
  .venv/bin/python -m pytest -q \
    tests/test_chemical_equilibrium.py -k "accessible_face"
  ```

  Expected: failure because the current compiler does not analyze or reduce
  support.

- [ ] **Step 3: Add accessible mapping evidence**

  Extend `CompiledReactionSystem` with:

  ```cpp
  HomogeneousSupportAnalysis support;
  std::vector<std::size_t> retained_species_indices;
  std::vector<std::size_t> removed_species_indices;
  DenseMatrix accessible_reaction_transform;
  std::vector<std::string> original_species_ids;
  std::vector<int> original_charges;
  std::vector<double> original_feed_amounts;
  ```

- [ ] **Step 4: Implement accessible reaction combinations**

  For removed columns \(N_{\mathrm{rem}}\), construct a basis of
  `null(transpose(N_rem))`. Form:

  ```text
  accessible_reactions = transpose(Z_acc) * reactions * P_acc
  accessible_ln_k = transpose(Z_acc) * ln_k
  ```

  Then restrict and rank-reduce invariant/charge rows, preserving mass, and
  verify:

  ```text
  rank([B_acc; z_acc]) + rank(nu_acc) == retained_species_count
  [B_acc; z_acc] * transpose(nu_acc) == 0
  ```

- [ ] **Step 5: Reconstruct accessible references and evidence**

  Build `g_ref` only after accessible-face recompilation. Record transform and
  reconstruction residuals. Reject a detected rank loss or inconsistent
  accessible constants; do not discard another reaction direction.

- [ ] **Step 6: Serialize accessible evidence**

  Return original and retained species IDs, maps, support classifications,
  accessible reaction matrix/constants, and accessible transform from
  `_chemical_compile_system`.

- [ ] **Step 7: Rebuild and verify GREEN**

  Run the cancellation test and the Task 1 compiler tests together.

### Task 4: Integrate Support and GREPE Levels Into the Existing Solver

**Files:**

- Modify: `tests/test_chemical_equilibrium.py`
- Modify: `cpp/src/chemical_equilibrium.hpp`
- Modify: `cpp/src/chemical_equilibrium_solver.cpp`
- Modify: `cpp/src/chemical_equilibrium_bindings.cpp`

**Interfaces:**

- Consumes: accessible `CompiledReactionSystem`.
- Produces: one full-species `ChemicalSolveResult` with
  `chemical_certification_level`, `boundary_status`, `support_qualifiers`, and
  restored exact structural zeros where supported.

- [ ] **Step 1: Write the failing result-level tests**

  Extend existing accepted manufactured and Belov assertions:

  ```python
  assert result["chemical_certification_level"] == "LOCAL_EQUILIBRIUM"
  assert result["boundary_status"] == "strict_interior"
  assert result["support_qualifiers"] == []
  ```

  Add the manufactured `(A, B, C, X)` structural-face solve and require:

  ```python
  assert result["amounts"][1] == 0.0
  assert result["amounts"][3] == 0.0
  assert result["boundary_status"] == "structural_face"
  assert result["chemical_certification_level"] == "LOCAL_EQUILIBRIUM"
  ```

  Add a focused Provider guard case using an exact installed component topology
  when one is available. If the installed catalog has no physically complete
  mass-conserving structural-face topology, exercise the guard through a
  private native helper over a compiled reduced system. Require fail-closed
  `BOUNDARY_DIRECTION_UNRESOLVED` before any Provider callback and without
  epsilon amounts.

- [ ] **Step 2: Run the result-level tests and verify RED**

  Run:

  ```bash
  .venv/bin/python -m pytest -q \
    tests/test_chemical_equilibrium.py \
    -k "certification_level or structural_face or boundary_direction"
  ```

  Expected: failure because the result has no GREPE-compatible evidence and
  the solver still assumes every original species is positive.

- [ ] **Step 3: Extend the private solve result**

  Add:

  ```cpp
  std::string chemical_certification_level = "FEASIBLE_ONLY";
  std::string boundary_status = "not_adjudicated";
  std::vector<std::string> support_qualifiers;
  std::vector<std::size_t> retained_species_indices;
  std::vector<std::size_t> structural_zero_species_indices;
  ```

- [ ] **Step 4: Solve the accessible manufactured ideal face**

  Run the existing amount chart, max-min initializer, Ipopt TNLP, KKT polish,
  and postsolve checks on the accessible system. Expand amounts to original
  order after acceptance and insert exact zeros only at validated structural
  zeros.

  Recompute original-space balances, charge, and named supplied-reaction
  affinities. A named reaction containing a removed intermediate is diagnostic
  only; require equilibrium for every retained accessible reaction
  combination.

- [ ] **Step 5: Fail Provider structural faces closed**

  Before Provider evaluation, check whether the compiled system removed any
  original component. Because native SDK v1 lacks a reduced-component topology
  contract, return:

  ```text
  chemical_certification_level = BOUNDARY_DIRECTION_UNRESOLVED
  boundary_status = boundary_direction_unresolved
  accepted = false
  ```

  Do not call Provider and do not create epsilon amounts.

- [ ] **Step 6: Assign primary levels from independent evidence**

  Assign `LOCAL_EQUILIBRIUM` only when current solver, numerical, physical,
  Provider-domain, local-minimum, trace, and support gates pass. Assign
  `INFEASIBLE_DECLARED_MODEL` only when exact support validation proves
  infeasibility. Otherwise retain `FEASIBLE_ONLY` or
  `BOUNDARY_DIRECTION_UNRESOLVED`.

  Never assign search-stable or global levels in this module.

- [ ] **Step 7: Serialize the result evidence**

  Add the new fields to `chemical_result()` without changing public package
  exports.

- [ ] **Step 8: Rebuild and verify GREEN**

  Run the focused result-level selection, then all
  `tests/test_chemical_equilibrium.py`.

### Task 5: Consolidate Documentation and Verify the Complete Migration

**Files:**

- Modify: `CONTEXT.md`
- Modify: `docs/designs/2026-07-21-private-reacting-phase-kernel.md`
- Modify: `docs/designs/2026-07-24-reactive-prework-handoff.md`
- Modify: `docs/designs/2026-07-27-grepe-homogeneous-chemical-layer.md`

**Interfaces:**

- Consumes: verified implementation and test evidence.
- Produces: one canonical description of the private GREPE homogeneous
  chemical layer and its remaining external prerequisites.

- [ ] **Step 1: Update canonical private-capability documentation**

  Record that the compiler now accepts redundant supplied reactions, validates
  cycles before reduction, owns homogeneous support, and recompiles accessible
  faces. Preserve the explicit limits:

  - no source-reference conversion without a real installed consumer;
  - no Provider reduced-component/boundary certificate;
  - no phase stability or coupled equilibrium;
  - no sensitivity or globality claim.

- [ ] **Step 2: Run formatting and static checks**

  ```bash
  git diff --check
  .venv/bin/python -m ruff check tests/test_chemical_equilibrium.py
  .venv/bin/python -m mypy src/epcsaft_equilibrium
  ```

- [ ] **Step 3: Run the focused chemistry suite**

  ```bash
  .venv/bin/python -m pytest -q tests/test_chemical_equilibrium.py
  ```

  Require zero failures with a freshly rebuilt extension.

- [ ] **Step 4: Run the full package suite**

  ```bash
  .venv/bin/python -m pytest -q
  ```

  Require zero failures. If a non-chemistry test fails, distinguish a real
  regression from an environment or stale-extension problem before changing
  code.

- [ ] **Step 5: Inspect the final code surface**

  ```bash
  git diff --stat
  git diff --check
  rg -n \
    "SEARCH_STABLE|CERTIFIED_EPS_GLOBAL|public.*chemical|eos_x_gamma|Reaktoro" \
    cpp src tests docs
  ```

  Confirm there is one compiler, one amount chart, and one Provider-backed
  chemical solver. Confirm forbidden claims or alternate implementations were
  not introduced.

- [ ] **Step 6: Run repository cleanup**

  ```bash
  bash "$HOME/.codex/hooks/codex-cleanup.sh" --repo-root .
  ```

  Remove only task-owned ignored artifacts through the hook's explicit removal
  mode, then verify Git status.

- [ ] **Step 7: Review requirements against the approved design**

  Check all eight acceptance items in
  `docs/designs/2026-07-27-grepe-homogeneous-chemical-layer.md`. Report exact
  remaining Provider, source-record, phase-pricing, boundary-oracle, and
  sensitivity limitations without weakening the implemented local claim.
