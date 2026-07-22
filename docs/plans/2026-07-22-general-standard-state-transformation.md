# General Standard-State Transformation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the retired water-only Provider reference tail with the installed Provider neutral-reference basis and transform complete charge-neutral source reaction records into the Provider Helmholtz basis without assigning single-ion references.

**Architecture:** `ProviderContext` owns exact SDK validation and returns the Provider neutral basis, its dimensionless log-fugacity contractions, salt-free reference composition, reference scalars, identities, and derivative-availability mask. `chemical_equilibrium.cpp` expresses each charge-neutral reaction in that basis, combines the Provider contractions with explicit source activity-scale factors, and returns Provider-basis `lnK` values plus residual and availability evidence. The existing private water-specific transformation is removed rather than retained as a compatibility path.

**Tech Stack:** C++17, pybind11 private bindings, installed `epcsaft` native SDK v1, pytest, Ipopt 3.11.9.

## Global Constraints

- Consume only the immutable Provider wheel with SHA-256 `a8b6376193301673429a8d8b648896b6881c6f693ca4d8fc3f6a9d2d16f3b39c`.
- Keep every new route private; add no public solve function, result family, selector, or compatibility shim.
- Accept only dimensionless source constants with exact component, charge, T/P, reference, and fingerprint identities.
- Transform only charge-neutral reactions representable by the complete Provider neutral basis; never construct single-ion standards.
- Report Provider derivative availability exactly; do not infer unavailable T/P/composition/parameter derivatives.
- Preserve local Ipopt/KKT/globality classifications; this leaf changes reference construction, not solver authority.
- Do not implement #46 unless an exact source-complete MEA Provider artifact and all ledger blockers are resolved.

---

### Task 1: Consume the neutral-reference SDK tail

**Files:**
- Modify: `cpp/src/provider.hpp`
- Modify: `cpp/src/provider.cpp`
- Test: `tests/test_chemical_equilibrium.py`

**Interfaces:**
- Consumes: `epcsaft_evaluate_neutral_reference_v1` from the installed SDK.
- Produces: `NeutralReferenceEvaluation ProviderContext::evaluate_neutral_reference(double temperature_k, double pressure_pa)`.

- [x] Write a failing installed-Provider test that requests the Figiel water/sodium/chloride neutral basis and checks exact component order, charge-neutral rows, salt-free composition, basis identity, fingerprint, T/P binding, convergence, and unavailable derivatives.
- [x] Run the focused test and confirm the current extension cannot build against the new SDK tail.
- [x] Replace the obsolete standard-reference ABI adapter with exact-capacity neutral-reference evaluation and fail-closed scalar, buffer, rank, identity, and derivative-mask validation.
- [x] Rebuild against the immutable wheel and run the focused test.

### Task 2: Transform arbitrary charge-neutral reaction records

**Files:**
- Modify: `cpp/src/chemical_equilibrium.hpp`
- Modify: `cpp/src/chemical_equilibrium.cpp`
- Modify: `cpp/src/chemical_equilibrium_bindings.cpp`
- Test: `tests/test_chemical_equilibrium.py`

**Interfaces:**
- Consumes: ordered reactions `nu`, dimensionless source `lnK`, explicit per-component log activity-scale factors, and `NeutralReferenceEvaluation`.
- Produces: Provider-basis `lnK`, standard offsets, reaction-to-neutral-basis coefficients, representation residual, and derivative availability.

- [x] Write one analytic multi-reaction charged test covering basis rescaling/reordering, species permutation, and gauge-free reconstruction.
- [x] Write one consolidated negative-path test for non-neutral reactions, incomplete basis span, source/reference identity mismatch, non-finite scale factors, and unsupported derivative claims.
- [x] Run both tests and confirm the private generic transform is absent.
- [x] Implement one rank-revealing full-column solve for reaction coordinates in the neutral basis, independently recompute the original-coordinate residual, and combine source scale and Provider contraction terms.
- [x] Replace the water-only private binding and caller with the generic transform; remove obsolete ABI test structures and special-case fields.
- [x] Run the focused chemical-equilibrium suite.

### Task 3: Adjudicate the MEA sentinel gate

**Files:**
- Read: `data/reference/mea_speciation_input_ledger.json`
- Read: installed Provider wheel manifest and catalog
- Update: GitHub issue #46 only with immutable blocker evidence if the gate remains closed

**Interfaces:**
- Consumes: exact Provider artifact identity and the machine-validated MEA ledger.
- Produces: either an authorized source-complete private sentinel input or a fail-closed blocker report.

- [x] Run the MEA ledger validator with `--require-executable` and retain its nonzero result.
- [x] Confirm the installed Provider catalog contains every ordered MEA true species before writing any sentinel code.
- [x] If either gate fails, do not add a chemistry-specific fixture or fake reference; report #46 blocked with exact missing records and leave it open.
- [ ] If both gates pass, start a separate TDD cycle for the one-state generic Provider solve.

### Task 4: Verify and publish #34

**Files:**
- Modify: `docs/phase-equilibrium.md` only if the internal reference contract description changed.

- [x] Build against the exact Provider wheel.
- [x] Run focused chemical-equilibrium tests, full pytest, Ruff, strict mypy, diff checks, and the repository cleanup audit.
- [x] Review the staged diff for one direct private path and no obsolete water-only compatibility surface.
- [ ] Commit and push #34, open its PR against `codex/chemical-equilibrium`, and merge only after the exact head is clean and verified.
