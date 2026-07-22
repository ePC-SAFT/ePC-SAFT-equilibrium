# Private Reacting-Phase Kernel Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the D-028 private homogeneous fixed-`T,P` reacting-phase kernel with exact compiled invariants, electroneutral coordinates, native Ipopt, certificates, and manufactured evidence.

**Architecture:** Keep one private C++ owner beside the existing native equilibrium code. A small underscored pybind test seam exercises typed C++ behavior without adding a public package export. The installed Provider remains the sole owner of the mechanical phase value and nonlinear tensors; manufactured analytic blocks supply independent closed-form oracles.

**Tech Stack:** C++17, pybind11, installed `epcsaft.native_sdk.v1`, Ipopt 3.11.9, Python 3.13, pytest.

## Global Constraints

- Start from clean Equilibrium commit `bb118bec4e6abdec417417e25a9b56d03491de4b` on `codex/chemical-equilibrium`.
- D-028 is sequencing authority only and has no capability-admission effect.
- Add no public solve function, export, result family, selector, compatibility shim, backend registry, EOS equation, receipt, release, or authority change.
- Do not import Provider source, copy Provider equations, or use the ordinary checkout.
- Treat every installed-Provider reaction case as manufactured and nonpredictive.
- Production derivatives are analytic/Provider exact; directional finite differences exist only in tests.
- `globality` remains `not_guaranteed`; missing evidence remains `not_adjudicated`.

---

### Task 1: Typed reaction-system compiler

**Files:**
- Create: `cpp/src/chemical_equilibrium.hpp`
- Create: `cpp/src/chemical_equilibrium.cpp`
- Create: `cpp/src/chemical_equilibrium_bindings.cpp`
- Modify: `cpp/src/module.cpp`
- Modify: `CMakeLists.txt`
- Test: `tests/test_chemical_equilibrium.py`

**Interfaces:**
- Consumes: ordered species/provider metadata, `B`, `nu`, feed, `lnK`, source/reference records, `T`, `P`, and completeness flags.
- Produces: private `CompiledReactionSystem`, minimum-norm `g_ref`, reduced independent balance rows, and underscored compiler evidence.

- [x] Write parameterized failing tests for valid `A <-> B`, gauge reconstruction, and every compiler rejection class: order, charge, fingerprint, rank, nonconservation, reaction dependence, incomplete `lnK`, reference binding, and non-neutral feed.
- [x] Run `pytest -q tests/test_chemical_equilibrium.py -k compiler` and confirm failure because the private seam is absent.
- [x] Implement dense rank/solve helpers and the immutable compiler with scale-aware residual checks.
- [x] Bind only underscored native test functions and add the new sources to the existing native target.
- [x] Re-run the focused compiler tests and retain only claim-bearing cases.
- [x] Commit the design/compiler checkpoint after fresh focused and existing-suite verification.

### Task 2: Exact positive electroneutral chart

**Files:**
- Modify: `cpp/src/chemical_equilibrium.hpp`
- Modify: `cpp/src/chemical_equilibrium.cpp`
- Modify: `cpp/src/chemical_equilibrium_bindings.cpp`
- Test: `tests/test_chemical_equilibrium.py`

**Interfaces:**
- Consumes: ordered integer charges and unconstrained chart coordinates.
- Produces: strictly positive amounts, exact charge, inverse coordinates, Jacobian, and one Hessian per physical amount.

- [x] Write failing tests for neutral logs, multi-cation/multi-anion softmax shares, exact charge, round trip, species permutation, amount scale, and trace-floor classification.
- [x] Add test-only centered directional checks for the Jacobian and all amount Hessians; confirm RED before implementation.
- [x] Implement reference-category softmax value/first/second derivatives and exact `Q`/neutral exponential chain rules.
- [x] Re-run chart tests at ordinary and trace scales, then the compiler tests.
- [x] Commit the chart checkpoint after fresh verification.

### Task 3: Phase block, max-min initializer, TNLP, and certificates

**Files:**
- Modify: `cpp/src/chemical_equilibrium.hpp`
- Modify: `cpp/src/chemical_equilibrium.cpp`
- Modify: `cpp/src/chemical_equilibrium_bindings.cpp`
- Test: `tests/test_chemical_equilibrium.py`

**Interfaces:**
- Consumes: `CompiledReactionSystem`, a differentiable phase block, pressure, volume/domain bounds, trace floor, and fixed solver profile.
- Produces: private solve evidence, exact KKT residual/Jacobian, reaction affinities, conservation/charge/pressure/domain certificates, and reduced-Hessian status.

- [x] Write failing analytic solve tests for `A <-> B`, `2A <-> B`, and `A <-> C+ + D-`, with expected compositions derived independently in the tests.
- [x] Write failing invariance tests for gauge, species order/permutation, amount scale, and reaction-basis transformations.
- [x] Write failing negative tests for max-min infeasibility, trace-bound terminals, callback/domain failures, non-success Ipopt statuses, indeterminate states, and synthetic false-success inputs.
- [x] Implement the analytic ideal phase block and exact physical-to-chart objective chain rule.
- [x] Implement the linear max-min TNLP and independently recompute strict-positive feasibility.
- [x] Implement one homogeneous TNLP, direct `lambda=1` attempt, bounded continuation fallback, and mandatory final `lambda=1` acceptance.
- [x] Implement independent postsolve recomputation, exact KKT assembly, and reduced-Hessian local-minimum adjudication.
- [x] Run focused analytic and negative tests, then all chemical-equilibrium tests.
- [x] Commit the native kernel checkpoint after fresh verification.

### Task 4: Installed Provider manufactured integration

**Files:**
- Modify: `cpp/src/chemical_equilibrium.cpp`
- Modify: `cpp/src/chemical_equilibrium_bindings.cpp`
- Test: `tests/test_chemical_equilibrium.py`

**Interfaces:**
- Consumes: exact installed Provider capsule, fingerprint, component order/charges, electrolyte phase callback, molar-volume callback, and packing callback.
- Produces: one private manufactured/nonpredictive reacting solve and exact transformed derivative/domain evidence.

- [x] Write a failing installed-artifact test using the immutable Perdomo Table-3 Provider wheel and synthetic reaction reference data.
- [x] Add test-only directional gradient/Hessian checks around one strict interior Provider state.
- [x] Implement the Provider phase adapter through `ProviderContext`, rejecting order/fingerprint/domain/packing mismatches before acceptance.
- [x] Verify the installed case is labeled `manufactured_nonpredictive`, consumes all declared Provider callbacks, and makes no source-backed chemical claim.
- [x] Commit the Provider-integration checkpoint after focused installed-artifact verification.

### Task 5: Documentation, static checks, and full evidence

**Files:**
- Modify: `CONTEXT.md`
- Modify: `ARCHITECTURE.yaml`
- Modify: `docs/phase-equilibrium.md`
- Modify: `docs/superpowers/plans/2026-07-21-private-reacting-phase-kernel.md`

**Interfaces:**
- Consumes: final tested implementation and D-028.
- Produces: truthful private-candidate documentation and exact verification record without promotion semantics.

- [x] Update package documentation to name the private foundation and preserve the accepted/public boundaries.
- [x] Run focused pytest, full pytest, Ruff, mypy, native build warnings, and the installed-artifact case from the exact Provider wheel.
- [x] Inspect `git diff`, staged diff, public exports, dependency edges, and generated artifacts.
- [x] Run `bash "$HOME/.codex/hooks/codex-cleanup.sh" --repo-root .` and remove only task-owned ignored artifacts if reported.
- [ ] Mark plan checkboxes from actual evidence, commit a stable checkpoint, push `codex/chemical-equilibrium`, and report commit/tree/artifact/test evidence plus unresolved scientific blockers.
