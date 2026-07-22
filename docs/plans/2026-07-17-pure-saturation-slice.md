# Pure-Component Saturation Slice Implementation Plan

Status: executed and accepted under `promotion-0018-equilibrium-pure-saturation-v1`; retained as implementation provenance

> **For agentic workers:** REQUIRED SUB-SKILL: use
> `superpowers:test-driven-development` for discrete contracts and the
> scientific-computing verification ladder for numerical claims.

**Goal:** Build and verify the first installed-artifact equilibrium slice for
local pure methane, ethane, and propane saturation.

**Architecture:** One Python package calls one native extension. The extension
holds one provider capsule, evaluates two explicit phase contexts, assembles a
three-variable Ipopt feasibility problem, and returns one immutable result
representation. The provider remains the only Helmholtz and nonlinear tensor
owner.

**Tech stack:** Python 3.13, scikit-build-core, CMake, pybind11, C++17, Ipopt,
pytest, Ruff, mypy, and the pinned provider wheel.

## Outcome Proof

**Intent:** Preserve the executable contract that created the accepted pure-component saturation capability without broadening it.
**Current Behavior:** `saturation` solves one local methane, ethane, or propane fixed-temperature boundary through the installed Provider and Ipopt.
**Expected Outcome:** The accepted operation, certificate, installed-artifact evidence, and receipt remain reproducible for their exact admitted inputs.
**Target Output:** One `SaturationResult` or typed `SaturationError` with separate solver, numerical, and physical evidence.
**Owner:** The Equilibrium package owns the saturation formulation and result; the installed Provider owns Helmholtz values and derivatives.
**Interface:** Public `saturation(model, temperature)` plus private native capsule and NLP evaluation seams.
**Cutover:** Completed by the implementation and promotion receipt identified in this repository.
**Replaced Path:** No predecessor equilibrium route was retained or redirected by this first capability slice.
**Evidence:** Focused tests, derivative checks, three scientific anchors, isolated-wheel inspection, and the accepted promotion receipt.
**Acceptance Proof:** The exact installed artifacts pass the receipt-bound pressure, potential, stability, domain, and confirmation certificates for all admitted components.
**Stop Criteria:** Stop on Provider mismatch, failed local certificate, source-domain violation, changed authority, or any request to widen the admitted scope.
**Avoid:** Do not add mixture equilibrium, phase discovery, globality claims, Provider equations, compatibility APIs, or release work.
**Risk:** Re-running this historical plan against a different Provider or silently widening its component/domain table would invalidate the accepted evidence.

## Implementation Boundaries

**Files To Create:** The executed slice created `pyproject.toml`, `CMakeLists.txt`, the native saturation owner, public API files, tests, anchors, validation script, architecture record, and candidate receipt.
**Files To Modify:** Only the package API, native module, saturation owner, tests, README, context, and architecture metadata named by the tasks.
**Files To Avoid:** Provider source, sibling repositories, unrelated equilibrium routes, release metadata, and user-owned IDE files.
**Source Of Truth:** The pure-saturation design, installed Provider SDK contract, retained source anchors, and `promotion-0018-equilibrium-pure-saturation-v1`.
**Read Path:** Public Provider header and capsule, retained anchor data, package design, and receipt-bound artifact records.
**Write Path:** Equilibrium source, tests, documentation, architecture metadata, and receipt locations explicitly named by this plan.
**Integration Points:** `epcsaft.native_sdk.v1`, the pybind11 extension, Ipopt, Python result conversion, and installed-wheel verification.
**Migration Or Cutover:** The accepted receipt completed cutover for only the named pure-component local boundary.
**Replaced Path Handling:** No displaced route exists; negative-space checks must continue to prevent duplicate Provider or saturation owners.
**Acceptance Proof Gate:** Reproduce the receipt-bound installed-artifact tests and scientific anchors without changing fingerprints, tolerances, or domains.

## Global constraints

- Provider commit: `02a03a440cf53f3b1d4304556a19a54c79c70179`.
- Provider wheel SHA-256:
  `e83f4b108d2df73888f7768f99f1c54b01ad86541b55afc1e43bbef2e1fb8f93`.
- Provider receipt SHA-256:
  `7c097c3ec8742d1714ea435973c605d5c18290e27729c1a12b3852926e61c399`.
- Build only against the installed public provider header and call only the
  public Python API and `epcsaft.native_sdk.v1` capsule.
- Keep one distribution, one Python package, one extension, and one native
  target.
- Do not add a backend selector, compatibility API, parameter table, provider
  kernel, workflow, release, or broader equilibrium route.

---

### Task 1: Package and capsule seam

**Files:**

- Create `pyproject.toml` and `CMakeLists.txt`.
- Create `src/epcsaft_equilibrium/__init__.py`.
- Create `cpp/src/module.cpp`.
- Create `tests/test_saturation.py`.

**Interfaces:**

```cpp
py::dict sdk_info(const py::capsule& capsule);
py::dict evaluate_nlp(
    const py::capsule& capsule,
    double temperature_k,
    const std::string& fingerprint,
    const std::array<double, 3>& variables,
    const std::array<double, 3>& multipliers
);
```

**Use Cases:**

- Establish one visible package/capsule seam while proving no duplicate Provider owner enters the wheel.

- [ ] Write a failing package test that imports the extension and asks for the
  capsule ABI, table size, result size, model context, and evaluation entry.
- [ ] Run the focused test and confirm that the extension is missing.
- [ ] Add the minimal scikit-build and CMake metadata. Locate the provider
  header from the installed `epcsaft` package and link only Ipopt and Python.
- [ ] Implement capsule-name and ABI validation, size checks, fingerprint
  checks, and one phase evaluation wrapper.
- [ ] Rebuild and run the focused test. Delete the Python model before a second
  capsule evaluation to prove capsule-owned lifetime.

### Task 2: Exact local formulation

**Files:**

- Create `cpp/src/saturation.hpp` and `cpp/src/saturation.cpp`.
- Modify `cpp/src/module.cpp`.
- Extend `tests/test_saturation.py`.

**Interfaces:**

```cpp
struct NlpEvaluation {
    double objective;
    std::array<double, 3> objective_gradient;
    std::array<double, 3> constraints;
    std::array<double, 9> jacobian;
    std::array<double, 6> lagrangian_hessian_lower;
};

NlpEvaluation evaluate_saturation_nlp(
    const ProviderContext& provider,
    double temperature_k,
    const std::array<double, 3>& variables,
    const std::array<double, 3>& multipliers
);
```

**Use Cases:**

- Produce derivative evidence that the exact local formulation, rather than a numerical surrogate, satisfies the acceptance basis.

- [ ] Write a failing parameterized derivative test at representative and
  boundary-near methane, ethane, and propane states.
- [ ] Confirm objective and gradient zeros, residual directional differences
  against the Jacobian, and finite-difference Jacobian-transpose directional
  differences against the exact Lagrangian Hessian.
- [ ] Implement only the `(n,V)` tensor transformations and saturation
  constraints documented in the design.
- [ ] Run the derivative test and confirm the declared tolerances without
  changing expected values to match the implementation.

### Task 3: Ipopt solve and public result

**Files:**

- Modify `cpp/src/saturation.hpp`, `cpp/src/saturation.cpp`, and
  `cpp/src/module.cpp`.
- Create `src/epcsaft_equilibrium/_api.py`.
- Modify `src/epcsaft_equilibrium/__init__.py`.
- Extend `tests/test_saturation.py`.

**Public interfaces:**

```python
@dataclass(frozen=True)
class PhaseState:
    amount_mol: float
    volume_m3: float
    molar_density_mol_m3: float
    pressure_pa: float
    chemical_potential_over_rt: float

@dataclass(frozen=True)
class SaturationDiagnostics:
    solver_converged: bool
    solver_status: str
    iterations: int
    numerical_converged: bool
    physical_accepted: bool
    pressure_relative_residual: float
    chemical_potential_absolute_residual: float
    phase_density_distance: float
    exact_derivatives: bool
    globality_certificate: bool

@dataclass(frozen=True)
class SaturationResult:
    temperature_k: float
    saturation_pressure_pa: float
    parameter_fingerprint: str
    vapor: PhaseState
    liquid: PhaseState
    diagnostics: SaturationDiagnostics

class SaturationError(RuntimeError):
    diagnostics: Mapping[str, object]

def saturation(
    model: epcsaft.EPCSAFT,
    temperature: pint.Quantity,
) -> SaturationResult:
    """Solve and certify one local pure-component saturation boundary."""
```

**Use Cases:**

- Return a typed result whose visible solver, numerical, and physical fields prevent Ipopt success from acting as acceptance proof.

- [ ] Write failing tests for input units, approved fingerprints, source
  domains, unsupported models, and topology collapse.
- [ ] Write a failing representative ethane solve test that requires separate
  solver, numerical, and physical acceptance plus exact derivatives and no
  globality certificate.
- [ ] Implement the TNLP with exact objective gradient, constraint Jacobian,
  and Lagrangian Hessian. Record Ipopt status, iteration count, constraint
  violation, bounds, seed, and callback errors.
- [ ] Implement deterministic seed discovery and two perturbed confirmation
  solves. Require agreement before setting numerical convergence.
- [ ] Implement the immutable Python result types and strict input conversion.
- [ ] Run the focused tests and inspect every status layer.

### Task 4: Scientific and installed-artifact evidence

**Files:**

- Create `data/reference/pure_saturation_anchors.csv` with one retained NIST
  row per admitted component, source URL, units, and lab comparison value.
- Create `scripts/validate_saturation.py` for the three-row scientific check.
- Extend `tests/test_saturation.py` with one parameterized compact anchor test.

**Use Cases:**

- Verify the installed artifact against retained scientific evidence and detect copied or duplicate Provider implementation content.

- [ ] Freeze methane, ethane, and propane expected rows from retained lab
  evidence before running the new implementation.
- [ ] Run the source-installed extension test against the pinned provider
  artifact.
- [ ] Build the equilibrium wheel in an isolated environment with no sibling
  source path, install both wheels into a second isolated environment, and run
  the public tests and validation script.
- [ ] Inspect wheel contents, build commands, dynamic links, and symbols. Prove
  that the wheel contains no provider source, implementation objects, private
  headers, or duplicated provider equation symbols.
- [ ] Record wheel SHA-256 and the exact verification commands.

### Task 5: Architecture, receipt, review, and completion

**Files:**

- Update `README.md` and `CONTEXT.md`.
- Create `ARCHITECTURE.yaml`.
- Create a candidate promotion receipt under `receipts/` after the
  implementation commit exists.

**Use Cases:**

- Complete the receipt-bound migration and leave capability acceptance visible without implying a broader cutover.

- [ ] Record the public exports, native target, dependencies, admitted
  fingerprints and domains, derivative order, local certificate, artifact
  hashes, evidence commands, and exclusions.
- [ ] Run Ruff, mypy, CMake build, focused tests, scientific validation, and
  isolated installed-wheel verification.
- [ ] Commit the implementation on local `main` without pushing.
- [ ] Ask an independent reviewer to inspect requirements, ownership,
  derivatives, physical gates, tests, and negative space. Fix accepted findings
  and rerun the affected verification.
- [ ] Create the candidate receipt against the implementation commit. Mark
  authority acceptance and user promotion approval pending.
- [ ] Run the full verification commands, repository cleanup audit, and Git
  status check. Commit the receipt and review record locally.
