# ePC-SAFT Equilibrium

`epcsaft-equilibrium` owns bounded equilibrium formulations over an installed
`epcsaft` provider. The accepted first capability solves a pure-component
saturation boundary for the provider-approved methane, ethane, or propane
model under promotion receipt `promotion-0018-equilibrium-pure-saturation-v1`.
The package also contains a non-authoritative local candidate for one bounded
neutral methane/ethane HELD `T,P,z` calculation.

```python
import epcsaft
from epcsaft_equilibrium import saturation

parameters = epcsaft.ParameterBundle.from_catalog(
    "gross-2001-methane-ethane", version=1
).select(("methane",))
model = epcsaft.EPCSAFT(parameters)

result = saturation(model, 150 * epcsaft.unit_registry.kelvin)
print(result.saturation_pressure_pa)
print(result.vapor.molar_density_mol_m3, result.liquid.molar_density_mol_m3)
```

The operation requires a Pint temperature quantity and accepts these source
domains: methane 97 to 300 K, ethane 90 to 305 K, and propane 85 to 523 K. A
temperature inside a parameter domain does not guarantee a two-phase boundary.
The operation raises `SaturationError` when Ipopt cannot produce a numerically
confirmed state that passes the local physical checks.

The result reports solver convergence, numerical confirmation, and physical
acceptance as separate fields. It checks pressure and chemical-potential
equality, positive distinct phase states, local mechanical stability, bounds,
finite values, and the provider fingerprint. It carries no globality or phase
discovery certificate.

## Local neutral HELD candidate

```python
import epcsaft
from epcsaft_equilibrium import tp_flash

parameters = epcsaft.ParameterBundle.from_catalog(
    "gross-2001-methane-ethane", version=1
).select(("methane", "ethane"))
model = epcsaft.EPCSAFT(parameters)

result = tp_flash(
    model,
    243.61 * epcsaft.unit_registry.kelvin,
    6.691 * epcsaft.unit_registry.megapascal,
    (0.5627, 0.4373),
)
print(result.phases, result.phase_fractions)
print(result.diagnostics.outcome, result.diagnostics.globality_certificate)
```

The controller admits only the reviewed binary fingerprint inside the
rectangular May et al. (2015), Table 5 source domain: 203.22--243.61 K,
2.124--6.885 MPa, and methane feed fraction 0.4661--0.66705. It returns one or
two phases after the declared finite Stage-I/II/III search. `FlashError`
retains invalid-input, provider, exhausted-search, scope, and indeterminate
diagnostics. Every result reports `globality_certificate="not_guaranteed"`;
the accepted pure-saturation authority above is unchanged. HELD diagnostics
classify solver, numerical, and physical evidence independently as `passed`,
`failed`, or `not_adjudicated`; those axes do not change the globality claim.

The candidate fingerprint is
`sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286`.
The retained Pereira source is the permanent-lab Markdown at commit
`13ce345b6dcc41d399bb2a4c7b9bedb18f74b45b`, blob
`dde7f02d4c93cce86804a8e6b62d37602990ac21`; it is provenance, not a build,
test, validation, or runtime dependency.

An archived non-production WIP subject extended the same public `tp_flash`
operation to qualifying installed Provider electrolyte SDKs and the Perdomo
HELD2 Stage-I/II/III controller. That dispatch was capability-driven rather
than case- or component-name-driven, but it is not part of current `main`.
Current `main` retains a private HELD2 development implementation. It includes
the deterministic pressure-root envelope, DIRECT-L Stage I, HiGHS Stage-II
Problem (64), and deterministic Stage-II basin discovery followed by
exact-Hessian Ipopt refinement. The experimental public runtime is preserved
by tag
`archive/held2-pre-strategy-2026-07-21` for strategy review.

The archived homogeneous reference search reports detected pressure-root
accounting separately as `root_completeness`. Its installed evidence reports
`root_completeness="not_proven"`: detecting and refining the retained roots is
not proof that every root exists in the finite domain. Root completeness is
independent of solver, numerical, physical, search-completeness, predictive,
and globality status.

The archived installed public Perdomo Table 3 evidence returns one accepted
homogeneous phase and a cross-EOS source-topology disagreement. It does not admit
electrolyte LLE or reproduce Perdomo's SAFT-gamma-Mie numerical endpoints. The
future gate remains one source-complete installed ePC-SAFT case that reaches
and passes Stage II and Stage III with two distinct liquids. It requires a new
bounded implementation assignment and exact corrected Provider artifact before
runtime work resumes.

The HELD2 solver strategy is canonized in
`docs/plans/2026-07-21-perdomo-held2-solver-strategy.md`. It assigns
deterministic pressure-root enumeration to homogeneous and trial-composition
density topology, NLopt DIRECT-L to the reduced Stage-I TPD search, HiGHS to
the Stage-II upper LP, global basin discovery plus exact-Hessian Ipopt to the
Stage-II lower problem, and exact-Hessian Ipopt to Stage III. SLSQP is not the
default replacement for Ipopt. Current `main` implements the Stage-I and
Stage-II assignments as private development routes and links pinned NLopt
2.11.0 and HiGHS 1.15.1 with Ipopt. Stage-III hardening, installed two-liquid
evidence, public electrolyte dispatch, and capability admission remain outside
that landed scope.

The earlier fixed-two-phase campaign remains `NON_ADMISSION` under the frozen
May `3*u_c` composition contract: 12 of 17 rows passed, rows 002/009/010/011
were solved model/data misses, and row 012 had no package-accepted local state.
HELD does not erase or reclassify those results. The later installed-artifact
HELD campaign is also stable `NON_ADMISSION`: it returned two one-phase
results, thirteen fail-closed third-candidate results, and three exhausted
searches, so no two-phase composition comparison was evaluable. Permanent-lab
review requires a focused controller-lifecycle redesign before another runtime
candidate. The current wheel remains immutable and unpromoted; the design-only
delta changes no runtime, tolerance, resource, public type, or globality claim.
That neutral redesign and the archived Perdomo HELD2 WIP remain provenance
while Stage-III hardening and installed two-liquid evidence await separate
assignments.

## Native boundary

The extension calls `epcsaft.native_sdk(model)` and retains one model-bound
`epcsaft.native_sdk.v1` capsule while Equilibrium evaluates phase contexts.
DIRECT-L explores the reduced HELD2 envelopes, HiGHS solves Problem (64), and
Ipopt refines smooth Stage-II and Stage-III NLPs. The HELD candidate consumes
the reviewed mixture value/gradient/Hessian tail;
the pure route continues to consume the accepted prefix. The extension
compiles against the declaration header installed by the provider wheel. It
does not link provider implementation symbols, compile provider sources, or
import private provider modules.

Source builds require Python 3.13, CMake, a C++17 compiler, pkg-config, Ipopt,
network or populated FetchContent caches for the pinned NLopt 2.11.0 and
HiGHS 1.15.1 archives, and the non-editable provider wheel installed in the
build environment. The
local candidate gate hashes the exact provider wheel before creating an
isolated build environment. Candidate wheels are retained as read-only files
under a commit-bound `artifacts/equilibrium-neutral-held-v1/<commit>/`
directory and are never overwritten by a correction. The exact candidate
commit, declared `SOURCE_DATE_EPOCH`, artifact hashes, review status, and
replayable commands are recorded in
`receipts/promotion/promotion-0002-neutral-held-v1-candidate.yaml`.
One clean build does not support an archive-byte reproducibility claim.
The current corrected subject is implementation commit `8318e755`, retained
under `artifacts/equilibrium-neutral-held-v1/8318e75/`; the rejected
`549162a3` artifact remains immutable provenance.

```text
uv run --isolated --no-project --python 3.13 \
  --with "$PINNED_PROVIDER_WHEEL" --with scikit-build-core --with pybind11 \
  -- uv build --no-build-isolation --wheel .
```

Run the compact package proof and scientific anchors with:

```text
pytest -q
python scripts/validate_saturation.py
```

The canonical documentation map is `docs/phase-equilibrium.md`. Detailed
formulation owners are
`docs/designs/2026-07-17-pure-saturation-slice.md`,
`docs/designs/2026-07-17-neutral-held-v1.md`, and
`docs/designs/2026-07-21-perdomo-held2.md`; the canonical HELD2 execution plan
and landed task record is
`docs/plans/2026-07-21-perdomo-held2-solver-strategy.md`. Migration receipt
`promotion-0018-equilibrium-pure-saturation-v1` makes this repository the
production owner of that exact local boundary capability. One local boundary
solve is not a phase-discovery or global-stability proof. The local HELD
candidate has authority effect `none` until separate review, validation,
provider-tail promotion, equilibrium promotion, and explicit user approval.
