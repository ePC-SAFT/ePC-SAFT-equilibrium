# ePC-SAFT Equilibrium

`epcsaft-equilibrium` owns bounded equilibrium formulations over an installed
`epcsaft` EOS. The accepted first capability solves a pure-component
saturation boundary for the EOS-approved methane, ethane, or propane
model under promotion receipt `promotion-0018-equilibrium-pure-saturation-v1`.
The package also contains a non-authoritative local candidate for one bounded
neutral methane/ethane HELD `T,P,z` calculation.

```python
import epcsaft
from epcsaft_equilibrium import saturation

parameters = epcsaft.Parameters.from_catalog(
    "gross-2001-methane-ethane", components=("methane",), version=1
)
mixture = epcsaft.Mixture(parameters)

result = saturation(mixture, 150 * epcsaft.unit_registry.kelvin)
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
finite values, and the EOS fingerprint. It carries no globality or phase
discovery certificate.

## Local neutral HELD candidate

```python
import epcsaft
from epcsaft_equilibrium import tp_flash

parameters = epcsaft.Parameters.from_catalog(
    "gross-2001-methane-ethane",
    components=("methane", "ethane"),
    version=1,
)
mixture = epcsaft.Mixture(parameters)

result = tp_flash(
    mixture,
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
retains invalid-input, EOS, exhausted-search, scope, and indeterminate
diagnostics. Every result reports `globality_certificate="not_guaranteed"`;
the accepted pure-saturation authority above is unchanged. HELD diagnostics
classify solver, numerical, and physical evidence independently as `passed`,
`failed`, or `not_adjudicated`; those axes do not change the globality claim.

The candidate fingerprint is
`sha256:3a840001adcb8b82f44e48307ad61e566f6a65d9b82d8312299a439dbce09195`.
The retained Pereira source is the permanent-lab Markdown at commit
`13ce345b6dcc41d399bb2a4c7b9bedb18f74b45b`, blob
`dde7f02d4c93cce86804a8e6b62d37602990ac21`; it is provenance, not a build,
test, validation, or runtime dependency.

The non-production development route extends the same public `tp_flash`
operation to qualifying installed EOS electrolyte SDKs and the Perdomo
HELD2 Steps 1--10 controller. Dispatch is capability-driven rather than case-
or component-name-driven. It includes
the deterministic pressure-root envelope, DIRECT-L Stage I, HiGHS Stage-II
Problem (64), and deterministic capped-multistart Ipopt refinement. One typed
Steps 1--10 workflow is shared by installed and manufactured problems; Stage
II has one major loop, one Eq. (66) decision owner, and one deterministic
bounded Step-5 start set per major. Replaced controllers and duplicate
manufactured orchestration are deleted rather than retained as compatibility
paths. The experimental public runtime is preserved by tag
`archive/held2-pre-strategy-2026-07-21` for strategy review.

The HELD2 coordinate contract is hybrid only where the source requires it:
Steps 1--9 use linear modified compositions and linear material balances.
The finite \(10^{-10}\) Step-1 floor regularizes those searches; it is not a
physical minimum. Step 10 alone refines a strictly positive charged trace
fraction in bounded \(\log_{10}\) coordinates, preserves every other Step-1
constraint, and then recertifies the reconstructed linear balances and
chemical potentials.

The homogeneous reference search reports detected pressure-root
accounting separately as `root_completeness`. Its installed evidence reports
`root_completeness="not_proven"`: detecting and refining the retained roots is
not proof that every root exists in the finite domain. Root completeness is
independent of solver, numerical, physical, search-completeness, predictive,
and globality status.

The installed public Perdomo Table 3 evidence returns one accepted
homogeneous phase and a cross-EOS source-topology disagreement. It does not admit
electrolyte LLE or reproduce Perdomo's SAFT-gamma-Mie numerical endpoints. The
admission gate remains one source-complete installed ePC-SAFT case that reaches
and passes Stage II and Stage III with two distinct liquids. It requires an
exact corrected EOS artifact and independent installed-artifact evidence.

The HELD2 algorithm is specified in
`docs/designs/2026-07-24-held2-paper-algorithm.md` and implemented by
`docs/plans/2026-07-24-held2-paper-rewrite.md`. It assigns
deterministic pressure-root enumeration to homogeneous and trial-composition
density topology, NLopt DIRECT-L to the reduced Stage-I TPD search, HiGHS to
the Stage-II upper LP, deterministic capped-multistart exact-Hessian Ipopt to
Step 5, and exact-Hessian Ipopt to Problem (67) in Step 8. SLSQP is not the
default replacement for Ipopt. The development route implements the Steps
1--10 assignments and links pinned NLopt 2.11.0 and HiGHS 1.15.1 with Ipopt.
Installed two-liquid evidence and capability admission remain outside the
landed scope.

The development-only native diagnostic is the live-progress route. After
exporting the installed EOS model configuration, run:

```bash
build/epcsaft-equilibrium-diagnostic \
  --model-config MODEL.json --temperature 298.15 --pressure 2508 \
  --feed 0.8321050353538131,0.08394748232309347,0.08394748232309347 \
  --trace
```

The trace shows reference roots, Stage-I evaluations, certificates, and exact
failure reasons as they occur. Enabling it does not change results, gates,
budgets, or the finite-search globality label.

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
compiles against the declaration header installed by the EOS wheel. It
does not link EOS implementation symbols, compile EOS sources, or
import private EOS modules.

Source builds require Python 3.13, CMake, a C++17 compiler, pkg-config, Ipopt,
network or populated FetchContent caches for the pinned NLopt 2.11.0 and
HiGHS 1.15.1 archives and the header-only Boost 1.88.0 archive, and the
non-editable EOS wheel installed in the build environment. The
local candidate gate hashes the exact EOS wheel before creating an
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
`docs/designs/2026-07-24-held2-paper-algorithm.md`; the canonical HELD2 execution plan
and landed task record is
`docs/plans/2026-07-24-held2-paper-rewrite.md`. The D-028 homogeneous
reacting-phase foundation is documented separately in
`docs/designs/2026-07-21-private-reacting-phase-kernel.md`. It currently has
Belov trace evidence, installed-EOS-manufactured evidence, and one
source-complete Held/IAPWS water self-ionization value case. The public typed
`chemical_equilibrium` operation reports a local fixed-`T,P` homogeneous value
and optionally returns exact conditioned derivatives with respect to compiled
balance totals, final EOS-basis `ln(K)`, and pressure. Value-only and
value-plus-Jacobian responses are explicit. Requested unsupported columns fail
closed. A source standard state's reference pressure remains immutable
provenance, while the EOS reference is re-evaluated and the transformed
EOS-basis record is bound at each actual trial pressure. Exact
source-reference pressure derivatives are included in the returned pressure
column when the installed derivative tail and its branch certificates are
available. Typed active-parameter requests consume only coordinates advertised
by the installed EOS and require one atomic callback to supply the
active-model Helmholtz, packing, pressure, chemical-potential, and
neutral-reference derivative blocks. Unsupported or incomplete requests fail
closed; no derivative is approximated. Results bind species/parameter order, units, chart topology,
EOS fingerprint, installed distribution RECORD fingerprints, and EOS
SDK ABI identity. The reported `condition_number_inf` is the deterministic
infinity-norm condition number of the four-pass row/column-equilibrated KKT
system and is the one numerical conditioning metric used by the sensitivity
gate; no second or unscaled condition metric is reported. The operation has no predictive admission, coupled
phase-equilibrium claim, or globality proof.
Accepted receipt
`promotion-0018-equilibrium-pure-saturation-v1` makes this repository the
production owner of that exact local boundary capability. One local boundary
solve is not a phase-discovery or global-stability proof. The local HELD
candidate has authority effect `none` until separate review, validation,
EOS-tail promotion, equilibrium promotion, and explicit user approval.
