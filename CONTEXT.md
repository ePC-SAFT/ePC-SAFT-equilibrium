# Equilibrium Repository Context

This repository is the production owner of the accepted pure-component
saturation slice recorded by migration receipt
`promotion-0018-equilibrium-pure-saturation-v1`. The accepted scope is one
local fixed-temperature methane, ethane, or propane saturation boundary over a
compatible installed provider artifact.

`governance_doctrine_revision: 2`

Canonical local doctrine: `../ePC-SAFT-organization/GOVERNANCE.md`.

The lab copy is non-authoritative provenance and research for this accepted
slice. The receipt does not imply phase discovery, global stability, mixture
equilibrium, release publication, or any capability outside its exact subject.

The repository now also contains the corrected complete local
`neutral-held-v1` candidate at implementation commit
`8318e755d4a8e490822fdf7bb2685d8c5af6436c`, tree
`3e8c98a13f1daca975b33c26fce3d143a1f34493`. It is limited to the reviewed
Gross--Sadowski methane/ethane fingerprint
`sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286`
and the audited May et al. (2015), Table 5 rectangle: 203.22--243.61 K,
2.124--6.885 MPa, and methane feed 0.4661--0.66705.

The sole public `tp_flash` operation runs the bounded Pereira HELD
Stage-I/II/III controller and returns one or two phases without a caller-
supplied phase count or phase guesses. Finite search cannot prove global
stability: every result retains
`globality_certificate="not_guaranteed"`, and search exhaustion, unsupported
trace or third-phase cases, provider failures, and indeterminate states fail
closed. HELD diagnostics report solver, numerical, and physical evidence on
independent `passed`, `failed`, or `not_adjudicated` axes. The accepted pure-
saturation route and authority are unchanged.

The retained scientific source is the Pereira Markdown in the permanent lab
at commit `13ce345b6dcc41d399bb2a4c7b9bedb18f74b45b`, blob
`dde7f02d4c93cce86804a8e6b62d37602990ac21`. The source is provenance only;
the clean package has no lab, migration, sibling-source, or provider-
implementation dependency.

The prior `neutral-two-phase-tp-flash-v1` campaign remains historical
`NON_ADMISSION` under the frozen `3*u_c` composition contract. Validation HEAD
`5a678beff38717478fd333c65e77f005cc2f6b15` and record SHA-256
`239c84788f75f8c66240c83e4f5874f112e1197dafad6273e1c8ec4efe994d24`
record 12/17 admitted rows, solved composition misses 002/009/010/011, row 012
without a package-accepted local state, and maximum material-balance error
1.11e-16. HELD does not erase or reclassify those model/data results.

The corrected HELD candidate is retained as one wheel built from implementation
commit `8318e755` with `SOURCE_DATE_EPOCH=1784573453` under
`artifacts/equilibrium-neutral-held-v1/8318e75/`. Its SHA-256 is
`8ecd70e0192b76b3a107629201c3e8bf34f2d945ca7c8192f824a0df7c9dde12`.
The earlier `549162a3` wheel with SHA-256
`be50837d73facbf0f2cb02cc3cfa7568f820b317b7b3984780ec1d992dbce76c`
remains immutable rejected-subject provenance and was not overwritten.
The candidate receipt is
`receipts/promotion/promotion-0002-neutral-held-v1-candidate.yaml` and remains
immutable evidence for that exact pre-validation subject. Permanent lab
approved the local candidate, while Validation commit
`93ff18541d2fe277a27671e4e6d12b6b009a58ed`, tree
`5aa2bc81941d1e807ba4e579231c4af9b7be15d7`, retained the installed-artifact
campaign as stable `NON_ADMISSION`: two one-phase results, thirteen
`scope_exceeded/third_candidate` results, three `search_exhausted` results,
and no evaluable two-phase composition comparisons. No promotion or authority
transfer has occurred.

Permanent-lab causal review accepted that evidence and requires a focused
binary controller redesign. Loose Stage-II cuts must remain separate from
provisional phase-candidate pairs; feasible pairs must be refined before
degenerate merging or genuine ambiguity is adjudicated; duplicate basins must
not stop the remaining 20 starts; and a complete nonprogressing pass must fail
closed promptly with exact reason `no_progress`. The design-only delta is
Section 13 of `docs/designs/2026-07-17-neutral-held-v1.md`. Runtime, tests,
artifacts, resources, tolerances, public types, binary scope, and globality are
unchanged. The user-authorized Perdomo HELD2 pivot defers this neutral redesign
as non-production provenance; it is no longer the active implementation gate.

`runtime_source_of_truth: true`

`accepted_capability: pure-component-saturation-v1`

`promotion_status: accepted`

`local_candidate: neutral-held-v1`

`local_candidate_authority_effect: none`

`local_candidate_review: approved-local-candidate-controller-redesign-required`

`validation_admission: NON_ADMISSION`

`local_candidate_promotion_status: not_requested`

The implemented design is `docs/designs/2026-07-17-neutral-held-v1.md`. It
replaced the unpromoted fixed-two-phase route and its duplicate public/result
surface without an alias. The complete controller still adds no accepted
capability or authority until independent review, installed-artifact
validation, ordered provider-tail and equilibrium promotion receipts, and
explicit user approval all complete.

`implemented_local_design: neutral-held-v1`

`implemented_local_design_authority_effect: none`

The executable plan is
`docs/plans/2026-07-17-neutral-held-v1-plan.md`. It freezes internal Stage I,
Stage II, Stage III, final public cutover, fixed-route deletion, isolated-wheel
proof, and review checkpoints. Tasks 1--7 remain frozen executed-v1 provenance;
Task 8 records the smallest post-validation controller correction, but the
Perdomo HELD2 pivot defers it. Its local status is
`controller_redesign_design_deferred`; the plan and
candidate receipt do not admit the provider mixture tail or create an accepted
equilibrium capability.

`implementation_plan_status: controller_redesign_design_deferred`
