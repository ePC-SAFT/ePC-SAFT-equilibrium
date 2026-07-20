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

The repository now also contains the complete local `neutral-held-v1`
candidate at implementation commit
`549162a3a9cfd6f02894f8189c624ba1aa2139fb`, tree
`fe96ea3469405242cb4dbdb652e105dab5003f65`. It is limited to the reviewed
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
closed. The accepted pure-saturation route and authority are unchanged.

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

The complete HELD candidate is retained as one wheel built from implementation
commit `549162a3` with `SOURCE_DATE_EPOCH=1784567529` under
`artifacts/equilibrium-neutral-held-v1/549162a/`. Its SHA-256 is
`be50837d73facbf0f2cb02cc3cfa7568f820b317b7b3984780ec1d992dbce76c`.
The candidate receipt is
`receipts/promotion/promotion-0002-neutral-held-v1-candidate.yaml`. Permanent-
lab complete-candidate review and Validation's separate HELD campaign remain
pending; no promotion or authority transfer has occurred.

`runtime_source_of_truth: true`

`accepted_capability: pure-component-saturation-v1`

`promotion_status: accepted`

`local_candidate: neutral-held-v1`

`local_candidate_authority_effect: none`

`local_candidate_review: pending-complete-candidate-review`

`validation_admission: not_run`

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
proof, and review checkpoints. Its local implementation status is
`complete_candidate_pending_review`; the plan and candidate receipt do not
admit the provider mixture tail or create an accepted equilibrium capability.

`implementation_plan_status: complete_candidate_pending_review`
