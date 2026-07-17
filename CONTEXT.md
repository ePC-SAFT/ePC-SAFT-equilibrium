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

The repository now also contains the locally implemented
`neutral-two-phase-tp-flash-v1` candidate. It is limited to the reviewed
Gross--Sadowski methane/ethane fingerprint and the audited May et al. (2015),
Table 5 rectangular input domain. It fixes the phase count at two and carries
only solver, numerical-confirmation, and local-physical evidence. It has no
promotion receipt and changes no runtime authority.

The audited row `may2015-ch4-c2h6-001` is the package source anchor. Permanent
lab review status is `APPROVE_LOCAL_CANDIDATE`. Validation repository final
evidence HEAD `5a678beff38717478fd333c65e77f005cc2f6b15` and record SHA-256
`239c84788f75f8c66240c83e4f5874f112e1197dafad6273e1c8ec4efe994d24`
returned `NON_ADMISSION` under the frozen `3*u_c` composition contract: 12/17
admitted, rows 002/009/010/011 solved but missed composition allowances, row
012 had no package-accepted local state; maximum material-balance error
1.11e-16. Package local-candidate approval remains distinct from validation
admission. Authority effect remains none and promotion is blocked.

`runtime_source_of_truth: true`

`accepted_capability: pure-component-saturation-v1`

`promotion_status: accepted`

`local_candidate: neutral-two-phase-tp-flash-v1`

`local_candidate_authority_effect: none`

`local_candidate_review: APPROVE_LOCAL_CANDIDATE`

`validation_admission: NON_ADMISSION`

`local_candidate_promotion_status: blocked`

The approved next design is `docs/designs/2026-07-17-neutral-held-v1.md`.
It replaces a separately promoted fixed-two-phase route with one bounded
Pereira HELD Stage-I/II/III controller and final `tp_flash` operation. The
design uses this repository's existing native module, Ipopt installation,
`ProviderContext`, and the installed provider Helmholtz value/gradient/Hessian
callback. It adds no current capability or authority until implementation,
independent review, installed-artifact validation, and a promotion receipt all
complete.

`next_approved_design: neutral-held-v1`

`next_design_authority_effect: none`
