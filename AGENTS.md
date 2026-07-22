# Equilibrium Repository Agent Contract

Repository Profile: scientific-computing

Accepted capability: pure-component-saturation-v1.

Before work, read doctrine revision 4 at
`../ePC-SAFT-organization/GOVERNANCE.md`. Its future published home is
`https://github.com/ePC-SAFT/.github/blob/main/GOVERNANCE.md`. Local policy may
narrow that doctrine but must not contradict it.

This repository is limited to equilibrium routes, phase discovery,
deterministic root and finite global exploration, stage-owned LPs and local
NLPs, Ipopt, certification, and equilibrium results. It consumes a compatible
installed provider and does not own provider internals or Ceres.

Canonical package-local formulation ownership and capability boundaries are
indexed in `docs/phase-equilibrium.md`. Historical permanent-lab M4 documents
are scientific and migration provenance, not current package authority.

## Agent home and Validation work

- The Equilibrium agent remains rooted in this repository and is the sole owner
  of equilibrium formulations, phase discovery, equilibrium-solver execution,
  certification, diagnostics, and equilibrium result contracts.
- When Migration assigns a bounded Equilibrium campaign, this same agent may use
  a clean checkout or worktree of the sibling Validation repository to author
  and execute it. Its task home does not move to Validation.
- Validation work must use exact immutable installed Provider and Equilibrium
  public artifacts. It must not import either source checkout, use private
  adapters, copy EOS or controller logic, or tune package behavior from the
  campaign.
- The Equilibrium agent owns the campaign execution and evidence commit, then
  returns to Equilibrium work. It does not create or wait on a resident
  Validation task, accept its own promotion, or turn finite search, local KKT,
  or a completed campaign into a globality proof.
- Migration serializes Validation writers and dispatches any required distinct
  review after the exact stable subject is available.

Accepted migration receipt `promotion-0018-equilibrium-pure-saturation-v1` is
the sole authority record for the fixed-temperature methane, ethane, and
propane saturation capability. The capability supplies one local boundary
solve with no phase-discovery or
global-stability certificate. Do not broaden it to mixtures, association,
electrolytes, flash, critical continuation, regression, workflows, or releases
without a separately admitted slice.
