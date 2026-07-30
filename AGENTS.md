# Equilibrium Repository Agent Contract

Repository Profile: scientific-computing

Accepted capability: pure-component-saturation-v1.

Before work, read doctrine revision 5 and the machine-readable access contract
from `../ePC-SAFT-management/GOVERNANCE.md` and
`../ePC-SAFT-management/manifests/management.json`. During the local checkout
transition, resolve both files from `../ePC-SAFT-organization/`. Their published home is
`https://github.com/ePC-SAFT/ePC-SAFT-management/blob/main/GOVERNANCE.md`.
Local policy may narrow that doctrine but must not contradict it.

This repository is limited to equilibrium routes, phase discovery,
deterministic root and finite global exploration, stage-owned LPs and local
NLPs, Ipopt, certification, and equilibrium results. It consumes a compatible
installed EOS and does not own EOS internals or Ceres.

Canonical package-local formulation ownership and capability boundaries are
indexed in `docs/phase-equilibrium.md`. Historical permanent-lab M4 documents
are scientific and transition provenance, not current package authority.

## Agent access and Validation work

- The Equilibrium agent remains rooted in this repository and is the sole owner
  of equilibrium formulations, phase discovery, equilibrium-solver execution,
  certification, diagnostics, and equilibrium result contracts.
- When Management assigns a bounded Equilibrium campaign, this same agent uses
  a target-owned Validation worktree created on demand under Validation's
  `.worktrees/` directory for a task-specific branch and path. Its task home
  does not move to Validation.
- Validation work must use exact immutable installed EOS and Equilibrium
  public artifacts. It must not import either source checkout, use private
  adapters, copy EOS or controller logic, or tune package behavior from the
  campaign.
- The Equilibrium agent owns the campaign execution and evidence commit, then
  returns to Equilibrium work. It does not create or wait on a resident
  Validation task, accept its own promotion, or turn finite search, local KKT,
  or a completed campaign into a globality proof.
- Management serializes Validation writers and dispatches any required distinct
  review after the exact stable subject is available.

Accepted receipt `promotion-0018-equilibrium-pure-saturation-v1` is
the sole authority record for the fixed-temperature methane, ethane, and
propane saturation capability. The capability supplies one local boundary
solve with no phase-discovery or
global-stability certificate. Do not broaden it to mixtures, association,
electrolytes, flash, critical continuation, regression, workflows, or releases
without a separately admitted slice.

## Agent skills

### Issue tracker

Issues and PRDs are tracked in this repository's GitHub Issues using `gh`.
See `docs/agents/issue-tracker.md`.

### Triage labels

Use the canonical Matt triage labels without aliases.
See `docs/agents/triage-labels.md`.

### Domain docs

Use the single-context layout: root `CONTEXT.md` plus repository-wide ADRs
under `docs/adr/`. See `docs/agents/domain.md`.
