# ePC-SAFT Equilibrium

- Native-extension package owning equilibrium formulations, phase discovery, deterministic root/finite global exploration, Ipopt/LP execution, certification, diagnostics, and equilibrium result contracts.
- Provider is an installed non-editable `epcsaft==0.1.0.dev0` distribution and remains sole EOS/derivative/domain owner; never copy Provider equations into this repo.
- Canonical formulation/capability index: `docs/phase-equilibrium.md`; HELD2 design and execution authority: `docs/designs/2026-07-21-perdomo-held2.md` and `docs/plans/2026-07-21-perdomo-held2-solver-strategy.md`.
- Repository policy: `AGENTS.md`; organization doctrine: sibling `../ePC-SAFT-organization/GOVERNANCE.md`.
- Current public accepted capability is only pure-component saturation v1. HELD/HELD2 work remains private/non-production until separately promoted.
- C++ implementation is under `cpp/src`; Python package/bindings under `src/epcsaft_equilibrium`; scientific contracts under `tests`.
- Read `mem:tech_stack` for build dependencies and `mem:task_completion` before declaring implementation complete.