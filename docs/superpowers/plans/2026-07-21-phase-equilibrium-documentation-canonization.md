# Phase-Equilibrium Documentation Canonization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Establish one current package-local phase-equilibrium authority map and one detailed Perdomo HELD2 design while correcting stale governing metadata without changing runtime behavior or capability authority.

**Architecture:** `docs/phase-equilibrium.md` is the compact navigation and authority owner. Existing pure-saturation and neutral-HELD designs remain detailed owners, while a new Perdomo design owns the implemented electrolyte formulation. Root package documents link to those owners and preserve receipts, Migration, Validation, papers, and lab M4 material in evidence/provenance roles.

**Tech Stack:** Markdown, YAML, Git, repository Python/static/test tooling

## Global Constraints

- Accepted capability remains exactly `pure-component-saturation-v1`.
- Perdomo HELD2 remains non-production, unpromoted, and `globality_certificate="not_guaranteed"`.
- No runtime source, tests, public symbols, native targets, artifacts, receipts, Provider contracts, scientific tolerances, resources, or seeds change.
- Do not copy EOS equations, Provider parameters, Validation case constants, old monorepo routes, archived private adapters, or M4 issue/queue history.
- Preserve the pure-saturation and neutral-HELD designs as their existing formulation owners.
- Perdomo and Ascani coordinate families remain explicitly non-equivalent.
- Keep unrelated `.idea` worktree changes untouched and out of every commit.

---

### Task 1: Add the package authority map and Perdomo design

**Files:**
- Create: `docs/phase-equilibrium.md`
- Create: `docs/designs/2026-07-21-perdomo-held2.md`

**Interfaces:**
- Consumes: Organization doctrine revision 3; existing pure/neutral designs; Perdomo and Pereira primary sources; live `held2` implementation/tests; Migration D-024 through D-026; Validation source ledgers.
- Produces: one package documentation entrypoint and one detailed current Perdomo formulation owner.

- [ ] **Step 1: Write the canonical authority map**

Create `docs/phase-equilibrium.md` with these exact top-level sections:

```markdown
# Phase-Equilibrium Documentation Authority

## Authority and capability status
## Public operations
## Formulation owners
## Shared package contract
## Evidence and claim axes
## Closed future formulations
## Historical and scientific provenance
```

The capability table must classify pure saturation as accepted, neutral HELD as frozen `NON_ADMISSION`, Perdomo HELD2 as non-production development, and Ascani/reactive formulations as closed.

- [ ] **Step 2: Check the authority map against current public exports and authority**

Run:

```bash
rg -n "saturation|tp_flash|pure-component-saturation-v1|NON_ADMISSION|not_guaranteed" \
  docs/phase-equilibrium.md README.md CONTEXT.md ARCHITECTURE.yaml \
  src/epcsaft_equilibrium/__init__.py
```

Expected: only `saturation` and `tp_flash` are described as current public operations; no electrolyte acceptance or globality claim appears.

- [ ] **Step 3: Write the Perdomo HELD2 design**

Create `docs/designs/2026-07-21-perdomo-held2.md` with these exact top-level sections:

```markdown
# Perdomo HELD2 Strong-Electrolyte Phase Equilibrium

## Status and authority
## Scope and nonclaims
## Scientific sources and notation
## Modified-mole coordinates
## Provider and derivative contract
## Homogeneous reference selection
## Stage I: stability
## Stage II: dual cutting and candidate discovery
## Stage III: direct total-free-energy refinement
## Certification and status axes
## Current evidence and next gate
## Provenance and deliberately excluded material
```

Include the eliminated-ion, modified-mole, modified-potential, and direct Stage-III equations from the approved spec. Describe the implemented `log(V)`, `q=log(eta)`, feasible-simplex, exact dual-pullback, reference-root, domain, complementarity, and fail-closed contracts. Do not embed case constants or artifact hashes as formulation doctrine.

- [ ] **Step 4: Run a scientific consistency scan**

Run:

```bash
rg -n "individual.ion|Galvani|electroneutral|modified potential|negative TPD|no_negative|root_completeness|globality|third derivative|A.*P.*V" \
  docs/designs/2026-07-21-perdomo-held2.md
```

Expected: individual-ion equality is conditional on an explicit Galvani convention; negative-witness and complete-no-negative asymmetry are explicit; direct `A+PV` Stage III and finite-search non-globality are explicit.

- [ ] **Step 5: Commit the canonical design owners**

```bash
git add -- docs/phase-equilibrium.md docs/designs/2026-07-21-perdomo-held2.md
git diff --cached --check
git commit -m "docs: canonize phase-equilibrium formulations"
```

### Task 2: Correct and link the governing package documents

**Files:**
- Modify: `AGENTS.md`
- Modify: `CONTEXT.md`
- Modify: `ARCHITECTURE.yaml`
- Modify: `README.md`

**Interfaces:**
- Consumes: canonical documents from Task 1.
- Produces: discoverable, revision-3-consistent package metadata without changing capability authority.

- [ ] **Step 1: Update `CONTEXT.md`**

Change doctrine revision 2 to 3. Add links to `docs/phase-equilibrium.md` and the Perdomo design. State that `runtime_source_of_truth` means implementation source, while accepted authority remains receipt-bound. Add the Perdomo primary-source identity, independent root-completeness semantics, and D-026 waiting boundary.

- [ ] **Step 2: Update `ARCHITECTURE.yaml`**

Change doctrine revision 2 to 3 and add:

```yaml
documentation:
  phase_equilibrium_authority: docs/phase-equilibrium.md
  perdomo_held2_design: docs/designs/2026-07-21-perdomo-held2.md

development_candidate:
  root_completeness: not_proven
  globality_certificate: not_guaranteed
  next_evidence_gate: source-complete-installed-two-liquid-stage-ii-iii-case
```

Integrate these fields into the existing structure rather than duplicating the `development_candidate` mapping. Preserve `authority_effect: none` and the existing accepted capability.

- [ ] **Step 3: Update `README.md` and `AGENTS.md`**

Add the canonical authority/design links to README and one compact discoverability pointer to AGENTS. State that public electrolyte dispatch is non-production and has not admitted electrolyte LLE. Explain that detected root accounting does not prove root completeness or globality.

- [ ] **Step 4: Parse and inspect the governing files**

Run:

```bash
python - <<'PY'
from pathlib import Path
import yaml

architecture = yaml.safe_load(Path("ARCHITECTURE.yaml").read_text())
assert architecture["authority"]["governance_doctrine_revision"] == 3
assert architecture["accepted_capability"]["name"] == "pure-component-saturation-v1"
assert architecture["development_candidate"]["authority_effect"] == "none"
assert architecture["development_candidate"]["globality_certificate"] == "not_guaranteed"
PY
```

Expected: exit 0.

- [ ] **Step 5: Commit the governing-file reconciliation**

```bash
git add -- AGENTS.md CONTEXT.md ARCHITECTURE.yaml README.md
git diff --cached --check
git commit -m "docs: reconcile equilibrium authority and status"
```

### Task 3: Verify canonization and minimum surface

**Files:**
- Verify only: all files changed by Tasks 1 and 2

**Interfaces:**
- Consumes: completed documentation changes.
- Produces: fresh schema, link, source-test, diff, cleanup, and status evidence.

- [ ] **Step 1: Check local Markdown and repository links**

Run a small read-only Python checker that extracts relative Markdown links from the six governing/design files, resolves them from each source file, and fails if any package-local relative target is missing. Explicit absolute provenance paths must be checked separately and labeled workstation-local.

- [ ] **Step 2: Search for stale governing claims**

```bash
rg -n "governance_doctrine_revision: 2|current_public_routes:|runtime_status: closed|stage_ii_status: not-executed|stage_iii_status: not-executed|packages/epcsaft-equilibrium" \
  AGENTS.md CONTEXT.md ARCHITECTURE.yaml README.md docs/phase-equilibrium.md \
  docs/designs/2026-07-21-perdomo-held2.md
```

Expected: no matches. Historical files outside this governing set may retain those terms.

- [ ] **Step 3: Run source and static verification**

Inspect `pyproject.toml` and repository commands, then run the current source suite and owned static checks. At minimum:

```bash
python -m pytest
python -m ruff check .
python -m ruff format --check .
python -m mypy --strict src tests
```

Use the repository's configured environment/runner if direct modules are unavailable. Expected: all commands exit 0.

- [ ] **Step 4: Review minimum surface and authority boundaries**

Confirm from `git diff --stat` and `git diff` that:

- only documentation and architecture metadata changed;
- no M4 archive, registry, Provider equation, case constant, or artifact was copied;
- there is one Perdomo design owner and one authority index;
- no public operation, result owner, tolerance, resource, or authority claim changed; and
- unrelated `.idea` changes remain unstaged and untouched.

- [ ] **Step 5: Run final repository checks**

```bash
git diff --check
bash "$HOME/.codex/hooks/codex-cleanup.sh" --repo-root .
git status --short --branch
git log -5 --oneline --decorate
```

Expected: documentation commits present; repository has no task-owned temporary artifacts; only pre-existing user-owned `.idea` changes may remain.

- [ ] **Step 6: Request review**

Send the exact final commit/tree, changed files, checks, remaining authority exclusions, and links to the canonical documents to the Migration/permanent reviewer when a review task is available. Do not push, promote, publish, create a receipt, or alter authority without a separate instruction.
