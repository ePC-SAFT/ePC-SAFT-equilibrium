# Domain Docs

The engineering skills use this repository's single-context domain
documentation.

## Before exploring

- Read the root `CONTEXT.md`.
- Read ADRs under `docs/adr/` that affect the work.
- If either location is absent, proceed without creating speculative
  documentation.

## Layout

```text
/
├── CONTEXT.md
├── docs/
│   └── adr/
└── src/
```

## Vocabulary

Use the domain terms established by `CONTEXT.md`. If a required concept is
missing, first determine whether the proposed term is unnecessary or exposes a
real modeling gap.

## ADR conflicts

Surface any proposed work that contradicts an existing ADR. Do not silently
override the recorded decision.

<!-- project-truss:setup:start -->
## Project Truss domain contract

Layout: `single-context`.

Use one root `CONTEXT.md` and root `docs/adr/`.

Read the nearest applicable context before shaping or implementation. Keep durable terminology, invariants, ownership, and accepted decisions there; do not mirror lifecycle state.
<!-- project-truss:setup:end -->
