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
