---
status: accepted
date: 2026-08-02
---

# Admit a typed Provider-model continuation input

## Decision

`chemical_equilibrium` may accept the optional typed
`ProviderModelContinuation` input defined by the Provider-model continuation
design. This narrowly supersedes ADR 0001's public-freeze sentence for this one
recovery input. It does not change the operation's owner, result claim, strict
certification gates, or the separate future `equilibrate` facade.

The input names a complete installed initial Provider endpoint, not an optimizer
seed or backend selector. Endpoint-bound reaction references are mandatory when
one source-standard contract cannot be transformed independently.

## Consequences

Downstream applications can request a reproducible model-family recovery path
without package-owned chemistry knowledge. A continuation terminal remains a
typed failure unless the exact target endpoint independently certifies. The
package public surface grows by one immutable input record and one optional
keyword; no compatibility alias or fallback route is introduced.
