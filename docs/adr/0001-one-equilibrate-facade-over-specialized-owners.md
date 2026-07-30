---
status: accepted
date: 2026-07-29
---

# Present one `equilibrate` facade over specialized equilibrium owners

The future public equilibrium surface will be one typed `equilibrate`
operation, but it will not force every physical problem through one numerical
algorithm. A homogeneous reacting problem delegates to the existing chemical
owner, a qualifying nonreactive phase-discovery problem delegates to its
existing HELD/HELD2 owner, and a reactive phase-discovery problem uses the
simultaneous GREPE controller. Delegation preserves each existing route's
qualification and rejection behavior; the facade does not admit a previously
unsupported nonreactive input. This preserves the proven limiting algorithms,
avoids sequential reaction-then-flash and copied controllers, and keeps each
result's scientific claim explicit.

During development, `tp_flash` and `chemical_equilibrium` remain frozen public
operations. The private `equilibrate` variants must reproduce those routes
through the same native owners before any public cutover is considered.
Removal of the older symbols is a later evidence-based decision, not part of
the initial GREPE implementation.
