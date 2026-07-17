# Independent Review: Pure Saturation v1

Reviewer: `codex-independent-reviewer:/root/independent_review`

Implementation owner: `codex:/root`

Final implementation commit: `eff707fbf4fb011627e8089bc13aad99b9a21dce`

Final verdict: `APPROVE_FOR_MANAGER_REVIEW`

## Scope

The review covered the first local fixed-temperature pure-component saturation
slice for the provider-approved methane, ethane, and propane fingerprints. It
checked doctrine and repository standards, provider ownership, exact derivative
use, Ipopt callbacks, physical acceptance, installed-artifact isolation,
scientific evidence, diagnostics, and excluded capabilities.

## Initial findings and closure

The first pass required fixes before manager review. The implementation owner
closed every Critical and Important item:

- The provider gate now locates a non-editable installed distribution, checks
  the exact provider wheel hash before the isolated build, proves both packages
  import from the isolated install, and retains a wheel/link/symbol audit.
- The ABI proof now evaluates after deleting the Python model, rejects capsule
  name, ABI, table-size, result-size, returned-struct-size, and fingerprint
  mismatches, and keeps the capsule alive across the GIL-released solve.
- Solver evidence now retains bounds and every search/confirmation seed,
  termination, iteration count, constraint violation, and callback error.
  Mechanically unstable seeds are rejected before Ipopt.
- Pint is a direct dependency; dimensional input errors are normalized; native
  failures become structured `SaturationError` results.
- Directional checks now cover representative and boundary-near states for all
  three components and assert the zero objective and gradient.
- The architecture manifest records every numerical and scientific tolerance.
  The anchor record now includes source commits, source hashes, transformation,
  use basis, and tolerance rationale.
- Native, bytecode, and IDE ignore patterns were restored.

The reviewer treated duplicate scope declarations and the single native source
module as non-blocking for this first bounded baseline. Splitting or introducing
a new registry layer now would enlarge the surface without improving the
admitted capability.

## Final evidence

The reviewer independently recomputed the final hashes and confirmed:

- equilibrium wheel: `f5105830a7cb30da1f17277285a530b7cb5ee5ca9173d4b9ea5d7d4bd1401e0b`
- installed-wheel JUnit XML: `e0d27c4505e4f8c2650abd8350a7fb18ca1df470497f682ecf2a4c88efa8e4cc`
- scientific validation log: `e5419987c68bc91274d6e5a48b75e51abf880a06390f23cea1b5c9ec763b3c5d`
- binary audit: `f1ee92b1ec0323b7ddd88e0de3dd5b24b39bfdaab8659d74a201746da3ce0cac`
- extracted native extension: `ca89139293608c0dd983ab31bcff3883a2b94a5f898b62f34adf609756ccea3a`

The XML contains 15 tests with zero failures or errors. The binary audit shows
only the equilibrium package in the wheel, one native target built from the two
equilibrium translation units, a direct Ipopt link, no provider implementation
material, and no provider/EOS symbols. `git diff --check 27c89bd..eff707f` is
clean.

## Authority boundary

This approval is for manager review, not promotion. The lab remains the
transitional runtime authority. No globality, phase-discovery, mixture,
association, electrolyte, flash, regression, release, or workflow claim was
reviewed or accepted.
