# Independent Review: Pure Saturation v1

Reviewer: `codex-independent-reviewer:/root/independent_review`

Implementation owner: `codex:/root`

Baseline implementation commit: `eff707fbf4fb011627e8089bc13aad99b9a21dce`

Correction implementation commit: `faba1f9ab9343db333e404dbd2bb1cb9345d798c`

Compatibility-test commit: `95c55705581808ccf38c7bb57891d16ac4b775aa`

Final verdict: `APPROVE_CORRECTIONS_FOR_MANAGER_REVIEW`

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
- Native and bytecode ignore patterns were restored. The user later directed
  the owner to track the existing PyCharm and CLion `.idea/` metadata while
  retaining ignores for generated builds, caches, binaries, and environments.

The reviewer treated duplicate scope declarations and the single native source
module as non-blocking for this first bounded baseline. Splitting or introducing
a new registry layer now would enlarge the surface without improving the
admitted capability.

## Correction verification

The focused correction review found that:

- Solver acceptance requires a converged Ipopt status, a finite constraint
  violation no greater than `1.0e-10`, and an empty callback error before the
  physical gates run. The accepted-case test checks the finite top-level
  violation and every recorded callback error.
- The consumer decodes provider error and parameter-fingerprint arrays through
  one length-bounded decoder. A synthetic malformed capsule test rejects a
  missing NUL in both paths, including the module fingerprint seam.
- Native-exception diagnostics report `exact_derivatives = false`. The anchor
  test and scientific validator require exactly methane, ethane, and propane.
- CMake provider discovery uses explicit failures instead of Python assertions.
  The isolated wheel build passed with `PYTHONOPTIMIZE=2`.

The existing accepted-case test proves distinct phases and inactive bounds.
A direct bound-failure injection would require a new native or public test seam,
so the owner retained the implementation rejection gate without expanding the
runtime surface.

## Final evidence

The reviewer independently recomputed the final hashes and confirmed:

- equilibrium wheel: `a7f119253dd3fd41b609105916b8184cbe24b64e67e37873426d22e43501f490`
- installed-wheel JUnit XML: `d4d58183487a2d8e01ef48473ea41a0582271dfaac2ede73f88ead5493f55ba7`
- scientific validation log: `e5419987c68bc91274d6e5a48b75e51abf880a06390f23cea1b5c9ec763b3c5d`
- binary audit: `b498c273c6377f2ea933ec754a098af8e14216656a0b55fe82491374f6dc546a`
- extracted native extension: `673367e6bfc50f249e998ea8110148a461ad9594cde4c7e34e8126b35663f4f9`

The XML contains 15 tests with zero failures or errors. The binary audit shows
only the equilibrium package in the wheel, one native target built from the two
equilibrium translation units, a direct Ipopt link, no provider implementation
material, and no provider/EOS symbols. `git diff --check 5f2b85f..faba1f9` is
clean.

## Final provider compatibility replay

The implementation owner replayed the unchanged equilibrium wheel against
provider commit `4b10cb899c94687cae734980285badb224dc95e6`. The provider extends the v1 SDK
table from the 40-byte fixed-state prefix to 56 bytes. The consumer accepted
the larger table, retained ABI version 1 and the exact 376-byte fixed-state
result, required the model context and fixed-state evaluator, and did not read
the parameterized tail.

The first replay exposed an exact-size assertion in the test contract. Commit
`95c55705581808ccf38c7bb57891d16ac4b775aa` changed that assertion to require a
table size of at least the fixed-state v1 prefix. It kept the result-size and
required-prefix checks exact. The corrected isolated replay passed all 15
tests. The scientific validator reproduced the original byte-identical log,
and the binary audit found no provider implementation material or EOS symbols.

Replay evidence:

- provider wheel: `f92f79c8d6f614660e5c201b7061c9b02b5cd1a25a4ed8c8fee0b59adaabf2bf`
- provider test receipt: `07447721abaca946c6e9221e7d49e431e13fcb8e6867944f67b6ba8337901480`
- installed-wheel JUnit XML: `f9241eace37708bd81bae8df95f2b38ded3bee6d372f5016194fecfc248e5916`
- scientific validation log: `e5419987c68bc91274d6e5a48b75e51abf880a06390f23cea1b5c9ec763b3c5d`
- binary audit: `a916f6f95ff5157b87c242bc7f41d3aa6b942db42a352cb6c8a467bbed7d79ae`
- extracted native extension: `673367e6bfc50f249e998ea8110148a461ad9594cde4c7e34e8126b35663f4f9`

This replay changes no candidate capability or authority decision.

## Authority boundary

This approval is for manager review, not promotion. The lab remains the
transitional runtime authority. No globality, phase-discovery, mixture,
association, electrolyte, flash, regression, release, or workflow claim was
reviewed or accepted.
