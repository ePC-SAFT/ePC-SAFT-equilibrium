# HELD2.0 Validated Working Baseline

## Purpose and status

This document records the first local HELD2.0 implementation state known to
complete the Perdomo electrolyte cases and the now-retired Figure-2 electrolyte
case in reasonable time. It is historical evidence retained for review,
cleanup, and performance provenance. It is not a capability admission, a claim
that the finite search is globally complete, or authority to tune the Provider
EOS or parameter bundles.

The Figure-2 consumer and its numerical gate were retired under Equilibrium
issue #50 so that the EOS case-study bundle can be removed. The records below
are not current validation requirements.

The normative scientific workflow remains
[`2026-07-24-held2-paper-algorithm.md`](2026-07-24-held2-paper-algorithm.md).
At the time of this checkpoint, the implementation was evaluated against two
independent numerical gates; those records are retained as historical
provenance only.
The relationship between the paper's necessary conditions and the package
evidence is recorded in
[`2026-07-28-held2-necessary-condition-map.md`](2026-07-28-held2-necessary-condition-map.md).

## Reproducible source identity

- Repository branch: `codex/held2-step5-best-qualified`
- Validated integration base: `8abbe02f0d79cad05202f8998c119a5e3083353f`
  (Provider frontend 0.2, PR #76)
- Tracked production/test working-diff SHA-256:
  `f24031fa81a6612d2befa5ee6709190822893fd1db31ada49df24c2d64a50735`
- Retired Figure-2 numerical-gate SHA-256:
  `25bfc55842052134d62a0c8c9a3448451231e734ee460acbc854f787913825fe`

The tracked-diff digest covers `cpp/src`, `cpp/tests`, and
`tests/test_perdomo_held2_trace.py` relative to the validated integration base.
The separate digest covers the former Figure-2 numerical gate.

### Merged recovery anchor

Commit `409d7464d8be98cd3e9d58b1458fa29111e5f8ff` is the immutable
pre-GREPE recovery anchor. Compare or recover through a new branch or worktree,
never by overwriting an active checkout. Immediately before the GREPE design
work, that clean anchor passed:

- native `held2-manufactured`: passed in 0.12 s; and
- the unchanged Perdomo plus Figure-2 public gates: `6 passed in 27.62s`.

Together with the fixed values below, these are the behavior freeze. GREPE
work may add equivalence evidence and design documentation, but may not change
the protected results, tolerances, budgets, transitions, or evidence meanings.

## Working workflow invariants

There is one production controller, `run_held2_algorithm`, implementing
Steps 1–10. Both the native diagnostic and the Python API call that controller,
and all thermodynamic values and derivatives come through the installed
Provider public SDK.

The working state depends on these invariants:

1. Step 2 retains a negative-TPD witness only for the multivariate search in
   which it supplies information not already present in Appendix C's complete
   one-dimensional interval search.
2. Candidate identity is stable. Duplicate removal preserves the surviving
   candidate indices; it does not reconstruct identity from nearest
   composition.
3. Step 5 uses named tolerances and treats numerically equivalent terminal
   points consistently.
4. Step 8 receives the actual feasible HiGHS phase amounts and compositions
   as its LP-to-NLP start. It does not replace them with equal phase fractions
   or raw candidate centers.
5. A Step-8 phase touching its local neighborhood is recentered conditionally;
   neighborhood radii are not enlarged globally.
6. Step 8 retires at most one KKT-certified inactive candidate and then
   re-solves the full Problem (67) on the retained neighborhoods. It does not
   freeze compositions and repair only phase fractions.
7. Step 8 may retry cold after a failed warm continuation, with both solves
   included in the work accounting.
8. A non-converged Ipopt return is never accepted even if Ipopt supplied a
   finite iterate. Because Step 8 is an acceleration, the failed evidence is
   retained and the controller returns to finite Stage II search. Malformed or
   uncertified feasibility evidence remains terminal and indeterminate. A
   failed reduced-active-set solve or balance reconstruction is likewise
   rejected acceleration evidence and returns to Stage II; it is not a phase
   result. Candidate retirements certified by preceding converged KKT audits
   remain valid acceleration state: the next attempt starts from the terminal
   retained IDs and adds only a genuinely unattempted newest eligible cut.
9. Steps 5 and 8 use Ipopt gradient-based NLP scaling. Physical phase amounts
   and modified compositions remain linear optimization variables; small,
   valid compositions are neither clipped nor retired merely because of their
   magnitude. Step 10 alone uses its documented logarithmic trace coordinate.
10. Independently certified Step-8 phases remain available as Step-9/10 cuts.
    Search exhaustion remains
    `globality_certificate = "not_guaranteed"`.

These are observable scientific and state-machine constraints, not a mandate
to retain every current helper, diagnostic field, branch, or allocation.

## Numerical gate A: Perdomo

`tests/test_perdomo_held2_trace.py` supplies the fast regression matrix:

- Table 3 NaCl/water at 298.15 K, 2508 Pa, and 5.6 mol NaCl/kg water returns
  `one_phase_no_negative_witness_detected`.
- Figure 1a at 298.15 K and 0.005 mol NaCl/kg water returns a two-phase
  equilibrium both below the boundary at 3181.454397 Pa and at the digitized
  source boundary of 3213.5903 Pa.

The public-result gates compare total reduced free energy, phase fractions,
phase compositions, and phase volumes. Their fixed tolerances are:

- free energy and feed: absolute/relative `1e-8`;
- phase fractions: absolute `5e-8`, relative `1e-8`;
- phase volumes: absolute `5e-8`, relative `1e-8`;
- phase compositions: absolute `5e-7`, relative `1e-8`.

The observer parity test additionally requires tracing to leave the numerical
result unchanged.

## Retired evidence record: Figure 2 case study

The former standalone public-API check supplied one deliberately bounded
case-study record:

- T = 293.15 K
- P = 100000 Pa
- components = water, ethanol, isobutanol, sodium-cation, chloride-anion
- parameter bundle = retired Figure-2 case-study bundle@1
- parameter fingerprint =
  `sha256:b43fac77754d9d5cca8b3db2cbe709892a786d97b756084b167ce126ab4c3007`
- feed =
  `(0.7005034356224062, 0.03286847384962303, 0.21275572005079824,
  0.026936185238586256, 0.026936185238586256)`

The accepted result has two phases, reduced free energy
`-8.156553792022185`, and phase fractions approximately
`(0.5928509426897598, 0.4071490573102402)`. The test records the complete
phase compositions and volumes rather than weakening the gate to organic-side
agreement alone. It also requires material balance, pressure stationarity,
KKT stationarity, passed numerical/physical diagnostics, no more than 16 upper
solves, and `globality_certificate = "not_guaranteed"`.

The checkpoint native execution completed in 24.69 s. Its dominant costs were
Step 2 (13.23 s, 6620 Provider evaluations) and Step 3 (3.23 s); Step 5 took
4.13 s and Step 8 took 4.05 s. This evidence does not support attributing the
remaining runtime primarily to Ipopt scaling.

### Step-2 joint-search performance replacement

The paper-aligned Step-2 redesign replaced the nested pressure envelope at
every DIRECT-L composition with a joint modified-composition/log-volume
search. The complete pressure envelope remains at the feed, every negative
state is independently re-evaluated, and only promising downstream seeds
receive a bounded local pressure-root refinement.

For the exact gate above, a controlled pre-change receipt recorded:

- wall/CPU: `14.328897813 s` / `14.323798 s`;
- Provider evaluations: `6620`;
- returned witness TPD: `-0.010218719562037012`;
- returned modified fractions:
  `(0.05555555573333337, 0.49999999986666654,
  0.0555555556, 0.38888888880000005)`; and
- returned volume: `4.851779848875743e-05 m3/mol`.

The final joint-search receipt recorded:

- wall/CPU: `5.09454932 s` / `5.093721 s`;
- Provider evaluations: `2819`;
- raw minimum joint TPD: `-0.12934921110306705`; and
- the identical downstream witness TPD, modified fractions, and volume listed
  above.

This is a 2.81-fold Step-2 wall-time improvement and a 57.4% reduction in
Provider calls. The final accepted free energy remained
`-8.156553792022187`, with phase fractions
`(0.5928510002832944, 0.4071489997167056)`, inside the unchanged public gate.
The raw minimum is deliberately distinct from the pressure-stationary witness
prepared for later steps; a more-negative off-pressure state is already
sufficient to establish instability.

After the lifecycle cleanup, a fresh receipt recorded `5.254560755 s` in
Step 2 and the same `2819` Provider evaluations, raw minimum, downstream
witness, final free energy, and phase fractions. The cleanup therefore changed
neither the deterministic work count nor the accepted numerical state.

## Validation evidence

At this checkpoint:

- the native manufactured Steps 1–10 workflow passed in about 0.1 s;
- the five focused Perdomo tests plus the former Figure-2 evidence check passed
  (`6 passed in 42.36s`);
- the full Python suite passed before the final fail-closed Step-8 adjustment
  (`88 passed in 63.21s`);
- the exact focused native, Perdomo, and former Figure-2 paths passed again
  after that
  adjustment;
- after resolving persistent-state identity from Appendix C, the native
  manufactured workflow passed in 0.11 s, the five Perdomo checks plus the
  former Figure-2 evidence check passed (`6 passed in 44.32s`), and the full
  suite passed (`88 passed in 105.99s`);
- after the code-surface cleanup, the native manufactured workflow passed in
  0.11 s and the full Python suite passed (`88 passed in 30.93s`);
- after integration with Provider frontend 0.2, the native manufactured
  workflow passed in 0.09 s and the complete Python suite, including all
  native/Python diagnostic parity cases, passed (`131 passed in 40.01s`);
- Ruff, mypy, native/Python diagnostics, and repository cleanup checks passed.

The eight-row exploratory Figure-2 campaign remains historical
`NON_ADMISSION`/`SOLVER_LIMIT` receipt. This single accepted gate does not
retroactively change that campaign or establish broad Figure-2 coverage.

## Change rule

Cleanup and optimization may remove incidental machinery only when:

1. the native manufactured workflow passes;
2. the complete Perdomo gate passes without widening its tolerances;
3. the former Figure-2 evidence record remains unchanged;
4. fail-closed and independent-evidence semantics remain intact; and
5. a measured runtime improvement is not purchased by reducing the declared
   finite search or weakening certification.

Changes to tolerances, search coverage, EOS calls, parameter data, or
certificate meanings are scientific changes and require separate evidence.

## Full-workflow review boundary

The post-checkpoint review classified the following behavior as essential to
the successful path rather than cleanup debris:

- append-only \(\mathcal M\) membership, explicit attempted/terminal Step-8
  candidate identities and neighborhood variables, and exact reuse when the
  ordered effective Problem-(67) state is unchanged;
- active-set retirement followed by a complete reduced Problem-(67) solve;
- conditional neighborhood recentering and one cold retry;
- Step-8 phase feedback after a Step-9/10 return to Stage II; and
- linear small-composition variables with logarithmic refinement confined to
  Step 10.

The stabilization cleanup may remove duplicate evaluation and bookkeeping
work around those decisions, but not the decisions themselves.

The same review identified two remaining scientific-contract improvements
that are deliberately not folded into behavior-preserving cleanup:

1. Step-5 local termination currently relies on the Ipopt terminal plus a
   finite-domain check; the normative independent KKT/dual certificate needs
   a focused implementation and manufactured failure tests.
2. A Step-8 HiGHS infeasible status is not yet accompanied by the documented
   independently validated Farkas evidence.

Those items must be resolved as coherent scientific changes, one at a time,
under both numerical gates. They are not justification for changing iteration
limits, clipping trace components, widening neighborhoods, or adding another
recovery path.

The identity question is resolved directly from the paper rather than by
adding data to Step 3:

- Appendix C stores each persistent member as
  \((V,\bar{\mathbf{x}}^{(EC)})\), so \(\mathcal M\) numerical identity uses
  modified composition and molar volume. Initial members do not require a
  packing-fraction cache.
- Step 6 owns the distinctness test in Eq. (66), where Provider packing
  fraction and modified composition are the two independent axes.
- Step 8 removes phases that have the same composition and volume. The
  implementation compares physical composition and log-volume difference,
  with log-volume serving only as a dimensionless relative-volume measure.

The associated identity tolerances are named implementation policy. They may
be adjusted when convergence evidence shows that the numerical gate is
impractical, provided the change is exercised by both public numerical gates
and does not merge materially distinct states or convert indeterminate
evidence into acceptance.
