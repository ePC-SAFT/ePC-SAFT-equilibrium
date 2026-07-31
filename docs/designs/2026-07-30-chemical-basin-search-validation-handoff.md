# Chemical basin-search validation handoff

This handoff binds the generic homogeneous chemical-equilibrium implementation
from [Equilibrium issue #96](https://github.com/ePC-SAFT/ePC-SAFT-equilibrium/issues/96)
to the immutable installed artifact used by downstream
[Validation issue #18](https://github.com/ePC-SAFT/ePC-SAFT-validation/issues/18)
and [Regression issue #16](https://github.com/ePC-SAFT/ePC-SAFT-regression/issues/16).
It does not transfer chemistry ownership, define an MEA campaign, or claim
global or phase stability.

## Immutable subject

| Item | Identity |
| --- | --- |
| Equilibrium implementation commit | `f433e4f706301c50e88c93572c3933e0bb702fef` |
| Equilibrium implementation tree | `9c8ade799876a01fec1460f56f71c7f93bca9f52` |
| Equilibrium wheel | `epcsaft_equilibrium-0.2.0.dev0-cp313-cp313-linux_x86_64.whl` |
| Equilibrium wheel SHA-256 | `adf241cfdc8a8e1aa186419440adf36f0405ce195f403a9c0461932ba345e198` |
| Installed-artifact `RECORD` SHA-256 | `235b8b8f111c98baff0273e2fb2f26be951554fa155c5da69b3c894640a7bcc1` |
| Provider wheel SHA-256 | `66b7ea8fb29e0a268b555cbdf401c3502517c088669a4157e8f64ab985b59ce9` |
| Ipopt / linear solver | `3.11.9` / MUMPS |

The wheel was built from the clean implementation commit and retained
read-only in the project artifact object store under its SHA-256. A fresh
Python 3.13 environment installed only the two wheels above plus declared
runtime dependencies. Its public ideal tracer returned one certified basin
from five ordered primary attempts. Its installed manufactured nonconvex
tracer retained two certified basins and selected the lower exact objective.

## Public receipt

`ChemicalEquilibriumDiagnostics.search` is a
`ChemicalEquilibriumSearch` containing:

- top-level outcome, continuation status, fixed primary budget, actual primary
  count, selected basin/value, and selection label;
- every primary and launched recovery attempt, with deterministic lineage,
  construction/retraction/domain status, solver/callback status, physical
  state, exact objective, residuals, KKT rank/conditioning, curvature class,
  terminal class, and duplicate-basin reference;
- every materially distinct certified basin and its canonical representative;
  and
- clipped nested budget projections at 1, 5, 13, and 25 starts.

The accepted selection label is
`lowest_observed_certified_local_value`. The existing globality field remains
`not_guaranteed`. Continuation is not used by this implementation and is
reported as `not_used`.

Regression evaluator v1 requires the native selected result to be accepted
before it emits an `OK` row. Rank-deficient, condition-limited, saddle,
inconclusive, boundary, domain-rejected, or exhausted searches therefore fail
the existing ABI gate without an ABI schema change.

## Frozen downstream sentinel

The application-owned source subject remains:

- MEA-Thermodynamics commit
  `269c954230b73bffe19d157137143a52d9c685f6`, tree
  `7d58f7b50b2d3e0682a862b92e8f1c998501cf80`;
- Hilliard `vle_obs_0137`: 313.15 K, 7326.7 Pa, unloaded MEA mass
  fraction 0.30, loading 0.466 mol CO2/mol MEA, and observed
  `pCO2 = 574 Pa`;
- Böttinger `cheq_canon_00194` at the same declared state, with observed
  `x_MEACOO- = 0.0502`;
- nine species and five source-bound reactions from
  `mea-nine-species-reaction-source-contract-v2`; and
- Provider regression-input parameter fingerprint
  `sha256:3773585e061b37643f5c7794e18424b83c86b82fa658983a0ee13fd8f1876fd6`
  and domain fingerprint
  `93510b66543e4e9e49c409a658b1bf7a01599ccd9ce3feef41bbab6b6eb668ab`.

The Equilibrium repository does not copy the application chemistry or source
rows. [Validation issue #18](https://github.com/ePC-SAFT/ePC-SAFT-validation/issues/18)
must consume those immutable application inputs and the two installed wheels
above through public package boundaries. It must retain the complete search
receipt even when the sentinel remains a saddle or no certified basin is
found.

## Interpretation boundary

A successful result certifies one observed strict local minimum of the declared
single homogeneous fixed-`T,P` problem. It does not certify the lowest
unobserved basin, phase stability, reactive phase equilibrium, global
thermodynamic stability, kinetics, predictive validity, or regression
readiness. Finite search failure never emits `infeasible_certified`; that
status remains reserved for a future independent feasibility certificate.
