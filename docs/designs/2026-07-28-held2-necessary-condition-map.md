# HELD2.0 Necessary-Condition Map

## Purpose and authority

This map connects three layers that must not be conflated:

1. the thermodynamic conditions a valid electrolyte phase equilibrium must
   satisfy;
2. the Perdomo HELD2.0 search and convergence conditions; and
3. the independent evidence required by this package before it reports an
   accepted result.

The normative workflow is
[`2026-07-24-held2-paper-algorithm.md`](2026-07-24-held2-paper-algorithm.md).
This companion document explains how that workflow is realized and where the
package deliberately adds evidence or numerical policy. It is not a
capability admission or a globality proof.

“Paper” below means Perdomo et al., *Computers & Chemical Engineering* 194
(2025) 108977, especially Algorithm 1, Eqs. (61)–(69), and Appendix C.
“Package” means the one native Equilibrium controller consuming the installed
Provider public SDK.

## Physical equilibrium conditions

| Necessary condition | Paper representation | Package evidence | Failure meaning |
|---|---|---|---|
| Valid phase state | Positive molar volume and admissible modified composition | Provider-domain state, finite thermodynamic values, and complete Step-1 polytope checks | Indeterminate or invalid input; never a phase |
| Composition normalization | Eliminated-ion modified-mole coordinate transformation | Reconstructed physical fractions sum to one | Physical certificate fails |
| Phase electroneutrality | Charge is eliminated by the modified coordinates | Reconstructed, scaled charge residual for every phase | Physical certificate fails |
| Material balance | Phase amounts and modified compositions in Problem (67) | Modified-coordinate and explicit physical-component balances | Stage III is not accepted |
| Mechanical equilibrium | Common specified \(T,P\); volume is a phase variable | Independent pressure-stationarity residual for every phase | Local state or Stage III is not accepted |
| Chemical/electrochemical equilibrium | Equality of modified potentials in Eq. (69) | Modified-potential residual plus original-coordinate KKT pullback | Stage III is not accepted |
| Local stationarity | Problems (65) and (67) terminate locally | Primal, stationarity, dual-sign, and complementarity audits owned by the relevant step | Solver success alone is insufficient |
| Stability against a trial phase | Tangent-plane search and lower/upper construction | Negative witness is independently evaluated; finite search retains `not_guaranteed` globality | No witness is an operational one-phase result, not proof of global stability |
| Equilibrium phase set | Total-Gibbs Problem (67) and Eqs. (68)–(69) | Certified nonlinear result, phase identity, balances, pressure, KKT, and paper convergence gates | Any missing evidence fails closed |

## Step-by-step condition map

| Step | Perdomo HELD2.0 condition | Package realization and evidence | Deliberate difference or policy |
|---|---|---|---|
| 1 | Define the \(C-2\) modified-composition domain and volume bounds, Eq. (61) | Exact eliminated-ion chart, complete polytope, immutable Provider pressure envelope | Finite positive search floors regularize exploration; they are not physical composition minima |
| 2 | Minimize the tangent-plane distance jointly over admissible modified composition and volume, Problem (63) | Complete feed pressure-root enumeration selects the reference; DIRECT-L searches \(C-2\) modified-composition coordinates plus normalized log volume; every strict-negative witness is independently re-evaluated | Pressure stationarity is not required to prove negative TPD. A promising witness may be locally refined to a stable pressure root only for downstream initialization. The finite search cannot certify global absence, so every result retains `globality_certificate="not_guaranteed"` |
| 3 | Initialize the upper problem and persistent \(\mathcal M\) using Appendix C | Reproduce the Appendix-C endpoints and feed state; store immutable IDs with \((V,\bar{\mathbf{x}}^{(EC)})\) | A Step-2 witness is additionally retained only when it adds multivariate information not already supplied by the complete one-dimensional construction |
| 4 | Solve upper LP (64) | HiGHS solve with independent primal, dual, and complementarity audit | Optimizer status is evidence input, not the certificate |
| 5 | Solve nonconvex lower problem (65) and stop when the best qualified lower value is no greater than the upper bound | Deterministic Ipopt multistart; retain Ipopt bound and polytope multipliers; independently lift the audited chart variables back to physical composition and audit that reconstruction, physical-volume KKT, the terminal composition's Provider volume bounds, normalization/charge gauge multipliers, pressure closure, primal/dual feasibility, complementarity, chart pullback, and active-set rank. Bind the terminal to the exact Step-1 coordinate snapshot and same-major Step-4 certificate, upper-solve identity, immutable cut snapshot, recomputed active-cut IDs, upper bound, and multipliers before it may enter \(\mathcal M\) | Numerically equivalent polytope and variable bounds are canonicalized once, with their dual force preserved. A qualifying terminal may receive bounded exact-derivative composition/pressure polishing before the independent audit; diagnostics keep raw solver variables separate from the audited terminal. Missing, stale, nonfinite, rank-deficient, or inconsistent evidence is indeterminate; packing fraction is not persistent-\(\mathcal M\) identity |
| 6 | Search all of \(\mathcal M\) using the gap, gradient, and distinctness conditions in Eq. (66) | Re-evaluate every eligible member; use actual Provider packing fraction and modified composition as the two distinctness axes | Packing fraction belongs here only. It is never replaced by log-volume |
| 7 | Advance the HELD2 major iteration | Persistent major counter and the next unused deterministic Step-5 start ordinal | Resource caps are declared finite-work policy, not paper convergence conditions |
| 8 | Solve total-Gibbs Problem (67); remove numerical duplicate phases with the same composition and volume | Exact perspective-form LP start followed by Ipopt; an infeasible HiGHS status is independently checked with a normalized Farkas ray against the original row matrix and row/column bounds; feasible NLP terminals retain independent balances, pressure, KKT, and phase-identity audits | Log-volume difference is a dimensionless relative-volume measure, not a packing proxy. Missing, malformed, sign-inconsistent, dual-infeasible, or marginal Farkas evidence is indeterminate |
| 9 | Apply free-energy and modified-potential convergence tests, Eqs. (68)–(69) | Evaluate both paper gates and preserve independent balance, pressure, domain, KKT, and identity axes | Passing paper ratios cannot erase a failed physical certificate |
| 10 | Refine trace electrolyte compositions and recheck the equilibrium | Refine only strictly positive charged trace coordinates in \(\log_{10}x_i\), reconstruct all phases, and repeat the full physical and numerical audit | This is the only logarithmic composition coordinate. Small accurate phase or component amounts are not retired |

## Identity map

| Decision | State variables | Tolerance owner | Why it is separate |
|---|---|---|---|
| Persistent \(\mathcal M\) numerical equivalence | Modified composition and molar volume | `m_representation_equivalent`; componentwise modified-composition difference and log-volume difference `<=1e-8` | Appendix C stores \((V,\bar{\mathbf{x}}^{(EC)})\); this prevents accumulation of numerical copies |
| Step-6 candidate distinctness | Modified composition and Provider packing fraction | Perdomo Eq. (66); either difference `>=1e-3` | This is a deliberately looser phase-discovery gate, not storage identity |
| Step-8 numerical phase merge | Physical composition and molar volume | `phase_merge`; maximum physical-composition or log-volume difference `<=1e-6` | The paper removes phases with the same composition and volume after Problem (67) |

The two volume identities use a log-volume difference because
\(\lvert\log V_1-\log V_2\rvert\) measures relative rather than
unit-dependent absolute volume. It does not estimate packing fraction.

## Necessary package additions

The package adds only evidence and recovery required to make the paper-style
finite algorithm reliable for its public use:

- named, purpose-owned tolerances with raw residual diagnostics;
- stable candidate IDs and append-only \(\mathcal M\);
- independently audited LP and physical/KKT evidence rather than solver status
  alone;
- the exact feasible Step-8 LP state as the NLP initial point;
- conditional neighborhood recentering when a phase reaches its local
  boundary;
- at most one KKT-supported retirement followed by a complete reduced
  Problem-(67) re-solve;
- one cold retry after a failed warm continuation;
- exact Step-8 reuse only when the ordered immutable candidate IDs and problem
  are unchanged; and
- independently certified Step-8 phase feedback when Step 9 or 10 returns to
  Stage II.

These policies do not change the EOS, copy Provider internals, create a second
controller, claim guaranteed globality, or permit a finite but uncertified
iterate to become an equilibrium.

## Tolerance doctrine

Perdomo does not prescribe every software tolerance. Package tolerances are
therefore practical numerical policy: they should be as loose as needed for
the relevant solve to converge reliably and as tight as needed to preserve
the scientific distinction owned by the gate. A tolerance change does not
need a claim of mathematical uniqueness, but it must remain named, retain its
failure meaning, and pass the native manufactured workflow plus the public
Perdomo numerical evidence.

Tolerance adjustment must not:

- merge materially different phases or persistent states;
- reinterpret a search limit as scientific convergence;
- turn missing independent evidence into a pass; or
- clip or retire a valid state because a composition is small.

## Closed focused evidence work

The Step-5 original-physical-coordinate KKT/dual-reconstruction audit and the
Step-8 Farkas infeasibility audit are implemented as separate fail-closed
certificates. Solver status is only an evidence input in both paths.

SAFT-\(\gamma\) Mie support is a stretch goal for Provider. Equilibrium should
consume it only through a separately admitted installed Provider capability;
it must not implement or copy that EOS here.
