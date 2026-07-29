# HELD2.0 Paper Algorithm Specification

## Status and authority

This document is the normative specification for the clean-room rewrite of
the charged-mixture HELD2.0 controller in this repository. The implementation
must be traceable to this document. Existing HELD2 controller behavior is not
authority and must not be copied merely for compatibility.

The primary scientific source is:

> F. A. Perdomo et al., “Phase stability criteria and fluid-phase equilibria
> in strong-electrolyte systems,” *Computers & Chemical Engineering* 194
> (2025) 108977, especially Algorithm 1, Sections 3.2, 3.5, 4.1–4.5, and
> Appendix C.

Equation numbers in this document refer to that paper. Where the paper leaves
a numerical or software choice open, this document labels the choice
**implementation policy**. Such a policy may not be presented as a rule from
Perdomo.

This document supersedes the Stage-II same-major search behavior described in
`docs/designs/2026-07-22-held2-installed-completion.md`. In particular, HELD2.0
does not require two new candidates from one major iteration.

The first implementation state demonstrated against both the Perdomo
electrolyte cases and one Khudaida Figure-2 electrolyte LLE case is frozen in
[`2026-07-28-held2-validated-working-baseline.md`](2026-07-28-held2-validated-working-baseline.md).
That evidence checkpoint is the regression guardrail for later cleanup and
performance work; it does not replace this specification.

The paper conditions, package certificates, and deliberate implementation
policies are cross-referenced in the
[`2026-07-28-held2-necessary-condition-map.md`](2026-07-28-held2-necessary-condition-map.md).
That map explains differences; it does not relax this specification.

The 2026-07-28 Khudaida Figure-2 investigation supersedes the earlier Step-8
prohibition on a
same-major active-set re-solve. Multiplier-only deletion followed by a
phase-fraction LP removed phases carrying material and froze the surviving
compositions, so the reduced state could not satisfy the unchanged balance
gate. Step 8 now retires at most one KKT-inactive candidate at a time and
re-solves Problem (67) on the retained candidate neighborhoods. A failed
reduced solve terminates indeterminate; it is not accepted through
phase-fraction recovery and is not converted into Stage-II feedback.

## Scope

The algorithm accepts:

- temperature \(T\);
- pressure \(P_o\); and
- total species mole fractions \(\mathbf{x}_o\).

It returns:

- the detected and retained active-phase count \(mp\);
- the species mole-fraction vector \(\mathbf{x}^m\) in every phase \(m\);
- the phase fraction \(\phi^m\) of every phase; and
- diagnostics that identify the completed step, transition, work count, and
  elapsed time for each of Steps 1–10.

The detected count is not a completeness proof. As Pereira et al. (2012)
note, when more stable phases exist than components, material balance can be
satisfied by a proper subset and HELD need not locate every stable phase.
Results therefore report
`phase_enumeration_certificate = completeness_not_guaranteed` independently
of the finite-search globality label.

The algorithm owns phase discovery and equilibrium. It does not own EOS
formulations or parameters. All thermodynamic values and derivatives come
through the installed Provider native SDK.

This rewrite covers the charged HELD2.0 route. Neutral mixtures continue to
enter through the common flash dispatch but use the existing neutral HELD
formulation. Perdomo describes HELD2.0 as applicable to molecular mixtures as
well, but replacing the accepted neutral implementation is outside this
rewrite. No completion statement in this document applies the new charged
coordinate state machine to the neutral route.

## Notation and persistent state

Let \(C\) be the species count. Provider order is public result order and is
never assumed to match the paper's notation. Step 1 constructs an immutable
internal permutation:

1. retained charged species, preserving Provider order;
2. the selected eliminated charged species \(E\);
3. molecular species other than the closure species, preserving Provider
   order; and
4. the selected molecular closure species \(C\).

The inverse permutation is applied to every public composition and
component-indexed diagnostic. In the internal paper order, species
\(1,\ldots,E\) are charged, species \(E\) is eliminated, and species \(C\) is
the molecular closure species. Let

\[
\mathcal C^{(EC)} = \{1,\ldots,E-1,E+1,\ldots,C-1\}
\]

be the \(C-2\) independent modified-composition indices. Species \(C\) is the
composition-closure species.

For each retained species \(i\ne E\), Eq. (23) first defines the modified mole
number

\[
\bar n_i^{(E)} =
\left(1-\frac{z_i}{z_E}\right)n_i,
\qquad i\in\mathcal C^{(E)}.
\tag{23}
\]

The transformation preserves the total mole number:
\(\sum_{i\in\mathcal C^{(E)}}\bar n_i^{(E)}=n_t\). Eq. (30) then defines the
\(C-1\) retained modified mole fractions

\[
\bar{x}^{(E)}_i =
\left(1-\frac{z_i}{z_E}\right)x_i,
\qquad i\in\mathcal C^{(E)}.
\tag{30}
\]

Define

\[
\alpha_i=1-\frac{z_i}{z_E}.
\]

Every retained factor \(\alpha_i\), \(i\ne E\), must be strictly nonzero.
These fractions sum to one. Eliminating the closure coordinate
\(\bar x_C^{(E)}\) gives the \(C-2\) vector
\(\bar{\mathbf x}^{(EC)}\) indexed by \(\mathcal C^{(EC)}\). The inverse
transformation is

\[
x_i=\frac{\bar x_i^{(E)}}{\alpha_i},
\quad i\in\mathcal C^{(EC)},\qquad
x_C=1-\sum_{i\in\mathcal C^{(EC)}}\bar x_i^{(E)},
\]

\[
x_E=-\frac{1}{z_E}
\sum_{i\in\mathcal C^{(EC)}}z_i x_i.
\]

Here \(z_C=0\), so \(\alpha_C=1\) and \(x_C=\bar x_C^{(E)}\). Every
state used by the algorithm must satisfy nonnegativity, normalization, and
electroneutrality after this inverse transformation.

The implementation stores independent coordinates compactly as
\(r=1,\ldots,C-2\). The compact-to-paper-index map is \(i=r\) for \(r<E\)
and \(i=r+1\) for \(r\geq E\). Equations written with paper indices and arrays
stored with compact indices may not be mixed implicitly.

For every retained species \(i\in\mathcal C^{(E)}\), Eq. (28) defines the
modified electrochemical potential

\[
\bar\mu_i^{el}
\equiv
\left(
\frac{\partial\underline{\bar G}^{el}}
{\partial\bar n_i^{(E)}}
\right)_{T,P,\bar n_{j\ne i}^{(E)}}
=
\frac{\mu_i^{el}}{\alpha_i}.
\tag{28}
\]

This includes the closure species, for which
\(\alpha_C=1\) and \(\bar\mu_C^{el}=\mu_C\). In terms of the intensive
modified Gibbs function and the \(C-2\) independent gradients, Eqs. (32)–(34)
are

\[
\bar\mu_C^{el}
=
\bar G^{el}
-\sum_{j\in\mathcal C^{(EC)}}
\bar G_{\bar x_j^{(EC)}}^{el}\bar x_j^{(EC)},
\tag{33}
\]

\[
\bar\mu_i^{el}
=
\bar G^{el}
+\bar G_{\bar x_i^{(EC)}}^{el}
-\sum_{j\in\mathcal C^{(EC)}}
\bar G_{\bar x_j^{(EC)}}^{el}\bar x_j^{(EC)}
=
\bar\mu_C^{el}+\bar G_{\bar x_i^{(EC)}}^{el},
\quad i\in\mathcal C^{(EC)}.
\tag{32, 34}
\]

Steps 9 and 10 compare all applicable \(C-1\) retained modified potentials,
not merely the \(C-2\) constrained-composition gradients. The eliminated
species has no separate modified potential in this coordinate system.

## Thermodynamic and tolerance basis

**Implementation policy:**

All algorithmic energies use the Provider's one-mole reduced basis:

\[
a=\frac{\bar A^{el}}{RT},\qquad
g=\frac{\bar G^{el}}{RT},\qquad
p_v=\frac{P_oV}{RT}.
\]

Thus every \(A+P_oV\) expression below is implemented as \(a+p_v\), every
\(\bar\lambda_i\) and modified chemical potential is reduced by \(RT\), and
all energy tolerances are dimensionless. Serialization must name the reduced
basis; dimensional and reduced quantities may not share a field.

The named effective-tolerance registry is:

| quantity | reduced value | role |
|---|---:|---|
| `tpd_negative` | \(-10^{-8}\) | strict Step-2 negative-witness threshold |
| \(\epsilon_b\) | \(10^{-2}\) | Step-6 upper/lower agreement |
| \(\epsilon_\lambda\) | \(0.5\) | Step-6 derivative agreement |
| \(\epsilon_\eta\) | \(10^{-3}\) | Step-6 packing-fraction distinctness |
| \(\epsilon_x\) | \(10^{-3}\) | Step-6 composition distinctness |
| \(\epsilon_g\) | \(10^{-4}\) | Step-9 reduced free-energy convergence |
| \(\epsilon_\mu\) | \(10^{-3}\) | Step-9 reduced-potential convergence |

Perdomo 2025 states \(10^{-6}\) defaults for \(\epsilon_g\) and
\(\epsilon_\mu\). The effective values above are repository implementation
policy validated against the installed Provider route; the Stage-II values
retain the Pereira et al. (2012) HELD settings. Every effective tolerance is
emitted in diagnostics. Independent Provider-domain, pressure, KKT, balance,
charge, phase-fraction, and physical-state-equivalence tolerances remain
separately named project certificates.

The following state persists across major iterations:

- \(k\): current paper major-iteration index, initialized to zero and
  incremented only on a return to Step 4;
- \(N_{\mathrm{upper}}\): diagnostic count of completed Step-4 LP solves;
- \(\mathcal M\): all retained upper-bounding states
  \((V^m,\bar{\mathbf{x}}^{(EC),m})\);
- \(UBD^V\): current upper bound on the dual optimum;
- \(\bar{\lambda}^k\): multipliers returned by the current Step-4 LP;
- \(\bar L^{V,k}\): best Step-5 lower-level value for the current multipliers;
- \(\mathcal M^*\): current Step-6 candidate subset of \(\mathcal M\); and
- \(N_{\mathrm{start}}\): next unused ordinal in the deterministic Step-5
  start stream; and
- Stage-III feedback, when Steps 8 or 9 return control to Step 4.

The set \(\mathcal M\) is persistent. Step 7 never discards it. Step 6 searches
the entire set under the current \(UBD^V\) and \(\bar{\lambda}^k\).

## Global invariants

1. There is one scientific implementation. The native diagnostic executable
   and the public Python API call the same controller and differ only in
   serialization and invocation.
2. No Step-5 state is required to share a major iteration with another
   candidate. Candidate multiplicity is determined only by Step 6 over
   persistent \(\mathcal M\).
3. A finite multistart search is not a global proof. Every result reports
   `globality_certificate = not_guaranteed`.
4. Resource exhaustion, Provider failure, incomplete pressure-root evidence,
   solver failure, or failed physical certification is indeterminate. It may
   not be converted to a stable one-phase result.
5. Mathematical-set identity in \(\mathcal M\) uses one named,
   representation-level physical equivalence rule. It may collapse numerical
   copies of one state, but it must be tighter than the Eq. (66) phase-identity
   rule and may not collapse genuinely distinct composition or density
   branches.
6. Solver diagnostics may strengthen evidence, but they may not silently
   change an Algorithm-1 stop condition. Any additional acceptance gate must
   be named as an implementation policy and reported separately.
7. Trace collection is observational. Enabling it must not alter starts,
   solver options, transitions, results, or status.
8. The final public result is constructed from the native result. Python
   postprocessing may not introduce component-count, charge, neutral-mixture,
   midpoint-feed, or route-specific scientific gates.

## Stage I — stability test and initialization

### Step 1 — define modified mole fractions and bounds

**Paper source:** Algorithm 1; Eqs. (23), (30), and (57)–(61); Section 4.2.

**Input:** \(T\), \(P_o\), \(\mathbf{x}_o\), and species charges.

**Procedure:**

1. Validate finite positive \(T\) and \(P_o\), a finite nonnegative
   \(\mathbf{x}_o\), unit normalization, and feed electroneutrality.
2. Require at least two charged species and at least one molecular species.
   Among charged species with maximum \(|z|\), retain only candidates for
   which \(1-z_i/z_E\ne0\) for every other species \(i\). Choose the last
   admissible candidate in Provider order, matching the paper convention that
   the eliminated ion is the final charged species. If no admissible candidate
   exists, terminate indeterminate with
   `unsupported_singular_charge_transformation`; a component-order tie may
   not silently select a singular map.
3. Choose the last molecular species in Provider order as closure species
   \(C\), construct the internal permutation defined above, and retain its
   inverse.
4. Apply Eq. (30) to obtain
   \(\bar{\mathbf{x}}_o^{(EC)}\).
5. Construct the physical upper bound for each retained charged species as

   \[
   x_i^{U}=
   1-\frac{|z_i|}
   {|z_i|+
   \max_{j\in\{1,\ldots,E\},\,\operatorname{sgn}(z_j)\ne
   \operatorname{sgn}(z_i)}|z_j|},
   \quad i=1,\ldots,E-1,
   \tag{59 corrected}
   \]

   and set \(x_i^U=1\) for every independent molecular species
   \(i=E+1,\ldots,C-1\). Apply Eq. (58):

   \[
   \bar x_i^{(EC),U}=\alpha_i x_i^U.
   \]

   **Paper erratum resolution:** printed Eq. (59) searches only
   \(j\in\{1,\ldots,E-1\}\), which excludes the eliminated counter-ion and
   yields an empty opposite-sign set for a binary salt. The corrected search
   includes all charged species \(1,\ldots,E\). Printed Eq. (60) uses compact
   indices \(i=E,\ldots,C-2\), whereas the surrounding equations use paper
   species indices. Under the explicit compact map above, those are exactly
   the molecular paper indices \(E+1,\ldots,C-1\). Tests must cover both a
   binary monovalent salt and an asymmetric salt.
6. Set the finite modified lower bound from Eq. (61):

   \[
   \bar{x}^{(EC),L}_i =
   10^{-10}\left(1-\frac{z_i}{z_E}\right).
   \]

   This is a finite-search regularization, not a physical lower bound on an
   equilibrium mole fraction. It keeps the multidimensional searches away
   from singular logarithms. A charged component may cross it only in the
   Step-10 logarithmic trace solve.

7. Reject any nonfinite, nonpositive, reversed, or feed-excluding coordinate
   interval after transformation.
8. Obtain the molar-volume domain corresponding to
   \(\eta\in[10^{-6},0.74]\) from the Provider.

**Output:** immutable coordinate metadata, transformed feed, composition
bounds, and volume-domain evaluator.

**Failure:** invalid input or an empty transformed domain terminates
indeterminate at Step 1.

### Hybrid coordinate contract

HELD2 uses one coordinate policy throughout the ten steps:

| Owner | Coordinates | Domain |
| --- | --- | --- |
| Steps 1--2 | linear independent modified fractions | complete Step-1 polytope, including the finite search floor |
| Steps 3--7 | linear independent modified fractions | complete Step-1 polytope |
| Step 8 | linear phase amounts and independent modified fractions | complete Step-1 polytope and exact linear balances |
| Step 9 | physical and modified residuals evaluated from the linear Step-8 state | no coordinate change |
| Step 10 | one charged trace fraction at a time in \(\log_{10}x_i\) | \(x_i\in[10^{-300},5\times10^{-10}]\), with every non-lower-bound Step-1 constraint still enforced |

The algorithm does not logarithmically transform material-balance variables.
At trace scale, differences decisive for chemical-potential equality are below
the useful resolution of a linear balance NLP. Step 10 therefore solves only
the scalar trace chemical-potential residual in logarithmic coordinates,
reconstructs the physical phase, and then re-evaluates the ordinary and
modified linear balances. Zero is never passed to a logarithmic Provider
state.

### Step 2 — tangent-plane stability test

**Paper source:** Algorithm 1; Eqs. (49), (52)–(55), and (62)–(63);
Section 4.2.

First solve the fixed-feed volume problem from Eq. (49). Enumerate every
Provider pressure root in the Step-1 volume interval, independently certify
its pressure residual and mechanical stability, and evaluate

\[
g_o(V)=
\bar A^{el}(T,V,\bar{\mathbf{x}}_o^{(EC)})+P_oV.
\]

The homogeneous reference \(V_o\) is the unique strict interior,
mechanically stable global minimizer among those roots. Evaluate both volume
boundaries as domain checks: if a boundary has a lower or unresolved-tied
objective, the fluid interval has truncated the minimizer and Step 2 is
indeterminate. A missing root, a marginally stable root, or two
representation-distinct roots tied within the named reference-objective
tolerance is likewise indeterminate; an arbitrary root may not define the
feed tangent.

At the selected reference, set

\[
\bar\lambda_{o,i}=
\left(
\frac{\partial\bar A^{el}
(T,V_o,\bar{\mathbf{x}}_o^{(EC)})}
{\partial\bar x_i^{(EC)}}
\right)_{T,V,\bar x_{j\ne i}},
\]

and define the complete Eq. (55) tangent in Helmholtz coordinates:

\[
\begin{aligned}
T_A(T,P_o,V,V_o,\bar{\mathbf{x}}^{(EC)},
\bar{\mathbf{x}}_o^{(EC)})
=\;&
\bar A^{el}(T,V_o,\bar{\mathbf{x}}_o^{(EC)})
+P_o(V_o-V)\\
&+\sum_i\bar\lambda_{o,i}
\left(\bar x_i^{(EC)}-\bar x_{o,i}^{(EC)}\right).
\end{aligned}
\tag{55 at the feed}
\]

Minimize the tangent-plane distance

\[
d(V,\bar{\mathbf{x}}^{(EC)};
T,P_o,\bar{\mathbf{x}}_o^{(EC)})
=
\bar A^{el}(T,V,\bar{\mathbf{x}}^{(EC)})
-T_A(T,P_o,V,V_o,\bar{\mathbf{x}}^{(EC)},
\bar{\mathbf{x}}_o^{(EC)})
\tag{62}
\]

over the Step-1 composition and volume domain.

Problem (63) is a joint minimization over the \(C-2\) independent modified
compositions and molar volume. It does not constrain a trial state to
\(P(T,V,\bar{\mathbf{x}}^{(EC)})=P_o\). Consequently, a finite admissible
state with independently re-evaluated \(d<-10^{-8}\) proves that the
homogeneous reference is unstable even when that trial volume is not a
pressure-stationary root. Pressure stationarity is required for a phase used
by the later equilibrium construction, not for the logical validity of the
negative witness itself.

**Paper erratum resolution:** the typeset Eq. (62) prints the reference-state
term
\(\bar A^{el}(T,V_o,\bar{\mathbf{x}}_o^{(EC)})-T_A\).
With Eq. (55), that expression is affine in every trial variable and cannot be
the nonconvex tangent-plane distance globally minimized in Problem (63).
The original HELD Eq. (9), the definition of tangent-plane distance, and the
requirement that the objective compare the trial free-energy surface with the
feed tangent all require the trial-state term
\(\bar A^{el}(T,V,\bar{\mathbf{x}}^{(EC)})\) used above. The implementation
must test zero TPD at the feed, exact first and second derivatives away from
the feed, and agreement with the original-coordinate Gibbs TPD.

The paper uses a tunnelling global search up to

\[
N_S=10C
\]

times, with a fixed random seed. Search stops immediately when an
independently certified witness satisfies
\(d<-10^{-8}\) on the reduced basis. A candidate at or above that strict
threshold is not a negative witness.

**Transition:**

- negative TPD witness: proceed to Step 3;
- no negative witness after all declared searches: return the operational
  one-phase outcome `no_negative_witness_detected`; or
- incomplete/failed search evidence: terminate indeterminate.

**Implementation policy:** a deterministic finite global method may replace
tunnelling, but it must retain the same early-negative semantics, declare its
work budget, and report `globality_certificate = not_guaranteed`. Absence of a
negative witness is an operational finite-search result, not certified global
stability. It is returned only when every declared Step-2 search completed
successfully.

The package implements Problem (63) with DIRECT-L in the joint coordinates

\[
\left(q_V,\mathbf q_x\right)\in[0,1]^{C-1},
\qquad
\log V =
\log V_L(\bar{\mathbf{x}}^{(EC)})
+q_V\left[
\log V_U(\bar{\mathbf{x}}^{(EC)})
-\log V_L(\bar{\mathbf{x}}^{(EC)})
\right],
\]

where the Step-1 map sends \(\mathbf q_x\) to the complete admissible
modified-composition polytope and the installed Provider supplies
composition-dependent \(V_L,V_U\). Each joint trial therefore consumes one
Provider bounds query and one state evaluation at the explicit mapped volume.
The deterministic search allowance is 6500 such Provider callbacks, rather
than 50 composition points that each trigger a new 64-interval pressure-root
enumeration. Complete pressure-root enumeration remains mandatory at the feed
because selecting the homogeneous reference is a different scientific
decision.

The first strict-negative joint state stops the instability search and is
freshly re-evaluated before it can become evidence. All Step-1 polytope,
ion-domain, finite-state, explicit-volume consistency, and strict-negative
checks are repeated. Failure of that independent check is indeterminate.
The minimum raw joint TPD is retained separately from any state prepared for
the later candidate construction.

A pressure-stationary state is useful to Steps 3–8 even though it is not
needed for the instability proof. After a certified negative witness exists,
the remaining search allowance may therefore prepare a downstream seed. A
composition DIRECT-L trial samples 65 deterministic log-volume points but
does not enumerate a pressure envelope. Only when the best sampled TPD is at
most 0.2 does it refine the single most-promising pressure-residual bracket,
using at most 64 safeguarded evaluations. The 0.2 value is a routing threshold
for local refinement; it cannot accept a witness. The refined state must be a
strict mechanically stable pressure root and must independently retain
\(d<-10^{-8}\). If the selected certified witness remains off-pressure, one
complete-envelope refinement at that composition is permitted for downstream
initialization. Thus complete envelopes are never nested inside every global
composition trial.

The joint trials have an exact two-callback charge and each downstream
preparation trial has an explicit worst-case charge of
\(1+65+64=130\) Provider callbacks. DIRECT-L receives only the number of
whole trials that fit the remaining 6500-callback allowance; zero permitted
trials is invalid rather than NLopt's unlimited-work convention. Mandatory
feed-reference work, charged trace-boundary evidence, independent
recertification, and the single targeted fallback are reported in the total
Step-2 Provider count but are not allowed to shrink the declared joint-search
coverage. Every finite result continues to report
`globality_certificate = not_guaranteed`.

### Step 3 — initialize the dual problem and \(\mathcal M\)

**Paper source:** Algorithm 1; Section 4.2; Appendix C.

Set \(k=0\), \(N_{\mathrm{upper}}=0\), and \(N_{\mathrm{start}}=0\). Evaluate
the homogeneous feed state and initialize

\[
UBD^V=\bar G_o^{el}.
\]

Initialize the lower value to
\(\bar L^V=-\infty\), matching Algorithm 1’s lower-bound initialization and
the cited Pereira implementation.

Initialize \(\mathcal M\) with the feed state and the Appendix-C bounding
states.

Appendix C's printed physical-coordinate formulas do not define the claimed
bracketing vectors. In particular, its \(1/[2(C-1)]\) upper prefactor places
both nominal bounding points below the feed for the four-component Table 5
mixtures. Replacing only that prefactor by \(1/2\) makes the neutral-solvent
upper point leave the electroneutral composition simplex. Clipping either
result would invalidate the corresponding dual cut.

The implementation therefore constructs the required opposing cuts directly
in the modified-fraction polytope introduced by Eq. (30). Modified fractions
already incorporate electroneutrality and are the coordinates of the
\(C-2\)-dimensional dual problem. For every independent coordinate \(i\),
find the exact distances \(\alpha_i^-\) and \(\alpha_i^+\) from the feed to
the Step-1 polytope boundary along the negative and positive coordinate axes:

\[
\alpha_i^\pm =
\max\left\{
\alpha\geq0:
\bar{\mathbf{x}}_o^{(EC)}
\alpha(\pm\mathbf e_i)\in\bar X^{(EC)}
\right\}.
\]

Use the interior midpoint on each side:

\[
\hat{\bar{\mathbf{x}}}^{(EC),i}
=\bar{\mathbf{x}}_o^{(EC)}
-\frac{\alpha_i^-}{2}\mathbf e_i,
\qquad
\tilde{\bar{\mathbf{x}}}^{(EC),i}
=\bar{\mathbf{x}}_o^{(EC)}
+\frac{\alpha_i^+}{2}\mathbf e_i.
\]

This retains Appendix C's operational requirement—one lower and one upper
dual cut per coordinate—while respecting every constraint that defines the
actual Problem (65) domain, including a Provider-declared ion ceiling. Before
any pressure solve, invert Eq. (30) and require the complete physical vector
to lie inside the Step-1 polytope, including
finite lower and upper bounds, closure-species bounds, reconstructed-ion
nonnegativity, normalization, and electroneutrality. Do not clip or project an
Appendix-C state. A violation terminates indeterminate with
`appendix_c_state_outside_domain`, because a cut outside Problem (65)'s domain
is not a valid bound on its dual function.

For every valid vector, solve
\(P(T,V,\bar{\mathbf{x}}^{(EC)})=P_o\), and insert the resulting
\((V,\bar{\mathbf{x}}^{(EC)})\) into \(\mathcal M\). Step 3 therefore
contributes the feed plus \(2(C-2)\) Appendix-C bounding states as the
paper-defined initialization.

**Implementation policy for multivariate finite search:** when \(C-2>1\),
also retain the independently certified negative-TPD state returned by Step 2
as one additional physical cut in \(\mathcal M\). The axis-aligned Appendix-C
states bracket each multiplier direction but do not represent the interior
multivariate witness that established instability. In one dimension the two
Appendix-C states already bracket the complete interval, so the witness is not
inserted. This is a declared finite-search acceleration; it does not replace
an Appendix-C state, change the lower-bound semantics, or upgrade
`globality_certificate="not_guaranteed"`.

**Paper erratum resolution:** the prose of Appendix C says
“for each \(i\in\mathcal C\),” while the equations, the \(C-2\)-dimensional
dual problem, and Section 4.2’s explicit count require
\(i\in\mathcal C^{(EC)}\). This specification uses the latter. Section 4.2
describes the Appendix-C construction as providing bounds on
\(\bar\lambda\), while Appendix C says the pressure-root pairs are stored in
\(\mathcal M\). Eq. (64) leaves \(\bar\lambda\in\mathbb R^{C-2}\), and neither
Section 4.2 nor Appendix C gives formulas for independent numerical box
bounds. The specification therefore preserves both statements as follows:
the \(2(C-2)\) states are inserted into \(\mathcal M\), and their opposing cut
slopes must make the first Eq. (64) LP bounded in every multiplier direction.
Those cuts provide implicit two-sided bounds on the optimizing multipliers;
one LP optimizer is not itself called a pair of numerical bounds. No separate
finite box may be invented without a cited formula and a revision to this
specification.

The rewrite must test lower/upper bracketing, exact simplex closure, the
inverse Eq. (30) lift, the one-dimensional no-extra-cut rule, and the
four-component Table 5 coordinate layout with its retained interior witness.

**Implementation policy:** if a constructed composition has multiple pressure
roots, select the unique lowest-objective strict-stable root. A missing root or
unresolved objective tie terminates indeterminate. Do not insert arbitrary
perturbations or synthetic points into \(\mathcal M\).

## Stage II — identify candidate stable phases

### Step 4 — solve the upper-level LP

**Paper source:** Algorithm 1; Eq. (64); Section 4.3.

Solve

\[
\begin{aligned}
UBD^V = \max_{v,\bar{\lambda}\in\mathbb R^{C-2}}\quad &v\\
\text{s.t.}\quad
v \leq&
\bar A^{el}(T,V^m,\bar{\mathbf{x}}^{(EC),m})+P_oV^m\\
&+\sum_i\bar\lambda_i
\left(\bar x_{o,i}^{(EC)}-\bar x_i^{(EC),m}\right),
\quad \forall m\in\mathcal M,\\
v\leq&\bar G_o^{el}.
\end{aligned}
\tag{64}
\]

Set \(\bar\lambda^k=\bar\lambda^*\) from this LP and retain the full LP
certificate before entering Step 5. Increment
\(N_{\mathrm{upper}}\leftarrow N_{\mathrm{upper}}+1\); do not increment \(k\)
in Step 4.

**Output:** current \(UBD^V\), \(\bar\lambda^k\), active cuts, residuals, and
solver status.

**Failure:** a nonoptimal, infeasible, unbounded, or uncertified LP terminates
indeterminate at Step 4.

### Step 5 — solve the nonconvex lower-level problem

**Paper source:** Algorithm 1; Eq. (65); Section 4.3.

For fixed \(\bar\lambda^k\), solve

\[
\begin{aligned}
\bar L^{V,k} =
\min_{V,\bar{\mathbf{x}}^{(EC)}}\quad&
\bar A^{el}(T,V,\bar{\mathbf{x}}^{(EC)})+P_oV\\
&+\sum_i\bar\lambda_i^k
\left(\bar x_{o,i}^{(EC)}-\bar x_i^{(EC)}\right)
\end{aligned}
\tag{65}
\]

over the Step-1 domain.

Run local solves from one persistent deterministic sequence generated with a
fixed seed. Start ordinal \(N_{\mathrm{start}}\) is consumed exactly once and
never reset at a major boundary or Stage-III return. The sequence includes the
paper's random interior starts and the retained HELD shifted-previous and
near-pure profiles; the exact schedule and resource budget are named
implementation policy and are emitted in diagnostics.

Track the best independently certified local solution. Certification requires
a finite in-domain physical state, pressure stationarity, primal feasibility,
original-coordinate first-order KKT stationarity, valid bound-multiplier
signs, complementarity, and coordinate-to-physical dual reconstruction. A
library success status alone is not a solved Problem (65).

Stop the multistart search as soon as the best certified local solution
satisfies

\[
\bar L^{V,k}\leq UBD^V.
\]

Add that best terminal
\((V^k,\bar{\mathbf{x}}^{(EC),k})\) to persistent \(\mathcal M\). If it is
representation-equivalent to an existing member, the mathematical set remains
unchanged and the insertion is recorded as `equivalent_member`; this may not
be reported as a new cut.

**Implementation policy:**

- solve Eq. (65) directly in the \(C-2\) modified-composition variables and
  one volume variable;
- use exact Provider derivatives and an exact-Hessian local NLP;
- use the repository’s named physical KKT and dual-reconstruction tolerances
  for the independent local-solution certificate;
- maintain one tight physical-state equivalence rule for \(\mathcal M\)
  membership, expressed in modified composition and molar volume, matching
  Appendix C's stored \((V,\bar{\mathbf{x}}^{(EC)})\) pair. Compare volume
  through log-volume difference as a dimensionless relative-volume measure;
  keep Eq. (66) Provider-packing-fraction distinctness as a separate, looser
  rule;
- if the qualifying solution is representation-equivalent, proceed through
  Steps 6–7 without claiming a new cut; the next major consumes the next
  unused start ordinal rather than replaying the same search; and
- if the declared multistart budget ends first, terminate indeterminate at
  Step 5.

Within one invocation, this step does not seek two candidates and does not
continue merely to find a different basin after its stop condition has been
met. Algorithm 1 can return through Step 7 with an unchanged \(\mathcal M\);
the persistent start ordinal makes that later invocation new work. The
Khudaida Figure-2 audit confirmed that forcing all 128 default starts after
the paper stop condition adds 238--285 seconds per feed and changes later
candidate selection; it is prohibited.

### Step 6 — search all of \(\mathcal M\)

**Paper source:** Algorithm 1; Eq. (66); Section 4.3.

For every point \(m\in\mathcal M\), evaluate its current lower function value

\[
\bar L^{V,m} =
\bar A^{el}(T,V^m,\bar{\mathbf{x}}^{(EC),m})+P_oV^m+
\sum_i\bar\lambda_i^k
\left(\bar x_{o,i}^{(EC)}-\bar x_i^{(EC),m}\right).
\]

Populate \(\mathcal M^*\) with points satisfying all applicable Eq. (66)
criteria:

1. upper/lower agreement:

   \[
   |UBD^V-\bar L^{V,m}|\leq\epsilon_b;
   \]

2. fixed-\((T,V,\bar x_{j\ne i})\) modified-composition derivative agreement
   for the selected component set \(\mathcal I^m\):

   \[
   \left|
   \frac{\partial\bar A^{el}}{\partial\bar x_i^{(EC)}}-
   \bar\lambda_i^k
   \right|
   \leq\epsilon_\lambda|\bar\lambda_i^k|;
   \]

3. pairwise phase distinctness: for each other candidate \(n\), either

   \[
   |\eta^m-\eta^n|\geq\epsilon_\eta
   \]

   or at least one independent modified composition satisfies

   \[
   |\bar x_i^{(EC),m}-\bar x_i^{(EC),n}|\geq\epsilon_x.
   \]

A coordinate at its Step-1 lower bound is not included in
\(\mathcal I^m\).

Every point in \(\mathcal M\) is eligible for this test, regardless of whether
it came from the feed, Appendix C, or Step 5. There is no same-major origin
gate, Step-5-KKT origin gate, or requirement for two fresh points.

Packing fraction \(\eta\) must come from the Provider packing-fraction
endpoint. Log-volume distance is not Eq. (66).

**Transition:**

- if \(mp=|\mathcal M^*|\geq2\), proceed to Step 8;
- otherwise proceed to Step 7.

**Implementation policy:** the paper leaves
\(\epsilon_b,\epsilon_\lambda,\epsilon_\eta,\epsilon_x\) user-defined. The
rewrite must use the repository’s named tolerance registry, expose the
effective values in diagnostics, and preserve the mathematical relation of
each criterion. Any additive floor used when \(\bar\lambda_i^k=0\) must be
named and reported rather than hidden in a generic comparison.

The paper says \(C_{pp}\) is “usually” \(C-2\) for \(C\leq5\) and may be
smaller for larger mixtures. The rewrite uses every non-bound independent
coordinate for every supported \(C\), so \(C_{pp}=|\mathcal I^m|\). A reduced
component subset is not introduced without a separate specification.

Construct \(\mathcal M^*\) deterministically. First evaluate the non-identity
criteria for every \(\mathcal M\) member, then visit eligible members by
stable insertion ID. Build a deterministic maximal pairwise-distinct subset:
retain the first eligible member, then retain a later member only when it
satisfies Eq. (66)'s packing-fraction-or-composition distinctness test against
every already-retained member. Tolerance proximity is not called an
equivalence relation and union-find clustering is forbidden because the
pairwise relation is not transitive. Retain every comparison and keep/drop
decision in diagnostics. A comparison in the unresolved tolerance margin is
indeterminate, not an arbitrary keep/drop decision.

### Step 7 — advance the major iteration

**Paper source:** Algorithm 1; Section 4.3.

Increment \(k\leftarrow k+1\) exactly once and return to Step 4 with the full
persistent \(\mathcal M\).

The reported “Iter.” count is \(N_{\mathrm{upper}}\), the number of Step-4
upper-level problems solved. Local NLP iterations and multistart attempts are
separate work counters.

**Implementation policy:** a configured major-iteration limit is a resource
limit, not part of HELD2.0. The native validation profile permits 80 upper
solves and 128 Step-5 starts. The upper-solve budget is shared across feeds
and accounts for the installed ePC-SAFT model rather than treating Table 5's
SAFT-\(\gamma\) Mie iteration counts as an equivalent resource bound.
Reaching either limit terminates indeterminate and reports the last completed
step.

## Stage III — acceleration and convergence

### Step 8 — minimize total Gibbs free energy

**Paper source:** Algorithm 1; Eq. (67); Section 4.4.

Starting from the \(mp\) members of \(\mathcal M^*\), solve

\[
\min_{\mathbf V',\bar{\mathbf x}',\boldsymbol\phi}
\sum_{m=1}^{mp}\phi^m
\left[\bar A^{el}(T,V'^m,\bar{\mathbf x}'^m)+P_oV'^m\right]
\tag{67}
\]

subject to:

- modified material balances
  \(\sum_m\phi^m\bar x_i'^m=\bar x_{o,i}^{(EC)}\);
- \(\sum_m\phi^m=1\);
- \(\phi^m\in[0,1]\);
- \(V'^m\in[V^L,V^U]\); and
- each modified composition is within \(10^{-3}\) of its Step-6 candidate,
  intersected with the Step-1 physical bounds.

Before invoking the nonlinear optimizer, solve an exact linear
material-balance feasibility problem. Let \(X_m\) be the Step-1 transformed
composition polytope intersected with candidate \(m\)'s \(10^{-3}\)
neighborhood, represented as \(B_m\bar{\mathbf x}\leq b_m\). Introduce
\(w_i^m=\phi^m\bar x_i^m\) and solve

\[
\sum_m\phi^m=1,\qquad
\sum_m w_i^m=\bar x_{o,i}^{(EC)},\qquad
B_m\mathbf w^m\leq b_m\phi^m,\qquad
\phi^m\geq0.
\tag{67-feasibility}
\]

This perspective formulation is an exact feasibility test for the bounded
composition neighborhoods; every Step-1 linear physical-domain inequality is
included. The volume intervals are checked separately for nonemptiness.
Only a certified LP infeasibility status with a validated Farkas certificate
is `certified_infeasible`. A feasible presolve supplies a deterministic
Problem-(67) start. Ipopt failure is never an infeasibility certificate.

**Transition:**

- certified infeasibility from Problem (67-feasibility): increment \(k\) and
  return to Step 4;
- feasible, independently certified nonlinear solve: finalize the active
  phase set as below, update \(mp\), and proceed to Step 9.

A Provider failure, optimizer failure, invalid terminal, or inability to
certify either feasibility or infeasibility terminates indeterminate at
Step 8. It is not Stage-II feedback.

Duplicate removal is part of Step 8. It may not be performed earlier to
manufacture a Step-6 candidate count.

**Implementation policy for duplicate phases:**

1. Visit solved phases by stable candidate ID. Build a greedy maximal
   pairwise-distinct set using the named `phase_merge` test in physical
   composition and molar volume, matching the paper's “same composition and
   volume” duplicate rule. Use log-volume difference only as a dimensionless
   relative-volume measure; do not apply a nontransitive equivalence-class
   algorithm. Retain the larger-weight numerical representative of a
   duplicate pair; break equal-weight ties by stable ID.
2. Sum a removed numerical duplicate's phase fraction into its retained
   representative and re-run the full balance and physical certificates.
3. Do not retire a distinct phase because its amount or any component
   fraction is small. Trace component fractions remain valid linear Step-8
   variables and do not trigger phase identity changes.
4. After a full KKT audit passes, a positive lower-bound phase-amount
   multiplier above the named retirement margin may nominate one inactive
   candidate. Remove at most that one candidate and re-solve Problem (67) on
   the retained neighborhoods, preserving their complete composition freedom.
5. Accept the reduced active set only after the reduced solve independently
   passes every nonlinear, KKT, balance, physical, and phase-identity gate.
   Preserve the reduced candidate-ID/continuation mapping. If the reduced solve
   fails, terminate indeterminate at Step 8.
6. A supplied active-set continuation is an optimizer initial condition, not
   scientific evidence. If that warm solve fails any Step-8 acceptance gate,
   retry the identical candidate problem once from the deterministic cold
   initialization. Combine resource accounting across both attempts and
   terminate indeterminate if the cold attempt also fails.

**Implementation policy for repeated Stage-III work:** persistent
\(\mathcal M\) is append-only, and a stable insertion ID identifies an
immutable physical cut. After a Stage-III return, initialize the next
Problem-(67) candidate problem from the previously retained active IDs plus
the newest currently eligible Step-6 ID. If that selection contains fewer
than two candidates, use the complete current \(\mathcal M^*\). Surviving
candidate neighborhoods may be recentered only under the boundary rule above.

If the resulting ordered stable-ID vector is exactly unchanged from the
previous Step-8 solve, reuse that certified Problem-(67) result instead of
replaying the same nonlinear problem. This is exact memoization, not scientific
acceptance: the feed and immutable candidate states are unchanged, and
Problem (67) contains neither the Step-4 multipliers nor \(UBD^V\). Step 9 is
still executed against the current Step-4 upper bound. Any new cut, changed
eligibility, changed stable-ID order, or uncertified previous result requires
a new Step-8 solve.

Step 8 remains in linear phase amounts and modified composition coordinates.
Logarithmic coordinates would remove the exact zero needed for inactive phase
amounts and would make the exact linear material balances nonlinear. Only the
Step-10 charged trace refinement uses the bounded logarithmic coordinate
specified by the hybrid coordinate contract.

### Step 9 — convergence tests

**Paper source:** Algorithm 1; Eqs. (68)–(69); Section 4.4.

Apply both tests to the Step-8 solution.

The reported Eq. (68) gap remains signed. Its numerical lower-bound test
allows negative roundoff up to \(\epsilon_g\); larger negative gaps fail.

Order active phases by stable candidate ID before evaluating adjacent pairs.
This ordering is deterministic and does not change the all-phases equality
represented by adjacent comparisons.

Free-energy convergence requires

\[
0\leq UBD^V-\underline{\bar G}^{el,*}\leq\epsilon_g,
\qquad \epsilon_g=10^{-4}\ \text{by implementation policy}.
\tag{68}
\]

Modified chemical-potential convergence requires, for every retained species
and adjacent phase pair,

\[
\left|
\frac{
\bar\mu_i^{el}(T,P_o,\bar{\mathbf x}^{*,m})-
\bar\mu_i^{el}(T,P_o,\bar{\mathbf x}^{*,m+1})
}{
\bar\mu_i^{el}(T,P_o,\bar{\mathbf x}^{*,m})
}
\right|
\leq\epsilon_\mu,
\qquad \epsilon_\mu=10^{-3}\ \text{by implementation policy}.
\tag{69}
\]

**Implementation policy for Eq. (69):** use the paper’s denominator without a
hidden scale floor. If the denominator is exactly zero, the comparison passes
only when the numerator is also exactly zero; otherwise it fails and returns
to Step 4. Record numerator, denominator, and ratio for every comparison. A
separate symmetric or absolute project certificate may be reported but may
not be substituted silently for Eq. (69).

**Transition:**

- Eq. (68), or Eq. (69) for a non-trace component, violated: increment \(k\)
  and return to Step 4 without discarding \(\mathcal M\). Retain the
  independently certified Step-8 phase states as cutting planes, and
  recenter the next Problem-(67) neighborhoods for surviving candidate
  identities on those refined states. The radius remains the paper's fixed
  \(10^{-3}\); feedback advances the local neighborhoods rather than
  widening them;
- both tests pass: proceed to Step 10;
- Eq. (69) fails only for charged components pinned to the finite Step-1
  search floor: proceed to Step 10 for bounded logarithmic refinement and
  require complete Step-9 recertification afterward.

**Implementation policy:** material balance, electroneutrality, pressure,
finite-state, solver KKT, and active-phase checks remain mandatory
fail-closed result certificates. They must be reported separately from the
paper’s Eqs. (68)–(69), so a stricter project certificate is never mislabeled
as a Perdomo transition criterion. Unavailable Provider evidence, solver
failure, or failure of an auxiliary project certificate terminates
indeterminate. These transitions are typed controller actions; reason strings
are diagnostic evidence and do not control the state machine.

### Step 10 — refine trace components

**Paper source:** Algorithm 1; Eqs. (70)–(72); Section 4.4.

Revert to original mole fractions. For each component \(i\) that remains at
the finite modified lower bound in a phase \(m\), define

\[
y_i^{(EC)}=\log_{10}x_i^{(EC)}.
\]

Hold all other transformed mole fractions at their Step-8 values:

\[
y_j^{(EC)}=
\log_{10}
\left(
\frac{\bar x_j^{(EC),m}}{1-z_j/z_E}
\right),
\qquad j\in\mathcal C^{(EC)},\ j\ne i.
\tag{71}
\]

**Paper erratum resolution:** printed Eq. (71) uses \(z_i\) in the
denominator for every fixed component \(j\). Inverting Eq. (30) for component
\(j\) requires its own factor \(1-z_j/z_E\), as written above. The rewrite
must test the inverse transformation component by component for unequal ion
charges; it may not copy the printed free-index mismatch.

Refine the trace coordinate over

\[
x_i\in[10^{-300},5\times10^{-10}]
\]

by solving the paper’s one-dimensional problem

\[
\min_{y_i^{(EC)}\in Y}
\left[
\bar\mu_i^{el}(T,P_o,\mathbf y^{(EC)})
-\bar\mu_i^{el,*}
\right],
\tag{72}
\]

so its modified chemical potential matches the value
\(\bar\mu_i^{el,*}\) in a phase where the component is not at the lower
bound.

Only the charged coordinate being refined may fall below its Step-1 finite
search floor. Upper bounds, closure-species positivity, eliminated-ion
nonnegativity, the Provider ion-domain ceiling, and every other Step-1
polytope constraint remain active. Molecular coordinates never receive this
exception.

After every refinement, reconstruct the full species vector and recheck
normalization, electroneutrality, material balance, pressure, and the Step-9
chemical-potential criterion. A refinement failure is indeterminate; the
unrefined lower-bound state is not silently accepted.

If no component is at the finite lower bound, Step 10 is a recorded no-op.

**Implementation policy:** Eq. (72) is written as a one-dimensional
minimization of the signed chemical-potential difference. The numerical owner
must solve the corresponding zero-residual condition within the stated
bounded log interval and report the final residual. It may not minimize the
signed residual toward a bound. When multiple phases contain component \(i\)
above the lower bound, use the active phase with the greatest physical mole
fraction of \(i\); resolve an exact tie by stable phase ID.

If no active phase contains component \(i\) above the finite lower bound, the
paper provides no reference \(\bar\mu_i^{el,*}\). Terminate indeterminate with
`trace_reference_absent`; do not invent a reference or silently retain the
finite bound.

Step 10 does not invent a material-balance correction absent from the paper.
After all trace refinements, reconstruct closure species and evaluate the
complete physical certificate. If the trace changes disturb material balance
or any other required invariant beyond its named tolerance, terminate
indeterminate and retain the residual evidence.

## Controller pseudocode

```text
held2(T, P0, x0):
    s1 = step1_define_coordinates_and_bounds(T, P0, x0)
    s2 = step2_stability_test(s1)
    if s2 == no_negative_witness_detected:
        return operational_homogeneous_result(s1, s2)
    if s2 != negative_witness:
        return indeterminate(step=2)

    state = step3_initialize_persistent_M(s1, s2)

    loop:
        s4 = step4_upper_lp(state.M, state.feed_gibbs)
        if not s4.certified:
            return indeterminate(step=4)
        state.UBD = s4.UBD
        state.lambda = s4.lambda
        state.upper_solve_count += 1

        s5 = step5_lower_multistart(state)
        if not s5.completed:
            return indeterminate(step=5)
        state.M.insert_by_representation(s5.best_terminal)

        s6 = step6_search_all_M(state)
        if size(s6.M_star) < 2:
            if major_budget_exhausted(state):
                return indeterminate(step=7, reason=resource_limit)
            step7_continue(state)
            continue

        s8 = step8_total_gibbs(state, s6.M_star)
        if s8 == certified_infeasible:
            step7_continue(state)
            continue
        if s8 != certified_feasible:
            return indeterminate(step=8)
        if s8.active_phase_count < 2:
            step7_continue(state)
            continue

        s9 = step9_convergence(state, s8)
        if s9 == paper_convergence_failed:
            state.M.insert_independently_certified(s8.active_phases)
            step7_continue(state)
            continue
        if s9 != certified:
            return indeterminate(step=9)

        s10 = step10_refine_trace_components(state, s8)
        if not s10.certified:
            if s10 requests_stage_ii:
                state.M.insert_independently_certified(s8.active_phases)
                step7_continue(state)
                continue
            return indeterminate(step=10)
        return accepted_multiphase_result(s10)
```

`step7_continue` performs the single
\(k\leftarrow k+1\) transition required by Step 7 or by Stage-III feedback.
The independent `upper_solve_count` is the Table-5 “Iter.” value.

## File and dependency design

The rewrite has one substantive C++ owner per paper step:

```text
held2_step1.{hpp,cpp}   coordinates and bounds
held2_step2.{hpp,cpp}   tangent-plane stability
held2_step3.{hpp,cpp}   persistent-M initialization
held2_step4.{hpp,cpp}   upper-level LP
held2_step5.{hpp,cpp}   lower-level multistart NLP
held2_step6.{hpp,cpp}   all-M Eq. (66) search
held2_step7.{hpp,cpp}   iteration/resource transition
held2_step8.{hpp,cpp}   total-Gibbs NLP and duplicate removal
held2_step8_nlp.cpp     Problem (67) Ipopt implementation
held2_step9.{hpp,cpp}   Eqs. (68)–(69) and physical certificates
held2_step10.{hpp,cpp}  trace-component refinement
held2_algorithm.{hpp,cpp}  thin state machine only
```

Shared low-level primitives may remain in focused modules:

- installed Provider model transport and lifetime;
- EOS state and derivative evaluation;
- modified-coordinate transformation;
- pressure-root enumeration;
- HiGHS and Ipopt adapters;
- progress events and result serialization; and
- named tolerance contracts.

No step file may call a later step directly. It returns a typed result to the
thin state machine. No pybind file owns scientific control flow.

The development-only native diagnostic and public Python API call the same
native Steps 1--10 implementation. The replaced controller is deleted; two
permanent scientific engines are forbidden.

## Diagnostics and timing

Every step records:

- invocation count;
- wall and CPU elapsed time;
- Provider evaluation count;
- optimizer solve count;
- optimizer iteration count;
- generated-start count where applicable;
- input/output state sizes;
- terminal status and reason; and
- next step.

Step 4 separately reports upper-level LP solves. This is the quantity reported
as “Iter.” in Perdomo Table 5. Step 5 reports local starts and local optimizer
iterations separately.

The native `--trace` stream emits start and completion events for every step
in real time. The final JSON contains the same structured timing records
whether or not the terminal trace was enabled.

No serializer may replay Provider evaluations or optimizer iterates merely to
construct diagnostics.

## Verification and cutover criteria

### Specification traceability

- Every controller transition has a test citing the corresponding section of
  this document.
- Every paper-derived equation in code cites its equation number.
- Every non-paper numerical choice is present in the implementation-policy
  sections above.
- A review must find no same-major candidate quota, candidate-origin gate,
  log-volume proxy for the Eq. (66) packing fraction \(\eta\), arbitrary
  Step-3 cuts, or Python-only scientific gate. Log-volume remains the
  relative-volume coordinate for identities whose paper variable is \(V\).

### Step tests

Each step has manufactured tests for:

- its successful output;
- each outgoing transition;
- invalid and incomplete evidence;
- deterministic replay;
- trace on/off result identity; and
- resource exhaustion where applicable.

Step 3 tests the exact \(1+2(C-2)\) Appendix-C states and the additional
multivariate Step-2 witness. Step 5 tests stopping at the first best value
satisfying \(\bar L^{V,k}\leq UBD^V\) and advancing the persistent start
ordinal by exactly the starts consumed. Step 6 tests reuse of old
\(\mathcal M\) members under new multipliers and packing-fraction
distinctness.

### End-to-end native tests

The native executable is the first acceptance route. It must:

1. load the exact installed Provider model through the Provider SDK;
2. expose the Provider fingerprint without parsing Provider model JSON;
3. emit real-time Steps 1–10 progress;
4. emit one JSON result from the shared native result structure; and
5. preserve `globality_certificate = not_guaranteed`.

### Table-5 input provenance gate

Perdomo Table 5 publishes three LiCl molalities, the two equilibrium phase
compositions and volumes, the upper-level solve count, and CPU time. It does
not publish the algorithm input \(\mathbf x_o\) or the equilibrium phase
fraction. Molality fixes the salt-to-water ratio but does not determine the
amount of 1-butanol.

The cited Al-Sahhaf and Kapetanovic experiment mixed equal volumes of
1-butanol and aqueous salt solution. Its paper reports phase mass fractions
but does not report every stock-solution density needed to convert that
preparation into a unique species mole-fraction vector. Perdomo also does not
state that its selected Table-5 algorithm inputs are those experimental
equal-volume preparations.

Consequently, no overall Table-5 feed vector is authoritative yet. The three
user-authorized vectors reconstructed by applying a lever rule to the rounded
Table-5 phase endpoints are useful synthetic diagnostics, but the unpublished
phase fraction remains an assumed quantity. They must be labelled
`reverse_engineered_synthetic`, not source-reported or uniquely reconstructed,
and cannot satisfy the product cutover gate.

Before a Table-5 run is accepted, its input record must contain:

- the exact source or derivation of \(\mathbf x_o\);
- the water molar mass and all density data used in a volume-to-mole
  conversion;
- the resulting LiCl molality recomputed from \(\mathbf x_o\);
- electroneutrality and normalization residuals; and
- a statement of whether the feed is source-reported, uniquely reconstructed,
  or deliberately synthetic.

The user-selected product cutover gate remains three source-qualified Table-5
inputs, one each for 4.58, 4.95, and 5.74 mol/kg, reaching accepted two-phase
native results through Steps 1–10 without hidden postprocessing. This is a
product-validation requirement, not proof that the algorithm is wrong if a
different EOS predicts different phase topology. Until the three
\(\mathbf x_o\) records exist, this gate is explicitly blocked rather than
satisfied with fabricated feeds.

For every run, per-step timings and work counts must show no unexplained
repeated search or diagnostic replay.

Perdomo’s published SAFT-\(\gamma\) Mie results are comparison references, not
numerical acceptance targets for the installed ePC-SAFT model:

| molality | paper upper-level solves | paper CPU / s |
|---:|---:|---:|
| 4.58 | 59 | 2.58 |
| 4.95 | 54 | 2.60 |
| 5.74 | 61 | 2.76 |

Different EOS formulations and parameter sets preclude exact phase-composition
parity. Algorithm-flow, work-count meaning, conservation, equilibrium
certificates, deterministic replay, and absence of wasted work remain valid
acceptance criteria.

The native Release-build audit of those three synthetic feeds at 298.15 K and
101300 Pa produced:

| molality label | native outcome | Step-4 solves | wall / s | Step 5 / s | Step 8 / s |
|---:|---|---:|---:|---:|---:|
| 4.58 | accepted two-phase | 63 | 14.33 | 6.70 | 6.36 |
| 4.95 | indeterminate: active-set balance | 60 | 16.03 | 6.50 | 8.24 |
| 5.74 | indeterminate: 64-solve resource limit | 64 | 13.76 | 7.06 | 5.40 |

These are diagnostic results for the installed, unfitted ePC-SAFT model, not
reproductions of the SAFT-\(\gamma\) Mie phase compositions. The dominant
measured cost is not repeated HELD2 control flow. Provider commit `e545ed2`
adds an append-only installed value endpoint that uses the same CppAD
zero-order tape as the derivative endpoint and returns before Jacobian,
site-sensitivity, and Hessian work. Step 8 uses it only for objective callbacks;
all gradients, Hessians, KKT checks, pressure polish, and Steps 9–10 retain the
full endpoint. This preserves the accepted 4.58 result bit for bit while
reducing Step 8 by 5.7–7.2 seconds. Step 5 still costs 6.5–7.1 seconds and the
remaining exact Step-8 derivative work costs 5.4–8.2 seconds, so the requested
sub-10-second wall target is not met. No unexplained controller replay or
failed shortcut remains.

### Public Python cutover

After the three native synthetic diagnostics established unchanged scientific
outcomes and the accepted 4.58 case preserved every reported phase and
certificate field:

1. route the public Python API to the same new native controller;
2. compare native and Python scientific payloads field-for-field;
3. remove the old public charged serializer and postprocessing fork; and
4. run installed-wheel validation as the public route.

The installed-wheel public route returns the same accepted 4.58 two-phase
result with 63 upper-level solves. The native and Python JSON contracts differ
only in nondeterministic wall and CPU measurements from independent runs. The
native payload owns phase pressure, chemical potentials, molar volume, and
total free energy, so Python performs typed conversion without a second
Provider evaluation or scientific reclassification.

Python may convert native structures to public Python objects. It may not
reclassify phases, reject feeds by component count or charge pattern, invent a
neutral/electrolyte postprocessing gate, or recompute scientific results.

## Explicitly rejected legacy behaviors

The rewrite must not contain:

- a target of two Step-6 candidates from one major iteration;
- continuation after Step 5 solely to find another basin;
- Step-6 eligibility restricted to Step-5 KKT-certified origins;
- Step-6 search restricted to current-major terminals;
- Stage-I witnesses inserted unconditionally or used to replace Appendix-C
  cuts;
- arbitrary coordinate perturbations inserted in place of Appendix C;
- infinite multiplier bounds made finite by unrelated synthetic cuts;
- basin-identity tolerances used to suppress nonidentical \(\mathcal M\)
  members;
- exact-bit comparison used to accumulate numerical copies of one
  \(\mathcal M\) member;
- log-volume distance substituted for Provider packing fraction in Eq. (66);
- Provider packing fraction substituted for Appendix C's persistent
  composition-and-volume identity;
- a 24-major scientific default presented as part of Perdomo HELD2.0;
- Provider or optimizer replay for JSON construction;
- an overall Table-5 feed fabricated from molality plus an arbitrary water,
  butanol, midpoint, or phase-fraction choice;
- a Python midpoint feed, binary-only gate, neutral-only gate, or charged
  postprocessing fork; or
- permanent dual implementations.

## Completion claim

Implementation is complete only when the old charged HELD2 scientific
controller has been removed, every charged production path uses the new
Steps 1–10 state machine, the native and Python payloads agree, the three
Table-5 native cases satisfy the cutover criteria, and the full repository
verification and cleanup audit pass. The accepted neutral HELD route remains
behind the same public flash dispatch and is not duplicated or rewritten by
this scope. No result from this finite search may claim guaranteed globality
or complete stable-phase enumeration.
