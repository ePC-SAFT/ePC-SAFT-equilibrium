# Reactive-LLE Benchmark Review

Status: source review and evidence selection; no runtime authority

Date: 2026-07-29

## Purpose

This review selects evidence subjects for the first GREPE implementation from
the local Zotero collection `Phase + Chemical Equilibrium`. A useful subject
must exercise simultaneous reaction and liquid-phase allocation and must
provide enough thermodynamic inputs and reference outputs to distinguish an
algorithm failure from missing model data.

The review does not copy paper parameters into Equilibrium. Provider parameter
bundles remain installed Provider artifacts, while reproduction data and
black-box results belong in Validation.

## Ranked subjects

### 1. Ascani--Senina pentyl-acetate reactive LLE

**Selected physical neutral-reactive-LLE subject.**

Ascani, Sadowski, and Held, “Simultaneous Predictions of Chemical and Phase
Equilibria in Systems with an Esterification Reaction Using PC-SAFT,”
*Molecules* 28 (2023) 1768,
[doi:10.3390/molecules28041768](https://doi.org/10.3390/molecules28041768),
defines

\[
\text{acetic acid}+\text{1-pentanol}
\rightleftharpoons
\text{pentyl acetate}+\text{water}
\]

at \(318.15\ \mathrm K\) and approximately \(1\ \mathrm{bar}\). The paper
reports the four pure-component PC-SAFT parameter records, all six binary
interaction records, and \(K_a=43.99\) for the modeled standard-state
convention. Ascani obtained that constant from one Senina homogeneous
equilibrium composition using PC-SAFT activity coefficients. That composition
is calibration evidence, not a validation case.

Senina et al., “Chemical equilibria in the quaternary reactive mixtures and
liquid phase splitting: a system with n-amyl acetate synthesis reaction at
318.15 K and 101.3 kPa,” *Journal of Molecular Liquids* 345 (2022) 118246,
[doi:10.1016/j.molliq.2021.118246](https://doi.org/10.1016/j.molliq.2021.118246),
provides the corresponding experimental record. Its Table 3 contains nine
reactive liquid-liquid tie-lines with water-rich and organic-rich mole
fractions and a stated standard composition uncertainty of \(0.002\).
Senina used less than \(2\ \mathrm{wt}\%\) aqueous HCl as catalyst. Ascani
omitted the catalyst from its four-species calculation as a stated
low-concentration approximation, so the calculation is not a source-identical
model of the complete experimental mixture.

Ascani plots the PC-SAFT-predicted tie-lines but does not tabulate their
calculated endpoints. Senina tabulates experimental endpoints but does not
provide phase amounts for a single flash feed. A campaign can construct a
deterministic overall conserved state from a declared convex combination of
each endpoint pair, solve that state, and compare the calculated endpoints
with experiment. It must report this constructed-feed provenance and may not
describe the result as reproduction of unpublished Ascani endpoints.

Ascani reports the modeled condition as \(1\ \mathrm{bar}\), while Senina
reports \(101.3\ \mathrm{kPa}\). The Provider parameter record and Validation
case must retain the exact pressure attached to each comparison rather than
silently treating the two values as identical.

Together the papers provide:

- one reaction and an unambiguous conserved-component matrix;
- fixed \(T,P\), source compositions, and experimental reactive tie-lines;
- a complete published neutral PC-SAFT parameter hypothesis for the
  four-species approximation;
- both homogeneous and heterogeneous chemical-equilibrium data; and
- a direct one-liquid-versus-two-liquid discovery challenge.

The subject is neutral. Passing it establishes neither charged coordinates nor
electrolyte reactive LLE. It also requires a separately admitted installed
Provider artifact using the exact conventions of Ascani's four-species model
hypothesis. The campaign must preserve both the equilibrium-constant
calibration provenance and the omitted-catalyst model discrepancy.

### 2. Ascani hypothetical \(A+B\rightleftharpoons C\)

**Selected source-derived manufactured topology subject.**

The same 2023 paper supplies three hypothetical PC-SAFT pure-component
records, all three binary interactions, and several reaction constants. At
\(K_a=2.25\), the reported topology contains disconnected chemical-equilibrium
branches and two reactive tie-lines. This is useful for testing phase
generation and topology changes, but the reference endpoints are graphical
rather than tabulated. An independent analytic manufactured phase block
remains necessary for exact numerical assertions.

### 3. Tsanas, Stenby, and Yan (2017)

“Calculation of simultaneous chemical and phase equilibrium by the method of
Lagrange multipliers,” *Chemical Engineering Science* 174 (2017) 112--126,
[doi:10.1016/j.ces.2017.08.033](https://doi.org/10.1016/j.ces.2017.08.033),
contains tabulated reactive VLE and VLLE cases. It is valuable independent
algorithm evidence, but its vapor phases and ideal, activity-coefficient, or
cubic-EOS models are outside the first installed-PC-SAFT two-liquid slice.

### 4. Tsanas et al. electrolyte multiphase cases

The 2019 and 2022 papers cover aqueous electrolyte reactions with vapor,
liquid, and solid phases and include tabulated results. They are valuable
future electrolyte and phase-family extensions, but Pitzer/PR models,
phase-incidence rules, vapor, and solids exceed the initial same-Provider
two-liquid contract.

### 5. Koulocheris and Coatléven reactive cases

These papers provide useful cross-method and reactive-azeotrope comparisons.
Their model families and phase scopes are less aligned with the first
installed-Provider reactive-LLE capability than the Ascani--Senina pair.

## Evidence ladder

The first coupled implementation uses three independent levels:

1. an analytic manufactured two-liquid reaction problem with exact expected
   balances, phase states, multipliers, and reduced costs;
2. the Ascani hypothetical PC-SAFT system for source-derived topology and
   narrow-basin behavior; and
3. the Ascani--Senina pentyl-acetate system for an installed-Provider neutral
   reactive-LLE campaign.

An electrolyte reactive-LLE subject is a later capability slice. It must add a
source-complete charged-species parameter bundle, electrolyte standard-state
transforms, phase-specific electroneutrality evidence, and tabulated or
otherwise reproducible reactive phase states. Neutral evidence cannot be
relabelled as electrolyte evidence.
