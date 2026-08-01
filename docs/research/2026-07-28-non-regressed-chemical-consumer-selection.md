# Non-regressed chemical-equilibrium consumer selection

Date: 2026-07-28

## Decision

Use aqueous water self-ionization at 298.15 K and 1 bar as the first real
installed-Provider reference-transformation consumer:

```text
2 H2O <=> H3O+ + OH-
```

This case is non-MEA, nonideal, charge balanced, and source complete when its
two source responsibilities remain explicit:

- Held, Cameretti, and Sadowski (2008) supply the ePC-SAFT mechanical model
  and ion parameters.
- IAPWS R11-07(2019), whose formulation is unchanged from the 2007 release,
  supplies the water ionization constant. The source standard state must be
  transformed explicitly into the Provider Helmholtz basis.

It requires no new parameter regression. It is one-reaction Provider
integration evidence, not a replacement for the existing multi-reaction
Belov--Aristova algorithm and trace benchmark. Those two tests establish
different claims and should be retained together.

## Installed Provider evidence

Data commit `c096285415d4d3198b9d00fc75af48b837dd1305` contains the canonical
`held-cameretti-sadowski-2008/1` packet. Equilibrium pins its complete packet
fingerprint and passes the packet's parameter bundle to the installed EOS with
`Parameters.from_bundle(...)`. The latest EOS SDK therefore remains the
runtime owner while reusable scientific data remains in Data; tests fail
closed when the packet or any required component, association, interaction,
or formulation record changes.

The catalog fixes:

- ordered species `water`, `hydronium-cation`, `hydroxide-anion`;
- charges `(0, +1, -1)`;
- the Held water segment, dispersion, temperature-dependent diameter, and 2B
  association records;
- published hydronium and hydroxide segment diameters and dispersion energies;
- water--ion dispersion rules and neglected ion--ion dispersion;
- the fully dissociated Debye--Huckel topology; and
- the concentration-independent solvent dielectric formulation.

The immutable Provider wheel at
`artifacts/equilibrium-perdomo-figure1a-v1/efc4c6c/epcsaft-0.1.0.dev0-cp313-cp313-linux_x86_64.whl`
has SHA-256
`9b05f4d202424f3db570dfaa80210e5aa54c9ee3586a4cbef4309f5881648d08`.
It is the first retained artifact that combines this catalog and
neutral-reference callback with the electrolyte phase-value SDK tail required
by current Equilibrium `main`. A black-box check using only the installed wheel
constructed the three-component model, evaluated its liquid reference at
298.15 K and 1 bar, and obtained a finite reference convergence error of
`1.0628611260443677e-05`. The wheel exports the native SDK capsule and its
neutral-reference callback.

The SDK neutral basis contains the pure-water row and the charge-neutral
hydronium-plus-hydroxide row. Therefore the self-ionization reference
contraction is reconstructed without assigning absolute single-ion standards:

```text
reference(H3O+ + OH-) - 2 reference(H2O)
```

## Source records

Primary mechanical source:

- C. Held, L. F. Cameretti, and G. Sadowski, “Modeling aqueous electrolyte
  solutions: Part 1. Fully dissociated electrolytes,” *Fluid Phase
  Equilibria* 270 (2008), 87--96,
  [doi:10.1016/j.fluid.2008.06.010](https://doi.org/10.1016/j.fluid.2008.06.010).

Primary reaction-constant source:

- International Association for the Properties of Water and Steam,
  [IAPWS R11-07(2019), Revised Release on the Ionization Constant of
  H2O](https://iapws.org/relguide/Ionization.pdf).

IAPWS records that the 2019 update changed explanatory text and references,
not the 2007 formulation. The 2024 release readjusts parameters and should not
silently replace the frozen benchmark. The implementation must bind the exact
release, equation, density input, logarithm convention, and molality standard
state used to calculate the 298.15 K value.

## Required Equilibrium test

Restore only the generic neutral-reference transformation previously proved by
Equilibrium PR #56, distilled into the existing singular chemical solver. The
test must:

1. load the exact installed Provider wheel and catalog through public
   Provider interfaces;
2. bind temperature, pressure, species order, charges, Provider fingerprint,
   source identity, source standard state, and Helmholtz basis identity;
3. transform the IAPWS source constant to the Provider basis;
4. solve through the existing fixed-`T,P` homogeneous owner;
5. independently certify water balance, hydrogen/oxygen conservation,
   electroneutrality, reaction affinity, positivity, pressure, Provider
   domain/packing, KKT stationarity, and reduced curvature;
6. preserve species-permutation, reaction-sign, and conservation-gauge
   invariance; and
7. fail closed for source-unit, fingerprint, order, temperature, pressure,
   basis-rank, or reference-convergence mismatch.

Acceptance establishes one source-bound local homogeneous liquid result. It
does not establish predictive chemistry, a general carbonate model,
sensitivities, coupled phase equilibrium, or globality.

## Why carbonate is later

The Ascani carbonate acid--base system remains the better later
multi-reaction nonideal validation. Its published mechanical parameters can be
transcribed without fitting, but the exact carbonate catalog, all source
interactions/correlations, and its installed artifact are not yet present.
Calcium complexes remain excluded because no source-complete SAFT mechanical
parameters were found for `CaHCO3+` or `CaCO3(aq)`.

The separate Provider catalog
`ascani-2022-case-study-2-electrolyte-lle` is a water/1-butanol/Na+/K+/Cl-
phase-equilibrium bundle. It is not the carbonate reacting system and must not
be relabeled as one.

## Regression readiness

Regression `origin/main` has a general typed parameter framework, while the
local `codex/multi-parameter-regression-core` branch adds ordered fitted
parameter vectors. Its fit-ready observations cover selected pure saturation,
fixed-composition VLE, MIAC, solvation-Gibbs, and dielectric targets.

It does not currently accept reactive-speciation observations or consume an
Equilibrium value/sensitivity contract. Its own architecture forbids copying
reaction equations, implementing a second equilibrium solver, or finite
differencing a black-box equilibrium solve. Cross-association identities,
temperature-dependent interaction coefficients, and several dielectric
correlation coefficients also remain not ready.

Consequently Regression cannot infer missing calcium-complex mechanical
parameters from PHREEQC formation constants, and it is neither required nor
useful for the selected water-self-ionization checkpoint. Regression becomes
relevant only if a later application deliberately selects unknown parameters
and supplies identifiable source-bound observations plus exact Provider and
Equilibrium sensitivities.
