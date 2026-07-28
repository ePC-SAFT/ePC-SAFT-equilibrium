# Search for SAFT parameters for aqueous calcium–carbonate complexes

Date: 2026-07-28

## Question

Do published ePC-SAFT, PC-SAFT, or electrolyte-SAFT mechanical parameter sets
exist specifically for the aqueous species `CaHCO3+` and `CaCO3(aq)`? A
qualifying record must assign the complex itself a documented SAFT mechanical
description, such as segment number, segment diameter, dispersion energy, and
any association or unlike-interaction parameters. Parameters for free `Ca2+`,
`HCO3-`, or `CO3^2-`, and formation constants alone, do not qualify.

## Conclusion

No source-complete, published SAFT mechanical parameter set was found for
either `CaHCO3+` or `CaCO3(aq)`.

- Published ePC-SAFT carbonate models provide parameters for free ions and
  molecular species, but do not represent these two calcium complexes as
  independently parameterized SAFT components.
- Geochemical databases and public implementations provide formation
  constants for both complexes. Those records determine speciation but do not
  supply SAFT segment, dispersion, or association parameters.
- Held and Sadowski provide a transferable *procedure* for constructing and
  fitting an ePC-SAFT ion-pair species, but they neither fit nor publish a
  parameter set for either calcium–carbonate complex.
- A public NeqSim database row named `CaCO3` is not usable evidence: its
  metadata are internally inconsistent, it has no traceable calcium-carbonate
  SAFT source, and NeqSim's actual aqueous-complex calculation uses formation
  constants and Davies activity coefficients rather than this row as a
  fitted ePC-SAFT aqueous-complex model.

Consequently, neither complex can presently be admitted to a source-backed
Provider bundle as an explicit ePC-SAFT species. Doing so would require a new,
documented regression or an explicitly approved inherited-parameter
construction, followed by validation. Such a construction would be a new
model assumption, not recovery of an existing published parameter set.

## Evidence

### Published ePC-SAFT carbonate work

1. **Pabsch, Held, and Sadowski (2020), “Modeling the CO2 Solubility in
   Aqueous Electrolyte Solutions Using ePC-SAFT.”**
   [Publisher record and supporting information](https://pubs.acs.org/doi/10.1021/acs.jced.0c00704).
   The paper publishes free-ion parameters including `Ca2+` and `HCO3-`.
   It explicitly states that, for its pH range, carbonic dissociation is
   neglected, only physical interactions are considered, and complete salt
   dissociation is assumed. It therefore supplies neither `CaHCO3+` nor
   `CaCO3(aq)` as a mechanical ePC-SAFT component.

2. **Ascani et al. (2022), “Prediction of pH in multiphase multicomponent
   systems with ePC-SAFT advanced.”**
   [Publisher article](https://pubs.rsc.org/en/content/articlehtml/2022/cc/d2cc02943j)
   and [official supporting information](https://www.rsc.org/suppdata/d2/cc/d2cc02943j/d2cc02943j1.pdf).
   The source-complete carbonate acid–base set contains water, carbon dioxide,
   hydronium, hydroxide, bicarbonate, carbonate, and spectator ions. It does
   not include calcium or either requested calcium complex.

3. **Schick et al. (2023), “Predicting CO2 solubility in aqueous and organic
   electrolyte solutions with ePC-SAFT advanced.”**
   [Publisher record](https://doi.org/10.1016/j.fluid.2022.113714).
   The model publishes a mechanical parameter for free `Ca2+` and includes
   carbonic-acid dissociation through free `HCO3-` and `CO3^2-`. Full-text
   inspection found no `CaHCO3+` or `CaCO3(aq)` component and no mechanical
   parameter table for calcium–carbonate complexes.

4. **Held and Sadowski (2009), “Modeling aqueous electrolyte solutions.
   Part 2. Weak electrolytes.”**
   [Publisher record](https://doi.org/10.1016/j.fluid.2009.02.015).
   This is the closest applicable ePC-SAFT method source. For an explicit ion
   pair, it requires segment diameter, segment number, and dispersion energy.
   Diameter and segment number are inherited from the constituent ions through
   stated combining rules; ion-pair dispersion energy and the ion-pairing
   constant are fitted to mean ionic activity data. The fitted examples are
   weak acids, acetates, cadmium halides, zinc halides, and related systems—not
   calcium bicarbonate or aqueous calcium carbonate. Applying its rules to the
   requested complexes would therefore create a new transferred model and
   still leave complex-specific fitted quantities unresolved.

5. **Yang et al. (2024), “Ion Pairing in ePPC-SAFT for Aqueous and
   Mixed-Solvent Alkali Halide Solutions.”**
   [Publisher article and supporting information](https://pubs.acs.org/doi/10.1021/acs.iecr.4c01579).
   This work provides explicit ion-pair association frameworks and parameters
   for alkali halides. It does not parameterize calcium, bicarbonate, carbonate,
   `CaHCO3+`, or `CaCO3(aq)`. Its association rules are methodology, not a
   transferable published parameter record for the requested species.

### First-party databases and repositories

6. **PHREEQC version 3.**
   [Official repository](https://github.com/phreeqc-dev/phreeqc3).
   PHREEQC databases and tests contain `CaHCO3+` and `CaCO3(aq)` formation
   reactions and thermodynamic constants. They use aqueous activity models,
   not SAFT mechanical component parameters, so they are suitable as
   independent speciation evidence but cannot populate an ePC-SAFT component
   record.

7. **NeqSim.**
   [Official repository](https://github.com/equinor/neqsim),
   [`COMP.csv`](https://github.com/equinor/neqsim/blob/master/src/main/resources/data/COMP.csv),
   [`REACTIONDATA.csv`](https://github.com/equinor/neqsim/blob/master/src/main/resources/data/REACTIONDATA.csv),
   and
   [`ScalePredictionCalculator`](https://github.com/equinor/neqsim/blob/master/src/main/java/neqsim/pvtsimulation/flowassurance/ScalePredictionCalculator.java).
   `REACTIONDATA.csv` contains `CaHCO3ionpair` and `CaCO3ionpair` equilibrium
   records. The scale calculator identifies both aqueous complexes but computes
   them from association constants and Davies activity corrections. No
   `CaHCO3+` mechanical PC-SAFT component was found. `COMP.csv` has a row named
   `CaCO3`, but it carries the methane CAS number `74-82-8`, a molar mass of
   `124.0`, and an implausible SAFT diameter field of `481.8`; no cited fit or
   applicability to `CaCO3(aq)` is attached. It is not a defensible published
   parameter set.

8. Exact code searches were also run over the public GitHub corpus and targeted
   at the first-party repositories for **FeOs**, **Clapeyron.jl**,
   **ThermoPack**, **teqp**, and public PC-SAFT implementations. Searches for
   `CaHCO3`, `CaHCO3+`, `CaCO3(aq)`, and combinations with `PC-SAFT`,
   `ePC-SAFT`, `ePPC-SAFT`, and `SAFT` found geochemical reaction databases,
   speciation examples, mineral-scale code, and molecular-simulation work, but
   no explicit SAFT mechanical parameter set for either aqueous complex.
   GitHub rate limiting interrupted the final repeated repository queries; the
   unrestricted corpus results obtained before that limit likewise contained
   no qualifying SAFT record.

### Other exact searches

Publisher, repository, dissertation, and general scholarly searches used these
exact species/model combinations:

- `"CaHCO3+" "PC-SAFT"`, `"CaHCO3+" "ePC-SAFT"`,
  `"calcium bicarbonate" "electrolyte SAFT"`, and
  `"CaHCO3" "SAFT" thermodynamic model`
- `"CaCO3(aq)" "PC-SAFT"`, `"calcium carbonate" ePC-SAFT`,
  `"CaCO3" "electrolyte PC-SAFT"`, and
  `"calcium" "carbonate" "ePPC-SAFT"`
- `PC-SAFT parameter database calcium bicarbonate`,
  `ePC-SAFT parameter database CaCO3 ion pair`, and
  `SAFT gamma Mie CaCO3 aqueous ion pair parameters`

The searches covered ACS Publications, Elsevier/ScienceDirect, RSC, Wiley,
TU Dortmund's Eldorado repository, institutional dissertation indexes,
GitHub, and the local Zotero full-text collection. Results about solid calcite,
free `Ca2+`, calcium salts treated as fully dissociated, PHREEQC/Pitzer/SIT
formation records, or molecular-dynamics binding free energies were excluded
because they do not provide the requested SAFT mechanical component
parameters.

## Scientific boundary

The negative finding is a literature-search result, not proof that no private
or unindexed parameterization exists. The defensible implementation choices
are presently:

1. validate homogeneous carbonate chemistry without explicit calcium complexes;
2. retain `CaHCO3+` and `CaCO3(aq)` only in an independent geochemical
   cross-model comparison; or
3. authorize a new Regression-owned fit with a declared source dataset,
   inherited-parameter rules, uncertainty assessment, and an immutable
   Provider artifact.

## Current Regression readiness

The latest local Regression work inspected was branch
`codex/general-parameter-regression` at commit `b3c1fb1`. It is substantially
ahead of public `main` (`11c664c`) and is not currently published on the
remote. The checkout also contains user-owned uncommitted test/data work, so
this assessment is read-only and does not treat the branch as admitted.

The branch has useful foundations:

- one typed Ceres engine with exact installed-Provider derivatives;
- source, units, bounds, transforms, partitions, rank, conditioning, active
  bounds, and confirmation-start diagnostics;
- scalar fit paths for `segment_count`, `segment_diameter`,
  `dispersion_energy_over_k`, Born diameter, solvation factor, relative
  permittivity, selected dielectric parameters, and selected `k_ij`/`l_ij`;
- one bounded joint five-coordinate pure 2B saturation fit.

It is not yet capable of fitting the requested calcium complexes:

- the reusable `RegressionProblem` explicitly accepts exactly one shared
  parameter;
- component `m`, `sigma`, and `epsilon/k` fits use neutral pure-saturation
  pressure/density observations, which are unavailable and physically
  inappropriate for dissolved `CaHCO3+` and `CaCO3(aq)`;
- direct electrolyte observations currently activate one advertised
  coordinate, not a new species' coupled mechanical parameter block;
- there is no reactive-speciation observation contract or exact
  Equilibrium implicit parameter-sensitivity path for fitting mechanical
  parameters against calcium-carbonate speciation;
- the special five-parameter association route is a neutral pure 2B
  saturation workflow, not a general aqueous ion-pair fitter.

Focused current-checkout verification passed 51 tests in 220.23 seconds,
covering the general contracts, scalar pure parameters, bounded association
derivatives, aqueous `k_ij`, Born diameters, and solvent/ionic permittivity.
This confirms the documented mechanics but does not broaden their domain.

The shortest scientifically defensible extension would follow Held and
Sadowski's ion-pair construction: fix inherited quantities under an explicit
source rule, identify only the unresolved complex-specific coordinates, and
fit them against source-bound aqueous activity/speciation data. That requires
new Provider exact derivatives and a new Regression observation contract,
likely consuming Equilibrium's exact reactive-state sensitivities. Formation
constants alone cannot identify SAFT mechanical parameters.
