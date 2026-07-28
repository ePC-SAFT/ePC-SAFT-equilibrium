# Carbonate Provider parameter sources

Date: 2026-07-28

## Decision

Use the source-coherent aqueous carbonate system from Ascani et al. (2022)
as the first non-ideal homogeneous chemical-equilibrium consumer:

`H2O`, `CO2(aq)`, `H3O+`, `OH-`, `HCO3-`, `CO3^2-`, and `Na+`.

Do not yet admit `CaHCO3+` or `CaCO3(aq)` as ePC-SAFT components. The
located sources provide free `Ca2+` parameters and PHREEQC complexation
constants, but not source-complete ePC-SAFT mechanical parameters for those
two complexes.

## Required Provider records

Each component needs its identity, charge, segment count, segment diameter,
and dispersion energy. Water and CO2 additionally need the published
association records. The model also needs a source-coherent dielectric/Born
formulation, all applicable unlike-pair `k_ij` and `l_ij` records, molar
masses used by the dielectric rule, and explicit temperature/pressure
applicability.

Reaction stoichiometry and source standard-state constants are not mechanical
ePC-SAFT component parameters. They remain chemical-equilibrium inputs and
must be transformed to the Provider reference basis without silently treating
concentration or PHREEQC constants as Provider-basis constants.

## Located source-complete carbonate core

Primary source:

- Moreno Ascani, Gabriele Sadowski, and Christoph Held, “Prediction of pH in
  multiphase multicomponent systems with ePC-SAFT advanced,” *Chemical
  Communications* 58 (2022), 8436–8439,
  [doi:10.1039/D2CC02943J](https://doi.org/10.1039/D2CC02943J).
- Official [supporting information](https://www.rsc.org/suppdata/d2/cc/d2cc02943j/d2cc02943j1.pdf),
  Tables S1–S5.

Published pure-component values:

| Species | m | sigma / Å | epsilon/k / K | Charge |
|---|---:|---:|---:|---:|
| H2O | 1.2047 | temperature correlation | 353.95 | 0 |
| CO2(aq) | 2.0729 | 2.7852 | 169.21 | 0 |
| Na+ | 1 | 2.8232 | 230.00 | +1 |
| HCO3- | 1 | 2.9296 | 70.00 | -1 |
| CO3^2- | 1 | 2.4422 | 249.26 | -2 |
| H3O+ (reported as H+) | 1 | 3.4656 | 500.00 | +1 |
| OH- | 1 | 2.0177 | 650.00 | -1 |

Water uses a 2B association scheme with association energy 2425.7 K and
association volume 0.04510. CO2 uses the source's induced-association
convention. The source treats aqueous `H+` as `H3O+` and uses
`H2O + CO2` as the carbonic-acid reaction participant; it does not introduce
a separate mechanical `H2CO3` component.

Relevant reported unlike-pair values include:

- water–ion `k_ij = 0`;
- `Na+`–`HCO3-`: `k_ij = -0.514`;
- `Na+`–`OH-`: `k_ij = 0.649`;
- water–CO2 temperature-correlation coefficients:
  `k_ij,a = 0.0122` and `k_ij,T = 0.0003016 K^-1`.

The exact equation convention for the temperature coefficient, every unlisted
pair, and the source pKa correlations must be resolved from the cited records
before constructing the bundle.
Unlisted values must not be invented.

The dielectric model is the ePC-SAFT advanced/original-Born formulation:
water relative permittivity `-105.2 ln(T/K) + 677.480`, ionic relative
permittivity 8, and ion Born radius tied to the source ion diameter.

The current Provider already implements the corresponding
`ionic-ascani-original-born` topology and generic neutral-reference SDK shape.
It does not yet contain this carbonate catalog bundle or its exact admitted
fingerprint/domain.

## Corroborating sources

Daniel Schick et al., “Predicting CO2 solubility in aqueous and organic
electrolyte solutions with ePC-SAFT advanced,” *Fluid Phase Equilibria* 567
(2023), 113714,
[doi:10.1016/j.fluid.2022.113714](https://doi.org/10.1016/j.fluid.2022.113714),
repeats the carbonate-ion values above and reports `Ca2+` as
`m = 1`, `sigma = 3.2648 Å`, and `epsilon/k = 1060 K`.

Martin Pabsch et al., “Modeling the CO2 Solubility in Aqueous Electrolyte
Solutions Using ePC-SAFT,” *Journal of Chemical & Engineering Data* (2020),
[doi:10.1021/acs.jced.0c00704](https://doi.org/10.1021/acs.jced.0c00704),
also provides free `Ca2+` and physical CO2/electrolyte interactions, but
explicitly neglects carbonate dissociation in its acidic scope. It therefore
does not close the reacting calcium-complex model.

M. Uyan et al., “Predicting CO2 solubility in aqueous
N-methyldiethanolamine solutions with ePC-SAFT,” *Fluid Phase Equilibria* 393
(2015), 91–100,
[doi:10.1016/j.fluid.2015.02.026](https://doi.org/10.1016/j.fluid.2015.02.026),
independently reports `H+`, `OH-`, `HCO3-`, and `CO3^2-` values consistent
with the Ascani carbonate set.

## Calcium-complex gap

PHREEQC 3.9.0 `wateq4f.dat` supplies reactions for:

- `Ca2+ + HCO3- = CaHCO3+`, log K = 1.106;
- `Ca2+ + CO3^2- = CaCO3(aq)`, log K = 3.224.

Those records are PHREEQC equilibrium constants, not ePC-SAFT mechanical
parameters. No located ePC-SAFT source provides the required component
diameter, dispersion energy, association treatment, and source-coherent pair
records for `CaHCO3+` and `CaCO3(aq)`. General ion-pair parameter-fitting
methods do not constitute published values for these species.

## Non-mixing rule

Do not combine the Ascani carbonate ions with the existing Held water
self-ionization bundle merely because species names overlap. For example, the
Held hydronium parameters and Ascani hydronium parameters are materially
different. Each installed bundle must preserve one complete source
formulation, reference convention, and applicability domain.
