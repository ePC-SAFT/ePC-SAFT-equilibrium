# 1-Butanol–LiCl interaction-parameter findings

## Direct ePC-SAFT evidence

Gabriel Grundl's dissertation, Table II-5, reports the ePC-SAFT binary
parameters used to model liquid-liquid equilibria at 298.15 K and 1 bar:

| Pair | k_ij | l_ij | Provenance |
| --- | ---: | ---: | --- |
| n-butanol/Li+ | 0 (table cell is blank) | -0.08 | fitted in that work |
| n-butanol/Cl- | 0.22 | 0.245 | prior reference 412 |

The accompanying text says the ion/n-butanol parameters were fitted to ternary
water/n-butanol/salt LLE data. The related 2016 Fluid Phase Equilibria paper
(DOI 10.1016/j.fluid.2016.05.001) reports successful ePC-SAFT modeling of the
water/1-butanol salt effect, including LiCl salting-in above 6 mol/kg water.

## Current Provider coverage

Provider supports both `k_ij` and `l_ij` and applies `l_ij` through the
cross-diameter and association combining rules. Its current Ascani 2022
water/1-butanol electrolyte bundle uses zero ion-neutral values and contains
Na+ and K+, not Li+. Its Figiel 2025 reference-electrolyte bundle includes Li+
and Cl- interactions with water, methanol, and ethanol, but no 1-butanol/Li+
pair.

Therefore the published n-butanol/LiCl parameter set is not represented by the
current bundle used for the campaign.

## Different-model evidence

Perdomo et al. (2025) use SAFT-gamma Mie group parameters, including unlike
energies for Li+–CH2OH and Cl-–CH2OH. Those values are relevant qualitative
evidence but cannot be inserted as ePC-SAFT molecular `k_ij` or `l_ij` values:
the model, units, and combining rules differ.

## Consequence

The organic-phase salt discrepancy should first be retested with the published
ePC-SAFT n-butanol/Li+ and n-butanol/Cl- parameters, together with the matching
water/ion, Li+/Cl-, pure-component, dielectric, and solvation conventions.
Changing HELD2 before this model-level correction would conflate equilibrium
search behavior with an incorrect thermodynamic parameterization.
