# 1-Butanol–LiCl Parameterization Review

Originating commit: `2a98e8b689616a8159c494934bc2aede251829ff`

Recovered path: `tmp/research/butanol-licl-kij.md`

Original Git blob: `71955f4246405d5b10a62e45f3cc2a88fb62fde7`

Original content SHA-256:
`ece5505857310d182289888cfad5fbc3b7eaa37dfc502e048996b5ff544964a1`

Record type: historical package-local scientific review

Authority effect: `none`

This review records the parameter-coverage state observed at the originating
commit. It is not a Validation receipt, promotion, capability admission, or
predictive claim.

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

## Provider coverage at the originating commit

Provider supported both `k_ij` and `l_ij` and applied `l_ij` through the
cross-diameter and association combining rules. Its Ascani 2022
water/1-butanol electrolyte bundle used zero ion-neutral values and contained
Na+ and K+, not Li+. Its Figiel 2025 reference-electrolyte bundle included Li+
and Cl- interactions with water, methanol, and ethanol, but no 1-butanol/Li+
pair.

Therefore the published n-butanol/LiCl parameter set was not represented by
the bundle used for the archived campaign.

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

## Archived Figure 3 campaign disposition

The related `tmp/perdomo-figure3-final-contract/` campaign is intentionally not
restored here or relocated to Validation. Its metadata retained a parameter
fingerprint but not immutable Provider and Equilibrium artifact hashes, used
absolute local paths, and did not supply the complete source, transformation,
license/use-basis, and tolerance provenance required for durable Validation
evidence. Its own summary recorded 10 accepted splits, 20 indeterminate cases,
and 3 one-phase outcomes against a published Perdomo SAFT-gamma Mie model
curve, not experimental data.

Those outputs are therefore disposable exploratory diagnostics rather than a
Validation campaign. Their summarized scientific limitation is retained above;
the bulk results, traces, and plots remain excluded from `HEAD`. This
disposition has authority effect `none` and does not accept, reject, or promote
an electrolyte equilibrium capability.
