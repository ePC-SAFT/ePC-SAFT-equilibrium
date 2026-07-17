# Neutral Two-Phase TP Flash Local Candidate

Stage: locally implemented candidate; independent validation and promotion pending

Authority effect: none

## Capability boundary

The candidate adds one operation:

```python
two_phase_flash(
    model: epcsaft.EPCSAFT,
    temperature: pint.Quantity,
    pressure: pint.Quantity,
    overall_mole_fractions: Sequence[float],
) -> TwoPhaseFlashResult
```

It admits only the reviewed Gross--Sadowski methane/ethane binary fingerprint
`sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286`
inside the rectangular domain audited from May et al. (2015), Table 5:
203.22--243.61 K, 2.124--6.885 MPa, and methane feed fraction
0.4661--0.66705. Inputs must be finite positive Pint quantities and a positive
two-component feed summing to one within `1e-12`; invalid feeds are not
normalized.

The phase count is fixed at two. A source-domain input can still have no
distinct local two-phase solution. The candidate carries no phase-discovery,
TPD, global-stability, critical, continuation, association, ionic, reactive,
or multiphase certificate.

## Ownership and provider boundary

The installed provider wheel owns the dimensionless extensive fixed-temperature
Helmholtz function and its exact gradient and Hessian in component amounts and
phase volume. Equilibrium owns phase variables, material balances, the
Helmholtz-plus-`PV` transform, Ipopt, initialization, confirmation, canonical
density ordering, local acceptance, and typed results.

The candidate was built only against the public header installed by provider
wheel SHA-256
`17c6735ee117469b13b76b6a669d0b4430071c8eaebdf1e405baefc9adcb838b`.
The reviewed public-header SHA-256 is
`b1dc4a666799603ea87fe515ff08226455a11e52908b7109ebe1841369cb92df`;
provider code commit is `0cd75b96f1dddba024c324df500493bd1564bf5c` and
candidate-records commit is `25e597289b044ff3e5af70b3572279293555b9f9`.
No provider source is compiled or imported, and the equilibrium extension
links no provider implementation symbol.

`checked_sdk()` remains the ABI validation owner in `cpp/src/module.cpp`.
Only the provider context/call ownership seam is shared in
`cpp/src/provider.cpp`. The Python capsule stays alive across the GIL-released
solve.

## Formulation and exact derivatives

For liquid and vapor candidate phases, equilibrium solves

```text
min  [Phi_L(T,n_L,V_L) + P V_L/(R T)]
   + [Phi_V(T,n_V,V_V) + P V_V/(R T)]

subject to n_L,i + n_V,i = z_i,  i in {CH4,C2H6}.
```

The normalized feed is exactly one mole. Solver coordinates are

```text
(n_L,CH4, n_L,C2H6, log V_L, n_V,CH4, n_V,C2H6, log V_V).
```

Component amounts have explicit positive bounds and phase volumes use log
coordinates. The two balances are linear and have zero constraint Hessian.
The objective gradient and Hessian use only provider first and second
derivatives, including the exact log-volume chain rules. Independent centered
directional differences test both the objective-gradient and
gradient-Hessian products; production uses no finite differences.

## Attempts and local acceptance

Two deterministic component-agnostic enrichment directions are tried with a
`0.5/0.5` phase split, liquid density `20000 mol/m3`, and an ideal-gas vapor
volume. The first locally acceptable solution receives one materially
perturbed confirmation solve. There is no fallback solver, continuation, or
paper-identity seed.

A successful result separately requires:

- acceptable Ipopt termination and constraint violation at most `1e-10`;
- one agreeing confirmation solve within `1e-7`;
- material-balance residual at most `1e-10`;
- relative specified-pressure stationarity at most `1e-8` in both phases;
- component chemical-potential difference at most `1e-8`;
- KKT stationarity at most `1e-8`;
- relative phase-density distance at least `1e-3`; and
- finite positive states with all bounds inactive by the named margin.

The native owner retains Ipopt equality and lower/upper bound multipliers. Its
independent diagnostic evaluates
`grad(f) + J^T lambda - z_L + z_U`, so the reported KKT residual is not a
copy of solver status. A collapsed or bound-active state is rejected with
structured `FlashError` diagnostics.

## Audited source evidence

The source contract is validation commit
`73a37f5935e919a34d1e4fa3af285951d6fac8e7`, CSV SHA-256
`5cd1e74925a3c6504f5106dcf911f2cae2d6e99a5133fccc20454d8991bdbc7f`,
metadata SHA-256
`d43433e93b354e01f96d330c760818a24b775026461ce795e45774cfb11ac94e`,
and tolerance-payload SHA-256
`ad744526678355be6ca47cf27ab9ff7ae66b7661c27e36ffe259c5b6295f1016`.
The controlling PDF SHA-256 is
`53fd1bdd55dc6807ec76cf88626438d8dfceb3ec09149d4405ea36cfbe6b842a`;
the locator is Table 5, page 3612 (PDF page 7). The frozen comparison is
absolute composition error at most `3*uc`, with no model-accuracy floor.
This source is evidence input, not a runtime dependency.

The required package anchor `may2015-ch4-c2h6-001` at 243.58 K and
3.949 MPa uses `x_CH4=0.3099`, `y_CH4=0.6664`, and midpoint feed
`z_CH4=0.48815`. The candidate returns:

```text
x_CH4 = 0.3025223259589743   |dx| = 0.0073776740410257 <= 0.0111
y_CH4 = 0.6703563353120439   |dy| = 0.0039563353120439 <= 0.0114
liquid fraction = 0.49534934421235416
material balance max abs = 1.1102230246251565e-16
pressure stationarity max relative = 1.452287772722744e-12
chemical-potential max abs = 6.757794324130373e-11
KKT stationarity max abs = 7.993605777301127e-15
confirmation max difference = 3.597815403046881e-15
```

A package-side, non-authoritative characterization of all 17 frozen rows found
16 locally certified states and 12 rows within both composition allowances.
Rows 002, 009, 010, and 011 missed at least one allowance. Row 012 converged
to indistinguishable phase densities and was correctly rejected as collapsed.
These are model/data observations, not automatic solver defects, and they do
not constitute the later validation-owned predictive campaign. No tolerance,
source row, seed, or scope was changed in response.

## Artifact evidence

The local candidate wheel
`epcsaft_equilibrium-0.1.0.dev0-cp313-cp313-linux_x86_64.whl` has SHA-256
`2259eacee30792b3b473b278923eeef8dbad057dcac330d6ec122f79bdc03c19`.
Package tests cover provider transport, capsule lifetime, exact directional
derivatives, input and rejection paths, retained multipliers, local
certification, the source anchor, and the unchanged accepted pure route.

The all-row concerns prevent a predictive-validation or promotion claim. This
document records local candidate evidence only and transfers no authority.
