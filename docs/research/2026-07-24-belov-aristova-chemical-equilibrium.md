# Belov-Aristova chemical-equilibrium benchmark contract

Date: 2026-07-24
Scope: primary-source extraction only; no code, tests, Reaktoro, commit, or push.

Sources read: the supplied PDF and its supplied `.mmd` extraction for G.V. Belov
and N.M. Aristova, “Calculation of Complex Chemical Equilibrium Using
Optimization Package Ipopt,” *International Journal of Thermodynamics* 26(4),
pp. 077-083 (2023), plus the authors'
[`Heterogeneous-Equilibrium`](https://github.com/gvbelov/Heterogeneous-Equilibrium)
repository at commit `c74b87545a3418262e8e38a7c7a2e31e1b12966a`. Page
references below are PDF pages followed by the printed page in parentheses. The
PDF rendering was checked for Eqs. (1)-(26), Figure 5, and Table 1. `[Fact]` is
stated by the paper or authors' repository; `[Inference]` is a consequence of
those facts; `[Blocker]` is missing evidence needed for a reproducible
source-faithful benchmark.

## Bottom line

[Fact] The paper gives a generic constrained Gibbs/Helmholtz minimization
template and two illustrative calculations, not fully specified benchmark
datasets. The companion repository nevertheless supplies a small, exact
`test1.dat` carbon/oxygen input that can support a derived homogeneous gas-only
algorithm benchmark. It does not publish an expected composition for that
restricted problem. The printed Julia snippets assign into `phase_moles` after
allocating `phase_mols`; correcting that name is necessary before they can run.

[Inference] Neither published application case nor the derived gas-only test can
be represented as an ePC-SAFT validation case by the accepted installed
Provider. The nuclear case needs a large source-specific chemical database and
phase inventory; the rocket case needs RP-1/combustion and nozzle
thermochemistry; the gas-only test uses ideal-gas standard Gibbs values. Belov
test 1 is therefore useful for the generic numerical formulation and
trace-species certification, not for Provider source-reference integration.

## 1. General executable formulation

### Fixed-\textit{T},\textit{p} Gibbs route

[Fact] The unknown vector is substance amounts `x_i` in mol. For `N` substances
and `m` chemical elements, the paper defines the formula matrix `nu_ji` as the
number of atoms of element `j` in substance `i`, and `b_j` as the supplied
element amounts. The fixed-\textit{T},\textit{p} problem is

```text
minimize       G(T,p,x),   x in R^N
subject to     T = const, p = const
               sum_i nu_ji*x_i = b_j,  j = 1..m
               x_i >= 0,              i = 1..N
```

This is Eq. (1), PDF p. 2 (078). The paper later calls the input matrix `A` and
the code enforces `A' * x == b` (Figure 1, PDF p. 3 (079)); therefore the
concrete storage orientation must be fixed when a benchmark is authored.

[Fact] For `N_c` pure condensed phases and `N_m` solution phases, the general
Gibbs model is Eqs. (2)-(3): pure-phase standard contributions plus
`x_i*(g_i + ln(a_i))` in each solution, with `g = G/(R*T)` and
`g_i = G_i/(R*T)`. The executable special case is Eq. (4), ideal gas + ideal
solution + zero condensed-phase volume:

```text
g(T,p,x) = sum_(i=1..Nc) x_i*g_i
         + sum_(j=1..Nm) [ sum_(i in I_j) x_i*(g_i + ln(x_i))
                           - y_j*ln(y_j) ]
y_j = sum_(i in I_j) x_i
```

For condensed species, Eq. (5) gives
`g_i = [H_i°(T) - T*S_i°(T)]/(R*T)`; for gas species Eq. (6) adds
`ln(p/p°)`. These are PDF p. 2 (078), Eqs. (2)-(6).

### Fixed-\textit{T},\textit{V} Helmholtz route

[Fact] The alternative problem is Eq. (7), PDF p. 2 (078): minimize `F(T,V,x)`
with fixed `T,V`, the same element balances, and `x_i >= 0`. The general form
is Eq. (8); after declaring mixture `j=1` to be gas, Eqs. (9)-(10) use the gas
term `x_i*(f_i + ln(x_i))` and condensed-solution terms
`x_i*(f_i° + ln(x_i)) - y_j*ln(y_j)`, with `f = F/(R*T)`. Eqs. (11)-(12) define

```text
f_i° = [H_i°(T) - T*S_i°(T)]/(R*T)                 condensed
f_i  = [H_i°(T) - T*S_i°(T)]/(R*T) + ln[R*T/(p°*V)] - 1   gas
```

[Fact] In the printed `calc_Helmholtz` code (Figure 1, PDF p. 3 (079)), `y[1]`
has no gas-phase amount constraint; only `j=2..ns` gets `sum(x_i)=y[j]`.
Consequently a returned `phase_mols[1]` is not defined by the shown model.
The code also uses `f[i]` in all terms although Eq. (10) distinguishes `f_i°`
from `f_i`; the intended array convention is therefore unresolved. [Blocker]

### Function inputs, bounds, coordinates, and outputs

[Fact] Both routines accept `(m, k, ns, nc, g_or_f, jx, A, b)`, where `k` is
the number of substances, `ns` the number of solution phases, `nc` the number
of pure condensed phases, `jx` stores the inclusive index interval for each
solution, `A` is the element matrix, and `b` the element feed. The component
ordering is: pure condensed substances first, gaseous substances second, then
condensed-solution components (Section 4, PDF p. 3 (079)).

[Fact] The JuMP variables are the direct amounts `x[1:k]` and phase amounts
`y[1:ns]`, each with lower bound `0` and start value `1.e-3` (Figure 1). The
Gibbs routine minimizes Eq. (4), enforces `sum(x[jx[1,j]:jx[2,j]]) == y[j]`
for every `j=1..ns`, and enforces `A' * x == b`. The Helmholtz routine uses the
gas interval `j=1`, condensed-solution intervals `j=2..ns`, and enforces only
the latter `y` relations plus `A' * x == b`. Both return `x`, `y`, and the
JuMP shadow prices `lam` for the element constraints.

[Inference] Because `log(x_i)` and `log(y_j)` occur while the declared bounds
permit zero, a faithful implementation needs a defined interior/trace-floor
policy. The paper supplies none. The source also says equilibrium amounts can
range from about `10^0` to about `10^-100` mol, below machine zero in practice
(Introduction, PDF p. 1 (077)); this is a numerical warning, not a declared
bound.

### Unknown temperature or pressure

[Fact] If temperature is not assigned, the outer equation is Eq. (13), PDF
p. 4 (080): `Z - sum_i x_i(T)*z_i(T) = 0`, where `Z` may be enthalpy, entropy,
or internal energy and `z_i` is the corresponding partial-molar property. The
bracketed Roots.jl path is exactly:

```text
T = findrootS(find_S, tmin, tmax, eps)
find_zero(f, [tmin,tmax], atol=eps, Order1())
```

The alternative is Newton's Eq. (14), with Eqs. (15)-(18) for `H,p`, `S,p`,
`U,V`, and `S,V`, respectively. If `T,S` are known, pressure uses Eq. (19).
The paper does not publish values for `tmin`, `tmax`, `eps`, initial `T` or
`p`, maximum iterations, safeguarding, or convergence tests. [Blocker]

## 2. Companion-repository test 1: executable homogeneous restriction

### Exact source identity and input

[Fact] The authors' repository is MIT licensed. At commit
`c74b87545a3418262e8e38a7c7a2e31e1b12966a`, Git blob
`5490caab3f7a3f8cdea0d1cb883cce4323902657`, `test1.dat` has SHA-256
`2f24904e5db64e96aa2ee96c53bd80863c9518b6c2ddd914e266181fe7459f84`.
It defines `T=2000 K`, `p=0.1 MPa`, and one mole of `C2O`, hence element totals
`C=2 mol` and `O=1 mol`. The complete source candidate list contains graphite,
diamond, and these eight gas species in this order:

| Species | `g_i°/(RT)` | C atoms | O atoms |
|---|---:|---:|---:|
| `O(g)` | -7.072468897041624 | 0 | 1 |
| `O2(g)` | -28.765987200068952 | 0 | 2 |
| `O3(g)` | -26.51300119486301 | 0 | 3 |
| `C(g)` | 21.44955454258845 | 1 | 0 |
| `C2(g)` | 21.26761422175228 | 2 | 0 |
| `CO(g)` | -34.35125913287487 | 1 | 1 |
| `CO2(g)` | -55.36542486137457 | 1 | 2 |
| `C2O(g)` | -17.10224518423043 | 2 | 1 |

[Inference] Removing graphite and diamond is a deliberate homogeneous
gas-phase restriction, not a reproduction of the source's full heterogeneous
answer. It supplies a private chemical-speciation benchmark while leaving phase
selection and phase equilibrium out of scope.

### Reference conversion and independent oracle

[Fact] The source convention is ideal gas with `p°=101325 Pa`. The private
manufactured ideal kernel uses `ln(n_i/V)` with a `1 mol/m3` concentration
reference. At fixed `T`, the exact per-species conversion is

```text
g_i,manufactured = g_i,source + ln(R*T/p°).
```

[Inference] For element multipliers `lambda_C` and `lambda_O`, the independent
ideal-gas KKT equations are

```text
n_i = n_total*(p°/p)
      * exp[-g_i,source - a_C,i*lambda_C - a_O,i*lambda_O],
sum_i a_C,i*n_i = 2,
sum_i a_O,i*n_i = 1,
sum_i n_i = n_total.
```

Solving those three equations with 80-decimal arithmetic, independently of the
Equilibrium implementation, gives

```text
lambda_C = -10.0776517753450153959143507476192439107332024990741
lambda_O =  44.8480568890243540815059574241188924984510845101041
n_total  =   1.49962490524610185240315148798558911129621759471888 mol
```

and amounts, in source gas-species order,

```text
O    5.96993828070211394402844371512880493560268644194452e-17
O2   5.24905958060066503006510156103611049344002112807903e-27
O3   1.83805491064889610693933315321367549351472938908179e-47
C    1.74962363295254877517350452410666696344037406210204e-5
C2   4.99607409054311182829527221270435265660344809626476e-1
CO   9.99232314211335204005658153304511473719499834407715e-1
CO2  4.45388559141274630810276404330965885463563769000304e-11
C2O  7.67685699587024466704103014174643928539777026177891e-4
```

[Fact, local diagnostic] The unmodified private log-coordinate/Ipopt kernel
conserves both elements and returns optimizer success, but pins multiple trace
species near `1e-12 mol` even when the declared trace floor is reduced to
`1e-50 mol`. Its independently recomputed reaction-affinity infinity norm
remains approximately `37` to `47`, so the existing acceptance logic correctly
rejects the result. [Inference] This case exposes transformed-coordinate
stationarity degeneracy and is suitable for the generic trace/boundary
readiness leaf; it must not be made green by relaxing the physical-affinity
gate.

## 3. Published benchmark A: nuclear-reactor emergency state

### Case definition

[Fact] Figure 5 and the surrounding text (PDF p. 5 (081)) define `T=2000 K`
and `p=1 bar = 0.1 MPa`. The printed initial feed is:

```text
1014.5 UO2 + 0.096 Np + 2.754 Pu + 0.824 Ce + 0.215 Y
+ 0.138% Te + 0.332 La + 1.442 Zr + 0.389 Ba + 0.899 Ru
+ 1.15% Mo + 0.265 Pr + 0.421 Sr + 0.0385 I2 + 0.859 Nd
+ 0.043% Nb + 0.0064 Am + 0.745 Cs + 0.166 Rh + 0.006 Sb
+ 0.025 Eu + 3725 H2O + 3725 H2
```

The percent signs on Te, Mo, and Nb are present in the PDF and `.mmd`; their
mass/mole meaning and their denominator are not explained. [Blocker] The
paper gives no reaction stoichiometry and no concrete `A`/`nu` matrix. [Fact]
It describes a heterogeneous system with 13 displayed pure condensed phases
and one displayed gas mixture; the listed concentrations are dominant species
only, not a complete candidate or equilibrium inventory.

### Published outputs

[Fact] Figure 5 prints the following properties exactly as shown: `p=0.1 MPa`,
`t=2000 K`, `v=1240.58 cub.m`, `s=1972.06 kJ/K`,
`h=-1.38049e+06 kJ`, `u=-1.50455e+06 kJ`, `Cp=428.554 kJ/K`, and
`Cv=366.525 kJ/K`. The figure labels the phase amounts as mol. Pure condensed
phase amounts are:

```text
Ru(c) 0.898991       Rh(c) 0.164611       Mo(c) 0.310964
Nb2O5(c) 0.0166375   Y2O3(c) 0.107237     Ce2O3(c) 0.411914
Pr2O3(c) 0.13243     Nd2O3(c) 0.429459     Eu2O3(c) 0.0123781
UO2(c) 1006.11       Np(c) 0.0831693       PuO2(c) 2.75149
BaZrO3(c) 0.333271
```

The dominant gas-mixture amounts are:

```text
H2(g) 3741.23         H2O(g) 3698.8          H(g) 8.59741
UO2OH(g) 8.00739      OH(g) 0.843116         CsOH(g) 0.593765
Sr(OH)2(g) 0.405405   LaO2(g) 0.331998       MoO3(g) 0.301278
MoO2(OH)2(g) 0.281665 MoOOH(g) 0.145607      Cs(g) 0.134998
Te(g) 0.1149          U2O6(g) 0.11274
```

The source itself warns that direct comparison with Pelton [33] is impossible
because that database is unavailable; it presents Figure 5 as an illustration,
not a reproduced reference result (PDF p. 5 (081)). The paper does not publish
the objective value, element residuals, KKT residuals, phase-rule calculation,
shadow prices, omitted trace species, or complete phase list. [Blocker]

### Numerical/model interpretation

[Fact] Section 7 says all cases used pure-substance properties from
IVTANTHERMO [29]. The algorithm checks mass balance, Kuhn-Tucker conditions,
and the Gibbs phase rule; scaling was sometimes done by recalculating element
content to `1 kg`. The most difficult described system had `22` chemical
elements, about `130` condensed substances, about `135` gaseous substances,
`2` solutions, and several dozen single-component phases (PDF p. 5 (081)).
[Inference] Figure 5 does not publish enough information to establish whether
all those candidate phases were loaded for this particular run; only the
displayed dominant phase set is executable from the article.

[Fact] The reported preparation time is “about several seconds” including
package/database loading and compilation. Measured calculation time was from
several hundredths to tenths of a second when temperature was given, and could
increase to several seconds when temperature was unknown (PDF p. 5 (081)).

## 4. Published benchmark B: LOX/RP-1 rocket performance

[Fact] The case is liquid oxygen plus kerosene RP-1, with oxidizer/fuel ratio
`3`, combustion-chamber pressure `200 bar`, and nozzle-exit pressure `0.4 bar`
(PDF pp. 5-6 (081-082)). The chamber equilibrium is computed at prescribed
pressure and enthalpy; throat and exit equilibrium are computed at prescribed
pressure and entropy. The individual RP-1 surrogate species, absolute feed
amount, enthalpy/entropy values, throat pressure, nozzle geometry/area ratio,
and phase/topology assumptions are not published. [Blocker]

The only numerical outputs are Table 1 (PDF p. 6 (082)):

| Program | `T_chamber` (K) | `I_vac` (m/s) | `C*` (m/s) |
|---|---:|---:|---:|
| CEA | 3867.2 | 3574.8 | 1781.1 |
| RENGINE | 3871.7 | 3575.3 | 1781.5 |

[Fact] These are performance metrics, not a published equilibrium composition.
The paper attributes discrepancies to different thermodynamic databases. The
outer unknown-property contract is Eqs. (13)-(19), but no root tolerances,
initial guesses, iteration limits, CEA settings, RENGINE revision, or complete
thermochemical input are supplied. [Blocker]

## 5. Standard state, activities, and numerical solver

[Fact] The paper's generalized activity convention is `a_i` in Eqs. (2)-(3);
the runnable special case substitutes ideal mole-fraction activities in each
solution and ideal-gas pressure dependence in Eq. (6). Condensed phases have
zero volume in Eq. (4). This is not an ePC-SAFT activity-coefficient or
Helmholtz-provider convention.

[Fact] IVTANTHERMO/reference-book standard-state data use `p°=1 atm=101325 Pa`
and the paper gives the polynomial variables `X=T/10000` in Eqs. (20)-(22),
including `Delta_f H°(298.15)`. It also records alternative NIST equations
(23)-(24), with `t=T/1000` and `p°=1 bar`, and NASA equations (25)-(26).
The coefficients for any of the reactor or rocket species are absent. [Blocker]

[Fact] The implementation is Julia + JuMP + Ipopt. JuMP automatic
differentiation supplies derivatives; Ipopt is described only as an
interior-point solver. The only explicit optimization settings are direct
amount coordinates, `x>=0`, `y>=0`, starts `1.e-3`, and the shown equality
constraints (Figure 1, PDF p. 3 (079)). No Ipopt version, linear solver,
`tol`, `max_iter`, print level, scaling policy, derivative-check setting,
termination acceptance rule, or multistart/global-search procedure is given.
The unknown-temperature bracket uses Roots.jl `Order1()` and `atol=tol`.

## 6. Installed ePC-SAFT Provider disposition

[Fact, local] The repository pins `epcsaft==0.1.0.dev0` (`pyproject.toml`). Its
accepted Equilibrium capability is only the fixed-temperature pure-component
local saturation boundary for provider-approved methane, ethane, and propane
(`README.md`, `ARCHITECTURE.yaml`, `docs/phase-equilibrium.md`). The local
provider catalog probe found Gross-Sadowski methane/ethane, Gross-Sadowski
propane, and a bounded Figiel electrolyte catalog containing water, methanol,
ethanol, Li+, Na+, K+, Cl-, and Br-; the latter is private/non-production for
Equilibrium and does not contain nuclear or RP-1 species.

[Inference] The nuclear Figure 5 case cannot be represented by the installed
Provider: its explicit species, element inventory, 13 condensed products, gas
species, candidate phase list, and source thermochemistry are outside the
catalog and outside the accepted capability. The rocket case likewise cannot be
represented: RP-1 chemistry, combustion species, chamber `H,p`/nozzle `S,p`
workflow, and performance calculation are not Provider inputs or an accepted
Equilibrium route. At most, the paper's abstract ideal-model shell could be
reimplemented with externally supplied `g/f`, `A`, `b`, and phase intervals;
that would be a separate ideal chemical-equilibrium implementation, not an
ePC-SAFT benchmark.

## 7. Complete blocker/prerequisite list

1. For the reactor and rocket application cases, exact source
   species/candidate-phase inventories, ordered component IDs, formulas, and
   the full concrete element matrix `A`/`nu`; neither published case supplies
   these. For the reactor case, resolve the `%Te`, `%Mo`, and `%Nb` semantics
   and all feed units. The companion `test1.dat` gas restriction is exempt from
   this blocker.
2. A complete source thermodynamic database for every admitted species,
   including `H°(T)`, `S°(T)`, `Delta_f H°(298.15)`, coefficient tables, source
   version, and one unambiguous standard pressure. The cited nuclear database
   is unavailable, and the rocket surrogate species are unspecified.
3. A declared activity/phase model. The paper's runnable equations are
   ideal-gas/ideal-solution/zero-condensed-volume; an ePC-SAFT reproduction
   would instead need Provider parameters, association/ionic model choices,
   binary interactions, domains, and a source-faithful standard-state bridge.
4. If reactions are to be represented rather than only element balances, the
   independent reaction stoichiometry, reaction rank, equilibrium constants,
   source standard states, activity scales, and exact transformation into the
   Provider basis are missing. The paper supplies none of these.
5. Complete phase-incidence data and a phase-discovery/stability policy. Figure
   5 reports dominant outputs only; the paper does not establish that the
   displayed phase set is complete.
6. A corrected, runnable implementation of Figure 1, including the
   `phase_moles`/`phase_mols` typo, the Helmholtz gas `y[1]` ambiguity, a
   positive-log boundary policy, and the intended `f_i` versus `f_i°` array
   convention.
7. Reproducible solver metadata: Julia/JuMP/Ipopt/Roots.jl versions, Ipopt
   options and linear solver, exact starts/brackets, tolerances, iteration
   limits, scaling, failure handling, and the code revision at repository [36].
8. Published acceptance evidence: objective values, complete `x/y` outputs,
   element and phase residuals, KKT/complementarity residuals, shadow-price
   convention, rank/conditioning evidence, and an independent reference output.
   Table 1 alone is not an equilibrium-composition oracle.
9. For the rocket case specifically: RP-1 surrogate composition, absolute
   propellant feed, chamber enthalpy, throat/exit entropy states, nozzle
   geometry or expansion ratio, CEA and RENGINE versions/settings, and the
   route from equilibrium composition to `I_vac` and `C*`.
10. A new, separately admitted Equilibrium capability and an exact installed
    Provider artifact containing the needed species/model fingerprints would be
    required before either case could become an ePC-SAFT validation subject.
