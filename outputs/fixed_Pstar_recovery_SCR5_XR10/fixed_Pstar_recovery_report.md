# Fixed-Pstar numerical recovery domain

Pstar=1 for EVERY trajectory, including limiter-active transient states. Same saved SCR5/XR10 radial benchmark; H=3.4, K=20, Kd=5, D=25, Imax=1.2. No EMT, metric optimization, or physical model changes. This is not a contraction certificate.

## Fixed-power equilibria

| gamma deg | limiter active | Iconv pu | Vc pu | Pe_gamma pu/rad | eigenvalues | classification |
|---|---|---|---|---|---|---|
| 13.810758941138 | 0 | 1.015307154 | 0.997477278 | 4.082037814 | -0.171259194103, -3.50521139413 | LOCALLY_STABLE |
| 56.301096020206 | 1 | 1.200000000 | 0.894388940 | -0.345863660 | 0.0137828751915, -3.69025346343 | SADDLE |

Nominal attractor: 13.810758941138 deg. Current switches: [-16.34092603758;16.34092603758] deg. Roots searched over [-180,180], separately across nonsmooth branches, including stationary-point cuts. Periodic copies are not nominal recovery.

## Recovery definition and validity

Angle error <=0.1 deg, |nu-1|<=1e-4, |gamma_dot|<=1e-4 and |nu_dot|<=1e-4 continuously for 10 normalized analytical time units. Recovery entry/exit uses adaptive ODE events on the maximum normalized residual; the tighter check divides all tolerances by 10. recovery_time is the start of the qualifying interval; recovery_confirmed_time is its end. Post-limiter recovery time ends at confirmation.

Each trajectory first runs to 100, with unresolved/tighter-ball and branch-departure cases extended to 200. Near-saddle boundary calculations additionally use 2000 because slow saddle departure can dominate finite-horizon classification. PHASE_SLIP_OR_ESCAPE means an observed |gamma-gamma_eq|>=180 deg, a numerical departure from the nominal angle branch. Data/episodes continue through the horizon to show the subsequent branch; a separate full-turn event at 360 deg is an observation bailout. nu<=0.1 leaves the analytical domain. These are stopping/observation criteria, not theorems. Voltage-invalid trajectories are reported separately, never automatically called unstable. No voltage limiter is added.

The vectorized ode45 mesh uses RelTol=2e-9, AbsTol=2e-11, MaxStep=0.25 (0.5 for extended studies), independent event functions per trajectory, no state reset/hysteresis. Saved electrical maxima are evaluated on the adaptive mesh plus dense 0.1-time samples and event times; they are numerical extrema, not certified global bounds. Full representative records include phasors, vector fields and episodes. Open episodes are right-censored: exit time NaN, dwell measured through observation end.

## Domain and fractions

Primary scan: gamma=0:0.25:90, nu=0.95:0.002:1.05. Adaptive extensions and the 0.01 deg angle slice are separate.

- Boundary intersects nu edges 0.950/1.050: expanded to 0.900/1.100.
- Boundary intersects nu edges 0.900/1.100: expanded to 0.800/1.200.
- Boundary intersects nu edges 0.800/1.200: expanded to 0.500/1.500.
- The separatrix remains open at the finite nu limits; this is a sampled window, not a globally exhausted basin.

primary: 18411 states; recovered 56.656347%; limiter-active episode 81.869535%; limited AND recovered 38.525881%. These are point-count fractions, not phase-space area estimates.

all sampled rows: 35244 states; recovered 56.591193%; limiter-active episode 82.890705%; limited AND recovered 39.481898%. These are point-count fractions, not phase-space area estimates.

| Classification | Count |
|---|---:|
| NOT_RECOVERED_WITHIN_HORIZON | 4195 |
| PHASE_SLIP_OR_ESCAPE | 11104 |
| RECOVERED | 19945 |

Recovered categories: inactive/never limited 6030; inactive/later limited 405; initially limited 13510.

## Angle-jump and boundary results

At nu0=1, T=200 recovery bracket: [51.001453400, 51.001539230] deg. Extended T=2000 bracket: [56.301069260, 56.301155090] deg; Delta gamma bracket [42.490310319, 42.490396149] deg. The lower endpoint recovers; the upper does not satisfy the finite-horizon recovery definition. Brackets are refined below 1e-4 deg; this bracket resolution is not a claim of an exact global stability boundary.

Positive limiter-active saddle: 56.301096020206 deg. Extended nu0=1 bracket midpoint minus saddle = 1.6154782e-05 deg. The 2-D refined boundary tilts with initial nu; the saddle local stable eigenvector is overlaid. This supports its organizing role but no global stable manifold has been rigorously computed.

Recovered-side extended near-boundary trajectory: initial gamma 56.301069260, limiter entry at t=0.000000000 and gamma=56.301069260; max gamma 56.301069260; total limiter dwell 1051.907180523; last exit 1051.907180523; recovery confirmation 1087.162423648. All times are normalized analytical time.

## Penetration and dwell

Among the finite-horizon scan recovered cases: farthest positive gamma 58.500000000 deg, from (gamma0,nu0)=(58.500000000,0.500000000); penetration beyond positive switch 42.159073962 deg. Largest absolute angle excursion 44.689241059 deg. Longest total limiter dwell 164.741065252, from (51.250000000,0.984000000); longest individual episode 164.741065252. These are sampled maxima, not universal allowable dwell times. Near-boundary extended trajectories can dwell much longer.

## Recovery mechanism

Between the positive switch and the positive saddle, Pe>1 even though Pe_gamma<0. At nu=1 this gives nu_dot<0, followed by nu<1 and gamma_dot<0: the state can exit current limitation and return to nominal. For initially nonunit nu, damping also matters: the sign of 1-Pe ALONE is not the sign of nu_dot. Above the saddle Pe<1 creates positive accelerating tendency at nu=1.

A transient at gamma=20 deg is not a new equilibrium with Pstar=Pe(20). Local instability of the Pe(gamma0)-parameterized equilibrium family does not prohibit recovery of transient states at fixed Pstar=1. The empirical recovery boundary separates outcomes at a stated horizon; it is not a contraction-certified boundary.

## Numerical checks

- savedClosureMaxError: 0
- representativeClassificationAgreement: 1
- maxRepresentativeFinalAngleDifference_deg: 7.95807864051e-13
- primaryVsTightBallDifferent: 409
- extendedBoundaryTightDisagreements: 0
- maxRepresentativeDwellDifference: 8.17708473733e-06
- maxRepresentativeExcursionDifference_deg: 6.3346522694e-08
- phaseDepartureThreshold_deg: 180

Boundary endpoints also carry tight-ball disagreement flags. Any flagged finite-horizon boundary depends on the settling tolerance and must not be presented as a physical separatrix. Representative classifications were repeated with RelTol=2e-11, AbsTol=2e-13 and MaxStep=0.125.

## Additional physical and branch checks

Maximum voltage command across the scan = 1.242006458247 pu < 1.3 pu; voltage-invalid count = 0; analytical-domain-invalid count = 0.

Largest tested recovered-side initial angle among the extended boundary brackets = 64.052009583 deg. Longest tested recovered-side boundary dwell = 1470.567014010 normalized analytical time. These extend the T<=200 scan results; neither is a universal bound.

| Representative | Initial gamma / nu | Pe initially | nu_dot initially | Minimum nu | Final unwrapped gamma | Classification |
|---|---|---|---|---|---|---|
| Inactive recovery | 14.810758941 / 1.000000000 | 1.071115344 | -0.010458139 | 0.997438578 | 13.810758941 | RECOVERED |
| Initially limited recovery | 20.000000000 / 1.000000000 | 1.168172974 | -0.024731320 | 0.992785449 | 13.810758941 | RECOVERED |
| Near boundary: recovered | 56.301069260 / 1.000000000 | 1.000000162 | -0.000000024 | 0.992785221 | 13.810758941 | RECOVERED |
| Near boundary: not recovered | 56.301155090 / 1.000000000 | 0.999999643 | 0.000000052 | 1.000000000 | 373.810758941 | PHASE_SLIP_OR_ESCAPE |
| Longest limited dwell (scan) | 51.250000000 / 0.984000000 | 1.029544609 | 0.054408087 | 0.984000000 | 13.810758941 | RECOVERED |
| Farthest penetration (scan) | 58.500000000 / 0.500000000 | 0.986549557 | 1.842191307 | 0.500000000 | 13.810758941 | RECOVERED |
| Inactive start, later limited | 16.000000000 / 1.500000000 | 1.155330512 | -1.853463776 | 0.992785435 | 13.810758941 | RECOVERED |
| Farthest extended-boundary recovery | 64.052009583 / 0.500000000 | 0.951039342 | 1.852635488 | 0.500000000 | 13.810758941 | RECOVERED |
| Longest extended-boundary dwell | 56.859998703 / 0.964000000 | 0.996614745 | 0.132869364 | 0.964000000 | 13.810758941 | RECOVERED |

The outer nu0=1 example approaches the next angular copy near 373.810758941 deg, not the original 13.810758941 deg equilibrium. Its branch departure is classified separately even though the reduced model later settles in that other branch.
