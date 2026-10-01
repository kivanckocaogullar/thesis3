# Q-priority common-metric study: SCR5, Xeq = 0.350 pu

The common constant metric succeeds on finite nonlinear domains. A degree-two state-dependent metric was therefore not necessary for this first enlargement. These are numerical enclosure results, not a rigorously validated uniform certificate or an EMT-controller certificate.

All frozen parameters and the network reconstruction are unchanged. Current-reference strategy is q-priority in the rotating VSM dq frame; the voltage-command limit remains 1.3 pu, constant-ratio. H=3.4, K=20, Kd=5, D=25, Imax=1.2. No current-loop/voltage-loop PI states are added. The raw reference is still the ideal q(gamma) proxy, not the actual PI reference.

## Before / after

| Branch | Original Theta: equilibrium screening, deg | New equilibrium family, deg | Nonlinear gamma domain, deg | nu domain | Metric |
|---|---|---|---|---|---|
| Inactive | (20.167968, 24.421813) | [0.040000, 24.381813345295] | [0.020000, 24.401813345295] | [0.99,1.01] | One common constant M, degree 0 |
| Active | (45.459623897403, 47.521506348664) | [29.457028775316, 47.481506348664] | [29.437028775316, 47.501506348664] | [0.99,1.01] | A separate common constant M, degree 0 |

The original transformation is unchanged: beta=-Pe_U/6.8, c=-Pe_gammaU/6.8, Theta=[-beta,0;-c,-beta]. Its entries depend on the reference equilibrium; its old interval is not a common-metric nonlinear certificate. Zero-tolerance screening boundaries are reproduced; the old 1e-9 tolerance moves the active lower bound to approximately 45.459623980416 deg.

| Branch | min eig(M) | Worst eig(Q), Q=J'M+MJ+2lambda M | lambda used | Numerical rate bound | Minimum current margin, pu | Minimum voltage margin, pu |
|---|---:|---:|---:|---:|---:|---:|
| Inactive | 0.0516552793828 | -0.166427815167 | 1e-5 | 0.0879651389904 | 0.000967803166945 | 0.0552534637565 |
| Active | 0.105647858147 | -1.75458423097e-5 | 1e-5 | 2.06327063139e-5 | 0 (roundoff -4.44e-16) | 0.000143574036696 |

Rates are per normalized analytical time. The rate bounds inherit the numerical-enclosure limitation; they are not formally validated bounds. Metrics have trace(M)=1. Full precision matrices and network data are saved in results.mat. Approximate matrices are:

Inactive: [[0.8849,0.2299],[0.2299,0.1151]].
Active: [[0.8244,0.2242],[0.2242,0.1756]].

Use full precision MAT matrices, not the rounded display, to reproduce inequalities.

## Boundary meanings

- Original inactive lower: original Sym(F) loses negative definiteness. Upper: current activation, 24.421813345295 deg.
- Original active lower: original Sym(F) loses negative definiteness. Upper: required voltage reaches 1.3 pu, 47.521506348664 deg.
- Physical active stable/voltage-feasible envelope remains (29.400102759326,47.521506348664) deg. The lower endpoint is Pe_gamma=0. No metric can certify a positive contraction rate at that zero-eigenvalue equilibrium.
- New inactive lower: chosen positive-angle search floor 0 deg plus two 0.02-deg guards. It is NOT a physical stability boundary or a claimed maximal interval.
- New inactive upper: activation minus 0.02-deg domain guard and 0.02-deg equilibrium-interior guard.
- New active lower: numerical common-metric feasibility search reaches approximately 29.417028775316 deg for the prescribed nu range and lambda; add a 0.02-deg slack guard and 0.02-deg equilibrium-interior guard. This is a sufficient-method boundary, not an intrinsic physical limit or globally optimal metric boundary.
- New active upper: hard voltage limit minus the same two guards. No voltage-limit relaxation was used.

The 0.02-degree guards deliberately trade a small interval loss for positive slack and a finite neighbourhood at every reported equilibrium. The fixed rate 1e-5 and nu halfwidth 0.01 are explicitly chosen search specifications, not claimed optimal choices. Other domain/rate tradeoffs were not exhausted.

## Nonlinear calculation and fixed Pstar

For each equilibrium gamma0, Pstar=Pe(gamma0), nu0=1. Pstar then stays constant:

f=[nu-1; ((Pstar-Pe(gamma))/nu -25*(nu-1))/6.8].

J=[[0,1],[-Pe_gamma/(6.8*nu),(-(Pstar-Pe)/nu^2-25)/6.8]].

Each common metric is tested for all independent gamma and nu in its rectangle and all Pstar in a conservative Pe range covering the equilibrium family. Pstar is never substituted by Pe(gamma) off equilibrium. Analytic J was checked by fixed-Pstar finite differences with maximum error below 3.3e-10.

J depends affinely on two scalar entries a,b. Numerical extrema of Pe and Pe_gamma enclose these entries in a rectangle; checking the four Jacobian vertices then covers that entire matrix polytope by convexity. The optimization maximizes the minimum negative LMI eigenvalue subject to trace(M)=1 and min eig(M)>=1e-6, at lambda=1e-5. fmincon is used for this small convex feasibility problem; lower domain boundaries are bisected. A second initial metric and independent 100001-angle electrical check are included.

Electrical extrema use 4001 samples, refinement of detected extrema and small floating pads. These are NOT directed-rounding interval bounds, and unsampled extrema are not mathematically excluded. Therefore the polytope argument is uniform CONDITIONAL ON the electrical enclosures; its present end-to-end status is numerical nonlinear-domain evidence, not a formal certificate. No SOS certificate is claimed. Exact trigonometric and square-root limiter dynamics were not falsely replaced by polynomials.

The full rectangle need not be invariant. Conditional on its inequalities, invariant local ellipsoids about each fixed equilibrium follow from V=(x-xd)'M(x-xd), using a sublevel contained in the convex rectangle. A common conservative choice for every family member is:

- Inactive rho=2.59354118354e-8, coordinate projections +/-0.0141421356237 deg and +/-0.000684478614084 in nu.
- Active rho=3.27808444651e-8, coordinate projections +/-0.0141421356237 deg and +/-0.000534802516741 in nu.

These projections describe ellipsoids, not full rectangular perturbation boxes. They are small sufficient neighbourhoods, not full basins or tested disturbance recovery limits. Interior equilibria can permit larger sublevels; that maximization was not performed.

For a constant M, Mdot=0 exactly. A quadratic metric with omitted Mdot was not used. No bridge between the inactive and active domains or contraction through a switching boundary is asserted. No fault or EMT simulation was run.

## Interpretation and next step

The previous small passing interval was substantially limited by the original metric choice. A common degree-zero metric already supplies much wider numerical nonlinear contraction evidence on the q-priority active branch. It does not change Pe_gamma, local stability, voltage feasibility, or the physical system. A rigorous follow-up should validate the electrical extrema/branch/voltage bounds with directed interval arithmetic and then validate the vertex LMIs; quadratic SOS optimization is not yet needed to demonstrate this enlargement.

Methodological reference: Aylward, Parrilo and Slotine, Stability and Robustness Analysis of Nonlinear Systems via Contraction Metrics and SOS Programming (2008), https://web.mit.edu/nsl/www/preprints/SOS_metrics.pdf . The present implementation uses their contraction inequality with constant M; it is not an implementation of a polynomial SOS certificate.
