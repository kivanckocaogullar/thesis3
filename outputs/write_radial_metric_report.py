from pathlib import Path
import json
p=Path('outputs/radial_quadratic_robust_20261001_161903');d=json.loads((p/'metric.json').read_text());ck=json.loads((p/'independent_checks.json').read_text())
report='''# Degree-two state-dependent contraction metric: nominal-voltage radial active transient

SUCCESS on the FULL original rectangle, for the specified two-state quasi-static Q-V model with fixed published network coefficients. One common quadratic metric with lambda=0.1 s^-1 passes both independent dense numerical checks and a directed interval-arithmetic uniform check. No smaller domain, new fault, parameter change, or independent pointwise metrics were used.

## Domain and physical-time vector field

Gamma in degrees: [24.423430419712,33.858410380369].
Nu: [.996174514620,1.005178107392].
Actual radian endpoints used and validated: [0.42626927545571497,0.5909407406288647]. This is the floating-point representation of the requested angular bounds; metric.json defines the numerical model precisely.

Nominal voltage ONLY, U=.98039215686274506 on the 20.4-kV base. External SCR5/XR10 and the physical additional .0858989549982-mH reactor are unchanged; reconstructed Zeq=.03020391991902692+j*.35000000000000009 pu. Vset=1 Qstar=0 Imax=1.2 umax=1.3 H=3.4 K=20 Kd=5 D=25; no virtual impedance. Pstar=1 is held fixed at every state. Fault-on data were NOT used to construct this vector field.

f1=100*pi*(nu-1),
f2=((1-Pe(gamma))/nu-25*(nu-1))/6.8,
J=[[0,100*pi],[-Pe_gamma/(6.8*nu),(-(1-Pe)/nu^2-25)/6.8]].

Electrical closure uses q=(Eref exp(j*(gamma+arg(Vth)))-Vth)/Zeq and Iconv=1.2*q/|q| on the active branch; Vc=Vth+Zeq Iconv, Vgrid=vgE Vc+vgU U, Igrid=igE Vc+igU U, Pe+jQgrid=Vgrid conj(Igrid), u_required=Vc+(Zconverter+.2)Iconv.

The quasi-static Q-V root solves Eref+.1Qgrid-1=0. Its TOTAL derivative is E_gamma=-.1Q_gamma/(1+.1Q_E), Pe_gamma=P_gamma+P_E E_gamma. No polynomial replacement of the trigonometric/root equations was made.

Recomputed current switch: 24.373430419712 deg. Recomputed fixed-Pstar saddle: 64.880551979056 deg, excluded from the domain. The instantaneous Jacobian has a positive real eigenvalue throughout the rectangle, numerically ranging from 1.3022749338 to 1.66171858275 s^-1. Thus a constant metric cannot satisfy a strict continuous-time contraction inequality there. The new metric uses its full Mdot contribution.

## Actual common polynomial

Define scaled coordinates (gamma in radians):

z1=(gamma-0.50860500804228981)/0.082335732586574867,
z2=(nu-1.000676311006)/0.0045017963860000676.

M=M0+M1*z1+M2*z2+M3*z1^2+M4*z1*z2+M5*z2^2.
Each Mi=[[a_i,b_i],[b_i,d_i]]. The entries below reproduce metric.json's full decimal coefficients:

| i | a_i | b_i | d_i |
|---|---:|---:|---:|
'''
for i,m in enumerate(d['coefficients']):report+=f'| {i} | {m[0][0]:.17g} | {m[0][1]:.17g} | {m[1][1]:.17g} |\n'
report+='''
This is M in the PHYSICAL state coordinates [gamma,nu], not a metric that still requires a coordinate congruence transformation. Only its polynomial arguments are scaled. The derivatives use the chain rule:

M_gamma=(M1+2*z1*M3+z2*M4)/0.082335732586574867,
M_nu=(M2+z1*M4+2*z2*M5)/0.0045017963860000676,
Mdot=M_gamma*f1+M_nu*f2.

trace(M(center))=1 to double-precision arithmetic. Directional finite-difference validation of Mdot gives maximum entrywise error 1.333e-8 on independent off-equilibrium test points. Use the full precision JSON/MAT coefficients when reproducing the result.

## Uniform and numerical bounds

Let Qc=J'M+MJ+Mdot+2*lambda*M with lambda=.1 s^-1.

| Quantity | Independent 301x301 numerical grid | Validated uniform bound |
|---|---:|---:|
| Minimum eigenvalue of M | 0.000100885494400996 | > 1e-5 |
| Maximum eigenvalue of M | 9.93193031701739 | < 10 |
| Worst maximum eigenvalue of Qc | -0.00107414116608405 | < -1e-5 |
| Contraction rate | .1 s^-1 | .1 s^-1 |
| Raw-current excess over 1.2 | >= .002451136739 pu | See interval_validation.json |
| Voltage-command headroom | >= .051276946768 pu | See interval_validation.json |
| Limited-current magnitude | 1.2 pu | Exactly 1.2 by the verified active radial formula |

Numerical extrema are grid extrema, not asserted exact optimizers. The rigorous bounds are deliberately weaker. The metric remains anisotropic: the dense global max/min ratio is approximately 9.84e4. Explicit positive lower and upper bounds prevent collapse, but Euclidean-distance conclusions must account for this conditioning.

## Optimization and validation method

A single common set of 17 independent coefficients was optimized (six symmetric matrices, one exact trace constraint). For fixed lambda the constraints are LMIs, solved with MATLAB LMI Lab feasp. A 9x9 grid was adaptively augmented with failing validation points. Complete Mdot was included in every constraint. A first rate-bisection experiment found a densely passing numerical metric at lambda=3.570375 s^-1, but its minimum eigenvalue was only about 2.814e-6. That high-rate metric is retained separately as a numerical experiment; it is NOT the metric certified here, and no global maximum rate is claimed. The final .1 s^-1 solve strengthened the positivity floor from 1e-6 to 1e-5 and provided wider validation margins. Some solver/grid iterates failed off-grid and were rejected/refined.

The uniform check uses mpmath directed interval arithmetic at 30 decimal digits, not just denser point sampling:
1. Cover the entire rectangle by adaptive interval boxes, with no untested gaps.
2. For each angular box, prove opposite signs of h at E=.98 and E=1.05 and positive h_E over that bracket. This establishes an existing unique root in the bracket for every angle in the box. Interval Newton contracts the root enclosure; the active limiter condition is then verified at the enclosed root.
3. Evaluate actual sin/cos, square roots, radial derivatives and implicit Q-V derivatives using intervals. The validation uses exact interval 100*pi for omega_b.
4. Check M-1e-5 I and 10 I-M positive definite by interval 2x2 Sylvester tests.
5. Enclose Qc directly or using a centered mean-value form. Interval first derivatives of Qc include analytical second electrical derivatives and the second implicit E derivative; no finite-difference derivative is used in the proof. Prove -Qc-1e-5 I positive definite.
6. Check raw-current activation and voltage-command feasibility over every box.

The interval verifier's complete result, leaf count and conservative physical margins are in interval_validation.json. The centered derivative implementation was independently checked by second-derivative finite differences (error below 1.7e-13) and 44 enclosure/corner checks. These sanity tests supplement, rather than replace, interval arithmetic.

Scope of validation: real-valued dynamics defined by the frozen DECIMAL network/metric constants in metric.json, with exact interval pi in the physical-time angular equation. It does not propagate component-identification uncertainty or prove exact equality between floating-point network reduction and unknown physical hardware. The chosen smooth Q-V root is certified within E=.98..1.05; global uniqueness outside that bracket is not claimed. This is a computer-assisted uniform inequality check, not a proof-assistant formalization or an SOS certificate.

## Existing trajectory containment and limitations

No new trajectory/fault was simulated. The previously saved primary nominal-voltage postfault trajectory enters this domain at clearing t=.1 s. Of its 265 sampled active points, 264 lie inside the rectangle; the remaining point is in the deliberately excluded switching layer. Interpolation of the saved trajectory places lower-bound exit at t=.428127071599 s. Actual saved limiter release is t=.428985231116 s, leaving about .00085816 s outside this interior before release. The interpolation is numerical containment evidence, not a validated ODE enclosure.

The uniform inequality certifies differential contraction while trajectories and the relevant comparison curves remain within this smooth active domain. The rectangle is not asserted forward-invariant. It does not cover the .05-degree switching layer, the inactive recovery region, or the fault-on dynamics. No metric matching, hybrid transition certificate, or final recovery guarantee has been established.

The underlying analytical model still uses instantaneous Q-V response and ideal current-reference tracking, with PLL frequency approximated by nominal and PI/filter/electrical transient states omitted. This certificate does not transfer automatically to the full dynamic Simulink controller.

Methodological reference: Aylward, Parrilo and Slotine (2008), Stability and Robustness Analysis of Nonlinear Systems via Contraction Metrics and SOS Programming, https://web.mit.edu/nsl/www/preprints/SOS_metrics.pdf . Their state-dependent contraction inequality is used; exact GFM dynamics remain trigonometric, piecewise and implicit. The polynomial object here is the metric, not a falsely assumed polynomial vector field.

Files: metric.json and metric.mat; console.txt; independent_checks.json; interval_math_checks.json; interval_validation.json. The search, robust solve and interval checker are separate new scripts in outputs. Frozen baseline scripts, original metrics and Simulink remain unchanged.
'''
(p/'report.md').write_text(report,encoding='utf-8')
print(p/'report.md')
