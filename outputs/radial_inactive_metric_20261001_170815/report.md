# Limiter-inactive recovery contraction analysis

SUCCESS: an independently optimized degree-two metric is uniformly interval-validated on the full rectangle, for the specified two-state quasi-static Q-V model. The original active metric is preserved.

| Quantity | Result |
|---|---|
| Gamma bounds | [16.426345885756, 24.323430419712] deg |
| Nu bounds | [0.99652713910517932, 1.0020607312196841] |
| Certified rate | lambda = 1 s^-1 |
| Uniform metric bounds | 1e-5 I < M < 10 I |
| Uniform contraction inequality | Qc < -1e-5 I |
| Dense metric eigenvalue extrema | 0.00106345455772, 9.90507398282 |
| Dense worst eigenvalue of Qc | -0.000656986669102 |
| Uniform current margin below 1.2 pu | > 0.0008986 pu |
| Uniform command-voltage margin below 1.3 pu | > 0.048374 pu |
| Dense current / voltage margins | 0.002428071661 / 0.051636547519 pu |
| Nominal equilibrium | 20.633904875622 deg, nu=1; contained |
| Saved primary inactive samples contained | 6060 / 6060 |
| Sampled maximum condition number | 1441.45 (physical coordinates) |

Uniform bounds are conservative certified bounds, not exact extrema. Dense values are numerical observations. Conditioning is anisotropic but finite; gamma uses radians and nu uses pu. Rate optimality is not claimed. Lambda=.1 first passed dense testing. A separate lambda=1 search passed adaptive refinement and 301-by-301 verification; its first coarse-grid solution failed off-grid, and violating nodes were added before acceptance.

## Equations and assumptions

All SCR5/XR10, physical Xeq=.350, H=3.4, K=20, Kd=5, D=25, Imax=1.2, umax=1.3, Vset=1, Qstar=0 settings are unchanged. Pstar=1 remains FIXED throughout this nonlinear domain. Nominal grid voltage only; no fault-on characteristic is used.

Inactive electrical closure:
q=(E exp(j*(gamma+arg(Vth)))-Vth)/Zeq; Iconv=q;
Vc=Vth+Zeq*q;
Vgrid=vgE*Vc+vgU*U; Igrid=igE*Vc+igU*U;
Pe+jQgrid=Vgrid*conj(Igrid);
u_required=Vc+(Zconverter+.2)*Iconv.

Solve E+.1*Qgrid(gamma,E)-1=0. Then E_gamma=-.1*Q_gamma/(1+.1*Q_E) and Pe_gamma_total=P_gamma+P_E*E_gamma.

Physical-time dynamics:
f=[100*pi*(nu-1); ((1-Pe)/nu-25*(nu-1))/6.8].
J=[[0,100*pi],[-Pe_gamma_total/(6.8*nu), (-(1-Pe)/nu^2-25)/6.8]].

The Q-V closure is quasi-static; it does not retain the actual Q filter or PI controller states. The ideal current-reference proxy, nominal-frequency electrical phasors and algebraic command-voltage feasibility are the existing benchmark approximations, not full EMT validation.

## Full polynomial metric

M=sum(Mi*bi), b=[1,z1,z2,z1^2,z1*z2,z2^2].
z=(x-center)./scale; x=[gamma(rad),nu]. The metric acts on physical-state variations; only its polynomial arguments are scaled. trace(M0)=1 fixes scaling.

center = [0.35560888299079796, 0.9992939351624317]

scale = [0.06891506321292276, 0.0027667960572523853]

M0:
```text
0.0019393908480797441  0.0048196741442657299
0.0048196741442657299  0.9980606091519203
```

M1:
```text
-0.0013341268114822709  0.035957661929351015
0.035957661929351015  -0.79210448926199917
```

M2:
```text
-0.005458055915879746  -0.013102188599069205
-0.013102188599069205  -1.0595536958512497
```

M3:
```text
0.0034356508337600902  -0.0015081223526455226
-0.0015081223526455226  4.8318438397706789
```

M4:
```text
0.002208966635259468  -0.10523765378969285
-0.10523765378969285  -0.80005104807768312
```

M5:
```text
0.0093858101268214696  0.057945075616615403
0.057945075616615403  3.0029704814201703
```

Mdot=sum(Mi*bd_i), with dz=f./scale and bd=[0,dz1,dz2,2*z1*dz1,z2*dz1+z1*dz2,2*z2*dz2]. Qc=J'M+MJ+Mdot+2*lambda*M includes the full Lie derivative.

## Verification

Directed mpmath interval arithmetic (30 decimal digits) validates actual trigonometric functions and the implicit Q-V root. Uniform root brackets, positive implicit derivative and interval Newton establish the electrical closure. Natural and derivative-based Taylor enclosures verify the inequalities. All 1294 terminal boxes passed (2587 total tests), with no pending boxes. The certificate concerns the fixed published decimal network and metric coefficients and exact omega=100*pi. It does not cover parameter-identification uncertainty or the full EMT controller. This is interval validation, not an SOS certificate.

Independent finite-difference Jacobian error: 1.836e-7. Directional Mdot error: 1.155e-8. Total second-derivative error: 1.914e-13; 44 Taylor-corner tests passed. Full machine-readable details are in the adjacent JSON files.

## Containment and remaining gap

Bounds came from the saved primary inactive recovery data: .05 deg lower padding, upper endpoint switch-.05 deg, and 2e-5 nu padding. No new simulation was run. All 6060 saved inactive samples and the nominal equilibrium are contained. The rectangle contains the observed minimum-angle state. Sample containment is not continuous-time containment: a thin segment just after release can lie between switch=24.373430419712 deg and the rectangle upper endpoint even if no saved inactive sample lies there. Recorded release time: 0.428985231116 s.

The preserved active domain starts at 24.423430419712 deg. The angular gap is therefore 0.1 deg, spanning (24.323430419712,24.423430419712). The inactive nu interval lies within the active nu interval [.996174514620,1.005178107392]; there is no further nu gap on their overlap.

## Certified small invariant subset

A separate constant quadratic Lyapunov function V=(x-xeq)'P(x-xeq), with P=M at the numerical nominal center, is checked by directed interval arithmetic. The center xeq is the exact equilibrium enclosed by root bracketing, not merely a rounded approximate root.

```json
{
  "validated": true,
  "equilibriumGammaBracket_rad": [
    0.36012957762191905,
    0.36012957762391906
  ],
  "P": [
    [
      0.001122016393126452,
      0.005840231645732315
    ],
    [
      0.005840231645732315,
      0.878671873769381
    ]
  ],
  "halfwidth_gamma_rad": ".001",
  "halfwidth_nu": ".00005",
  "rho": "5.415991846817519e-10",
  "definition": "(x-xeq)^T P (x-xeq) <= rho; xeq is exact unique implicit equilibrium, nu=1",
  "proof": "J^T P+P J negative throughout containing box; mean-value integral establishes Vdot<0; ellipsoid contained strictly in box."
}
```
The checker validates J'P+PJ<0 throughout the containing box with half-widths .001 rad and .00005 pu, including current/voltage feasibility. The mean-value integral for f(x)-f(xeq) gives Vdot<0 off equilibrium. The ellipsoid V<=rho is strictly contained in the box, as verified using P inverse, and is therefore forward-invariant and locally attracting. Rounded displayed root endpoints are illustrative; directed decimal endpoints are used internally.

The entire large contraction rectangle is NOT asserted forward-invariant. Entry of a switched recovery trajectory into the small ellipsoid has not been certified. Metric matching, switching-gap treatment, continuous-time containment and connection to terminal recovery remain necessary. Two separately contracting regions do not yet constitute a global recovery certificate.

## Files and preservation

New scripts: radial_inactive_metric_search.m, radial_inactive_metric_rate_trial.m, check_radial_inactive_metric.m, validate_radial_inactive_intervals.py, check_inactive_interval_math.py, check_inactive_invariant_subset.py. Network and trajectory data are read from the prior saved experiment. No new faults, EMT, physical tuning, baseline modifications or active-metric reoptimization. All 59 protected baseline hashes remain unchanged.

Method reference: Aylward, Parrilo and Slotine (2008), https://web.mit.edu/nsl/www/preprints/SOS_metrics.pdf. Polynomial-metric synthesis with the full Lie derivative; here LMIs provide candidates and directed interval bounds provide the independent uniform check.
