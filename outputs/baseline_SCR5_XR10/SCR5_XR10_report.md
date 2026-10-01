# Limiter-free SCR=5, X/R=10 baseline

Generated from saved GFMVSG1 model; normalized analytical time. Frozen analytical files were not changed.

## Network and scope

Bridge -> converter RL -> capacitor node (shunt series RC) -> transformer T equivalent -> closed breaker -> external RL -> ideal source.

Old PI line and old nonideal source removed; fault block unconnected; no local load.

No_Transformer grid inductor is inactive; active damping is not passive series resistance.

Sbase=2000000 VA; capacitor Vbase=690 V line-line; transformer HV base=20400 V (1.02 times measurement nominal). Physical ideal bus=20000 V, referred U=0.980392156863 pu. Transformer D1 phase shift is removed on referral. Rgrid=3.98014876 ohm; Lgrid=0.12669207 H.

E is the physical capacitor-node fundamental amplitude under ideal inner-loop tracking and unity measurement transfer approximation. Constant E=0.997477277664 is calibrated at Pstar=1 using E=1-0.1 Qgrid. The secondary algebraic Q-V curve is sensitivity only; it does not change primary derivatives or certificates. Dynamic Q-V, voltage/current control and measurement states are omitted. This is not an exact EMT reduction.

Transformer z1=0.005+j0.02, z2=0.005+j0.02; Ym=0.002+j-0.002; series-RC shunt Yc=0.000143519441+j0.0456895492. Converter RL is upstream of E and excluded from downstream equivalent.

Nodal reduction: Y=[1/z1+1/z2+Ym,-1/z2;-1/z2,1/z2+1/(zg+zb)]. Nodes=Y\[Vcap/z1;U/(zg+zb)]. Ig=(node2-U)/(zg+zb), Vg=node2-zb Ig, Ic=(Vcap-node1)/z1+Yc Vcap. Zeq=1/icE, Vth=-icU/icE U.

Zeq=0.0298789620445+j0.233596600585 pu; Vth=0.990384568064+j-0.00174104618759 pu. Gamma=arg(Vcap)-arg(Vth); Vcap=E exp(j(gamma+phi)), phi=-0.100722990788 deg.

Power is Re(Vg conj(Ig)) at the grid measurement, NOT capacitor source power. With actual referred bus magnitude U, Pe=a E^2+b U^2+E U(c cos(gamma)+d sin(gamma)); [a,b,c,d]=[0.351963430883615 -0.184913782926074 -0.174685441407252 4.25553863553776].

Pe_gamma=EU(-c sin(gamma)+d cos(gamma)); Pe_U=2bU+E(c cos(gamma)+d sin(gamma)); Pe_gammaU=E(-c sin(gamma)+d cos(gamma)).

Current block paths (prefix GFMVSG1/Grid_Forming_Converter):

- /Converter_Control/Outer_Control/Virtual_Synchronous_Machine
- /Converter_Control/Operation_Mode_Selector
- /Converter_Control/Inner_Control/Cascaded_PI_controller/v_cap_Controller
- /Electrical_Part/LC_harmonic_filter/V_cap_abc
- /Converter_Control/measurement_transformations/Power_Calculation_(Positive_Sequence)

The VSM Q-V amplitude command feeds the mode selector (Lvirt=0), then the capacitor voltage controller; physical V_cap is the chosen E node. The measured grid positive-sequence voltage/current feed active power. Bridge voltage is not the fixed E source.

LV bases: I_RMS=1673.47904113 A, Z=0.23805 ohm. HV transformer bases: I_RMS=56.6029675676 A, Z=208.08 ohm. P and Q use Sbase. Peak dq and RMS phasors have identical pu magnitude on compatible bases. The Thevenin equivalent is exact for this linear positive-sequence network at 50 Hz; the control/measurement quasi-steady assumptions are the approximation.

## Dynamics and unchanged metric construction

gamma_dot=nu-1; nu_dot=((Pstar-Pe)/nu-D_eff(nu-1))/(2H). H=3.4, K=20, Kd=5, D_eff=25. Each equilibrium has its own fixed Pstar=Pe(gamma0), nu0=1.

J=[0,1;-Pe_gamma/(2H nu),(-(Pstar-Pe)/nu^2-D_eff)/(2H)]. beta=-Pe_U/(2H), c_metric=-Pe_gammaU/(2H), Theta=[-beta,0;-c_metric,-beta], M=Theta'Theta. F=Theta J/Theta with fixed reference Theta (Theta_dot=0). No optimized metric.

Full-period scan: -180..180 deg, 0.0005 deg spacing, refined scalar roots. Certificate: min eig(M)>1e-9 and max eig(Sym(F))<-1e-9. CSV distinguishes tolerance boundaries and mathematical zeros. Local stability and sufficient contraction are different tests. Singular Theta points are excluded.

## Results

Nominal angle=13.8107589411 deg. Nominal is locally stable; certificate=0. All Pstar=1 roots in one period are in nominal_equilibria.csv; copies repeat by 360 deg. The equilibrium family always sets Pstar=Pe; a fixed-Pstar root can disappear at a power extremum.

- R1: 22.3523784001 to 29.217483849 deg.
- R2: 50.7057644612 to 84.7903294107 deg.

Boundary causes and structural points are exported separately. All primary changes arise from electrical network, power measurement mapping and nominal E calibration; controller/metric method stays fixed. These changes are not attributable to X alone.

### Old/new numerical comparison

| Quantity | Old | New | Change |
|---|---:|---:|---:|
| E | 1 | 0.997477277664 | -0.00252272233582 |
| Vth magnitude | 1 | 0.990386098399 | -0.00961390160133 |
| Req | 0 | 0.0298789620445 | 0.0298789620445 |
| Xeq | 0.3 | 0.233596600585 | -0.0664033994155 |
| Nominal gamma | 17.4576031237 | 13.8107589411 | -3.64684418258 |
| H | 3.4 | 3.4 | 0 |
| D_eff | 25 | 25 | 0 |
| R1 boundary 1 | 15.3225991806 | 22.3523784001 | 7.02977921957 |
| R1 boundary 2 | 21.0916682548 | 29.217483849 | 8.12581559428 |
| R2 boundary 1 | 48.6263618306 | 50.7057644612 | 2.07940263062 |
| R2 boundary 2 | 81.4159628061 | 84.7903294107 | 3.37436660465 |

Local stability intervals (deg):
- -87.649387093 to 92.3506122192

Nominal Pstar=1 mathematical equilibria (modulo 360 deg):
- gamma0=13.8107589411 deg; J eigenvalues -0.171259194103, -3.50521139413; Locally stable.
- gamma0=170.890466185 deg; J eigenvalues 0.156610232499, -3.83308082073; Locally unstable.

Full matrices for every selected point are stored in MAT pointMatrices and printed in the console log. The second certified region survives.

## Trajectories and voltage sags

Representative trajectories use our common +1 deg test (not a paper-specified magnitude); unstable main panels show early departure. Full diagnostics retain unwrapped angles and L relative to ORIGINAL equilibrium, so another branch does not yield L=0. A neighboring trajectory starts at +1.01 deg.

Sag references include the new nominal and midpoints of all positive certified regions. Nominal is explicitly labelled uncertified. Initial states are exact equilibria, no angle jump. Retained voltage 0.9..0.1 with duration 1; duration sweep 1,2,3,4,5,7,10 at 0.1 pu. Pstar is fixed during sag; U is restored at clearing.

Recovery uses the last 10 percent of elapsed trajectory: angle error <0.1 deg and |nu-1|<1e-4. Other-branch test also checks power residual <1e-4. nu=0.5 terminates integration rather than modifying dynamics. Tests use a fixed reference metric: sampled instantaneous test, not a regional proof. Discrete recovering tests are not critical voltages or CCT.

## Verification

- powerDerivative: 7.62888e-11
- Jacobian: 3.58954e-11
- betaDerivative: 9.14124e-12
- cDerivative: 1.12019e-11
- matrixSymmetry: 0
- matrixFormula: 1.01436e-15
- equilibriumResidual: 0
- phasorPower: 1.11022e-15
- thevenin: 1.79018e-15
- nominalAngle: 1.47105e-15

Tighter independent sag integrations are checked and saved. PNG and editable FIG outputs include power/derivative, three matrix spectra, local/certificate regions, boundary causes, fixed-reference instantaneous regions, trajectories, Figure 5A/5B equivalents and sag magnitude/duration diagnostics. H-D retuning/sweep and unrelated legacy SG-only figures are not part of this fixed-tuning baseline.

The baseline is ready as an internally verified reference for a subsequent reduced current-limiter study within the stated approximations. No current-limiter study was started.
