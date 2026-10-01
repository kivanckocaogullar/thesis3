# Limiter-active transient-domain preparation — radial SCR5/XR10 Xeq=.350

This is numerical transient-domain identification, NOT a contraction certificate or a recovery certificate. No M/Theta optimization, controller tuning, Simulink edits or frozen-result overwrites were performed. All 59 protected hashes remain unchanged.

## Model/time audit

Fixed physical network: the .0858989549982 mH additional LV interfacing reactor, original transformer/shunt/converter filter, external SCR5/XR10 RL, 2-MVA bases. H=3.4 K=20 Kd=5 D=25; Vset=1 Qstar=0; radial reference Imax=1.2 and original radial voltage-command maximum 1.3. Pstar=1 CONSTANT throughout every experiment, including the Jacobian.

Physical time t is seconds:
  gamma_dot = omega_b*(nu-1), omega_b=2*pi*50;
  nu_dot = ((1-Pe(gamma;U))/nu -25*(nu-1))/6.8.
The VSG frequency equation retains its physical-time H convention. This is NOT a uniform time rescaling of the earlier normalized-angle model (which used gamma_dot=nu-1 while retaining the same frequency-equation coefficients). Do not compare the old normalized duration 1 directly with .1 s here.

PLL frequency is approximated by nominal 1; all network impedance phasors are evaluated at 50 Hz. The actual dynamic PLL, measurement filters, inner PI/current tracking and electrical transients remain omitted, as in the ideal reference benchmark. No implication of full PI/EMT validation is made.

Q-V primary approximation: its actual .2-s filter is replaced by the instantaneous equation Eref+.1 Qgrid(gamma,Eref,U)-1=0 at every evaluation. This creates a two-state quasi-static Q-V reduction. It is NOT a frozen Eref model and not the actual dynamic filter. Eref jumps at sag application and clearing; an actual filter state would remain continuous. Since the fault lasts .1 s, shorter than .2 s, this approximation is material and is not established as quantitatively accurate for the physical transient. A dynamic Q-V investigation would need an additional state and a different Jacobian/domain.

At sag, the physical external source is reduced to .6 of nominal (20 kV -> 12 kV). On the existing 20.4-kV HV base, U changes from .980392156863 to .588235294118. The nodal reduction is rebuilt with this source; Vth scales by .6 and all voltages/currents/Pe/Q and implicit Eref are recalculated. Zeq and passive nodal coefficients remain unchanged because none of the linear passive components change. The pre-fault Pe(gamma) curve is NOT reused during the sag.

## Electrical vector field

theta=gamma+arg(Vth), q=(Eref exp(j theta)-Vth)/Zeq.
Iconv=q when |q|<=1.2; Iconv=1.2 q/|q| on the active branch.
Radial limiting commutes with the VSM dq/network-frame rotation. q remains an ideal raw PI-reference proxy.
Vc=Vth+Zeq Iconv; Vgrid=vgE Vc+vgU U; Igrid=igE Vc+igU U;
Pe+jQgrid=Vgrid conj(Igrid);
u_required=Vc+(Zconverter+.2)Iconv.
The actual saturated capacitor voltage is reconstructed, not imposed equal to Eref exp(j theta). The 1.3 limit is checked; no inconsistent clipping of u while retaining ideal current tracking is used. Trajectories are set to stop if voltage feasibility is lost. No such stop occurred.

On each smooth branch, with partial derivatives at fixed Eref,
E_gamma=-.1 Q_gamma/(1+.1 Q_E),
P_gamma_total=P_gamma+P_E E_gamma.

The physical-time fixed-Pstar Jacobian is
J=[[0,314.159265358979],[-P_gamma_total/(6.8 nu), (-(1-Pe)/nu^2-25)/6.8]].
The off-equilibrium term -(1-Pe)/nu^2 is retained. Pstar is NEVER replaced by Pe(gamma) inside the vector field/Jacobian.

## Nominal and disturbance

Nominal gamma=20.633904875622 deg, nu=1, Eref=1.000698736667, Iconv=1.017813065253 pu, u_required=1.210039788936 pu, physical Vc=1.000698736667 pu and PCC=.978868334506 pu. Current limitation is inactive.

Adaptation tested retained voltage .8 and .7 first, at .1 s duration; neither activated the current limiter during the fault. Retained voltage .6 did activate it. No further sag-depth sweep was made.

Primary fault: sag from t=0 to t=.100000 s, starting exactly at the nominal state. No angle or frequency jump is applied to this primary experiment.

| Primary trajectory quantity | Value |
|---|---|
| Limiter activation | t=0+, immediately from source-voltage step |
| Clearing state | gamma=25.027638262623 deg, nu=1.004481785586 |
| Current limiter at clearing, both sides | Active |
| Limiter release | t=.428985231116 s |
| Total continuous active duration | .428985231116 s |
| Gamma range over 8 s | [16.476345885756,30.257330295649] deg, unwrapped |
| Nu range | [.996547139105,1.004481785586] |
| Pe range | [.615813546219,1.173082735657] pu |
| Raw current-reference range | [.814044239081,1.487596434689] pu |
| Active limited-current magnitude | 1.2 pu |
| Required voltage-command range | [1.106505379054,1.248870881579] pu |
| State at 8 s | gamma=20.633900997646 deg, nu=.999999998099 |
| Final power residual Pe-1 | -1.81647e-7 pu |

At t=0+ Pe=.615813546219 and Eref=.928331072129. At clearing, Pe changes from .663207024726 to 1.171470353387 and Eref from .936286998135 to 1.001782541255; gamma and nu remain continuous.

27 numerical runs cover initial gamma offsets [-.1,0,+.1] deg, initial nu offsets [-1e-4,0,+1e-4], and fault durations [.08,.10,.12] s. These small initial offsets are solely the requested neighbourhood study; the central run starts exactly nominal. Every sampled run recovered under the same model. Recovery means maximum original-equilibrium angle error over t=7..8 s <.01 deg, maximum |nu-1|<1e-5, and final |Pe-1|<1e-4. This is observed recovery for 27 runs, not a guarantee for the continuous initial set/duration interval.

Some .08-s runs release instantaneously at clearing, re-enter current limitation afterward and release again. Exact per-mode active intervals including t=0 and clearing jumps are recorded in active_time_intervals.csv; ensemble_summary.csv additionally retains smooth event times. Do not interpret an empty event list during the fault as absence of activation: all these faults start on the active branch.

## Two-dimensional candidate regions

Different source voltages give different vector fields. Keep the two regions separate for any subsequent optimization.

| Domain | gamma bounds, deg | nu bounds | Total Pg range, pu/rad | Min raw-current excess above 1.2 | Min voltage headroom to 1.3 |
|---|---|---|---|---|---|
| Fault-on Uret=.6 | [20.482353937532,27.047756563847] | [.999880000000,1.005178107392] | [.440653848565,.734843795575] | .005857013449 pu | .174535336238 pu |
| Postfault active interior | [24.423430419712,33.858410380369] | [.996174514620,1.005178107392] | [-.190495953490,-.140353448980] | .002451136739 pu | .051276946768 pu |

The fault-on rectangle encloses the sampled active cloud plus .05-deg and 2e-5-nu padding. The initial postfault rectangle crossed the current switch and was REJECTED. Its corrected lower gamma bound is the postfault switch 24.373430419712 deg plus .05 deg. Consequently the postfault interior deliberately excludes a thin activation/release layer; it is not an enclosure of the complete trajectory all the way to release. That layer needs separate nonsmooth treatment. The upper gamma and nu bounds retain the ensemble/padding limits.

Both boxes were sampled on 501x101 grids with implicit Q-V solved throughout. All sampled points lie on the current-active branch with |Iconv|=1.2 and |u_required|<1.3. Analytic off-equilibrium Jacobians were independently finite-difference checked at 77 points per box with fixed Pstar. These are finite two-dimensional numerical domains, not just curve samples, but not interval-arithmetic validation or forward-invariant sets.

Jacobian entry ranges (J12=314.159265359, J11=0):
- Fault-on: J21=[-.108078233443,-.064468213222], J22=[-3.733266464935,-3.723032293160]. Maximum real eigenvalue at a point is in [-1.866633232467,-1.861516146580] s^-1.
- Postfault: J21=[.020533886416,.028121689921], J22=[-3.655264332100,-3.650838419159]. Maximum real eigenvalue is in [1.302274933802,1.661718582752] s^-1.

Fault-on Pe in the proposed box is only .613881.. .680093, so no Pstar=1 equilibrium lies there. Positive total Pg during the sag is compatible with the earlier negative postfault active-family result: the external voltage and hence electrical characteristic are different.

## Fixed-Pstar saddles and obstacles

For the exact pre/postfault model at Pstar=1 the active saddle is gamma=64.880551979056 deg, nu=1, u_required=1.092193280060 pu; physical-time eigenvalues are +2.55690231355 and -6.23337290179 s^-1. It is outside the postfault box (upper bound 33.858410380369 deg). It is not a reused old-SCR5 angle. The nominal equilibrium eigenvalues are -1.83823529412 +/- j10.9823028726 s^-1.

The fault-on model has no Pstar=1 root on the scanned positive principal-angle interval 0..180 deg; no saddle from that interval lies in its candidate box. Periodic angular copies do not intersect either small box.

Potential for next optimization, without claiming success:
- Fault-on pointwise Hurwitz Jacobians and positive Pg make this a sensible first finite-domain metric-search candidate, but do not prove one common contraction metric exists.
- The postfault active interior has positive instantaneous Jacobian eigenvalues, so a constant positive-definite metric cannot satisfy the continuous-time contraction LMI there. A state-dependent/time-dependent metric could include Mdot on a non-equilibrium transient domain; that possibility has NOT been tested.
- Neither rectangle is asserted to be invariant; trajectories enter/leave, and the fault/clearing vector field changes. Metric matching/switching and the limiter boundary require separate checks. Positive-rate contraction on a finite transient domain would not by itself prove eventual recovery.
- The instantaneous Q-V approximation is particularly important on a .1-s sag. Retaining the actual .2-s Q-V filter, PLL and PI dynamics changes the problem and its state dimension.

Files: fault_trajectory.csv and postfault_trajectory.csv preserve both sides of clearing; ensemble_summary.csv; active_time_intervals.csv; candidate_domains.csv/.mat; results.mat; console.txt and final_checks.txt. No figures generated.
