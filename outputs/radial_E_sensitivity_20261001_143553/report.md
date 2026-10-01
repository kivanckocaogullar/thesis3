# Radial limiter E sensitivity, physical Xeq=0.350 pu

Frozen controls/network: SCR5/XR10 external RL; reconstructed transformer, shunt and added reactor; H=3.4 K=20 Kd=5 D=25, Imax=1.2, voltage-command radial limit 1.3, active damping 0.2. E alone is varied as a frozen target in an ideal reference/current-tracking benchmark. Saved model and 59 protected file hashes remain unchanged.

Results: see comparison.csv. Positive-angle principal branch 0<gamma<=180 deg is screened; negative-angle importing equilibria are not covered. Smooth roots of Pe_gamma, Pe and |u|-1.3 delimit candidate intervals. Analytic Pe_gamma is checked by finite differences and the reconstructed converter current by network KCL. Numerical root results are not interval-arithmetic certificates.

q=(E exp(j(gamma+arg(Vth)))-Vth)/Zeq.
Radial s=min(1,Imax/|q|) q. Rotation into the VSM dq frame followed by radial scaling and inverse rotation gives the same result (radial limiting commutes with rotation).
Vc=Vth+Zeq*s; Vgrid=vgE*Vc+vgU*U; Igrid=igE*Vc+igU*U.
Pe+jQgrid=Vgrid*conj(Igrid); u_required=Vc+(Zconverter+0.2)*s.
During saturation Vc is reconstructed and is NOT fixed at E exp(jtheta).

Pe_gamma is the analytic fixed-E derivative on each smooth branch. At nu0=1, Pstar=Pe(gamma0), J=[0,1;-Pe_gamma/6.8,-25/6.8]; trace=-25/6.8 and determinant=Pe_gamma/6.8. Thus every reported positive-slope interval is locally asymptotically stable in this two-state benchmark. This is not a computed contraction region or a common-metric result. I=1.2 exactly on the active branch. Every reported interval also has Pe>0 and |u_required|<=1.3.

Examples preserving a feasible UNSATURATED nominal Pstar=1:
E=.875: active interval (24.994574206146,28.493771) deg, Pstar about 1.021989..1.02454.
E=.89: active interval (25.040119361334,26.628001) deg, Pstar about 1.044865..1.04537.
E=.9: active interval (25.050727608626,25.313809) deg, Pstar about 1.059386..1.05940.
Lower boundary: current activation; upper boundary: Pe_gamma=0. These are open intervals. Required voltage is not the limiting boundary in these low-E cases.

Numerically refined low-E window for coexistence of a strictly unsaturated nominal and a positive-width active stable interval:
0.861336941872 < E < 0.901943039019.
Lower endpoint has Pe(gamma_switch)=1 (nominal exactly at activation); upper endpoint has Pe_gamma(gamma_switch+)=0 and the active stable interval collapses. This identifies the connected low-E candidate window, not a proof excluding all conceivable E values/controllers.
E=.86 retains a feasible Pstar=1 point, but it is ALREADY current-limited (gamma=25.905022425963 deg, I=1.2, u=1.092193280060).
E=.85 cannot reach Pstar=1 on the tested positive-angle branch (peak Pe approximately .990262).
High E=1.15,1.2,1.3 gives positive slope at activation, but no intersection with voltage-feasible stable exporting active equilibria. Their switching voltage requirements are 1.39397,1.43745,1.51283 pu respectively; their nominal Pstar=1 points also violate 1.3. Merely increasing E therefore does not solve this constrained problem in the tested range.

## Actual Q-V compatibility

Existing model inspection: GFMVSG1/Grid_Forming_Converter/Converter_Control/Outer_Control/Virtual_Synchronous_Machine; normal-mode OuterCtrl_Lvirt=0. Vouter=v_set+yQ, with yQ a filtered 0.1(q_set-Qgrid), TQ=.2 s. DC relation: Eref=v_set+0.1(Qstar-Qgrid). Capacitor-voltage PI integrators, freeze logic, inner current PI and their states are not represented by q(gamma).

For Qstar=0, the required setpoint at an individual frozen-E equilibrium is v_set_required=E+0.1Qgrid. Nominal values:
E=.875 -> v_set=.840735476598.
E=.89 -> v_set=.859466803141.
E=.9 -> v_set=.871995762156.
These require substantial voltage-setpoint changes; they do NOT preserve the actual original v_set=1,Qstar=0 operating condition. At v_set=1 the equivalent required Qstar values are approximately -1.59265,-1.40533,-1.28004 pu respectively, not proposed controller settings.

Even the required v_set varies along each fixed-E active equilibrium family (ranges in CSV). Setting v_set to its nominal calibrated value once does not hold E exactly constant at every gamma. With actual Q-V feedback, E changes with Q and the total Pe derivative includes that dependence; the Q filter adds a state if retained. Thus the fixed-E positive-slope intervals are candidates for a next controller-consistent analysis, not verified intervals of the full PI/Q-V controller.

Solving the unchanged quasi-steady Q-V nominal condition jointly with Pe=1 in the .350 network gives E=1.000698736668, gamma=20.633904873783 deg. This is different from the frozen old E=.997477277664, and illustrates why that scalar is not an exact controller setpoint.

Recommendation: retain the radial algorithm. E=.875 and E=.89 are useful mathematical candidates, subject to acceptable voltage operation and a Q-V/PI-consistent equilibrium calculation. Do not claim contraction or tune M until this compatibility is resolved. If original voltage/Q setpoints must remain near nominal, investigate an explicitly modeled virtual-impedance/control-reference modification instead of assigning an arbitrary lower physical E; no success for that untested modification is claimed here. No model edits, fault simulations, EMT, metric optimization or figures were performed.
