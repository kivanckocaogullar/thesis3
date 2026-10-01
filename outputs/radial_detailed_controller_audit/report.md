# Detailed radial PI/Q-V investigation — partial audit, NOT completed stability verification

Status: source inspection and isolated physical variant completed. Detailed-controller equilibrium search and full local spectrum NOT verified because model update/compilation did not return. Do not classify the actual controller as A (saddle), B (stable different characteristic), or C (no persistent active equilibrium) on this evidence. No numerical equilibrium candidate is certified by this audit.

## Physical variant and unchanged sources

The saved GFMVSG1.slx contains the external SCR5/XR10 RL and the original converter filter/transformer, but NOT the additional .0858989549982-mH LV reactor used by the Xeq=.350 analytical variant.

New isolated model: GFMVSG1_Xeq0350_audit.slx.
New isolated electrical referenced subsystem: GFM_Electrical_Xeq0350.slx.
Added block: GFMVSG1_Xeq0350_audit/Grid_Forming_Converter/Electrical_Part/Additional_Interfacing_Reactor.
BranchType=L, L=8.58989549982e-5 H; connected phase-by-phase between LC_harmonic_filter grid-side ports and Transformer_Model LV ports. This is an actual additional block, not a change to external SCR or transformer winding tuning. The main converter library link is disabled only in the COPY to redirect its electrical subsystem reference. All controller subsystems retain their original references and parameters. The copied model has been saved but NOT successfully compiled or simulated.

Original GFMVSG1.slx and frozen baselines were not saved or overwritten. The copied electrical subsystem's .slx.r2022b file is Simulink's automatic version backup inside this new output directory.

## Verified controller data and equilibrium equations

Current mask values re-read from the saved model: H=3.4, K=20, Kd=5, Kp_v=.8, Ki_v=.5, Kp_i=3.3, Ki_i=16.5, Imax=1.2, umax=1.3, AD=.2, Lvirt=0, Q-V droop=.1, TQ=.2 s, T_s_ctrl=T_s_PWM_meas=T_s_SPS=1e-4 s. Measurement low-pass enabled at 2500 Hz, positive-sequence filter cutoff 20 Hz. Average bridge model. This corrects any obsolete values in the September-30 inspection report.

Paths relative to C=GFMVSG1/Grid_Forming_Converter/Converter_Control:
- Outer_Control/Virtual_Synchronous_Machine
- Operation_Mode_Selector
- Inner_Control/Cascaded_PI_controller/v_cap_Controller
- Inner_Control/Cascaded_PI_controller/i*_Limiter/dq_lim/dq_Limiter_i_conv*
- Inner_Control/Cascaded_PI_controller/i_conv_Controller
- Inner_Control/Cascaded_PI_controller/dq_Limiter_v_conv*
- Inner_Control/Cascaded_PI_controller/v_cap_Controller/Integrator_Freeze_Logic
- Inner_Control/Cascaded_PI_controller/i_conv_Controller/Integrator_Freeze_Logic
- measurement_transformations/Power_Calculation_(Positive_Sequence)
- Grid_Synchronization/Synchronization_Unit

All following dq quantities with subscript m are MEASURED/filter outputs, not silently replaced by physical phasors. Let J2=[0,-1;1,0]. In normal VSM mode:

vstar=[1+yQ;0],
ev=vstar-vc_m,
r=.8*ev+xi_v+if_m+nu*.04569*J2*vc_m,
s=sat_radial(r,1.2), aI=(norm(r)>1.2),
ei=s-ic_m,
u_cmd=3.3*ei+xi_i+vc_m+nu*.05097*J2*ic_m,
w=sat_radial(u_cmd,1.3), aV=(norm(u_cmd)>1.3),
v_bridge_command=w-.2*ic_dq0_m.

Here u_cmd denotes voltage command, NOT the frequency nu. The current limiter acts on the actual PI expression r. It is not q=(Eref exp(jtheta)-Vth)/Zeq. Active damping occurs AFTER the 1.3 limiter.

For each a=d,q, normal-mode F_force=0 and delayed steady flags satisfy:
F_va=aI AND [sign(r_a)==sign(ev_a)],
F_ia=aV AND [sign(u_cmd,a)==sign(ei_a)].
The implementation includes one-sample delay memories in these tests, which must equal their steady input values in an equilibrium.

xi_va[k+1]-xi_va[k]=Ts*.5*(1-F_va)*ev_a,
xi_ia[k+1]-xi_ia[k]=Ts*16.5*(1-F_ia)*ei_a.
Every such increment must vanish. Reset must be inactive. With norm(u_cmd)<1.3, aV=0 and the current PI requires ei=0 in steady state; this equates LIMITED reference to MEASURED converter current, not automatically unfiltered physical current.

VSM equilibrium balance:
0=(Pstar-Pgrid_m)/max(nu,.5)+20*(1-nu)-5*(nu-f_PLL).
At locked nominal frequency nu=f_PLL=1, Pstar=Pgrid_m. This equality must use actual filtered measured power from the solved network/controller point.
Physical-time angle dynamics are theta_dot=2*pi*50*nu; relative synchronous gamma is constant at nu=1. Do not splice the old normalized two-state Jacobian into PI dynamics measured in seconds.

The actual Q-V block is a discrete filter:
yQ[k]-alphaQ*yQ[k-1]=(1-alphaQ)/2 * (vQ[k]+vQ[k-1]),
vQ=.1*(0-Qgrid_m).
Thus yQ=-.1Qgrid_m in DC steady state; its state must still remain in the dynamic linearization. Its exact filter numerator/denominator, and measurement filters, are in model_inspection.txt. The raw abc measurement low-pass has numerator [(1-alpha_meas_LPF)/2,(1-alpha_meas_LPF)/2], denominator [1,-alpha_meas_LPF], Ts=T_s_PWM_meas. DDSRF decoupling/positive-sequence filter states and PLL states have not been eliminated.

Physical network steady-state requirements include Vbridge=Vc+Zconverter*Iconv, Vc=Vth+Zeq*Iconv and reconstructed transformer nodes/grid phasors for the .350 network. These fundamental relations alone do not enforce all sampled/periodic electrical and measurement states of the EMT realization. Balanced abc sinusoidal steady state is a periodic orbit, not a constant full abc Simulink state vector; a complete check needs consistent synchronous-frame states or periodic trimming/Floquet multipliers. No full PLL/measurement/network state residual is claimed here.

## Structural neutral-mode finding (conditional, not a full numerical spectrum)

For a constant synchronous-dq, strictly current-limited equilibrium, let r0 be its raw PI reference with norm(r0)>1.2. If a nearby equilibrium retains the same normal-mode freeze branch, replacing
xi_v -> xi_v + epsilon*r0
makes r -> (1+epsilon)*r0 while sat_radial(r,1.2) stays exactly unchanged. Physical outputs, errors and current-loop signals can remain unchanged. The signs of r also remain unchanged for small epsilon.

For each voltage integrator, either it is frozen or its equilibrium error is zero, so its steady update still vanishes. This produces a local internal-memory equilibrium freedom under the stated constant-dq assumptions. The limiter derivative Ds=(1.2/norm(r0))*(I-rhat*rhat') has a radial null direction.

Consequences:
- Any integrator strictly frozen throughout the local mode has an identity update row, hence a discrete unit eigenvalue in that mode. Its continuous surrogate has a zero eigenvalue.
- Such a neutral mode prevents local exponential/asymptotic convergence to one isolated FULL-STATE equilibrium and prevents a strict positive-rate full-state contraction certificate in that smooth mode.
- A unit multiplier does NOT by itself establish instability, nor does it prove that all other modes are stable. Electrical outputs might be stable transverse to an equilibrium manifold.
- At ev_a=0, r_a=0, limiter equality, or release boundaries the hybrid sign logic is nonsmooth; a single frozen-branch Jacobian is insufficient.
- For an actual sampled periodic orbit with ripple and changing freeze flags, the simple constant-dq scaling argument is not automatically an exact invariance. Its period map must be evaluated.

Therefore this structural observation cannot substitute for the requested complete detailed-controller eigenvalue calculation. It also cannot be described as confirmation of the old two-state saddle mechanism.

## Requested numerical candidate table — not yet available

| Quantity | Verified detailed-controller candidates |
|---|---|
| gamma, Pstar, measured Pe/Q | Not solved |
| Raw/limited current, physical Vc/PCC, actual u_cmd | Not solved |
| PI/Q-V/PLL/measurement/network residuals | Not verified |
| Freeze/release consistency of a solved point | Not verified |
| Full continuous/discrete local spectrum | Not computed |
| Structural neutral mode | Conditional result above, not a complete spectrum |

Zero verified candidates is NOT a finding that zero equilibria exist.

## Execution limitation and next unresolved work

The original model's SimulationCommand='update' and repeated update attempts on the isolated .350 copy failed to return a completion or a MATLAB exception. The batch calls were interrupted after waiting; no root cause was established. The variant_compile.txt log ends after the model was loaded and the initial watchdog message; no COMPILE SUCCESS or state inventory was obtained. The source loads, can be inspected, and the physical copy saves successfully. This is an unresolved compilation/execution limitation, not a demonstrated model-equation inconsistency or missing license.

Remaining work: resolve model compilation, construct a synchronous/periodic operating point with the actual states, search PI/freeze-consistent active equilibria and assess the full discrete/continuous or period-map spectrum. Until then A/B/C remains unresolved. No faults, metric optimization, controller retuning or original-model modifications were performed.
