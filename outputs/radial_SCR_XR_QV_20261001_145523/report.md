# Final SCR-X/R sensitivity: frozen radial limiter and Vset=1

RESULT: None of the 16 tested external-grid combinations creates a positive-width, exporting, locally stable limiter-active equilibrium interval in the specified network-consistent quasi-static Q-V / ideal current-reference benchmark. All 16 retain a feasible, locally stable, limiter-INACTIVE Pstar=1 equilibrium. This result concerns this benchmark and tested grid only, not every radial-limited GFM or the full EMT controller.

The failure is the sign of TOTAL Pe_gamma, not an invented undervoltage limit. Across every exporting limiter-active branch, the numerically determined maximum total derivative remains negative; its supremum occurs at activation. The least-negative case is SCR10/XR15, -0.0964159756051 pu/rad. Immediately at activation every case still has positive voltage-command headroom (about .04264.. .05841 pu). Thus relaxing the 1.3 limit would not fix the observed negative synchronizing slope there.

## Full comparison

Intervals are open and angles are degrees; derivatives are pu/rad. `--` means there is no interval over which to define its voltage margin.

| SCR | X/R | Reconstructed Xeq | gamma switch | Total Pg after switch | Feasible active interval | Width | Min voltage margin in interval | Nominal Pstar=1 |
|---:|---:|---:|---:|---:|---|---:|---|---|
| 3 | 3 | 466111335980878.000000000 | 333184547412626.000000000 | -266395671759838.000000000 | none | 0 | -- | feasible, stable, inactive |
| 3 | 5 | 477147057693967.000000000 | 336621292881053.000000000 | -220856323768587.000000000 | none | 0 | -- | feasible, stable, inactive |
| 3 | 10 | 482190895083943.000000000 | 338407540651268.000000000 | -186952301923511.000000000 | none | 0 | -- | feasible, stable, inactive |
| 3 | 15 | 483162188194413.000000000 | 338828136962054.000000000 | -175821188733858.000000000 | none | 0 | -- | feasible, stable, inactive |
| 5 | 3 | 340583783275419.000000000 | 240096258376183.000000000 | -209966208097276.000000000 | none | 0 | -- | feasible, stable, inactive |
| 5 | 5 | 347050698593099.000000000 | 24253554766733.000000000 | -169496876237476.000000000 | none | 0 | -- | feasible, stable, inactive |
| 5 | 10 | 349999999853243.000000000 | 243734304093614.000000000 | -140081398856864.000000000 | none | 0 | -- | feasible, stable, inactive |
| 5 | 15 | 350566253058811.000000000 | 243996043230913.000000000 | -130570818387517.000000000 | none | 0 | -- | feasible, stable, inactive |
| 7 | 3 | 287137948604549.000000000 | 20117713086961.000000000 | -180605107784392.000000000 | none | 0 | -- | feasible, stable, inactive |
| 7 | 5 | 291710684892469.000000000 | 203128103039835.000000000 | -145250038875877.000000000 | none | 0 | -- | feasible, stable, inactive |
| 7 | 10 | 293794189669668.000000000 | 204063235087875.000000000 | -119636791719602.000000000 | none | 0 | -- | feasible, stable, inactive |
| 7 | 15 | 294193695439396.000000000 | 204259907759599.000000000 | -11137874651079.000000000 | none | 0 | -- | feasible, stable, inactive |
| 10 | 3 | 247191195012604.000000000 | 17235744565268.000000000 | -153607056791066.000000000 | none | 0 | -- | feasible, stable, inactive |
| 10 | 5 | 250367924107347.000000000 | 173879691071797.000000000 | -124489371584322.000000000 | none | 0 | -- | feasible, stable, inactive |
| 10 | 10 | 251814346736462.000000000 | 17459673552015.000000000 | -103274036788991.000000000 | none | 0 | -- | feasible, stable, inactive |
| 10 | 15 | 252091422652436.000000000 | 17474337359485.000000000 | -964159756051268.000000000 | none | 0 | -- | feasible, stable, inactive |
## Actual nominal voltages and currents

These are reconstructed physical capacitor-node and grid-measurement/PCC voltages, not Eref substituted during saturation. Voltage bases remain LV 690 V and HV 20.4 kV. No minimum voltage limit is imposed.

| SCR | X/R | Nominal capacitor pu | Nominal PCC pu | Converter current pu | Required command pu |
|---:|---:|---:|---:|---:|---:|
| 3 | 3 | 100941751683208.000000 | 100166696543026.000000 | 100296945028646.000000 | 121274445469471.000000 |
| 3 | 5 | 100226540829583.000000 | 983015658680647.000000 | 101445870793333.000000 | 121048453009079.000000 |
| 3 | 10 | 996606184110206.000000 | 967923385692772.000000 | 102806582513477.000000 | 120896923887881.000000 |
| 3 | 15 | 994672604066509.000000 | 962694361144879.000000 | 103368025217956.000000 | 120851127103285.000000 |
| 5 | 3 | 101032713678729.000000 | 100400732747974.000000 | 100194103344415.000000 | 121305756828875.000000 |
| 5 | 5 | 100493614239096.000000 | 990033482533166.000000 | 10094465210381.000000 | 121128532836307.000000 |
| 5 | 10 | 100069873667408.000000 | 978868334526168.000000 | 101781306523572.000000 | 121003978893764.000000 |
| 5 | 15 | 999260940714494.000000 | 975041718192674.000000 | 102116677332563.000000 | 120964840464408.000000 |
| 7 | 3 | 100957082888013.000000 | 100206190817581.000000 | 100278947594321.000000 | 121279683719903.000000 |
| 7 | 5 | 100522951050034.000000 | 99080043060317.000000 | 100894927256613.000000 | 121137649175341.000000 |
| 7 | 10 | 100181350315224.000000 | 981821723461039.000000 | 101539446867662.000000 | 121035431542237.000000 |
| 7 | 15 | 100065540830522.000000 | 978753306658286.000000 | 101791025699573.000000 | 121002775786285.000000 |
| 10 | 3 | 100829011275981.000000 | 998756586275816.000000 | 100437631167615.000000 | 121236417442475.000000 |
| 10 | 5 | 100492896679002.000000 | 990014713836815.000000 | 100945881428942.000000 | 121128310638818.000000 |
| 10 | 10 | 100227651163645.000000 | 983044969873693.000000 | 101443603884273.000000 | 121048774899996.000000 |
| 10 | 15 | 100137693410918.000000 | 980666497193225.000000 | 101632289630775.000000 | 121023000003058.000000 |
Additional Pe=1 roots on the active branch are saved in nominal_points.csv. They can be electrically feasible but are unstable (negative total Pg); they are not silently discarded as nonexistent equilibria. Full reconstructed Vc, Vgrid, Igrid, Eref, Pe, Qgrid, current and command arrays for all 16 positive-angle families are in results.mat. No candidate active stable interval exists for an additional interval-detail table.

## Frozen physical network and baseline reproduction

The additional LV interfacing reactor was calibrated ONCE using the previous SCR5/XR10 physical network: Ladd=0.0000858989549982 H = .0858989549982 mH. It is inserted between the capacitor node and transformer primary branch. This physical inductance, transformer primary/secondary series branches, magnetizing admittance, capacitor-node series-RC shunt, converter filter, breaker and active damping are retained for all 16 cases. The added reactor is NOT retuned to keep Xeq=.350.

External grid only:
Zmag_ohm=(20000 V)^2/(2e6 VA * SCR),
Rgrid=Zmag/sqrt(1+(X/R)^2), Lgrid=(X/R)*Rgrid/(2*pi*50),
Zgrid_pu=(Rgrid+j*2*pi*50*Lgrid)/(20400^2/2e6).
This respects the 20-kV source used for SCR and the existing 20.4-kV transformer HV voltage base. All nodal coefficients, Zeq and Vth are rebuilt; LV per-unit impedance uses the consistent 690-V/2-MVA base. Complete physical R/L and complex equivalent data are saved in networks.csv.

At SCR5/XR10 the rebuilt reference has Xeq=.349999999853 pu, |delta Zeq|=1.47e-10 and |delta Vth|=6.3e-12 relative to the saved reconstructed network. Tiny differences come from stored external R/L decimal precision. Fixed-Vset nominal Eref=1.000698736674 and gamma=20.633904866585 deg reproduce the prior Q-V-consistent nominal. Switching is now 24.373430409361 deg; it is not the old frozen-E switching angle 24.421813345295 deg.

## Equations, branch tests and local stability

Frozen controller data: radial current-reference Imax=1.2; original radial voltage-command maximum=1.3; Vset=1; Qstar=0; droop=.1; H=3.4 K=20 Kd=5 D=25; original PI gains and active damping .2. No virtual impedance was added.

q=(Eref*exp(j*(gamma+arg(Vth)))-Vth)/Zeq,
Iconv=q if |q|<=1.2, otherwise Iconv=1.2*q/|q|.
Radial scaling commutes with rotating VSM dq/network transformations. This q remains an IDEAL raw-reference proxy, not a reconstructed physical PI output.

Vc=Vth+Zeq*Iconv;
Vgrid=vgE*Vc+vgU*U; Igrid=igE*Vc+igU*U;
Pe+jQgrid=Vgrid*conj(Igrid);
u_required=Vc+(Zconverter+.2)*Iconv.
No fixed Eref target is imposed as the actual saturated capacitor voltage.

Solve h=Eref+.1Qgrid-1=0 at every gamma. On smooth branches,
Eref_gamma=-.1Q_gamma/(1+.1Q_E),
Pe_total_gamma=P_gamma+P_E*Eref_gamma.

At each equilibrium nu0=1 and Pstar=Pe(gamma0), held FIXED in its Jacobian:
J=[0,1;-Pe_total_gamma/6.8,-25/6.8].
Negative total Pg gives det(J)<0 and a saddle, not an asymptotically stable operating point. A positive-definite contraction metric cannot turn such an unstable equilibrium into a contracting equilibrium. No metric optimization or certificate was attempted.

Scans cover 0..180 deg, explicitly split at activation before derivative-root refinement. Roots of total Pg, Pe and |u|-1.3 define all candidate partitions. Each partition is checked for active current, positive export, positive total derivative, current and voltage feasibility. One-sided total derivatives are evaluated separately at switching. Local maxima of the exporting active-branch derivative are also refined to verify the observed absence of a positive-slope interval. Numerical root scans are not rigorous interval-arithmetic proofs of global completeness.

Independent central and one-sided finite differences agree with analytic derivatives to below 6.7e-9 pu/rad; network KCL errors are below 1.8e-15. The Q-V implicit denominator remains above .9992 in all angle scans. Additional root scans at 37 angles per network find one positive Q-V root within Eref=.05..2 pu; global uniqueness outside that tested range is not claimed.

## Electrical interpretation and scope

Increasing SCR lowers the rebuilt interfacing reactance and moves activation to a smaller angle. Increasing X/R makes the post-activation slope less negative in this tested grid, but does not reverse its sign. Once radial saturation fixes converter-current magnitude, further reference-angle movement changes current direction and reconstructed node voltages. In these networks with Vset=1, the exported-power characteristic falls with angle on the active exporting branch; the implicit Q-V adjustment is insufficient to restore a positive synchronizing slope.

The Q-V filter is reduced to its DC relation, and current tracking is ideal. Voltage/current PI integrator dynamics, conditional anti-windup/freeze, measurements and EMT electrical dynamics are omitted. Therefore this is a physically interconnected electrical/quasi-static feasibility study, not proof of exact PI equilibria or EMT stability. These qualifications do not change the negative result inside the specified benchmark.

Only external grid strength and X/R were varied. No low-E setpoint, gain change, priority algorithm, virtual impedance, minimum voltage threshold, fault simulation, figure or M/Theta optimization was introduced. Original Simulink and protected baselines were left unchanged.
