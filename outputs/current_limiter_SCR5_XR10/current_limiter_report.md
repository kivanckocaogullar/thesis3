# SCR5/XR10 current-limiter comparison

**Network-consistent two-state ideal radial-current benchmark.** Not an exact PI/freeze controller model. All electrical coefficients are read from the saved SCR5/XR10 baseline. H=3.4, K=20, Kd=5, D_eff=25; normalized analytical time and baseline metric method unchanged.

## Electrical closure

Vtarget=E exp(j(gamma+phi)); Vth=tU; q=(Vtarget-tU)/Zeq. Inactive: Iconv=q, Vc=Vtarget. Active: Iconv=1.2 q/abs(q), Vc=tU+Zeq Iconv. Grid voltage/current and transformer nodes use baseline nodal coefficients; Pe=real(Vg conj(Ig)). Gamma is the target angle during saturation, not actual capacitor angle. No rho-times-unlimited-power shortcut is used.

Zf is excluded from downstream Zeq and used only in Vbridge=Vc+Zf Iconv; u_required=Vbridge+0.2 Iconv. Points exceeding abs(u_required)=1.3 are marked voltage-invalid, with no added voltage saturation. The primary certificates require voltage validity.

## Derivatives and metric

All branches use analytic Pe_gamma, Pe_U and Pe_gammaU from voltage/current reconstruction. For r=abs(q), r_g=real(conj(q)q_g)/r, r_U=real(conj(q)q_U)/r, r_gU=real(conj(q_U)q_g)/r-r_g r_U/r. The active mixed current derivative is 1.2[-q_g r_U/r^2-q_U r_g/r^2-q r_gU/r^2+2q r_g r_U/r^3]. q_g=j Vtarget/Zeq and q_U=-t/Zeq.

J=[0,1;-Pe_gamma/(6.8nu),(-(Pstar-Pe)/nu^2-25)/6.8]. At each new equilibrium Pstar=Pe(gamma0), nu0=1. beta=-Pe_U/6.8; c=-Pe_gammaU/6.8; Theta=[-beta,0;-c,-beta]; M=Theta transpose times Theta; F=Theta J/Theta with fixed reference Theta_dot=0. No metric optimization.

Scan covers -180..180 deg at 0.001 deg spacing, partitioned at all refined current switches. Each interval endpoint is refined by a scalar root solve or is an exact switching/domain endpoint. Numerical certificate tolerances are 1e-9, matching baseline. Exact switches have only one-sided derivatives and no ordinary smooth certificate. Adjacent branches are not joined into a switched-system proof.

## Results

Nominal gamma=13.8107589411 deg, Pe=1, current=1.01530715382 pu, current margin=0.184692846183 pu, u_required=1.20360014576 pu. Inactive=1, voltage-valid=1, local-stable=1, contraction-certified=0. Eigenvalues -0.171259194103 and -3.50521139413; max eig SymF=20.5140979568. Nominal remains unchanged.

Current-limit switching angles (deg): -16.3409260376 16.3409260376 

Voltage-capability boundaries (deg): none in scanned domain

| Limited certified region | Lower deg | Upper deg | Mode | Lower cause | Upper cause |
|---|---:|---:|---|---|---|
| None | | | | | |

Baseline region coverage (widths in degrees; same target-angle domains, limited equilibria have new Pstar):

| Region | Baseline interval | Current feasible | Current limited | Limited voltage valid | Limited locally stable | Limited certified | Lost certified width |
|---|---|---:|---:|---:|---:|---:|---:|
| R1 | 22.352378400..29.217483849 | 0.000000000 | 6.865105449 | 6.865105449 | 0.000000000 | 0.000000000 | 6.865105449 |
| R2 | 50.705764461..84.790329411 | 0.000000000 | 34.084564950 | 34.084564950 | 0.000000000 | 0.000000000 | 34.084564950 |

Limiter-active locally stable interval width=0.747958727 deg; limiter-active voltage-valid certified width=0 deg. Active does not automatically imply unstable. See all_classification_intervals.csv for exact intervals and current_limiter_regions.csv for boundary causes.

- Limiter-free R1 disappears as a limited branchwise contraction-certified angle interval. Current-feasible here means the original prescribed-voltage equilibrium does not require current limitation.
- Limiter-free R2 disappears as a limited branchwise contraction-certified angle interval. Current-feasible here means the original prescribed-voltage equilibrium does not require current limitation.

All electrical network, E target, measurement mapping, tuning and metric methodology are identical; changes in power, equilibrium family, derivatives and certificates come from the ideal current projection and its necessary voltage reconstruction. Voltage-capability exclusion is separately reported, not voltage-saturation dynamics.

## Trajectories and optional sags

+1 deg is our test choice; the narrow stable limiter-active import branch uses +0.1 deg to stay on that branch initially. The neighboring initial angle is always a further +0.01 deg. Baseline/limited trajectory pairs share target equilibrium angle but generally have DIFFERENT equilibrium Pstar; tables state both. For metric comparisons the SAME limited-reference Theta is used for both curves, with distance to their original common target-angle equilibrium. Invalid voltage paths remain mathematical diagnostics and are flagged.

Optional nominal sag comparisons use identical Pstar=1 and initial equilibrium, retained voltages 0.5/0.3/0.1, duration 1, with zero angle jump. All requested electrical quantities are retained in MAT. Only nominal is used because old R1/R2 powers are not shared feasible equilibria. These representative sag comparisons are not a full duration/magnitude boundary sweep.

## Numerical verification

- inactiveBaseline: 2.220446e-16
- powerContinuity: 3.056444e-12
- voltageContinuity: 7.1400444e-13
- currentLimit: 2.220446e-16
- Pe_gamma: 2.7862068e-10
- Pe_U: 1.8145861e-09
- Pe_gammaU: 1.0893371e-08
- projectionDerivative: 1.0889795e-10
- Jacobian: 5.7849725e-11
- matrixFormula: 8.8817842e-16
- matrixSymmetry: 0
- networkResidual: 9.4246426e-16
- equilibriumResidual: 0
- nominalMatch: 7.1764816e-13

Analytical formulation is ready as an ideal benchmark for later EMT validation within its stated voltage-capability domain. It omits current/voltage PI states, delayed freeze and release behavior, measurement/Q-V/PLL dynamics and physical current overshoot. A smooth-branch test is not a global switched-system certificate. No EMT validation was started.
