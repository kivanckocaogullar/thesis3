# ConstantRatioLimiter: baseline-style numerical disturbance validation

**All 10 experiments recover numerically to the same nominal equilibrium.** Current limitation is encountered in every experiment, and no required voltage command exceeds 1.3 pu. Only C/30% and D/30% remain within the certified piecewise domains at every saved postfault point. All ten numerically enter Omega and the terminal ellipsoid. These observations do not establish rigorous fault-to-recovery certificates.

## Baseline replication audit

Before creating any figures, inspected `baseline_gfm_contraction_SCR5_XR10.m` and all four corresponding PNG exports in `baseline_SCR5_XR10`: matrix_eigenvalues, 07_representative_trajectories, figure5A_operating_points and figure5B_temporal_gamma_with_sag_examples.

| Original methodology | Limiter-enabled replication |
|---|---|
| Figure 3: three stacked spectra, J/M/Sym(F), zero lines and blue boundary guides | Same quantities and arrangement, with the full Mdot correction; fixed-Pstar state slice at nu=1, not an equilibrium family |
| Figure 4: phase plane above L(t), one case per column | Two rows x five cases for each sag separately; gamma unwrapped and nu-1; equilibrium/initial/clearing/final markers |
| Figure 5A: horizontal region bars and angle markers | Same bar-map layout for the five initial angles; zoomed to 10-42 degrees for readable labels |
| Figure 5B: full-width temporal top panel, four lower examples, early-time inset | Same 3x2 tiled structure; separate 30%/50% exports; top shows all five cases, lower panels A/C/D/E, with B present in the top and its own Figure 4 column |

Case colours remain A blue, B orange, C yellow, D purple, E green throughout. Baseline figure numbering is retained with ConstantRatioLimiter suffixes. Supplementary frequency/current figures do not replace these four figure types. Every figure is exported as 300-dpi PNG, editable FIG and vector PDF. All nine PDF files were inspected: one page each, zero embedded raster images.

## Necessary mathematical and timing adaptations

The earlier baseline varied Pstar=Pe(gamma0) for different equilibria and used a +1-degree trajectory test. Here Pstar=1 for every state and experiment, with the exact requested five initial angles and no extra +1-degree offset. Only B is equilibrium. Figure 3 is explicitly labelled a fixed-Pstar state slice; its Jacobian eigenvalues away from B are not equilibrium-stability classifications. M and Sym(F) are left blank outside their validated domains and split at switching. J is shown over the original 0-180-degree angular axis as a mathematical reduced-model diagnostic, not an assertion of feasibility or stability of all those states.

Old `Usag=.3` meant retained .3, i.e. a 70% reduction. The new 30% and 50% reductions use retained source factors .70 and .50. They are not 0.3 and 0.5 retained-voltage tests.

The old baseline used fault onset 0, duration 1 in normalized analytical time, adaptive long horizons and gamma_dot=nu-1. That system cannot consistently be converted by relabelling its time axis. The authorized fallback is therefore used: fault-on 0 s, clearing .100 s, end 8 s. Initial angle displacement is specified at 0-minus, before the source jump at 0-plus; no artificial prefault dwell or clearing reset is added. Duplicate time rows at onset and clearing preserve both electrical sides.

Physical-time dynamics are gamma_dot=100*pi*(nu-1) and nu_dot=((1-Pe)/nu-25*(nu-1))/6.8. Numerical settings follow the compatible physical-time benchmark: ode45, RelTol=1e-9, AbsTol=[1e-11,1e-12], MaxStep=.002 s during fault and .005 s after clearing. Every experiment is independently re-integrated with RelTol=1e-11, AbsTol=[1e-13,1e-14] and half MaxStep. Saved outputs use 501 fault points, 6001 postfault points and exact detected event times. Recovery retains baseline tail criteria: last 10% of the 8-s horizon, |angle error|<.1 degree, |nu-1|<1e-4, with final |Pe-1|<1e-4.

## Physical model and validity

External SCR5/XR10, physical Xeq=.350, added LV reactor .0858989549982 mH, Vset=1, Qstar=0, H=3.4, K=20, Kd=5, radial Imax=1.2, radial umax=1.3 are unchanged. The saved physical network is reconstructed for each source magnitude with the reactor included once. Passive coefficients and Zeq stay fixed; Vth changes consistently with the source. Measured PCC voltage is not simply multiplied by a sag factor. Q-V closure remains Eref=1-.1*Qgrid algebraically; actual filter and PI states are omitted.

The voltage-command limiter is a hard validity condition of the ideal-current-tracking benchmark, not simulated by silently clipping voltage while keeping ideal current tracking. The implementation stops at |u_required|=1.3 or the nu=.5 guard, including a check immediately after source restoration. None of these ten tests triggers either stop.

### Pre-fault initial quantities

| Case | gamma0 (deg) | Raw I | Limited I | Pe (pu) | Required u | Limiter |
|---|---:|---:|---:|---:|---:|---|
| A | 14.000000 | 0.692191 | 0.692191 | 0.683938 | 1.140583 | inactive |
| B | 20.633905 | 1.017813 | 1.017813 | 1.000000 | 1.210040 | inactive |
| C | 24.373430 | 1.200000 | 1.200000 | 1.173090 | 1.248880 | boundary |
| D | 27.000000 | 1.328609 | 1.200000 | 1.166343 | 1.240434 | active |
| E | 36.000000 | 1.766110 | 1.200000 | 1.138371 | 1.209140 | active |

Case C is checked with both smooth equations: the Q-V root and Pe agree at the exact switching root. Dynamics use the physically continuous min(1,Imax/|q|) radial law; the reporting tolerance never changes the dynamics. D and E are already limiter-active before fault application, so their event logs need not contain a new activation at t=0; the initial-condition table records this explicitly.

## Ten-case comparison

All cases: limiter encountered = Y, voltage feasibility = Y, sampled Omega entry = Y, sampled terminal entry = Y, numerical recovery = Y, model-validity violation = none. Final angles are within 0.000015 degrees of 20.633904875622; final nu errors are below 8e-9. Detailed final values, extrema, errors, entries/exits and margins are in comparison_summary.csv.

| Case | Initial angle | Sag depth | Clearing gamma | Clearing nu | Active duration (s) | Clearing domain | All postfault samples inside? | Final gamma | Final nu |
|---|---:|---:|---:|---:|---:|---|---|---:|---:|
| A | 14.000000 | 30% | 19.457642 | 1.005299276 | 0.263286 | outside | N | 20.633901533 | 0.9999999983 |
| A | 14.000000 | 50% | 21.210337 | 1.007291538 | 0.539155 | outside | N | 20.633909322 | 0.9999999956 |
| B | 20.633905 | 30% | 23.819639 | 1.003101588 | 0.196753 | outside | N | 20.633901642 | 1.0000000010 |
| B | 20.633905 | 50% | 26.819433 | 1.006312097 | 0.599506 | outside | N | 20.633912042 | 0.9999999961 |
| C | 24.373430 | 30% | 26.418826 | 1.002112510 | 0.310279 | R2 | Y | 20.633901393 | 1.0000000013 |
| C | 24.373430 | 50% | 30.107639 | 1.005878050 | 0.645106 | outside | N | 20.633914754 | 0.9999999981 |
| D | 27.000000 | 30% | 28.882353 | 1.001951088 | 0.355266 | R2 | Y | 20.633900190 | 1.0000000006 |
| D | 27.000000 | 50% | 32.466668 | 1.005619936 | 0.684645 | outside | N | 20.633915456 | 1.0000000007 |
| E | 36.000000 | 30% | 37.617553 | 1.001690325 | 0.512893 | outside | N | 20.633904687 | 0.9999999953 |
| E | 36.000000 | 50% | 40.814825 | 1.004989897 | 0.895418 | outside | N | 20.633890628 | 1.0000000072 |

For A/50%, the event sequence explicitly captures re-entry: activation at the initial source jump, release at clearing (.1 s), reactivation at .125789455851 s and release at .564944752915 s. Total active duration is .539155297064 s. Durations are summed over event-partitioned solution intervals, not inferred from a coarse activity plot.

## Metrics: what the displacement plots mean

Full-precision existing coefficient files are loaded without optimization. For x_d=[gamma_eq,1], the baseline algebraic norm is retained locally as

L_i(x)=||Theta_i(x)(x-x_d)||=sqrt((x-x_d)' M_i(x)(x-x_d)), Theta_i=chol(M_i).

With a state-dependent metric this is a **local weighted chord/displacement**, not a geodesic distance or an independently proven decaying Lyapunov function. In particular, the active domain does not contain x_d, and contraction of tangent vectors does not automatically imply monotonic decay of this finite chord. Fault-on and outside-domain samples are blank. Active and inactive curves use different line styles and NaN breaks; they are never silently joined into one global distance. One-sided chord values at compatible postfault switch events are saved separately and marked on Figure 4. No rescaling is used to hide metric jumps.

The prior validated quadratic-form gains are <=256 for active-to-inactive and <=2 in reverse on the full compatible event strip; on the smaller Omega strip, active-to-inactive gain is <=4. Square roots apply to lengths. Event files labelled activeToInactive report that directional comparison; outside the compatible event strip no validated gain is asserted. These bounds are conditional on the previously verified nominal-voltage/transversal-switch assumptions and do not cover source jumps.

At each certified postfault sample, the full Qc=J'M+MJ+Mdot+2*lambda*M is evaluated with lambda1=1 or lambda2=.1. Figure 3 uses Sym(F)=.5*Theta^(-T)*(J'M+MJ+Mdot)*Theta^(-1), equivalent to the symmetric transformed Jacobian including Theta_dot. Dropping Mdot would not reproduce the certified analysis.

### Three distinct conclusions

1. **Numerical recovery:** observed in all ten simulations.
2. **Contraction validity on portions:** the appropriate full inequality passes at every saved sample inside its certified branch/domain; only C/30 and D/30 have all saved postfault samples inside the union.
3. **Formal fault-to-recovery guarantee:** not newly established for any of these numerical disturbance trajectories. The earlier Omega certificate remains valid, but sampled/interpolated Omega entry is not a validated ODE enclosure. Entry times are explicitly labelled sampled. No outside-domain point is automatically classified unstable.

Omega membership uses the prior energy formula with numerical quadrature/interpolation; terminal membership uses the saved Lyapunov ellipsoid. Neither calculation is promoted to a rigorous enclosure of the integrated trajectory. No new metric optimization, interval certificate or EMT validation is claimed.

## Verification and files

- Maximum tighter-run final-state difference: 5.7e-13.
- Maximum limiter-event time difference: 1.49e-06 s; event counts match.
- Independent direct nodal power error: 8.88e-16 pu.
- Physical gamma_dot check error: 1.82e-08 rad/s.
- Maximum limiter-event residual: 1.87e-14 pu; clearing state discontinuity is zero.
- Highest required command across all runs: 1.248880365 pu, below 1.3.
- All current magnitudes obey 1.2 pu, and full metric inequalities pass where applicable. Baseline and certified metric hashes are unchanged (including the 59-file protected baseline manifest).
- `results.mat`: all ten trajectories, solver solutions, complete settings, network, metrics and initial summary. `audit_results.mat` contains the augmented summary and independent checks.
- `case_*_trajectory.csv`: all requested state/electrical/limiter/metric/membership histories. Event, domain-transition and metric-transition CSVs accompany each run.
- `comparison_summary.csv`, `initial_conditions.csv`, `integration_verification.csv`, `model_consistency_audit.csv`, `settings.json` and `export_verification.json` provide machine-readable results.
- `run_study.m` implements the experiments, calls `plot_baseline_figures.m`, then `audit_study.m`. To rerun without overwriting, copy these scripts into a new sibling result directory under outputs; the main script refuses to overwrite an existing results.mat. To reproduce figures alone, call plot_baseline_figures with this result directory. `write_report.py` produces this report and checks exports/hashes.

## Figure links

- Figure3_ConstantRatioLimiter: [PNG](../radial_baseline_replication_20261001/Figure3_ConstantRatioLimiter.png) | [vector PDF](../radial_baseline_replication_20261001/Figure3_ConstantRatioLimiter.pdf)
- Figure4_ConstantRatioLimiter_Sag30: [PNG](../radial_baseline_replication_20261001/Figure4_ConstantRatioLimiter_Sag30.png) | [vector PDF](../radial_baseline_replication_20261001/Figure4_ConstantRatioLimiter_Sag30.pdf)
- Figure4_ConstantRatioLimiter_Sag50: [PNG](../radial_baseline_replication_20261001/Figure4_ConstantRatioLimiter_Sag50.png) | [vector PDF](../radial_baseline_replication_20261001/Figure4_ConstantRatioLimiter_Sag50.pdf)
- Figure5A_ConstantRatioLimiter: [PNG](../radial_baseline_replication_20261001/Figure5A_ConstantRatioLimiter.png) | [vector PDF](../radial_baseline_replication_20261001/Figure5A_ConstantRatioLimiter.pdf)
- Figure5B_ConstantRatioLimiter_Sag30: [PNG](../radial_baseline_replication_20261001/Figure5B_ConstantRatioLimiter_Sag30.png) | [vector PDF](../radial_baseline_replication_20261001/Figure5B_ConstantRatioLimiter_Sag30.pdf)
- Figure5B_ConstantRatioLimiter_Sag50: [PNG](../radial_baseline_replication_20261001/Figure5B_ConstantRatioLimiter_Sag50.png) | [vector PDF](../radial_baseline_replication_20261001/Figure5B_ConstantRatioLimiter_Sag50.pdf)
- Supplement_Frequency_ConstantRatioLimiter: [PNG](../radial_baseline_replication_20261001/Supplement_Frequency_ConstantRatioLimiter.png) | [vector PDF](../radial_baseline_replication_20261001/Supplement_Frequency_ConstantRatioLimiter.pdf)
- Supplement_Current_ConstantRatioLimiter_Sag30: [PNG](../radial_baseline_replication_20261001/Supplement_Current_ConstantRatioLimiter_Sag30.png) | [vector PDF](../radial_baseline_replication_20261001/Supplement_Current_ConstantRatioLimiter_Sag30.pdf)
- Supplement_Current_ConstantRatioLimiter_Sag50: [PNG](../radial_baseline_replication_20261001/Supplement_Current_ConstantRatioLimiter_Sag50.png) | [vector PDF](../radial_baseline_replication_20261001/Supplement_Current_ConstantRatioLimiter_Sag50.pdf)
