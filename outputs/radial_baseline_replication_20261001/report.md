# ConstantRatioLimiter: baseline-style numerical disturbance validation

This update evaluates 12 initial angles at 50% and 70% source-voltage reduction (24 experiments). The original A-E points remain unchanged; F=16, G=17, H=23, I=24.2, J=24.55, K=32 and L=34.2 degrees are supplementary numerical cases. All figures display ONLY A-E. Numerical outcomes are listed below and are not formal fault-to-recovery certificates.

## Baseline replication audit

Before creating any figures, inspected `baseline_gfm_contraction_SCR5_XR10.m` and all four corresponding PNG exports in `baseline_SCR5_XR10`: matrix_eigenvalues, 07_representative_trajectories, figure5A_operating_points and figure5B_temporal_gamma_with_sag_examples.

| Original methodology | Limiter-enabled replication |
|---|---|
| Figure 3: three stacked spectra, J/M/Sym(F), zero lines and blue boundary guides | Same quantities and arrangement, with the full Mdot correction; fixed-Pstar state slice at nu=1, not an equilibrium family |
| Figure 4: phase plane above L(t), one case per column | Two rows x five cases for each sag separately; gamma unwrapped and nu-1; equilibrium/initial/clearing/final markers |
| Figure 5A: horizontal region bars and angle markers | Same bar-map layout for the five initial angles; zoomed to 10-42 degrees for readable labels |
| Figure 5B: full-width temporal top panel, four lower examples, early-time inset | Same 3x2 tiled structure; separate 50%/70% exports; top shows A-E, lower panels A/C/D/E, with B present in the top and its own Figure 4 column |

Case colours remain A blue, B orange, C yellow, D purple, E green throughout. Baseline figure numbering is retained with ConstantRatioLimiter suffixes. Supplementary frequency/current figures do not replace these four figure types. Every figure is exported as 300-dpi PNG, editable FIG and vector PDF. All nine PDF files were inspected: one page each, zero embedded raster images.

## Necessary mathematical and timing adaptations

The earlier baseline varied Pstar=Pe(gamma0) for different equilibria and used a +1-degree trajectory test. Here Pstar=1 and initial nu=1 for all 12 requested angles, without an additional angle offset. Only B is equilibrium. Figure 3 is explicitly labelled a fixed-Pstar state slice; its Jacobian eigenvalues away from B are not equilibrium-stability classifications. M and Sym(F) are left blank outside their validated domains and split at switching. J is shown over the original 0-180-degree angular axis as a mathematical reduced-model diagnostic, not an assertion of feasibility or stability of all those states.

The 30% reduction was removed. The two retained source factors are now .50 and .30, corresponding to 50% and 70% reductions. Timing and integration settings are unchanged from the previous ConstantRatioLimiter study.

The old baseline used fault onset 0, duration 1 in normalized analytical time, adaptive long horizons and gamma_dot=nu-1. That system cannot consistently be converted by relabelling its time axis. The authorized fallback is therefore used: fault-on 0 s, clearing .100 s, end 8 s. Initial angle displacement is specified at 0-minus, before the source jump at 0-plus; no artificial prefault dwell or clearing reset is added. Duplicate time rows at onset and clearing preserve both electrical sides.

Physical-time dynamics are gamma_dot=100*pi*(nu-1) and nu_dot=((1-Pe)/nu-25*(nu-1))/6.8. Numerical settings follow the compatible physical-time benchmark: ode45, RelTol=1e-9, AbsTol=[1e-11,1e-12], MaxStep=.002 s during fault and .005 s after clearing. Every experiment is independently re-integrated with RelTol=1e-11, AbsTol=[1e-13,1e-14] and half MaxStep. Saved outputs use 501 fault points, 6001 postfault points and exact detected event times. Recovery retains baseline tail criteria: last 10% of the 8-s horizon, |angle error|<.1 degree, |nu-1|<1e-4, with final |Pe-1|<1e-4.

## Physical model and validity

External SCR5/XR10, physical Xeq=.350, added LV reactor .0858989549982 mH, Vset=1, Qstar=0, H=3.4, K=20, Kd=5, radial Imax=1.2, radial umax=1.3 are unchanged. The saved physical network is reconstructed for each source magnitude with the reactor included once. Passive coefficients and Zeq stay fixed; Vth changes consistently with the source. Measured PCC voltage is not simply multiplied by a sag factor. Q-V closure remains Eref=1-.1*Qgrid algebraically; actual filter and PI states are omitted.

The voltage-command limiter is a hard validity condition of the ideal-current-tracking benchmark, not simulated by silently clipping voltage while keeping ideal current tracking. The implementation stops at |u_required|=1.3 or the nu=.5 guard, including a check immediately after source restoration. Stopped trajectories are not assigned recovery beyond their validity interval.

### Pre-fault initial quantities

| Case | gamma0 (deg) | Raw I | Limited I | Pe (pu) | Required u | Limiter |
|---|---:|---:|---:|---:|---:|---|
| A | 14.000000 | 0.692191 | 0.692191 | 0.683938 | 1.140583 | inactive |
| B | 20.633905 | 1.017813 | 1.017813 | 1.000000 | 1.210040 | inactive |
| C | 24.373430 | 1.200000 | 1.200000 | 1.173090 | 1.248880 | boundary |
| D | 27.000000 | 1.328609 | 1.200000 | 1.166343 | 1.240434 | active |
| E | 36.000000 | 1.766110 | 1.200000 | 1.138371 | 1.209140 | active |
| F | 16.000000 | 0.790630 | 0.790630 | 0.780233 | 1.161558 | inactive |
| G | 17.000000 | 0.839770 | 0.839770 | 0.828083 | 1.172040 | inactive |
| H | 23.000000 | 1.133223 | 1.133223 | 1.110020 | 1.234656 | inactive |
| I | 24.200000 | 1.191577 | 1.191577 | 1.165160 | 1.247087 | inactive |
| J | 24.550000 | 1.208655 | 1.200000 | 1.172657 | 1.248324 | active |
| K | 32.000000 | 1.572378 | 1.200000 | 1.151721 | 1.223449 | active |
| L | 34.200000 | 1.679091 | 1.200000 | 1.144559 | 1.215651 | active |

Case C is checked with both smooth equations: the Q-V root and Pe agree at the exact switching root. Dynamics use the physically continuous min(1,Imax/|q|) radial law; the reporting tolerance never changes the dynamics. D and E are already limiter-active before fault application, so their event logs need not contain a new activation at t=0; the initial-condition table records this explicitly.

## 24-case comparison

All previous CSV columns are retained, including angle/frequency extrema, equilibrium errors, membership, exit/re-entry times and voltage feasibility. The table below summarizes all 24 runs; the following table adds extrema and equilibrium errors. Supplementary F-L cases appear only in data and tables, never in figures.

| Case | Initial angle | Sag depth | Clearing gamma | Clearing nu | Limiter? / active duration (s) | Clearing domain | All postfault samples inside? | Voltage feasible? | Numerical recovery | Final gamma | Final nu |
|---|---:|---:|---:|---:|---|---|---|---|---|---:|---:|
| A | 14.000000 | 50% | 21.210337 | 1.007291538 | Y / 0.539155 | outside | N | Y | Y | 20.633909322 | 0.9999999956 |
| A | 14.000000 | 70% | 23.433575 | 1.009693410 | Y / 0.812101 | outside | N | Y | Y | 20.633902705 | 1.0000000089 |
| B | 20.633905 | 50% | 26.819433 | 1.006312097 | Y / 0.599506 | outside | N | Y | Y | 20.633912042 | 0.9999999961 |
| B | 20.633905 | 70% | 29.506845 | 1.009139181 | Y / 0.985673 | outside | N | Y | Y | 20.633886426 | 0.9999999970 |
| C | 24.373430 | 50% | 30.107639 | 1.005878050 | Y / 0.645106 | outside | N | Y | Y | 20.633914754 | 0.9999999981 |
| C | 24.373430 | 70% | 32.964313 | 1.008861442 | Y / 1.137812 | outside | N | Y | Y | 20.633917879 | 0.9999999843 |
| D | 27.000000 | 50% | 32.466668 | 1.005619936 | Y / 0.684645 | outside | N | Y | Y | 20.633915456 | 1.0000000007 |
| D | 27.000000 | 70% | 35.408465 | 1.008682230 | Y / 1.304033 | outside | N | Y | Y | 20.633931454 | 1.0000000126 |
| E | 36.000000 | 50% | 40.814825 | 1.004989897 | Y / 0.895418 | outside | N | Y | Y | 20.633890628 | 1.0000000072 |
| E | 36.000000 | 70% | 43.884698 | 1.008170035 | Y / 7.363016 | outside | N | Y | N | 4590.334361219 | 1.0422831356 |
| F | 16.000000 | 50% | 22.869892 | 1.006967000 | Y / 0.549924 | outside | N | Y | Y | 20.633909121 | 0.9999999956 |
| F | 16.000000 | 70% | 25.256940 | 1.009518476 | Y / 0.858743 | outside | N | Y | Y | 20.633896085 | 1.0000000088 |
| G | 17.000000 | 50% | 23.709976 | 1.006814293 | Y / 0.560769 | outside | N | Y | Y | 20.633909432 | 0.9999999956 |
| G | 17.000000 | 70% | 26.171005 | 1.009433476 | Y / 0.882491 | outside | N | Y | Y | 20.633892494 | 1.0000000078 |
| H | 23.000000 | 50% | 28.890006 | 1.006028086 | Y / 0.627016 | outside | N | Y | Y | 20.633913867 | 0.9999999972 |
| H | 23.000000 | 70% | 31.691471 | 1.008960404 | Y / 1.074110 | outside | N | Y | Y | 20.633900309 | 0.9999999862 |
| I | 24.200000 | 50% | 29.953264 | 1.005896411 | Y / 0.642729 | outside | N | Y | Y | 20.633914655 | 0.9999999980 |
| I | 24.200000 | 70% | 32.803389 | 1.008873741 | Y / 1.129114 | outside | N | Y | Y | 20.633915408 | 0.9999999841 |
| J | 24.550000 | 50% | 30.264990 | 1.005859527 | Y / 0.647553 | outside | N | Y | Y | 20.633914850 | 0.9999999982 |
| J | 24.550000 | 70% | 33.128208 | 1.008848979 | Y / 1.146890 | outside | N | Y | Y | 20.633920412 | 0.9999999847 |
| K | 32.000000 | 50% | 37.058096 | 1.005224947 | Y / 0.783025 | outside | N | Y | Y | 20.633908048 | 1.0000000078 |
| K | 32.000000 | 70% | 40.098030 | 1.008378138 | Y / 7.468812 | outside | N | Y | N | 3824.902353641 | 1.0356611748 |
| L | 34.200000 | 50% | 39.115852 | 1.005087430 | Y / 0.839853 | outside | N | Y | Y | 20.633899160 | 1.0000000092 |
| L | 34.200000 | 70% | 42.176856 | 1.008259826 | Y / 7.363016 | outside | N | Y | N | 4457.944444931 | 1.0299146742 |

| Case/sag | min gamma | max gamma | min f (Hz) | max f (Hz) | final gamma error (deg) | final nu error | max required u |
|---|---:|---:|---:|---:|---:|---:|---:|
| A/50% | 14.000000 | 33.705693 | 49.803512 | 50.364577 | 4.44687e-06 | -4.35326e-09 | 1.24888037 |
| A/70% | 14.000000 | 43.379412 | 49.769642 | 50.484670 | -2.17074e-06 | 8.85035e-09 | 1.24888037 |
| B/50% | 15.568473 | 36.598577 | 49.789494 | 50.315605 | 7.1667e-06 | -3.87728e-09 | 1.24888037 |
| B/70% | 14.898022 | 49.368516 | 49.761491 | 50.456959 | -1.84496e-05 | -3.01154e-09 | 1.24888037 |
| C/50% | 15.345115 | 39.179029 | 49.780170 | 50.293902 | 9.87832e-06 | -1.89864e-09 | 1.24888037 |
| C/70% | 14.838674 | 53.285197 | 49.759016 | 50.443072 | 1.30031e-05 | -1.5671e-08 | 1.24888037 |
| D/50% | 15.209118 | 41.197673 | 49.774493 | 50.280997 | 1.05801e-05 | 7.35106e-10 | 1.24888037 |
| D/70% | 14.815742 | 56.404022 | 49.758060 | 50.434111 | 2.65781e-05 | 1.26371e-08 | 1.24888037 |
| E/50% | 14.897617 | 49.385841 | 49.761477 | 50.249495 | -1.42476e-05 | 7.1627e-09 | 1.24888037 |
| E/70% | 36.000000 | 4590.334361 | 50.000000 | 52.616709 | 4569.7 | 0.0422831 | 1.24888037 |
| F/50% | 15.853791 | 34.088137 | 49.801401 | 50.348350 | 4.24535e-06 | -4.44348e-09 | 1.24888037 |
| F/70% | 15.022566 | 45.046944 | 49.766700 | 50.475924 | -8.7902e-06 | 8.7706e-09 | 1.24888037 |
| G/50% | 15.804645 | 34.475553 | 49.799344 | 50.340715 | 4.55639e-06 | -4.44271e-09 | 1.24888037 |
| G/70% | 14.990458 | 45.940692 | 49.765357 | 50.471674 | -1.23818e-05 | 7.82188e-09 | 1.24888037 |
| H/50% | 15.423490 | 38.190446 | 49.783438 | 50.301404 | 8.99165e-06 | -2.8413e-09 | 1.24888037 |
| H/70% | 14.856636 | 51.789557 | 49.759759 | 50.448020 | -4.56683e-06 | -1.38023e-08 | 1.24888037 |
| I/50% | 15.354765 | 39.051637 | 49.780571 | 50.294821 | 9.77941e-06 | -2.03335e-09 | 1.24888037 |
| I/70% | 14.840758 | 53.092030 | 49.759097 | 50.443687 | 1.05326e-05 | -1.58657e-08 | 1.24888037 |
| J/50% | 15.335420 | 39.309502 | 49.779761 | 50.292976 | 9.9742e-06 | -1.75667e-09 | 1.24888037 |
| J/70% | 14.836691 | 53.483254 | 49.758930 | 50.442449 | 1.55365e-05 | -1.53017e-08 | 1.24888037 |
| K/50% | 15.006121 | 45.495736 | 49.766003 | 50.261247 | 3.17218e-06 | 7.75835e-09 | 1.24888037 |
| K/70% | 32.000000 | 3824.902354 | 50.000000 | 52.616737 | 3804.27 | 0.0356612 | 1.24888037 |
| L/50% | 14.940590 | 47.581690 | 49.763271 | 50.254371 | -5.71548e-06 | 9.23173e-09 | 1.24888037 |
| L/70% | 34.200000 | 4457.944445 | 50.000000 | 52.616735 | 4437.31 | 0.0299147 | 1.24888037 |

**Observed outcomes:** 21/24 numerical recoveries; 24/24 voltage-feasible trajectories.

Non-recovering rows: E (36 deg), 70% sag, K (32 deg), 70% sag, L (34.2 deg), 70% sag. These trajectories show sustained angle drift within the 8-s horizon, not a voltage-validity stop. No critical sag or exact stability boundary is claimed.

Complete postfault sampled containment: 0/24. Clearing-state membership, Omega entry and terminal entry remain distinct from recovery.

The unchanged absolute state-comparison threshold of 2e-6 was exceeded for K/70%, L/70%, E/70%. Integration tolerances were NOT relaxed; flags are retained in both verification and comparison CSVs. Both integrations agree on non-recovery and voltage feasibility in these drifting cases. Exact final-angle digits should not be interpreted beyond the reported integration discrepancy.

For A/50%, the event sequence explicitly captures re-entry: activation at the initial source jump, release at clearing (.1 s), reactivation at .125789455851 s and release at .564944752915 s. Total active duration is .539155297064 s. Durations are summed over event-partitioned solution intervals, not inferred from a coarse activity plot.

## Metrics: what the displacement plots mean

Full-precision existing coefficient files are loaded without optimization. For x_d=[gamma_eq,1], the baseline algebraic norm is retained locally as

L_i(x)=||Theta_i(x)(x-x_d)||=sqrt((x-x_d)' M_i(x)(x-x_d)), Theta_i=chol(M_i).

With a state-dependent metric this is a **local weighted chord/displacement**, not a geodesic distance or an independently proven decaying Lyapunov function. In particular, the active domain does not contain x_d, and contraction of tangent vectors does not automatically imply monotonic decay of this finite chord. Fault-on and outside-domain samples are blank. Active and inactive curves use different line styles and NaN breaks; they are never silently joined into one global distance. One-sided chord values at compatible postfault switch events are saved separately and marked on Figure 4. No rescaling is used to hide metric jumps.

The prior validated quadratic-form gains are <=256 for active-to-inactive and <=2 in reverse on the full compatible event strip; on the smaller Omega strip, active-to-inactive gain is <=4. Square roots apply to lengths. Event files labelled activeToInactive report that directional comparison; outside the compatible event strip no validated gain is asserted. These bounds are conditional on the previously verified nominal-voltage/transversal-switch assumptions and do not cover source jumps.

At each certified postfault sample, the full Qc=J'M+MJ+Mdot+2*lambda*M is evaluated with lambda1=1 or lambda2=.1. Figure 3 uses Sym(F)=.5*Theta^(-T)*(J'M+MJ+Mdot)*Theta^(-1), equivalent to the symmetric transformed Jacobian including Theta_dot. Dropping Mdot would not reproduce the certified analysis.

### Three distinct conclusions

1. **Numerical recovery:** assessed separately for each row, only within a valid reduced-model trajectory and the unchanged final-tail criteria.
2. **Contraction validity on portions:** the appropriate full inequality is checked at every saved sample inside its certified branch/domain. Complete sampled containment is reported separately for each row.
3. **Formal fault-to-recovery guarantee:** not newly established for any of these numerical disturbance trajectories. The earlier Omega certificate remains valid, but sampled/interpolated Omega entry is not a validated ODE enclosure. Entry times are explicitly labelled sampled. No outside-domain point is automatically classified unstable.

Omega membership uses the prior energy formula with numerical quadrature/interpolation; terminal membership uses the saved Lyapunov ellipsoid. Neither calculation is promoted to a rigorous enclosure of the integrated trajectory. No new metric optimization, interval certificate or EMT validation is claimed.

## Verification and files

- Maximum tighter-run final-state difference: 2.61e-05.
- Maximum limiter-event time difference: 2.94e-06 s; event counts match.
- Independent direct nodal power error: 8.88e-16 pu.
- Physical gamma_dot check error: 9.97e-07 rad/s.
- Maximum limiter-event residual: 4.41e-12 pu; clearing state discontinuity is zero.
- Highest required command across all runs: 1.248880365 pu, below 1.3.
- All current magnitudes obey 1.2 pu, and full metric inequalities pass where applicable. Baseline and certified metric hashes are unchanged (including the 59-file protected baseline manifest).
- `results.mat`: all 24 trajectories, solver solutions, complete settings, network, metrics and initial summary. `audit_results.mat` contains the augmented summary and independent checks.
- `case_*_trajectory.csv`: all requested state/electrical/limiter/metric/membership histories. Event, domain-transition and metric-transition CSVs accompany each run.
- `comparison_summary.csv`, `initial_conditions.csv`, `integration_verification.csv`, `model_consistency_audit.csv`, `settings.json` and `export_verification.json` provide machine-readable results.
- `run_study.m` implements the experiments, calls `plot_baseline_figures.m`, then `audit_study.m`. The same scripts were modified directly. To reproduce figures alone, call plot_baseline_figures with this result directory. `write_report.py` produces this report and checks exports/hashes. No new metric optimization or interval certification was performed.
- The complete prior 10-case directory, including its scripts and all outputs, was copied and hash-checked before the changes: `../radial_baseline_replication_20261001_archive_10cases_20261001_184741`. The current script refuses to overwrite results.mat; archive results first for another run.

## Figure links

- Figure3_ConstantRatioLimiter: [PNG](Figure3_ConstantRatioLimiter.png) | [vector PDF](Figure3_ConstantRatioLimiter.pdf)
- Figure4_ConstantRatioLimiter_Sag50: [PNG](Figure4_ConstantRatioLimiter_Sag50.png) | [vector PDF](Figure4_ConstantRatioLimiter_Sag50.pdf)
- Figure4_ConstantRatioLimiter_Sag70: [PNG](Figure4_ConstantRatioLimiter_Sag70.png) | [vector PDF](Figure4_ConstantRatioLimiter_Sag70.pdf)
- Figure5A_ConstantRatioLimiter: [PNG](Figure5A_ConstantRatioLimiter.png) | [vector PDF](Figure5A_ConstantRatioLimiter.pdf)
- Figure5B_ConstantRatioLimiter_Sag50: [PNG](Figure5B_ConstantRatioLimiter_Sag50.png) | [vector PDF](Figure5B_ConstantRatioLimiter_Sag50.pdf)
- Figure5B_ConstantRatioLimiter_Sag70: [PNG](Figure5B_ConstantRatioLimiter_Sag70.png) | [vector PDF](Figure5B_ConstantRatioLimiter_Sag70.pdf)
- Supplement_Frequency_ConstantRatioLimiter: [PNG](Supplement_Frequency_ConstantRatioLimiter.png) | [vector PDF](Supplement_Frequency_ConstantRatioLimiter.pdf)
- Supplement_Current_ConstantRatioLimiter_Sag50: [PNG](Supplement_Current_ConstantRatioLimiter_Sag50.png) | [vector PDF](Supplement_Current_ConstantRatioLimiter_Sag50.pdf)
- Supplement_Current_ConstantRatioLimiter_Sag70: [PNG](Supplement_Current_ConstantRatioLimiter_Sag70.png) | [vector PDF](Supplement_Current_ConstantRatioLimiter_Sag70.pdf)
