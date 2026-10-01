from pathlib import Path
import csv,json,hashlib,re
from pypdf import PdfReader
from PIL import Image
p=Path(__file__).parent
read=lambda name:list(csv.DictReader((p/name).open(encoding='utf-8-sig')))
rows=read('comparison_summary.csv');initial=read('initial_conditions.csv');checks=read('integration_verification.csv');audit=read('model_consistency_audit.csv')
assert len(rows)==24 and len(initial)==12 and len(checks)==24 and len(audit)==24
assert len(list(p.glob('case_*_trajectory.csv')))==24
archive=p.parent/'radial_baseline_replication_20261001_archive_10cases_20261001_184741'
oldrows=list(csv.DictReader((archive/'comparison_summary.csv').open(encoding='utf-8-sig')))
assert set(oldrows[0]).issubset(rows[0])
truth=lambda x:x.lower() in ['true','1']
qa=[]
for f in sorted(p.glob('*.pdf')):
 reader=PdfReader(f);assert len(reader.pages)==1 and len(reader.pages[0].images)==0
 assert not re.search(r'\b[F-L]\s*:\s*\d',reader.pages[0].extract_text()),'Supplementary case label found in figure'
 png=f.with_suffix('.png');im=Image.open(png);assert min(im.size)>1000
 qa.append({'name':f.name,'pages':1,'embedded_raster_images':0,'png_pixels':im.size})
assert len(qa)==9
(p/'export_verification.json').write_text(json.dumps(qa,indent=2))
hashes=read('protected_before.csv')
for r in hashes:assert hashlib.sha256(Path(r['Path']).read_bytes()).hexdigest().upper()==r['Hash'].upper()
for r in csv.DictReader((p.parent/'current_limiter_SCR5_XR10/protected_before.csv').open(encoding='utf-8-sig')):
 assert hashlib.sha256(Path(r['Path']).read_bytes()).hexdigest().upper()==r['Hash'].upper()
r='''# ConstantRatioLimiter: baseline-style numerical disturbance validation

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
'''
for x in initial:r+=f"| {x['Case']} | {float(x['gamma0_deg']):.6f} | {float(x['rawCurrent']):.6f} | {float(x['limitedCurrent']):.6f} | {float(x['Pe']):.6f} | {float(x['requiredVoltage']):.6f} | {x['limiterState']} |\n"
r+='''
Case C is checked with both smooth equations: the Q-V root and Pe agree at the exact switching root. Dynamics use the physically continuous min(1,Imax/|q|) radial law; the reporting tolerance never changes the dynamics. D and E are already limiter-active before fault application, so their event logs need not contain a new activation at t=0; the initial-condition table records this explicitly.

## 24-case comparison

All previous CSV columns are retained, including angle/frequency extrema, equilibrium errors, membership, exit/re-entry times and voltage feasibility. The table below summarizes all 24 runs; the following table adds extrema and equilibrium errors. Supplementary F-L cases appear only in data and tables, never in figures.

| Case | Initial angle | Sag depth | Clearing gamma | Clearing nu | Limiter? / active duration (s) | Clearing domain | All postfault samples inside? | Voltage feasible? | Numerical recovery | Final gamma | Final nu |
|---|---:|---:|---:|---:|---|---|---|---|---|---:|---:|
'''
for x in rows:
 r+=f"| {x['Case']} | {float(x['gamma0_deg']):.6f} | {float(x['sagDepth_percent']):.0f}% | {float(x['clearGamma_deg']):.6f} | {float(x['clearNu']):.9f} | {'Y' if truth(x['limiterActivated']) else 'N'} / {float(x['activeDuration_s']):.6f} | {x['clearRegion']} | {'Y' if truth(x['postfaultContained_sampled']) else 'N'} | {'Y' if truth(x['voltageFeasible']) else 'N'} | {x['recovery']} | {float(x['finalGamma_deg']):.9f} | {float(x['finalNu']):.10f} |\n"
r+='\n| Case/sag | min gamma | max gamma | min f (Hz) | max f (Hz) | final gamma error (deg) | final nu error | max required u |\n|---|---:|---:|---:|---:|---:|---:|---:|\n'
for x in rows:r+=f"| {x['Case']}/{float(x['sagDepth_percent']):.0f}% | {float(x['minGamma_deg']):.6f} | {float(x['maxGamma_deg']):.6f} | {float(x['minFrequency_Hz']):.6f} | {float(x['maxFrequency_Hz']):.6f} | {float(x['finalGammaError_deg']):.6g} | {float(x['finalNuError']):.6g} | {float(x['maxRequiredVoltage']):.8f} |\n"
r+='\n**Observed outcomes:** '+str(sum(x['recovery']=='Y' for x in rows))+'/24 numerical recoveries; '+str(sum(truth(x['voltageFeasible']) for x in rows))+'/24 voltage-feasible trajectories.\n\n'
r+='Non-recovering rows: '+', '.join(x['Case']+' ('+x['gamma0_deg']+' deg), '+str(round(float(x['sagDepth_percent'])))+'% sag' for x in rows if x['recovery']!='Y')+'. These trajectories show sustained angle drift within the 8-s horizon, not a voltage-validity stop. No critical sag or exact stability boundary is claimed.\n\n'
r+='Complete postfault sampled containment: '+str(sum(truth(x['postfaultContained_sampled']) for x in rows))+'/24. Clearing-state membership, Omega entry and terminal entry remain distinct from recovery.\n\n'
r+='The unchanged absolute state-comparison threshold of 2e-6 was exceeded for '+', '.join(x['Case']+'/'+str(round(float(x['sagDepth_percent'])))+'%' for x in checks if not truth(x['withinOriginalStateCheck']))+'. Integration tolerances were NOT relaxed; flags are retained in both verification and comparison CSVs. Both integrations agree on non-recovery and voltage feasibility in these drifting cases. Exact final-angle digits should not be interpreted beyond the reported integration discrepancy.\n'
r+='''
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

'''
r+=f"- Maximum tighter-run final-state difference: {max(float(x['stateDifference']) for x in checks):.3g}.\n"
r+=f"- Maximum limiter-event time difference: {max(float(x['currentEventTimeDifference_s']) for x in checks):.3g} s; event counts match.\n"
r+=f"- Independent direct nodal power error: {max(float(x['directNodalPowerError']) for x in audit):.3g} pu.\n"
r+=f"- Physical gamma_dot check error: {max(float(x['physicalTimeDerivativeError']) for x in audit):.3g} rad/s.\n"
r+=f"- Maximum limiter-event residual: {max(float(x['limiterEventResidual']) for x in audit):.3g} pu; clearing state discontinuity is zero.\n"
r+=f"- Highest required command across all runs: {max(float(x['maxRequiredVoltage']) for x in rows):.9f} pu, below 1.3.\n"
r+='''- All current magnitudes obey 1.2 pu, and full metric inequalities pass where applicable. Baseline and certified metric hashes are unchanged (including the 59-file protected baseline manifest).
- `results.mat`: all 24 trajectories, solver solutions, complete settings, network, metrics and initial summary. `audit_results.mat` contains the augmented summary and independent checks.
- `case_*_trajectory.csv`: all requested state/electrical/limiter/metric/membership histories. Event, domain-transition and metric-transition CSVs accompany each run.
- `comparison_summary.csv`, `initial_conditions.csv`, `integration_verification.csv`, `model_consistency_audit.csv`, `settings.json` and `export_verification.json` provide machine-readable results.
- `run_study.m` implements the experiments, calls `plot_baseline_figures.m`, then `audit_study.m`. The same scripts were modified directly. To reproduce figures alone, call plot_baseline_figures with this result directory. `write_report.py` produces this report and checks exports/hashes. No new metric optimization or interval certification was performed.
- The complete prior 10-case directory, including its scripts and all outputs, was copied and hash-checked before the changes: `../radial_baseline_replication_20261001_archive_10cases_20261001_184741`. The current script refuses to overwrite results.mat; archive results first for another run.

## Figure links

'''
for name in ['Figure3_ConstantRatioLimiter','Figure4_ConstantRatioLimiter_Sag50','Figure4_ConstantRatioLimiter_Sag70','Figure5A_ConstantRatioLimiter','Figure5B_ConstantRatioLimiter_Sag50','Figure5B_ConstantRatioLimiter_Sag70','Supplement_Frequency_ConstantRatioLimiter','Supplement_Current_ConstantRatioLimiter_Sag50','Supplement_Current_ConstantRatioLimiter_Sag70']:
 r+=f"- {name}: [PNG]({(p/name).resolve().as_posix()}.png) | [vector PDF]({(p/name).resolve().as_posix()}.pdf)\n"
(p/'report.md').write_text(r,encoding='utf-8');print(p/'report.md')
