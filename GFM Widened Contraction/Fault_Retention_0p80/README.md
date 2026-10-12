# Five-FCT study: voltage retention 0.80 pu

Voltage drop: 20%. Selected FCTs: [.10 .25 .50 .75 1.00] s. Fault starts at t=0, returns to the unchanged nominal voltage at the selected FCT. The observation ends at t=8 s from inception, exactly as in the reference; it is not eight additional seconds after clearing.

## Running

In MATLAB, add this folder's `scripts` directory to the path and call `run_fault_retention_0p80`. Saved numerical results are reused if present. To regenerate figures only: `run_fault_retention_0p80('plots')`. For an explicit fresh five-run reproduction without overwriting this study: `run_fault_retention_0p80('all', fullfile(pwd,'reproduction_0p80'))`, using a new output directory. No CCT search, optimization, interval verification or additional tolerance runs are performed.

Shared code and immutable input snapshots are in the sibling `Fault_Retention_0p60/scripts`; keep the two folders together. `fct_retention_engine.m` holds the original model, ode45/event integration and unchanged recovery functions (copied verbatim from `analyze_FCT_60pct.m`). `plot_retention_FCT.m` and `plot_retention_containment.m` adapt the existing widened-region plotters. `compare_retention_results` loads the saved original 0.40 study and produces the combined table. `validate_retention_results` checks saved trajectories without rerunning simulations. Entry points and shared functions do not modify reference files.

## Configuration and references

SCR=5; X/R=10; Pstar=1; mass=6.8; D=25; Imax=1.2; Vmax=1.3; omega=100*pi. Initial gamma=20.633904875622 deg (the full precision equilibrium is saved in MAT), nu=1. The reference model's nominal U is 0.980392156862745; retention multiplies that same nominal U, so fault U=0.784313725490196. It is restored exactly to the reference nominal U on clearing. Do not replace U with 1: that would change the reference model.

ode45: RelTol=1e-9, AbsTol=[1e-11 1e-12], fault MaxStep=0.002 s, postfault MaxStep=0.005 s. All original event definitions and saved-sample settings are unchanged. `cfg.selectedFCT` records the sweep; the unused baseline scalar cfg.tClear is removed to avoid a misleading single-duration label.

Recovery: on all saved postfault samples t>=7.2 s, |gamma-gammaEq|<0.1 degrees and |nu-1|<1e-4, and final |Pe-1|<1e-4. Original classification strings are preserved. NOT_RECOVERED means no recovery by 8 s, not proven instability. Validity stops and unresolved errors retain their original distinct labels; all ten new runs are valid.

Reference sources (exact paths and SHA256 in `Fault_Retention_0p60/scripts/inputs/source_hashes.json`):

- Original `FCT_sensitivity_60pct/analyze_FCT_60pct.m` and saved five-FCT MAT from the 2026-10-11 study; the 0.40 case is never simulated here.
- `GFM Widened Contraction/CCT_analysis_20261011_181410/baseline_results.mat`, model inputs and refined CCT CSV.
- `GFM Widened Contraction/contraction region/metrics/inactive` and `active`. Loaded coefficients and bounds match `GFM Contraction Analsys/limiter region/widened_nu_20261012_130655` exactly.
- Original `FCT_60pct_widened/plot_saved_widened_FCT.m` and `plot_widened_containment.m`.

Inactive rectangle: gamma [10,24.3234] degrees, nu [0.98,1.02], lambda=1. Active: gamma [24.4234,43], nu [0.98,1.02], lambda=0.1. M>0 and M<10I; saved metrics and strict interval certificates reused. The 0.1-degree uncertified gap is excluded. Physical switching boundary is 24.373430419712 degrees (vertical dashed line in region figures).

## Files

- `results/fct_analysis.mat`: full original solver solutions, dense saved samples, p/n/pf/nf/cfg, FCTs, metrics, certificates, detailed and summary tables.
- `case_details.csv`: all requested extrema, membership, recovery and sampled exit quantities. Gamma extrema are unwrapped and postfault; nu/current/voltage extrema cover fault-on plus postfault samples. Sampled extrema are not continuous-time certified bounds.
- `five_FCT_summary.csv`, `assessment_details.csv`, `limiter_active_intervals.csv`, per-FCT `trajectory_*.csv`, `numerical_configuration.json`, `cct_reference.csv`, `validation.csv`.
- `combined_040_060_080.csv/.mat`: original 0.40 numbers loaded from MAT unchanged, with all three studies assessed using the same widened rectangles.
- `figures/`: 01_gamma_vs_time, 02_nu_vs_time, 03_phase_plane, 04_contraction_zoom, trajectory_containment, 05_current_limiter, each 300-dpi PNG, vector PDF and editable FIG.

## Interpretation and numerical limits

Membership is mode-aware. `LeavesClearingModeRectangle` refers to the rectangle applicable at clearing; `LeavesCertifiedUnion` includes the gap. For initially outside cases these flags are false and exit time is NaN (not applicable), while `FirstOutsideUnion_s_after_clear`=0. The legacy containment table names its first-outside quantity `First_Exit_s_after_clearing`; it too uses 0 for initial outside. NaN means no sampled outside point. Exit times are the first outside saved sample, not a continuous-time root. Limiter intervals use the existing ode45 event boundaries and midpoint mode test. Time-inside fractions use the reference trapezoidal estimate.

The strict active certificate is VERIFIED, but it does not certify a uniform 1e-5 epsilon margin over every box. No new stronger certificate is claimed. Metric eigenvalues/condition numbers in case_details are clearing-state diagnostics; evaluations outside a rectangle are extrapolations, not certificates. Region membership alone guarantees neither recovery nor invariance, and leaving a region does not prove instability. Separate metrics do not establish switching certification. The 0.80 CCT record is a tested lower-bound result through 1 s, not an actual CCT of 1 s.

## Results

|FCT (s)|gamma clear (deg)|nu clear|Recovery by 8 s|Limiter active at clearing|Rectangle|
|---|---|---|---|---|---|
|0.1|22.779333|1.002062486|RECOVERED|0|INACTIVE|
|0.25|28.518478|1.001849631|RECOVERED|1|ACTIVE|
|0.5|35.064564|1.001179442|RECOVERED|1|ACTIVE|
|0.75|39.954482|1.001060189|RECOVERED|1|ACTIVE|
|1|45.076371|1.001261382|RECOVERED|1|OUTSIDE|
