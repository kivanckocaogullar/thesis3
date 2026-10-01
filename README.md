# Thesis backup: current-limited grid-forming inverter analysis

This repository is a preservation copy of the recent MATLAB and numerical-analysis work on a grid-forming inverter (GFM) with current limitation. The source snapshots and computed results are retained under `outputs/`; the selected earlier baseline figures are retained under `work/`. Original file and directory names are preserved wherever practical. Numerical data, mathematical certificates, model equations, simulation settings, and limiter behavior have not been changed for this backup.

## Contents

- `outputs/baseline_SCR5_XR10/` — limiter-free SCR 5 / X/R 10 baseline, original MATLAB analysis, network data, disturbance results, and editable/exported figures.
- `outputs/current_limiter_SCR5_XR10/` and root-level `gfm_current_limiter_*.m` — original constant-radial-current-limiter benchmark, branchwise analysis, boundary and sensitivity data, and disturbance simulations.
- `outputs/current_limited_*` — reduced-order current-limited VSG model and contraction-region work.
- `outputs/radial_inactive_metric_20261001_170815/` — certified state-dependent inactive metric `M_inactive` on R1, with full-precision coefficients and directed interval checks.
- `outputs/radial_quadratic_robust_20261001_161903/` — certified state-dependent active metric `M_active` on R2, with full-precision coefficients and directed interval checks.
- `outputs/radial_switch_connection_20261001/` — one-sided switching-boundary analysis, metric-domain extensions, transition-gain verification, and the recovery-set analysis using the preserved metrics.
- `outputs/active_metric_R1_crossvalidation_20261001/` — latest cross-validation of the unchanged `M_active` on inactive R1; this is a numerical cross-validation and records its contraction counterexample.
- `outputs/radial_baseline_replication_20261001/` — latest 24-case study, with initial angles A–L at 14, 20.633904875622, 24.373430419712, 27, 36, 16, 17, 23, 24.2, 24.55, 32, and 34.2 degrees. Each is simulated at 50% and 70% source-voltage sag. The folder includes saved trajectories, current events, metric transitions, model settings, figures, comparison tables, and integration audits.
- Other `outputs/` folders retain the preceding disturbance studies, diagnostics, controller/model inspections, switching studies, sensitivity studies, and archived 10-case results.
- `work/previous_gfm_outputs/` and `work/superseded_sg_validation/` retain selected earlier numerical figures and result snapshots.
- The remaining `work/` files preserve intermediate source snapshots, pre-edit hash manifests, and preliminary result data from the same analysis sessions.
- `model/GFMVSG1.slx` and `model/tdmsquared-master/tdmsquared-master/Library/` preserve the source GFM Simulink model and its referenced controller library files.

## Models and analysis

The work starts with the limiter-free SCR5/XR10 electrical model and baseline contraction analysis. It then constructs a two-state current-limited reduced-order VSG benchmark with the **Constant Ratio / radial limiter**. The limiter applies radial projection to the converter current at the stated 1.2 pu limit; the model also checks the 1.3 pu voltage-command limit. The network reconstruction and switching branches are documented in the accompanying reports. This is a reduced model with algebraic, instantaneous Q–V closure and ideal current-reference tracking; it omits the dynamic PI, filter, and PLL states of the detailed controller.

The certified state-dependent metrics are branch-specific:

- `M_inactive` (`R1`) is interval-validated at contraction rate 1 s⁻¹ on its reported inactive rectangle. Its polynomial coefficients, center, scales, interval checks, and report are in `outputs/radial_inactive_metric_20261001_170815/`.
- `M_active` (`R2`) is interval-validated at contraction rate 0.1 s⁻¹ on its reported active rectangle. Its full-precision coefficients, center, scales, interval checks, and report are in `outputs/radial_quadratic_robust_20261001_161903/`.
- `outputs/radial_switch_connection_20261001/` preserves the one-sided switch calculations, validated extensions of the original metrics, conservative transition gains, and a separate energy/Lyapunov recovery subset. It does not claim a gain-free hybrid-contraction or full EMT certificate.
- `outputs/active_metric_R1_crossvalidation_20261001/` tests the unchanged active metric on R1 and records that positivity passed the sampled grid while the full contraction inequality failed there. It is not an R1 certificate.

The detailed reports are the authoritative source for certificate scope, interval methods, coefficients, assumptions, numerical-versus-uniform bounds, and limitations. MATLAB verification scripts and Python interval-validation helpers are saved alongside the data.

## Reproducing the latest study

The latest run script is `outputs/radial_baseline_replication_20261001/run_study.m`. It uses saved model and certificate inputs from sibling folders in `outputs/`, and refuses to overwrite an existing `results.mat`. To reproduce it, make a **disposable copy of the complete `outputs/` tree**, remove `radial_baseline_replication_20261001/results.mat` from that copy, and run the copied script in MATLAB, for example:

```matlab
addpath('path/to/disposable-copy/outputs/radial_baseline_replication_20261001');
run_study
```

The script calls the companion plotting and audit scripts after the integrations. Keep the checked-in results intact and use the disposable copy for any rerun. Settings for the 24 cases are recorded in `settings.json`; the script also contains assertions for the saved physical model, metric coefficients, switching point, and integration outcomes.

To inspect or repeat the metric checks, use the `check_*.m`, `validate_*.py`, and `check_*.py` scripts beside their reports. The Python interval routines use `mpmath`; the supplied interval-validation dependency bundle is intentionally not vendored. Install `mpmath` in the Python environment used for those scripts. Some detailed controller-audit folders contain Simulink `.slx` files and require MATLAB/Simulink to open.

The external publication PDF used by the paper-reading helper is not included. Set `GFM_REFERENCE_PDF` to a local copy before running `read_paper.py` or `plain.py`. Model inspection scripts resolve the bundled `model/` assets relative to the repository.

## Data and file formats

- `.m` — MATLAB analysis, simulation, verification, and plotting source.
- `.py` — interval arithmetic and report/validation helpers.
- `.mat` — saved model states and numerical results.
- `.csv` / `.json` — case data, settings, coefficients, checks, and machine-readable results.
- `.fig`, `.png`, `.pdf` — editable MATLAB figures and view/export copies.
- `.slx` — saved Simulink models used for controller/model inspection.

Large `.mat`, `.csv`, and `.fig` files are stored with Git LFS. Install Git LFS before cloning if you need to download their contents.
