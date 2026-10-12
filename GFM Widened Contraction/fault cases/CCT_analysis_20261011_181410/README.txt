CCT analysis: same nominal GFM model, independent clearing-time searches

Run calculate_CCT_SCR5_XR10.m. A status-aware figure is exported after calculation. Each execution creates a
new timestamped results_* directory. Existing scripts/results are never run or overwritten.
The script is standalone with local functions; keep model_inputs.mat and baseline_results.mat
beside it. baseline_source_snapshot.txt records the inspected original source.

Inspected baseline: fault cases/Tclear_0p25/fixed_Pstar_recovery_SCR5_XR10_batch.m,
verify_affected_cases.m, model_inputs.mat, fixed_Pstar_sag_sweep_results.mat.
The latest saved cfg and source use tClear=0.10 s, despite the folder name and old
summary prose saying 0.25 s. This discrepancy is retained in the original files.
The CCT script first reproduces all 12 saved cases at the actual saved cfg.tClear.
Physics and segment integration functions were copied verbatim. No contraction
analysis, metric optimization, H/D/Imax scan, or controller change is performed.

Pstar=1; gamma0=20.633904875622 deg; nu0=1; omega=100*pi; mass=6.8; D=25;
Imax=1.2; converter voltage-command limit=1.3. All other network/controller
values come from the copied baseline model inputs, checked against saved results.

Recovery test (unchanged): full model-valid integration to tEnd=8 s from fault onset;
ALL saved postfault samples at absolute time t>=7.2 s must satisfy
abs(gamma-gammaEq)<0.1 degree and abs(nu-1)<1e-4; final abs(Pe-1)<1e-4.
Gamma is unwrapped, not compared modulo 360 degrees.
Samples include adaptive nodes, events and 6001 dense postfault samples;
fault samples include adaptive nodes, events and 501 dense samples.
Default RelTol=1e-9; AbsTol=[1e-11,1e-12]; MaxStep fault=0.002 s, postfault=0.005 s.
All limiter event equations and voltage/nu terminal event safeguards are unchanged.
Numerical integration stopping early without a physical event is UNRESOLVED.
Physical voltage-limit or nu<=0.5 violations are NOT_ASSESSED_VALIDITY_STOP.

Search: [0.001,1.0] s; 0.05 s coarse spacing, plus exactly 0.10,0.25,1.0 s.
Only one sampled recovered-to-not-recovered transition permits bisection.
Multiple/reversed transitions or unknown classifications are explicitly reported.
Bisection bracket width <=0.001 s. Estimate is bracket midpoint, not an exact value.
Both endpoints are repeated at RelTol=1e-11, AbsTol=[1e-13,1e-14], half maximum steps,
and at both 8 s and 16 s absolute horizons. At 16 s the recovery window starts at14.4 s.
Any changed endpoint classification prevents reporting a CCT estimate.
Coarse sampling cannot mathematically exclude unsampled recovery islands.

RECOVERED: above finite-horizon test passed.
NOT_RECOVERED: valid full run but that test failed; NOT proven instability.
NOT_ASSESSED_VALIDITY_STOP: outside the physical validity limits.
UNRESOLVED: numerical exception, premature numerical termination or ambiguity.
NUMERICAL_CCT: operational simulated recovery threshold with verified endpoints,
not a rigorous asymptotic stability CCT. NO_TRANSITION_UP_TO_MAX means all sampled
times recovered; it is not a computed CCT or proof for all longer/untested times.

Outputs: cct_results.csv (summary and bounds); cct_trials.csv (all search and
verification evaluations); baseline_regression.csv; cct_results.mat (model,
settings, classifications and full tightly checked endpoint trajectories).

CCT VISUALIZATION
The existing script now calls its local plot_cct_results function after saving results.
For plotting without recalculation, set cctPlotResultsDir to an existing results folder
and run the same calculate_CCT_SCR5_XR10.m script. Existing figure files cause an
explicit no-overwrite error; prior run figures and calculation results are preserved.
Outputs in the selected run: CCT_vs_Voltage_Sag.png (300 dpi), .fig, .pdf (vector).
All status labels and brackets come from the results table. The unresolved diamond
is located at the saved lower bracket endpoint, NOT a calculated CCT midpoint.
The numerical line is broken at unresolved or lower-bound cases. FCT=0.25 s is a
reference to the earlier fixed-clearing-time simulations.
Plot-only execution verified: 9 numerical, 2 lower-bound, 1 unresolved cases.
Previous result file hashes, original calculation block and physics functions unchanged.

Extended verification: set cctRefineResultsDir to the original results folder and run the same script. New results are saved in a fresh results_refined_* directory. Latest: results_refined_20261011_184900_867. Uret=0.45 CCT=0.15656640625 s, endpoints agree at 32/64 s under two tight tolerance sets; see refinement_notes.txt.
