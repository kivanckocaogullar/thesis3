# H–D contraction-bound sweep (GFM/VSG, current-limited)

How the numerically feasible contraction rectangles of the current-limited GFM/VSG
model change with virtual inertia *H* and damping *D*. The contraction method is the
baseline method, unchanged; only *H* and *K*<sub>d</sub> are varied.

> **All sweep results are numerically grid-checked (301×301) — not interval-certified.**
> A feasible rectangle is not an invariant set, a region of attraction or a recovery
> guarantee. The two limiter branches use separate metrics; no common switching
> certificate exists.

## How to run

```matlab
addpath('scripts')
run_HD_contraction_sweep                                  % 4x4 grid, tag 'grid4x4_v2'
run_HD_contraction_sweep('Kd',0:30,'tag','dense_Kd')      % denser Kd grid (D = 20..50)
run_HD_contraction_sweep('H',1:8,'Kd',0:5:30,'tag','g8x7')% any user grid
```

Completed (*H*, *K*<sub>d</sub>) points are checkpointed in
`results/diagnostics/checkpoints_<tag>/` and skipped on restart (a checkpoint is reused
only if its configuration signature matches). Every output file carries the tag, so runs
never overwrite each other. Plots only: `plot_HD_gamma_surfaces(root,tag)` etc.

## Model facts verified from the code

| Item | Finding |
|---|---|
| Reduced model | γ̇ = ω<sub>b</sub>(ν−1), ν̇ = ((P\*−P<sub>e</sub>)/ν − D(ν−1))/m |
| *H* | enters only through m = 2H (VSM block gain `1/(2*H)`) |
| *D* | D = K + K<sub>d</sub> (VSM: K(ω<sub>ref</sub>−ω) − K<sub>d</sub>(ω−ω<sub>PLL</sub>), ω<sub>ref</sub> = 1, ω<sub>PLL</sub> ≈ 1). Baseline K = 20, K<sub>d</sub> = 5, D = 25 |
| Sweep convention | **K = 20 fixed, K<sub>d</sub> varied**, so D = 20 + K<sub>d</sub> ≥ 20. In this reduced model only the sum D matters |
| Independent of H, D | equilibrium (P<sub>e</sub>(γ) = 1, γ<sub>eq</sub> = 20.634°), P<sub>e</sub>(γ), limiter switch (γ<sub>sw</sub> = 24.373°), voltage-command limit, mode validity, active-branch saddle γ<sub>u</sub> = 64.88° |
| Baseline | H = 3.4 s, K = 20, K<sub>d</sub> = 5, D = 25, m = 6.8, SCR 5, X/R 10, P\* = 1, I<sub>max</sub> = 1.2, V<sub>max</sub> = 1.3 |

## Unchanged method

λ<sub>inactive</sub> = 1 s⁻¹, λ<sub>active</sub> = 0.1 s⁻¹ (fixed); ε = 1e-5; M ≻ εI, M ≺ 10I,
JᵀM + MJ + Ṁ + 2λM ≺ −εI; quadratic metric in normalized coordinates; Ṁ from the
nonlinear vector field; 9×9 initial LMI grid, `feasp`, 101×101 adaptive refinement
(≤ 8 rounds, ≤ 30 nodes per round), 301×301 dense check. The core functions in
`scripts/hd_fit_metric.m`, `hd_verify_metric.m`, `hd_metric_poly.m`, `hd_coupled.m` are
verbatim copies of the local functions of
`GFM Contraction Analsys/limiter region/widen_nu_contraction.m` (names changed only).
Metric coefficients are re-optimized from scratch for every candidate rectangle; the
baseline metrics are **not** used as warm starts.

Documented procedural extensions (acceptance criterion unchanged):
* `feasp` t<sub>min</sub> < 0 accepted as in the original code; t<sub>min</sub> ∈ [0, 1e-7] inconclusive.
* If the 301×301 check fails, its worst violating nodes are added to the LMI grid for at
  most 4 extra rounds (`cfg.maxDenseRounds`).
* Edge searches start from the baseline edges and move outward (or inward if the baseline
  rectangle fails), avoiding very small anchor rectangles.

The first run without these extensions (tag `grid4x4`) is kept in `results/diagnostics`
for traceability; its figures are in `figures/_superseded_run1_grid4x4`. See
`results/diagnostics/HD_sweep_report_grid4x4_v2.md`, section 3.

Added checks (γ-only, independent of H, D): limiter-mode consistency
(|q| < I<sub>max</sub> inactive, > I<sub>max</sub> active) and |u| < V<sub>max</sub> on a 2001-point γ grid.

## Boundary definitions

**Study A (angle bounds), ν ∈ [0.98, 1.02] fixed.**
Inactive: upper edge fixed at 24.3234° (switch − 0.05°), lower edge searched from the
baseline edge 10° down to the inactive mode-validity floor −24.51° (or up to 20° if the
baseline rectangle fails). Active: lower edge fixed at 24.4234° (switch + 0.05°), upper
edge searched from the baseline edge 43° up to 64.83° (or down to 26°)
(0.05° below the saddle γ<sub>u</sub>; no contraction rectangle can contain a non-attracting
equilibrium). Bisection tolerance 0.1°. The fixed edges are set by the limiter switch,
not by optimization.

**Study B (frequency bounds), γ fixed at the baseline intervals**
(inactive [10°, 24.3234°], active [24.4234°, 43°]). Search domain ν ∈ [0.95, 1.05]
(±2.5 Hz at 50 Hz; covers the 47.5–51.5 Hz grid-code withstand band). Rule:
(1) largest symmetric band [1−δ, 1+δ], searched from the baseline δ = 0.02 outward to
0.05 (or inward to 0.001); (2) ν<sub>min</sub> searched with
ν<sub>max</sub> = 1+δ\* held; (3) ν<sub>max</sub> searched with ν<sub>min</sub> = 1−δ\* held; (4) the
combined rectangle [ν<sub>min</sub>, ν<sub>max</sub>] is tested once and its status reported.
Tolerance 0.0005 pu. Edges that reach the domain limit are flagged.

Bisection is justified because feasibility is monotone under domain inclusion for the
exact problem; the numerical procedure may violate this, so every probe is stored and
inconclusive probes inside a bracket are flagged.

**Outcome classes per candidate:** `feasible`; `infeasible_lmi` (feasp t<sub>min</sub> > 1e-7 on
the LMI nodes — no metric of this structure satisfies the constraints even on the finite
grid); `mode_or_voltage_invalid`; `inconclusive_lmi` (|t<sub>min</sub>| ≤ 1e-7);
`inconclusive_refinement`; `inconclusive_dense`; `solver_failed`.

**Analytic necessary condition.** At the equilibrium Ṁ = 0, so tr Sym(F) = tr J = −D/(2H).
Both eigenvalues ≤ −λ requires **D ≥ 4λH** for every rectangle containing the
equilibrium (all inactive rectangles).

## Runs

| Tag | Grid | Report |
|---|---|---|
| `grid4x4_v2` | H ∈ {1, 3.4, 5.5, 8}, K<sub>d</sub> ∈ {0, 5, 15, 30} | `results/diagnostics/HD_sweep_report_grid4x4_v2.md` |
| `dense_Kd` | H ∈ {1, 3.4, 5.5, 8}, K<sub>d</sub> = 0:30 (D = 20…50), 124 points | `results/diagnostics/HD_sweep_report_dense_Kd.md` |
| `grid4x4` | first run, superseded (see report v2, section 3) | — |

For dense D grids the clearest figure is `figures/region_comparisons/profiles_vs_D_<tag>`
(`plot_HD_line_profiles.m`). `plot_HD_nu_surfaces.m` exports both the one-sided extents
(`nu_bounds_<tag>`) and the largest symmetric band (`nu_bounds_symmetric_<tag>`); in the
active branch the lower and upper frequency edges trade off strongly, so the symmetric band
is the jointly feasible quantity. `hd_robustness_check.m` reproduces the initial-grid
robustness test of the dense report.

## Folder layout

```
scripts/        entry point, searches, verbatim core functions, plotting
results/angle_bounds/       angle_bounds_<tag>.csv/.mat, metrics_json_<tag>/ (interval-checker format)
results/frequency_bounds/   frequency_bounds_<tag>.csv/.mat, metrics_json_<tag>/
results/diagnostics/        summary_table_<tag>.csv, all_probes_<tag>.csv, sweep_results_<tag>.mat,
                            baseline_check.csv, domain_limits.json, init_sensitivity.csv,
                            run_log_<tag>.txt, checkpoints_<tag>/, HD_sweep_report_<tag>.md
figures/...                 PNG (300 dpi), PDF (vector), FIG; plotted data as CSV
```

Interval verification of a selected point:
`python "GFM Widened Contraction/contraction region/widened_interval_checker.py" results/angle_bounds/metrics_json_<tag>/<point>/<branch> <branch>`
