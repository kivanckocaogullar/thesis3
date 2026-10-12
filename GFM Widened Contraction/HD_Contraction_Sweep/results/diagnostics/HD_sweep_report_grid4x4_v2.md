# H–K<sub>d</sub> contraction-bound sweep — preliminary 4×4 results (tag `grid4x4_v2`)

**Status: numerically grid-checked (301×301) — not interval-certified.** A feasible
rectangle is not an invariant set, a region of attraction or a recovery guarantee; the
two limiter branches use separate metrics and no switching certificate exists. No fault,
FCT, CCT or EMT computations were made.

## 1. Model facts verified from the code

* Reduced model (physical time): γ̇ = ω<sub>b</sub>(ν−1), ν̇ = ((P\*−P<sub>e</sub>)/ν − D(ν−1))/m.
* VSM block (`Subsys_GFormConv_OuterCtrl_VSM.slx`):
  2H ω̇ = (P\*−P)/ω + K(ω<sub>ref</sub>−ω) − K<sub>d</sub>(ω−ω<sub>PLL</sub>). With ω<sub>ref</sub> = 1 and
  ω<sub>PLL</sub> ≈ 1 (the existing approximation): **m = 2H, D = K + K<sub>d</sub>**. Stored values
  H = 3.4, K = 20, K<sub>d</sub> = 5, D = 25, m = 6.8 (asserted in `hd_set_HD.m`).
* H enters only through m; D only through the damping term. As requested, **K = 20 is
  fixed and K<sub>d</sub> is varied** (K<sub>d</sub> ∈ {0, 5, 15, 30} → D ∈ {20, 25, 35, 50}). In this
  reduced model only the sum D matters; the split matters for later EMT studies.
* The electrical closure (`hd_coupled`) uses neither m nor D, so the equilibrium
  (γ<sub>eq</sub> = 20.634°), P<sub>e</sub>(γ), the limiter switch (γ<sub>sw</sub> = 24.373°), the
  voltage-command limit, the inactive mode-validity floor (−24.51°, where |q| reaches
  I<sub>max</sub> again) and the active-branch saddle (γ<sub>u</sub> = 64.88°) are independent of H, K<sub>d</sub>.

## 2. Baseline consistency

* Saved interval-verified metrics (`widened_nu_20261012_130655`) re-evaluated with the
  copied functions: inactive dense min eig M = 9.88109e-4, max eig Q = −9.22243e-4;
  active 1.27997e-5, −1.08629e-5 — identical to the saved values (`baseline_check.csv`).
* Re-fitting the baseline rectangles with the sweep code reproduces the same numbers.
* The baseline point of the sweep (H = 3.4, K<sub>d</sub> = 5) contains both baseline rectangles.

## 3. Method and documented deviations

Unchanged: λ<sub>inactive</sub> = 1, λ<sub>active</sub> = 0.1 s⁻¹, ε = 1e-5, M ≺ 10I, quadratic metric,
Ṁ from the vector field, 9×9 → `feasp` → 101×101 refinement (≤ 8 rounds) → 301×301
acceptance. Metrics re-optimized from scratch per candidate (no warm start).

A first run (tag `grid4x4`, kept in `results/diagnostics`, figures moved to
`figures/_superseded_run1_grid4x4`) gave many inconclusive outcomes, including at the
baseline. Causes and corrections in `grid4x4_v2`:

1. **Dense-grid feedback (documented extension).** Metrics passing the 101×101 test
   often missed the ε margin on the 301×301 grid by ~1e-6. The worst 301×301 violators
   are now added to the LMI grid for at most 4 extra rounds. The acceptance criterion is
   unchanged. Of 244 feasible probes, 134 needed no extra round, 70 one, 34 two, 6 three.
2. **feasp rule.** t<sub>min</sub> < 0 is accepted as in the original code (run 1 wrongly marked
   |t<sub>min</sub>| < 1e-7 inconclusive); t<sub>min</sub> ∈ [0, 1e-7] inconclusive; > 1e-7 infeasible.
3. **Baseline-anchored searches.** Run 1 started from very small rectangles (e.g.
   γ ∈ [24.42°, 26°]); the coordinate normalization then makes Ṁ large and the finite-grid
   fit unreliable. Searches now start at the baseline edges (10°, 43°, δ = 0.02) and move
   outward, or inward if the baseline rectangle fails.

With these changes the initial-grid sensitivity test (7×7, 9×9, 11×11 at the baseline)
gives identical bounds (−24.51°, 45.047°); in run 1 it did not.

## 4. Results

### Study A — angle bounds (ν ∈ [0.98, 1.02] fixed)

Fixed edges (limiter switch, independent of H, D): γ<sub>max,inactive</sub> = 24.3234°,
γ<sub>min,active</sub> = 24.4234°.

γ<sub>min,inactive</sub>: **−24.51° at every feasible point** — the search reached the inactive
mode-validity floor, i.e. the bound is set by the current limiter, not by contraction.
No feasible rectangle at (H, D) = (5.5, 20), (8, 20), (8, 25).

γ<sub>max,active</sub> (deg), limited by contraction feasibility (tolerance 0.1°):

| H \ D | 20 | 25 | 35 | 50 |
|---|---|---|---|---|
| 1 | 48.97 | 48.88 | 48.37 | 47.52 |
| 3.4 | 44.45 | **45.05** | 45.64 | 46.24 |
| 5.5 | 41.80 | 42.67 | 43.94 | 44.79 |
| 8 | 39.15 | 40.14 | 41.74 | 43.34 |

### Study B — frequency bounds (γ fixed at the baseline intervals)

Search domain ν ∈ [0.95, 1.05], tolerance 0.0005. "≥1.05" / "≤0.95": domain limit reached.

Inactive (γ ∈ [10°, 24.32°]): at every feasible point ν<sub>min</sub> ≤ 0.95 and ν<sub>max</sub> ≥ 1.05
(one-sided). Symmetric band 1 ± 0.05 jointly feasible except (5.5, 35): 1 ± 0.0458 and
(8, 50): 1 ± 0.0453, where the combined [0.95, 1.05] rectangle was inconclusive
(refinement budget, not infeasibility). No feasible rectangle at (5.5, 20), (8, 20), (8, 25).

Active (γ ∈ [24.42°, 43°]), ν<sub>max</sub>:

| H \ D | 20 | 25 | 35 | 50 |
|---|---|---|---|---|
| 1 | ≥1.05 | ≥1.05 | ≥1.05 | ≥1.05 |
| 3.4 | ≥1.05 | **≥1.05** | ≥1.05 | ≥1.05 |
| 5.5 | 1.0126 | 1.0170 | 1.0365\* | ≥1.05 |
| 8 | 1.0060 | 1.0075 | 1.0117 | 1.0256 |

ν<sub>min,active</sub> ≤ 0.95 everywhere except (8, 20): 0.9898. \*(5.5, 35): one-sided extents
not jointly feasible; jointly feasible symmetric band 1 ± 0.0345.

## 5. Trends (only what the 4×4 data support)

* **Inactive branch:** feasibility is exactly separated by the analytic necessary
  condition **D > 4λ<sub>inactive</sub>H = 4H** on this grid (infeasible: D/H = 3.6, 2.5, 3.1;
  feasible down to D/H = 4.375 at (8, 35) and 4.55 at (5.5, 25)). Where feasible, both
  the angle and the frequency bounds are at their limits (mode floor / search domain), so
  H and D act as an on/off switch for this branch at λ = 1. The fixed rate λ<sub>inactive</sub> = 1
  is what excludes the high-H, low-D corner.
* **Active branch, angle:** γ<sub>max</sub> decreases monotonically with H at every D
  (≈ 9.8° from H = 1 to 8 at D = 20; 4.2° at D = 50). With D it increases for
  H ≥ 3.4 but decreases slightly for H = 1 — non-monotone across H.
* **Active branch, frequency:** the upper bound ν<sub>max</sub> is the binding edge; it shrinks
  with H and grows with D wherever it is not domain-limited. ν<sub>min</sub> binds only at (8, 20).
* The selected rates prevented feasibility only for the inactive branch (three points).
  No active-branch search failed.

## 6. Solver diagnostics

430 candidate rectangles: 244 feasible, 160 `infeasible_lmi`, 14 `inconclusive_lmi`,
12 `inconclusive_refinement`, 0 solver failures. Inconclusive probes occurred inside
bisection brackets in 12 of the 64 study/branch searches (flagged in `summary_table_grid4x4_v2.csv`,
column `warnings`); treating them as "not shown feasible" makes the reported bounds
conservative by at most the bracket width. The non-monotone inactive δ\* values at
(5.5, 35) and (8, 50) coincide with `inconclusive_refinement`, so they are numerical
rather than evidence of infeasibility. Total probe time 1.68 h (41 min wall, 4 workers).

## 7. Limitations

Numerical 301×301 checks only; ε-margin acceptance; quadratic metric; bisection assumes
monotone feasibility (true for the exact problem, not guaranteed numerically); search
domains (mode floor, γ<sub>u</sub>, ν ∈ [0.95, 1.05]) cap several reported bounds; the reduced
model lumps K and K<sub>d</sub>; separate branch metrics without a switching certificate.

## 8. Recommendations

* **Denser grid:** worthwhile for the active branch (γ<sub>max</sub>, ν<sub>max</sub> vary smoothly with
  H and D) and near the inactive threshold D ≈ 4H. A K<sub>d</sub> step of 1 (D = 20…50) with
  the four H values needs ≈ 3–4 h (≈ 6.3 CPU-min per point, 4 workers). A denser H grid
  (e.g. 1:1:8) would make the surfaces more informative along H.
* **Interval verification candidates** (JSON already exported):
  baseline (3.4, 25); (1, 20) largest active angle; (8, 20) most restrictive active
  frequency band; (8, 35) and (5.5, 25) just above the inactive threshold.
* **Later fault-recovery tests:** the same points plus (8, 25), where no inactive
  contraction rectangle exists at λ = 1, as a contrast case.
