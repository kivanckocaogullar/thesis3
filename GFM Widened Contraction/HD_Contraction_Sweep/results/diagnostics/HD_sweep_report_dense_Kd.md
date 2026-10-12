# H–K<sub>d</sub> contraction-bound sweep — dense K<sub>d</sub> grid (tag `dense_Kd`)

**Status: numerically grid-checked (301×301) — not interval-certified.** Not an invariant
set, region of attraction or recovery guarantee; separate branch metrics, no switching
certificate. No fault, FCT, CCT or EMT computations.

Grid: H ∈ {1, 3.4, 5.5, 8} s × K<sub>d</sub> ∈ {0, 1, …, 30} with K = 20 fixed, i.e.
D ∈ {20, 21, …, 50}; 124 points. Method and settings identical to `grid4x4_v2` (see
`HD_sweep_report_grid4x4_v2.md`); the 16 shared points were taken from the `grid4x4_v2`
checkpoints (identical configuration signature), the other 108 were computed
(17:13–20:01, 4 workers). 3202 candidate rectangles in total: 1871 feasible,
1146 `infeasible_lmi`, 123 `inconclusive_lmi`, 61 `inconclusive_refinement`,
1 `inconclusive_dense`, 0 solver failures (12.5 CPU-h).

Main figure for the dense grid: `figures/region_comparisons/profiles_vs_D_dense_Kd`
(bounds versus D, one line per H). The 3D surfaces, the symmetric-band surfaces and the
feasibility map are also exported.

## 1. Limiter-inactive branch (λ = 1 s⁻¹)

* **Feasibility coincides with the analytic necessary condition D > 4H at all 124 points**
  (0 mismatches). No feasible rectangle for H = 5.5, D ≤ 22 and H = 8, D ≤ 32 (16 points);
  the first feasible values are D = 23 (D/H = 4.18) and D = 33 (D/H = 4.13). Within the
  1-unit D resolution the necessary condition is therefore also numerically sufficient
  for this metric family. Infeasibility for D ≤ 4H is a mathematical result (trace
  argument), not only a numerical one.
* Where feasible, γ<sub>min</sub> = −24.51° (inactive mode-validity floor) and the symmetric
  frequency band reaches the search limit 1 ± 0.05. Eight isolated values below 0.05
  (down to 0.0345, e.g. H = 8, D = 44 and 46) coincide with `inconclusive_refinement` at the
  limit probe and are numerical, not a trend.
* Conclusion: at λ<sub>inactive</sub> = 1 the inactive branch is an on/off property governed by
  D/H; inside the feasible set its rectangle is limited by the current limiter (angle) and
  by the search domain (frequency), not by H or D.

## 2. Limiter-active branch (λ = 0.1 s⁻¹)

Feasible at all 124 points in both studies.

**γ<sub>max</sub> (Study A, ν ∈ [0.98, 1.02]):**

| H (s) | D = 20 | D = 25 | D = 35 | D = 50 | trend in D |
|---|---|---|---|---|---|
| 1 | 48.97 | 48.88 | 48.37 | 47.52 | decreasing, ±0.5° scatter (7 sign changes) |
| 3.4 | 44.45 | 45.05 | 45.64 | 46.24 | non-decreasing |
| 5.5 | 41.80 | 42.67 | 43.94 | 44.79 | non-decreasing |
| 8 | 39.15 | 40.14 | 41.74 | 43.34 | strictly increasing |

At every D, γ<sub>max</sub> decreases with H. The H = 1 curve is the only non-monotone one;
its scatter exceeds the 0.1° bisection tolerance and several of its brackets contain
inconclusive probes, so the fine structure there is not resolved.

**Frequency (Study B, γ ∈ [24.42°, 43°]) — largest symmetric band 1 ± δ\*:**

* H = 1 and 3.4: δ\* at the search limit 0.05 (except isolated values 0.0402 at
  H = 3.4, D = 27 and 0.0495 at D = 29, 38).
* H = 5.5: δ\* rises smoothly 0.0126 (D = 20) → 0.0345 (D = 35), then reaches the limit
  from D = 36 (two isolated lower values, 0.0448 at D = 38 and 0.0439 at D = 41).
* H = 8: δ\* rises monotonically at every step, 0.0060 (D = 20) → 0.0256 (D = 50).

**One-sided extents ν<sub>min</sub>, ν<sub>max</sub>:** these are conditional on the opposite edge held
at 1 ∓ δ\*(D), and the active branch shows a strong trade-off between the two edges.
Example (H = 8): the upper search tests [1 − δ\*, 1.05]; it is feasible at D = 41 and 43 but
not at D = 40, 42, 44, which makes ν<sub>max</sub> jump between ≈ 1.015 and 1.05. A separate
robustness test (`robustness_active_numax_H8.csv`) fitted the *same* rectangle
ν ∈ [0.9845, 1.05] with initial grids 7, 9, 11 and 13: infeasible for all grids at D = 40,
feasible for all grids at D = 41–44. The fit is therefore insensitive to the
initialization; the jumps come from the lower edge 1 − δ\*(D) moving down as D grows
(0.9836 at D = 42 is already too low for ν<sub>max</sub> = 1.05). Use δ\* as the jointly feasible
frequency measure; the one-sided values show the trade-off, not independent bounds.

## 3. Numerical caveats

* Near the boundary the deciding `feasp` t<sub>min</sub> values are of order 1e-9 – 1e-6 while LMI
  entries are of order 10³ (ω<sub>b</sub> = 314 in J). "Feasible" results are independently
  confirmed on the 301×301 grid; "infeasible_lmi" results close to the boundary are not
  infeasibility certificates for the continuous problem. Reported bounds are the largest
  rectangles found feasible.
* Grey crosses in the profile figure mark bounds whose final bracket contained an
  inconclusive probe (bound possibly conservative by at most the bracket width).
* Search-domain limits (inactive mode floor, active saddle γ<sub>u</sub>, ν ∈ [0.95, 1.05]) cap
  many reported values (open markers).

## 4. Suggested next steps

1. Interval-verify representative metrics (JSON exported per point and branch), e.g.
   baseline (3.4, 25); (8, 33) and (5.5, 23) just above the inactive threshold; (8, 20)
   narrowest active band; (1, 20) largest active angle.
2. If smoother surfaces along H are wanted, add H values (e.g. 2, 4.5, 6.5) with the same
   script; the threshold D = 4H predicts the inactive feasibility boundary in advance.
3. The ν search domain (±0.05) saturates most points; widen it only if frequencies beyond
   ±2.5 Hz are physically relevant.
