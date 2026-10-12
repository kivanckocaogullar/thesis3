# Five FCT cases with the widened contraction regions

Run `analyze_saved_widened_FCT.m` in MATLAB to reproduce the analysis and all five figures in `results`. It reads the included saved FCT trajectories and the sibling CCT/contraction-region folders. It does not execute an ODE solver, a CCT search or metric optimization.

## Reuse validation
The CCT script is byte-identical to the script used by the original five-case study. The cached trajectories' p, n, fault voltage, rebuilt fault network, initial state, recovery horizon and solver settings match the copied CCT inputs and widened-region model. The widened metric MAT and JSON coefficients and bounds agree, network parameters agree, and both supplied interval certificates report VERIFIED with no unresolved/counterexample boxes. Existing project files are unchanged (SHA-256 audit).

60% sag means Uret=0.40, not Uret=0.60. FCT values are 0.05, 0.10, 0.14, 0.15 and 0.25 s; the simulation horizon remains 8 s from fault inception. Gamma is unwrapped; nu=omega/(100*pi). The unchanged recovery test uses all saved postfault states at t>=7.2 s: angle error <0.1 degree, |nu-1|<1e-4, and final |Pe-1|<1e-4.

## Regions and certificates
Inactive gamma [10,24.3234] degrees, active gamma [24.4234,43] degrees; both nu [0.98,1.02]. Rates remain 1/s and 0.1/s, respectively, and 0<M<10I is retained. The original switching gap is excluded.
The supplied active certificate establishes strict M>0 and Q<0 but does not certify the full 1e-5 optimization margin uniformly. We preserve that distinction; no stronger certificate is claimed here.

## Results
| FCT (s) | gamma_clear (deg) | nu_clear | Inside widened rectangle | Saved numerical recovery |
|---:|---:|---:|---|---|
| 0.05 | 22.700162 | 1.004427792 | Yes, inactive | RECOVERED |
| 0.10 | 28.302907 | 1.007857540 | Yes, active | RECOVERED |
| 0.14 | 34.737983 | 1.009925119 | Yes, active | RECOVERED |
| 0.15 | 36.564048 | 1.010359418 | Yes, active | NOT_RECOVERED within 8 s |
| 0.25 | 58.251734 | 1.013414523 | No | NOT_RECOVERED within 8 s |

All five were outside the original small rectangles. The first four now also pass the bounded pointwise contraction-rate test with their new branch metric. At FCT=0.25, metric extrapolation has max eig(M)=12.9519>10 and max eig(Sym(F))=0.3969>0.
The saved numerical CCT is 0.141234375 s, bracket [0.14084375,0.141625] s. Existing 8/16-second endpoint checks are retained. This is a numerical finite-horizon threshold, not proof of global stability/instability.

FCT=0.15 demonstrates why clearing-point membership does not imply recovery: the trajectory need not remain in its certified rectangle, and the two metrics do not establish a common switching certificate. The containment figure uses solid links only when both stored endpoints are inside the correct-mode rectangle; other links are dashed. Percent time inside is the original trapezoidal estimate on saved, nonuniform samples, not a continuous-time containment proof.

## Figures and tables
- 01_gamma_vs_time, 02_nu_vs_time: original trajectories and recovery labels.
- 03_phase_plane: original fault-on and postfault trajectories, with an inset.
- 04_contraction_zoom: clearing states and widened rectangles.
- trajectory_containment: mode-aware solid/dashed containment and full-trajectory inset.
Every figure is exported as PNG, PDF and MATLAB FIG.
- clearing_state_comparison.csv includes gamma in radians/degrees, nu, omega in rad/s, Hz, old/new membership, metric eigenvalues/conditioning and original numerical recovery.
- containment_summary.csv and containment_data.mat retain sampled-membership details.
- cct_reference.csv, region_bounds.csv and validation/audit files document the unchanged inputs and selected results.

The plotting and containment routines were adapted from the existing October 11 workflows. The metric_at, coupled, raw and rebuild functions are reused from the original FCT script. The main CCT and contraction-region scripts were not modified.
