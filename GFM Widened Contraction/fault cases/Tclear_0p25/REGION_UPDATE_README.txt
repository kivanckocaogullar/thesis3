TCLEAR=0.25 s — WIDENED CONTRACTION REGIONS

The simulation is configured for Tclear=0.25 s. All 12 retained-voltage cases are recalculated using the existing model, limiter, integration settings and recovery criteria. The standard tight-tolerance comparison and affected-case extra-tight verification are retained.

Region source: Documents/MATLAB2/GFM Widened Contraction/contraction region/metrics.
Inactive: gamma [10,24.3234] deg, nu [0.98,1.02].
Active: gamma [24.4234,43] deg, nu [0.98,1.02].
M<10I and original contraction rates remain unchanged. The active certificate establishes strict inequalities, not a uniform 1e-5 margin.

inside_Cpost: active rectangle with active postfault limiter mode.
inside_inactive_Cpost: inactive rectangle with inactive mode.
inside_any_Cpost: union of the two correct-mode memberships, excluding the switching gap. Membership does not guarantee recovery.

CSV and MAT results, trajectories, events, verification tables, summaries and all 18 PNG/FIG plots use Tclear=0.25 s. The plots read the actual duration from the data.

To rerun simulations: fixed_Pstar_recovery_SCR5_XR10_batch.m, then verify_affected_cases.m.
To update saved region classification and regenerate plots: update_saved_region_membership.m.
To regenerate plots only: plot_fixed_Pstar_nominal_sag_sweep_SCR5_XR10().
The previous working state was backed up outside this folder.
See RUN_STATUS.txt and widened_region_update_audit.json for final verification status.
