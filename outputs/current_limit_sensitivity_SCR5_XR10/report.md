# Imax sensitivity: SCR5/XR10

Network-consistent two-state ideal radial-current benchmark. Hypothetical analytical ratings; saved Simulink setting is not changed. H=3.4, K=20, Kd=5, D=25, E/network/metric unchanged. No metric optimization, sag or trajectory study is included. Full period -180..180 deg scanned at 0.002 deg followed by root refinement. Exact switches are excluded from ordinary smooth certificates.

| Imax pu | Positive switch deg | Positive local-stable width deg | Positive certified width deg | Full-period certified intervals | Max required voltage pu |
|---|---:|---:|---:|---:|---:|
| 1.2 | 16.340926038 | 16.340926038 | 0.000000000 | 0 | 1.242006458 |
| 1.3 | 17.714020567 | 17.714020567 | 0.000000000 | 0 | 1.262890687 |
| 1.4 | 19.089626698 | 19.089626698 | 0.000000000 | 0 | 1.283825609 |

Increasing Imax delays saturation and enlarges the positive-angle locally stable interval in these tests. Certification is evaluated separately, both with and without the 1.3-pu voltage-capability check. See all_intervals.csv for every connected interval and switching_and_voltage_boundaries.csv for refined switching/voltage boundaries. The voltage limit is not simulated; any exceedance invalidates the current-only realization assumption.

The unchanged baseline R1 starts at 22.352378400 deg and requires 1.636365934 pu converter current there. Raising Imax to 1.3/1.4 alone need not expose this R1: the unchanged metric must still be tested on each branch. No current rating is recommended as a hardware change by this analytical experiment.

## Validation

Imax=1.2

- inactiveBaseline: 2.22044605e-16
- powerContinuity: 0
- voltageContinuity: 5.55111512e-17
- currentLimit: 2.22044605e-16
- Pe_gamma: 2.78620682e-10
- Pe_U: 1.81458612e-09
- Pe_gammaU: 1.08933711e-08
- projectionDerivative: 1.08897946e-10
- Jacobian: 5.7849725e-11
- matrixFormula: 8.8817842e-16
- matrixSymmetry: 0
- networkResidual: 9.42464265e-16
- equilibriumResidual: 0
- nominalMatch: 7.17648163e-13

Imax=1.3

- inactiveBaseline: 2.22044605e-16
- powerContinuity: 1.98641104e-12
- voltageContinuity: 4.6386889e-13
- currentLimit: 2.22044605e-16
- Pe_gamma: 2.78620682e-10
- Pe_U: 1.94572308e-09
- Pe_gammaU: 1.18113288e-08
- projectionDerivative: 1.17972465e-10
- Jacobian: 5.7849725e-11
- matrixFormula: 7.77156117e-16
- matrixSymmetry: 0
- networkResidual: 9.99200722e-16
- equilibriumResidual: 0
- nominalMatch: 7.17648163e-13

Imax=1.4

- inactiveBaseline: 2.22044605e-16
- powerContinuity: 4.4408921e-16
- voltageContinuity: 1.11022302e-16
- currentLimit: 2.22044605e-16
- Pe_gamma: 2.78620682e-10
- Pe_U: 1.56859548e-09
- Pe_gammaU: 6.45613096e-09
- projectionDerivative: 8.32613284e-11
- Jacobian: 5.7849725e-11
- matrixFormula: 7.77156117e-16
- matrixSymmetry: 0
- networkResidual: 9.29811783e-16
- equilibriumResidual: 0
- nominalMatch: 7.17648163e-13

