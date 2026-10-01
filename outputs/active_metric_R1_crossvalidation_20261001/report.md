# Exact active metric cross-validation on extended inactive R1

The unchanged, full-precision coefficients, centre and scaling were loaded directly from `../radial_quadratic_robust_20261001_161903/metric.json`. The saved nominal network was loaded from the existing transient study. No previous files were edited. Pstar=1 and Vset=1 are fixed; the nominal external-source value in the saved per-unit convention is U=0.98039215686274506. This value is not a newly imposed voltage sag.

The inactive quasi-static Q-V equation E+0.1 Qgrid=1 was solved at each angle. Its implicit derivative E_gamma=-0.1 Q_gamma/(1+0.1 Q_E) was included in the total Pe derivative. Physical-time dynamics are f1=100*pi*(nu-1), f2=((1-Pe)/nu-25*(nu-1))/6.8. Both spatial metric derivatives and their complete contribution Mdot=M_gamma*f1+M_nu*f2 were included. Finite-difference cross-check errors: Mdot 1.14e-6; Jacobian 1.84e-7.

| Check | Original R2, active | Extended R1, inactive |
|---|---|---|
| Same M positive definite? | Uniformly validated, M > 1e-5 I | All 1,002,001 samples pass; uniform positivity not independently certified |
| Contraction at lambda=0.1 /s? | Uniformly validated, Qc < -1e-5 I | Fails; sampled worst eigenvalue +0.281560808388 |
| Contraction at lambda=0? | Follows from the positive-rate certificate | Fails; sampled worst eigenvalue +0.277628936340 |
| Uniform full-domain validation? | Existing directed-interval certificate passed | No; numerical counterexamples, so no full-domain interval attempt |

The tested R1 rectangle is gamma=[16.426345885756,24.373430419712] degrees, nu=[0.99652713910517932,1.0025], with a 1001-by-1001 grid including all endpoints. The sampled minimum eigenvalue of M is 6.474509390975e-5, at (22.895272696396187 degrees, 0.99652713910517932). No positivity failure was found. This sampled minimum is not asserted to be an exact global bound.

For both rates an explicit contraction counterexample, also the worst sampled point, is (gamma,nu)=(16.426345885756 degrees,0.99652713910517932). At lambda=0.1, Qc=[[0.146374181591375,-4.973410379050178],[-4.973410379050178,-182.68632132962028]]. Its largest eigenvalue is positive. 635123 grid points fail at lambda=0.1 and 630933 fail at lambda=0. These are grid counts, not area estimates or uniform proofs.

## One-sided switching-boundary check

The switch is at gamma=24.3734304197121 degrees. Pe=1.17309010862117 on both sides, while total Pe_gamma is 2.61821792763416 on the inactive side and -0.140081398881124 on the active side. The vector field is continuous; the Jacobian is not. For 10001 nu samples across the requested R1 strip:

| One-sided branch | Worst eig(Qc), lambda=0.1 | Worst eig(Qc), lambda=0 |
|---|---:|---:|
| Inactive | +0.273102936342 | +0.260722485430 |
| Active | -0.00190389543744 | -0.00193047548055 |

All four extrema occur at nu=0.99652713910517932. The same metric remains positive on the sampled switching line, minimum eigenvalue 9.867080593384e-5. Thus using active equations at the boundary would wrongly hide the inactive-side failure. The boundary results here are numerical one-sided checks, not a new switched-system certificate.

Raw current spans 0.811587095124287 to 1.199999999999995 pu, consistent with the inactive branch including its closure. Maximum required voltage command is 1.248880365016437 pu, below 1.3 pu by 0.051119634983563 pu.

**Conclusion:** This exact M_active cannot certify both requested branch domains. The demonstrated obstruction is failure of the complete contraction inequality even at lambda=0, not an observed loss of positive definiteness. Uniform R1 positivity remains unproved because the requested full-domain interval stage was conditional on successful numerical screening. No coefficients, model parameters or existing certificates were changed, and no metric optimization was performed.
