# Radial-current limiter: instability diagnosis

Same saved SCR5/XR10 network and parameters; no redesign, tuning, new angle parameter, metric optimization or EMT. This is the ideal two-state radial benchmark, not the complete PI/freeze controller.

Current switch: **16.340926037578 deg**. At this point Pe and Vc are continuous, but their derivatives need not be. No smooth Pe_gamma=0 root was found on either positive-angle branch over 0..180 deg. The slope changes SIGN BY A JUMP at the switch, not by continuously crossing zero. Therefore there is no positive-width stable active interval.

## Exact one-sided values

All magnitude/power derivatives are per RADIAN of target gamma. Phase derivatives are rad/rad (numerically also deg/deg). Global phasors use the referred physical grid bus as zero angle; gamma itself references Vth.

| Quantity | Inactive side | Active side |
|---|---:|---:|
| Pe | 1.17939641482 | 1.17939641482 |
| Qe | 0.0338889230818 | 0.0338889230818 |
| Pe_gamma | 4.0415260692 | -0.16830696933 |
| Qe_gamma | 0.257170325271 | -0.690153574 |
| Iconv_mag_gamma | 4.1763396411 | -2.5442610981e-16 |
| Iconv_angle_gamma | 0.588253445882 | 0.588253445882 |
| Vc_mag_gamma | 5.56515446259e-17 | -0.163915250859 |
| Vc_angle_gamma | 1 | 0.0277757802035 |
| Vgrid_mag_gamma | -0.0604937298024 | -0.136478909543 |
| Igrid_mag_gamma | 4.1917481825 | -0.0246441589027 |
| delta_VI_gamma | 0.119488458089 | -0.580595333826 |
| C_V | -0.100874796149 | -0.164214005983 |
| C_I | 4.14240086535 | -0.00409296334672 |
| C_Vmag | -0.0725909122877 | -0.163771164121 |
| C_Imag | 4.11816631665 | -0.0242115558179 |
| C_phase | -0.00404933516533 | 0.0196757506097 |
| max_real_eigJ | -0.169473199768 | 0.00671999570409 |

## What removes synchronizing power?

Before limiting both the magnitude and direction of converter current can change. The prescribed capacitor voltage has fixed magnitude E and follows the target angle with unit angle derivative. At saturation the radial component of the current derivative is removed: dIconv/dgamma=(Imax/abs(q)) [q_gamma-q real(conj(q)q_gamma)/abs(q)^2]. The surviving derivative is tangent to the current circle. Current direction is exactly that of q; there is no independent limiter-angle control.

The capacitor voltage must then move according to Vc=Vth+Zeq Iconv rather than remain equal to Vtarget. Its magnitude decreases and its angle barely follows the target near switching. Grid current is not the converter current because shunt/magnetizing currents and transformer referral remain in the network.

At the active-side limit the COMPLEX decomposition is C_V=-0.164214006 and C_I=-0.00409296335, totaling -0.168306969 pu/rad. The POLAR decomposition is voltage magnitude -0.163771164, grid-current magnitude -0.0242115558, and V-I phase contribution 0.0196757506 pu/rad. The signed terms, rather than current magnitude clamping alone, explain the negative synchronizing slope.

Reactive-power slope changes from 0.257170325 to -0.690153574 pu/rad while active-power slope changes from 4.04152607 to -0.168306969. Positive Q is verified against Q=v_gq*i_gd-v_gd*i_gq. The stored P/Q curves distinguish reactive increase from changes in apparent-power magnitude; no conservation of P+Q or fixed grid-current magnitude is assumed.

**Near the switch the dominant negative term is grid-voltage magnitude reduction**, not V-I phase rotation: the phase contribution is initially positive. The current-magnitude growth term that previously supported P is removed and the reconstructed Vc/grid voltage can sag. At larger angles the phase term becomes negative and its importance changes, as quantified below.

| gamma deg | Pe | Qe | Voltage-magnitude contribution | Grid-current-magnitude contribution | Phase contribution |
|---|---:|---:|---:|---:|---:|
| 16.340926038 | 1.17939641 | 0.0338889231 | -0.163771164 | -0.0242115558 | 0.0196757506 |
| 20.000000000 | 1.16817297 | -0.00878216503 | -0.155949432 | -0.0225852899 | -0.00489384287 |
| 40.000000000 | 1.08835817 | -0.216780859 | -0.139151775 | -0.0180660053 | -0.117162645 |
| 60.000000000 | 0.977172794 | -0.399766109 | -0.12317021 | -0.0142707036 | -0.223821744 |
| 90.000000000 | 0.75794722 | -0.627578394 | -0.0870130545 | -0.00834120787 | -0.376078657 |

On the active branch Q first DECREASES from positive injection, crosses zero at 19.229481506893 deg, and subsequently becomes negative (reactive absorption under this convention). Thus the initial loss of P slope is NOT explained by an increase in positive reactive-power injection. The power-factor cosine initially improves; farther along the active branch growing phase separation contributes to declining P.

## Stability versus certificate

At each point Pstar=Pe and nu=1. trace(J)=-25/6.8 and det(J)=Pe_gamma/6.8. Thus a negative slope gives a saddle regardless of the chosen contraction metric. The ordinary Jacobian is undefined at the exact switch; the one-sided stable and saddle classifications meet there. On the positive principal interval 0..180 deg: inactive locally stable 0<=gamma<16.340926037578; active locally stable EMPTY; active saddle 16.340926037578<gamma<=180. Displayed plots cover 0..90 and 10..40 deg.

The unchanged metric already fails to certify the locally stable pre-switch branch (including nominal): its positive maximum symmetric eigenvalue must not be confused with local instability. At the switch genuine local instability occurs in addition. The negative limited power slope is a property of this electrical equilibrium family; the earlier pre-switch certificate failure is metric conservativeness.

## Validation

- activeCurrentMagnitude: 4.4408921e-16
- radialAngleDifference_rad: 2.00928071e-16
- complexPower: 0
- polarPower: 4.4408921e-16
- complexDecomposition: 2.66453526e-15
- polarDecomposition: 3.83157048e-15
- dqQSign: 0
- activeCurrentMagnitudeDerivative: 2.15105711e-15
- complexDerivatives: 1.07774893e-09
- magnitudeDerivatives: 5.22532861e-09
- phaseDerivatives: 2.07879857e-07
- PeDerivative: 4.71746142e-10
- QeDerivative: 6.33649577e-10
- Jacobian: 1.67739749e-10
- equilibriumResidual: 0
- switchPowerContinuity: 2.44249065e-15
- switchVoltageContinuity: 4.99600361e-16
- oneSidedPeDerivative: 1.42912987e-10

CSV includes Cartesian phasors and derivatives, continuous unwrapped global angles, relative V-I angle, polar derivatives, both power-slope decompositions, grid P/Q, voltage capability, Jacobian spectrum and metric context. Exact one-sided diagnostics are separately exported; smooth derivative columns at the exact switch are marked NaN. MAT retains all figure source data.
