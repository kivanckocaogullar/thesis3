# Fixed-Vset Q-V consistency of low-E radial-limiter candidates

The limiter-active, locally asymptotically stable interval survives for BOTH candidates in the network-consistent, quasi-static Q-V ideal-current-reference benchmark. This is not a full PI-controller equilibrium/stability validation and not a contraction certificate.

Frozen: external SCR5/XR10, physically reconstructed Xeq=.350 pu downstream network, H=3.4 K=20 Kd=5 D=25, radial current-reference Imax=1.2, radial voltage-command maximum=1.3, converter impedance and active damping unchanged. The voltage limiter is checked for feasibility, not applied while retaining an inconsistent ideal current source. No old files were overwritten; all 59 protected hashes remain unchanged.

## Compact comparison

| Nominal E candidate | Fixed Vset | Old fixed-E stable active interval, deg | New fixed-Vset/Q-V stable active interval, deg | Pe at activation | Total Pe_gamma just after activation, pu/rad |
|---|---|---|---|---|---|
| .875 | .840735476598092 | (24.994574206146,28.493771) | (24.994829900309,28.525361699) | 1.022082859154 | .071536240557 |
| .890 | .859466803140713 | (25.040119361334,26.628001) | (25.040346077641,26.620648610) | 1.045048915123 | .029842125422 |

In both rows the lower endpoint is radial current activation. The upper endpoint is TOTAL Pe_gamma=0, where the slow eigenvalue becomes zero. Both are excluded from the strictly smooth/stable interval. The limiter remains active beyond these upper endpoints up to 180 deg in the principal positive-angle scan, but those points are not included as stable feasible exporting operating conditions. No bridge across the nonsmooth activation boundary is certified.

## Voltages and nominal feasibility

| Quantity | .875 candidate | .890 candidate |
|---|---|---|
| Nominal Pstar | 1 pu | 1 pu |
| Nominal gamma, deg | 24.440460105138 | 23.930462925388 |
| Nominal converter current, pu | 1.175971650606 | 1.150972545726 |
| Nominal required command, pu | 1.102896541595 | 1.114881190872 |
| Nominal capacitor magnitude, pu | .875000000000 | .890000000000 |
| Nominal PCC magnitude, pu | .904352701456 | .913338952897 |
| Active Eref range, pu | [.874737992,.875060114] | [.890057253,.890124086] |
| Active physical capacitor magnitude, pu | [.875060114,.876396079] | [.890124086,.890395701] |
| Active PCC magnitude, pu | [.903512394,.904232759] | [.911662498,.911809999] |
| Active required command, pu | [1.108554268,1.110134039] | [1.126238483,1.126554790] |
| Active converter current | 1.2 pu | 1.2 pu |
| Active equilibrium Pstar range, pu | (1.022082859,1.024173361) | (1.045048915,1.045450831) |

Nominal points remain unsaturated and satisfy current/command limits. Active current margin is zero by construction; active voltage margins are approximately .189866 and .173445 pu. Capacitor-node and PCC voltages are distinct reconstructed physical voltages; Eref is not forced onto the capacitor during saturation. The retained bases are LV 690 V and HV 20.4 kV. These low voltage operating levels are not a claim of acceptable practical voltage regulation: no minimum capacitor/PCC voltage requirement was supplied. A requirement |Vc|>=.9 pu would exclude both candidates.

## Equations and derivative

Let theta=gamma+arg(Vth). The ideal raw PI-reference proxy is q=(Eref*exp(j theta)-Vth)/Zeq. Apply s=q if |q|<=1.2 and s=1.2*q/|q| otherwise. Radial scaling commutes with the VSM dq/network-frame rotation, so this exactly represents the same radial rule within the proxy. It does not replace the actual PI output with a verified EMT measurement.

Reconstruct Vc=Vth+Zeq*s, Vgrid=vgE*Vc+vgU*U, Igrid=igE*Vc+igU*U, Pe+jQgrid=Vgrid*conj(Igrid), u_required=Vc+(Zconverter+.2)*s.

For each fixed Vset solve h(Eref,gamma)=Eref+.1*Qgrid(Eref,gamma)-Vset=0. Qstar=0. The target Eref, including during saturation, is the unknown; no artificial constant-E continuation is used.

On each smooth branch:
Eref_gamma=-.1*Q_gamma/(1+.1*Q_E),
Pe_total_gamma=P_gamma+P_E*Eref_gamma.

All partial derivatives are analytic derivatives of the reconstructed network and radial limiter. Across 0..180 deg the minimum implicit denominator was .99988249 and .99988593; maximum Q-V residual was 2.22e-16 pu. A separate root scan at 181 angles found exactly one positive algebraic solution in Eref=.05..2 pu for each angle. This is numerical branch/uniqueness evidence over that tested range, not a global formal uniqueness proof. Total power derivatives were independently checked against finite differences of the coupled solution (maximum error 1.16e-8 and 3.08e-10 pu/rad). Network KCL was checked to 1e-11.

The normalized reduced dynamics remain gamma_dot=nu-1, nu_dot=((Pstar-Pe(gamma;Vset))/nu-25*(nu-1))/6.8. At each equilibrium nu0=1, Pstar=Pe(gamma0;Vset). Pstar remains fixed when linearizing; its family parametrization is NOT differentiated inside the dynamics.

J0=[[0,1],[-Pe_total_gamma/6.8,-25/6.8]].

| Candidate | Midpoint gamma, deg | Fixed-E partial Pg | Eref_gamma | TOTAL Pg | Jacobian eigenvalues |
|---|---|---|---|---|---|
| .875 | 26.760095799786 | .040474791839 | -.005095447264 | .033008459660 | -.001320812903, -3.67514977533 |
| .890 | 25.830497343630 | .017838667804 | -.002395021704 | .014397548159 | -.000575992167, -3.67589459607 |

Both eigenvalues are negative throughout each OPEN reported interval. The slow eigenvalue tends to zero at the upper boundary. The post-switch total derivatives are smaller than the old fixed-E switch derivatives (.088910761 and .037583663), although the boundary shifts are small.

## What is and is not physically verified

Vset is held fixed at the value calibrated from each old nominal point. It is not adjusted again at each gamma. This removes the previous DC Q-V inconsistency in the ideal benchmark. Because saturation is radial, the physical capacitor voltage phase/magnitude need not equal the rotating target; gamma continues to mean the VSM reference angle relative to Vth.

The actual controller contains a Q-V filter (TQ=.2 s), capacitor-voltage/current PI integrators, feedforward, delayed conditional freeze/release, measurement filters, and network electrical dynamics. Here the Q-V filter has been replaced by its DC algebraic relation and the raw PI reference by q. Thus the calculated two-state eigenvalues do NOT certify the full dynamic PI/Q-V controller. A real steady PI/freeze state satisfying its equations has not been constructed. In particular, a total equilibrium-family derivative is the Jacobian stiffness only of this instantaneous algebraic reduction, not automatically of a model retaining the Q-V state.

Conclusion: both candidates survive the requested fixed-Vset quasi-static consistency test and remain potential starting points for reduced-model contraction analysis, subject to explicit acceptance of the low voltages and model approximation. Full controller consistency still requires PI/freeze equilibria and the Q-V dynamic state; it cannot be inferred from this two-state result. No M optimization, faults, EMT or saved Simulink modification was performed.
