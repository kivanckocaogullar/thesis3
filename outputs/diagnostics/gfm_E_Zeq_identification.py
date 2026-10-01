"""First-principles balanced 50-Hz identification; never edits contraction files.
Parameters transcribed from GFMVSG1.slx and its resolved tDMSquared libraries.
No EMT curve fitting, parameter optimization, or contraction calculation.
"""
from pathlib import Path
import hashlib
import json
import numpy as np

OUT = Path(__file__).resolve().parent.parent
S, VL, VH, F, TS = 2e6, 690., 20400., 50., 1e-4
W = 2*np.pi*F
ZBL, ZBH = VL**2/S, VH**2/S
U = 20000/VH
ZT = .005+.02j
ZS = (4/np.sqrt(101)+1j*40/np.sqrt(101))/ZBH
ZL = (2*.501+1j*W*2*2.279e-3)/ZBH
ZBR = (.01*1e6/(.01+1e6))/ZBH  # closed breaker, parallel snubber
YM = 1/500+1/(1j*500)
YC = 1/(.06875+1/(1j*.04569))
YH = 1j*W*(2*151.2e-9)*ZBH/2
YFAULT = ZBH/1e6  # open fault phase-switch snubber, positive sequence
ZF = .001+.05097j

def network(v, u=U, shunts=True):
    """Nodes: core, HV before breaker, line sending end, line receiving end.
    All quantities referred to the 690/20400 V, 2 MVA transformer bases.
    Ideal transformer ratio and D1 phase shift removed by referral.
    """
    y=np.zeros((4,4),complex)
    for i,j,z in [(0,1,ZT),(1,2,ZBR),(2,3,ZL)]:
        adm=1/z
        y[i,i]+=adm; y[j,j]+=adm; y[i,j]-=adm; y[j,i]-=adm
    y[0,0]+=1/ZT+(YM if shunts else 0)
    y[2,2]+=YH if shunts else 0
    y[3,3]+=1/ZS+(YH+YFAULT if shunts else 0)
    b=np.zeros(4,complex); b[0]=v/ZT; b[3]=u/ZS
    n=np.linalg.solve(y,b)
    ip=(v-n[0])/ZT
    ig=(n[0]-n[1])/ZT
    ic=ip+(YC*v if shunts else 0)
    return dict(ic=ic, ip=ip, ig=ig, vg=n[2], sg=n[2]*np.conj(ig), n=n,
                nodal_residual=float(max(abs(y@n-b))))

def solve2(fun, guess):
    x=np.array(guess,float)
    for _ in range(50):
        r=np.asarray(fun(x))
        if max(abs(r))<2e-11:
            return x, float(max(abs(r)))
        h=1e-6
        jac=np.column_stack([(np.asarray(fun(x+h*e))-np.asarray(fun(x-h*e)))/(2*h)
                             for e in np.eye(2)])
        x-=np.linalg.solve(jac,r)
    raise RuntimeError('Phasor equilibrium did not converge')

def cplx(z):
    return dict(real=float(np.real(z)),imag=float(np.imag(z)),
                magnitude=float(abs(z)),angle_deg=float(np.angle(z,deg=True)))

a=network(1,0); b=network(0,1)
ZTH=1/a['ic']; UTH=-b['ic']/a['ic']*U
# Independent ladder/Thevenin reduction; compare against nodal elimination.
za=1/(1/ZS+YH+YFAULT); ua=U/(1+ZS*(YH+YFAULT))
zb=1/(YH+1/(ZL+za)); ub=ua/(1+YH*(ZL+za))
zm=1/(YM+1/(ZT+ZBR+zb)); um=ub/(1+YM*(ZT+ZBR+zb))
zp=ZT+zm
zcheck=zp/(1+YC*zp); ucheck=um/(1+YC*zp)
assert abs(ZTH-zcheck)<1e-11 and abs(UTH-ucheck)<1e-11
alpha=(1-TS*np.pi*2500)/(1+TS*np.pi*2500)
zminus=np.exp(-1j*W*TS); kh=np.tan(np.pi*2*TS)
G=(1-zminus)/(1+kh+(kh-1)*zminus)*(1-alpha)/2*(1+zminus)/(1-alpha*zminus)

def operating_point(pstar, gain):
    k=abs(gain)
    def residual(x):
        e,g=x; n=network(e*np.exp(1j*g))
        return [k*k*n['sg'].real-pstar, k*e-(1-.1*k*k*n['sg'].imag)]
    (e,g),res=solve2(residual,[1.,pstar*.068])
    v=e*np.exp(1j*g); n=network(v)
    bridge=v+ZF*n['ic']
    uraw=bridge+.2*gain*n['ic']
    # Thevenin gamma is measured relative to UTH, not raw HV source phase.
    d=g-np.angle(UTH)
    arg=(e*e+abs(UTH)**2-(1.2/k)**2*abs(ZTH)**2)/(2*e*abs(UTH))
    assert -1<arg<1
    glim=np.arccos(arg)
    energy_in=(v*np.conj(n['ic'])).real
    lossgrid=energy_in-n['sg'].real
    simple=(v-U)/(2*ZT+ZBR+ZL+ZS)
    return dict(Pstar_pu=pstar,E_pu=float(e),Eref_pu=float(k*e),
        gamma_bus_deg=float(np.rad2deg(g)),gamma_thevenin_deg=float(np.rad2deg(d)),
        gamma_vsg_referred_bus_deg=float(np.rad2deg(g)+np.angle(gain,deg=True)),
        Ireq=cplx(n['ic']),Iref_pu=float(k*abs(n['ic'])),grid_power=cplx(n['sg']),
        measured_Pe_pu=float(k*k*n['sg'].real),measured_Qe_pu=float(k*k*n['sg'].imag),
        capacitor_injection_power_pu=float(energy_in),passive_power_loss_pu=float(lossgrid),
        bridge_voltage_pu=float(abs(bridge)),raw_voltage_reference_pu=float(abs(uraw)),
        gamma_lim_frozen_E_deg=float(np.rad2deg(glim)),arccos_argument=float(arg),
        residual=res,naive_series_I_pu=float(abs(simple)),
        naive_series_complex_current_error_pct=float(100*abs(simple-n['ic'])/abs(n['ic'])))

nominal=operating_point(1.,1.)
corrected=operating_point(1.,G)
def threshold(gain):
    k=abs(gain)
    def fun(x):
        e,g=x; n=network(e*np.exp(1j*g))
        return [k*abs(n['ic'])-1.2,k*e-(1-.1*k*k*n['sg'].imag)]
    (e,g),res=solve2(fun,[1.,.08])
    n=network(e*np.exp(1j*g)); v=e*np.exp(1j*g)
    return dict(E_pu=float(e),gamma_bus_deg=float(np.rad2deg(g)),
        gamma_lim_deg=float(np.rad2deg(g-np.angle(UTH))),
        gamma_vsg_referred_bus_deg=float(np.rad2deg(g)+np.angle(gain,deg=True)),
        Pstar_pu=float(k*k*n['sg'].real),
        raw_voltage_reference_pu=float(abs(v+ZF*n['ic']+.2*gain*n['ic'])),residual=res)

validation=0.
for e in [.98,1.,1.02]:
    for deg in [-10,0,4,6,20]:
        v=e*np.exp(1j*np.deg2rad(deg)); n=network(v)
        validation=max(validation,abs(n['ic']-(v-UTH)/ZTH),n['nodal_residual'])
        assert (v*np.conj(n['ic'])).real-n['sg'].real>=-1e-10
assert validation<1e-10
summary=dict(
    assumptions='Balanced positive sequence, 50 Hz, normal closed breaker, inactive limits; quasi-steady control tracking and algebraic Q-V droop. Not EMT validation.',
    selected_E_node='Physical capacitor/filter node voltage; measured-reference tracking approximation stated separately.',
    E_pu=nominal['E_pu'],U1_pu=float(abs(UTH)),Req_pu=float(ZTH.real),Xeq_pu=float(ZTH.imag),
    Zeq_mag_pu=float(abs(ZTH)),Zeq_angle_deg=float(np.angle(ZTH,deg=True)),Imax_pu=1.2,
    gamma_lim_deg=nominal['gamma_lim_frozen_E_deg'],
    gamma_lim_definition='Frozen nominal E, gamma relative to Thevenin voltage, unity measurement approximation. See threshold_algebraic_QV for the variable-E closure.',
    common_Sbase=S,common_Vbase=VL,common_Ibase=S/(np.sqrt(3)*VL),common_Zbase=ZBL,
    U1_physical_bus_referred_pu=U,Thevenin_source= cplx(UTH),Zeq=cplx(ZTH),
    passive_series_sum=cplx(2*ZT+ZBR+ZL+ZS),measurement_gain=cplx(G),
    physical_Imax_with_measurement_gain=1.2/abs(G),
    base_HV_transformer=dict(V=VH,I=S/(np.sqrt(3)*VH),Z=ZBH),
    base_HV_measurement=dict(V=20000.,I=S/(np.sqrt(3)*20000),Z=20000**2/S),
    nominal_unity_measurement=nominal,nominal_static_measurement_correction=corrected,
    threshold_algebraic_QV=threshold(1.),threshold_static_measurement_correction=threshold(G),
    saved_model_setpoints=[operating_point(ps,1.) for ps in [.1,.35,.4]],
    components=dict(ZT1=cplx(ZT),ZT2=cplx(ZT),ZS=cplx(ZS),ZL=cplx(ZL),ZBR=cplx(ZBR),
        YM=cplx(YM),YC=cplx(YC),Yline_half=cplx(YH),Yfault_pu=float(YFAULT),ZF=cplx(ZF)),
    grid_voltage_coeff_E=cplx(a['vg']),grid_voltage_coeff_U=cplx(b['vg']),
    grid_current_coeff_E=cplx(a['ig']),grid_current_coeff_U=cplx(b['ig']),
    max_nodal_and_thevenin_residual=float(validation),
    script_hashes={name:hashlib.sha256((OUT/name).read_bytes()).hexdigest()
                   for name in ['baseline_gfm_contraction.m','gfm_current_limiter_analysis.m']})
(OUT/'gfm_E_Zeq_identification_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
console='Previously assumed: E=1 pu; Xeq=0.3 pu; Imax=1.2 pu; gamma_lim=20.739519610955 deg\n'
console+=f"Identified 50-Hz reduction: E_nom={nominal['E_pu']:.12f}; Zeq={ZTH.real:.12f}+j{ZTH.imag:.12f} pu\n"
console+=f"Base: 2 MVA, 690 V LL; Imax=1.2 pu of converter current reference. Ubus_referred={U:.12f}; Uth={abs(UTH):.12f} pu\n"
console+=f"gamma_lim: frozen nominal E={nominal['gamma_lim_frozen_E_deg']:.12f}; algebraic Q-V={threshold(1.)['gamma_lim_deg']:.12f} deg\n"
console+=f"Static measurement correction: E_nom={corrected['E_pu']:.12f}; gamma_lim_QV={threshold(G)['gamma_lim_deg']:.12f} deg\n"
console+='Recommendation: structurally revise electrical algebraic closure and power measurement, then rerun; do not change only Xeq.\n'
console+='No contraction calculation rerun; both analysis scripts untouched.\n'
(OUT/'gfm_E_Zeq_identification_console.txt').write_text(console,encoding='utf-8')
print(console)
print(json.dumps(summary,indent=2))
