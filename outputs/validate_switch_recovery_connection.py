"""Directed interval switching checks and independent energy recovery argument.
No ODE integrations or new faults. All coefficients are loaded at full precision.
"""
import sys,json,math,time
from pathlib import Path
out=Path('outputs/radial_switch_connection_20261001')
def namespace(kind):
 script=out/(kind+'_checker.py');sys.argv=['v',str(out/('inactive_recovery_extension' if kind=='inactive' else 'active'))]
 ns={'__file__':str(script)};exec(script.read_text().split('start=time.time();')[0],ns);return ns
A=namespace('active');N=namespace('inactive');iv=N['iv'];I=N['I'];pos=N['pos'];pd=N['pd']
gs=I('0.42539661082971958');gap=I('1e-12');sw=iv.mpf([(gs-gap).a,(gs+gap).b]);ge=I('0.36012957762291908');eq=iv.mpf([(ge-gap).a,(ge+gap).b])
def enc(ns,g):
 e,reason=ns['enclosure'](g);assert e is not None,reason;return e
assert enc(N,sw.a)[2].b<I('1.2') and enc(N,sw.b)[2].a>I('1.2')
assert enc(N,eq.a)[0].b<1 and enc(N,eq.b)[0].a>1
def raw_derivative(ns,g):
 er=enc(ns,g);E=er[4];rot=iv.mpc(iv.cos(g),iv.sin(g))*ns['rot0'];v=E*rot;q=(v-ns['T'])/ns['Z'];S,Sg,SE,_,_=ns['electrical'](g,E)
 Eg=-I('.1')*Sg.imag/(1+I('.1')*SE.imag);qg=(1j*v+rot*Eg)/ns['Z'];return (ns['conj'](q)*qg).real/iv.sqrt(ns['abs2'](q))
sg1=raw_derivative(N,sw);sg2=raw_derivative(A,sw);assert pos(sg1) and pos(sg2)
def matrix(ns,g,v):return ns['matpoly']((g-ns['cen'][0])/ns['sc'][0],(v-ns['cen'][1])/ns['sc'][1],I(0),I(0))[0]
def transition(k,reverse=False):
 lo=I('0.99652713910517932');hi=I('1.0025')
 for j in range(2048):
  v=iv.mpf([(lo+(hi-lo)*j/2048).a,(lo+(hi-lo)*(j+1)/2048).b]);m1=matrix(N,sw,v);m2=matrix(A,sw,v)
  if reverse:m1,m2=m2,m1
  q=[[I(k)*m2[i][j]-m1[i][j] for j in range(2)] for i in range(2)]
  if not pd(q[0][0],q[0][1],q[1][1]):return False
 return True
gains=[]
for rev in [False,True]:
 for k in [2,4,8,16,32,64,128,256,512,1024]:
  if transition(k,rev):gains.append(k);break
 else:raise AssertionError('gain failed')
print('Validated quadratic-form gains',gains,flush=True)
# Electrical sign/feasibility certification on compact recovery box. Strict
# branch activation is imposed by the exact switch, not interval overestimation.
glo=I('0.28669381977787517');ghi=I('0.5909407406288647')
def partition(ns,left,right,count,pg_positive=False,p_positive=False):
 minpg=1e9;minvm=1e9;minrawder=1e9
 for j in range(count):
  g=iv.mpf([(left+(right-left)*j/count).a,(left+(right-left)*(j+1)/count).b]);P,Pg,r,u,E=enc(ns,g)
  if pg_positive:assert pos(Pg)
  if p_positive:assert pos(P-1)
  rd=raw_derivative(ns,g);assert pos(rd)
  vm=I('1.3')-iv.sqrt(ns['abs2'](u));assert pos(vm)
  minpg=min(minpg,float(Pg.a));minvm=min(minvm,float(vm.a));minrawder=min(minrawder,float(rd.a))
 return {'min_Pgamma':minpg,'voltage_margin':minvm,'raw_current_derivative_lower':minrawder}
phys1=partition(N,glo,sw.b,512,pg_positive=True);phys2=partition(A,sw.a,ghi,512,p_positive=True)
# Validated Riemann sums for W(g)=integral_eq^g (Pe(a)-1) da. Endpoint
# root uncertainty is bounded using a tiny interval correction.
def integral_piece(ns,a,b,count=1024):
 total=I(0);h=(b-a)/count
 for j in range(count):
  z=iv.mpf([(a+j*h).a,(a+(j+1)*h).b]);total+=h*(enc(ns,z)[0]-1)
 return total
def W(g):
 if g.b<ge.a:value=-integral_piece(N,g,ge)
 elif g.a>sw.b:value=integral_piece(N,ge,sw.a)+integral_piece(A,sw.b,g)+iv.mpf(['-1e-10','1e-10'])
 else:value=integral_piece(N,ge,g)
 return value+iv.mpf(['-1e-10','1e-10'])
wl=W(glo);wh=W(ghi);ws=W(sw.a)
mass=I('6.8');omega=100*iv.pi
def T(v):w=v-1;return mass*omega*(w*w/2+w*w*w/3)
nuLo=I('0.9975');nuHi=I('1.0025');h0=I('.0063')
assert pos(wl-h0) and pos(wh-h0) and pos(T(nuLo)-h0) and pos(T(nuHi)-h0)
# Nontrivial limiter-active clearing set, fixed decimal radian endpoints.
cLo=I('0.4267322646126135');cHi=I('0.4276056667386107') # near 24.45 / 24.50 deg
vcLo=I('.99995');vcHi=I('1.00005');wc=W(cHi);energyUpper=wc+max(T(vcLo).b,T(vcHi).b)
assert pos(cLo-sw.b) and pos(h0-energyUpper)
result={'switch_bracket_rad':[str(sw.a),str(sw.b)],'switch_raw_derivative_inactive':[float(sg1.a),float(sg1.b)],'switch_raw_derivative_active':[float(sg2.a),float(sg2.b)],'gain_active_to_inactive_quadratic':gains[0],'gain_inactive_to_active_quadratic':gains[1],'gain_nu_range':[.99652713910517932,1.0025],'physical_inactive':phys1,'physical_active':phys2,'W_lower_wall':[float(wl.a),float(wl.b)],'W_upper_wall':[float(wh.a),float(wh.b)],'W_switch':[float(ws.a),float(ws.b)],'energy_level':.0063,'kinetic_lower_wall':[float(T(nuLo).a),float(T(nuLo).b)],'kinetic_upper_wall':[float(T(nuHi).a),float(T(nuHi).b)],'recovery_box_rad':[[str(glo.a),str(nuLo.a)],[str(ghi.b),str(nuHi.b)]],'clearing_subset_rad':[[str(cLo.a),str(vcLo.a)],[str(cHi.b),str(vcHi.b)]],'clearing_subset_max_energy_upper':float(energyUpper.b),'complete':True}
(out/'connection_validation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
