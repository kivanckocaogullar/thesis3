"""Directed interval validation of a saved widened contraction rectangle.
Based on original_scripts/{inactive,active}_checker.py (same closure, same
implicit Q-V root enclosure, same mean-value fallback for Q). Changes:
 * one file for both branches (argv[2] = inactive | active);
 * the limiter mode is ENFORCED per box (|q|<Imax inactive, |q|>Imax active);
 * every condition is first tried with the epsilon margin used in the
   optimization and then with the strict (zero-margin) condition;
 * certified eigenvalue bounds are accumulated over passed boxes;
 * failing boxes are tested at their centre to separate genuine counterexamples
   from interval overestimation.
Usage: python widened_interval_checker.py metrics/<branch> <branch> [box budget]
Writes interval_verification.json into that folder (overwrites only that file).
"""
import sys,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'deps'))
from mpmath import iv
iv.dps=30
base=Path(sys.argv[1]);branch=sys.argv[2];assert branch in ('inactive','active')
ACTIVE=branch=='active'
data=json.loads((base/'metric.json').read_text(),parse_float=str)
I=lambda x:iv.mpf(str(x))
def pos(x):return bool(x.a>0)
def neg(x):return bool(x.b<0)
def conj(x):return iv.mpc(x.real,-x.imag)
def abs2(x):return x.real*x.real+x.imag*x.imag
def cx(pair):return iv.mpc(I(pair[0]),I(pair[1]))
net={k:cx(v) for k,v in data['network'].items()}
C=[[[I(v) for v in row] for row in mat] for mat in data['coefficients']]
cen=list(map(I,data['center']));sc=list(map(I,data['scale']));bounds=[[I(v) for v in row] for row in data['bounds']]
U=I(data['U']);Z=net['Zeq'];T=net['Vth'];tnorm=iv.sqrt(abs2(T));rot0=T/tnorm
mass=I(data['mass']);D=I(data['D']);om=100*iv.pi;ilim=I(data['Imax']);vlim=I(data['umax']);AD=I(data['AD']);lam=I(data['lambda']);eps=I(data['epsilon']);up=I(data['upper'])

def electrical(g,E):
 rot=iv.mpc(iv.cos(g),iv.sin(g))*rot0;v=E*rot;q=(v-T)/Z;rg=1j*v/Z;rE=rot/Z;r=iv.sqrt(abs2(q))
 if ACTIVE:
  ic=ilim*q/r
  ig=ilim*(rg/r-q*(conj(q)*rg).real/(r*r*r));iE=ilim*(rE/r-q*(conj(q)*rE).real/(r*r*r))
 else:
  ic=q;ig=rg;iE=rE
 vc=T+Z*ic;vg=net['vgE']*vc+net['vgU']*U;ii=net['igE']*vc+net['igU']*U
 S=vg*conj(ii);Sg=net['vgE']*Z*ig*conj(ii)+vg*conj(net['igE']*Z*ig);SE=net['vgE']*Z*iE*conj(ii)+vg*conj(net['igE']*Z*iE)
 return S,Sg,SE,r,vc+(net['Zconverter']+AD)*ic

_enc={}
def enclosure(g):
 key=g._mpi_
 if key not in _enc:_enc[key]=_enclosure(g)
 return _enc[key]
def _enclosure(g):
 E=iv.mpf(['.98','1.05']);S0=electrical(g,I('.98'))[0];S1=electrical(g,I('1.05'))[0]
 if not (neg(I('.98')+I('.1')*S0.imag-1) and pos(I('1.05')+I('.1')*S1.imag-1)):return None,'root bracket'
 hE=1+I('.1')*electrical(g,E)[2].imag
 if not pos(hE):return None,'root uniqueness'
 for _ in range(10):
  mid=(E.a+E.b)/2;S=electrical(g,mid)[0];he=1+I('.1')*electrical(g,E)[2].imag
  if not pos(he):return None,'root derivative'
  N=mid-(mid+I('.1')*S.imag-1)/he
  aa=max(E.a,N.a);bb=min(E.b,N.b)
  if aa>bb:return None,'empty root'
  new=iv.mpf([aa,bb]);oldwidth=float(E.delta.b);E=new
  if float(E.delta.b)>.995*oldwidth:break
 S,Sg,SE,r,u=electrical(g,E);he=1+I('.1')*SE.imag
 if not pos(he):return None,'implicit derivative'
 Pg=Sg.real-SE.real*I('.1')*Sg.imag/he
 return (S.real,Pg,r,u,E),None

def matpoly(z1,z2,dz1,dz2):
 b=[I(1),z1,z2,z1*z1,z1*z2,z2*z2]
 bd=[I(0),dz1,dz2,2*z1*dz1,z2*dz1+z1*dz2,2*z2*dz2]
 M=[[sum(C[k][i][j]*b[k] for k in range(6)) for j in range(2)] for i in range(2)]
 MD=[[sum(C[k][i][j]*bd[k] for k in range(6)) for j in range(2)] for i in range(2)]
 return M,MD

import math
def dn(x):return math.nextafter(float(x),-math.inf)
def upf(x):return math.nextafter(float(x),math.inf)
def pd(a,b,c):return pos(a) and pos(a*c-b*b)
def eigs(a,b,c):
 """Enclosures of (lambda_min, lambda_max) of [[a,b],[b,c]] over the interval matrix."""
 h=(a+c)/2;w=((a-c)/2)**2+b*b
 # The discriminant is a sum of squares; clip only the spurious negative part of its enclosure.
 if w.a<0:w=iv.mpf([0,w.b])
 s=iv.sqrt(w);return h-s,h+s

_sec={}
def second(g,E):
 key=(g._mpi_,E._mpi_)
 if key not in _sec:_sec[key]=_second(g,E)
 return _sec[key]
def _second(g,E):
 rot=iv.mpc(iv.cos(g),iv.sin(g))*rot0;v=E*rot;q=(v-T)/Z;r=iv.sqrt(abs2(q));qi=[1j*v/Z,rot/Z];qij=[[-v/Z,1j*rot/Z],[1j*rot/Z,iv.mpc(0)]]
 ri=[(conj(q)*u).real/r for u in qi]
 if ACTIVE:ci=[ilim*(qi[i]/r-q*ri[i]/(r*r)) for i in range(2)];c=ilim*q/r
 else:ci=qi;c=q
 vc=T+Z*c;vg=net['vgE']*vc+net['vgU']*U;ig=net['igE']*vc+net['igU']*U
 si=[net['vgE']*Z*ci[i]*conj(ig)+vg*conj(net['igE']*Z*ci[i]) for i in range(2)]
 sij=[[None,None],[None,None]]
 for i,j in [(0,0),(0,1),(1,1)]:
  rij=((conj(qi[i])*qi[j]+conj(q)*qij[i][j]).real-ri[i]*ri[j])/r
  if ACTIVE:cij=ilim*(qij[i][j]/r-qi[i]*ri[j]/(r*r)-qi[j]*ri[i]/(r*r)-q*rij/(r*r)+2*q*ri[i]*ri[j]/(r*r*r))
  else:cij=qij[i][j]
  val=net['vgE']*Z*cij*conj(ig)+net['vgE']*Z*ci[i]*conj(net['igE']*Z*ci[j])+net['vgE']*Z*ci[j]*conj(net['igE']*Z*ci[i])+vg*conj(net['igE']*Z*cij)
  sij[i][j]=val;sij[j][i]=val
 ep=-I('.1')*si[0].imag/(1+I('.1')*si[1].imag)
 epp=-I('.1')*(sij[0][0].imag+2*sij[0][1].imag*ep+sij[1][1].imag*ep*ep)/(1+I('.1')*si[1].imag)
 return sij[0][0].real+2*sij[0][1].real*ep+sij[1][1].real*ep*ep+si[1].real*epp

def qvalues(g,nu,er):
 P,Pg,r,u,E=er;z1=(g-cen[0])/sc[0];z2=(nu-cen[1])/sc[1];f1=om*(nu-1);f2=((1-P)/nu-D*(nu-1))/mass;a=-Pg/(mass*nu);b=(-(1-P)/(nu*nu)-D)/mass
 M,MD=matpoly(z1,z2,f1/sc[0],f2/sc[1]);m11=M[0][0];m12=M[0][1];m22=M[1][1]
 return M,[2*a*m12+MD[0][0]+2*lam*m11,om*m11+b*m12+a*m22+MD[0][1]+2*lam*m12,2*om*m12+2*b*m22+MD[1][1]+2*lam*m22]

def taylorQ(g,nu,er):
 P,Pg,r,u,E=er;Pgg=second(g,E);z1=(g-cen[0])/sc[0];z2=(nu-cen[1])/sc[1]
 f1=om*(nu-1);f2=((1-P)/nu-D*(nu-1))/mass;a=-Pg/(mass*nu);b=(-(1-P)/(nu*nu)-D)/mass
 ag=-Pgg/(mass*nu);bg=Pg/(mass*nu*nu);av=Pg/(mass*nu*nu);bv=2*(1-P)/(mass*nu*nu*nu)
 M,_=matpoly(z1,z2,I(0),I(0))
 mg=[[(C[1][i][j]+2*z1*C[3][i][j]+z2*C[4][i][j])/sc[0] for j in range(2)] for i in range(2)]
 mv=[[(C[2][i][j]+z1*C[4][i][j]+2*z2*C[5][i][j])/sc[1] for j in range(2)] for i in range(2)]
 mgg=[[2*C[3][i][j]/(sc[0]*sc[0]) for j in range(2)] for i in range(2)]
 mgv=[[C[4][i][j]/(sc[0]*sc[1]) for j in range(2)] for i in range(2)]
 mvv=[[2*C[5][i][j]/(sc[1]*sc[1]) for j in range(2)] for i in range(2)]
 mdg=[[mgg[i][j]*f1+mgv[i][j]*f2+mv[i][j]*a for j in range(2)] for i in range(2)]
 mdv=[[mgv[i][j]*f1+mvv[i][j]*f2+mg[i][j]*om+mv[i][j]*b for j in range(2)] for i in range(2)]
 dg=[2*ag*M[0][1]+2*a*mg[0][1]+mdg[0][0]+2*lam*mg[0][0],om*mg[0][0]+bg*M[0][1]+b*mg[0][1]+ag*M[1][1]+a*mg[1][1]+mdg[0][1]+2*lam*mg[0][1],2*om*mg[0][1]+2*bg*M[1][1]+2*b*mg[1][1]+mdg[1][1]+2*lam*mg[1][1]]
 dv=[2*av*M[0][1]+2*a*mv[0][1]+mdv[0][0]+2*lam*mv[0][0],om*mv[0][0]+bv*M[0][1]+b*mv[0][1]+av*M[1][1]+a*mv[1][1]+mdv[0][1]+2*lam*mv[0][1],2*om*mv[0][1]+2*bv*M[1][1]+2*b*mv[1][1]+mdv[1][1]+2*lam*mv[1][1]]
 gc=iv.mpf((g.a+g.b)/2);vc=iv.mpf((nu.a+nu.b)/2);ec,reason=enclosure(gc)
 if ec is None:return None
 qc=qvalues(gc,vc,ec)[1]
 return [qc[i]+dg[i]*(g-gc)+dv[i]*(nu-vc) for i in range(3)]

def mode_ok(r):return pos(r-ilim) if ACTIVE else neg(r-ilim)
def mode_violated(r):return neg(r-ilim) if ACTIVE else pos(r-ilim)

def check(g,nu):
 """Returns (ok, reason, info). info holds certified bounds for a passed box."""
 er,reason=enclosure(g)
 if er is None:return False,reason,None
 P,Pg,r,u,E=er
 if not mode_ok(r):return False,'mode',None
 if not pos(vlim*vlim-abs2(u)):return False,'voltage',None
 M,q=qvalues(g,nu,er);m11=M[0][0];m12=M[0][1];m22=M[1][1]
 # Margin levels: 2 = epsilon margin holds, 1 = only the strict condition holds.
 if pd(m11-eps,m12,m22-eps):ml=2
 elif pd(m11,m12,m22):ml=1
 else:return False,'M lower',None
 if pd(up-m11-eps,-m12,up-m22-eps):mu=2
 elif pd(up-m11,-m12,up-m22):mu=1
 else:return False,'M upper',None
 qforms=[q];qm=0
 if pd(-q[0]-eps,-q[1],-q[2]-eps):qm=2
 elif pd(-q[0],-q[1],-q[2]):qm=1
 if qm<2:
  qt=taylorQ(g,nu,er)
  if qt is not None:
   qforms.append(qt)
   if pd(-qt[0]-eps,-qt[1],-qt[2]-eps):qm=2
   elif qm==0 and pd(-qt[0],-qt[1],-qt[2]):qm=1
 if qm==0:return False,'Q',None
 e=float(eps.a)
 lmin=dn(eigs(m11,m12,m22)[0].a);umin=dn(eigs(up-m11,-m12,up-m22)[0].a)
 qmax=min(upf(eigs(*f)[1].b) for f in qforms)
 info=dict(minM=max(lmin,e if ml==2 else 0.0),min10IminusM=max(umin,e if mu==2 else 0.0),maxQ=min(qmax,-e if qm==2 else 0.0),
  eps_margin=bool(ml==2 and mu==2 and qm==2),
  current_margin=dn(((r-ilim) if ACTIVE else (ilim-r)).a),voltage_margin=dn((vlim-iv.sqrt(abs2(u))).a))
 return True,'pass',info

def point_counterexample(g,nu):
 """Evaluate at the box centre with degenerate intervals. Returns a reason if a
 condition is CERTIFIABLY violated there (genuine counterexample), else None."""
 gc=iv.mpf((g.a+g.b)/2);vc=iv.mpf((nu.a+nu.b)/2);er,reason=enclosure(gc)
 if er is None:return None
 P,Pg,r,u,E=er
 if mode_violated(r):return 'mode'
 if pos(abs2(u)-vlim*vlim):return 'voltage'
 M,q=qvalues(gc,vc,er)
 if neg(eigs(M[0][0],M[0][1],M[1][1])[0]):return 'M lower'
 if neg(eigs(up-M[0][0],-M[0][1],up-M[1][1])[0]):return 'M upper'
 if pos(eigs(*q)[1]):return 'Q'
 return None

start=time.time();stack=[(bounds[0][0],bounds[1][0],bounds[0][1],bounds[1][1],0)]
limit=int(sys.argv[3]) if len(sys.argv)>3 else 400000;MAXDEP=40
# Boxes that satisfy only the strict (zero-margin) conditions are split further up to
# SOFTDEP to try to certify the epsilon margin; beyond that the strict pass is accepted.
SOFTDEP=22
passed=0;tested=0;failures={};depthmax=0;unresolved=[];counter=[]
agg=dict(minM=None,min10IminusM=None,maxQ=None);epsall=True;minc=1e99;minv=1e99
while stack and tested<limit:
 ga,gb,va,vb,dep=stack.pop();g=iv.mpf([ga.a,gb.b]);nu=iv.mpf([va.a,vb.b]);tested+=1
 ok,reason,info=check(g,nu)
 if ok and not info['eps_margin'] and dep<SOFTDEP:ok,reason=False,'eps margin only'
 if ok:
  passed+=1;epsall&=info['eps_margin'];minc=min(minc,info['current_margin']);minv=min(minv,info['voltage_margin'])
  for k,fn in (('minM',min),('min10IminusM',min),('maxQ',max)):agg[k]=info[k] if agg[k] is None else fn(agg[k],info[k])
 else:
  failures[reason]=failures.get(reason,0)+1
  box=[float(g.a)*180/float(iv.pi.a),float(g.b)*180/float(iv.pi.a),float(nu.a),float(nu.b)]
  if dep>=6 and reason!='eps margin only':
   ce=point_counterexample(g,nu)
   if ce is not None:counter.append(dict(reason=ce,box_deg_nu=box,depth=dep));continue
  if dep>=MAXDEP:unresolved.append(dict(reason=reason,box_deg_nu=box,depth=dep));continue
  wg=float(((gb-ga)/sc[0]).b);wv=float(((vb-va)/sc[1]).b)
  if wg>=wv:mid=(ga+gb)/2;stack.extend([(ga,mid,va,vb,dep+1),(mid,gb,va,vb,dep+1)])
  else:mid=(va+vb)/2;stack.extend([(ga,gb,va,mid,dep+1),(ga,gb,mid,vb,dep+1)])
  depthmax=max(depthmax,dep+1)
 if tested%2000==0:print(json.dumps({'tested':tested,'passed':passed,'pending':len(stack),'depth':depthmax,'counterexamples':len(counter),'seconds':round(time.time()-start,1),'last':reason}),flush=True)
complete=not stack and not unresolved and not counter
status='VERIFIED' if complete else ('FAILED' if counter else 'INCONCLUSIVE')
f=lambda x:None if x is None else float(x)
result={'branch':branch,'status':status,
 'gamma_bounds_deg':[float(bounds[0][0].a)*180/float(iv.pi.a),float(bounds[1][0].a)*180/float(iv.pi.a)],'nu_bounds':[float(bounds[0][1].a),float(bounds[1][1].a)],
 'lambda':data['lambda'],'epsilon':data['epsilon'],'upper':data['upper'],
 'certified_min_eig_M_lower_bound':f(agg['minM']),'certified_min_eig_10I_minus_M_lower_bound':f(agg['min10IminusM']),
 'certified_max_eig_contraction_matrix_upper_bound':f(agg['maxQ']),
 'epsilon_margin_holds_in_all_boxes':bool(epsall and complete),
 'verified_boxes':passed,'unresolved_boxes':len(unresolved)+len(stack),'counterexample_boxes':len(counter),'tested_boxes':tested,'max_depth':depthmax,
 'budget_exhausted':bool(stack),'failures_before_subdivision':failures,
 'min_certified_current_margin_pu':minc if passed else None,'min_certified_voltage_margin_pu':minv if passed else None,
 'mode_or_voltage_violations':[c for c in counter if c['reason'] in ('mode','voltage')],
 'counterexamples':counter[:50],'unresolved':unresolved[:50],
 'precision_decimal_digits':iv.dps,'elapsed_s':time.time()-start,
 'notes':'Directed-rounding interval arithmetic (mpmath iv, 30 digits) on fixed decimal coefficients from metric.json; omega=100*pi exact; Mdot = dM/dz * f/scale with the nonlinear vector field of the branch.'}
(base/'interval_verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=1),flush=True)
