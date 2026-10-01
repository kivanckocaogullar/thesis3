"""Directed interval validation for fixed published decimal network/metric.
Not SOS; actual sin/cos, square roots, and an implicit Q-V root enclosure.
"""
import sys,json,time,heapq
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent.parent/'radial_interval_validation_deps'))
from mpmath import iv
iv.dps=30
base=Path(sys.argv[1]);data=json.loads((base/'metric.json').read_text(),parse_float=str)
I=lambda x:iv.mpf(str(x))
lo=lambda x:x._mpi_[0]
# Use interval endpoints directly, never round inward through float.
def lower(x):return x.a
def upper(x):return x.b
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
 ic=q
 ig=rg;iE=rE
 vc=T+Z*ic;vg=net['vgE']*vc+net['vgU']*U;ii=net['igE']*vc+net['igU']*U
 S=vg*conj(ii);Sg=net['vgE']*Z*ig*conj(ii)+vg*conj(net['igE']*Z*ig);SE=net['vgE']*Z*iE*conj(ii)+vg*conj(net['igE']*Z*iE)
 return S,Sg,SE,r,vc+(net['Zconverter']+AD)*ic

_enclosure_cache={}
def enclosure(g):
 key=g._mpi_
 if key not in _enclosure_cache:_enclosure_cache[key]=_enclosure(g)
 return _enclosure_cache[key]
def _enclosure(g):
 E=iv.mpf(['.98','1.05']);S0=electrical(g,I('.98'))[0];S1=electrical(g,I('1.05'))[0]
 if not (neg(I('.98')+I('.1')*S0.imag-1) and pos(I('1.05')+I('.1')*S1.imag-1)):return None,'root bracket'
 hE=1+I('.1')*electrical(g,E)[2].imag
 if not pos(hE):return None,'root uniqueness'
 # Uniform bracket + positive derivative establish one root for every g.
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

def pd(a,b,c):return pos(a) and pos(a*c-b*b)
_second_cache={}
def second(g,E):
 key=(g._mpi_,E._mpi_)
 if key not in _second_cache:_second_cache[key]=_second(g,E)
 return _second_cache[key]
def _second(g,E):
 rot=iv.mpc(iv.cos(g),iv.sin(g))*rot0;v=E*rot;q=(v-T)/Z;r=iv.sqrt(abs2(q));qi=[1j*v/Z,rot/Z];qij=[[-v/Z,1j*rot/Z],[1j*rot/Z,iv.mpc(0)]]
 ri=[(conj(q)*u).real/r for u in qi];ci=qi;c=q
 vc=T+Z*c;vg=net['vgE']*vc+net['vgU']*U;ig=net['igE']*vc+net['igU']*U
 si=[net['vgE']*Z*ci[i]*conj(ig)+vg*conj(net['igE']*Z*ci[i]) for i in range(2)]
 sij=[[None,None],[None,None]]
 for i,j in [(0,0),(0,1),(1,1)]:
  rij=((conj(qi[i])*qi[j]+conj(q)*qij[i][j]).real-ri[i]*ri[j])/r
  cij=qij[i][j]
  val=net['vgE']*Z*cij*conj(ig)+net['vgE']*Z*ci[i]*conj(net['igE']*Z*ci[j])+net['vgE']*Z*ci[j]*conj(net['igE']*Z*ci[i])+vg*conj(net['igE']*Z*cij)
  sij[i][j]=val;sij[j][i]=val
 ep=-I('.1')*si[0].imag/(1+I('.1')*si[1].imag)
 epp=-I('.1')*(sij[0][0].imag+2*sij[0][1].imag*ep+sij[1][1].imag*ep*ep)/(1+I('.1')*si[1].imag)
 return sij[0][0].real+2*sij[0][1].real*ep+sij[1][1].real*ep*ep+si[1].real*epp

def qvalues(g,nu,er):
 P,Pg,r,u,E=er;z1=(g-cen[0])/sc[0];z2=(nu-cen[1])/sc[1];f1=om*(nu-1);f2=((1-P)/nu-D*(nu-1))/mass;a=-Pg/(mass*nu);b=(-(1-P)/(nu*nu)-D)/mass
 M,MD=matpoly(z1,z2,f1/sc[0],f2/sc[1]);m11=M[0][0];m12=M[0][1];m22=M[1][1]
 return [2*a*m12+MD[0][0]+2*lam*m11,om*m11+b*m12+a*m22+MD[0][1]+2*lam*m12,2*om*m12+2*b*m22+MD[1][1]+2*lam*m22]

def taylorQ(g,nu,er):
 P,Pg,r,u,E=er;Pgg=second(g,E);z1=(g-cen[0])/sc[0];z2=(nu-cen[1])/sc[1]
 f1=om*(nu-1);f2=((1-P)/nu-D*(nu-1))/mass;a=-Pg/(mass*nu);b=(-(1-P)/(nu*nu)-D)/mass
 ag=-Pgg/(mass*nu);bg=Pg/(mass*nu*nu);av=Pg/(mass*nu*nu);bv=2*(1-P)/(mass*nu*nu*nu)
 M,_=matpoly(z1,z2,I(0),I(0))
 mg=[[ (C[1][i][j]+2*z1*C[3][i][j]+z2*C[4][i][j])/sc[0] for j in range(2)] for i in range(2)]
 mv=[[ (C[2][i][j]+z1*C[4][i][j]+2*z2*C[5][i][j])/sc[1] for j in range(2)] for i in range(2)]
 mgg=[[2*C[3][i][j]/(sc[0]*sc[0]) for j in range(2)] for i in range(2)]
 mgv=[[C[4][i][j]/(sc[0]*sc[1]) for j in range(2)] for i in range(2)]
 mvv=[[2*C[5][i][j]/(sc[1]*sc[1]) for j in range(2)] for i in range(2)]
 mdg=[[mgg[i][j]*f1+mgv[i][j]*f2+mv[i][j]*a for j in range(2)] for i in range(2)]
 mdv=[[mgv[i][j]*f1+mvv[i][j]*f2+mg[i][j]*om+mv[i][j]*b for j in range(2)] for i in range(2)]
 dg=[2*ag*M[0][1]+2*a*mg[0][1]+mdg[0][0]+2*lam*mg[0][0],om*mg[0][0]+bg*M[0][1]+b*mg[0][1]+ag*M[1][1]+a*mg[1][1]+mdg[0][1]+2*lam*mg[0][1],2*om*mg[0][1]+2*bg*M[1][1]+2*b*mg[1][1]+mdg[1][1]+2*lam*mg[1][1]]
 dv=[2*av*M[0][1]+2*a*mv[0][1]+mdv[0][0]+2*lam*mv[0][0],om*mv[0][0]+bv*M[0][1]+b*mv[0][1]+av*M[1][1]+a*mv[1][1]+mdv[0][1]+2*lam*mv[0][1],2*om*mv[0][1]+2*bv*M[1][1]+2*b*mv[1][1]+mdv[1][1]+2*lam*mv[1][1]]
 gc=(g.a+g.b)/2;vc=(nu.a+nu.b)/2;ec,reason=enclosure(gc)
 if ec is None:return None
 qc=qvalues(gc,vc,ec)
 return [qc[i]+dg[i]*(g-gc)+dv[i]*(nu-vc) for i in range(3)]

def check(g,nu):
 er,reason=enclosure(g)
 if er is None:return False,reason,None
 P,Pg,r,u,E=er
 # Branch restriction is imposed separately at the validated switch root.
 if not pos(vlim*vlim-abs2(u)):return False,'voltage',None
 z1=(g-cen[0])/sc[0];z2=(nu-cen[1])/sc[1]
 f1=om*(nu-1);f2=((1-P)/nu-D*(nu-1))/mass;a=-Pg/(mass*nu);b=(-(1-P)/(nu*nu)-D)/mass
 M,MD=matpoly(z1,z2,f1/sc[0],f2/sc[1]);m11=M[0][0];m12=M[0][1];m22=M[1][1]
 if not pd(m11-eps,m12,m22-eps):return False,'M lower',None
 if not pd(up-m11,-m12,up-m22):return False,'M upper',None
 q11=2*a*m12+MD[0][0]+2*lam*m11
 q12=om*m11+b*m12+a*m22+MD[0][1]+2*lam*m12
 q22=2*om*m12+2*b*m22+MD[1][1]+2*lam*m22
 if not pd(-q11-eps,-q12,-q22-eps):
  qt=taylorQ(g,nu,er)
  if qt is None or not pd(-qt[0]-eps,-qt[1],-qt[2]-eps):return False,'Q',None
 return True,'pass',(float((ilim.a-r.b).a),float((vlim.a-iv.sqrt(abs2(u)).b).a))

start=time.time();stack=[(bounds[0][0],bounds[1][0],bounds[0][1],bounds[1][1],0)];passed=0;tested=0;failures={};depthmax=0;minbranch=1e99;minvoltage=1e99
limit=int(sys.argv[2]) if len(sys.argv)>2 else 20000
while stack and tested<limit:
 ga,gb,va,vb,dep=stack.pop();g=iv.mpf([ga.a,gb.b]);nu=iv.mpf([va.a,vb.b]);tested+=1
 ok,reason,marg=check(g,nu)
 if ok:passed+=1;minbranch=min(minbranch,marg[0]);minvoltage=min(minvoltage,marg[1])
 else:
  failures[reason]=failures.get(reason,0)+1
  if dep>=28:stack.append((ga,gb,va,vb,dep));break
  wg=float(((gb-ga)/sc[0]).b);wv=float(((vb-va)/sc[1]).b)
  if wg>=wv:
   mid=(ga+gb)/2;stack.extend([(ga,mid,va,vb,dep+1),(mid,gb,va,vb,dep+1)])
  else:
   mid=(va+vb)/2;stack.extend([(ga,gb,va,mid,dep+1),(ga,gb,mid,vb,dep+1)])
  depthmax=max(depthmax,dep+1)
 if tested%500==0:print(json.dumps({'tested':tested,'passed':passed,'pending':len(stack),'depth':depthmax,'seconds':round(time.time()-start,1),'last':reason}),flush=True)
result={'complete':not stack,'tested':tested,'passed_boxes':passed,'pending_boxes':len(stack),'max_depth':depthmax,'precision_decimal_digits':iv.dps,'failures_to_subdivide':failures,'elapsed_s':time.time()-start,'min_branch_margin_passed':minbranch,'min_voltage_margin_passed':minvoltage,'epsilon':data['epsilon'],'lambda':data['lambda'],'status':'validated full rectangle for fixed decimal coefficients' if not stack else 'incomplete interval attempt; no full certificate'}
(base/'interval_validation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)






