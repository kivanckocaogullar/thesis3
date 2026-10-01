import sys,json
from pathlib import Path
sys.argv=['check','outputs/radial_inactive_metric_20261001_170815']
exec(Path('outputs/validate_radial_inactive_intervals.py').read_text(encoding='utf-8-sig').split('start=time.time();')[0])
gc=I('0.36012957762291908');rg=I('1e-12');gl=gc-rg;gh=gc+rg
el,_=enclosure(gl);eh,_=enclosure(gh)
assert neg(el[0]-1) and pos(eh[0]-1)
# P is frozen at a decimal numerical center; Lyapunov function is centered
# on the unique exact equilibrium enclosed by [gl,gh].
MM,_=matpoly((gc-cen[0])/sc[0],(I(1)-cen[1])/sc[1],I(0),I(0))
P=[[I(format(float(MM[i][j].mid),'.17g')) for j in range(2)] for i in range(2)]
p11,p12,p22=P[0][0],P[0][1],P[1][1]
assert pd(p11,p12,p22)
for rad in ['.001','.0005','.0001','.00005']:
 hg=I(rad);hv=I('.00005');g=iv.mpf([(gl-hg).a,(gh+hg).b]);v=iv.mpf([(1-hv).a,(1+hv).b]);er,reason=enclosure(g)
 if er is None:continue
 power,pg,r,u,E=er;a=-pg/(mass*v);b=(-(1-power)/(v*v)-D)/mass
 q11=2*a*p12;q12=om*p11+b*p12+a*p22;q22=2*om*p12+2*b*p22
 if pd(-q11,-q12,-q22) and pos(ilim-r) and pos(vlim*vlim-abs2(u)):
  det=p11*p22-p12*p12;inv11=p22/det;inv22=p11/det
  rho=I(format(.5*min(float((hg*hg/inv11).a),float((hv*hv/inv22).a)),'.15g'))
  assert pos(hg*hg-rho*inv11) and pos(hv*hv-rho*inv22)
  assert pos(pg)
  result={'validated':True,'equilibriumGammaBracket_rad':[float(gl.a),float(gh.b)],'P':[[float(t.a) for t in row] for row in P],'halfwidth_gamma_rad':rad,'halfwidth_nu':'.00005','rho':str(float(rho.a)),'definition':'(x-xeq)^T P (x-xeq) <= rho; xeq is exact unique implicit equilibrium, nu=1','proof':'J^T P+P J negative throughout containing box; mean-value integral establishes Vdot<0; ellipsoid contained strictly in box.'}
  break
else:result={'validated':False}
(base/'invariant_subset.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
