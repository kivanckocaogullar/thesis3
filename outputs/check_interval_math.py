import sys,json
from pathlib import Path
sys.argv=['test', 'outputs/radial_quadratic_robust_20261001_161903']
src=Path('outputs/validate_radial_polynomial_intervals_taylor.py').read_text(encoding='utf-8-sig')
exec(src.split('start=time.time();')[0])
worst=0.;checked=0
for i in range(11):
 gc=cen[0]+sc[0]*I(str(-.9+.18*i));vc=cen[1]+sc[1]*I(str(-.8+.16*i));h=I('0.000001')
 e,_=enclosure(gc);ep,_=enclosure(gc+h);em,_=enclosure(gc-h)
 analytic=second(gc,e[4]);fd=(ep[1]-em[1])/(2*h)
 difference=analytic-fd;worst=max(worst,abs(float(difference.a)),abs(float(difference.b)))
 g=iv.mpf([(gc-I('.0001')).a,(gc+I('.0001')).b]);v=iv.mpf([(vc-I('.000001')).a,(vc+I('.000001')).b]);er,_=enclosure(g);Q=taylorQ(g,v,er)
 for gg in [g.a,g.b]:
  ee,_=enclosure(gg)
  for vv in [v.a,v.b]:
   qq=qvalues(gg,vv,ee)
   for k in range(3):assert qq[k].a>=Q[k].a and qq[k].b<=Q[k].b
   checked+=1
assert worst<1e-7
result={'secondDerivativeFiniteDifferenceMaxError':worst,'TaylorCornerChecks':checked,'passed':True}
(base/'interval_math_checks.json').write_text(json.dumps(result,indent=2));print(result)
