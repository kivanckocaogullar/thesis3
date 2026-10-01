from pathlib import Path
import json
s=Path('outputs/validate_switch_recovery_connection.py').read_text();exec(s.split('gains=[]')[0])
v=json.loads((out/'connection_validation.json').read_text());h0=I('.0063');wl=I('0.0056626116124433')
def kinetic(v):w=v-1;return I('6.8')*100*iv.pi*(w*w/2+w*w*w/3)
assert pos(wl+kinetic(I('.9992'))-h0) and pos(wl+kinetic(I('1.0008'))-h0)
# At any switching event in the energy sublevel, nu lies in this smaller strip.
src=s[s.index('def transition('):s.index('gains=[]')].replace("lo=I('0.99652713910517932');hi=I('1.0025')","lo=I('.9992');hi=I('1.0008')")
exec(src)
gains=[]
for rev in [False,True]:
 for k in [1,2,3,4,6,8,12,16,32,64,128,256]:
  if transition(k,rev):gains.append(k);break
assert len(gains)==2
# Explicit checks for the small omitted integration intervals in W().
swp=enc(N,sw)[0]-1;eqp=enc(N,eq)[0]-1
assert max(abs(swp.a),abs(swp.b))*(sw.b-sw.a)<I('1e-10')
assert max(abs(eqp.a),abs(eqp.b))*(eq.b-eq.a)<I('1e-10')
result={'nu_switch_bounds':[.9992,1.0008],'gain_active_to_inactive_quadratic':gains[0],'gain_inactive_to_active_quadratic':gains[1],'scope':'energy recovery sublevel H<=.0063 only','small_integral_corrections_verified':True}
(out/'recovery_transition_gain.json').write_text(json.dumps(result,indent=2));print(result)
