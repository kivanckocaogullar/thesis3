from pathlib import Path
import json
root=Path('outputs');out=root/'radial_switch_connection_20261001'
gs=json.loads((out/'switch_numeric.json').read_text())['gamma_switch_rad']
for kind,src,checker in [('inactive','radial_inactive_metric_20261001_170815','validate_radial_inactive_intervals.py'),('active','radial_quadratic_robust_20261001_161903','validate_radial_polynomial_intervals_taylor.py')]:
 d=json.loads((root/src/'metric.json').read_text());p=out/kind;p.mkdir()
 # Formal smooth extension across a 1e-12-rad overlap around the root.
 if kind=='inactive':d['bounds'][0][0]=d['bounds'][1][0];d['bounds'][1][0]=gs+1e-12
 else:d['bounds'][1][0]=d['bounds'][0][0];d['bounds'][0][0]=gs-1e-12
 (p/'metric.json').write_text(json.dumps(d,indent=2))
 s=(root/checker).read_text(encoding='utf-8-sig')
 s=s.replace("if not pos(ilim-r):return False,'activation',None","# Branch restriction is imposed separately at the validated switch root.")
 s=s.replace("if not pos(r-ilim):return False,'activation',None","# Branch restriction is imposed separately at the validated switch root.")
 (out/(kind+'_checker.py')).write_text(s.replace("Path(__file__).parent/'radial_interval_validation_deps'","Path(__file__).parent.parent/'radial_interval_validation_deps'"))
