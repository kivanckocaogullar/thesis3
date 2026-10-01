import json
from pathlib import Path
p=Path('outputs/radial_switch_connection_20261001');d=json.loads(Path('outputs/radial_inactive_metric_20261001_170815/metric.json').read_text());d['bounds'][1]=[0.42539661083071958,1.0025];d['upper']=15
o=p/'inactive_recovery_extension';o.mkdir();(o/'metric.json').write_text(json.dumps(d,indent=2))
# Verify that prior extension serialization did not alter coefficients.
e=json.loads((p/'inactive_upper12/metric.json').read_text(encoding='utf-8-sig'));assert e['coefficients']==d['coefficients']
