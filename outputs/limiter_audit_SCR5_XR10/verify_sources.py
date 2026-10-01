from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET
import hashlib
repo=Path(__file__).resolve().parents[2]
p=repo/'model'/'GFMVSG1.slx'
with ZipFile(p) as z:
 for fn in z.namelist():
  if fn.startswith('simulink/systems/') and fn.endswith('.xml'):
   root=ET.fromstring(z.read(fn))
   for b in root.findall('Block'):
    if b.attrib.get('SID') in ['40','41']:
     print(b.attrib)
     for x in b.findall('.//P'):
      if x.attrib.get('Name') in ['PositiveSequence','Resistance','Inductance','BranchType','VariationEntity']:print(x.attrib.get('Name'),x.text)
     for x in b.findall('.//MaskParameter'):
      if x.attrib.get('Name') in ['PositiveSequence','Resistance','Inductance','BranchType','VariationEntity']:print(ET.tostring(x,encoding='unicode'))
paths=[p,Path('outputs/baseline_gfm_contraction.m'),Path('outputs/gfm_current_limiter_analysis.m'),Path('outputs/baseline_gfm_contraction_SCR5_XR10.m')]
Path('outputs/limiter_audit_SCR5_XR10/source_hashes.txt').write_text('\n'.join(str(f)+' '+hashlib.sha256(f.read_bytes()).hexdigest() for f in paths),encoding='utf-8')
f=Path('outputs/limiter_audit_SCR5_XR10/audit_report.md');s=f.read_text(encoding='utf-8-sig');s=s.replace('Required |I_conv| (pu) | 1.2-|I_conv| (pu) | Required pre-active-damping |u| (pu)', 'Required current magnitude (pu) | Current margin (pu) | Pre-active-damping voltage magnitude (pu)');f.write_text(s,encoding='utf-8')
