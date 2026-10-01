from zipfile import ZipFile
from pathlib import Path
import xml.etree.ElementTree as ET
repo=Path(__file__).resolve().parents[2]
lib=repo/'model'/'tdmsquared-master'/'tdmsquared-master'/'Library'
p=next(lib.rglob('Subsys_Conv_Electrical_Part_SPS.slx'))
with ZipFile(p) as z:
 for name in z.namelist():
  if name.startswith('simulink/systems/') and name.endswith('.xml'):
   root=ET.fromstring(z.read(name))
   blocks=root.findall('Block')
   if any(b.get('Name')=='LC_harmonic_filter' for b in blocks):
    print('SYSTEM',name)
    for b in blocks:print(b.attrib,[(v.get('Name'),v.text) for v in b.findall('P') if v.get('Name') in ['Side','Port','SourceBlock']])
    for l in root.findall('Line'):print(ET.tostring(l,encoding='unicode'))
