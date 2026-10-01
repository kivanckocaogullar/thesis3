from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET
repo=Path(__file__).resolve().parents[2]
lib=repo/'model'/'tdmsquared-master'/'tdmsquared-master'/'Library'
names=['Subsys_GFormConv_InnerCtrl_cascadedPI.slx','dq_Limiter.slx','PI_freeze.slx','Subsys_GFormConv_Control.slx','Subsys_Conv_Ctrl_measurement_transformations.slx']
out=repo/'outputs'/'limiter_audit_SCR5_XR10'/'xml_control_trace.txt'
with out.open('w',encoding='utf-8') as f:
 for name in names:
  path=next(lib.rglob(name));f.write('\nFILE '+path.relative_to(repo).as_posix()+'\n')
  with ZipFile(path) as z:
   for entry in z.namelist():
    if not(entry.startswith('simulink/systems/') and entry.endswith('.xml')):continue
    root=ET.fromstring(z.read(entry));f.write('\nSYSTEM '+entry+'\n')
    for b in root.findall('Block'):
     f.write('BLOCK '+repr(b.attrib)+' '+repr({p.attrib.get('Name'):p.text for p in b.findall('P')})+'\n')
     for s in b.findall('System'):f.write(' SUBSYSTEM '+repr(s.attrib)+'\n')
    for line in root.findall('Line'):f.write('LINE '+ET.tostring(line,encoding='unicode')+'\n')
    for s in root.findall('.//Initialization'):f.write('INIT '+str(s.text)+'\n')
print(out)
