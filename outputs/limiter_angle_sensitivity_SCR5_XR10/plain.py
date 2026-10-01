from pypdf import PdfReader
from pathlib import Path
import os
repo=Path(__file__).resolve().parents[2]
p=Path(os.environ.get('GFM_REFERENCE_PDF','')).expanduser()
if not p.is_file(): raise SystemExit('Set GFM_REFERENCE_PDF to a local copy of the reference PDF.')
r=PdfReader(p)
(repo/'outputs/limiter_angle_sensitivity_SCR5_XR10/paper_plain.txt').write_text('\n'.join('\nPAGE '+str(i+1)+'\n'+pg.extract_text() for i,pg in enumerate(r.pages)),encoding='utf-8')
