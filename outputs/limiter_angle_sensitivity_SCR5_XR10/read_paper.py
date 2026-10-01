from pathlib import Path
import os
from pypdf import PdfReader
repo=Path(__file__).resolve().parents[2]
p=Path(os.environ.get('GFM_REFERENCE_PDF','')).expanduser()
if not p.is_file(): raise SystemExit('Set GFM_REFERENCE_PDF to a local copy of the reference PDF.')
d=PdfReader(p);o=repo/'outputs'/'limiter_angle_sensitivity_SCR5_XR10'
with (o/'paper_extracted.txt').open('w',encoding='utf-8') as f:
 for k,page in enumerate(d.pages):f.write('\n\n===== PDF PAGE '+str(k+1)+' =====\n'+page.extract_text(extraction_mode='layout'))
print(len(d.pages))
