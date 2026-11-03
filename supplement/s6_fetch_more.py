from pathlib import Path
from urllib.request import Request,urlopen
from concurrent.futures import ThreadPoolExecutor
import hashlib,json
O=Path('work/S6_sources')
items={
 'nhsr021_stacks.pdf':'https://stacks.cdc.gov/view/cdc/5756/cdc_5756_DS1.pdf',
 'PH_B.htm':'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/PH_B.htm',
 'PH_B.xpt':'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/PH_B.xpt',
 'lab_variablelist.html':'https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Laboratory&Cycle=2001-2002',
 'labmethods.html':'https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/labmethods.aspx?BeginYear=2001',
 'seconddaylabdoc_b.pdf':'https://wwwn.cdc.gov/nchs/data/nhanes/nhanes_01_02/seconddaylabdoc_b.pdf',
}
def fetch(it):
 name,url=it;p=O/name
 try:
  if not p.exists():
   with urlopen(Request(url,headers={'User-Agent':'Mozilla/5.0'}),timeout=90) as f:b=f.read()
   p.write_bytes(b)
  b=p.read_bytes()
  return dict(file=name,url=url,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),status='OK')
 except Exception as e:return dict(file=name,url=url,status='NOT_RUN',reason=str(e))
with ThreadPoolExecutor(max_workers=4) as ex:res=list(ex.map(fetch,items.items()))
(O/'downloads_more.json').write_text(json.dumps(res,indent=2))
for r in res:print(r['file'],r['status'],r.get('bytes',r.get('reason')))
