from pathlib import Path
from urllib.request import Request,urlopen
from concurrent.futures import ThreadPoolExecutor
import hashlib,json
root=Path('work/S_sources');root.mkdir(exist_ok=True)
items={
 'Liu2018_S1_table.docx':'https://journals.plos.org/plosmedicine/article/file?id=10.1371/journal.pmed.1002718.s005&type=supplementary',
 'NHANES_2001_2002_MORT_2019_PUBLIC.dat':'https://ftp.cdc.gov/pub/health_statistics/nchs/datalinkage/linked_mortality/NHANES_2001_2002_MORT_2019_PUBLIC.dat'}
def fetch(pair):
 n,u=pair;p=root/n
 try:
  if not p.exists():
   with urlopen(Request(u,headers={'User-Agent':'BIOAGE-REL-01 authorized writing supplements'}),timeout=90) as r:p.write_bytes(r.read())
  b=p.read_bytes();result=dict(file=n,url=u,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),status='OK');print(result,flush=True);return result
 except Exception as e:result=dict(file=n,url=u,status='NOT_RUN',reason=str(e));print(result,flush=True);return result
with ThreadPoolExecutor(max_workers=2) as ex:rows=list(ex.map(fetch,items.items()))
(root/'downloads.json').write_text(json.dumps(rows,indent=2)+'\n')
