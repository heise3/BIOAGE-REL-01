from fetch_sources import fetch, ROOT
from concurrent.futures import ThreadPoolExecutor
import json
base='https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1988/DataFiles/'
items={n+'.'+ext:base+n+'.'+ext for n in ['SSNH3CYS','SSNH3BTP'] for ext in ['xpt','htm']}
items['SSNH3CYS.xpt']='https://wwwn.cdc.gov/nchs/data/nhanes3/37a/ssnh3cys.xpt'
items['SSNH3BTP.xpt']='https://wwwn.cdc.gov/nchs/data/nhanes3/36a/SSNH3BTP.xpt'
items['sscystat.xpt']='https://wwwn.cdc.gov/nchs/data/nhanes3/27a/sscystat.xpt'
with ThreadPoolExecutor(max_workers=4) as ex: out=list(ex.map(fetch,items.items()))
(ROOT/'downloads_kidney.json').write_text(json.dumps(out,indent=2)+'\n')
