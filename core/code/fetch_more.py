from fetch_sources import fetch, ROOT
from concurrent.futures import ThreadPoolExecutor
import json

base='https://ftp.cdc.gov/pub/health_statistics/nchs/datalinkage/linked_mortality/'
items={n:base+n for n in ['NHANES_III_MORT_2019_PUBLIC.dat','R_ReadInProgramAllSurveys.R']}
pin='https://raw.githubusercontent.com/dayoonkwon/BioAge/b1f9fc02f086cd4aa74185f2335ab1366082e7fe/'
items['NHANES3.rda']=pin+'data/NHANES3.rda'
items['NHANES3_HDTrain.rda']=pin+'data/NHANES3_HDTrain.rda'
items['bioage_nhanes_all.R']=pin+'data-raw/nhanes_all.R'
for n in ['phenoage_calc','kdm_calc','kdm_nhanes','hd_calc','hd_nhanes']:
    items[n+'.R']=pin+'R/'+n+'.R'
with ThreadPoolExecutor(max_workers=5) as ex: out=list(ex.map(fetch,items.items()))
(ROOT/'downloads_more.json').write_text(json.dumps(out,indent=2)+'\n')
