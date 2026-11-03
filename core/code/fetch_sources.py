from pathlib import Path
from urllib.request import urlopen, Request
from concurrent.futures import ThreadPoolExecutor
import hashlib, json, datetime, time

ROOT=Path(__file__).resolve().parent/'source'
ROOT.mkdir(parents=True,exist_ok=True)
base='https://wwwn.cdc.gov/nchs/data/nhanes3/'
items={f'{n}.{ext}':base+f'3a/{n}.{ext}' for n in ['labse','lab2se','examse'] for ext in ['dat','sas']}
items.update({f'{n}-acc.pdf':base+f'3a/{n.upper()}-acc.pdf' for n in ['labse','lab2se','examse']})
items.update({'readme_3a.txt':base+'3a/readme.txt','lab.dat':base+'1a/lab.dat','lab.sas':base+'1a/lab.sas','lab-acc.pdf':base+'1a/lab-acc.pdf'})
for n in ['L11_2_B','L25_2_B','L40_2_B','DEMO_B','L11_B','L25_B','L40_B']:
    for ext in ['htm','xpt']:
        items[n+'.'+ext]='https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/'+n+'.'+ext

def fetch(pair):
    name,url=pair; p=ROOT/name
    try:
        if not p.exists():
            for attempt in range(3):
                try:
                    with urlopen(Request(url,headers={'User-Agent':'BIOAGE-REL-01 public-data feasibility audit'}),timeout=180) as r:
                        length=r.headers.get('Content-Length');part=p.with_suffix(p.suffix+'.part')
                        with part.open('wb') as f:
                            while True:
                                chunk=r.read(1024*1024)
                                if not chunk: break
                                f.write(chunk)
                        if length and part.stat().st_size!=int(length):raise ValueError('Incomplete HTTP body')
                        b=part.read_bytes()
                        if not name.endswith('.htm') and b[:100].lstrip().lower().startswith((b'<!doctype html',b'<html')):
                            raise ValueError('HTML received instead of data')
                        part.replace(p)
                    break
                except Exception:
                    if attempt==2: raise
                    time.sleep(1)
        b=p.read_bytes()
        d={'file':name,'url':url,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'downloaded_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'OK'}
    except Exception as e:
        d={'file':name,'url':url,'status':'ERROR','error':str(e)}
    print(json.dumps(d),flush=True)
    return d

if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=5) as ex: out=list(ex.map(fetch,items.items()))
    (ROOT/'downloads.json').write_text(json.dumps(out,indent=2)+'\n')
