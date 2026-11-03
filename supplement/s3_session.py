from pathlib import Path
import re,csv,json,hashlib
O=Path('outputs/BIOAGE-REL-01/S');SRC=Path('work/source')
def schema(n):
 s=(SRC/(n+'.sas')).read_text();block=re.search(r'\bINPUT\b(.*?);',s,re.S).group(1)
 return {m[0]:(int(m[1]),int(m[2] or m[1])) for m in re.findall(r'^\s*(\w+)\s+(\d+)(?:-(\d+))?\s*$',block,re.M)}
pdef=schema('lab');rdef=schema('labse');assert 'MXPSESSR' not in rdef
def rows(n,cols,defs):
 result={}
 for line in (SRC/(n+'.dat')).read_text().splitlines():
  r={c:line[defs[c][0]-1:defs[c][1]].strip() for c in cols};assert r['SEQN'] not in result;result[r['SEQN']]=r
 return result
p=rows('lab',['SEQN','HSAGEIR','MXPSESSR'],pdef);r=rows('labse',['SEQN','MXRSESSR'],rdef)
adult=[(p[k]['MXPSESSR'],v['MXRSESSR']) for k,v in r.items() if k in p and p[k]['HSAGEIR'] and int(p[k]['HSAGEIR'])>=20]
assert len(adult)==2160
valid=[x for x in adult if x[0] in ['1','2','3'] and x[1] in ['1','2','3']]
same=sum(x==y for x,y in valid)
cross=[dict(label='写作期预定补充分析',population='all LABSE adults age>=20; not restricted to main complete cases',MXPSESSR=i,MXRSESSR=j,n_people=sum(x==str(i) and y==str(j) for x,y in adult),total_adults=len(adult),total_both_valid=len(valid),row_source='LAB primary MXPSESSR',column_source='LABSE repeat MXRSESSR') for i in [1,2,3] for j in [1,2,3]]
with (O/'S3_session_crosstab.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(cross[0]));w.writeheader();w.writerows(cross)
positions=[]
for variable,n,defs in [('MXPSESSR','lab',pdef),('MXRSESSR','labse',rdef)]:
 start,end=defs[variable];positions.append(dict(variable=variable,SAS_file=n+'.sas',data_file=n+'.dat',start_position_1based=start,end_position_1based=end,source_sha256=hashlib.sha256((SRC/(n+'.sas')).read_bytes()).hexdigest(),same_name_in_other_file=variable in (rdef if n=='lab' else pdef)))
with (O/'S3_variable_positions.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(positions[0]));w.writeheader();w.writerows(positions)
doc=(SRC/'labse-acc.txt').read_text();at=doc.index('The second exams were scheduled');page=int(re.findall(r'=== PDF page (\d+) ===',doc[:at])[-1])
summary=dict(status='COMPUTED',n_all_adults=len(adult),n_both_valid=len(valid),n_invalid_or_missing=len(adult)-len(valid),n_same=same,n_different=len(valid)-same,all_adults_100_percent_same=(same==len(adult)),interpretation='按设计同时段安排' if same==len(adult) else 'Investigate discrepancy versus main sample',documentation_qualification='NCHS says when possible at the same time of day; equal session categories do not prove identical clock time of venipuncture.',source='https://wwwn.cdc.gov/nchs/data/nhanes3/3a/LABSE-acc.pdf',source_pdf_page_1based=page)
(O/'S3_session_audit.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n');print(summary)
