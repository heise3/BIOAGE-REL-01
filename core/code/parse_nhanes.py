from pathlib import Path
import json,csv,re
root=Path(__file__).resolve().parent
src=root/'source'
choices={
 'lab':['SEQN','HSAGEIR','HSSEX','DMARETHN','DMARACER','MXPSESSR','PHPFAST','AMP','AMPSI','APPSI','CEP','CEPSI','SGP','SGPSI','CRP','WCP','LMPPCNT','MVPSI','RWP','BUP','UAP','GHP','TCP','CHP'],
 'labse':['SEQN','MXRSESSR','MXPRDAYS','PHRFAST','AMR','AMRSI','APRSI','CER','CERSI','SGR','SGRSI','CRR','WCR','LMRPCNT','MVRSI','RWR','BUR','UAR','GHR','TCR','CHR'],
 'examse':['SEQN','MARF12','MARF12R','PERPREG','PERMNK1R','SPRFEV1'],
}
if __name__=='__main__':
 for name,wanted in choices.items():
    sas=(src/(name+'.sas')).read_text()
    block=re.search(r'\bINPUT\b(.*?);',sas,re.S).group(1)
    schema={m[0]:[int(m[1]),int(m[2] or m[1])] for m in re.findall(r'^\s*(\w+)\s+(\d+)(?:-(\d+))?\s*$',block,re.M)}
    (src/(name+'_schema.json')).write_text(json.dumps(schema,indent=2))
    wanted=[x for x in wanted if x in schema]
    lines=(src/(name+'.dat')).read_text().splitlines()
    rows=[]
    for line in lines:
        row={}
        for k in wanted:
            a,b=schema[k];v=line[a-1:b].strip()
            if k not in ['SEQN','HSAGEIR','HSSEX','DMARETHN','DMARACER','MXPSESSR','MXRSESSR','MXPRDAYS','MARF12','MARF12R','PERPREG'] and re.fullmatch('8{2,}',v):v=''
            row[k]=v
        rows.append(row)
    assert len({r['SEQN'] for r in rows})==len(rows)
    with (root/(name+'_selected.csv')).open('w') as f:
        w=csv.DictWriter(f,fieldnames=wanted);w.writeheader();w.writerows(rows)
    print(name,len(rows),min(map(len,lines)),max(map(len,lines)))
