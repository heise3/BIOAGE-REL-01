from pathlib import Path
import csv,json,hashlib,shutil,datetime,re,zipfile,ast

ROOT=Path.cwd();OUT=ROOT/'outputs/BIOAGE-REL-01';SRC=ROOT/'work/source'
def read(n):return list(csv.DictReader((OUT/n).open()))
def write(n,rows,fields=None):
    with (OUT/n).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]));w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def table(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','|'+'|'.join(['---']*len(headers))+'|']+['| '+' | '.join(map(str,r))+' |' for r in rows])
def fmtrow(rows,metric,group=None):
    r=next(r for r in rows if r['metric']==metric and (group is None or r.get('group')==group))
    return f"{float(r['value']):.3f}（{float(r['ci_lower']):.3f}–{float(r['ci_upper']):.3f}）" if r['value'] else 'NOT_RUN'

assert (OUT/'contract_v1.0.md').exists()
assert sha(OUT/'contract_v1.0.md')=='c7cba996dde4a84323445a6f876711e7561723077e39ce9d8ec8d33f2f27a29d'
qa=json.loads((OUT/'QA.json').read_text());assert qa['status']=='PASS' and all(x['pass'] for x in qa['checks'].values())
rel=read('B_reliability.csv');rep=read('D_replication.csv');mort=read('C_mortality.csv');sec=read('B_secondary.csv');sens=read('B_sensitivity.csv');comp=read('B_components.csv');counts={r['metric']:int(float(r['value'])) for r in read('A1_counts.csv')};cal=read('A2_calibration.csv')

# Exact column names, documented missing codes, and observed ranges.
spec=[('albumin','AMP','AMR','g/dL','888','', 'g/dL ×10 → g/L'),('alp','APPSI','APRSI','U/L','8888','','unchanged'),('creatinine','CEP','CER','mg/dL','8888','','(0.960×mg/dL−0.184)×88.4017 → µmol/L'),('glucose','SGP','SGR','mg/dL','888','','serum biochemistry; ×0.0555 → mmol/L'),('crp','CRP','CRR','mg/dL','88888','8888','log(1+CRP); retain published 0.21 nondetect fill'),('wbc','WCP','WCR','10^3 cells/µL','88888','88888','unchanged'),('lymphocyte_percent','LMPPCNT','LMRPCNT','%','88888','88888','percentage points; not fraction'),('mcv','MVPSI','MVRSI','fL','88888','88888','unchanged'),('rdw','RWP','RWR','%','88888','88888','percentage points; not RWPSI/RWRSI fraction')]
schemas={n:json.loads((SRC/(n+'_schema.json')).read_text()) for n in ['lab','labse']}
raw={n:(SRC/(n+'.dat')).read_text().splitlines() for n in ['lab','labse']}
book=[]
for marker,v1,v2,unit,fill1,fill2,transform in spec:
    for visit,n,v,fill in [(1,'lab',v1,fill1),(2,'labse',v2,fill2)]:
        start,end=schemas[n][v];vv=[line[start-1:end].strip() for line in raw[n]];valid=[float(x) for x in vv if x and x!=fill]
        book.append(dict(cohort='NHANES_III',marker=marker,visit=visit,file=n.upper(),variable=v,positions=f'{start}-{end}',unit=unit,missing_codes=('blank; '+fill if fill else 'blank only in released codebook'),missing_code_action='convert to NA; no imputation',LOD_mgdL=.3 if marker=='crp' else '',LOD_code='0.21 = below 0.3 mg/dL; rounded LOD/sqrt(2)' if marker=='crp' else 'No separate nondetect code listed for this field; retain numeric results',n_valid_all_people=len(valid),n_blank=sum(not x for x in vv),n_fill=sum(x==fill for x in vv) if fill else 0,observed_min=min(valid),observed_max=max(valid),transformation=transform,source_url='https://wwwn.cdc.gov/nchs/data/nhanes3/'+('1a/lab-acc.pdf' if visit==1 else '3a/LABSE-acc.pdf'),status='VERIFIED'))
repmap=[('albumin','LBXSAL','LB2SAL','g/dL','L40_B','L40_2_B'),('alp','LBDSAPSI','LB2SAPSI','U/L','L40_B','L40_2_B'),('creatinine','LBDSCR','LB2SCR','mg/dL','L40_B','L40_2_B'),('glucose','LBXSGL','LB2SGL','mg/dL','L40_B','L40_2_B'),('crp','LBXCRP','LB2CRP','mg/dL','L11_B','L11_2_B'),('wbc','LBXWBCSI','LB2WBCSI','10^3 cells/µL','L25_B','L25_2_B'),('lymphocyte_percent','LBXLYPCT','LB2LYPCT','%','L25_B','L25_2_B'),('mcv','LBXMCVSI','LB2MCVSI','fL','L25_B','L25_2_B'),('rdw','LBXRDW','LB2RDW','%','L25_B','L25_2_B')]
for marker,v1,v2,unit,f1,f2 in repmap:
    for visit,v,f in [(1,v1,f1),(2,v2,f2)]:
        row={k:'' for k in book[0]};row.update(cohort='NHANES_2001_2002',marker=marker,visit=visit,file=f,variable=v,unit=unit,missing_codes='SAS missing / dot',missing_code_action='convert to NA; no imputation',LOD_mgdL=.02 if marker=='crp' else '',LOD_code='0.01 low-value fill proxy; approximately 0.02 mg/dL reportable limit' if marker=='crp' else 'No separate nondetect code used',transformation='No NHANES III creatinine recalibration in 2001–2002; otherwise same score units',source_url=f'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/{f}.htm',status='VERIFIED');book.append(row)
write('A1_codebook_check.csv',book)
crp_counts=read('B_CRP_LOD.csv')

# Preserve the actual bounded supplemental search, including its access limits.
queries=[
 'site.medrxiv.org "phenotypic age" ("reliability" OR "repeatability" OR "test-retest" OR "within-person") after:2025-11-01 before:2026-11-03',
 'site.biorxiv.org "phenotypic age" ("reliability" OR "repeatability" OR "test-retest" OR "within-person") after:2025-11-01 before:2026-11-03',
 'site.medrxiv.org PhenoAge NHANES repeat second exam',
 'site.biorxiv.org PhenoAge NHANES repeat second exam',
 'site.medrxiv.org "phenotypic age" "test-retest"',
 'site.biorxiv.org "phenotypic age" "within-person"',
 'site.medrxiv.org "phenotypic age" "repeatability"',
 'site.biorxiv.org "phenotypic age" "reliability" "2026"',
 'site.medrxiv.org "phenotypic age" reliability after:2025-11-02 before:2026-11-03',
 'site.biorxiv.org "phenotypic age" repeatability after:2025-11-02 before:2026-11-03']
write('duplicate_search_log.csv',[dict(id=i+1,searched_on='2026-11-02',window_reviewed='2025-11-02 through 2026-11-02',query=q,classification='NONE_DETECTED_EXACT',scope='Returned indexed records screened for same repeat-exam datasets and same comparison; date reviewed manually; not an exhaustive registry search',limitation='Search returned older and off-domain records too; two 2026 preprint primary pages returned HTTP 403') for i,q in enumerate(queries)])

# Final source inventory is based on actual completed files, not failed first attempts.
urls={}
for fn in ['downloads.json','downloads_more.json','downloads_kidney.json']:
    for r in json.loads((SRC/fn).read_text()):urls[r['file']]=r['url']
urls['nhanes3_datafiles.html']='https://wwwn.cdc.gov/nchs/nhanes/nhanes3/datafiles.aspx'
urls['bioage_nhanes_all.R']='https://raw.githubusercontent.com/dayoonkwon/BioAge/b1f9fc02f086cd4aa74185f2335ab1366082e7fe/data-raw/nhanes_all.R'
manifest=[]
for fn,u in sorted(urls.items()):
    p=SRC/fn
    if not p.exists():continue
    pin=u.replace('/BioAge/master/','/BioAge/b1f9fc02f086cd4aa74185f2335ab1366082e7fe/')
    manifest.append(dict(file=fn,url_observed=u,url_for_reproduction=pin,bytes=p.stat().st_size,sha256=sha(p),status='PRESENT_VERIFIED',local_cache='work/source/'+fn))
write('source_manifest.csv',manifest)
for fn,expected in [('NHANES3.rda','5b88f3fa73de31d40dde4168a07de00c7c3441b5'),('NHANES3_HDTrain.rda','f3e0a8fc193a46fc75777da6fc137fe9ed4d9db1'),('phenoage_calc.R','468e95b32b23a48566a8f3e266a4a103330c25ff')]:
    bb=(SRC/fn).read_bytes();assert hashlib.sha1(b'blob '+str(len(bb)).encode()+b'\0'+bb).hexdigest()==expected

mainicc=next(r for r in rel if r['group']=='all' and r['metric']=='ICC');half=(float(mainicc['ci_upper'])-float(mainicc['ci_lower']))/2
status=json.loads((OUT/'status.json').read_text());assert half<=.03 and all(r['status']=='PASS' for r in cal)
status.update(layer='FINAL',verdict='GO',terminal=True,reason='A1 and A2 passed; primary overall ICC bootstrap 95% CI half-width <=0.03; C and eligible D completed; no verified EXACT duplicate found in bounded supplemental search.',completed_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),contract_amended=False,effects_computed=True,implementation='BioAge orig branch and default data preprocessing, pinned b1f9fc02f086cd4aa74185f2335ab1366082e7fe',duplicate_screen={'classification':'NONE_DETECTED_EXACT','exhaustive':False,'full_text_access_limits':'Two 2026 preprint pages returned HTTP 403'},failure_layers={'engineering':'RESOLVED; 12 QA checks PASS','data_executability':'PASS','measurement':'A2_PASS; B_PRECISION_PASS','scientific_result':'COMPUTED_WITH_PH_AND_GENERALIZABILITY_LIMITATIONS','incremental_prediction':'NO_CLEAR_C_INDEX_GAIN'},D={'status':'COMPUTED','N_second_day_nine_complete_all_target_ages':544,'N_rep_second_day_complete_adult_first_pregnancy_excluded':415,'N_rep_paired':407,'population':'age>=20; exclude baseline RIDEXPRG=1; repeat-visit pregnancy not observed; implementation scope disclosed','bootstrap_replicates':2000},not_estimable=['different-session reliability: N=0'],scientific_limitations=['BioAge default log(1+CRP) path is not the direct original ln(CRP) score','Convenience samples; unweighted','Unknown second-exam pregnancy in 300 women age 20–59 retained under exclude-known-pregnancy rule','Exposure PH tests positive; Cox HRs are follow-up summary associations','Mortality training-participant overlap persists with exam2 exposure','Regression dilution correction assumes classical measurement error; approximate for Cox','Public mortality time is month-scale from first exam; not exact second-exam landmark time','Replication adult sample excludes baseline RIDEXPRG=1; repeat-visit pregnancy status unavailable'])
status['B'].update(ICC=float(mainicc['value']),ci_lower=float(mainicc['ci_lower']),ci_upper=float(mainicc['ci_upper']),half_width=half)
status['A2'].update(PhenoAge_r=float(cal[0]['value']),PhenoAge_MAD_years=float(cal[1]['value']),creatinine_CVw_percent=float(cal[2]['value']))
(OUT/'status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2)+'\n')

brief=table(['门 / 指标','结果'],[
 ['A1 主样本 / 死亡','1,801 人 / 829 例；死亡链接合格率100%'],['A1 复现资格','第二天九项完整成人415人；两次完整配对407人'],['A2 分数校准','N=8,924；r=1.000000；MAD=0.000001079年'],['A2 肌酐校准','821→795人；CVw=7.5870%，目标[6.1,9.1]%'],['B 主ICC精度',f'半宽={half:.6f} ≤0.03，PASS'],['C / D','均已执行；各bootstrap 2,000次'],['终局','GO（v1.0，无合同修改）']])
summary=table(['指标','发现 NHANES III（N=1,801）','独立复现 2001–2002（N=407）'],[[label,fmtrow(rel,m,'all'),fmtrow(rep,m)] for m,label in [('ICC','ICC(2,1)'),('within_SD','个体内SD（年）'),('MDC95','MDC95（年）'),('BA_bias','检查2−检查1偏差（年）'),('BA_lower','Bland–Altman下界（年）'),('BA_upper','Bland–Altman上界（年）'),('top20_retention','加速前20%保留比例')]])
subgroups=table(['分层','N（人）','ICC（95%CI）','MDC95年（95%CI）'],[[g,next(r['n'] for r in rel if r['group']==g),fmtrow(rel,'ICC',g),fmtrow(rel,'MDC95',g)] for g in ['20-39','40-59','60+','male','female']])
senstable=table(['分析层','N','ICC（95%CI）','MDC95年（95%CI）'],[[g,next(r['n'] for r in sens if r['group']==g),fmtrow(sens,'ICC',g),fmtrow(sens,'MDC95',g)] for g in ['CRP_both_below_1mgdL','fast_both_ge8h','fast_both_lt8h','fast_discordant_8h','fast_unknown','glucose_fixed_to_pair_mean','interval_known','interval_unknown','session_different']])
comprows=sorted([r for r in comp if r['metric']=='MDC_drop_percent'],key=lambda r:-float(r['value']))
comptable=table(['固定为两次均值的成分','MDC降幅%（95%CI）','系数×个体内SD排序'],[[r['component'],f"{float(r['value']):.2f}（{float(r['ci_lower']):.2f}–{float(r['ci_upper']):.2f}）",r['delta_contribution_rank']] for r in comprows])
morttable=table(['C指标','估计（95%CI）'],[[r['metric'],f"{float(r['value']):.4f}（{float(r['ci_lower']):.4f}–{float(r['ci_upper']):.4f}）"] for r in mort])
secondary=table(['次指标','N','ICC（95%CI）','MDC95（95%CI）'],[[alg,next(r['n'] for r in sec if r['algorithm']==alg),fmtrow([r for r in sec if r['algorithm']==alg],'ICC'),fmtrow([r for r in sec if r['algorithm']==alg],'MDC95')] for alg in ['KDM_BA_Accel','HD']])
report=f'''# BIOAGE-REL-01 v1.0 筛选终局报告

执行日期：2026-11-02（America/Los_Angeles）。终局 **GO**。这是按合同允许的 BioAge 原始分数分支及其默认输入预处理完成的可执行性与精度筛选；不把 GO 解释为临床有效性、原始 ln(CRP) 版的等价验证、因果结论或明确预测增益。

{brief}

## 工程

用户原文逐字节冻结；SHA-256：`c7cba996dde4a84323445a6f876711e7561723077e39ce9d8ec8d33f2f27a29d`。未产生v1.1或v1.2，未改科学阈值、算法集合、人群或终局标准。公式/单位选择在B效应前写入 implementation_lock.md；后续工程修正见 implementation_audit.md。

仅取得公开NCHS记录、公开2019死亡链接和公开BioAge参考数据。未使用2000年RDC第二天数据，未登录门户或依赖外部回信。输入文件最终SHA-256及URL见 source_manifest.csv。BioAge提交固定为 `b1f9fc02f086cd4aa74185f2335ab1366082e7fe`，两个RDA及分数源代码的Git blob SHA也与该仓库相符。

已修复的工程问题：首次大文件下载不完整并重取；PDF读取依赖缺失；HD合并时同名空列冲突；2001–2002首检生化变量实际为 LBDSCR/LBDSAPSI；肾功能校准极值规则的初始读法错误。均在相关结果定稿前修复，未改门槛。12项独立单位、字典计数、校准、两向ANOVA、bootstrap完成度检查全部通过，详见 QA.json。

NHANES III生化葡萄糖采用SGP/SGR，以复现BioAge输入；不以OGTT血浆值代替。NCHS一般实验室用途偏好专门检测值，但本合同选择的是已发表分数的既定输入。分数采用稳定的Gompertz反解，避免高风险记录出现1−概率的浮点消失；不截尾原始完整成人。全部bootstrap按人配对抽取2,000次，百分位95%CI。残差模型只在首检15,635名完整成人拟合一次，再用于两次检查和复现两次检查。

## 数据

LABSE、EXAMSE各实际2,596人，首检LAB 29,314人；重复检查中成人2,160人。成人九项两次完整1,814人，排除已知第二次怀孕13人，得到主样本1,801人。主样本死亡链接全部合格，死亡829例。300名20–59岁女性的第二次孕情不明确，按“排除已知怀孕”保留并单列，不能将其描述为已确认未孕。

九项重复变量依次为 AMR、APRSI、CER、SGR、CRR、WCR、LMRPCNT、MVRSI、RWR；精确位置、单位、各字段填充值和非缺失计数见 A1_codebook_check.csv。CRP低于0.3 mg/dL的发布值0.21保留，不再补一次LOD/√2，也不当缺失；888、8888、88888等仅按字段定义转NA。

主样本两次会话编码全相同，`session_different`结构性N=0，记NOT_RUN。这只能说明指定编码没有可用对照层，不能据此估计时段差异效应。两次间隔已知1,161人，缺失640人；间隔缺失未排除主ICC。

复现的三份第二天文件九项交集完整544人（含16–19岁者）；按≥20岁并排除DEMO_B的首检RIDEXPRG=1留下415人，再与首检九项完整记录配对得到407人。RIDEXPRG是首检孕情，并非复检孕情；本次D仅描述这个成人实现样本，不能声称已严格核实第二次非孕条件。A1的N_rep=415和D的实际配对N=407是不同计数阶段，不混写。

代表性比较使用全部有首检LAB记录的成人，参与第二次检查2,160人、未参与14,870人。参与者首检平均年龄49.03岁、未参与48.76岁；有完整分数者分别2,031和13,604人，其PhenoAgeAccel均值分别−0.212和0.032年。性别、四类种族/族裔比例见 B_representativeness.csv。这些小的均值差不构成总体代表性证明；所有模型均不加权。

CRP低值编码比例：NHANES III首检60.41%、复检62.69%（0.21编码）；2001–2002首检2.70%、复检3.44%（0.01低值编码、约0.02 mg/dL最低可报告水平）。两队列人群与检测条件均不同，不能把可靠性差异全归于检测下限。

## 科学结果

### 主指标与独立复现

{summary}

主ICC区间半宽为{half:.6f}，通过冻结的总体主指标精度门。ICC较高与MDC约7年可以同时成立：ICC取决于人与人之间的分布跨度，MDC反映本样本中两次测量差值的绝对尺度。复现ICC较高并不表示个体重测差值更小。MDC是两次重测总波动的统计阈值，含分析、短期生物变化及采样条件影响，不代表“真实衰老速率”或个人治疗反应阈值。

### 年龄与性别分层

{subgroups}

### 冻结竞争解释与敏感性

{senstable}

8小时用于空腹分层，未知单列；不按结果切换空腹门槛。两次同会话编码层与主样本相同，点估计相同；额外bootstrap导致CI有限Monte Carlo差异。CRP任一次≥1.0 mg/dL排除后的MDC为6.37年；这不能把残余变异分解为感染、其他真实变化与分析误差。将葡萄糖固定为两次均值后的MDC约4.88年，是预定的计算反事实，不是用新算法替换主指标。

{comptable}

固定成分后MDC降幅不可相加，因为成分间协方差仍在。精确变换后“系数×个体内SD”排序与MDC降幅排序不同；B_components.csv同时给出两者及其bootstrap CI，未据此删除成分。

### 次指标

{secondary}

KDM采用BioAge默认九项（FEV1、SBP、总胆固醇、HbA1c、白蛋白、肌酐、log1p CRP、ALP、BUN）及首检NHANES III按性别30–75岁、非孕参考拟合，参数固定后投射两次；只分析两次九项皆完整1,044人。其残差模型亦仅在首检数据拟合。HD固定BioAge内置按性别的年轻健康参考人群，两次记录合并投射以保持相同缩放；其MDC为分数单位而非年。不同次指标的样本、尺度不同，不按ICC数值单独排算法优劣，也不将次指标CI宽度当新增主门槛。

### 死亡与回归稀释

三种暴露在同一1,801人、829例死亡中比较；年龄、性别、四类种族/族裔调整。共同1 SD=首检PhenoAgeAccel的6.77464年，以保持单次与均值HR可比。

{morttable}

λ定义为 Cov(检查1,检查2)/Var(检查2)，校正检查2的logHR除以λ；bootstrap同时重估λ和模型。它依赖经典测量误差、重复误差结构等假设，是Cox情境的近似校正，不能从两次重测验证这些假设。均值模型相对单次模型的表观C-index增益CI均跨0，未建立明确判别增益。

PH检验：检查1暴露p=0.000573，检查2主暴露p=0.02755，均值暴露p=0.004346；全局p分别0.01060、0.12814、0.04607。不能把主模型的全局p>0.05解释为暴露满足PH。单一HR仅作为随访期间汇总关联，未据此另换模型求更好结果。完整检验见 C_PH_tests.csv。

检查2血样未参与原分数训练，但原训练中的受试者及死亡结局重叠并未因更换血样消失，因此死亡关联不是独立外部验证。公开随访为月尺度PERMTH_EXM，起点仍为原检查；间隔缺失使精确检查2 landmark分析不可执行。本合同的C比较保持统一起点，不将其声称为精确landmark死亡估计。

## 是否继续投入

**GO**。A1、A2、B总体主指标精度门通过，C与符合资格的D均完成。可靠性高低和C-index未显示明确增益均不是冻结合同的否决条件。

可以继续围绕“较高排序可靠性与约7年的绝对重测变化阈值并存”整理稿件；对死亡回归稀释仅作假设依赖的补充。当前没有临床可用性、干预响应或普适MDC的证据，也没有把C-index微小点增益解释为预测改进。

失败分层：工程问题已解决；数据可执行；主测量与精度门通过；科学估计已算但保留PH、选择与误差模型限制；增量预测未见明确增益。空时段层NOT_RUN属于结构性无对照，不是阴性科学结果。

近12个月bioRxiv/medRxiv补检未核实到“同重复数据、同重测/稀释比较”的EXACT命中，检索式及范围见 duplicate_search_log.csv。检索返回一些旧文及偏离站点的记录，两份2026预印本原始网页403；这里只能认定本轮限定检索未发现EXACT，不能据此声称穷尽查新或确定首创。未把用户提供的既往桌面检索伪称为本次重新全文审查。

## 与合同的偏离

1. 未修改合同或启用v1.1；未更换主算法来改善效应。采用第1.7条允许的BioAge实现，其默认CRP输入为log(1+CRP)、肌酐为已发表标准化值。它与合同风险文字中直接称“ln(CRP)”存在实质口径差异，已在效应前锁定并披露；本报告仅支持该BioAge路径，**不能被引述为原始直接ln(CRP)版Levine分数的结果**。若未来要求那个版本，需另行明确科学合同，不能把本次同版本重跑成另一指标。
2. A2内置对照只能在有phenoage0的8,924人比较，因为包本身只为特定年龄/空腹子集存储该列；已在全部首检完整成人复算，未将A2对照子集当作主样本。没有为通过校准把主分析限制为空腹或≤84岁。
3. 年龄采用首检访谈年龄并固定用于两次，避免生日与缺失间隔混入重测比较；复现也使用其首检年龄。复现投射同一个首检NHANES III残差模型，未在复现结果上重新优化。
4. HD必须用无量纲分数单位报告，不能套用“年”；两次共同缩放保证绝对一致性比较。KDM/HD不是主精度门的替代品。
5. 时段不同层N=0而不可估计；主样本孕情未知与公开死亡时间分辨率限制均如实保留。D的人群口径在效应前按成人和首检孕情实现：第二天全目标年龄九项完整544人，成人且首检未报告怀孕415人，两次完整407人。合同D未单独细化孕情口径，现有RIDEXPRG不能确认复检孕情；该复现范围限制不改变主GO门，也不把首检孕情伪称复检孕情。未以擅自排除/补齐或另换变量消除限制。
6. 肾功能校准初始极值实现错误已更正为原文的校准后差值，完整重现821→795与各标志物剔除计数，阈值未改；首次HD列名冲突与复现变量名错误均为工程修复。所有最终结果取修复后输出。

来源： [NCHS NHANES III目录](https://wwwn.cdc.gov/nchs/nhanes/nhanes3/datafiles.aspx)、[LABSE字典](https://wwwn.cdc.gov/nchs/data/nhanes3/3a/LABSE-acc.pdf)、[2001–2002首检生化字典](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/L40_B.htm)、[第二天血常规字典](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/L25_2_B.htm)、[第二天生化字典](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/L40_2_B.htm)、[第二天CRP字典](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/L11_2_B.htm)、[NCHS死亡链接说明](https://www.cdc.gov/nchs/linked-data/mortality-files/)、[肌酐校准研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC3628297/)、[cystatin C校准](https://pmc.ncbi.nlm.nih.gov/articles/PMC3703771/)、[BioAge固定提交](https://github.com/dayoonkwon/BioAge/tree/b1f9fc02f086cd4aa74185f2335ab1366082e7fe)。来源仅支持相应方法/字典；本报告的估计值来自本次本地计算。
'''
(OUT/'REPORT.md').write_text(report)

# Copy reviewable code; the cache remains local and is not duplicated into the aggregate archive.
code=OUT/'code';code.mkdir(exist_ok=True)
for n in ['fetch_sources.py','fetch_more.py','fetch_kidney.py','parse_nhanes.py','a1.R','a2.R','b_analysis.R','c_mortality.R','d_replication.R','qa.R','finalize.py']:
    shutil.copy2(ROOT/'work'/n,code/n)
for p in (ROOT/'work').glob('*.py'):ast.parse(p.read_text())
(OUT/'README.md').write_text('''# BIOAGE-REL-01 冻结筛选包

首先阅读 REPORT.md 和 status.json。contract_v1.0.md 与 implementation_lock.md 保留效应前冻结记录，implementation_audit.md说明工程校准修正。CSV均为按人统计的聚合输出；本包不含参加者逐行记录。

code/含完整下载、固定宽度读入、A1/A2/B/C/D和QA源代码；R依赖版本见sessionInfo.txt。输入公共文件缓存保留在本聊天工作区work/source，最终URL/字节数/哈希见source_manifest.csv。BioAge输入按固定Git提交下载。

可在任意空目录运行 `bash /本包绝对路径/code/run_reproduction.sh`；该脚本会在当前目录创建work/和outputs/BIOAGE-REL-01/，下载公开输入并复算。需要Python 3、R及haven/readr/dplyr/jsonlite/survey/survival。PDF文本已用于本次人工核对，数值复算不依赖PDF抽取。重取输入后先按source_manifest.csv校验；哈希漂移应停留在工程核查，不可把新版本输入默认为同一冻结数据。

manifest.sha256为本包所有交付文件的校验清单；ZIP是同一目录的便携副本。不得把GO、一次校准或探索性死亡关联引述为临床验证。
''')
runner='''#!/bin/bash
set -euo pipefail
BIOAGE_PACKAGE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p work/source outputs/BIOAGE-REL-01
cp "$BIOAGE_PACKAGE_DIR"/code/*.py "$BIOAGE_PACKAGE_DIR"/code/*.R work/
cp "$BIOAGE_PACKAGE_DIR"/contract_v1.0.md "$BIOAGE_PACKAGE_DIR"/contract_v1.0.sha256 "$BIOAGE_PACKAGE_DIR"/implementation_lock.md "$BIOAGE_PACKAGE_DIR"/implementation_lock.md.sha256 "$BIOAGE_PACKAGE_DIR"/implementation_audit.md outputs/BIOAGE-REL-01/
python3 - <<'PY'
from pathlib import Path
import json,hashlib
p=Path('outputs/BIOAGE-REL-01');h=hashlib.sha256((p/'contract_v1.0.md').read_bytes()).hexdigest()
(p/'status.json').write_text(json.dumps(dict(project='BIOAGE-REL-01',contract_version='1.0',contract_sha256=h,layer='A1',verdict='NOT_RUN',effects_computed=False)))
PY
python3 work/fetch_sources.py
python3 work/fetch_more.py
python3 work/fetch_kidney.py
python3 - "$BIOAGE_PACKAGE_DIR/source_manifest.csv" <<'PY'
import csv,hashlib,sys
from pathlib import Path
for r in csv.DictReader(open(sys.argv[1])):
 p=Path('work/source')/r['file']
 if r['file']=='nhanes3_datafiles.html':continue
 assert p.exists(),r['file']+' missing'
 assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],r['file']+' source changed: inspect before proceeding'
PY
python3 work/parse_nhanes.py
Rscript work/a1.R
Rscript work/a2.R
Rscript work/b_analysis.R
Rscript work/c_mortality.R
Rscript work/d_replication.R
Rscript work/qa.R
python3 work/finalize.py
'''
(code/'run_reproduction.sh').write_text(runner)
filelist=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.sha256')
(OUT/'manifest.sha256').write_text(''.join(f'{sha(p)}  {p.relative_to(OUT)}\n' for p in filelist))
zp=OUT.parent/'BIOAGE-REL-01_v1.0_GO.zip'
with zipfile.ZipFile(zp,'w',compression=zipfile.ZIP_DEFLATED) as z:
    for p in sorted(OUT.rglob('*')):
        if p.is_file():z.write(p,Path('BIOAGE-REL-01')/p.relative_to(OUT))
print(json.dumps(dict(verdict='GO',deliverables=len(filelist)+1,zip_bytes=zp.stat().st_size,report=OUT.joinpath('REPORT.md').as_posix(),primary_ICC=float(mainicc['value']),CI_half_width=half),ensure_ascii=False))
