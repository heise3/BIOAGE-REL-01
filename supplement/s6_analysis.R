suppressPackageStartupMessages({library(haven);library(jsonlite)})
O<-'outputs/BIOAGE-REL-01/S/S6';dir.create('work/S6',showWarnings=FALSE)
Rboot<-2000L;set.seed(202611026)
s<-readRDS('work/S/s1_results.rds');b<-readRDS('work/b_results.rds');cal<-readRDS('work/a2.rds');rr<-readRDS('work/d_results.rds');a<-readRDS('work/a1.rds')
stopifnot(nrow(rr$d)==407,nrow(a$main)==1801,b$base_n==15635,abs(coef(b$resmod)[1]+5.1991589906396616)<1e-12,abs(coef(s$model_ln)[1]+6.907959472007664)<1e-12)
icc21<-b$icc21;pheno<-cal$pheno
bioms<-c('albumin_gL','alp','creat_umol','glucose_mmol','lncrp','wbc','lymph','mcv','rdw')
coefs<-c(-.03359355,.001868778,.009506491,.1953192,.09536762,.05542406,-.01199984,.02676401,.3306156)/.090165
names(coefs)<-bioms
labels<-c(albumin_gL='白蛋白',alp='ALP',creat_umol='肌酐',glucose_mmol='葡萄糖',lncrp='CRP（按分数变换）',wbc='白细胞',lymph='淋巴细胞百分比',mcv='MCV',rdw='RDW')
ci<-function(v){if(!any(is.finite(v)))return(c(NA_real_,NA_real_));quantile(v[is.finite(v)],c(.025,.975),names=FALSE)}
bootindices<-lapply(s$cohorts,function(z)replicate(Rboot,sample.int(nrow(z$z1),nrow(z$z1),TRUE)))
component_rows<-list();kk<-1L;identity_errors<-list()
for(cohort in names(s$cohorts))for(version in c('log1p','ln')){
 zz<-s$cohorts[[cohort]];z1<-zz$z1;z2<-zz$z2
 if(version=='ln'){z1$lncrp<-log(z1$crp);z2$lncrp<-log(z2$crp)}
 xs<-s$scores[[cohort]][[version]];dd<-xs$x2-xs$x1;n<-length(dd);idx<-bootindices[[cohort]]
 contrib<-sweep(as.matrix(z2[,bioms]-z1[,bioms]),2,coefs,'*')
 identity_errors[[paste(cohort,version)]]<-max(abs(rowSums(contrib)-dd));stopifnot(max(abs(rowSums(contrib)-dd))<1e-10)
 weights<-apply(contrib,2,sd)/sqrt(2);ranks<-rank(-weights,ties.method='min')
 for(j in seq_along(bioms)){
  cc<-contrib[,j];ff<-dd-cc
  calc<-function(ii){d0<-dd[ii];c0<-cc[ii];f0<-ff[ii];v0<-var(d0);c(MDC95_original=1.96*sd(d0),MDC95_fixed=1.96*sd(f0),MDC_drop_years=1.96*(sd(d0)-sd(f0)),MDC_drop_percent=100*(1-sd(f0)/sd(d0)),coefficient_times_within_SD=sd(c0)/sqrt(2),covariance_allocated_variance=cov(c0,d0)/2,difference_variance_drop_percent=100*(1-var(f0)/v0))}
  obs<-calc(seq_len(n));bt<-apply(idx,2,calc);ci0<-t(apply(bt,1,ci));row<-data.frame(cohort=cohort,score=version,component=bioms[j],component_zh=unname(labels[j]),n=n,delta_contribution_rank=unname(ranks[j]),bootstrap_replicates=Rboot,status='COMPUTED')
  for(m in names(obs)){row[[m]]<-unname(obs[m]);row[[paste0(m,'_ci_lower')]]<-ci0[m,1];row[[paste0(m,'_ci_upper')]]<-ci0[m,2]}
  # Verify the algebraic fixed-component result against rescoring identical means.
  q1<-z1;q2<-z2;q1[[bioms[j]]]<-q2[[bioms[j]]]<-(z1[[bioms[j]]]+z2[[bioms[j]]])/2
  model<-if(version=='ln')s$model_ln else b$resmod
  fixed1<-pheno(q1)-predict(model,q1);fixed2<-pheno(q2)-predict(model,q2)
  stopifnot(max(abs((fixed2-fixed1)-ff))<1e-10)
  component_rows[[kk]]<-row;kk<-kk+1L
 }
}
comp<-do.call(rbind,component_rows)
idcols<-c('score','component','component_zh')
left<-comp[comp$cohort=='NHANES_III',setdiff(names(comp),'cohort')]
right<-comp[comp$cohort=='NHANES_2001_2002',setdiff(names(comp),'cohort')]
names(left)[!names(left)%in%idcols]<-paste0('NHANES_III_',names(left)[!names(left)%in%idcols])
names(right)[!names(right)%in%idcols]<-paste0('NHANES_2001_2002_',names(right)[!names(right)%in%idcols])
wide<-merge(left,right,by=idcols,sort=FALSE);wide<-wide[order(match(wide$score,c('log1p','ln')),wide$NHANES_2001_2002_delta_contribution_rank),]
write.csv(wide,file.path(O,'S6a_components_replication.csv'),row.names=FALSE,na='')
old<-read.csv('outputs/BIOAGE-REL-01/B_components.csv');old<-old[old$metric!='',]
old_err<-0
for(i in seq_len(nrow(old))){z<-comp[comp$cohort=='NHANES_III'&comp$score=='log1p'&comp$component==old$component[i],];old_err<-max(old_err,abs(z[[old$metric[i]]]-old$value[i]))}
stopifnot(old_err<1e-9)
cat('S6a complete, original B maximum point-estimate difference:',old_err,'\n');flush.console()

# S6b: raw concentrations; each analyte has its own paired Tukey filter.
z1<-s$cohorts$NHANES_2001_2002$z1;z2<-s$cohorts$NHANES_2001_2002$z2
refs<-data.frame(analyte=c('creat','glucose','albumin','alp'),unit=c('mg/dL','mg/dL','g/dL','U/L'),ref_cvw=c(6.7,7.5,3.5,8.1),ref_cva=c(4.6,1.5,1.3,2.1),ref_n=c(823L,789L,831L,828L))
fence<-function(x,y){qs<-quantile(c(x,y),c(.25,.75),type=7,names=FALSE);c(lower=qs[1]-3*diff(qs),upper=qs[2]+3*diff(qs))}
cvcalc<-function(x,y,cva,screen=FALSE){
 if(screen){f<-fence(x,y);keep<-x>=f[1]&x<=f[2]&y>=f[1]&y<=f[2];x<-x[keep];y<-y[keep]}
 avg<-mean((x+y)/2);w<-sqrt(mean((y-x)^2)/2);cvt<-100*w/avg;v<-cvt^2-cva^2
 c(mean_concentration=avg,total_within_RMS=w,CVt_percent=cvt,CVw_percent=if(v>=0)sqrt(v)else NA_real_,CVt_squared_minus_CVa_squared=v)
}
cvrows<-list();kk<-1L
for(j in seq_len(nrow(refs))){v<-refs$analyte[j];x<-z1[[v]];y<-z2[[v]];f<-fence(x,y);keep<-x>=f[1]&x<=f[2]&y>=f[1]&y<=f[2]
 for(screen in c(FALSE,TRUE)){
  obs<-cvcalc(x,y,refs$ref_cva[j],screen);bt<-apply(bootindices$NHANES_2001_2002,2,function(ii)cvcalc(x[ii],y[ii],refs$ref_cva[j],screen));ci0<-t(apply(bt,1,ci))
  row<-data.frame(analyte=v,unit=refs$unit[j],screen=if(screen)'after_Tukey_3IQR'else'before_Tukey',n_original=length(x),n=if(screen)sum(keep)else length(x),n_removed=if(screen)sum(!keep)else 0L,tukey_lower=f[1],tukey_upper=f[2],NHSR21_CVw_percent=refs$ref_cvw[j],NHSR21_CVa_percent=refs$ref_cva[j],NHSR21_analyte_n=refs$ref_n[j],bootstrap_replicates=Rboot,bootstrap_nonestimable=sum(!is.finite(bt['CVw_percent',])),status=if(is.finite(obs['CVw_percent']))'COMPUTED'else'NOT_ESTIMABLE',reason=if(is.finite(obs['CVw_percent']))''else'CVt squared is smaller than fixed literature CVa squared; biological CV is not identifiable under this subtraction')
  for(m in names(obs)){row[[m]]<-unname(obs[m]);row[[paste0(m,'_ci_lower')]]<-ci0[m,1];row[[paste0(m,'_ci_upper')]]<-ci0[m,2]}
  row$difference_CVw_percentage_points<-unname(obs['CVw_percent'])-refs$ref_cvw[j]
  row$classification<-if(!is.finite(row$difference_CVw_percentage_points))'不可估计'else if(abs(row$difference_CVw_percentage_points)<=1.5)'一致'else'不一致'
  row$source_url<-'https://www.cdc.gov/nchs/data/nhsr/nhsr021.pdf'
  row$method_note<-'One-way within-person RMS includes visit bias; pooled-day concentration Tukey fences, remove whole pair; fences recomputed in paired bootstrap; CVa fixed from NHSR21, not locally estimated.'
  cvrows[[kk]]<-row;kk<-kk+1L
 }
}
cvres<-do.call(rbind,cvrows);write.csv(cvres,file.path(O,'S6b_cvw_vs_nhsr21.csv'),row.names=FALSE,na='')
cat('S6b complete\n');flush.console()

# S6c: aggregate raw paired structure and requested score sensitivity.
x<-z1$creat;y<-z2$creat;dd<-y-x;n<-length(dd)
deming<-function(x,y){sxx<-var(x);syy<-var(y);sxy<-cov(x,y);if(!is.finite(sxy)||abs(sxy)<1e-14)return(c(intercept=NA_real_,slope=NA_real_));sl<-(syy-sxx+sqrt((syy-sxx)^2+4*sxy^2))/(2*sxy);c(intercept=mean(y)-sl*mean(x),slope=sl)}
calc_pair<-function(ii)c(day1_mean=mean(x[ii]),day2_mean=mean(y[ii]),mean_difference=mean(dd[ii]),SD_difference=sd(dd[ii]),Deming_intercept=unname(deming(x[ii],y[ii])[1]),Deming_slope=unname(deming(x[ii],y[ii])[2]),proportion_abs_difference_gt0_3=mean(abs(dd[ii])>.3))
pairrows<-list();kk<-1L
addrows<-function(group,score,n0,obs,bt,extra=list()){
 ci0<-t(apply(bt,1,ci));df<-data.frame(group=group,score=score,metric=names(obs),n=n0,value=as.numeric(obs),ci_lower=ci0[,1],ci_upper=ci0[,2],bootstrap_valid=rowSums(is.finite(bt)),bootstrap_replicates=Rboot,status='COMPUTED',reason='')
 for(v in names(extra))df[[v]]<-extra[[v]]
 df
}
obs<-calc_pair(seq_len(n));bt<-apply(bootindices$NHANES_2001_2002,2,calc_pair)
pairrows[[kk]]<-addrows('all_407','raw_creatinine',n,obs,bt);kk<-kk+1L
big<-abs(dd)>.3
pairrows[[kk]]<-data.frame(group='all_407',score='raw_creatinine',metric='n_abs_difference_gt0_3',n=n,value=sum(big),ci_lower=NA_real_,ci_upper=NA_real_,bootstrap_valid=NA_integer_,bootstrap_replicates=0L,status='COMPUTED',reason='');kk<-kk+1L
keep<-which(!big)
relcalc<-function(v1,v2)c(ICC=icc21(v1,v2),within_SD=sd(v2-v1)/sqrt(2),MDC95=1.96*sd(v2-v1))
for(version in c('log1p','ln')){
 xs<-s$scores$NHANES_2001_2002[[version]]
 for(g in c('all_407','exclude_abs_creatinine_difference_gt0_3')){
  ii<-if(g=='all_407')seq_len(n)else keep;obs<-relcalc(xs$x1[ii],xs$x2[ii]);idx<-if(g=='all_407')bootindices$NHANES_2001_2002 else replicate(Rboot,sample(ii,length(ii),TRUE));bt<-apply(idx,2,function(jj)relcalc(xs$x1[jj],xs$x2[jj]))
  pairrows[[kk]]<-addrows(g,version,length(ii),obs,bt);kk<-kk+1L
 }
}
days<-rr$d$LB2DAY;stopifnot(length(days)==n,all(is.finite(days)),all(days>=0))
q<-quantile(days,c(1/3,2/3),type=7,names=FALSE)
group<-cut(days,c(-Inf,q,Inf),labels=c('interval_tertile_1','interval_tertile_2','interval_tertile_3'),right=TRUE)
pairrows[[kk]]<-addrows('interval_relation_all_407','raw_creatinine',n,c(Spearman_rho=cor(days,dd,method='spearman')),matrix(apply(bootindices$NHANES_2001_2002,2,function(ii)cor(days[ii],dd[ii],method='spearman')),nrow=1,dimnames=list('Spearman_rho',NULL)));kk<-kk+1L
for(g in levels(group)){
 ii<-which(group==g);calc<-function(jj)c(mean_difference=mean(dd[jj]),SD_difference=sd(dd[jj]),LB2DAY_min=min(days[jj]),LB2DAY_max=max(days[jj]),LB2DAY_mean=mean(days[jj]));obs<-calc(ii);bt<-replicate(Rboot,calc(sample(ii,length(ii),TRUE)));pairrows[[kk]]<-addrows(g,'raw_creatinine',length(ii),obs,bt);kk<-kk+1L
}
pairs<-do.call(rbind,pairrows);pairs$unit<-ifelse(pairs$metric%in%c('ICC','Deming_slope','Spearman_rho'),'dimensionless',ifelse(pairs$metric=='proportion_abs_difference_gt0_3','proportion',ifelse(pairs$metric=='n_abs_difference_gt0_3','people',ifelse(grepl('LB2DAY',pairs$metric),'days',ifelse(pairs$score%in%c('log1p','ln'),'years','mg/dL')))))
pairs$deming_error_variance_ratio<-ifelse(grepl('Deming',pairs$metric),1,NA_real_)
pairs$interval_q1<-q[1];pairs$interval_q2<-q[2]
write.csv(pairs,file.path(O,'S6c_creatinine_pairs.csv'),row.names=FALSE,na='')
cat('S6c complete\n');flush.console()

# S6d: no second-exam fasting duration in the public inventory or paired files.
ph<-as.data.frame(read_xpt('work/S6_sources/PH_B.xpt'));stopifnot(!anyDuplicated(ph$SEQN),all(c('PHAFSTHR','PHAFSTMN')%in%names(ph)))
ph<-ph[match(rr$d$SEQN,ph$SEQN),];first_known<-sum(is.finite(ph$PHAFSTHR)&is.finite(ph$PHAFSTMN));first_ge8<-sum(ph$PHAFSTHR+ph$PHAFSTMN/60>=8,na.rm=TRUE)
inventory<-list()
for(f in c('L11_2_B','L25_2_B','L40_2_B')){z<-read_xpt(paste0('work/source/',f,'.xpt'));inventory[[f]]<-names(z);stopifnot(!any(grepl('FAST|FST|^PH',names(z),ignore.case=TRUE)))}
fast<-expand.grid(score=c('log1p','ln'),group=c('both_ge8h','both_lt8h','discordant_8h','unknown_pair_fasting'),metric=c('ICC','MDC95'),stringsAsFactors=FALSE)
fast$n_replication<-n;fast$n_stratum<-NA_integer_;fast$value<-NA_real_;fast$ci_lower<-NA_real_;fast$ci_upper<-NA_real_;fast$status<-'NOT_RUN';fast$bootstrap_replicates<-0L
fast$reason<-'Only first-exam PH_B PHAFSTHR/PHAFSTMN found; no second-exam food-fasting duration in complete laboratory variable inventory, L11_2_B, L25_2_B or L40_2_B. LB2DAY is inter-exam days, not fasting duration. Pair strata cannot be assigned.'
fast$files_checked<-'PH_B.htm/xpt; L11_2_B.htm/xpt; L25_2_B.htm/xpt; L40_2_B.htm/xpt; full 2001-2002 Laboratory Variable List; general second-day PDF retrieval attempted (404)'
fast$variables_checked<-'PHAFSTHR; PHAFSTMN; PHDSESN; LB2DAY; all repeat-file names searched FAST/FST/PH; complete inventory labels searched fasting/food fast/phlebotomy'
fast$first_exam_fasting_duration_available_n<-first_known;fast$first_exam_ge8h_n<-first_ge8
write.csv(fast,file.path(O,'S6d_fasting_replication.csv'),row.names=FALSE,na='')
write_json(list(first_fasting_known_n=first_known,first_ge8h_n=first_ge8,second_exam_variable_names=inventory,not_run_reason=unique(fast$reason)),file.path(O,'code/S6_fasting_inventory.json'),auto_unbox=TRUE,pretty=TRUE)
qa<-list(n_main=nrow(a$main),n_replication=n,bootstrap_replicates=Rboot,residual_log1p=coef(b$resmod),residual_ln=coef(s$model_ln),component_sum_max_errors=identity_errors,original_B_point_max_error=old_err,interval_tertile_cutpoints=q,interval_tertile_n=as.list(table(group)),creatinine_extreme_n=sum(big),S6b_tukey_counts=as.list(setNames(cvres$n_removed[cvres$screen=='after_Tukey_3IQR'],refs$analyte)),S6d_status='NOT_RUN',session_info=capture.output(sessionInfo()))
write_json(qa,file.path(O,'S6_QA.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
saveRDS(list(components=comp,component_wide=wide,cvw=cvres,pairs=pairs,fasting=fast,qa=qa), 'work/S6/results.rds')
print(comp[comp$component%in%c('creat_umol','glucose_mmol'),c('cohort','score','component','delta_contribution_rank','MDC95_fixed','MDC_drop_percent','difference_variance_drop_percent')],row.names=FALSE)
print(cvres[,c('analyte','screen','n','CVt_percent','CVw_percent','difference_CVw_percentage_points','classification')],row.names=FALSE)
print(pairs[pairs$metric%in%c('mean_difference','SD_difference','Deming_intercept','Deming_slope','n_abs_difference_gt0_3','Spearman_rho','ICC','MDC95'),c('group','score','metric','n','value','ci_lower','ci_upper')],row.names=FALSE)
