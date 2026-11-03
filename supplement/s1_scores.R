suppressPackageStartupMessages({library(haven);library(jsonlite)})
O<-'outputs/BIOAGE-REL-01/S';dir.create('work/S',showWarnings=FALSE)
st<-fromJSON(file.path(O,'S_status.json'));stopifnot(st$S5$classifications$EXACT==0,st$S5$classifications$PARENT==0)
a<-readRDS('work/a1.rds');cal<-readRDS('work/a2.rds');b<-readRDS('work/b_results.rds');d<-readRDS('work/d_results.rds')
pheno<-cal$pheno;project_input<-cal$project_input;metrics<-b$metrics;icc21<-b$icc21;Rboot<-2000L;set.seed(202611021)
stopifnot(nrow(a$main)==1801,sum(a$main$MORTSTAT==1)==829,nrow(d$d)==407,b$base_n==15635)
bm<-c('albumin','alp','creat','glucose','crp','wbc','lymph','mcv','rdw','age')
base<-cal$p[cal$p$age>=20&complete.cases(cal$p[,bm]),];stopifnot(nrow(base)==15635,all(base$crp>0))
base$lncrp<-log(base$crp);base$phenoage_ln<-pheno(base)
model_ln<-lm(phenoage_ln~age,data=base);saveRDS(model_ln,'work/S/ln_residual_model.rds')
write_json(list(label='写作期预定补充分析',score='v1.0 coefficients; natural-log CRP only',n=15635,intercept=unname(coef(model_ln)[1]),age_slope=unname(coef(model_ln)[2]),v1_intercept=unname(coef(b$resmod)[1]),v1_age_slope=unname(coef(b$resmod)[2]),refit_during_bootstrap=FALSE),file.path(O,'S1_ln_residual_model.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
z1<-b$z1;z2<-b$z2
rd<-d$d;first<-c('LBXSAL','LBDSAPSI','LBDSCR','LBXSGL','LBXCRP','LBXWBCSI','LBXLYPCT','LBXMCVSI','LBXRDW');second<-c('LB2SAL','LB2SAPSI','LB2SCR','LB2SGL','LB2CRP','LB2WBCSI','LB2LYPCT','LB2MCVSI','LB2RDW')
project_rep<-function(fields){z<-as.data.frame(rd[,fields]);names(z)<-bm[1:9];z$albumin_gL<-z$albumin*10;z$glucose_mmol<-z$glucose*.0555;z$creat_umol<-z$creat*88.4017;z$lncrp<-log1p(z$crp);z$age<-rd$RIDAGEYR;z$SEQN<-rd$SEQN;z}
r1<-project_rep(first);r2<-project_rep(second)
accel<-function(z,version){if(version=='ln'){z$lncrp<-log(z$crp);pheno(z)-predict(model_ln,z)}else pheno(z)-predict(b$resmod,z)}
cohorts<-list(NHANES_III=list(z1=z1,z2=z2),NHANES_2001_2002=list(z1=r1,z2=r2));scores<-list();rows<-list();k<-1L
for(cohort in names(cohorts)){
 zz<-cohorts[[cohort]];n<-nrow(zz$z1);boot_idx<-replicate(Rboot,sample.int(n,n,replace=TRUE))
 scores[[cohort]]<-list()
 for(version in c('log1p','ln')){
  x1<-accel(zz$z1,version);x2<-accel(zz$z2,version);f1<-zz$z1;f2<-zz$z2;crpmean<-(f1$crp+f2$crp)/2;f1$crp<-f2$crp<-crpmean
  f1$lncrp<-f2$lncrp<-if(version=='ln')log(crpmean)else log1p(crpmean)
  m<-if(version=='ln')model_ln else b$resmod;xf1<-pheno(f1)-predict(m,f1);xf2<-pheno(f2)-predict(m,f2)
  calc<-function(ii){m0<-metrics(x1[ii],x2[ii]);mdcf<-1.96*sd(xf2[ii]-xf1[ii]);c(m0,CRP_fixed_MDC95=mdcf,CRP_fixed_MDC_drop_years=unname(m0['MDC95'])-mdcf,CRP_fixed_MDC_drop_percent=100*(1-mdcf/unname(m0['MDC95'])))}
  obs<-calc(seq_len(n));boot<-apply(boot_idx,2,calc);ci<-t(apply(boot,1,quantile,c(.025,.975),na.rm=TRUE))
  rows[[k]]<-data.frame(label='写作期预定补充分析',supplement='S1',cohort=cohort,score=version,metric=names(obs),n=n,value=as.numeric(obs),ci_lower=ci[,1],ci_upper=ci[,2],unit=ifelse(names(obs)%in%c('ICC','top20_retention'),'proportion',ifelse(names(obs)=='CRP_fixed_MDC_drop_percent','percent','years')),bootstrap_replicates=Rboot,status='COMPUTED');k<-k+1L
  scores[[cohort]][[version]]<-list(x1=x1,x2=x2,fixed1=xf1,fixed2=xf2,boot=boot)
 }
}
res<-do.call(rbind,rows);write.csv(res,file.path(O,'S1_lncrp_reliability.csv'),row.names=FALSE,na='')
stopifnot(max(abs(scores$NHANES_III$log1p$x1-b$y1))<1e-10,max(abs(scores$NHANES_III$log1p$x2-b$y2))<1e-10,max(abs(scores$NHANES_2001_2002$log1p$x1-d$x1))<1e-10,max(abs(scores$NHANES_2001_2002$log1p$x2-d$x2))<1e-10)
# Plausibility: full 2001–2002 first-exam adult complete-case laboratory sample.
demo<-as.data.frame(read_xpt('work/source/DEMO_B.xpt'))[,c('SEQN','RIDAGEYR')]
labs<-list(as.data.frame(read_xpt('work/source/L11_B.xpt'))[,c('SEQN','LBXCRP')],as.data.frame(read_xpt('work/source/L25_B.xpt'))[,c('SEQN','LBXWBCSI','LBXLYPCT','LBXMCVSI','LBXRDW')],as.data.frame(read_xpt('work/source/L40_B.xpt'))[,c('SEQN','LBXSAL','LBDSAPSI','LBDSCR','LBXSGL')])
full<-Reduce(function(x,y)merge(x,y,by='SEQN'),c(list(demo),labs));full<-full[full$RIDAGEYR>=20&complete.cases(full[,first]),]
z<-as.data.frame(full[,first]);names(z)<-bm[1:9];z$albumin_gL<-z$albumin*10;z$glucose_mmol<-z$glucose*.0555;z$creat_umol<-z$creat*88.4017;z$age<-full$RIDAGEYR
pl<-list();k<-1L
for(version in c('ln','log1p')){
 z$lncrp<-if(version=='ln')log(z$crp)else log1p(z$crp);age<-pheno(z);n<-length(age);obs<-c(mean_PhenoAge=mean(age),SD_PhenoAge=sd(age));bt<-replicate(Rboot,{ii<-sample.int(n,n,TRUE);c(mean(age[ii]),sd(age[ii]))});ci<-t(apply(bt,1,quantile,c(.025,.975)))
 pl[[k]]<-data.frame(label='写作期预定补充分析',score=version,cohort='NHANES_2001_2002_first_exam_all_adults_complete',n=n,metric=names(obs),value=as.numeric(obs),ci_lower=ci[,1],ci_upper=ci[,2],Liu2018_reference_N=11432,Liu2018_reference_value=NA_real_,difference_from_Liu=NA_real_,local_status='COMPUTED',comparison_status='NOT_RUN',reason='Liu 2018 main text and S1 Table do not report raw Phenotypic Age mean or SD; S1 Table lists chronological age and individual biomarkers. The standardized acceleration SD=1 in Figure 3 is not raw-score SD.',source_doi='10.1371/journal.pmed.1002718',source_table_url='https://journals.plos.org/plosmedicine/article/file?id=10.1371/journal.pmed.1002718.s005&type=supplementary');k<-k+1L
}
write.csv(do.call(rbind,pl),file.path(O,'S1_plausibility.csv'),row.names=FALSE,na='')
saveRDS(list(scores=scores,cohorts=cohorts,model_ln=model_ln,reliability=res,full_first_adult_N=nrow(full)), 'work/S/s1_results.rds')
print(res[res$metric%in%c('ICC','MDC95','CRP_fixed_MDC_drop_percent'),],row.names=FALSE);cat('ln residual coefficients:',coef(model_ln),'\n')
