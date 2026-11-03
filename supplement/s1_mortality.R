suppressPackageStartupMessages({library(survival);library(readr);library(haven);library(jsonlite)})
O<-'outputs/BIOAGE-REL-01/S';a<-readRDS('work/a1.rds');b<-readRDS('work/b_results.rds');d<-readRDS('work/d_results.rds');s1<-readRDS('work/S/s1_results.rds');old<-read.csv('outputs/BIOAGE-REL-01/C_mortality.csv');set.seed(202611023);Rboot<-2000L
main<-data.frame(SEQN=a$main$SEQN,age=a$main$HSAGEIR,sex=factor(a$main$HSSEX),race=factor(a$main$DMARETHN),ELIGSTAT=a$main$ELIGSTAT,MORTSTAT=a$main$MORTSTAT,PERMTH_EXM=a$main$PERMTH_EXM)
repdata<-data.frame(SEQN=d$d$SEQN,age=d$d$RIDAGEYR,sex=factor(d$d$RIAGENDR));rep_reason<-''
path<-'work/S_sources/NHANES_2001_2002_MORT_2019_PUBLIC.dat'
if(file.exists(path)){
 m<-as.data.frame(read_fwf(path,fwf_cols(SEQN=c(1,6),ELIGSTAT=c(15,15),MORTSTAT=c(16,16),PERMTH_EXM=c(46,48)),col_types='dddd',na=c('','.')))
 de<-as.data.frame(read_xpt('work/source/DEMO_B.xpt'))[,c('SEQN','RIDRETH1')]
 mm<-m[match(repdata$SEQN,m$SEQN),];race0<-de$RIDRETH1[match(repdata$SEQN,de$SEQN)];repdata$race<-factor(ifelse(race0==3,1,ifelse(race0==4,2,ifelse(race0==1,3,ifelse(race0%in%c(2,5),4,NA)))))
 repdata$ELIGSTAT<-mm$ELIGSTAT;repdata$MORTSTAT<-mm$MORTSTAT;repdata$PERMTH_EXM<-mm$PERMTH_EXM
}else rep_reason<-'Public 2001–2002 mortality file unavailable; fixed 407-person replication sample retained.'
cohorts<-list(NHANES_III=main,NHANES_2001_2002=repdata);rows<-list();phrows<-list();k<-1L;j<-1L;boots<-list()
for(cohort in names(cohorts)){
 dd<-cohorts[[cohort]];n<-nrow(dd);scorecohort<-s1$scores[[cohort]];scale_sd<-sd(scorecohort$log1p$x1)
 set.seed(if(cohort=='NHANES_III')20261103L else 202611024L)
 bootstrap_indices<-replicate(Rboot,sample.int(n,n,TRUE))
 executable<-all(c('ELIGSTAT','MORTSTAT','PERMTH_EXM','race')%in%names(dd))
 if(executable)executable<-all(complete.cases(dd))&&all(dd$ELIGSTAT==1)&&all(dd$PERMTH_EXM>0)
 for(version in c('log1p','ln')){
  if(!executable){rows[[k]]<-data.frame(label='写作期预定补充分析',cohort=cohort,score=version,exposure=c('exam1','exam2_primary','mean'),n=n,deaths=NA,HR=NA,ci_lower=NA,ci_upper=NA,common_SD_years=scale_sd,bootstrap_requested=Rboot,bootstrap_success=0,status='NOT_RUN',reason=ifelse(nchar(rep_reason)>0,rep_reason,'At least one member of the fixed cohort lacks eligible mortality linkage, time, or covariates; no new subset selected.'),source='S1');k<-k+1L;next}
  if(cohort=='NHANES_III'&&version=='log1p'){
   rr<-old[match(c('HR_exam1','HR_exam2','HR_mean'),old$metric),];stopifnot(all(rr$n==1801),all(rr$deaths==829))
   rows[[k]]<-data.frame(label='写作期预定补充分析',cohort=cohort,score=version,exposure=c('exam1','exam2_primary','mean'),n=n,deaths=sum(dd$MORTSTAT==1),HR=rr$value,ci_lower=rr$ci_lower,ci_upper=rr$ci_upper,common_SD_years=scale_sd,bootstrap_requested=Rboot,bootstrap_success=rr$bootstrap_success,status='REUSED_FROZEN_V1_0',reason='',source='v1.0 C_mortality.csv');k<-k+1L;next
  }
  sc<-scorecohort[[version]];dd$x1<-sc$x1;dd$x2<-sc$x2;dd$xm<-(sc$x1+sc$x2)/2
  calc<-function(x,keep=FALSE){fits<-lapply(c('x1','x2','xm'),function(v){x$exposure<-x[[v]]/scale_sd;coxph(Surv(PERMTH_EXM,MORTSTAT)~exposure+age+sex+race,data=x,ties='efron',x=keep)});vals<-sapply(fits,function(f)coef(f)['exposure']);if(keep)list(values=vals,fits=fits)else vals}
  obs<-calc(dd,TRUE);bt<-matrix(NA_real_,Rboot,3)
  for(i in seq_len(Rboot)){
   ii<-bootstrap_indices[,i];bt[i,]<-tryCatch(calc(dd[ii,]),error=function(e)rep(NA_real_,3))
   if(i%%500==0){cat(cohort,version,'Cox bootstrap',i,'/2000\n');flush.console()}
  }
  ci<-t(apply(bt,2,quantile,c(.025,.975),na.rm=TRUE));rows[[k]]<-data.frame(label='写作期预定补充分析',cohort=cohort,score=version,exposure=c('exam1','exam2_primary','mean'),n=n,deaths=sum(dd$MORTSTAT==1),HR=exp(obs$values),ci_lower=exp(ci[,1]),ci_upper=exp(ci[,2]),common_SD_years=scale_sd,bootstrap_requested=Rboot,bootstrap_success=colSums(is.finite(bt)),status='COMPUTED',reason='',source='S1 fixed cohort; public 2019 LMF; age+sex+four-category race/ethnicity')
  k<-k+1L;boots[[paste(cohort,version)]]<-bt
  for(v in 1:3){z<-cox.zph(obs$fits[[v]]);phrows[[j]]<-data.frame(cohort=cohort,score=version,exposure=c('exam1','exam2_primary','mean')[v],term=rownames(z$table),z$table,n=n);j<-j+1L}
 }
}
res<-do.call(rbind,rows);rownames(res)<-NULL;write.csv(res,file.path(O,'S1_lncrp_mortality.csv'),row.names=FALSE,na='')
if(length(phrows))write.csv(do.call(rbind,phrows),file.path(O,'S1_PH_tests.csv'),row.names=FALSE,na='')
saveRDS(list(results=res,bootstrap=boots),'work/S/s1_mortality.rds');print(res,row.names=FALSE)
