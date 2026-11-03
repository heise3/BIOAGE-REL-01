suppressPackageStartupMessages({library(haven);library(readr);library(dplyr);library(jsonlite)})
out <- 'outputs/BIOAGE-REL-01'
p <- read.csv('work/lab_selected.csv'); r <- read.csv('work/labse_selected.csv'); e <- read.csv('work/examse_selected.csv')
stopifnot(!anyDuplicated(p$SEQN),!anyDuplicated(r$SEQN),!anyDuplicated(e$SEQN))
bm1 <- c('AMP','APPSI','CEP','SGP','CRP','WCP','LMPPCNT','MVPSI','RWP')
bm2 <- c('AMR','APRSI','CER','SGR','CRR','WCR','LMRPCNT','MVRSI','RWR')
d <- merge(merge(r,p,by='SEQN',all.x=TRUE),e,by='SEQN',all.x=TRUE)
complete <- complete.cases(d[,c(bm1,bm2)])
adult <- !is.na(d$HSAGEIR)&d$HSAGEIR>=20
pregnant <- !is.na(d$MARF12R)&d$MARF12R==1
main <- d[adult&complete&!pregnant,]
m <- read_fwf('work/source/NHANES_III_MORT_2019_PUBLIC.dat',fwf_cols(SEQN=c(1,6),ELIGSTAT=c(15,15),MORTSTAT=c(16,16),PERMTH_INT=c(43,45),PERMTH_EXM=c(46,48)),col_types='ddddd',na=c('','.'))
main <- merge(main,m,by='SEQN',all.x=TRUE)
cnt <- function(metric,value,denom=NA_real_,note='') data.frame(metric,value,denominator=denom,status='COMPUTED',note)
counts <- rbind(cnt('N_labse',nrow(r)),cnt('N_exam_se',nrow(e)),cnt('N_adult_second',sum(adult)),cnt('N_adult_18_second',sum(d$HSAGEIR>=18,na.rm=TRUE)),cnt('N_complete_adult_before_pregnancy',sum(adult&complete)),cnt('N_pregnant_complete',sum(adult&complete&pregnant)),cnt('N_main',nrow(main)),cnt('N_main_pregnancy_unknown_women20_59',sum(main$HSSEX==2&main$HSAGEIR<60 & (is.na(main$MARF12R)|!main$MARF12R%in%c(1,2)))),cnt('N_mort_eligible',sum(main$ELIGSTAT==1,na.rm=TRUE),nrow(main)),cnt('D',sum(main$ELIGSTAT==1&main$MORTSTAT==1,na.rm=TRUE)),cnt('N_same_session',sum(main$MXPSESSR==main$MXRSESSR,na.rm=TRUE)),cnt('N_sessions_both_known',sum(complete.cases(main[,c('MXPSESSR','MXRSESSR')]))),cnt('N_interval_known',sum(!is.na(main$MXPRDAYS))))
x1 <- read_xpt('work/source/L11_2_B.xpt');x2 <- read_xpt('work/source/L25_2_B.xpt');x3 <- read_xpt('work/source/L40_2_B.xpt');demo <- read_xpt('work/source/DEMO_B.xpt')
stopifnot(!anyDuplicated(x1$SEQN),!anyDuplicated(x2$SEQN),!anyDuplicated(x3$SEQN))
repbm <- c('LB2CRP','LB2WBCSI','LB2LYPCT','LB2MCVSI','LB2RDW','LB2SAL','LB2SAPSI','LB2SCR','LB2SGL')
repd <- Reduce(function(a,b) merge(a,b,by='SEQN'),list(x1,x2[,setdiff(names(x2),'LB2DAY')],x3[,setdiff(names(x3),'LB2DAY')],demo[,c('SEQN','RIDAGEYR','RIAGENDR','RIDEXPRG')]))
repn_all <- sum(complete.cases(repd[,repbm])); repd <- repd[complete.cases(repd[,repbm])&repd$RIDAGEYR>=20 & !(repd$RIDEXPRG%in%1),]
counts <- rbind(counts,cnt('N_rep_all_ages_complete',repn_all),cnt('N_rep',nrow(repd),NA,'Age >=20; exclude RIDEXPRG=1; complete nine second-exam markers'))
write.csv(counts,file.path(out,'A1_counts.csv'),row.names=FALSE,na='')
saveRDS(list(p=p,r=r,e=e,d=d,main=main,repd=repd,bm1=bm1,bm2=bm2), 'work/a1.rds')
st<-fromJSON(file.path(out,'status.json'));a1pass<-nrow(main)>=600&&sum(main$ELIGSTAT==1&main$MORTSTAT==1,na.rm=TRUE)>=300
st$A1<-list(verdict=ifelse(a1pass,'PASS','STOP'),N_main=nrow(main),D=sum(main$ELIGSTAT==1&main$MORTSTAT==1,na.rm=TRUE),N_rep=nrow(repd),N_rep_all_target_ages_complete=repn_all)
st$layer<-'A1';st$verdict<-ifelse(a1pass,'PASS','STOP');st$reason<-ifelse(a1pass,'A1 schema and counts passed; A2 pending.','STOP_DATA: A1 count gate not passed.');write_json(st,file.path(out,'status.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
stopifnot(a1pass)
print(counts,row.names=FALSE)
load('work/source/NHANES3.rda');ref <- NHANES3;ref$SEQN <- as.numeric(sub('.*_','',ref$sampleID));z <- merge(p,ref,by='SEQN');
cat('Raw versus BioAge preprocessing diagnostics:\n')
for(pair in list(c('AMP','albumin'),c('CEP','creat'),c('SGP','glucose'),c('LMPPCNT','lymph'),c('MVPSI','mcv'),c('RWP','rdw'),c('APPSI','alp'),c('WCP','wbc'),c('CRP','crp'))){a=z[[pair[1]]];b=z[[pair[2]]];ok=is.finite(a)&is.finite(b);cat(pair,paste('n',sum(ok),'max_abs',max(abs(a[ok]-b[ok])),'mean_abs',mean(abs(a[ok]-b[ok]))),'\n')}
cat('Calibration input preprocessing identities: lncrp log1p error=',max(abs(ref$lncrp-log1p(ref$crp)),na.rm=TRUE),'\n')
