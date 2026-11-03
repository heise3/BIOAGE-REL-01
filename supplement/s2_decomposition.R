suppressPackageStartupMessages(library(jsonlite))
O<-'outputs/BIOAGE-REL-01/S';a<-readRDS('work/a1.rds');b<-readRDS('work/b_results.rds');d<-readRDS('work/d_results.rds');s1<-readRDS('work/S/s1_results.rds');cal<-readRDS('work/a2.rds');pheno<-cal$pheno;icc21<-b$icc21
set.seed(202611022);Rboot<-2000L
calc<-function(x,y){sw<-sd(y-x)/sqrt(2);v1<-var(x)-sw^2;v2<-var(y)-sw^2;c(mean1=mean(x),mean2=mean(y),totalSD1=sd(x),totalSD2=sd(y),betweenSD1=if(v1>=0)sqrt(v1)else NA_real_,betweenSD2=if(v2>=0)sqrt(v2)else NA_real_,between_variance1=v1,between_variance2=v2,withinSD=sw,ICC=icc21(x,y),MDC95=1.96*sd(y-x),pair_mean_SD=sd((x+y)/2))}
analyses<-list();add<-function(cohort,scope,x,y){analyses[[length(analyses)+1]]<<-list(cohort=cohort,scope=scope,x=x,y=y)}
add('NHANES_III','v1_full',b$y1,b$y2);add('NHANES_2001_2002','v1_full',d$x1,d$x2)
ii<-a$main$HSAGEIR>=20&a$main$HSAGEIR<=69;add('NHANES_III','age_20_69',b$y1[ii],b$y2[ii]);ii<-d$d$RIDAGEYR>=20&d$d$RIDAGEYR<=69;add('NHANES_2001_2002','age_20_69',d$x1[ii],d$x2[ii])
z1<-s1$cohorts$NHANES_2001_2002$z1;z2<-s1$cohorts$NHANES_2001_2002$z2;avg<-(z1$creat_umol+z2$creat_umol)/2;z1$creat_umol<-z2$creat_umol<-avg
fx1<-pheno(z1)-predict(b$resmod,z1);fx2<-pheno(z2)-predict(b$resmod,z2);add('NHANES_2001_2002','creatinine_pair_mean',fx1,fx2)
stopifnot(abs(sd((fx1+fx2)/2)-sd((d$x1+d$x2)/2))<1e-10)
rows<-list();k<-1L
for(g in analyses){n<-length(g$x);obs<-calc(g$x,g$y);bt<-replicate(Rboot,{ii<-sample.int(n,n,TRUE);calc(g$x[ii],g$y[ii])});ci<-t(apply(bt,1,quantile,c(.025,.975),na.rm=TRUE))
 for(visit in 1:2){
  suffix<-as.character(visit)
  rows[[k]]<-data.frame(label='写作期预定补充分析',score='v1.0_log1p',cohort=g$cohort,scope=g$scope,visit=visit,n=n,mean_accel=obs[paste0('mean',suffix)],mean_ci_lower=ci[paste0('mean',suffix),1],mean_ci_upper=ci[paste0('mean',suffix),2],mean_CI_excludes_zero=ci[paste0('mean',suffix),1]>0||ci[paste0('mean',suffix),2]<0,total_SD=obs[paste0('totalSD',suffix)],total_SD_ci_lower=ci[paste0('totalSD',suffix),1],total_SD_ci_upper=ci[paste0('totalSD',suffix),2],between_SD=obs[paste0('betweenSD',suffix)],between_SD_ci_lower=ci[paste0('betweenSD',suffix),1],between_SD_ci_upper=ci[paste0('betweenSD',suffix),2],between_variance_raw=obs[paste0('between_variance',suffix)],within_SD=obs['withinSD'],within_SD_ci_lower=ci['withinSD',1],within_SD_ci_upper=ci['withinSD',2],ICC=obs['ICC'],ICC_ci_lower=ci['ICC',1],ICC_ci_upper=ci['ICC',2],MDC95=obs['MDC95'],MDC95_ci_lower=ci['MDC95',1],MDC95_ci_upper=ci['MDC95',2],pair_mean_SD=obs['pair_mean_SD'],bootstrap_replicates=Rboot,between_SD_bootstrap_valid=sum(is.finite(bt[paste0('betweenSD',suffix),])),status=ifelse(is.finite(obs[paste0('betweenSD',suffix)]),'COMPUTED','NOT_RUN_NEGATIVE_BETWEEN_VARIANCE'),note=ifelse(g$scope=='creatinine_pair_mean','Fixing paired creatinine mean removes within-person fluctuations while retaining between-person level differences; paired-score mean SD is invariant.','Descriptive only; fixed v1.0 residual model; between SD=sqrt(total variance - within variance)'))
  k<-k+1L
 }
}
res<-do.call(rbind,rows);rownames(res)<-NULL;write.csv(res,file.path(O,'S2_icc_decomposition.csv'),row.names=FALSE,na='');saveRDS(list(results=res,creatinine_fixed_x1=fx1,creatinine_fixed_x2=fx2),'work/S/s2_results.rds');print(res[,c('cohort','scope','visit','n','mean_accel','total_SD','between_SD','within_SD','ICC','MDC95')],row.names=FALSE)
