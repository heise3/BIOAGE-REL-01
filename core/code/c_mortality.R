suppressPackageStartupMessages({library(survival);library(jsonlite)})
a<-readRDS('work/a1.rds');b<-readRDS('work/b_results.rds');out<-'outputs/BIOAGE-REL-01';st<-fromJSON(file.path(out,'status.json'));stopifnot(st$B$verdict=='PASS')
set.seed(20261103);Rboot<-2000L;d<-a$main;d$x1<-b$y1;d$x2<-b$y2;d$xm<-(d$x1+d$x2)/2
d$age<-d$HSAGEIR;d$sex<-factor(d$HSSEX);d$race<-factor(d$DMARETHN)
ok<-d$ELIGSTAT==1&is.finite(d$PERMTH_EXM)&d$PERMTH_EXM>0&complete.cases(d[,c('x1','x2','age','sex','race','MORTSTAT')]);d<-d[ok,];scale_sd<-sd(d$x1)
calc<-function(dd,savefits=FALSE){
 fits<-lapply(c('x1','x2','xm'),function(v){dd$x<-dd[[v]]/scale_sd;coxph(Surv(PERMTH_EXM,MORTSTAT)~x+age+sex+race,data=dd,ties='efron',x=savefits,model=savefits)})
 beta<-sapply(fits,function(f)coef(f)['x']);cc<-sapply(fits,function(f)summary(f)$concordance['C']);lam<-cov(dd$x1,dd$x2)/var(dd$x2)
 vals<-c(logHR_exam1=unname(beta[1]),logHR_exam2=unname(beta[2]),logHR_mean=unname(beta[3]),lambda=lam,logHR_exam2_corrected=unname(beta[2])/lam,C_exam1=unname(cc[1]),C_exam2=unname(cc[2]),C_mean=unname(cc[3]),delta_C_mean_minus_exam1=unname(cc[3]-cc[1]),delta_C_mean_minus_exam2=unname(cc[3]-cc[2]))
 if(savefits)list(values=vals,fits=fits) else vals
}
obs<-calc(d,TRUE);bt<-matrix(NA_real_,Rboot,length(obs$values),dimnames=list(NULL,names(obs$values)))
for(i in seq_len(Rboot)){
 ii<-sample.int(nrow(d),nrow(d),TRUE);bt[i,]<-tryCatch(calc(d[ii,]),error=function(e)rep(NA_real_,ncol(bt)))
 if(i%%200==0){cat('C bootstrap',i,'of',Rboot,'\n');flush.console()}
}
ci<-t(apply(bt,2,quantile,c(.025,.975),na.rm=TRUE));res<-data.frame(metric=names(obs$values),n=nrow(d),deaths=sum(d$MORTSTAT==1),value=unname(obs$values),ci_lower=ci[,1],ci_upper=ci[,2],unit='index',status='COMPUTED',bootstrap_success=colSums(is.finite(bt)))
hr<-grepl('logHR',res$metric);res$metric[hr]<-sub('logHR','HR',res$metric[hr]);res$value[hr]<-exp(res$value[hr]);res$ci_lower[hr]<-exp(res$ci_lower[hr]);res$ci_upper[hr]<-exp(res$ci_upper[hr]);res$unit[hr]<-'HR per common exam1 SD';res$common_SD_years<-scale_sd
write.csv(res,file.path(out,'C_mortality.csv'),row.names=FALSE,na='')
ph<-do.call(rbind,lapply(seq_along(obs$fits),function(i){z<-cox.zph(obs$fits[[i]]);data.frame(exposure=c('exam1','exam2_primary','mean')[i],term=rownames(z$table),z$table,n=nrow(d),status='COMPUTED')}))
write.csv(ph,file.path(out,'C_PH_tests.csv'),row.names=FALSE,na='')
saveRDS(list(results=res,ph=ph,boot=bt,common_SD_years=scale_sd,fits=obs$fits),'work/c_results.rds')
st$C<-list(verdict='PASS',status='COMPUTED',n=nrow(d),deaths=sum(d$MORTSTAT==1),bootstrap_replicates=Rboot,bootstrap_success_min=min(res$bootstrap_success));st$layer<-'C';st$verdict<-'PASS';st$reason<-'A1, A2 and B precision passed; C complete; replication follows.';write_json(st,file.path(out,'status.json'),pretty=T,auto_unbox=T)
print(res,row.names=FALSE);print(ph,row.names=FALSE)
