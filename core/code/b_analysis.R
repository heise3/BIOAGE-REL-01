suppressPackageStartupMessages({library(dplyr);library(jsonlite)})
a<-readRDS('work/a1.rds');b<-readRDS('work/a2.rds');stopifnot(all(b$cal$status=='PASS'))
pheno<-b$pheno;project_input<-b$project_input;out<-'outputs/BIOAGE-REL-01';set.seed(20261102);Rboot<-2000L
base<-b$p[b$p$age>=20&complete.cases(b$p[,c('albumin','alp','creat','glucose','crp','wbc','lymph','mcv','rdw','age')]),]
resmod<-lm(phenoage~age,data=base);z1<-project_input(a$main);z2<-project_input(a$main,TRUE)
y1<-pheno(z1)-predict(resmod,z1);y2<-pheno(z2)-predict(resmod,z2)
icc21<-function(x,y){n<-length(x);rr<-(x+y)/2;gm<-mean(rr);msr<-2*sum((rr-gm)^2)/(n-1);msc<-n*((mean(x)-gm)^2+(mean(y)-gm)^2);mse<-sum(((x-y)-mean(x-y))^2)/(2*(n-1));(msr-mse)/(msr+mse+2*(msc-mse)/n)}
metrics<-function(x,y){dd<-y-x;sdiff<-sd(dd);q1<-quantile(x,.8,names=FALSE);q2<-quantile(y,.8,names=FALSE);c(ICC=icc21(x,y),within_SD=sdiff/sqrt(2),MDC95=1.96*sdiff,BA_bias=mean(dd),BA_lower=mean(dd)-1.96*sdiff,BA_upper=mean(dd)+1.96*sdiff,top20_retention=mean(y[x>=q1]>=q2))}
reliability<-function(x,y,label='PhenoAgeAccel',group='all',unit='years',status='COMPUTED'){
 ok<-is.finite(x)&is.finite(y);x<-x[ok];y<-y[ok];n<-length(x)
 if(n<3){mn<-names(metrics(1:4,2:5));return(data.frame(algorithm=label,group=group,metric=mn,n=n,value=NA,ci_lower=NA,ci_upper=NA,unit=ifelse(mn%in%c('ICC','top20_retention'),'proportion',unit),status='NOT_RUN',reason='Fewer than 3 paired people; stratum contrast is not identifiable'))}
 obs<-metrics(x,y);boot<-replicate(Rboot,{ii<-sample.int(n,n,replace=TRUE);metrics(x[ii],y[ii])});ci<-t(apply(boot,1,quantile,c(.025,.975),na.rm=TRUE))
 data.frame(algorithm=label,group=group,metric=names(obs),n=n,value=as.numeric(obs),ci_lower=ci[,1],ci_upper=ci[,2],unit=ifelse(names(obs)%in%c('ICC','top20_retention'),'proportion',unit),status=status,reason='')
}
rel<-reliability(y1,y2)
cat('Primary reliability computed\n');flush.console()
for(g in c('20-39','40-59','60+','male','female')){
 ii<-switch(g,'20-39'=which(a$main$HSAGEIR<40),'40-59'=which(a$main$HSAGEIR>=40&a$main$HSAGEIR<60),'60+'=which(a$main$HSAGEIR>=60),'male'=which(a$main$HSSEX==1),'female'=which(a$main$HSSEX==2))
 rel<-rbind(rel,reliability(y1[ii],y2[ii],group=g))
}
write.csv(rel,file.path(out,'B_reliability.csv'),row.names=FALSE,na='')
sens<-reliability(y1[a$main$CRP<1&a$main$CRR<1],y2[a$main$CRP<1&a$main$CRR<1],group='CRP_both_below_1mgdL')
fast1<-a$main$PHPFAST;fast2<-a$main$PHRFAST
groups<-list(session_same=which(a$main$MXPSESSR==a$main$MXRSESSR),session_different=which(a$main$MXPSESSR!=a$main$MXRSESSR),fast_both_ge8h=which(fast1>=8&fast2>=8),fast_both_lt8h=which(fast1<8&fast2<8),fast_discordant_8h=which((fast1>=8)!=(fast2>=8)),fast_unknown=which(is.na(fast1)|is.na(fast2)),interval_known=which(!is.na(a$main$MXPRDAYS)),interval_unknown=which(is.na(a$main$MXPRDAYS)))
for(g in names(groups)){ii<-groups[[g]];sens<-rbind(sens,reliability(y1[ii],y2[ii],group=g))}
# Fixed-component counterfactuals use identical raw arithmetic means at both visits.
bioms<-c('albumin_gL','alp','creat_umol','glucose_mmol','lncrp','wbc','lymph','mcv','rdw')
coefs<-c(-.03359355,.001868778,.009506491,.1953192,.09536762,.05542406,-.01199984,.02676401,.3306156)/.090165
dif<-y2-y1;components<-list();mdc0<-1.96*sd(dif)
for(j in seq_along(bioms)){
 bm<-bioms[j];contrib<-coefs[j]*(z2[[bm]]-z1[[bm]]);fixed<-dif-contrib
 obs<-c(MDC95_fixed=1.96*sd(fixed),MDC_drop_years=mdc0-1.96*sd(fixed),MDC_drop_percent=100*(1-sd(fixed)/sd(dif)),coefficient_times_within_SD=abs(coefs[j])*sd(z2[[bm]]-z1[[bm]])/sqrt(2),covariance_allocated_variance=cov(contrib,dif)/2)
 bt<-replicate(Rboot,{ii<-sample.int(length(dif),length(dif),TRUE);dd<-dif[ii];ff<-fixed[ii];cc<-contrib[ii];c(1.96*sd(ff),1.96*(sd(dd)-sd(ff)),100*(1-sd(ff)/sd(dd)),sd(cc)/sqrt(2),cov(cc,dd)/2)})
 ci<-t(apply(bt,1,quantile,c(.025,.975)));components[[j]]<-data.frame(component=bm,metric=names(obs),n=length(dif),value=as.numeric(obs),ci_lower=ci[,1],ci_upper=ci[,2],status='COMPUTED')
 if(bm=='glucose_mmol'){
  q1<-z1;q2<-z2;avg<-(q1$glucose_mmol+q2$glucose_mmol)/2;q1$glucose_mmol<-avg;q2$glucose_mmol<-avg
  sens<-rbind(sens,reliability(pheno(q1)-predict(resmod,q1),pheno(q2)-predict(resmod,q2),group='glucose_fixed_to_pair_mean'))
 }
}
comp<-do.call(rbind,components);rankval<-comp[comp$metric=='coefficient_times_within_SD',];comp$delta_contribution_rank<-match(comp$component,rankval$component[order(-rankval$value)])
write.csv(comp,file.path(out,'B_components.csv'),row.names=FALSE,na='');write.csv(sens,file.path(out,'B_sensitivity.csv'),row.names=FALSE,na='')
cat('Sensitivity and components computed\n');flush.console()
# Comparison of all second-exam participants with other first-exam adults.
pp<-b$p[b$p$age>=20,];pp$participated<-pp$SEQN%in%a$r$SEQN;pp$accel<-pp$phenoage-predict(resmod,pp);orig<-a$p[match(pp$SEQN,a$p$SEQN),]
represent<-list();k<-1L
for(g in c(FALSE,TRUE)){
 ii<-which(pp$participated==g)
 for(v in c('age','accel')){x<-pp[[v]][ii];represent[[k]]<-data.frame(participated=g,variable=v,category='continuous',n_total=length(ii),n_observed=sum(is.finite(x)),mean=mean(x,na.rm=T),sd=sd(x,na.rm=T),n_category=NA,proportion=NA);k<-k+1L}
 for(v in c('HSSEX','DMARETHN'))for(val in sort(unique(orig[[v]]))){n<-sum(orig[[v]][ii]==val,na.rm=T);represent[[k]]<-data.frame(participated=g,variable=v,category=as.character(val),n_total=length(ii),n_observed=sum(!is.na(orig[[v]][ii])),mean=NA,sd=NA,n_category=n,proportion=n/length(ii));k<-k+1L}
}
write.csv(do.call(rbind,represent),file.path(out,'B_representativeness.csv'),row.names=F,na='')
lod<-data.frame(cohort='NHANES_III',visit=c(1,2),n=length(y1),below_lod=c(sum(a$main$CRP==.21),sum(a$main$CRR==.21)),lod_mgdL=.3)
lod$proportion<-lod$below_lod/lod$n;write.csv(lod,file.path(out,'B_CRP_LOD.csv'),row.names=F)
# HD uses the package's frozen sex-specific reference population.
load('work/source/NHANES3_HDTrain.rda');source('work/source/hd_calc.R')
long<-rbind(transform(z1,sampleID=paste0(SEQN,'_1')),transform(z2,sampleID=paste0(SEQN,'_2')))
long$hd<-NA_real_;long$hd_log<-NA_real_
for(sx in 1:2){ii<-which(long$gender==sx);rr<-NHANES3_HDTrain[NHANES3_HDTrain$gender==sx,];h<-hd_calc(long[ii,setdiff(names(long),c('hd','hd_log'))],rr,bioms)$data;long$hd[ii]<-h$hd;long$hd_log[ii]<-h$hd_log}
h1<-long$hd[seq_along(y1)];h2<-long$hd[length(y1)+seq_along(y1)];secondary<-reliability(h1,h2,label='HD',unit='pooled sex-specific SD of distance')
# KDM coefficients are the default first-exam NHANES III fits, fixed before projection.
load('work/source/NHANES3.rda');source('work/source/kdm_calc.R')
kbm<-c('fev','sbp','totchol','hba1c','albumin','creat','lncrp','alp','bun')
k1<-z1;k2<-z2;firstref<-NHANES3[match(paste0('1991_',a$main$SEQN),NHANES3$sampleID),]
for(v in c('fev','sbp','totchol','hba1c','bun')) k1[[v]]<-firstref[[v]]
k2$fev<-a$main$SPRFEV1;k2$sbp<-a$main$PERMNK1R;k2$totchol<-a$main$TCR;k2$hba1c<-a$main$GHR;k2$bun<-a$main$BUR
k1$kdm<-k2$kdm<-NA_real_;kfit<-list()
for(sx in 1:2){tr<-NHANES3[NHANES3$age>=30&NHANES3$age<=75&NHANES3$pregnant==0&NHANES3$gender==sx,];fit<-kdm_calc(tr,kbm)$fit;kfit[[as.character(sx)]]<-fit
 ii<-which(k1$gender==sx&complete.cases(k1[,kbm])&complete.cases(k2[,kbm]));if(length(ii)>2){k1$kdm[ii]<-kdm_calc(k1[ii,],kbm,fit,fit$s_ba2)$data$kdm;k2$kdm[ii]<-kdm_calc(k2[ii,],kbm,fit,fit$s_ba2)$data$kdm}}
ko<-is.finite(k1$kdm)&is.finite(k2$kdm);km<-lm(kdm0~age,data=NHANES3);kacc1<-k1$kdm-predict(km,k1);kacc2<-k2$kdm-predict(km,k2)
secondary<-rbind(secondary,reliability(kacc1,kacc2,label='KDM_BA_Accel'));write.csv(secondary,file.path(out,'B_secondary.csv'),row.names=F,na='')
saveRDS(list(resmod=resmod,base_n=nrow(base),z1=z1,z2=z2,y1=y1,y2=y2,icc21=icc21,metrics=metrics,reliability=reliability,rel=rel,secondary=secondary,kdm_fit=kfit), 'work/b_results.rds')
writeLines(c(capture.output(summary(resmod)),paste('first exam complete adult N=',nrow(base))),file.path(out,'B_residual_model.txt'))
ci<-rel[rel$group=='all'&rel$metric=='ICC',];half<-(ci$ci_upper-ci$ci_lower)/2
status<-fromJSON(file.path(out,'status.json'));status$B<-list(verdict=ifelse(half<=.03,'PASS','STOP'),ICC=ci$value,ci_lower=ci$ci_lower,ci_upper=ci$ci_upper,half_width=half,bootstrap_replicates=Rboot);status$layer<-'B';status$verdict<-ifelse(half<=.03,'PASS','STOP');status$reason<-ifelse(half<=.03,'B precision gate passed; C pending.','STOP_PRECISION: ICC 95% CI half-width exceeds 0.03.');status$effects_computed<-TRUE;write_json(status,file.path(out,'status.json'),pretty=T,auto_unbox=T)
print(rel[rel$group=='all',],row.names=FALSE);cat('precision_half_width',half,'\n');print(secondary[secondary$metric=='ICC',],row.names=F)
