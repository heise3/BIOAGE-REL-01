# From S7/: Rscript code/run.R /path/to/read-only/REL-workspace
suppressPackageStartupMessages({library(survival);library(haven);library(jsonlite)})
root<-normalizePath(commandArgs(TRUE)[1]); work<-file.path(root,'work'); src<-file.path(work,'source'); rel<-file.path(root,'outputs/BIOAGE-REL-01')
stopifnot(!file.exists('A3_calibration.tsv'),!file.exists('tableS7_reliability.tsv'))
out<-function(d,n) write.table(d,n,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
a<-readRDS(file.path(work,'a1.rds')); b<-readRDS(file.path(work,'b_results.rds')); r<-readRDS(file.path(work,'d_results.rds'))
cali<-readRDS(file.path(work,'a2.rds')); sup<-readRDS(file.path(work,'S/s1_results.rds'))
icc<-b$icc21; m<-a$main
stopifnot(!anyDuplicated(m$SEQN),!anyDuplicated(r$d$SEQN),identical(as.numeric(m$SEQN),as.numeric(b$z1$SEQN)))
fields3<-list(c('GRP','LMP','MOP','PLPSI','TGP','G1P','PHPFAST'),c('GRR','LMR','MOR','PLRSI','TGR','G1R','PHRFAST'))
fields2<-list(c('LBDNENO','LBDLYMNO','LBDMONO','LBXPLTSI','LBXSTR','LBXSGL'),c('LB2NENO','LB2LYMNO','LB2MONO','LB2PLTSI','LB2STR','LB2SGL'))
read3<-function(v){
 nm<-c('lab','labse')[v];schema<-fromJSON(file.path(src,paste0(nm,'_schema.json')))
 lines<-readLines(file.path(src,paste0(nm,'.dat'))); pos<-schema$SEQN; ids<-as.integer(substr(lines,pos[1],pos[2]));stopifnot(!anyDuplicated(ids),all(m$SEQN%in%ids))
 lines<-lines[match(m$SEQN,ids)]; lapply(fields3[[v]],function(k){
  p<-schema[[k]]; stopifnot(length(p)==2);raw<-trimws(substr(lines,p[1],p[2])); raw[raw==strrep('8',p[2]-p[1]+1)]<-NA_character_;as.numeric(raw)
 })
}
read2<-function(v){
 nm<-c('L25_B','L25_2_B','L40_B','L40_2_B')
 z<-lapply(nm[c(v,v+2)],function(n){d<-as.data.frame(read_xpt(file.path(src,paste0(n,'.xpt'))));stopifnot(!anyDuplicated(d$SEQN),all(r$d$SEQN%in%d$SEQN));d[match(r$d$SEQN,d$SEQN),]})
 h<-lapply(seq_len(6),function(i)as.numeric(z[[if(i<=4)1 else 2]][[fields2[[v]][i]]]))
 ph<-read_xpt(file.path(work,'S6_sources/PH_B.xpt'))
 h[[7]]<-if(v==1)(ph$PHAFSTHR+ph$PHAFSTMN/60)[match(r$d$SEQN,ph$SEQN)] else rep(NA_real_,nrow(r$d));h
}
derive<-function(z,ref,nameset){
 names(z)<-c('N','L','M','P','T','G','F')
 fast<-is.finite(z$F)&z$F>=8
 g<-ifelse(fast,z$G,NA_real_);ty<-log(z$T*g/2)
 vals<-list(ref,z$N/z$L,z$P/z$L,z$P*z$N/z$L,z$N*z$M/z$L,ty,z$N,z$L,z$M,z$P,z$T,g)
 setNames(vals,nameset)
}
co<-c('NHANES_III','NHANES_2001_2002')
labs<-list(c('PhenoAgeAccel','GLR','PLR','SII_g','SIRI_g','TyG','Granulocytes','Lymphocytes','Mononuclear','Platelets','TG','Fasting_glucose'),c('PhenoAgeAccel','NLR','PLR','SII','SIRI','TyG','Neutrophils','Lymphocytes','Monocytes','Platelets','TG','Fasting_glucose'))
# Reproject frozen score coefficients and age residual model onto the same people.
scores<-lapply(co,function(id){ss<-sup$cohorts[[id]];lapply(ss[c('z1','z2')],function(z)cali$pheno(z)-predict(b$resmod,z))})
stopifnot(identical(as.numeric(r$d$SEQN),as.numeric(sup$cohorts[[co[2]]]$z1$SEQN)),identical(as.numeric(m$SEQN),as.numeric(sup$cohorts[[co[1]]]$z1$SEQN)))
pairs<-lapply(1:2,function(j)lapply(1:2,function(v)derive(if(j==1)read3(v)else read2(v),scores[[j]][[v]],labs[[j]])))
target<-function(n,j,k){if(j==2 && k%in%c(6,12))'NOT_RUN_A1' else if(n<c(900,270)[j])'NOT_RUN_A2' else 'PASS'}
counts<-list(); jobs<-list()
for(j in 1:2)for(k in seq_along(labs[[j]])){
 x<-pairs[[j]][[1]][[k]];y<-pairs[[j]][[2]][[k]]
 # Log of each positive conventional component or index, including the pre-defined TyG index.
 ok1<-is.finite(x)&(k==1|x>0);ok2<-is.finite(y)&(k==1|y>0);ii<-which(ok1&ok2);gate<-target(length(ii),j,k)
 counts[[length(counts)+1]]<-data.frame(dataset=co[j],analyte=labs[[j]][k],N_base=length(x),visit1_valid=sum(ok1),visit2_valid=if(gate=='NOT_RUN_A1')NA_integer_ else sum(ok2),both_valid=if(gate=='NOT_RUN_A1')NA_integer_ else length(ii),verified_pairs=length(ii),threshold=c(900,270)[j],status=gate)
 jobs[[length(jobs)+1]]<-list(j=j,k=k,i=ii,x=x[ii],y=y[ii],gate=gate)
}
out(do.call(rbind,counts),'A2_counts.tsv')
d<-data.frame(age=m$HSAGEIR,sex=factor(m$HSSEX),race=factor(m$DMARETHN),time=m$PERMTH_EXM,event=m$MORTSTAT)
mortok<-which(m$ELIGSTAT==1&d$time>0&complete.cases(d));stopifnot(length(mortok)==nrow(m))
design<-model.matrix(~age+sex+race,d)[,-1,drop=FALSE];Y<-Surv(d$time,d$event)
model<-function(x,y,idx){
 sd0<-sd(x);list(XX=lapply(list(x,y,(x+y)/2),function(z)cbind(x=z/sd0,design[idx,,drop=FALSE])),Y=Y[idx,,drop=FALSE],sd0=sd0)
}
betas<-function(M,ii,with_p=FALSE){
 fit<-function(X){
 f<-coxph.fit(X[ii,,drop=FALSE],M$Y[ii,,drop=FALSE],strata=NULL,offset=NULL,init=NULL,control=coxph.control(),weights=NULL,method='efron',rownames=NULL,resid=FALSE,nocenter=c(-1,0,1))
 beta<-unname(f$coefficients[1]);if(with_p)c(beta,2*pnorm(-abs(beta/sqrt(f$var[1,1]))))else beta
 }
 sapply(M$XX,fit)
}
frozenICC<-sapply(c('B_reliability.csv','D_replication.csv'),function(fn){z<-read.csv(file.path(rel,fn));z$value[which(z$metric=='ICC')[1]]})
frozenHR<-read.csv(file.path(rel,'C_mortality.csv'));frozenHR<-frozenHR[match(c('HR_exam1','HR_exam2','HR_mean'),frozenHR$metric),]
calrows<-list()
for(j in 1:2){
 e<-scores[[j]];obs<-icc(e[[1]],e[[2]]);n<-length(e[[1]])
 z<-read.csv(file.path(rel,c('B_reliability.csv','D_replication.csv')[j]));nn<-z$n[which(z$metric=='ICC')[1]]
 calrows[[j]]<-data.frame(dataset=co[j],metric='ICC',n=n,n_frozen=nn,deaths=NA,deaths_frozen=NA,value=obs,frozen=frozenICC[j],abs_difference=abs(obs-frozenICC[j]),tolerance=.005,pass=n==nn&&abs(obs-frozenICC[j])<=.005)
}
i<-mortok; M<-model(scores[[1]][[1]][i],scores[[1]][[2]][i],i);v<-betas(M,seq_along(i))
for(h in 1:3){z<-frozenHR[h,];dif<-abs(v[h]-log(z$value));calrows[[2+h]]<-data.frame(dataset=co[1],metric=paste0('log_',z$metric),n=length(i),n_frozen=z$n,deaths=sum(d$event[i]),deaths_frozen=z$deaths,value=v[h],frozen=log(z$value),abs_difference=dif,tolerance=.01,pass=length(i)==z$n&&sum(d$event[i])==z$deaths&&dif<=.01)}
cal<-do.call(rbind,calrows);out(cal,'A3_calibration.tsv');stopifnot(all(cal$pass))
q4<-function(x,y){t1<-quantile(x,.75,names=FALSE);t2<-quantile(y,.75,names=FALSE);mean(y[x>=t1]<t2)}
ci<-function(x)quantile(x,c(.025,.975),names=FALSE,na.rm=TRUE)
analyze<-function(job){
 j<-job$j;k<-job$k;idx<-job$i;name<-labs[[j]][k];reason<-switch(job$gate,NOT_RUN_A1='Second-visit fasting duration unavailable',NOT_RUN_A2='Below frozen paired-N threshold',PASS='')
 row<-data.frame(dataset=co[j],analyte=name,N=if(job$gate=='NOT_RUN_A1')NA_integer_ else length(idx),status=if(reason=='')'COMPUTED' else 'NOT_RUN',reason,scale=if(k==1)'years' else 'natural log',ICC=NA_real_,ICC_lo=NA_real_,ICC_hi=NA_real_,CVw_percent=NA_real_,MDC95=NA_real_,Q4_first_n=NA_integer_,Q4_exit_n=NA_integer_,Q4_exit=NA_real_,Q4_lo=NA_real_,Q4_hi=NA_real_,reference_ICC_same_people=NA_real_,delta_ICC=NA_real_,delta_lo=NA_real_,delta_hi=NA_real_,B_success=0L)
 mr<-data.frame(dataset=co[j],analyte=name,N=row$N,deaths=NA_integer_,status=row$status,reason,common_SD=NA_real_,HR1=NA_real_,lo1=NA_real_,hi1=NA_real_,HR2=NA_real_,lo2=NA_real_,hi2=NA_real_,HRmean=NA_real_,lomean=NA_real_,himean=NA_real_,logHR2_minus_1=NA_real_,diff_lo=NA_real_,diff_hi=NA_real_,p1=NA_real_,p2=NA_real_,p_flip=NA,B_success=0L)
 if(job$gate!='PASS')return(list(row,mr))
 x0<-job$x;y0<-job$y;x<-if(k==1)x0 else log(x0);y<-if(k==1)y0 else log(y0)
 u<-scores[[j]][[1]][idx];v<-scores[[j]][[2]][idx];n<-length(x)
 stats<-function(ii){ic<-icc(x[ii],y[ii]);c(ic,q4(x[ii],y[ii]),ic-icc(u[ii],v[ii]))}
 set.seed(c(20261102,20261104)[j]);bt<-replicate(2000,stats(sample.int(n,n,replace=TRUE)))
 s<-stats(seq_len(n));row$ICC<-s[1];row[1,c('ICC_lo','ICC_hi')]<-ci(bt[1,]);row$Q4_exit<-s[2];row[1,c('Q4_lo','Q4_hi')]<-ci(bt[2,])
 top<-x>=quantile(x,.75);row$Q4_first_n<-sum(top);row$Q4_exit_n<-sum(top&y<quantile(y,.75))
 row$delta_ICC<-s[3];row[1,c('delta_lo','delta_hi')]<-ci(bt[3,]);row$reference_ICC_same_people<-icc(u,v);row$MDC95<-1.96*sd(y-x);row$B_success<-sum(apply(bt,2,function(z)all(is.finite(z))))
 if(k!=1)row$CVw_percent<-100*sd(y0-x0)/sqrt(2)/mean(c(x0,y0))
 if(j==1){
  M<-model(x,y,idx);v0<-betas(M,seq_len(n),TRUE);mr$deaths<-sum(d$event[idx]);mr$common_SD<-M$sd0
  set.seed(20261103);bb<-replicate(2000,{ii<-sample.int(n,n,TRUE);tryCatch(betas(M,ii),error=function(e)rep(NA_real_,3))})
  for(h in 1:3){suf<-c('1','2','mean')[h];mr[[paste0('HR',suf)]]<-exp(v0[1,h]);mr[1,paste0(c('lo','hi'),suf)]<-exp(ci(bb[h,]))}
  mr$logHR2_minus_1<-v0[1,2]-v0[1,1];mr[1,c('diff_lo','diff_hi')]<-ci(bb[2,]-bb[1,]);mr$p1<-v0[2,1];mr$p2<-v0[2,2];mr$p_flip<-(mr$p1<.05)!=(mr$p2<.05);mr$B_success<-sum(is.finite(colSums(bb)))
 }
 list(row,mr)
}
ans<-parallel::mclapply(jobs,analyze,mc.cores=4,mc.set.seed=FALSE,mc.preschedule=FALSE)
stopifnot(!any(vapply(ans,inherits,logical(1),'try-error')))
out(do.call(rbind,lapply(ans,'[[',1)),'tableS7_reliability.tsv')
out(do.call(rbind,lapply(ans[1:12],'[[',2)),'tableS7_mortality.tsv')
