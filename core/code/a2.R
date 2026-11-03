suppressPackageStartupMessages({library(haven);library(dplyr);library(jsonlite)})
out <- 'outputs/BIOAGE-REL-01';a<-readRDS('work/a1.rds')
pheno <- function(x) {
 xb <- -19.90667-.03359355*x$albumin_gL+.009506491*x$creat_umol+.1953192*x$glucose_mmol+.09536762*x$lncrp-.01199984*x$lymph+.02676401*x$mcv+.3306156*x$rdw+.001868778*x$alp+.05542406*x$wbc+.08035356*x$age
 141.50225+(log(.0055305*1.51714/.007692696)+xb)/.090165
}
project_input <- function(d,second=FALSE) {
 bm<-if(second)a$bm2 else a$bm1
 z<-as.data.frame(d[,bm]);names(z)<-c('albumin','alp','creat','glucose','crp','wbc','lymph','mcv','rdw')
 z$age<-d$HSAGEIR;z$gender<-d$HSSEX;z$SEQN<-d$SEQN;z$sampleID<-as.character(d$SEQN)
 z$creat<-.960*z$creat-.184;z$creat_umol<-z$creat*88.4017;z$glucose_mmol<-z$glucose*.0555;z$albumin_gL<-z$albumin*10;z$lncrp<-log1p(z$crp)
 z
}
p<-project_input(a$p);p$phenoage<-pheno(p)
load('work/source/NHANES3.rda');ref<-NHANES3;ref$SEQN<-as.numeric(sub('.*_','',ref$sampleID))
z<-merge(p[,c('SEQN','phenoage')],ref[,c('SEQN','phenoage0')],by='SEQN')
ok<-is.finite(z$phenoage)&is.finite(z$phenoage0)
cal<-data.frame(check=c('PhenoAge_r','PhenoAge_MAD_years'),n=sum(ok),value=c(cor(z$phenoage[ok],z$phenoage0[ok]),mean(abs(z$phenoage[ok]-z$phenoage0[ok]))),target=c('>=0.99','<=0.5'),status='PASS',note='Raw first-exam reconstruction; BioAge published implementation and preprocessing; comparison on available finite phenoage0 overlap')
if(cal$value[1]<.99||cal$value[2]>.5)cal$status<-'CALIBRATION_FAIL'
# Rebuild the kidney paper calibration population independently of the main study.
cys1<-as.data.frame(read_xpt('work/source/sscystat.xpt'))[,c('SEQN','SSCYSTAT')]
k1<-as.data.frame(read_xpt('work/source/SSNH3BTP.xpt'))[,c('SEQN','SSB2M','SSBTP')]
k2<-as.data.frame(read_xpt('work/source/SSNH3CYS.xpt'))
kd<-Reduce(function(x,y)merge(x,y,by='SEQN'),list(a$d,cys1,k1,k2))
fields<-c('CEP','CER','SSCYSTAT','SSCYSCR','SSB2M','SSB2MR','SSBTP','SSBTPR')
kd<-kd[kd$HSAGEIR>=18&complete.cases(kd[,fields]),]
kd$CEP<-.960*kd$CEP-.184;kd$CER<-.960*kd$CER-.184
kd$SSCYSTAT<-1.12*(.022+.800*kd$SSCYSTAT);kd$SSCYSCR<-1.12*kd$SSCYSCR
keep<-rep(TRUE,nrow(kd));bad_counts<-c()
for(pair in list(c('CEP','CER'),c('SSCYSTAT','SSCYSCR'),c('SSBTP','SSBTPR'),c('SSB2M','SSB2MR'))){
 delta<-kd[[pair[2]]]-kd[[pair[1]]];bad<-abs(delta-mean(delta))>3*sd(delta);keep<-keep&!bad;bad_counts<-c(bad_counts,sum(bad))
}
cat('Kidney initial complete N',nrow(kd),'after calibrated per-marker difference >3SD exclusion',sum(keep),'counts',bad_counts,'\n')
stopifnot(nrow(kd)==821,sum(keep)==795,all(bad_counts==c(8,14,10,8)))
kk<-kd[keep,];cv<-100*sd(kk$CER-kk$CEP)/sqrt(2)/mean(c(kk$CEP,kk$CER))
cal<-rbind(cal,data.frame(check='Creatinine_CVw_percent',n=nrow(kk),value=cv,target='[6.1,9.1]',status=ifelse(cv>=6.1&&cv<=9.1,'PASS','CALIBRATION_FAIL'),note=paste('Reconstructed study-specific serum-panel sample; >=18; >3SD calibrated differences by marker; cystatin drift+IFCC standardization per PMC3703771; initial N',nrow(kd),'; excluded counts creat/cys/BTP/B2M=8/14/10/8; unique excluded 26; final n=795; reference CVw=7.6%.')))
write.csv(cal,file.path(out,'A2_calibration.csv'),row.names=FALSE,na='')
saveRDS(list(p=p,pheno=pheno,project_input=project_input,kidney_n_initial=nrow(kd),kidney_n=nrow(kk),cal=cal), 'work/a2.rds')
print(cal,row.names=FALSE)
st<-fromJSON(file.path(out,'status.json'));st$A2<-list(verdict=ifelse(all(cal$status=='PASS'),'PASS','STOP'),PhenoAge_r=cal$value[1],PhenoAge_MAD_years=cal$value[2],creatinine_CVw_percent=cv,kidney_N=nrow(kk));st$layer<-'A2';st$verdict<-ifelse(all(cal$status=='PASS'),'PASS','STOP');st$reason<-ifelse(all(cal$status=='PASS'),'A2 passed; B pending.','CALIBRATION_FAIL');write_json(st,file.path(out,'status.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
stopifnot(all(cal$status=='PASS'))
