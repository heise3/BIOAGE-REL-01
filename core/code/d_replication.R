suppressPackageStartupMessages({library(haven);library(jsonlite)})
a<-readRDS('work/a1.rds');b<-readRDS('work/b_results.rds');cal<-readRDS('work/a2.rds');out<-'outputs/BIOAGE-REL-01';st<-fromJSON(file.path(out,'status.json'));stopifnot(st$B$verdict=='PASS',nrow(a$repd)>=200)
Rboot<-2000L;set.seed(20261104);icc21<-b$icc21;metrics<-b$metrics;reliability<-b$reliability;pheno<-cal$pheno
p1<-as.data.frame(read_xpt('work/source/L11_B.xpt'))[,c('SEQN','LBXCRP')];p2<-as.data.frame(read_xpt('work/source/L25_B.xpt'))[,c('SEQN','LBXWBCSI','LBXLYPCT','LBXMCVSI','LBXRDW')];p3<-as.data.frame(read_xpt('work/source/L40_B.xpt'))[,c('SEQN','LBXSAL','LBDSAPSI','LBDSCR','LBXSGL')]
d<-Reduce(function(x,y)merge(x,y,by='SEQN'),list(a$repd,p1,p2,p3));first<-c('LBXSAL','LBDSAPSI','LBDSCR','LBXSGL','LBXCRP','LBXWBCSI','LBXLYPCT','LBXMCVSI','LBXRDW');second<-c('LB2SAL','LB2SAPSI','LB2SCR','LB2SGL','LB2CRP','LB2WBCSI','LB2LYPCT','LB2MCVSI','LB2RDW');d<-d[complete.cases(d[,c(first,second)]),]
project<-function(bm){z<-d[,bm];names(z)<-c('albumin','alp','creat','glucose','crp','wbc','lymph','mcv','rdw');z$albumin_gL<-z$albumin*10;z$glucose_mmol<-z$glucose*.0555;z$creat_umol<-z$creat*88.4017;z$lncrp<-log1p(z$crp);z$age<-d$RIDAGEYR;z}
z1<-project(first);z2<-project(second);x1<-pheno(z1)-predict(b$resmod,z1);x2<-pheno(z2)-predict(b$resmod,z2)
res<-reliability(x1,x2,label='PhenoAgeAccel',group='NHANES_2001_2002',unit='years');res$note<-'Fixed NHANES III first-exam residual model projected to both replication visits; second-exam complete adult N=415; paired N reported per row.'
write.csv(res,file.path(out,'D_replication.csv'),row.names=F,na='')
lod<-data.frame(cohort='NHANES_2001_2002',visit=c(1,2),n=nrow(d),below_lod=c(sum(d$LBXCRP==.01),sum(d$LB2CRP==.01)),lod_mgdL=.02);lod$proportion<-lod$below_lod/lod$n
old<-read.csv(file.path(out,'B_CRP_LOD.csv'));write.csv(rbind(old,lod),file.path(out,'B_CRP_LOD.csv'),row.names=F,na='')
saveRDS(list(d=d,x1=x1,x2=x2,results=res),'work/d_results.rds');print(res,row.names=F)
