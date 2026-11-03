O<-'outputs/BIOAGE-REL-01/S';r<-read.csv('outputs/BIOAGE-REL-01/B_reliability.csv');sen<-read.csv('outputs/BIOAGE-REL-01/B_sensitivity.csv');rep<-read.csv('outputs/BIOAGE-REL-01/D_replication.csv');s1<-read.csv(file.path(O,'S1_lncrp_reliability.csv'))
scenarios<-list(
 list(cohort='NHANES_III',scope='main_full',score='log1p',row=r[r$group=='all'&r$metric=='MDC95',],source='v1.0 B_reliability.csv'),
 list(cohort='NHANES_III',scope='both_fasting_ge8h',score='log1p',row=sen[sen$group=='fast_both_ge8h'&sen$metric=='MDC95',],source='v1.0 B_sensitivity.csv'),
 list(cohort='NHANES_2001_2002',scope='replication_fixed407',score='log1p',row=rep[rep$metric=='MDC95',],source='v1.0 D_replication.csv'),
 list(cohort='NHANES_III',scope='main_full',score='ln',row=s1[s1$cohort=='NHANES_III'&s1$score=='ln'&s1$metric=='MDC95',],source='S1_lncrp_reliability.csv'),
 list(cohort='NHANES_2001_2002',scope='replication_fixed407',score='ln',row=s1[s1$cohort=='NHANES_2001_2002'&s1$score=='ln'&s1$metric=='MDC95',],source='S1_lncrp_reliability.csv'))
rows<-list();k<-1L
for(sc in scenarios){z<-sc$row;stopifnot(nrow(z)==1);sdiff<-z$value/1.96;sdlo<-z$ci_lower/1.96;sdhi<-z$ci_upper/1.96
 for(delta in c(1,2,3)){
  raw<-2*(1.96+.84)^2*sdiff^2/delta^2;rawlo<-2*(1.96+.84)^2*sdlo^2/delta^2;rawhi<-2*(1.96+.84)^2*sdhi^2/delta^2
  rows[[k]]<-data.frame(label='写作期预定补充分析',supplement='S4',cohort=sc$cohort,scope=sc$scope,score=sc$score,N_source_people=z$n,delta_years=delta,alpha_two_sided=.05,power=.80,z_alpha=1.96,z_power=.84,SD_diff_years=sdiff,SD_diff_ci_lower=sdlo,SD_diff_ci_upper=sdhi,n_per_group_unrounded=raw,n_per_group=ceiling(raw),n_total=2*ceiling(raw),n_per_group_lower_from_SD_CI=ceiling(rawlo),n_per_group_upper_from_SD_CI=ceiling(rawhi),individual_MDC95_years=z$value,individual_MDC95_ci_lower=z$ci_lower,individual_MDC95_ci_upper=z$ci_upper,attrition_inflation=0,source=sc$source,status='DERIVED_FROM_2000_PERSON_BOOTSTRAP_RESULTS',note='Two independent parallel groups; change-score endpoint; published user-specified normal approximation. Variance scenario based on observational short-interval paired differences; not a fully specified intervention trial.')
  k<-k+1L
 }
}
res<-do.call(rbind,rows);write.csv(res,file.path(O,'S4_sample_size.csv'),row.names=FALSE,na='');print(res[,c('cohort','scope','score','delta_years','SD_diff_years','n_per_group','individual_MDC95_years')],row.names=FALSE)
