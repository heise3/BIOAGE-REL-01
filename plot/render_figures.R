# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# v0.3: scatter/BA, histogram, CI, signed stacks, forest and curve panels -> parameter inheritance
# Inspected PCA, BarAblation, LineTrend/plot_sweep.py and BarComposition/plot_rewriting.py + companion PNGs.
# Assets use unrelated PCA/performance/proportion computations; no native asset computations are run.
# PCA/ablation assets perform different computations and cannot be run on frozen score summaries.
# Only the inspected assets' visual principles are inherited; no PCA, tests or score fits are run.

# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
library(ggplot2)
theme_cns <- theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7, color = "#333333"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    strip.text = element_text(size = 8, face = "bold"),
    panel.grid = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank()
  )
# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
categorical <- c("#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666")
categorical_extended <- c(
  "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
  "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999"
)
diverging  <- c("#2166AC", "#F7F7F7", "#B2182B")
sequential <- c("#F7FBFF", "#6BAED6", "#08306B")
accent_red <- "#B2182B"
grey       <- "#999999"
black      <- "#222222"
# Academic Figure Skill Export Baseline — COPY VERBATIM
save_cns_figure <- function(plot, filename, width_mm = 183, height_mm = NULL) {
  ggsave(paste0(filename, ".pdf"), plot, device = cairo_pdf,
         width = width_mm, height = height_mm, units = "mm", dpi = 300)
  png(paste0(filename, ".png"), width = width_mm, height = height_mm,
      units = "mm", res = 300, type = "cairo")
  print(plot)
  dev.off()
}

# Explicit v0.3 user overrides: Okabe–Ito semantics, 8 pt floor, 600 dpi TIFF.
# Official Medical Sciences guide checked: 170 mm maximum double-column width.
suppressPackageStartupMessages({library(patchwork);library(ragg);library(jsonlite)})
bundle_dir<-if(length(commandArgs(trailingOnly=TRUE)))commandArgs(trailingOnly=TRUE)[1]else getwd()
O<-file.path(bundle_dir,'figures')
col_co<-c('NHANES III'='#0072B2','NHANES 2001–2002'='#E69F00')
col_dir<-c('Younger'='#009E73','Middle'='#BDBDBD','Older'='#D55E00')
col_comp<-c('Glucose [G]'='#CC79A7','Creatinine [C]'='#56B4E9','RDW [R]'='#F0E442','Other six'='#999999')
pt_size<-8/.pt
th<-theme_cns+theme(axis.text=element_text(size=8),legend.text=element_text(size=8),panel.border=element_blank(),panel.grid.major.y=element_line(colour='#E6E6E6',linewidth=.2),panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour='#444444',linewidth=.21),axis.ticks=element_line(linewidth=.21),axis.ticks.length=grid::unit(1,'mm'),legend.position='bottom',legend.title=element_blank(),legend.key.size=grid::unit(3,'mm'),plot.tag=element_text(size=9,face='bold'),plot.tag.position=c(0,1),plot.margin=margin(5,5,5,5),plot.caption=element_text(size=8,hjust=0))
tag<-function(p,t)p+labs(tag=t)+th
export_figure<-function(p,stem,h){save_cns_figure(p,file.path(O,stem),170,h);agg_tiff(file.path(O,paste0(stem,'.tiff')),width=170,height=h,units='mm',res=600,compression='lzw');print(p);dev.off()}
coords<-read.csv(file.path(O,'paired_scores.csv'))
rel<-read.csv(file.path(O,'plot_reliability.csv'));directions<-read.csv(file.path(O,'plot_direction.csv'))
hr<-read.csv(file.path(O,'plot_hr.csv'));scales<-read.csv(file.path(O,'plot_scales.csv'))
co_labels<-c('NHANES_III'='NHANES III','NHANES_2001_2002'='NHANES 2001–2002')
getrel<-function(co,m,group='all'){
 z<-rel[rel$cohort==co&rel$score=='log1p'&rel$group==group&rel$metric==m,];stopifnot(nrow(z)==1);z
}
annotations<-list();add_annotation<-function(figure,panel,label,value,source,source_col='',operation='Existing estimate, display formatting'){
 annotations[[length(annotations)+1]]<<-data.frame(figure=figure,panel=panel,label=label,value=as.character(value),source=source,source_column=source_col,operation=operation)
}
hist_qa<-list()
scatter<-function(co){
 z<-coords[coords$cohort==co,];lab<-unname(co_labels[co]);mdc<-getrel(lab,'MDC95')$value;icc<-getrel(lab,'ICC')$value
 lim<-range(c(z$draw1,z$draw2));lim<-lim+c(-1,1)*diff(lim)*.055
 band<-data.frame(x=lim,lower=lim-mdc,upper=lim+mdc)
 add_annotation(if(co=='NHANES_III')'Figure 1'else'Figure S2','A','ICC',sprintf('%.2f',icc),'plot_reliability.csv','value')
 ggplot(z,aes(draw1,draw2))+geom_ribbon(data=band,aes(x=x,ymin=lower,ymax=upper),inherit.aes=FALSE,fill='#EEEEEE')+geom_abline(slope=1,intercept=0,colour='#555555',linewidth=.35)+geom_point(colour=col_co[[lab]],alpha=.3,size=.65,stroke=0)+annotate('text',x=lim[2]-.05*diff(lim),y=lim[1]+.075*diff(lim),hjust=1,vjust=0,label=paste('ICC',sprintf('%.2f',icc)),size=pt_size,family='Arial')+coord_equal(xlim=lim,ylim=lim,expand=FALSE)+labs(x='PhenoAgeAccel, draw 1 (years)',y='PhenoAgeAccel, draw 2 (years)')+th
}
histogram<-function(co){
 z<-coords[coords$cohort==co,];lab<-unname(co_labels[co]);d<-z$difference;mdc<-getrel(lab,'MDC95')$value
 # Bin and overflow counts are plotting operations on the existing vector, not new estimators.
 left<-floor(d/.5)*.5;dr<-ifelse(d<= -2,'Younger',ifelse(d>=2,'Older','Middle'))
 h<-aggregate(rep(1,length(d)),list(left=left,direction=dr),sum);names(h)[3]<-'count';h$right<-h$left+.5
 h$direction<-factor(h$direction,levels=names(col_dir));h$cohort<-co
 write.csv(h,file.path(O,paste0('histogram_bins_',co,'.csv')),row.names=FALSE)
 ymax<-max(tapply(h$count,h$left,sum)[as.numeric(names(tapply(h$count,h$left,sum)))>= -15&as.numeric(names(tapply(h$count,h$left,sum)))<15])
 dn<-directions[directions$cohort==co&directions$score=='log1p'&directions$stratum=='all'&directions$threshold_years==2,]
 young<-dn$proportion[dn$direction=='decrease']*100;old<-dn$proportion[dn$direction=='increase']*100
 nleft<-sum(d< -15);nright<-sum(d>15);nshown<-sum(d>= -15&d<=15)
 stopifnot(sum(h$count)==nrow(z),nleft+nright+nshown==nrow(z))
 hist_qa[[co]]<<-data.frame(cohort=co,n=nrow(z),displayed_window_n=nshown,left_beyond_n=nleft,right_beyond_n=nright,bin_width_years=.5,all_points_accounted_for=TRUE)
 for(a in c('young','old'))add_annotation(if(co=='NHANES_III')'Figure 1'else'Figure S2','B',a,sprintf('%.1f%%',if(a=='young')young else old),'plot_direction.csv','proportion')
 for(a in c('left','right'))add_annotation(if(co=='NHANES_III')'Figure 1'else'Figure S2','B',paste(a,'beyond ±15'),if(a=='left')nleft else nright,'paired_scores.csv','difference','Explicit plot-window overflow count')
 ggplot(h,aes((left+right)/2,count,fill=direction))+geom_col(width=.5,colour='#555555',linewidth=.12)+scale_fill_manual(values=col_dir,guide='none')+annotate('segment',x=c(-mdc,mdc),xend=c(-mdc,mdc),y=0,yend=ymax*1.13,linetype='dashed',colour='#555555',linewidth=.35)+annotate('text',x=-14.8,y=ymax*1.35,hjust=0,vjust=1,label=sprintf('Appeared ≥2 y younger:\n%.1f%%',young),colour=col_dir[['Younger']],size=pt_size,family='Arial')+annotate('text',x=14.8,y=ymax*1.35,hjust=1,vjust=1,label=sprintf('Appeared ≥2 y older:\n%.1f%%',old),colour=col_dir[['Older']],size=pt_size,family='Arial')+annotate('segment',x=-13.4,xend=-14.9,y=ymax*.85,yend=ymax*.85,arrow=grid::arrow(length=grid::unit(1.2,'mm')),linewidth=.3)+annotate('text',x=-14.7,y=ymax*.91,hjust=0,label=paste(nleft,'below −15'),size=pt_size,family='Arial')+annotate('segment',x=13.4,xend=14.9,y=ymax*.85,yend=ymax*.85,arrow=grid::arrow(length=grid::unit(1.2,'mm')),linewidth=.3)+annotate('text',x=14.7,y=ymax*.91,hjust=1,label=paste(nright,'above 15'),size=pt_size,family='Arial')+annotate('text',x=0,y=ymax*1.06,label=sprintf('MDC95 = %.1f y',mdc),size=pt_size,family='Arial')+coord_cartesian(xlim=c(-15,15),ylim=c(0,ymax*1.43),expand=FALSE)+scale_x_continuous(breaks=seq(-15,15,5))+labs(x='Difference, draw 2 − draw 1 (years)',y='Participants (n)')+th
}
subgroups<-data.frame(cohort=c(rep('NHANES III',6),'NHANES 2001–2002'),group=c('all','20-39','40-59','60+','fast_both_ge8h','fast_both_lt8h','all'),label=c('Overall','Age 20–39 y','Age 40–59 y','Age ≥60 y','Both fasting ≥8 h','Both fasting <8 h','Replication'))
for(i in 1:nrow(subgroups)){r<-getrel(subgroups$cohort[i],'MDC95',subgroups$group[i]);subgroups$value[i]<-r$value;subgroups$lower[i]<-r$lower;subgroups$upper[i]<-r$upper}
subgroups$label<-factor(subgroups$label,levels=rev(subgroups$label));subgroups$cohort<-factor(subgroups$cohort,levels=names(col_co))
pc<-ggplot(subgroups,aes(value,label,colour=cohort,shape=cohort))+annotate('segment',x=c(1,2),xend=c(1,2),y=.6,yend=7.35,colour='#888888',linetype='dotted',linewidth=.35)+geom_errorbar(aes(xmin=lower,xmax=upper),orientation='y',width=.15,linewidth=.45)+geom_point(size=1.8)+scale_colour_manual(values=col_co,guide='none')+scale_shape_manual(values=c(16,17),guide='none')+annotate('text',x=c(1,2),y=7.75,label=c('1 y','2 y'),size=pt_size,family='Arial')+scale_x_continuous(limits=c(0,12),breaks=seq(0,12,2))+scale_y_discrete(expand=expansion(add=c(.4,1.2)))+labs(x='Minimal detectable change, MDC95 (years; 95% CI)',y=NULL)+th
fig1<-wrap_plots(tag(scatter('NHANES_III'),'A'),tag(histogram('NHANES_III'),'B'),tag(pc,'C'),design='AB\nCC',heights=c(1.15,1))
export_figure(fig1,'figure1_ranking_change',166)

comp<-read.csv(file.path(O,'plot_components.csv'));stopifnot(nrow(comp)==18)
ord<-comp$label[comp$cohort=='NHANES III'][order(comp$value[comp$cohort=='NHANES III'],decreasing=TRUE)]
comp$label<-factor(comp$label,levels=rev(ord));comp$cohort<-factor(comp$cohort,levels=names(col_co));dodge<-position_dodge(.5)
pa<-ggplot(comp,aes(value,label,colour=cohort,shape=cohort))+geom_vline(xintercept=0,colour='#777777',linetype='dotted',linewidth=.35)+geom_errorbar(aes(xmin=lower,xmax=upper),orientation='y',position=dodge,width=.15,linewidth=.45)+geom_point(position=dodge,size=1.8)+scale_colour_manual(values=col_co)+scale_shape_manual(values=c(16,17))+labs(x='Reduction in MDC95 after fixing\ncomponent to paired mean (%)',y=NULL,caption='95% CIs; reductions are not additive')+th+guides(colour=guide_legend(nrow=2),shape=guide_legend(nrow=2))
ext<-read.csv(file.path(O,'plot_extremes.csv'));stopifnot(nrow(ext)==216)
people<-unique(ext[,c('cohort','person','rank_abs_d')]);people<-people[order(match(people$cohort,c('NHANES_III','NHANES_2001_2002')),people$rank_abs_d),];people$y<-24:1
ext$y<-people$y[match(ext$person,people$person)];ext$fill<-ifelse(ext$analyte=='glucose_mmol','Glucose [G]',ifelse(ext$analyte=='creat_umol','Creatinine [C]',ifelse(ext$analyte=='rdw','RDW [R]','Other six')))
ext$fill<-factor(ext$fill,levels=names(col_comp));ext$xmin<-ext$xmax<-0
for(person in people$person){ii<-which(ext$person==person);pos<-neg<-0
 for(j in ii){v<-ext$contribution_years[j];if(v>=0){ext$xmin[j]<-pos;pos<-pos+v;ext$xmax[j]<-pos}else{ext$xmax[j]<-neg;neg<-neg+v;ext$xmin[j]<-neg}}
}
ext$code<-ifelse(ext$fill=='Glucose [G]','G',ifelse(ext$fill=='Creatinine [C]','C',ifelse(ext$fill=='RDW [R]','R','')))
pb<-ggplot(ext)+geom_rect(aes(xmin=xmin,xmax=xmax,ymin=y-.34,ymax=y+.34,fill=fill),colour='#333333',linewidth=.2)+geom_vline(xintercept=0,colour='#444444',linewidth=.35)+geom_hline(yintercept=5.5,colour='#555555',linewidth=.4)+geom_text(data=ext[ext$code!=''&abs(ext$contribution_years)>=6,],aes(x=(xmin+xmax)/2,y=y,label=code),size=pt_size,family='Arial')+scale_fill_manual(values=col_comp,drop=FALSE)+scale_y_continuous(breaks=people$y,labels=people$person,expand=expansion(add=.7))+labs(x='Signed component contribution to\ndraw 2 − draw 1 (years)',y='Participant (P: NHANES III; R: replication)')+th+guides(fill=guide_legend(nrow=2,byrow=TRUE))+theme(panel.grid.major.y=element_blank())
scales$cohort<-factor(scales$cohort,levels=names(col_co));scales$row<-ifelse(scales$cohort=='NHANES III',2,1);scales$method<-factor(scales$method,levels=c('Conventional','MAD-based'))
pcc<-ggplot(scales,aes(value,row,colour=cohort))+geom_line(aes(group=cohort),colour='#999999',linewidth=.5)+geom_point(aes(shape=method),size=2.2,stroke=.7)+geom_text(aes(label=sprintf('%.1f',value)),nudge_y=.19,size=pt_size,family='Arial',show.legend=FALSE)+scale_colour_manual(values=col_co,guide='none')+scale_shape_manual(values=c(16,1))+scale_y_continuous(breaks=c(1,2),labels=c('2001–2002','NHANES III'),limits=c(.55,2.55))+scale_x_continuous(limits=c(3.8,8),breaks=c(4,6,8))+labs(x='Normal-calibrated MDC scale (years)',y=NULL)+th+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())+guides(shape=guide_legend(nrow=2))
fig2<-wrap_plots(tag(pa,'A'),tag(pb,'B'),tag(pcc,'C'),design='AB\nCB',widths=c(1.12,1),heights=c(1.6,1))
export_figure(fig2,'figure2_noise_budget',183)
write.csv(ext,file.path(O,'component_stack_geometry.csv'),row.names=FALSE)

forest<-function(co){
 z<-hr[hr$cohort==co,];z$y<-nrow(z):1
 ggplot(z,aes(value,y))+geom_vline(xintercept=1,linetype='dotted',colour='#777777',linewidth=.35)+geom_errorbar(aes(xmin=lower,xmax=upper),orientation='y',width=.15,colour=col_co[[co]],linewidth=.5)+geom_point(colour=col_co[[co]],shape=ifelse(z$label=='Corrected draw 2',18,16),size=2)+geom_text(aes(x=1.03,y=y+.3,label=sprintf('%.2f (%.2f–%.2f)',value,lower,upper)),hjust=0,size=pt_size,family='Arial')+scale_y_continuous(breaks=z$y,labels=z$label,limits=c(.5,nrow(z)+.6))+scale_x_continuous(limits=c(.95,max(1.6,z$upper+.04)),breaks=c(1,1.25,1.5,2))+labs(x='Mortality HR per 1 SD (6.8 years; 95% CI)',y=NULL)+th
}
side_text<-function(text){ggplot()+annotate('text',x=0,y=1,label=text,hjust=0,vjust=1,size=pt_size,family='Arial',lineheight=1.2)+xlim(0,1)+ylim(0,1)+theme_void()+theme(plot.margin=margin(8,3,3,1))}
fa<-wrap_plots(tag(forest('NHANES III'),'A'),side_text('λ = 0.83\n\nΔC-index\n(mean vs draw 1):\n0.0009\n(−0.0005 to 0.0021)\n\nApparent discrimination'),widths=c(2.25,1))
curves<-read.csv(file.path(O,'sample_size_curve.csv'));curves$scenario<-factor(curves$scenario,levels=c('All participants','Both fasting ≥8 h'))
pp<-curves[curves$delta_years==1,]
fc<-ggplot(curves,aes(delta_years,n_unrounded,linetype=scenario))+geom_line(colour=col_co[[1]],linewidth=.6)+geom_point(data=pp,aes(y=n_ceiling,shape=scenario),colour=col_co[[1]],size=2.2)+geom_text(data=pp,aes(y=n_ceiling,label=n_ceiling),nudge_x=c(.11,-.13)[match(pp$scenario,levels(curves$scenario))],hjust=ifelse(pp$scenario=='All participants',0,1),nudge_y=ifelse(pp$scenario=='All participants',.1,-.10),size=pt_size,family='Arial',show.legend=FALSE)+geom_vline(xintercept=1,colour='#BBBBBB',linetype='dotted',linewidth=.3)+annotate('text',x=.5,y=1600,label='A single person: change must exceed 7.1 y (MDC95)\nto exceed the normal-calibrated retest scale.',hjust=0,vjust=1,size=pt_size,family='Arial')+scale_y_log10(limits=c(10,1800),breaks=c(10,20,50,100,200,500,1000))+scale_x_continuous(breaks=seq(.5,3,.5),limits=c(.5,3))+scale_linetype_manual(values=c('solid','longdash'))+scale_shape_manual(values=c(16,17))+labs(x='True between-group difference, Δ (years)',y='Participants required per group (log scale)',caption='Two-sided α = 0.05; 80% power. Difference SD: 3.63 y overall; 3.07 y with paired fasting.')+th+guides(linetype=guide_legend(nrow=1),shape=guide_legend(nrow=1))
fig3<-wrap_plots(fa,tag(fc,'B'),ncol=1,heights=c(1,1.28))
export_figure(fig3,'figure3_implications',171)

# Supplementary Figure S1: full-range Bland–Altman, all participants retained.
ylim<-range(coords$difference);ylim<-ylim+c(-1,1)*diff(ylim)*.06
ba<-lapply(names(co_labels),function(co){z<-coords[coords$cohort==co,];lab<-unname(co_labels[co]);bias<-getrel(lab,'BA_bias')$value;lim<-c(getrel(lab,'BA_lower')$value,getrel(lab,'BA_upper')$value)
 ggplot(z,aes(paired_mean,difference,colour=lab))+geom_point(alpha=.3,size=.65,stroke=0)+geom_hline(yintercept=0,colour='#BBBBBB',linewidth=.25)+geom_hline(yintercept=bias,colour='#555555',linewidth=.4)+geom_hline(yintercept=lim,colour='#555555',linewidth=.4,linetype='dashed')+scale_colour_manual(values=col_co)+coord_cartesian(ylim=ylim)+labs(x='Mean PhenoAgeAccel across draws (years)',y='Difference, draw 2 − draw 1 (years)')+th
})
export_figure((wrap_plots(tag(ba[[1]],'A'),tag(ba[[2]],'B'),guides='collect') & theme(legend.position='bottom')),'figureS1_bland_altman',91)
fs2c<-wrap_plots(tag(forest('NHANES 2001–2002'),'C'),side_text('76 deaths\nExploratory associations\n\nNo corrected HR or λ\nwas previously estimated;\nnone is added here.'),widths=c(2.25,1))
fs2<-wrap_plots(tag(scatter('NHANES_2001_2002'),'A'),tag(histogram('NHANES_2001_2002'),'B'),fs2c,design='AB\nCC',heights=c(1.25,1))
export_figure(fs2,'figureS2_replication',157)
write.csv(do.call(rbind,hist_qa),file.path(O,'histogram_window_QA.csv'),row.names=FALSE)
write.csv(do.call(rbind,annotations),file.path(O,'annotation_render_log.csv'),row.names=FALSE)
write_json(list(no_score_or_model_refit=TRUE,new_analysis='Only sample-size formula curve; binning, overflow counts and signed-stack geometry are rendering operations',all_scatter_points=as.list(table(coords$cohort)),all_signed_contributions=nrow(ext),all_extreme_people=nrow(people),components=18,width_mm=170,font_floor_pt=8,tag_pt=9,tiff_dpi=600,histogram_window='−15 to +15 years, overflows explicitly labeled',replication_corrected_HR='NOT_RUN: no prior estimate; no new analysis authorized'),file.path(O,'render_QA.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(O,'../code/render_sessionInfo.txt'))

