# Reconstruct participant graphics locally, from the frozen S1 intermediate.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2)
s <- readRDS(args[1])
out <- args[2]
cf <- c(albumin_gL=-.03359355, creat_umol=.009506491, glucose_mmol=.1953192,
        lncrp=.09536762, lymph=-.01199984, mcv=.02676401, rdw=.3306156,
        alp=.001868778, wbc=.05542406)
points <- extremes <- list()
for (co in names(s$scores)) {
  z <- s$scores[[co]]$log1p
  d <- z$x2-z$x1
  points[[co]] <- data.frame(cohort=co, person_row=seq_along(d), draw1=z$x1,
                             draw2=z$x2, paired_mean=(z$x1+z$x2)/2, difference=d)
  a <- s$cohorts[[co]]$z1
  b <- s$cohorts[[co]]$z2
  a$lncrp <- log1p(a$crp)
  b$lncrp <- log1p(b$crp)
  top <- head(order(-abs(d), a$SEQN), ceiling(length(d)*.01))
  extremes[[co]] <- do.call(rbind, lapply(seq_along(top), function(rank) {
    i <- top[rank]
    changes <- as.numeric(b[i,names(cf)]-a[i,names(cf)])*cf/.090165
    stopifnot(abs(sum(changes)-d[i]) < 1e-9)
    data.frame(cohort=co, person=sprintf('%s%02d',if(co=='NHANES_III')'P' else 'R',rank),
               rank_abs_d=rank, analyte=names(cf), contribution_years=changes, difference_years=d[i])
  }))
}
write.csv(do.call(rbind,points),file.path(out,'paired_scores.csv'),row.names=FALSE)
write.csv(do.call(rbind,extremes),file.path(out,'plot_extremes.csv'),row.names=FALSE)
