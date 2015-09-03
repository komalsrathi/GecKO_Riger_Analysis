# dat.null is shuffled data 
# dat.exp is real
calc.emp.pval <- function(dat.null,dat.exp){
  
  # dat.exp$mean.zscore.pos -> real z.scores for positive regulators
  dat.exp$rank.in.null.pos <- rank(-c(dat.exp$mean.zscore.pos,dat.null$mean.zscore.pos))[1:length(dat.exp$mean.zscore.pos)] - rank(-dat.exp$mean.zscore.pos)
  dat.exp$empirical.pval.pos <- dat.exp$rank.in.null.pos/length(dat.null$mean.zscore.pos)
  
  # dat.exp$mean.zscore.neg -> real z.scores for negative regulators
  dat.exp$rank.in.null.neg <- rank(c(dat.exp$mean.zscore.neg,dat.null$mean.zscore.neg))[1:length(dat.exp$mean.zscore.neg)] - rank(dat.exp$mean.zscore.neg)
  dat.exp$empirical.pval.neg <- dat.exp$rank.in.null.neg/length(dat.null$mean.zscore.neg)
  
  return(dat.exp)  
}
