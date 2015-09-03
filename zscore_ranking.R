library(data.table)
library(plyr)
library(dplyr)
library(preprocessCore)
source('../bin/calc_empirical_pval.R') # calculate empirical pvalues
source('../bin/quantile_normalize.R') # quantile normalize counts
source('../bin/calc_exp_zscore.R') # calculate expected zscore
source('../bin/calc_null_zscore.R') # calculate null 

# directory where counts are stored
setwd('~/Desktop/komalrclust/projects/Investigators/Zolt/gecko/GSNAP')

# read in sample info
sample <- read.csv('../sample_info_mitophagy.csv',stringsAsFactors = F)

# quantile normalize counts per library
norm.lib <- ddply(sample,.variables = c('library'),.fun = function(x) quantile.normalize(x))

# remove non targeting genes & predicted genes
norm.lib <- norm.lib[-grep('NonTargeting|Gm[0-9]+[b]*$',norm.lib$gene_id),]

# geometric means of lows & highs
norm.lib$meanLow = apply(norm.lib[,grep('Low',colnames(norm.lib))],MARGIN = 1,FUN = function(x) exp(mean(log(x))))
norm.lib$meanHigh = apply(norm.lib[,grep('High',colnames(norm.lib))],MARGIN = 1,FUN = function(x) exp(mean(log(x))))

# calculate arithmetic mean of meanLow & meanHigh
norm.lib$mean <- apply(norm.lib[grep('mean',colnames(norm.lib))],MARGIN = 1,mean)

# calculate fold change meanlow vs meanhigh
norm.lib$fc <- norm.lib$meanLow/norm.lib$meanHigh

# calculate expected zscores
z.exp <- calc.exp.zscore(norm.lib)

# calculate null zscore distribution
z.null <- do.call(rbind,lapply(1:100,function(y) calc.null.zscore(norm.dat = norm.lib)))

# calculate empirical pvalues
results <- calc.emp.pval(dat.null = z.null, dat.exp = z.exp)

# write.csv(z.exp,'../results/rankings.csv',quote=F,row.names=F)
