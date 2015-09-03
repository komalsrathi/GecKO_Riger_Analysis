library(ggplot2)

# distribution of read counts
setwd('/fujfs/d1/projects/Investigators/Zolt/gecko/GSNAP')

# get info file
summary <- read.csv('../bin/samples_fasta.csv',header=F,stringsAsFactors = F,col.names = c('Fasta','Barcode','Reference','Library'))

# get all normalized counts in one data frame
lib.total <- ddply(.data = summary, .variables = 'Barcode', .fun = function(x) concat.files(x))

concat.files <- function(x)
{
  lib <- read.csv(paste(x$Barcode,'_UID_mappings_normalized.txt',sep=""))
  lib$Library <- x$Library
  return(lib)
}

# distribution
ggplot(data=lib.total, aes(x=log2(normcounts),color=Barcode,linetype=Library)) + 
  geom_density() + xlab('\n log2(Normalized Counts)') + ggtitle('Distribution of Read Counts\n') + theme_bw(base_size = 18) + facet_wrap(~Barcode) + 
  guides(color=FALSE)

