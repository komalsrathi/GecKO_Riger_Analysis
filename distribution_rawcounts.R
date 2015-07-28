# distribution of read counts
setwd('/fujfs/d1/projects/Investigators/Zolt/gecko/GSNAP')

list.f <- list.files(pattern="*_UID_mappings.txt")
names <- sub('_UID_mappings.txt','',list.f)
summary <- read.csv('../bin/samples_fasta.csv',header=F,stringsAsFactors = F,col.names = c('Fasta','Barcode','Reference','Library'))
summary <- summary[,c(2,4)]
lib.total <- data.frame()

for(i in 1:20)
{
  lib <- read.csv(list.f[i])
  lib$name <- names[i]
  lib.total <- rbind(lib.total,lib)
}

lib.total <- merge(lib.total, summary, by.x='name', by.y='Barcode')

ggplot(data=lib.total, aes(x=log2(normcounts),color=name,linetype=Library)) + 
  geom_density() + xlab('\n log2(Normalized Counts)') + ggtitle('Distribution of Read Counts\n') + theme_bw(base_size = 18) + facet_wrap(~name) + 
  guides(color=FALSE)

