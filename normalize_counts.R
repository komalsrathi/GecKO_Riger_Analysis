library(reshape2)

setwd('/fujfs/d1/projects/Investigators/Zolt/gecko/GSNAP')

sample <- read.csv('../bin/samples_fasta.csv',header=F,stringsAsFactors = F)

for(i in 1:nrow(sample))
{
  fname <- paste(sample[i,2],'UID_mappings.txt',sep='_')
  uid <- read.csv(file = fname, header=F, stringsAsFactors = F,col.names = c('ID','rawcounts'))
  uid <- cbind(uid,colsplit(string = as.character(uid$ID),pattern = "\\|",names = c('UID','gene_id')))
  
  if(sample[i,4]=='A'){
    lib <- read.csv(file = '../mouse_gecko_libraries/mouse_geckov2_library_a_2.csv')
  } else if(sample[i,4]=='B'){
    lib <- read.csv(file = '../mouse_gecko_libraries/mouse_geckov2_library_b_1.csv')
  }
  
  uid.lib <- merge(lib,uid[,2:ncol(uid)],by=c('gene_id','UID'),all.x=T)
  uid.lib$rawcounts[is.na(uid.lib$rawcounts)] <- 0
  total <- sum(uid.lib$rawcounts)
  uid.lib$normcounts <- ((uid.lib$rawcounts/total)*10e6+1)
  
  write.csv(uid.lib, file = fname, quote=F, row.names = F)
}