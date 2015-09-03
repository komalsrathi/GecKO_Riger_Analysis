library(reshape2)

setwd('/fujfs/d1/projects/Investigators/Zolt/gecko/GSNAP')

sample <- read.csv('../bin/samples_fasta.csv',header=F,stringsAsFactors = F)

# apply to each sample
d_ply(.data = sample,.variables = 'V2',.fun = function(x) norm.by.total.reads(x))

norm.by.total.reads <- function(x){
  y <- paste(x$V2,'UID_mappings.txt',sep='_')
  uid <- read.csv(file = y, header=F, stringsAsFactors = F,col.names = c('ID','rawcounts')) # get raw counts
  uid <- cbind(uid, colsplit(string = as.character(uid$ID), pattern = "\\|", names = c('UID','gene_id'))) # split ID into UID and gene_id
  
  # get corresponding library of guide RNAs to add guide RNAs not reported by htseq-count
  if(length(grep('MGLibA',uid$UID))){
    lib <- read.csv(file = '../mouse_gecko_libraries/mouse_geckov2_library_a_2.csv')
  } else if(length(grep('MGLibB',uid$UID))){
    lib <- read.csv(file = '../mouse_gecko_libraries/mouse_geckov2_library_b_1.csv')
  }
  
  # normalize by dividing by total read counts & multiplying by 10^6
  uid <- merge(lib,uid[,2:ncol(uid)],by=c('gene_id','UID'),all.x=T)
  uid$rawcounts[is.na(uid$rawcounts)] <- 0
  total <- sum(uid$rawcounts)
  uid$normcounts <- ((uid$rawcounts/total)*10^6+1)
  
  # write out the normalized counts
  fileout = paste(x$V2,'UID_mappings_normalized.txt',sep='_')
  write.csv(uid, file = fileout, quote=F, row.names = F)
}
