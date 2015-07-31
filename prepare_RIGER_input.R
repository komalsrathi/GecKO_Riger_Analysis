setwd('Desktop/komalrclust//projects/Investigators/Zolt/gecko/GSNAP')

# sample info
sample <- read.csv('../sample_info.csv',stringsAsFactors = F)

fname <- paste('FGC1134_s_5_',sample$barcode,'_UID_mappings_normalized.txt',sep='')

names <- unique(paste('m',sample$samp,sample$name,sep=''))

liba <- read.csv('../mouse_gecko_libraries/mouse_geckov2_library_a_2.csv')
libb <- read.csv('../mouse_gecko_libraries/mouse_geckov2_library_b_1.csv')
lib <- rbind(liba,libb)

j = 1 
for(i in seq(from = 1,to = 8,by = 2))
{
  dat <- read.csv(fname[i])
  dat2 <- read.csv(fname[i+1])
  dat <- rbind(dat,dat2)
  lib <- merge(lib,dat[,c(1:2,5)],by=c('gene_id','UID'))
  colnames(lib)[ncol(lib)] <- names[j]
  j <- j+1
}
lib <- lib[,-3]

# counts
gene.count <- plyr::count(lib,vars='gene_id')

# low pass
low.pass <- lib
low.pass[low.pass==1] <- NA
low.pass <- na.omit(low.pass)

# leave 4
leave.four <- lib
leave.four$mean <- apply(leave.four[,3:6],MARGIN = 1,mean)
leave.four <- plyr::arrange(leave.four,gene_id,plyr::desc(mean))

leave.four <- within(leave.four, {
  s <- ave(mean, gene_id, FUN = seq_along)
})

leave.four <- leave.four[-which(leave.four$s>4),-c(7:8)]

new.gene.counts <- plyr::count(leave.four,vars='gene_id')

# make gct files
make.gct <- function(dt,nsample){
  nr <- nrow(dt)
  dat <- data.frame(V1=c('#1.2',nr),V2=c('',nsample),V3=c('',''),V4=c('',''),V5=c('',''),V6=c('',''),stringsAsFactors = F)
  colnames(dt)[1:2] <- c('DESCRIPTION','NAME')
  dt <- dt[,c(2,1,3:6)]
  colnames(dat) <- colnames(dt)
  dat <- rbind(dat,colnames(dt))
  dat <- rbind(dat,dt)
  return(dat)
}

leave.four.gct <- make.gct(dt = leave.four, nsample = 4)
low.pass.gct <- make.gct(dt = low.pass, nsample = 4)
lib.gct <- make.gct(dt = lib, nsample = 4)

write.table(leave.four.gct,'../results/leave_four.gct',sep='\t',quote=F,row.names = F,col.names = F)
write.table(low.pass.gct,'../results/low_pass.gct',sep='\t',quote=F,row.names = F,col.names = F)
write.table(lib.gct,'../results/no_filter.gct',sep='\t',quote=F,row.names = F,col.names = F)
