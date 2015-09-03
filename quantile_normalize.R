# quantile normalize raw counts
quantile.normalize <- function(x)
{
  files <- list.files(pattern = 'UID_mappings_normalized.txt') # get all files with UID counts
  y <- files[sapply(x$barcode, function(y) grep(y,files))] # get filenames matching only the current library
  dat <- as.data.frame(do.call("cbind", lapply(y, function(fn) data.frame(read.csv(fn,stringsAsFactors=FALSE))))) # read in files
  dat <- dat[,c(1,2,grep('rawcounts',colnames(dat)))] # retain first two columns and columns with counts
  names <- paste('m',x$sample,x$run,sep = '') # names to set for count columns based on sample & run
  dat[,3:ncol(dat)] <- dat[,3:ncol(dat)]+1 # add 1 to raw counts
  dat <- cbind(dat[,1:2],as.data.frame(normalize.quantiles(as.matrix(dat[,3:ncol(dat)])))) # quantile normalize 
  colnames(dat)[3:ncol(dat)] <- names # set names to count columns
  return(dat)
}
