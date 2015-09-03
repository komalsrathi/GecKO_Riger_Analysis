# this function does random shuffling of gene to guide RNAs
# input is a dataframe of UID, gene_id, quantile normalized counts for each sample along with means and fold change

calc.null.zscore <- function(norm.dat){
  # shuffle gene id
  # order data by mean(normalized counts) 
  dt <- transform(norm.dat, gene_id = sample(gene_id)) 
  dt <- dt[order(dt$mean),] 

  # calculate zscore
  # split in 12 bins, apply scale function on fold change (z-score) and combine the bins in dt2
  dt %>% mutate(grp=cut(1:nrow(dt),12)) %>%
    group_by(grp) %>%
    mutate(z.score=c(scale(fc))) %>%
    ungroup() %>%        
    select(-grp) %>%
    as.data.frame -> dt2 
  
  # fc_scaled is scaled(fc) i.e. z-score
  # get top 4 highest z-score per gene and calculate mean
  # this will give us a dataframe with gene_id & mean(z-score) for positive regulators
  dt2 %>%
    arrange(gene_id,desc(z.score)) %>% # sort by gene_id, then by z-score - highest to lowest
    mutate(ranks = ave(z.score, gene_id, FUN = seq_along)) %>% # for each gene, add ranks to z-score. Assigns ranks 1-6 for highest to lowest z-score per gene.
    group_by(gene_id) %>% # now group by gene_id
    filter(ranks %in% c(1:4)) %>% # get top 4 highest ranking z-scores per gene
    summarize(mean.zscore.pos=mean(z.score)) %>% # calculate mean of those z-scores -> collapsed to gene level
    as.data.frame -> dat.pos
  
  # similarly get top 4 lowest z-score per gene and calculate mean
  # this will give us a dataframe with gene_id & mean(z-score) for negative regulators
  dt2 %>%
    arrange(gene_id,z.score) %>% 
    mutate(ranks = ave(z.score, gene_id, FUN = seq_along)) %>% # for each gene, add ranks to z-score. Assigns ranks 1-6 for lowest to highest z-score per gene.
    group_by(gene_id) %>% 
    filter(ranks %in% c(1:4)) %>% 
    summarize(mean.zscore.neg=mean(z.score)) %>% 
    as.data.frame -> dat.neg
  
  # merge dat.pos and dat.neg
  res <- merge(dat.pos,dat.neg,by='gene_id')
  
  return(res)
}
