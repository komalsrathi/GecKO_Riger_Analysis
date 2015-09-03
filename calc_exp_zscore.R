calc.exp.zscore <- function(norm.dat){
  
  # sort by mean
  norm.dat <- norm.dat[order(norm.dat$mean),]
  
  # calculate zscore
  # split, apply and combine
  norm.dat %>% mutate(grp=cut(1:nrow(norm.dat),12)) %>%
    group_by(grp) %>%
    mutate(z.score=c(scale(fc))) %>%
    ungroup() %>%        
    select(-grp) %>%
    as.data.frame -> dt2 
  
  # positive regulators
  # top ranked
  dt2 %>%
    arrange(gene_id,desc(z.score)) %>% # sort by gene_id, then by fc_scaled - highest to lowest
    mutate(ranks = ave(z.score, gene_id, FUN = seq_along)) %>% # add ranks to fc_scaled
    group_by(gene_id) %>% # group by gene_id
    filter(ranks %in% c(1:4)) %>% # get top 4 highest ranks
    summarize(mean.zscore.pos=mean(z.score)) %>%
    as.data.frame -> z.exp.top
  
  # negative regulators
  # bottom ranked
  dt2 %>%
    arrange(gene_id,z.score) %>% # sort by gene_id, then by fc_scaled - lowest to highest
    mutate(ranks = ave(z.score, gene_id, FUN = seq_along)) %>% # add ranks to fc_scaled 
    group_by(gene_id) %>% # group by gene_id
    filter(ranks %in% c(1:4)) %>% # get top 4 lowest ranks
    summarize(mean.zscore.neg=mean(z.score)) %>%
    as.data.frame -> z.exp.bottom
  
  # positive regulators obtained with real data are stored in z.exp.top
  # negative regulators obtained with real data are stored in z.exp.bottom
  z.exp <- merge(z.exp.top,z.exp.bottom,by='gene_id')
  
  return(z.exp)
}
