# read in results from z-score analysis (zscore_ranking.R)
z.exp <- read.csv('../results/rankings.csv')

# make ranked list of positive regulators
ranked.list.pos <- z.exp[order(z.exp$mean.zscore.pos,decreasing = T),c(1,2)]
ranked.list.pos$gene_id <- toupper(ranked.list.pos$gene_id)
write.table(ranked.list.pos,'../results/GSEA/ranked_list_pos.rnk',quote=F,row.names = F,col.names = F,sep='\t')

# make ranked list of negative regulators
ranked.list.neg <- z.exp[order(z.exp$mean.zscore.neg,decreasing = T),c(1,3)]
ranked.list.neg$gene_id <- toupper(ranked.list.neg$gene_id)
write.table(ranked.list.neg,'../results/GSEA/ranked_list_neg.rnk',quote=F,row.names = F,col.names = F,sep='\t')
