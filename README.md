# GecKO_Riger_Analysis

# RIGER analysis
normalize_counts_bytotalreads.R - normalize reads by total count (formula - rawcounts/totalreads*10^6+1)

distribution_normalized_counts_bytotalreads.R - get distribution of normalized counts (formula - rawcounts/totalreads*10^6+1)

prepare_RIGER_input.R - process normalized counts, make three files : nofilter, leave four and low pass. Then prepare GCT files for RIGER input.

# Adapter cutting & alignment using GSNAP
gecKO_batch.pl - gsnap pipeline for gecKO

# Z-score calculation pipeline
quantile_normalize.R - quantile normalize raw counts

calc_exp_zscore.R - function to calculate real zscores

calc_null_zscore.R - function to calculate null distribution of zscores

calc_empirical_pval.R - function to calculate empirical pvalues

zscore_ranking.R - full script using various functions to calculate zscores & emp. pvalues

# using the ranking file obtained in z-score calculation pipeline make GSEA ranking lists
make_ranked_list_GSEA.R - make ranked list of positive & negative regulators for GSEA

# converts fastq to fasta
fastq_to_fasta.pl
