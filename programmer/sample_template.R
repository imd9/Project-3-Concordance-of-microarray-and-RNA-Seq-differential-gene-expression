#R

library(DESeq2)

# load counts
cnts <- read.csv('groups/group_EX_rna_counts.csv',row.names=1)

# filter out rows that have any zeros for funzies
cnts <- subset(cnts,rowSums(cnts==0)==0)

# sample information
info <- read.csv('groups/group_EX_rna_info.csv')

# create the DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = cnts,
  colData = info,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds$mode_of_action <- relevel(dds$mode_of_action, ref='Control')

# run DESeq
dds <- DESeq(dds)
res <- results(dds, contrast=c('mode_of_action','HMGCOA','Control'))
res <- lfcShrink(dds, coef=2)

# write out DE results
write.csv(res,'example_deseq_results.csv')

# write out matrix of normalized counts
write.csv(counts(dds,normalized=TRUE),'example_deseq_norm_counts.csv')
