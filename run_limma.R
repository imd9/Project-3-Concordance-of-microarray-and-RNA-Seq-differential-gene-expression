# Author: Monica Roberts

# Purpose: Part 5 of project. Run limma on microarray expression matrix to determine differentially expressed
# genes. Plot logFC distribution and relationship between p-value and log FC per treatment.
# Write results out to CSV.

library(limma)
library(dplyr)
library(ggplot2)

# sample info dataframe with array_id and chemical columns
samples <- read.csv('group_4_mic_info.csv',as.is=TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('/project/bf528/project_3/samples/liver-normalization-rma.txt',
  sep='\t',
  as.is=TRUE,
  header=TRUE,
  row.names=1
)

# subset the full expression matrix to just those in this comparison
rma.subset <- rma[paste0('X',samples$array_id)]

# construct a design matrix modeling treatment vs control for use by limma
design <- model.matrix(
  ~factor(
    samples$chemical,
    levels=c('Control', 'BETA-ESTRADIOL', 'BEZAFIBRATE', 'N-NITROSODIETHYLAMINE')
  )
)
colnames(design) <- c('Intercept','BETA-ESTRADIOL', 'BEZAFIBRATE', 'N-NITROSODIETHYLAMINE')

# run limma
fit <- lmFit(rma.subset, design)
fit <- eBayes(fit)
t1 <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH') %>% arrange(adj.P.Val)
t2 <- topTable(fit, coef=3, n=nrow(rma.subset), adjust='BH') %>% arrange(adj.P.Val)
t3 <- topTable(fit, coef=4, n=nrow(rma.subset), adjust='BH') %>% arrange(adj.P.Val)

# write out the results to file
write.csv(t1,'group4_beta_limma_results.csv')
write.csv(t2,'group4_bezafibrate_limma_results.csv')
write.csv(t3,'group4_nitro_limma_results.csv')

#histogram of fold change values from significant DE
t1_sig <- t1 %>% filter(adj.P.Val<0.05)
t2_sig <- t2 %>% filter(adj.P.Val<0.05)
t3_sig <- t3 %>% filter(adj.P.Val<0.05)

t1_sig %>% ggplot(aes(x=logFC)) + 
  geom_histogram(fill='coral4', bins = 30, color='white') + 
  theme_bw() +
  ggtitle('Distribution of Log Fold Change Values for Beta-Estradiol')

t2_sig %>% ggplot(aes(x=logFC)) + 
  geom_histogram(fill='rosybrown', bins = 30, color='white') + 
  theme_bw() +
  ggtitle('Distribution of Log Fold Change Values for Bezafibrate') +
  xlim(-3,5)

t3_sig %>% ggplot(aes(x=logFC)) + 
  geom_histogram(fill='lightsalmon1', bins = 20, color='white') + 
  theme_bw() +
  ggtitle('Distribution of Log Fold Change Values for N-Nitrosodiethylamine')

t1_sig %>% ggplot(aes(y=-log10(P.Value), x=logFC)) + 
  geom_point(color='coral4') +
  theme_bw() +
  ggtitle('Nominal P-value vs. Log FC for Beta-Estradiol DE Genes')

t2_sig %>% ggplot(aes(y=-log10(P.Value), x=logFC)) + 
  geom_point(color='rosybrown') +
  theme_bw() +
  ggtitle('Nominal P-value vs. Log FC for Bezafibrate DE Genes')

t3_sig %>% ggplot(aes(y=-log10(P.Value), x=logFC)) + 
  geom_point(color='lightsalmon1') +
  theme_bw() +
  ggtitle('Nominal P-value vs. Log FC for Bezafibrate DE Genes')
