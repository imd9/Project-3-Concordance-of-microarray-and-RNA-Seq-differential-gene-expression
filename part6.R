#Author: Monica Roberts
#Purpose: part 6 of project 3. Computes concordance between microarray DE and RNA-Seq DE. Produces plots describing the results.

library(ggplot2)
library(dplyr)
library(biomaRt)
library(ggpubr)
library(RColorBrewer)

#probe id to refseq mapping csv
probe_id_mapping <- read.csv('/project/bf528/project_3/refseq_affy_map.csv') %>% rename(probeid=PROBEID, refseq=REFSEQ)

#read in files from part 5 and merge with mapping file
limma_beta <- read.csv('group4_beta_limma_results.csv', as.is=TRUE, header=TRUE) %>% 
  rename(probeid=X) 
limma_beta_mapped <- merge(limma_beta, probe_id_mapping, by='probeid') %>% arrange(SYMBOL, P.Value)

limma_bezafibrate <- read.csv('group4_bezafibrate_limma_results.csv', as.is=TRUE, header=TRUE) %>% 
  rename(probeid=X)
limma_beza_mapped <- merge(limma_bezafibrate, probe_id_mapping, by='probeid') %>% arrange(SYMBOL, P.Value)

limma_nitro <- read.csv('group4_nitro_limma_results.csv', as.is=TRUE, header=TRUE) %>% 
  rename(probeid=X)
limma_nitro_mapped <- merge(limma_nitro, probe_id_mapping, by='probeid') %>% arrange(SYMBOL, P.Value)

# only keep most significant mapped probes for each gene symbol then filter out significant probes at padj<0.05
limma_beta_mappedu <- limma_beta_mapped[match(unique(limma_beta_mapped$SYMBOL), limma_beta_mapped$SYMBOL),]
limma_beta_DE <- limma_beta_mappedu %>% filter(adj.P.Val<0.05) %>% arrange(P.Value)

limma_beza_mappedu <- limma_beza_mapped[match(unique(limma_beza_mapped$SYMBOL), limma_beza_mapped$SYMBOL),]
limma_beza_DE <- limma_beza_mappedu %>% filter(adj.P.Val<0.05) %>% arrange(P.Value)

limma_nitro_mappedu <- limma_nitro_mapped[match(unique(limma_nitro_mapped$SYMBOL), limma_nitro_mapped$SYMBOL),]
limma_nitro_DE <- limma_nitro_mappedu %>% filter(adj.P.Val<0.05) %>% arrange(P.Value)

#load in deseq results and use mapping file to map symbols, keep unique genes, filter for significant
deseq_beta_results <- read.csv('/projectnb/bf528/users/swiss_cheese_2022/project_3/programmer/part4/DESeq_ER_results.csv') %>% 
  rename(refseq=X)
deseq_beta_mapped <- merge(deseq_beta_results, probe_id_mapping, by='refseq') %>% arrange(SYMBOL, pvalue)
deseq_beta_mappedu<- deseq_beta_mapped[match(unique(deseq_beta_mapped$SYMBOL), deseq_beta_mapped$SYMBOL),]
deseq_beta_DE <- deseq_beta_mappedu %>% filter(padj<0.05) 

deseq_beza_results <- read.csv('/projectnb/bf528/users/swiss_cheese_2022/project_3/programmer/part4/DESeq_PPARA_results.csv') %>% 
  rename(refseq=X)
deseq_beza_mapped <- merge(deseq_beza_results, probe_id_mapping, by='refseq') %>% arrange(SYMBOL, pvalue)
deseq_beza_mappedu<- deseq_beza_mapped[match(unique(deseq_beza_mapped$SYMBOL), deseq_beza_mapped$SYMBOL),]
deseq_beza_DE <- deseq_beza_mappedu %>% filter(padj<0.05) 

deseq_nitro_results <- read.csv('/projectnb/bf528/users/swiss_cheese_2022/project_3/programmer/part4/DESeq_DNA_Damage_Rresults.csv') %>% 
  rename(refseq=X)
deseq_nitro_mapped <- merge(deseq_nitro_results, probe_id_mapping, by='refseq') %>% arrange(SYMBOL, pvalue)
deseq_nitro_mappedu<- deseq_nitro_mapped[match(unique(deseq_nitro_mapped$SYMBOL), deseq_nitro_mapped$SYMBOL),]
deseq_nitro_DE <- deseq_nitro_mappedu %>% filter(padj<0.05) 

# create list of unique genes analysed in whole comparison to compute N for correction below
beta_all_genes <- merge(deseq_beta_mappedu, limma_beta_mappedu, all.x=TRUE, all.y=TRUE) %>% distinct(SYMBOL)
beza_all_genes <- merge(deseq_beza_mappedu, limma_beza_mappedu, all.x=TRUE, all.y=TRUE) %>% distinct(SYMBOL)
nitro_all_genes <- merge(deseq_nitro_mappedu, limma_nitro_mappedu, all.x=TRUE, all.y=TRUE) %>% distinct(SYMBOL)

#create function to calculate corrected intersection
background_correct_intersect <- function(deseq_results, limma_results, total_genes) {
  n0 <- nrow(inner_join(limma_results, deseq_results, by='SYMBOL'))
  N = nrow(total_genes)
  n1 = nrow(limma_results)
  n2 = nrow(deseq_results)
  x = ((n0*N)-(n1*n2))/(n0+N-n1-n2)
  return(x)
}

#calculate corrected intersection for each group
beta_true_intersect <- background_correct_intersect(deseq_beta_DE, limma_beta_DE, beta_all_genes)
beza_true_intersect <- background_correct_intersect(deseq_beza_DE, limma_beza_DE, beza_all_genes)
nitro_true_intersect <- background_correct_intersect(deseq_nitro_DE, limma_nitro_DE, nitro_all_genes)

#concordnace calculation
concordance <- function(intersect, DE_microarray, DE_RNASeq) {
  concordance <- (2*intersect) / (nrow(DE_microarray) + nrow(DE_RNASeq))
  return(concordance)
}

# calculate concordance for each treatment group
beta_concordance <- concordance(beta_true_intersect, limma_beta_DE, deseq_beta_DE)
beza_concordance <- concordance(beza_true_intersect, limma_beza_DE, deseq_beza_DE)
nitro_concordance <- concordance(nitro_true_intersect, limma_nitro_DE, deseq_nitro_DE)

# plot for concordance and # DE genes for microarray
microarray_concordance <- tibble(chemical=c('Beta', 'Beza', 'Nitro'),
                                    concordance = c(beta_concordance, beza_concordance, nitro_concordance),
                                    DE_genes=c(nrow(limma_beta_DE), nrow(limma_beza_DE), nrow(limma_nitro_DE))
                                    )

ggplot(microarray_concordance, aes(x=DE_genes, y=concordance, label=chemical)) + 
  geom_point(size=3, shape=21) +
  theme_bw() +
  geom_text(hjust=0.6, vjust=1.15) +
  xlab('Treatment effect (# of DEGs from microarray)') +
  ylab('Concordance of DEGs') + 
  geom_smooth(method='glm', se=FALSE, linetype='dashed', color='black') +
  stat_regline_equation(label.y=0.5, aes(label=..rr.label..))

# plot for concordance and # DE genes for RNAseq
RNASeq_concordance <- tibble(chemical=c('Beta', 'Beza', 'Nitro'),
                                 concordance = c(beta_concordance, beza_concordance, nitro_concordance),
                                 DE_genes=c(nrow(deseq_beta_DE), nrow(deseq_beza_DE), nrow(deseq_nitro_DE))
)

ggplot(RNASeq_concordance, aes(x=DE_genes, y=concordance, label=chemical)) + 
  geom_point(size=3, shape=21) +
  theme_bw() +
  geom_text(hjust=0.6, vjust=1.5) +
  xlab('Treatment effect (# of DEGs from RNASeq)') +
  ylab('Concordance of DEGs') + 
  geom_smooth(method='glm', se=FALSE, linetype='dashed', color='black') +
  stat_regline_equation(label.y=0.5, label.x=3900, aes(label=..rr.label..))
  
 
#separate above and below median genes for each group
limma_beta_above <- limma_beta_DE %>% filter(AveExpr > median(AveExpr))
deseq_beta_above <- deseq_beta_DE %>% filter(baseMean > median(baseMean))

limma_beta_below <- limma_beta_DE %>% filter(AveExpr < median(AveExpr))
deseq_beta_below <- deseq_beta_DE %>% filter(baseMean < median(baseMean))


limma_beza_above <- limma_beza_DE %>% filter(AveExpr > median(AveExpr))
deseq_beza_above <- deseq_beza_DE %>% filter(baseMean > median(baseMean))

limma_beza_below <- limma_beza_DE %>% filter(AveExpr < median(AveExpr))
deseq_beza_below <- deseq_beza_DE %>% filter(baseMean < median(baseMean))


limma_nitro_above <- limma_nitro_DE %>% filter(AveExpr > median(AveExpr))
deseq_nitro_above <- deseq_nitro_DE %>% filter(baseMean > median(baseMean))

limma_nitro_below <- limma_nitro_DE %>% filter(AveExpr < median(AveExpr))
deseq_nitro_below <- deseq_nitro_DE %>% filter(baseMean < median(baseMean))

#re-calculate true intersection for each group
split_total_genes <- beta_all_genes %>% slice_head(n=6374)

beta_up_intersect <- background_correct_intersect(deseq_beta_above, limma_beta_above, split_total_genes)
beta_down_intersect <- background_correct_intersect(deseq_beta_below, limma_beta_below, split_total_genes)

beza_up_intersect <- background_correct_intersect(deseq_beza_above, limma_beza_above, split_total_genes)
beza_down_intersect <- background_correct_intersect(deseq_beza_below, limma_beza_below, split_total_genes)

nitro_up_intersect <- background_correct_intersect(deseq_nitro_above, limma_nitro_above, split_total_genes)
nitro_down_intersect <- background_correct_intersect(deseq_nitro_below, limma_nitro_below, split_total_genes)

#recalculate concordance for each treatment for above and below median
beta_above_concordance <- concordance(beta_up_intersect, limma_beta_above, deseq_beta_above)
beta_below_concordance <- concordance(beta_down_intersect, limma_beta_below, deseq_beta_below)

beza_above_concordance <- concordance(beza_up_intersect, limma_beza_above, deseq_beza_above)
beza_below_concordance <- concordance(beza_down_intersect, limma_beza_below, deseq_beza_below)

nitro_above_concordance <- concordance(nitro_up_intersect, limma_nitro_above, deseq_nitro_above)
nitro_below_concordance <- concordance(nitro_down_intersect, limma_nitro_below, deseq_nitro_below)

#create dataset combining all concordance scores
all_concordance <- tibble(chemical=c('Beta', 'Beza', 'Nitro', 'Beta', 'Beza', 'Nitro', 'Beta', 'Beza', 'Nitro'),
                          group=c('all', 'all', 'all', 'above', 'above', 'above', 'below', 'below', 'below'),
                          concordance=c(beta_concordance, beza_concordance, nitro_concordance,
                                        beta_above_concordance, beza_above_concordance, nitro_above_concordance,
                                        beta_below_concordance, beza_below_concordance, nitro_below_concordance))
all_concordance$group <- factor(all_concordance$group)

ggplot(all_concordance, aes(x=chemical, y=concordance, fill=group)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() 
