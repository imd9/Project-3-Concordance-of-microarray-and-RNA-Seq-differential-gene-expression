###################################################
#PART 4 ######

if (FALSE){ 
  install.packages("BiocManager")
  BiocManager::install("DESeq2")
}# To install packages switch to TRUE if needed....

library(tidyverse)
library(DESeq2)
library(stringr)
library("apeglm")
library(hardhat)
library("ggplot2")
setwd("/projectnb/bf528/users/swiss_cheese_2022/project_3/programmer/part4")

######4.1 ###############################################
#controls
controls <- read.table('/project/bf528/project_3/samples/control_counts.csv', header = TRUE, sep = ',')
info <- read.csv('/project/bf528/project_3/groups/group_4_rna_info.csv')
controls_subset <- controls[,info$Run[info$mode_of_action=='Control']]
controls_subset$Geneid <- controls$Geneid
# treatment counts
treatment_counts <- read.csv('/projectnb/bf528/users/swiss_cheese_2022/project_3/programmer/part3/featurecounts/featureCounts_combined.csv')
#combine controls with treatments
treat_cont_comb <- merge(controls_subset, treatment_counts, by='Geneid') %>% column_to_rownames(var='Geneid')


######4.2 ###############################################
#run deseq
# filter out rows that have any zeros for funzies
cts_treat_control<- subset(treat_cont_comb,rowSums(treat_cont_comb==0)==0)
colnames(cts_treat_control) <- gsub("\\.", "", colnames(cts_treat_control))
# sample information
info1 <- info[info$vehicle == 'SALINE_100_%',] %>% 
  subset(mode_of_action == 'DNA_Damage' | mode_of_action == 'Control')
dna_subset <- cts_treat_control[,info1$Run]

info2 <- info[info$vehicle == 'CORN_OIL_100_%',] %>%
  subset(mode_of_action == 'ER' | mode_of_action == 'Control')
er_subset<- cts_treat_control[,info2$Run]

info3 <- info[info$vehicle == 'CORN_OIL_100_%',] %>% 
  subset(mode_of_action == 'PPARA' | mode_of_action == 'Control')
ppara_subset<- cts_treat_control[,info3$Run]

#cts_tc_1 <- cts_treat_control %>% select(dput(as.character(info1$Run))) 
#cts_tc_2 <- cts_treat_control %>% select(dput(as.character(info2$Run)))                                                   
#cts_tc_3 <- cts_treat_control %>% select(dput(as.character(info3$Run))) 

#4.3####################################
# create the DESeq object with info1 

dds1 <- DESeqDataSetFromMatrix(
  countData = dna_subset,
  colData = info1,
  design= ~ mode_of_action
)

dds2 <- DESeqDataSetFromMatrix(
  countData = er_subset,
  colData = info2,
  design= ~ mode_of_action
)

dds3 <- DESeqDataSetFromMatrix(
  countData = ppara_subset,
  colData = info3,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds1$mode_of_action <- relevel(dds1$mode_of_action, ref='Control')
dds2$mode_of_action <- relevel(dds2$mode_of_action, ref='Control')
dds3$mode_of_action <- relevel(dds3$mode_of_action, ref='Control')


# run DESeq
dds1 <- DESeq(dds1) #saline_100_% vehicle
dds2 <- DESeq(dds2) #CORN_OIL_100_% vehicle                       
dds3 <- DESeq(dds3)
#lfcShrink(dds, 2)

#4.4 ################################
#against dnadamage mode of action
res_DNA_Damage <- results(dds1, contrast=c('mode_of_action','DNA_Damage','Control'))
res_DNA_Damage <- lfcShrink(dds1, 2)
write.csv(res_DNA_Damage, 'DESeq_DNA_Damage_Rresults.csv')
write.csv(res_DNA_Damage[order(res_DNA_Damage$padj),], 'DNA_Damage_resultsPADJ.csv') #padj sorted

# lenght of the DE p-ajusted value:
length(which(res_DNA_Damage$padj < 0.05)) #number of significant genes at threshold #4369

####to get the top 10 #########
#DNA_Damage_sig <- subset(as_tibble(res_DNA_Damage, rownames = 'genes'), padj < 0.05) %>%
  #column_to_rownames('genes')
#hist(DNA_Damage_sig$log2FoldChange)
#ggplot(DNA_Damage_sig, aes(log2FoldChange, -log10(pvalue))) + geom_point()
#write.csv(top_n(res_DNA_Damage, 10))

#against ER mode of action
res_er <- results(dds2, contrast=c('mode_of_action','ER','Control'))
res_er <- lfcShrink(dds2, 2)
write.csv(res_er, 'DESeq_ER_results.csv')
write.csv(res_er[order(res_er$padj),], 'ER_resultsPADJ.csv') #padj sorted

# lenght of the DE p-ajusted value:
length(which(res_er$padj < 0.05))
#3252

#against PPARA mode of action
res_PPARA <- results(dds3, contrast=c('mode_of_action','PPARA','Control'))
res_PPARA <- lfcShrink(dds3, 2)
write.csv(res_PPARA, 'DESeq_PPARA_results.csv')
write.csv(res_PPARA[order(res_PPARA$padj),], 'PPARA_resultsPADJ.csv')

# lenght of the DE p-ajusted value:
length(which(res_PPARA$padj < 0.05))
#2834


######################################
#4.5 Report the top 10 DE genes from each analysis by p-value in a table

res1_sorted = res_DNA_Damage[order(res_DNA_Damage$padj),]
res1_sorted = res1_sorted[complete.cases(res1_sorted), ]
res2_sorted = res_er[order(res_er$padj),]
res2_sorted = res2_sorted[complete.cases(res2_sorted), ]
res3_sorted = res_PPARA [order(res_PPARA$padj),]
res3_sorted = res3_sorted[complete.cases(res3_sorted), ]

DNA_Damage_sig_genecounts = dim(res1_sorted[res1_sorted$padj<0.05,])[1]
DNA_Damage_top10 = row.names(res1_sorted)[1:10]
ER_sig_genecounts = dim(res2_sorted[res2_sorted$padj<0.05,])[1]
ER_top10 = row.names(res2_sorted)[1:10]
PPARA_sig_genecounts = dim(res3_sorted[res3_sorted$padj<0.05,])[1]
PPARA_top10 = row.names(res3_sorted)[1:10]

Result=data.frame(col1=c(DNA_Damage_sig_genecounts,DNA_Damage_top10),col2=c(ER_sig_genecounts,ER_top10),col3=c(PPARA_sig_genecounts,PPARA_top10))
names(Result) = c("DNA_Damage", "ER", "PPARA")
rownames(Result) = c("Significant genes count", "1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th", "9th", "10th")
write.csv(Result, "Toxgroup significant gene count and top 10 most significant genes.csv")

################################################################
#4.6 Histogram
# Histograms of FC Values
graph_df <- as.data.frame(res1_sorted) %>% mutate(Group = "DNA_Damage") %>% 
  bind_rows(as.data.frame(res2_sorted) %>% mutate(Group = "ER"), 
            as.data.frame(res3_sorted) %>% mutate(Group = "PPARA")) %>%
  filter(padj < 0.05)

graph_df %>%
  ggplot() + 
  geom_histogram(aes(x = log2FoldChange, fill = Group), bins = 50) + 
  facet_wrap(~ Group) + 
  labs(title = "Log2 Fold Change of Significant DE Genes", 
       x = "Log2 Fold Change", y = "Count") + 
  theme(legend.position="none")

ggsave("/projectnb/bf528/users/swiss_cheese_2022/project_3/programmer/part4/DEhist.png")

graph_df %>%
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = padj, color = Group), alpha = 0.1) + 
  facet_wrap(~ Group) + 
  labs(title = "Log2 Fold Change vs. Adjusted p-value", x = "Log2 Fold Change", 
       y = "Adjusted p-value") + 
  theme(legend.position="none")

ggsave("/projectnb/bf528/users/swiss_cheese_2022/project_3/programmer/part4/log2foldchange.png")  


##########################################################
# Volcano Plots
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

png("/projectnb/bf528/users/swiss_cheese_2022/project_3/programmer/part4/volcano_DNA_Damage.png", width = 800, height = 480)
EnhancedVolcano(as.data.frame(res1_sorted),
                lab = rownames(res1_sorted),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'DNA_Damage Differentially Expressed Genes',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                legendPosition = 'right')
dev.off()

png("/projectnb/bf528/users/swiss_cheese_2022/project_3/programmer/part4/volcano_ER.png", width = 800, height = 480)
EnhancedVolcano(as.data.frame(res2_sorted),
                lab = rownames(res2_sorted),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'ER Differentially Expressed Genes',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                legendPosition = 'right')
dev.off()

png("/projectnb/bf528/users/swiss_cheese_2022/project_3/programmer/part4/volcano_ppara.png", width = 800, height = 480)
EnhancedVolcano(as.data.frame(res3_sorted),
                lab = rownames(res3_sorted),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'PPARA Differentially Expressed Genes',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                legendPosition = 'right')
dev.off()

# 4.7 ################################
# Extract the normalized counts out of the DESeq2 object, and save the normalized counts matrix to a file.
write.csv(counts(dds1,normalized=TRUE),'DNA_Damage_deseq_norm_counts.csv')
write.csv(counts(dds2,normalized=TRUE),'ER_deseq_norm_counts.csv')
write.csv(counts(dds3,normalized=TRUE),'PPARA_deseq_norm_counts.csv')


##########################################################################
