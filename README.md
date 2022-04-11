# Project 3: Concordance of microarray and RNA-Seq differential gene expression

Microarrays have been at the forefront of analyzing transcriptomes for evaluating drug safety. It helps in identifying Differentially Expressed Genes (DEGs) and predicting patient/toxicity outcomes based on gene-expression data. With the advent of newer technologies like high throughput sequencing technologies provide new methods for whole-transcriptome analyses of gene expression like RNA-seq. Some studies go on to compare the technical reproducibility of the results obtained from microarrays and RNA-seq experiments. Some have reported that RNA-seq reports lower precision for weakly expressed genes owing to the nature of sampling or higher sensitivity of RNA-seq for gene detection. Based on these varied conclusions, Wang et al conducted a comprehensive study to evaluate RNA-seq in its differences and similarities to microarrays in terms of identifying DEGs and developing predictive models. The design of this study was to generate Illumina RNA-seq and Affymetrix microarray data from the same set of liver samples of rats under varying degrees of perturbation by 27 chemicals representing multiple modes of action (MOA). 

Reference:
Wang, Charles, Binsheng Gong, Pierre R. Bushel, Jean Thierry-Mieg, Danielle Thierry-Mieg, Joshua Xu, Hong Fang, et al. 2014. “A comprehensive study design reveals treatment- and transcript abundance–dependent concordance between RNA-seq and microarray data” Nature Biotechnology 32 (9): 926–32.

# Contributors

Monica Roberts: Analyst (monicapr@bu.edu)

Preshita Dave: Data curator (preshita@bu.edu)

Italo Duran: Programmer (duran01@bu.edu)

# Repository Contents

1. STAR.qsub - written by Preshita Dave. Executes the alignment used by the STAR module and outputs the aligned bam files and alignment statistics. 
2. multiqc.qsub - written by Preshita Dave. Inputs the MultiQC module which summarizes the FastQC and STAR output files to generate an html file that displays overall statistics. 
3. run_limma.R - written by Monica Roberts. Executes limma differential expression analysis for each treatment group on the RMA expression matrix provided by the authors. Writes the output to csv files that are used in part 6 for concordance calculations. Also plots statistics from limma analysis.
4. part6.R - written by Monica Roberts. Takes limma DE results and DESeq2 results and does analysis on concordance between the two methods. Calculates true intersection, concordance score, and generates plots from the analysis. 
