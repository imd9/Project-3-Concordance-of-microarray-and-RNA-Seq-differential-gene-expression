In this folder for the programmer you will find all the deliverables and the plots that are needed to understand this project.

## Part 3 of project,

1st and 2nd part for the programmer was to take the bam files and use the use the featurecounts.qsub file that was buit to get reports and in a new formart.

3rd part, use the built multiqc.qsub file to get a report from the 9 files from the feature counts.

4th,Combine the counts files from each sample into a single comma-delimited text file where the first column is the gene name and 
subsequent columns are the counts taken from each file for each sample. 
Create a header row in the combined file with “gene_id” in the first column and the corresponding sample name in subsequent columns. 
You might consider doing this using the read.csv function and dataframes in R.

5th, Create box plots for each of your samples showing the distribution of counts. 
Report your observations about distributional differences between your samples.

## Part 4 of project,

1st, Create a new counts matrix that includes control counts. Your sample metadata file group_N_rna_info.csv contains 
the sample IDs of the control samples for vehicles that correspond to your treatment samples. Identify the columns of 
your control samples in the sample matrix provided at /project/bf528/project_3/samples/control_counts.csv. 
You will need to subset columns out of the control matrix and merge them with your treatment counts matrix. 
NB: Be sure you are matching rows of the treatment and control samples so that the genes match.

2nd, Install the DESeq2 package from Bioconductor, read the vignette to learn how to use it.
We have provided an example script illustrating DESeq2 usage in /project/bf528/project_3/scripts/run_deseq.R to help you.

3rd, Write an R script to run DESeq2 comparing each group of your samples to the controls in the combined data matrix from 4.1. 
This will produce three separate DE gene lists, one for each condition with a corresponding set of controls. 
You will need to modify your script for each comparison; you may consider writing separate scripts for each analysis or put them all in the same script. 
You will need to subset your counts matrix to include only the relevant samples for a given comparison before running DESeq2. 
The appropriate controls for a given treatment are the ones that have the same vehicle value.

4th, Write out the differential expression results to files sorted by adjusted p-value. Report the number of genes significant at p-adjust < 0.05.

5th, Report the top 10 DE genes from each analysis by p-value in a table.

6th, Create histograms of fold change values from the significant DE genes for each analysis.
Also create scatter plots of fold change vs nominal p-value.

7th, 4 person groups only: Read the DESeq2 manual to identify how to extract the normalized counts out of the DESeq2 object, 
and save the normalized counts matrix to a file. You will use this matrix for clustering in part 7.


