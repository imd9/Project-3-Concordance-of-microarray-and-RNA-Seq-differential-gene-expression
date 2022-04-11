# Creates boxplot of count distribution for samples
# Clean R session:
rm(list=ls())
# Import libraries:
library(tidyverse)

# Read in csv as tibbles:
ds1 <- as_tibble(read.csv("featureCounts_combined.csv", header=TRUE))

# Create boxplot:
# Create boxplot, with labels horizontal (1) and increased left margin
jpeg("boxplot_panel.jpg")
boxplot(ds1[,c(2:10)], horizontal=TRUE, ylim=c(0, 1000), las=1, 
        pars=list(par(mar=c(4,7,4,4))),
        col=c("blue", "blue", "red", "red", "red", "green", "green", "green", "blue"),
        main="Distribution of Gene Counts", xlab="Gene Counts") 
legend("bottomright", title = "Mode of action" , legend=c("ER", "DNA DAMAGE", "PPARA"), 
       fill=c("green", "red", "blue"), cex=0.75)
dev.off()


# Produce boxplot on log-scale
# Add 0.0001 to each cell so log scale will work:
ds_shifted <- ds1
ds_shifted[,c(2:10)] <- ds1[,c(2:10)] + 0.0001

# Create boxplot logscale
jpeg("boxplot_panel_log.jpg")
boxplot(ds_shifted[,c(2:10)], horizontal=TRUE, las=1, log="x", 
        ylim=c(0.9, 1000000), pars=list(par(mar=c(4,7,4,4))),
        col=c("blue", "blue", "red", "red", "red", "green", "green", "green", "blue"),
        main="Distribution of Gene Counts", xlab="Gene Counts (log)")
dev.off()
