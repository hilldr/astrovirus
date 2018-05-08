#------------------------------------------------------------
#                RNA seq analysis workflow
#                    4_summary_table
#                   Ryan Berger 5-7-18
#
#   RNA seq workflow starting with kallisto aligned reads
#
#  Purpose: Make summary of gene differential expression
#  Input: Differential expression results dataframe (res.df)
#  Output: Table summary of gene changes
#
#
# Source: https://www.bioconductor.org/help/workflows/rnaseqGene/
#------------------------------------------------------------

#------------------------------------------------------------
#                     summary_table

## Function to make RNA seq summary table
## Input is: res dataframe, p value cutoff, log2fold change cutoff
GeneTable <- function(x, p = 0.05, fc = 1){
  library(dplyr)
  library(magrittr)
  rowIDs <- c('Total genes with counts',
              paste('Adjusted p value <', p),
              paste('Adj. p <', p, '& log2 Fold Change >', fc))
  Total <- c(nrow(x), ## Total column
             nrow(filter(x, padj < p)), 
             nrow(filter(x, padj < p & abs(log2FoldChange) > fc)))
  Increasing <- c(nrow(filter(x, log2FoldChange > 0)), ## Increasing column
                  nrow(filter(x, padj < p & log2FoldChange > 0)), 
                  nrow(filter(x, padj < p & log2FoldChange > fc)))
  Decreasing <- c(nrow(filter(x, log2FoldChange < 0)), ## Decreasing column
                  nrow(filter(x, padj < p & log2FoldChange < 0)), 
                  nrow(filter(x, padj < p & log2FoldChange < -fc)))
  ## Make table
  subsetTable <- data_frame(Total, Increasing, Decreasing) %>% 
    as.data.frame()
  rownames(subsetTable) <- rowIDs
  return(subsetTable)
}


## Make results tables for each differential expression file
results.dir <- "../results/DESeq2/"
library(gridExtra)

## 0h results
res <- read.csv(file = file.path(results.dir, 'V0_over_M0_Wald-test.csv'))
GeneSummary <- GeneTable(res, p = 0.05, fc =  1)
## Save png of table
png(filename = file.path(results.dir, "V0_over_M0_Gene-Summary.png"),
    height = 1.5, width = 6, units = 'in', res = 500)
grid.arrange(tableGrob(GeneSummary), 
             top = 'Virus over Mock 0h')
dev.off()


## 12h results
res <- read.csv(file = file.path(results.dir, 'V12_over_M12_Wald-test.csv'))
GeneSummary <- GeneTable(res, p = 0.05, fc = 1)
## Save png of table
png(filename = file.path(results.dir, "V12_over_M12_Gene-Summary.png"),
    height = 1.5, width = 6, units = 'in', res = 500)
grid.arrange(tableGrob(GeneSummary), 
             top = 'Virus over Mock 12h')
dev.off()


## 24h results
res <- read.csv(file = file.path(results.dir, 'V24_over_M24_Wald-test.csv'))
GeneSummary <- GeneTable(res, p = 0.05, fc = 1)
## Save png of table
png(filename = file.path(results.dir, "V24_over_M24_Gene-Summary.png"),
    height = 1.5, width = 6, units = 'in', res = 500)
grid.arrange(tableGrob(GeneSummary), 
             top = 'Virus over Mock 24h')
dev.off()