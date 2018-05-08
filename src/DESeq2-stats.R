## R script to generate Wald and LRT test results from DESeq2 output
## David R. Hill
## -----------------------------------------------------------------------------

## enable parallel processes
library("BiocParallel")
library('magrittr')
register(MulticoreParam(4))
results.dir <- "../results/DESeq2/"

## load output of DESeq-export-counts.R
load(file = file.path(results.dir, "dds.Rdata"))

## setup multifactor design
DESeq2::design(dds) <- ~ group

## Run DESeq
dds <- DESeq2::DESeq(dds)

## Add annotation - symbol, entrezID, and gene name
annotateRes <- function(res){
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  res$symbol <- rownames(res)
  res$entrez <- mapIds(org.Hs.eg.db,
                       keys = rownames(res),
                       column = 'ENTREZID',
                       keytype = 'SYMBOL',
                       multiVals = 'first')
  res$name <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = 'GENENAME',
                     keytype = 'SYMBOL',
                     multiVals = 'first')
  return(res)
}

## Apply Wald test for specific comparisons
## 'test = "Wald"' is default; specified here for clarity
res <- DESeq2::results(dds, test = "Wald",
                       contrast = c("group", "V0", "M0")) %>% annotateRes
write.csv(res, file = file.path(results.dir, "V0_over_M0_Wald-test.csv"))
res <- DESeq2::results(dds, test = "Wald",
                       contrast = c("group", "V12", "M12")) %>% annotateRes
write.csv(res, file = file.path(results.dir, "V12_over_M12_Wald-test.csv"))
res <- DESeq2::results(dds, test = "Wald",
                       contrast = c("group", "V24", "M24")) %>% annotateRes
write.csv(res, file = file.path(results.dir, "V24_over_M24_Wald-test.csv"))
