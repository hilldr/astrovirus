## R script to generate Wald and LRT test results from DESeq2 output
## David R. Hill
## -----------------------------------------------------------------------------

## enable parallel processes
library("BiocParallel")
register(MulticoreParam(4))
results.dir <- "../results/DESeq2/"

## load output of DESeq-export-counts.R
load(file = file.path(results.dir, "dds.Rdata"))

## setup multifactor design
DESeq2::design(dds) <- ~ group

## Run DESeq
dds <- DESeq2::DESeq(dds)

## Apply Wald test for specific comparisons
## 'test = "Wald"' is default; specified here for clarity
res <- DESeq2::results(dds, test = "Wald",
                       contrast = c("group", "V0", "M0"))
write.csv(res, file = file.path(results.dir, "V0_over_M0_Wald-test.csv"))
res <- DESeq2::results(dds, test = "Wald",
                       contrast = c("group", "V12", "M12"))
write.csv(res, file = file.path(results.dir, "V12_over_M12_Wald-test.csv"))
res <- DESeq2::results(dds, test = "Wald",
                       contrast = c("group", "V24", "M24"))
write.csv(res, file = file.path(results.dir, "V24_over_M24_Wald-test.csv"))
