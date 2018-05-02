## Wald test results
results/DESeq2/V0_over_M0_Wald-test.csv \
results/DESeq2/V12_over_M12_Wald-test.csv \
results/DESeq2/V24_over_M24_Wald-test.csv : results/DESeq2/dds.Rdata
	R -e "setwd('./src/'); source('DEseq2-stats.R')"
