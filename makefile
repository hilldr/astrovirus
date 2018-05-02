## Wald test results
results/DESeq2/V0_over_M0_Wald-test.csv \
results/DESeq2/V12_over_M12_Wald-test.csv \
results/DESeq2/V24_over_M24_Wald-test.csv : results/DESeq2/dds.Rdata \
	src/DESeq2-stats.R
	R -e "setwd('./src/'); source('DEseq2-stats.R')"

## volcano plot
img/volcano-plot.png : results/DESeq2/V0_over_M0_Wald-test.csv \
	results/DESeq2/V12_over_M12_Wald-test.csv \
	results/DESeq2/V24_over_M24_Wald-test.csv \
	src/volcano-plot.R \
	src/ggplot2-themes.R
	R -e "setwd('./src/'); source('volcano-plot.R')"

## principle component plot
img/pca.png : results/DESeq2/complete-dataset_DESeq2-normalized-counts.csv \
	data/Run_2127/Run_2127_wobus.csv \
	src/ggplot2-themes.R
	R -e "setwd('./src/'); source('pca.R')"

## rule to make all images
## add new images here for automatic generation
img : img/volcano-plot.png \
	img/pca.png

