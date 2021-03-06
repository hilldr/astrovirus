## normalized counts
## update this to also require kallisto output
results/DESeq2/complete-dataset_DESeq2-normalized-counts.csv : src/DESeq-export-counts.R
	R -e "setwd('./src/'); source('DESeq2-export-counts.R')"

## Wald test results
results/DESeq2/V0_over_M0_Wald-test.csv \
results/DESeq2/V12_over_M12_Wald-test.csv \
results/DESeq2/V24_over_M24_Wald-test.csv : results/DESeq2/dds.Rdata \
	src/DESeq2-stats.R
	R -e "setwd('./src/'); source('DESeq2-stats.R')"

## volcano plot
img/volcano-plot.png : results/DESeq2/V0_over_M0_Wald-test.csv \
	results/DESeq2/V12_over_M12_Wald-test.csv \
	results/DESeq2/V24_over_M24_Wald-test.csv \
	src/volcano-plot.R \
	src/ggplot2-themes.R
	R -e "setwd('./src/'); source('volcano-plot.R')"

## Summary table
results/DESeq2/V0_over_M0_Gene-Summary.png \
results/DESeq2/V12_over_M12_Gene-Summary.png \
results/DESeq2/V24_over_M24_Gene-Summary.png : results/DESeq2/V0_over_M0_Wald-test.csv \
	results/DESeq2/V12_over_M12_Wald-test.csv \
	results/DESeq2/V24_over_M24_Wald-test.csv
	R -e "setwd('./src/'); source('summary-table.R')"

## timecourse scatterplot
img/timecourse-scatterplot.png : results/DESeq2/V0_over_M0_Wald-test.csv \
	results/DESeq2/V12_over_M12_Wald-test.csv \
	results/DESeq2/V24_over_M24_Wald-test.csv \
	src/timecourse-scatterplot.R \
	src/ggplot2-themes.R
	R -e "setwd('./src/'); source('timecourse-scatterplot.R')"

## principle component plot
img/pca.png : results/DESeq2/complete-dataset_DESeq2-normalized-counts.csv \
	data/Run_2127/Run_2127_wobus.csv \
	src/pca.R \
	src/ggplot2-themes.R
	R -e "setwd('./src/'); source('pca.R')"

## GSEA analysis
results/GSEA/combined_GSEA_results.csv : results/DESeq2/V0_over_M0_Wald-test.csv \
	results/DESeq2/V12_over_M12_Wald-test.csv \
	results/DESeq2/V24_over_M24_Wald-test.csv \
	src/gsea.R
	R -e "setwd('./src/'); source('gsea.R')"

## GSEA heatmap plot
img/gsea-reactome-heatmap.png : results/GSEA/combined_GSEA_results.csv \
	src/gsea-reactome-heatmap.R
	R -e "setwd('./src/'); source('gsea-reactome-heatmap.R')"

## titer correlation analysis
## output correlation scores
results/counts_metadata_titer-correlation.csv : data/Sample-titers-for-RNAseq.csv \
	data/Run_2127/Run_2127_wobus.csv \
	results/DESeq2/complete-dataset_DESeq2-normalized-counts.csv \
	src/titer-correlation.R
	R -e "setwd('./src/'); source('titer-correlation.R')"

## single gene titer plots
img/titer-gene-correlation-plot_NEGATIVE.png \
img/titer-gene-correlation-plot_POSITIVE.png \
img/titer-gene-correlation-plot_gene_list.png : src/titer-plot.R \
	results/counts_metadata_titer-correlation.csv \
	results/gene_list.csv
	R -e "setwd('./src/'); source('titer-plot.R')"

## GSEA of titer results
results/GSEA/GSEA_titer-corr_GO.csv \
results/GSEA/GSEA_titer-corr_REACTOME.csv : results/counts_metadata_titer-correlation.csv \
	src/titer-GSEA.R
	R -e "setwd('./src/'); source('titer-GSEA.R')"

## rule to make all images
## add new images here for automatic generation
img : img/volcano-plot.png \
	img/pca.png \
	img/timecourse-scatterplot.png \
	results/DESeq2/V0_over_M0_Gene-Summary.png \
	results/DESeq2/V12_over_M12_Gene-Summary.png \
	results/DESeq2/V24_over_M24_Gene-Summary.png \
	img/gsea-reactome-heatmap.png

