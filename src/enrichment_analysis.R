## pathway enrichment analysis for top up- and down-regulated genes
## David R. Hill 2018-08-01
## -----------------------------------------------------------------------------

## load prerequisites
library(magrittr)

## load expression data
data <- readr::read_csv(file = "../results/DESeq2/V24_over_M24_Wald-test.csv",
                        col_names = TRUE)

## top 100 up-regulated genes
top_up <- subset(data,data$pvalue < 0.05) %>%
    .[order(-.$log2FoldChange),] %>%
    .$symbol %>%
    .[1:100]

## top 100 down-regulated genes
top_down <- subset(data,data$pvalue < 0.05) %>%
    .[order(-.$log2FoldChange),] %>%
    .$symbol %>%
    .[1:100]

## pathway analysis of clusters
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

## function to generate reactome enrichment scores
reactome.oea <- function(genes){ 
    ## create list of gene names from list
    ids <- bitr(genes,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db")

    ## scan REACTOME
    reactome <- enrichPathway(gene = ids$ENTREZID,
                              qvalueCutoff = 0.05,
                              readable = TRUE)
    out <- as.data.frame(reactome@result)
    ## write out results
    return(out)
}

## execute function and write out results
dir.create("../results/enrichment_analysis")
top_up_paths <- reactome.oea(top_up) %>%
    write.csv(file = "../results/enrichment_analysis/top_100_UP_REACTOME.csv")
top_down_paths <- reactome.oea(top_down) %>%
    write.csv(file = "../results/enrichment_analysis/top_100_DOWN_REACTOME.csv")
