## R script for GSEA of titer correlation data
## David R. Hill
## -----------------------------------------------------------------------------

library(magrittr)

## load data and reformat
data <- readr::read.csv("../results/counts_metadata_titer-correlation.csv")

## GSEA based on correlation results
## create function for GO & REACTOME GSEA
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

## convert gene symbols to entrez IDs
input <- bitr(data$SYMBOL,
              fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = "org.Hs.eg.db") %>%
    dplyr::left_join(., data, by = 'SYMBOL')

## sort by correlation
input <- input[order(-input$cor),]
up.list<- sort(input$cor, decreasing = TRUE)
names(up.list) <- input$ENTREZID

## GO gsea
go.gsea <- gseGO(geneList     = up.list,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 nPerm        = 1000,
                 minGSSize    = 50,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = TRUE)
save(go.gsea, file = "../results/GSEA/GSEA_titer-corr_GO.GO.RData")
go.gsea@result$database <- "GO"
write.csv(go.gsea@result, file = "../results/GSEA/GSEA_titer-corr_GO.csv", row.names = FALSE)
	  


## REACTOME GSEA
reactome.gsea <- gsePathway(geneList     = up.list,
                            nPerm        = 1000,
                            minGSSize    = 10,
                            pvalueCutoff = 0.05,
                            verbose      = TRUE)
save(reactome.gsea, file = "../results/GSEA/GSEA_titer-corr_REACTOME.REACTOME.RData")
reactome.gsea@result$database <- "REACTOME"
write.csv(reactome.gsea@result,
          file = "../results/GSEA/GSEA_titer-corr_REACTOME.csv", row.names = FALSE)
	  





