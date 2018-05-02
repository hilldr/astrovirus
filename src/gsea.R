## R script for Gene Set Enrichment Analysis
## David R. Hill
## -----------------------------------------------------------------------------

## load prerequisites
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

## create functions for GO & REACTOME GSEA
gsea.go <- function(input = input,
                    file = file,
                    data.dir = data.dir,
                    comparison = comparison){
    
   input <- readr::read_csv(file = input) %>%
        dplyr::rename(SYMBOL = X1)
    input <- bitr(input[input$baseMean >= 1,]$SYMBOL,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db") %>%
        dplyr::left_join(.,input[input$baseMean >= 1,], by = 'SYMBOL')
    
    input <- input[order(-input$log2FoldChange),]
    up.list<- sort(input$log2FoldChange, decreasing = TRUE)
    names(up.list) <- input$ENTREZID
    ## GO GSEA
    go.gsea <- gseGO(geneList     = up.list,
                     OrgDb        = org.Hs.eg.db,
                     ont          = "BP",
                     nPerm        = 1000,
                     minGSSize    = 50,
                     maxGSSize    = 500,
                     pvalueCutoff = 1,
                     verbose      = TRUE)
    save(go.gsea, file = file.path(data.dir,paste0(file,".GO.RData")))
    go.gsea@result$database <- "GO"
    go.gsea@result$comparison <- comparison
    write.csv(go.gsea@result,
              file = file.path(data.dir,paste0(file,"_GO.csv")), row.names = FALSE)
    return(go.gsea@result)
    
}

gsea.reactome <- function(input = input,
                          file = file,
                          data.dir = data.dir,
                          comparison = comparison){

    input <- readr::read_csv(file = input) %>%
        dplyr::rename(SYMBOL = X1)
    input <- bitr(input[input$baseMean >= 1,]$SYMBOL,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db") %>%
        dplyr::left_join(.,input[input$baseMean >= 1,], by = 'SYMBOL')
    
    input <- input[order(-input$log2FoldChange),]
    up.list<- sort(input$log2FoldChange, decreasing = TRUE)
    names(up.list) <- input$ENTREZID
    ## REACTOME GSEA
    reactome.gsea <- gsePathway(geneList  = up.list,
                 nPerm        = 1000,
                 minGSSize    = 10,
                 pvalueCutoff = 1,
                 verbose      = TRUE)
    save(reactome.gsea, file = file.path(data.dir,paste0(file,".REACTOME.RData")))
    reactome.gsea@result$database <- "REACTOME"
    reactome.gsea@result$comparison <- comparison
    write.csv(reactome.gsea@result,
              file = file.path(data.dir,paste0(file,"_REACTOME.csv")), row.names = FALSE)
    return(reactome.gsea@result)
    
}

## use mapply to iterate over all files in directory
data <- mapply(gsea.reactome,
               input = file.path("../results/DESeq2",
                                 list.files(path = "../results/DESeq2/")[grep("Wald", list.files(path = "../results/DESeq2/"))]),
               data.dir = "../results/GSEA",
               file = paste0("results", substr(list.files(path = "../results/DESeq2/")[grep("Wald",list.files(path = "../results/DESeq2/"))], 2, 3)),
               comparison = gsub("_", "", substr(list.files(path = "../results/DESeq2/")[grep("Wald",list.files(path = "../results/DESeq2/"))], 2, 3)),
	       SIMPLIFY = FALSE)

## combine list output into a single data.frame
df.reactome <- do.call("rbind", data)

## use mapply to iterate over all files in directory
data <- mapply(gsea.go,
               input = file.path("../results/DESeq2",
                                 list.files(path = "../results/DESeq2/")[grep("Wald", list.files(path = "../results/DESeq2/"))]),
               data.dir = "../results/GSEA",
               file = paste0("results", substr(list.files(path = "../results/DESeq2/")[grep("Wald",list.files(path = "../results/DESeq2/"))], 2, 3)),
               comparison = gsub("_", "", substr(list.files(path = "../results/DESeq2/")[grep("Wald",list.files(path = "../results/DESeq2/"))], 2, 3)),
	       SIMPLIFY = FALSE)

## combine list into a single data.frame
df.go <- do.call("rbind", data)

## combine REACTOME and GO results into a single dataframe
write.csv(rbind(df.reactome, df.go), file = "../results/GSEA/combined_GSEA_results.csv")

## clear working memory
rm(list = ls())
