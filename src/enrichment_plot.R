## Plot results of enrichment analysis
## David R. Hill 2018-08-01
## -----------------------------------------------------------------------------

## load prerequisites
library(magrittr)
library(ggplot2)
library(ggstance)
        
## load data from file
data_up <- readr::read_csv(file = "../results/enrichment_analysis/top_100_UP_REACTOME.csv",
                           col_names = TRUE)[1:10,]
data_down <- readr::read_csv(file = "../results/enrichment_analysis/top_100_DOWN_REACTOME.csv",
                             col_names = TRUE)[1:10,]
data_up$status <- "up-regulated"
data_down$status <- "down-regulated"
data <- rbind(data_up, data_down)

## sort by pvalue
data <- data[order(-data$pvalue),]
data$Description <- factor(data$Description, levels = unique(data$Description))

plot_oea <- ggplot(data =  data,
                  aes(x = -log10(pvalue),
                      y = Description)) +
    geom_barh(stat = "identity") +
    facet_grid(status~., scales = "free_y") +
    ggtitle("Top 100 up- and down-regulated genes") +
    ylab("REACTOME pathway") +
    annotate(geom = "segment",
             x = -log10(0.05),
             xend = -log10(0.05),
             y = 0,
             yend = 11,
             linetype = "dashed",
             color = "red") +
    theme(axis.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.position = "bottom",
          panel.background = element_blank(),
          panel.border = element_rect(color = "grey30",
                                      fill = NA),
          axis.title = element_text(size = 24),
          strip.text = element_text(size = 24),
          strip.text.y = element_text(size = 22))

png(filename = "../img/enrichment_plot.png",
    width = 1000, height = 600)
print(plot_oea)
dev.off()
