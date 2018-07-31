## Plot heatmap of TPMs for genes correlated with viral titer
## David R. Hill 2018-07-31
## -----------------------------------------------------------------------------

## load prerequisites
library(magrittr)

## load data
data <- readr::read_csv(file = "../results/counts_metadata_titer-correlation.csv")

## subset to top 10 positively correlated
top.corr <- c(unique(data$SYMBOL)[1:15],
              unique(data$SYMBOL)[(length(unique(data$SYMBOL)) - 14):length(unique(data$SYMBOL))])
data.top <- data[which(data$SYMBOL %in% top.corr),]

## mean by condition
df.tidy.mean <- data.top %>%
    dplyr::group_by(hpi, virus, SYMBOL) %>%
    dplyr::summarise(mean = mean(count),
                     var_sd = sd(count),
                     num = n())

scale_this <- function(x) as.vector(scale(x))

data.scaled <- df.tidy.mean %>% dplyr::group_by(SYMBOL) %>%
    dplyr::mutate(zscore = scale_this(mean))

## spread matrix for hierarchical clustering
data.scaled$samples <- paste0(data.scaled$virus, "_", data.scaled$hpi)
df.scaled.sprd <- data.scaled %>%
    dplyr::select(-var_sd, -num, -mean, -hpi, -virus) %>%
    tidyr::spread(samples, zscore)

## calculate distance between genes
ord <- hclust(dist(df.scaled.sprd[,2:7],
                   method = "euclidean"),
              method = "ward.D")$order
## impose order on genes
data.scaled$SYMBOL <- factor(data.scaled$SYMBOL,
                             levels = unique(data.scaled$SYMBOL)[ord])


## setup plot ------------------------------------------------------------------
library(ggplot2)

plot4 <- ggplot(data = data.scaled,
                aes(x = as.factor(hpi), y = SYMBOL)) +
    geom_tile(stat = "identity",
              color = "grey70",
              aes(fill = zscore)) + 
    xlab("Time post-infection (h)") + ylab("") +
    facet_grid(.~virus) +
    scale_fill_distiller(name = "Z-score ", palette = "RdYlBu") +
    theme(
        strip.text = element_text(size = 24),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA,
                                    color = NA,
                                    size = 1),
        plot.title = element_text(size = 45,
                                  face = "bold",
                                  hjust = 0),
        )

png(filename = "../img/heatmap_titer_correlation.png",
    width = 500, height = 1200, res = 1/1200)
print(plot4)
dev.off()
