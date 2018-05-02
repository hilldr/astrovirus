## R Script to generate principle component analysis of HIE response
## David R. Hill
## -----------------------------------------------------------------------------

library(magrittr)
## load dataset from file
dataset <- "../results/DESeq2/complete-dataset_DESeq2-normalized-counts.csv"
df <- readr::read_csv(file = dataset) %>% dplyr::rename(SYMBOL = X1)

## read in table with sample metadata
## load using the UM core provided sample submission form
samples <- readr::read_csv(file = "../data/Run_2127/Run_2127_wobus.csv", skip = 18)

## extract and format experimental design labels from Description string
samples$group <- gsub('.{2}$', '', samples$Description)
samples$treatment <- substring(samples$Description, 1, 1) %>%
    gsub('V', 'Virus', .) %>%
    gsub('M', 'Mock', .)
samples$hr <- gsub("[^0-9.]", "", samples$Description) %>% as.numeric()

## calculate variance by row (gene)
var <- dplyr::select(df, -SYMBOL) %>% apply(1, sd, na.rm = TRUE)
## adjust cut off according to variance percentile
## effectively, base PCA on only the most highly variable genes
pca.data <- dplyr::select(df, -SYMBOL)[var > quantile(var, 0.9) & var != 0,]
## calculate principle components
pca <- prcomp(t(pca.data), scale = TRUE, center = TRUE)

## put PCA results into dataframe with sample info
scores <- data.frame(Description = colnames(pca.data),
                     pca$x[,1:ncol(pca$x)]) %>%
    dplyr::left_join(samples, by = 'Description')

## PCA plot
library(ggplot2)
library(RColorBrewer)
source("ggplot2-themes.R")

## function to format decimals as precentage
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

plot <- ggplot( data = scores, aes(x = PC1, y = PC2)) +
    geom_point(shape = 21, stroke = 3, aes(fill = as.factor(hr), color = treatment), size = 12) +
    theme1 +
    scale_fill_brewer(palette = "Reds") +
    scale_color_brewer(palette = "Set1", direction = -1) +    
    theme(legend.position = "bottom",
          legend.background = element_rect(colour = "white"),
          legend.key = element_rect(color = "white",
                                    fill = "white")) +
    geom_hline(yintercept = 0,
               size = 1, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0,
               size = 1, linetype = "dashed", color = "grey70") +
    coord_fixed(ratio = 1) +
    xlab(paste("PC1 (",percent(round(summary(pca)$importance[2,1],4)),")",sep = "")) +
    ylab(paste("PC2 (",percent(round(summary(pca)$importance[2,2],4)),")",sep = ""))

## print to plotting device
png(filename = "../img/pca.png", width = 1200, height = 1000)
print(plot)
dev.off()

## alternate PCA
plot <- ggplot( data = scores, aes(x = PC1, y = PC3)) +
    geom_point(shape = 21, stroke = 3, aes(fill = as.factor(hr), color = treatment), size = 12) +
    theme1 +
    scale_fill_brewer(palette = "Reds") +
    scale_color_brewer(palette = "Set1", direction = -1) +    
    theme(legend.position = "bottom",
          legend.background = element_rect(colour = "white"),
          legend.key = element_rect(color = "white",
                                    fill = "white")) +
    geom_hline(yintercept = 0,
               size = 1, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0,
               size = 1, linetype = "dashed", color = "grey70") +
    coord_fixed(ratio = 1) +
    xlab(paste("PC1 (",percent(round(summary(pca)$importance[2,1],4)),")",sep = "")) +
    ylab(paste("PC3 (",percent(round(summary(pca)$importance[2,3],4)),")",sep = ""))

## print to plotting device
png(filename = "../img/pca2.png", width = 1200, height = 1000)
print(plot)
dev.off()

plot <- ggplot( data = scores, aes(x = PC2, y = PC3)) +
    geom_point(shape = 21, stroke = 3, aes(fill = as.factor(hr), color = treatment), size = 12) +
    theme1 +
    scale_fill_brewer(palette = "Reds") +
    scale_color_brewer(palette = "Set1", direction = -1) +    
    theme(legend.position = "bottom",
          legend.background = element_rect(colour = "white"),
          legend.key = element_rect(color = "white",
                                    fill = "white")) +
    geom_hline(yintercept = 0,
               size = 1, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0,
               size = 1, linetype = "dashed", color = "grey70") +
    coord_fixed(ratio = 1) +
    xlab(paste("PC2 (",percent(round(summary(pca)$importance[2,2],4)),")",sep = "")) +
    ylab(paste("PC3 (",percent(round(summary(pca)$importance[2,3],4)),")",sep = ""))

## print to plotting device
png(filename = "../img/pca3.png", width = 1200, height = 1000)
print(plot)
dev.off()
