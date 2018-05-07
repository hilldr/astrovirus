## R Script to generate principle component analysis of HIE response
## David R. Hill
##  Modified by Ryan Berger:
##   - Use the dds object instead of counts.csv to make PCA plot
##   - rlog transform data prior to making plot
##
## -----------------------------------------------------------------------------


# Load packages
library(magrittr)
library(ggplot2)
library(RColorBrewer)
source("ggplot2-themes.R")

## Load dds object
load('../results/DESeq2/dds.Rdata')

## Perform variance stabilizing transformations on dds using rlog
## use argument blind = FALSE when multiple replicates are present
rld <- DESeq2::rlog(dds, blind = FALSE)

## Rename V and M to Virus and Mock
colData(rld)
colData(rld)$treatment <- gsub('V', 'Virus', colData(rld)$treatment) %>% 
  gsub('M', 'Mock', .)


# PCA plot (use just only for PC variance estimates)
pca <- plotPCA(rld, intgroup = c('hr', 'treatment'))

# Get PCA data
pca.df <- plotPCA(rld, intgroup = c('group','hr','treatment'), returnData = TRUE)


plot <- ggplot(data = pca.df, aes(x = PC1, y = PC2))+ 
geom_point(shape = 21, stroke = 3, aes(fill = as.factor(hr), color = treatment), size = 12) +
  theme1 + 
  scale_fill_brewer(palette = "Reds") +
  scale_color_brewer(palette = "Set1", direction = -1) +    
  theme(legend.position = "right") +
  geom_hline(yintercept = 0,
             size = 1, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0,
             size = 1, linetype = "dashed", color = "grey70") +
  coord_fixed(ratio = 1) +
  xlab(pca$labels$x) + #pull variance estimates from al. plotPCA call
  ylab(pca$labels$y)

## open plot device
png(filename = "../img/pca.png",
    width = 1200, height = 1000)
print(plot)
dev.off()


