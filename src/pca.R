## R Script to generate principle component analysis of HIE response
## David R. Hill
##  Modified by Ryan Berger:
##   - Use the dds object instead of counts.csv to make PCA plot
##   - rlog transform data prior to making plot
##
## To Do: fix legend so that sample IDs aren't in it
## -----------------------------------------------------------------------------

# Working directory should be ../astrovirus/src
getwd()

# Load packages
library(magrittr)
library(ggplot2)
library(RColorBrewer)
source("ggplot2-themes.R")

# Load dds object
load('../results/DESeq2/dds.Rdata')

# Perform variance stabilizing transformations on dds using rlog
# use argument blind = FALSE when multiple replicates are present
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

# Rename V and M to Virus and Mock
colData(rld)
colData(rld)$treatment <- gsub('V', 'Virus', colData(rld)$treatment) %>% 
  gsub('M', 'Mock', .)


# PCA plot
pca <- plotPCA(rld, intgroup = c('hr', 'treatment'))
pca +
  ggtitle('PCA Plot of Enteroids') +
  geom_point(shape = 21, stroke = 3, aes(fill = hr, color = treatment, size = 12)) + 
  theme1 +
  scale_fill_brewer(palette = "Reds") +
  scale_color_brewer(palette = "Set1", direction = -1)  +  
  theme(legend.position = "bottom",
      legend.background = element_rect(colour = "white"),
      legend.key = element_rect(color = "white",
                                fill = "white")) +
  geom_hline(yintercept = 0,
             size = 1, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0,
             size = 1, linetype = "dashed", color = "grey70") +
  coord_fixed(ratio = 1)  + 
  guides(hr = 'legend', treatment = 'legend', group = 'none', size = 'none') # This works for size but not group

# To do: clean up legend so that the group IDs (e.g. 12:Mock) aren't in there


## Other way to do it
# plotPCA returnData = TRUE gives the plotting values as a dataframe
# Can make ggplot but doesn't give the percentage of variance values

# Get PCA data
pca.df <- plotPCA(rld, intgroup = c('group','hr','treatment'), returnData = TRUE)
pca.df


plot <- ggplot( data = pca.df, aes(x = PC1, y = PC2))+ 
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
  coord_fixed(ratio = 1) 

plot  
## Doesn't have the variance percentages to put on the plot!