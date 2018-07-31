## R script to generate a plot of GSEA results over time
## David R. Hill
## -----------------------------------------------------------------------------

library(magrittr)
## load in parent-child lists
rpr <- readr::read_delim(file = "../data/ReactomePathwaysRelation.txt",
                         delim = "\t",
                         col_names = FALSE) %>%
    dplyr::rename(parent = X1, child = X2)
rpr <- rpr[grep("R-HSA", rpr$parent),]

## immune system
rpr1 <- subset(rpr, rpr$parent == "R-HSA-168256")

## make a graph from the reactome pathways table
library(igraph)
net <- graph_from_data_frame(d = rpr, directed = TRUE) 

## find all nodes that connect to given REACTOME term
path.net <- induced_subgraph(net,
                            subcomponent(net,
                                         V(net) [ match("R-HSA-168256", V(net)$name) ],
                                         mode = "out"))

## annotate nodes according to mid-tier reactome hierarchy
## search for all nodes connected to each of the 1st tier nodes in innate immunity
index <- lapply(rpr1$child, function(x) {
    out <- subcomponent(net,
                        V(net) [ match(x, V(net)$name) ],
                        mode = "out") %>%
        induced_subgraph(net, .) 
        return(V(out)$name)
}
)
## apply names to lists
names(index) <- rpr1$child
## convert to stacked dataframe
stacked <- stack(index) %>%
    dplyr::rename(ID = values, top.tier = ind)

## join to most distal nodes
categories <- data.frame(ID = V(path.net)$name) %>%
    dplyr::left_join(stacked, by = 'ID')
categories$top.tier <- as.character(categories$top.tier)

## descriptive names for categories
## load reactome pathway list
rp <- readr::read_delim(file = "../data/ReactomePathways.txt",
                         delim = "\t",
                         col_names = FALSE) %>%
    dplyr::rename(top.tier = X1, categories = X2, Species = X3)
rp <- rp[grep("Homo sapiens", rp$Species),]

#path.names <- path.names[which(rp$top.tier %in% categories$top.tier),]

## join to categories index
categories %<>% dplyr::left_join(rp, by = 'top.tier')

## load data file
data <- readr::read_csv(file = "../results/GSEA/combined_GSEA_results.csv")

## subset to Innate immune pathways
data <- data[which(data$ID %in% V(path.net)$name),]

library(magrittr)

## subset here, after clustering is done
data <- subset(data, data$pvalue < 0.05)

## hiararchical clustering
dat <- data %>%
    dplyr::select(ID, NES, comparison) %>%
    tidyr::spread(key = comparison, value = NES)

## treat NA (not significant) as 0 
dat[is.na(dat)] <- 0
    
## determine order for axis clustering
ord <- hclust(dist(dat[,2:4], method = "euclidean"), method = "ward.D")$order
ord2 <- hclust(dist(t(dat[,2:4]), method = "euclidean"), method = "ward.D")$order

## hierarchical clustering of pathways for more coherent plotting
data$Description <- factor(data$Description, levels = unique(data$Description)[ord])
#data$comparison <- factor(data$comparison, levels = unique(data$comparison)[ord2])

## bind category IDs
data <- dplyr::left_join(data, categories, by = 'ID')

## setup the plot
library(ggplot2)
source("ggplot2-themes.R")
plot <- ggplot(data = data[!is.na(data$categories),],
               aes(x = as.factor(comparison),
                   y = Description,
                   fill = NES)) +
    scale_fill_distiller(name = "NES ",
                         palette = "RdYlBu",
                         na.value = "#2C6CAD") +
    geom_tile() +
    facet_grid(categories ~ ., scales = "free_y", space = "free_y") +
    xlab("") + ylab("") +
  #  ggtitle("REACTOME: Immunity") +
    ## theme1 +
    theme(axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 24),
          strip.text.y = element_text(size = 20,
                                      angle = 360),
          strip.background = element_rect(fill = "grey90"),
          legend.position = "bottom",
          legend.text = element_text(size = 24),
          legend.title = element_text(size = 24),
	  panel.spacing = unit(1, "lines"),
          plot.title = element_text(size = 45,
                                    face = "bold",
                                    hjust = 1),
          panel.background = element_rect(fill = "white",
                                          color = "grey80"))


## open device for plotting
png(filename = "../img/gsea-reactome-heatmap_Immunity.png", height = 600, width = 1100)
print(plot)
dev.off()
