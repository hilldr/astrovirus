## R script for plotting the results of correlation analysis
## between viral titer and RNA-seq gene counts
## David R. Hill
## -----------------------------------------------------------------------------

library(magrittr)

## load data
data <- readr::read_csv(file = "../results/counts_metadata_titer-correlation.csv")
## calculate position for 'r' on facet plots
#data.position <- data %>% dplyr::group_by(SYMBOL) %>%
#    dplyr::summarise(position = max(count))
#data <- data %>% dplyr::left_join(data.position, by = 'SYMBOL')

## subset to top 10 positively correlated
data.pos <- data[which(data$SYMBOL %in% unique(data$SYMBOL)[1:12]),]
## calculate position for 'r' on facet plots
data.position <- data.pos %>% dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(position = max(count))
data.pos <- data.pos %>% dplyr::left_join(data.position, by = 'SYMBOL')
## top 10 negatively correlated
data.neg <- data[which(data$SYMBOL %in% unique(data$SYMBOL)[(length(unique(data$SYMBOL)) - 11):length(unique(data$SYMBOL))]),]
data.position <- data.neg %>% dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(position = min(count))
data.neg <- data.neg %>% dplyr::left_join(data.position, by = 'SYMBOL')


## Plots -----------------------------------------------------------------------
## Positive correlations
library(ggplot2)
library(scales)
source("ggplot2-themes.R")

plot1 <- ggplot(data = data.pos[data.pos$virus == "VA1",], aes(x = PFU_well, y = count)) +
    geom_smooth(
        colour = "grey70",
        fill = "grey80",
        linetype = "dashed",
        size = 1,
        method = "lm",
        formula = y ~ x,
        level = 0.95
    ) +
    geom_point(shape = 21, size = 5, aes(fill = as.factor(hpi))) +
    facet_wrap(~SYMBOL, scales = "free_y") +
    scale_x_log10(
        labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "b", size = 1,
                        short = unit(.75,"mm"),
                        mid = unit(1,"mm"),
                        long = unit(2,"mm")) +
    xlab("PFU/well") + ylab("TPM") +
    scale_fill_brewer("HPI", palette = "Reds") +
    geom_text(
        size = 5,
        aes(x = 2e3,
            y = position,
            label = paste0("r = ", round(titer_correlation, digits = 3)))) + 
    theme(
        strip.text = element_text(size = 24),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA,
                                    color = "grey70",
                                    size = 1),
        plot.title = element_text(size = 45,
                                  face = "bold",
                                  hjust = 0),
        )

## print to open PNG device    
png(filename = "../img/titer-gene-correlation-plot_POSITIVE.png", width = 1200, height = 900)
print(plot1)
dev.off()


plot2 <- ggplot(data = data.neg[data.neg$virus == "VA1",], aes(x = PFU_well, y = count)) +
    geom_smooth(
        colour = "grey70",
        fill = "grey80",
        linetype = "dashed",
        size = 1,
        method = "lm",
        formula = y ~ x,
        level = 0.95
    ) +
    geom_point(shape = 21, size = 5, aes(fill = as.factor(hpi))) +
    facet_wrap(~SYMBOL, scales = "free_y") +
    scale_x_log10(
        labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "b", size = 1,
                        short = unit(.75,"mm"),
                        mid = unit(1,"mm"),
                        long = unit(2,"mm")) +
    xlab("PFU/well") + ylab("TPM") +
    scale_fill_brewer("HPI", palette = "Reds") +
    geom_text(
        size = 5,
        aes(x = 2e3,
            y = position,
            label = paste0("r = ", round(titer_correlation, digits = 3)))) + 
    theme(
        strip.text = element_text(size = 24),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA,
                                    color = "grey70",
                                    size = 1),
        plot.title = element_text(size = 45,
                                  face = "bold",
                                  hjust = 0),
        )

## print to open PNG device    
png(filename = "../img/titer-gene-correlation-plot_NEGATIVE.png", width = 1200, height = 900)
print(plot2)
dev.off()

