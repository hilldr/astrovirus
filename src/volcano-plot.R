## R script to generate volcano plot figure by time point
## David R. Hill
## -----------------------------------------------------------------------------
library(magrittr)
## data import
data.dir <- "../results/DESeq2/"

## read in individual time points
hr.0 <- readr::read_csv(file = file.path(data.dir,
                                         "V0_over_M0_Wald-test.csv"))
hr.0$hr <- 0
hr.12 <- readr::read_csv(file =  file.path(data.dir,
                                           "V12_over_M12_Wald-test.csv"))
hr.12$hr <- 12
hr.24 <- readr::read_csv(file =  file.path(data.dir,
                                           "V24_over_M24_Wald-test.csv"))
hr.24$hr <- 24

## combine into single dataframe
data <- rbind(hr.0, hr.12, hr.24) %>%
    dplyr::rename(SYMBOL = X1)

## create status catergory for assigning colors
data$status <- ifelse(data$pvalue > 0.05 | is.na(data$padj), "a",
                    ifelse(data$log2FoldChange > 0, "b", "c"))
## sort by status
data <- data[order(data$status),]


## pre-reqs
library(ggplot2)
source("ggplot2-themes.R")

## set up facet plot
plot <- ggplot(data = data, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(shape = 21, size = 2, aes(fill = status, color = status)) +
    facet_grid(. ~ hr) +
    ## set axis limits and labels
    xlim(c(-2,2)) + ylim(c(0,10)) +
    xlab(expression(paste("log"[2],"FC(Virus/Mock)"))) +
    ylab(expression(paste("-log"[10],"(P-value)"))) +
    ## colors and theme
    scale_fill_manual(values = c("grey70", color.set[1], color.set[2])) +
    scale_color_manual(values = c("grey70", color.set[1], color.set[2])) +
    theme(
        strip.text = element_text(size = 48),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 24),
	legend.position = "none",
	panel.background = element_blank(),
        panel.border = element_rect(fill = NA,
                                    color = "grey85",
                                    size = 1)				    
    )

## open plot device
png(filename = "../img/volcano-plot.png",
    width = 1000, height = 1000)
print(plot)
dev.off()
                               
