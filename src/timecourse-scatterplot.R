## R script to generate scatterplot of all gene expression at each time point
## David R. Hill
## -----------------------------------------------------------------------------

library(magrittr)

## set data directory
data.dir <- "../results/DESeq2/"
## load data
v0 <- readr::read_csv(file = file.path(data.dir, "V0_over_M0_Wald-test.csv"))
v0$hr <- 0
v12 <- readr::read_csv(file = file.path(data.dir, "V12_over_M12_Wald-test.csv"))
v12$hr <- 12
v24 <- readr::read_csv(file = file.path(data.dir, "V24_over_M24_Wald-test.csv"))
v24$hr <- 24
## one dataframe to bind them all
data <- rbind(v0, v12, v24)
## make a new status column that will
data$status <- ifelse(data$pvalue > 0.05 | is.na(data$pvalue), "a",
                    ifelse(data$log2FoldChange > 0, "b", "c"))
## set order so that blue and red are plotted on top of grey
data <- data[order(data$status),]

## plot
library(ggplot2)
source("ggplot2-themes.R")
multi.volcano <- ggplot(data = data, aes(y = log2FoldChange, x = factor(hr))) +
    geom_point(position = position_jitter(w = 0.33),
               aes(fill = status, color = status),
               shape = 21, size = 2) +
    scale_fill_manual(values = c("grey70", color.set[1], color.set[2])) +
    scale_color_manual(values = c("grey70", color.set[1], color.set[2])) +
    ylim(c(-5, 5)) +
    xlab("Time post-microinjection (h)") +
    ylab(expression(paste("log"[2],"FC over PBS"))) +
    theme1 + 
    theme(axis.text.x = element_text(size = 24))

## print to plotting device
png(filename = "../img/timecourse-scatterplot.png", width = 1000, height = 800)
print(multi.volcano)
dev.off()


