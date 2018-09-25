## plot viral rna-seq data
## David R. Hill 2018-09-24
## -----------------------------------------------------------------------------

## load prerequisites
library(magrittr)

## load qpcr data
## load in viral qpcr data
qpcr <- readr::read_csv(file = "../data/Sample-titers-for-RNAseq.csv", skip = 4)
## load using the UM core provided sample submission form
samples <- readr::read_csv(file = "../data/Run_2127/Run_2127_wobus.csv",
                           skip = 18)
## Create Sample ID column in 'qpcr' to match 'samples'
qpcr$Sample_ID <- paste0("Sample_",qpcr$code_seq)

## load RNA-seq data
seq <- readr::read_csv(file = "../results/VA1_alignment_stats.csv", col_names = TRUE)

## join datasets
data <- dplyr::left_join(seq, qpcr, by = 'Sample_ID')

## calculate viral TPM
data$tpm <- (data$n_pseudoaligned/data$n_processed)*1e6

## correlate TPM and genome copies (qpcr)
fit <- lm(data$tpm ~ data$PFU_well)

## generate fit plot
library(ggplot2)
library(scales)
plot <- ggplot(data = data,
               aes(x = PFU_well,
                   y = tpm)) +
    geom_smooth(method='lm',formula=y~x) +
    geom_point(shape = 21,
               size = 7,
               stroke = 1,
               aes(fill = virus)) +
    scale_y_log10(
        breaks = trans_breaks("log10", function(x) 10^x)(c(0.001, 1e4)),
        labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(
        breaks = trans_breaks("log10", function(x) 10^x)(c(1, 1e7)),
                  labels = trans_format("log10", math_format(10^.x))) +
    ylab(latex2exp::TeX("\\frac{viral transcripts}{1\\times{}10$^{6}$ RNA-seq reads}")) +
    xlab(latex2exp::TeX("genome copies \\times{}HIE$^{-1}$ (RT-qPCR)")) +
    scale_fill_discrete(name = "Treatment")+
    scale_fill_brewer(palette = "Set1", direction = -1) + 
    annotate(geom = "text",
             x = 2.5, y = 2e3,
             label = paste0("P = ", 2.5e-14),
             size = 6) +
    annotate(geom = "text",
             x = 2.5, y = 1e3,
             label = latex2exp::TeX("\\textit{r}$^{2}$ = 0.98"),
             size = 6) +
    annotation_logticks(sides = "lb", size = 0.5,
                        short = unit(0.5,"mm"),
                        mid = unit(1,"mm"),
                        long = unit(2,"mm"))+  
    theme(axis.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.text=element_text(size = 18),
          legend.position = "right",
          panel.background = element_blank(),
          panel.border = element_rect(color = "grey30",
                                      fill = NA),
          axis.title = element_text(size = 24),
          strip.text = element_text(size = 24))

#png(filename = "../img/seq_pcr_correlation.png",
 #   width = 800, height = 600)
print(plot)
#dev.off()

png(filename = "../img/seq_pcr_correlation.png",
    width = 800, height = 600)
print(plot)
dev.off()
