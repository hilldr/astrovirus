## R script to correlate astrovirus titer with gene expression
## David R. Hill
## -----------------------------------------------------------------------------

library(magrittr)

## load data and reformat
## load in viral titer data
titer <- readr::read_csv(file = "../data/Sample-titers-for-RNAseq.csv", skip = 4)
## load using the UM core provided sample submission form
samples <- readr::read_csv(file = "../data/Run_2127/Run_2127_wobus.csv",
                           skip = 18)
## Create Sample ID column in 'titer' to match 'samples'
titer$Sample_ID <- paste0("Sample_",titer$code_seq)

## Load RNA-seq dataset from file
df <- readr::read_csv(file = "../results/DESeq2/complete-dataset_DESeq2-normalized-counts.csv") %>%
    dplyr::rename(SYMBOL = X1) %>%
    tidyr::gather(Description, count, -SYMBOL)

## scaling function
scale_this <- function(x) {as.vector(scale(x))}

## combine datasets
data <- samples %>%
    dplyr::left_join(titer, by = 'Sample_ID') %>%
    dplyr::left_join(df, by = 'Description') %>%
    dplyr::group_by(SYMBOL) %>%
    ## corelate (unscaled)
    dplyr::mutate(cor = cor(count, PFU_well))

## sort based on correlation
data <- data[order(-data$cor),]

## write out correlation data to file
dplyr::rename(data, titer_correlation = cor) %>%
    write.csv(file = "../results/counts_metadata_titer-correlation.csv")
