## summarize viral alignment stats from kallisto
## David R. Hill 2018-09-24
## -----------------------------------------------------------------------------

## load prerequisites
library(magrittr)

## function to extract run info data from a directory
info_df <- function(dir) {
    info.files <- list.files(path = dir,
                         pattern = "run_info.json",
                         full.names = TRUE,
                         recursive = TRUE)
    df <- mapply(
        FUN = function(x) {
            ## load library for json files
            data <- jsonlite::fromJSON(x, flatten = TRUE) %>%
                as.data.frame()
            data$Sample_ID <- strsplit(x, split = "/")[[1]][6]
            data$Sample_ID <- strsplit(data$Sample_ID, split = "_")[[1]][1] %>%
                paste0("Sample_",.)
            data$idx <- strsplit(as.vector(data$call), split = " ")[[1]][4] %>%
                gsub("../data/genomes/", "", x = .)
            return(data)
        },
        x = info.files,
        SIMPLIFY = FALSE)
    return(df)
}

## extract data from run info files, concatenate into dataframe, & write out
VA1.info <- info_df(dir = "../results/astrovirus_va1/kallisto/") %>%
    do.call("rbind", .) %>%
    dplyr::left_join(
               readr::read_csv(file = "../data/Run_2127/Run_2127_wobus.csv",
                               skip = 18,
                               col_names = TRUE),
               by = 'Sample_ID') %>%
    dplyr::left_join(
               readr::read_csv(file = "../data/genome_index.csv",
                               col_names = TRUE),
               by = 'idx')

## write out results    
write.csv(VA1.info, file = "../results/VA1_alignment_stats.csv")

## extract data from run info files, concatenate into dataframe, & write out
Hs.info <- info_df(dir = "../results/Run_2127/") %>%
    do.call("rbind", .) %>%
    dplyr::left_join(
               readr::read_csv(file = "../data/Run_2127/Run_2127_wobus.csv",
                               skip = 18,
                               col_names = TRUE),
               by = 'Sample_ID') %>%
    dplyr::left_join(
               readr::read_csv(file = "../data/genome_index.csv",
                               col_names = TRUE),
               by = 'idx')

## write out results
write.csv(Hs.info, file = "../results/Hs_alignment_stats.csv")

## combine all alignment stats
aln.stats <- rbind(VA1.info,
                   Hs.info)
write.csv(aln.stats, file = "../results/ALL_alignment_stats.csv")
