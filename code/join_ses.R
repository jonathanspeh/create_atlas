library(SummarizedExperiment)
source("code/functions.R")


files <- list.files(here::here("data", "ses"), full.names = TRUE)
files <- files[grepl("_se.RDS", files)]

counts_list <- lapply(files, readRDS)

combined_se <- combine_se(counts_list)



