library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)


# TODO - figure out what to do with techniocal? replicates - average them? 
curl::curl_download("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE190680&format=file",
                    destfile = here::here("data", "GSE190680.tar")
                    )

untar(here::here("data", "GSE190680.tar"), 
      exdir = here::here("data", "GSE190680")
      )

sm <- getGEO("GSE190680", destdir = here::here("data", "GSE190680"))

meta <- pData(sm$GSE190680_series_matrix.txt.gz) |> 
  select(accession = geo_accession,
         sample_name = title,
         sample_type =  characteristics_ch1.2,
         age = "age:ch1",
         sex = "gender:ch1"
         ) |> 
  mutate(variant = sub("^(.*?)_S.*", "\\1", sample_name),
         individual = sub("^.*_(S\\d+).*", "\\1", sample_name),
         tech_replicate = sub("^.*?(\\d+)\\D*$", "\\1", sample_name),
         disease = "COVID-19",
         dataset = "GSE1906080"
         )
         
files <- list.files("data/GSE190680", full.names = FALSE)

read_counts <- function(name){
  sample <- stringr::str_remove(name, "_.*.txt.gz")
  file <- fread(paste0("data/GSE190680/", name)) |>
    filter(!grepl("^__", V1)) 
  colnames(file) <- c("gene", substitute(sample))
  file
}


counts <- lapply(files, read_counts)

counts_raw <- purrr::reduce(counts, full_join, by = "gene") 


se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts_raw[,-1])),
  colData = meta
)

rownames(se) <- counts_raw$gene
colnames(se)


prcomp(t(assay(se, "counts")))$x |>
  as_tibble(rownames = "accession") |>
  left_join(meta) |>
  ggplot(aes(x = PC1, y = PC2, colour = sample_type)) +
  geom_point()


sm
