library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)


# Different covid 19 variants


#curl::curl_download("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE190680&format=file",
#                    destfile = here::here("data", "GSE190680.tar")
#                    )

#untar(here::here("data", "GSE190680.tar"), 
#      exdir = here::here("data", "GSE190680"))

sm <- getGEO("GSE190680", destdir = here::here("data", "GSE190680"))

meta <- pData(sm$GSE190680_series_matrix.txt.gz) |> 
  dplyr::select(id = geo_accession,
         sample_name = title,
         sample_type =  characteristics_ch1.2,
         age = "age:ch1",
         sex = "gender:ch1"
         ) |> 
  mutate(variant = sub("^(.*?)_S.*", "\\1", sample_name),
         individual = sub("^.*_(S\\d+).*", "\\1", sample_name),
         sampling_point = sub("^.*?(\\d+)\\D*$", "\\1", sample_name),
         disease = "COVID-19",
         dataset = "GSE1906080"
         )
         
dplyr::filter(meta, 
              variant == "Alpha") |>
  group_by(individual) |>
  count() |>
  arrange(desc(n))


files <- list.files("data/GSE190680", full.names = FALSE)
files <- files[grepl("GSM", files)]

read_counts <- function(name){
  sample <- stringr::str_remove(name, "_.*.txt.gz")
  file <- fread(paste0("data/GSE190680/", name)) |>
    dplyr::filter(!grepl("^__", V1)) 
  colnames(file) <- c("gene", substitute(sample))
  file
}

counts <- lapply(files, read_counts)

counts_raw <- purrr::reduce(counts, full_join, by = "gene") 
meta_point1 <- dplyr::filter(meta, sampling_point == 1)

counts_point1 <-  dplyr::select(counts_raw, c(gene, meta_point1$id))


se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts_point1[,-1])),
  colData = meta_point1,
  metadata = list(notes = "assay counts contains counts only for sampling point 1, full data are saved in metadata",
                  full_counts = counts_raw,
                  full_meta = meta
                  )
)


rownames(se) <- mapIds(EnsDb.Hsapiens.v86,
                       keys = counts_raw$gene,
                       keytype = "SYMBOL",
                       columns = "ENSEMBL")

se <- se[!is.na(rownames(se)),]


saveRDS(se, "data/ses/GSE190680_se.RDS")



