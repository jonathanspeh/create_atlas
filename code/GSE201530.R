library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)


GEO_accs <- "GSE201530"
file <- paste0(GEO_accs, ".tar")
dir_path <- paste0(here::here("data", GEO_accs), "/")

# curl::curl_download(paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=", 
#                            GEO_accs, 
#                            "&format=file"),
#                     destfile = here::here("data", file))
# 
# untar(here::here("data", file), 
#       exdir = dir_path)

sm <- getGEO(GEO_accs, destdir = dir_path)


meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::select(id = geo_accession,
                sample_name = title,
                sample_type =  characteristics_ch1.2,
                age = "age:ch1",
                sex = "gender:ch1",
                group = "disease state:ch1",
                source = "cell type:ch1", 
                processing_info) |> 
  mutate(individual = sub("^.*_(\\d+)_.*$", "\\1", sample_name),
         disease = case_when(grepl("control", sample_type) ~ "healthy",
                             TRUE ~ "COVID-19"),
         variant = "Omicron",
         dataset = GEO_accs,
         pediatric = FALSE)



files <- list.files(dir_path, full.names = FALSE)
files <- files[grepl("GSM", files)]

read_counts <- function(name){
  sample <- stringr::str_remove(name, "_.*.txt.gz")
  file <- fread(paste0(dir_path, name)) |>
    dplyr::filter(!grepl("^__", V1)) 
  colnames(file) <- c("gene", substitute(sample))
  file
}

counts <- lapply(files, read_counts)

counts_raw <- purrr::reduce(counts, full_join, by = "gene") 



se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts_raw[,-1])),
  colData = meta,
  metadata = list(notes = "covid 19 (Omicron) patienst with ans without previous challenge (infection / vaccine) and healthy controls")
)

rownames(se) <- mapIds(EnsDb.Hsapiens.v86,
                       keys = counts_raw$gene,
                       keytype = "SYMBOL",
                       columns = "ENSEMBL")

se <- se[!is.na(rownames(se)),]
saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))



