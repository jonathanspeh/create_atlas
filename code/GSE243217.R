library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)

# Failed - normalised counts
GEO_accs <- "GSE243217"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE243217&format=file"
dir_path <- here::here("data", GEO_accs)
tar_path <- here::here("data", paste0(GEO_accs, ".tar"))


if(!file.exists(tar_path)){
  dir.create(dir_path)
  curl::curl_download(download,
                      destfile = tar_path)
  untar(tar_path, 
        exdir = dir_path)
}

# continue here, couldn't load sm file
sm <- getGEO(GEO_accs, destdir = dir_path) 

pData(sm[[1]]) |> dplyr::select(starts_with("char")) 
pData(sm[[1]]) |>
  as_tibble() |>
  dplyr::filter(`infection type:ch1` == "Noninfection") 
colnames(pData(sm[[1]]))
glimpse(pData(sm[[1]]))

unique(pData(sm[[1]])$ `characteristics_ch1.1` )




meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = title) |>
  mutate(disease = ifelse(grepl("healthy", `disease state:ch1`), "healthy", `disease state:ch1`),
         age = NA,
         sex = NA,
         dataset = GEO_accs, 
         pediatric = FALSE, 
         source = "whole blood"
         ) |>
  dplyr::select(id, sample_name, pediatric,
                disease, 
                processing_info, source, dataset, age, sex) 





files <- list.files(dir_path, full.names = FALSE)
files <- files[grepl("GSM", files)]

read_counts <- function(name){
  sample <- stringr::str_remove(name, "_.*.txt.gz")
  file <- data.table::fread(here::here(dir_path, name)) 
  colnames(file) <- c("gene", substitute(sample))
  file
}

counts_list <- lapply(files, read_counts)
counts <- purrr::reduce(counts_list, full_join, by = "gene") 


rownames(meta) <- meta$id

all(rownames(meta) == colnames(counts)[-1])

count_mat <- as.matrix(counts[,-1])
counts

meta |>
  count(disease)

se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts_filtered[,-1])),
  colData = meta_filtered,
  metadata = list(notes = "Patients with staph aureus infected foot ulcer, only befoore treatment in assay, after treatment in metadata",
                  full_counts = counts,
                  full_meta = meta))


rownames(se) <- stringr::str_remove(counts_filtered$gene, "\\.[0-9]*") # gene names have some suffices

se <- se[!is.na(rownames(se)),]
# saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))




