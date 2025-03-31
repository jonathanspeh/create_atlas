library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)


GEO_accs <- "GSE221234"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE221234&format=file"
dir_path <- here::here("data", GEO_accs)
tar_path <- here::here("data", paste0(GEO_accs, ".tar"))


if(!file.exists(tar_path)){
  dir.create(dir_path)
  curl::curl_download(download,
                      destfile = tar_path)
  untar(tar_path, 
        exdir = dir_path)
}


sm <- getGEO(GEO_accs, destdir = dir_path) 

meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = title,
                source = `source_name_ch1`,
                disease = `disease state:ch1`,
                timepoint = `sampling group:ch1`
                ) |>
  mutate(group = paste(disease, timepoint, sep = "_"),
         age = NA,
         sex = NA,
         dataset = GEO_accs, 
         pediatric = FALSE,
         ) |>
  dplyr::select(id, sample_name, pediatric,
                disease, 
                processing_info, source, dataset, group, age, sex, timepoint) 

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

meta_filtered <- meta |>
  dplyr::filter(timepoint == "early")

meta_filtered$timepoint


counts_filtered <- counts |>
  dplyr::select(gene, meta_filtered$id)

rownames(meta_filtered) <- meta_filtered$id

all(rownames(meta_filtered) == colnames(counts_filtered)[-1])

se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts_filtered[,-1])),
  colData = meta_filtered,
  metadata = list(notes = "Patients with covid-10, only early sampled (d1-5) in assay, full data in metadata",
                  full_counts = counts,
                  full_meta = meta))


rownames(se) <- counts_filtered$gene

se <- se[!is.na(rownames(se)),]
saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))


meta |>
  count(group)
