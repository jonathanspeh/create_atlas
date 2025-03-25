library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)
library(stringr)

GEO_accs <- "GSE255053"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE255053&format=file"
dir_path <- here::here("data", GEO_accs)
tar_path <- here::here("data", paste0(GEO_accs, ".tar"))


if(!file.exists(tar_path)){
  dir.create(dir_path)
  curl::curl_download(download,
                      destfile = tar_path)
  untar(tar_path, 
        exdir = dir_path)
}


sm <- getGEO("GSE255053", destdir = dir_path) # can't download smf?


meta <- pData(sm[[1]]) |>
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                source = "tissue:ch1",) |> 
  mutate(individual = stringr::str_match(title, "HC[0-9]*"),
         timepoint = case_when(
           grepl("baseline", title) ~ "baseline",
           grepl("during", title) ~ "during infection",
           grepl("recovery", title) ~ "after malaria"),
         disease = case_when(
           grepl("baseline", title) ~ "healthy",
           grepl("during", title) ~ "Malaria",
           grepl("recovery", title) ~ "convalescent"),
         dataset = GEO_accs, 
         pediatric = TRUE,
         sample_name = title,
         age = NA,
         sex = NA
         ) |>
  dplyr::select(id, sample_name, pediatric,
                disease, age, sex,
                processing_info, source, dataset)






files <- list.files(dir_path, full.names = FALSE)
files <- files[grepl("GSM", files)]


read_counts <- function(name){
  sample <- stringr::str_remove(name, "_.*.txt.gz")
  file <- data.table::fread(here::here(dir_path, name)) |>
    dplyr::filter(!grepl("^__", V1)) 
  colnames(file) <- c("gene", substitute(sample))
  file
}

counts_list <- lapply(files, read_counts)
counts <- purrr::reduce(counts_list, full_join, by = "gene") 


meta_active = filter(meta, disease == "Malaria")
counts_active = select(counts, all_of(c("gene", meta_active$id)))


se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts_active[,-1])),
  colData = meta_active,
  metadata = list(notes = "Pediatric patients Malaria, only active malaria in assay, baseline and reconvery in metadat",
                  full_counts = counts,
                  full_meta = meta))



rownames(se) <- str_remove(counts_active$gene, "\\.[0-9]*") # gene names have some suffices

se <- se[!is.na(rownames(se)),]
saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))



