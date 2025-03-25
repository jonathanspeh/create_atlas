library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)
library(stringr)


GEO_accs <- "GSE158163"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE158163&format=file"
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
                sex = "gender:ch1",
                age = "age:ch1",
                group = "group:ch1"
  ) |> 
  mutate(disease = case_when(grepl("(sGAS)", description) ~ "group a streptococcal infection",
                             grepl("(vCtrl)", description) ~ "asymptomatic infection",
                             grepl("(sGAS)", description) ~ "group a streptococcal colonisation",
                             grepl("(nCtrl)", description) ~ "healthy",
                             grepl("(sAdV)", description) ~ "adenovirus infection",
                             grepl("(sOP)", description) ~ "pharyngitis unknown pathogen",
                             TRUE ~ description),
         sample_name = stringr::str_match(title, "\\[(.*)\\]")[,2] ,
         dataset = GEO_accs, 
         pediatric = as.numeric(age) < 18,
         source = "whole blood"
         ) |>
  dplyr::select(id, sample_name, age, pediatric, sex,
                disease, 
                processing_info, source, dataset)



files <- list.files(dir_path, full.names = FALSE)
files <- files[grepl("GSM", files)]


read_counts <- function(name){
  sample <- stringr::str_remove(name, "_.*.txt.gz")
  file <- data.table::fread(here::here(dir_path, name)) #|>
    #dplyr::filter(!grepl("^__", V1)) 
  colnames(file) <- c("gene", substitute(sample))
  file
}

counts_list <- lapply(files, read_counts)
counts <- purrr::reduce(counts_list, full_join, by = "gene") 


se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts[,-1])),
  colData = meta,
  metadata = list(notes = "Pediatric patients with different causes of pharyngitis and healthy controls")
)


rownames(se) <- counts$gene
colnames(se) %in% meta$id
se <- se[!is.na(rownames(se)),]
saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))





