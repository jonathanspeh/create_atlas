library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)


# Failed - only normalised counts
GEO_accs <- "GSE211567"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE211567&format=file&file=GSE211567%5FnormData%5Fdiscovery%5F2021MAR24%2Etxt%2Egz"
dir_path <- here::here("data", GEO_accs)
file_path <- here::here(dir_path, paste0(GEO_accs, ".csv.gz"))


if(!file.exists(file_path)){
  dir.create(dir_path)
  curl::curl_download(download,
                      destfile = file_path)
  }



sm <- getGEO(GEO_accs, destdir = dir_path) 



pData(sm[[1]]) |> dplyr::select(starts_with("char")) 
pData(sm[[1]]) |>
  as_tibble() |>
  dplyr::filter(`infection type:ch1` == "Noninfection") 
colnames(pData(sm[[1]]))
glimpse(pData(sm[[1]]))

unique(pData(sm[[1]])$`pathogen:ch1`)



meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = title,
                sex = "gender:ch1",
                pathogen = `pathogen:ch1`,
                group = `infection type:ch1`
  ) |>
  mutate(disease = case_when(pathogen == "Noninfection" ~ "non-infectious disease",
                             TRUE ~ paste(pathogen, "infection", sep = "_")), 
         age = floor(as.numeric(`age:ch1`)),
         dataset = GEO_accs, 
         pediatric = age < 18,
         source = "whole blood") |>
  dplyr::select(id, sample_name, pediatric,
                disease, 
                processing_info, source, dataset, pathogen, group, age, sex)



counts <- fread(file_path) |> rename("gene" = V1)


rownames(meta) <- meta$id


se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts_selected[,-1])),
  colData = meta_filtered,
  metadata = list(notes = "atients with different causes of respiratory diseases"))

all(floor(assay(se)) == assay(se), na.rm = TRUE)

rownames(se) <- counts_selected$`Ensembl Gene ID`

se <- se[!is.na(rownames(se)),]
#saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))





