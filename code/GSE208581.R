library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)


GEO_accs <- "GSE208581"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE208581&format=file&file=GSE208581%5Fraw%5Fcounts%2Etsv%2Egz"

dir_path <- here::here("data", GEO_accs)
file_path <- here::here(dir_path, paste0(GEO_accs, ".tsv.gz"))



if(!file.exists(file_path)){
  dir.create(dir_path)
  curl::curl_download(download,
                      destfile = file_path)
  }


sm <- getGEO(GEO_accs, destdir = dir_path) 

glimpse(pData(sm[[1]]))
unique(pData(sm[[1]])$`condition:ch1`)

as.data.frame(pData(sm[[1]])) |>
  count(`condition:ch1`)

meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = `description`,
                sex = `Sex:ch1`,
                source = `tissue:ch1`
                ) |>
  mutate(disease = case_when(`condition:ch1` == "SIRS+" ~ "non-infectious_sirs",
                             `condition:ch1` == "SIRS-" ~ "uncomplicated-post-op",
                             `condition:ch1` == "Sepsis+" ~ "post-op_sepsis",
                             `condition:ch1` == "UInf+" ~ "simple_post-op_infection"),
         age = NA,
         dataset = GEO_accs, 
         pediatric = FALSE, 
         source = "whole_blood") |>

  dplyr::select(id, sample_name, pediatric,
                disease, 
                processing_info, source, dataset, age, sex) 

counts <- fread(file_path, data.table = FALSE) |> rename("gene" = "Ensembl_Gene_ID")


ids <- meta$sample_name
names(ids) <- meta$id

counts <- counts |> 
  dplyr::select(gene, all_of(meta$sample_name)) |>
  rename(all_of(ids))

all(colnames(counts)[-1] == meta$id)


se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts[,-1])),
  colData = meta,
  metadata = list(notes = "Patients with post-operative infection, post-op sepsis, non-inf SIRS  and normal post-op"))

all(floor(assay(se)) == assay(se), na.rm = TRUE)

rownames(se) <- counts$gene

se <- se[!is.na(rownames(se)),]
saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))
