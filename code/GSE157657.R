library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)



# Failed - couldn't make sens of metadat
GEO_accs <- "GSE157657"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE157657&format=file&file=GSE157657%5Fnorm%2Edata%2Etxt%2Egz"

dir_path <- here::here("data", GEO_accs)
file_path <- here::here(dir_path, paste0(GEO_accs, ".tsv.gz"))



if(!file.exists(file_path)){
  dir.create(dir_path)
  curl::curl_download(download,
                      destfile = file_path)
  }


sm <- getGEO(GEO_accs, destdir = dir_path) 

glimpse(pData(sm[[1]]))
unique(pData(sm[[1]])$`days_from_att:ch1`)

as.data.frame(pData(sm[[1]])) |>
  count(`patient id:ch1`, `characteristics_ch1.2`, `characteristics_ch1.1`,)

NA < 1

as.data.frame(pData(sm[[1]])) |> 
  filter(`days_from_att:ch1` < 1 | is.na(`days_from_att:ch1`)) |>
  count(`characteristics_ch1.2`, `characteristics_ch1.1`)

meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = `description`,
                sex = `gender:ch1`,
                age = `age:ch1`
                source = `tissue:ch1`,
                patient = `patient id:ch1`
                ) |>
  mutate(disease = case_when(`condition:ch1` == "SIRS+" ~ "non-infectious_sirs",
                             `condition:ch1` == "SIRS-" ~ "uncomplicated-post-op",
                             `condition:ch1` == "Sepsis+" ~ "post-op_sepsis",
                             `condition:ch1` == "UInf+" ~ "simple_post-op_infection"),
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
