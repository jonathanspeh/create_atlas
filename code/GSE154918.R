library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)

# Failed - no normalised data
GEO_accs <- "GSE154918"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE154918&format=file&file=GSE154918%5FSchughart%5FSepsis%5F200320%2Etxt%2Egz"

dir_path <- here::here("data", GEO_accs)
file_path <- here::here(dir_path, paste0(GEO_accs, ".txt.gz"))



if(!file.exists(file_path)){
  dir.create(dir_path)
  curl::curl_download(download,
                      destfile = file_path)
  }


sm <- getGEO(GEO_accs, destdir = dir_path) 

glimpse(pData(sm[[1]]))
unique(pData(sm[[1]])$`status:ch1`)


meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = `title`,
                sex = `Sex:ch1`) |>
  mutate(disease = case_when(`status:ch1` == "Hlty" ~ "healthy",
                             `status:ch1` == "Inf1_P" ~ "uncomplicated_infection",
                              grepl("Seps", `status:ch1`) ~ "Sepsis", 
                              grepl("Shock", `status:ch1`) ~ "septic_shock"),
         timepoint = ifelse(grepl("FU", `status:ch1`), "follow_up", "acute"),
         age = NA,
         dataset = GEO_accs, 
         pediatric = FALSE, 
         source = "whole_blood") |>

  dplyr::select(id, sample_name, pediatric,
                disease, 
                processing_info, source, dataset, age, sex, timepoint) 

count(meta, disease, timepoint)
meta$sample_name
counts <- fread(file_path, data.table = FALSE) |> rename("gene" = "ENSEMBL_gene_ID")

ids <- meta$sample_name
names(ids) <- meta$id

counts <- counts |> 
  dplyr::select(gene, all_of(meta$sample_name)) |>
  rename(all_of(ids))

all(colnames(counts)[-1] == meta$id)


se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts[,-1])),
  colData = meta,
  metadata = list(notes = "Patients with spesis and healthy control"))

all(floor(assay(se)) == assay(se), na.rm = TRUE)

rownames(se) <- counts$gene

se <- se[!is.na(rownames(se)),]
#saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))
