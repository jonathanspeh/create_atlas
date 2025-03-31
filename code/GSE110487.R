library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)


GEO_accs <- "GSE110487"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE110487&format=file&file=GSE110487%5Frawcounts%2Exlsx"

dir_path <- here::here("data", GEO_accs)
file_path <- here::here(dir_path, paste0(GEO_accs, ".xlsx"))



if(!file.exists(file_path)){
  dir.create(dir_path)
  curl::curl_download(download,
                      destfile = file_path)
  }


sm <- getGEO(GEO_accs, destdir = dir_path) 


meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = `title`,
                timepoint = `timepoint:ch1`,
                patient = `patient:ch1`) |>
  mutate(disease = "septic_shock",
         age = NA,
         sex = NA,
         dataset = GEO_accs, 
         pediatric = FALSE, 
         source = "whole_blood") |>

  dplyr::select(id, sample_name, pediatric,
                disease, 
                processing_info, source, dataset, age, sex, timepoint) 


meta_2 <- pData(sm[[2]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = `title`,
                timepoint = `timepoint:ch1`,
                patient = `patient:ch1`) |>
  mutate(disease = "septic_shock",
         age = NA,
         sex = NA,
         dataset = GEO_accs, 
         pediatric = FALSE, 
         source = "whole_blood") |>
  
  dplyr::select(id, sample_name, pediatric,
                disease, 
                processing_info, source, dataset, age, sex, timepoint) 


meta <- rbind(meta, meta_2)


counts <- readxl::read_xlsx(file_path) |> rename("gene" = "Geneid")

ids <- meta$sample_name
names(ids) <- meta$id


counts <- counts |> 
  dplyr::select(gene, all_of(meta$sample_name)) |>
  rename(all_of(ids))

all(colnames(counts)[-1] == meta$id)

meta_t1 <- meta |>
  filter(timepoint == "T1") 

counts_t1 <- counts |>
  dplyr::select(gene, all_of(meta_t1$id))
  
rownames(meta_t1) <- meta_t1$id

all(rownames(meta_t1) == colnames(counts_t1)[-1])

se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts_t1[,-1])),
  colData = meta_t1,
  metadata = list(notes = "Patients with Septic shock, only Timepoint 1 in assay",
                  full_meta = counts,
                  full_meta = meta))

all(floor(assay(se)) == assay(se), na.rm = TRUE)

rownames(se) <- counts_t1$gene

se <- se[!is.na(rownames(se)),]
saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))
