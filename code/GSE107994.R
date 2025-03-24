library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)


GEO_accs <- "GSE107994"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE107994&format=file&file=GSE107994%5FRaw%5Fcounts%5FLeicester%5Fwith%5Fprogressor%5Flongitudinal%2Exlsx"
file <- paste0(GEO_accs, ".xlsx")
dir_path <- paste0(here::here("data", GEO_accs), "/")

dir.create(dir_path)
curl::curl_download(download,
                    destfile = here::here(dir_path, file))

sm <- getGEO(GEO_accs, destdir = dir_path)


unique(pData(sm[[1]])$characteristics_ch1.1)


meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = title,
                group =  `characteristics_ch1.1`,
                timepoint = `characteristics_ch1.12`,
                source = "tissue:ch1") |> 
  mutate(individual = stringr::str_remove(`characteristics_ch1.3`, "patient_id: "),
         disease = case_when(
           grepl("Active", group) ~ "Tuberculosis",
           grepl("Progressor", group) ~ "progressive latent tuberculosis",
           grepl("LTBI", group) ~ "latent Tuberculosis",
           grepl("Control", group) ~ "healthy"),
         tb_type = stringr::str_remove(`characteristics_ch1.4`, "tb_disease_type: "),
         dataset = GEO_accs,
         pediatric = FALSE,
         sample_type = paste(source_name_ch1, timepoint, collapse = "_"),
         age = stringr::str_remove(`characteristics_ch1.10`, "age_at_baseline_visit: "),
         sex = stringr::str_remove(`characteristics_ch1.9`, "gender: ")) |>
  dplyr::select(id, individual, sample_name, sample_type, age, pediatric,sex, group, 
                disease, processing_info, source, dataset, timepoint)


counts_raw <- readxl::read_excel(here::here(dir_path, file))

meta$timepoint

# Filter on baseline
meta_bl <- dplyr::filter(meta, grepl("Baseline", timepoint))
counts_bl <-  dplyr::select(counts_raw, c(Genes, meta_bl$sample_name))


se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts_bl[,-1])),
  colData = meta_bl,
  metadata = list(notes = "Latent and active TB and controls, latent progressors",
                  full_counts = counts_raw,
                  full_meta = meta))
  


rownames(se) <- counts_raw$Genes
se <- se[!is.na(rownames(se)),]
saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))

