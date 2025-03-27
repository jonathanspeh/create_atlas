library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)


GEO_accs <- "GSE107991"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE107991&format=file&file=GSE107991%5FRaw%5Fcounts%5FBerry%5FLondon%2Exlsx"
file <- paste0(GEO_accs, ".xlsx")
dir_path <- paste0(here::here("data", GEO_accs), "/")

dir.create(dir_path)
curl::curl_download(download,
                    destfile = here::here(dir_path, file))

sm <- getGEO(GEO_accs, destdir = dir_path)
#sm <- getGEO(filename = here::here(dir_path, "GSE107991_series_matrix.txt.gz"))


meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = title,
                sample_type =  characteristics_ch1.1,
                group = "group:ch1",
                source = "tissue:ch1") |> 
  mutate(individual = stringr::str_remove(sample_name, "Berry_London_Sample"),
         disease = case_when(
           grepl("Active", sample_type) ~ "Tuberculosis",
           grepl("LTBI", sample_type) ~ "latent Tuberculosis",
           grepl("Control", sample_type) ~ "healthy"),
         dataset = GEO_accs,
         pediatric = FALSE,
         age = NA,
         sex = NA) |>
  dplyr::select(id, individual, sample_name, sample_type, age, pediatric,sex, group, 
                disease, processing_info, source, dataset)


counts_raw <- readxl::read_excel(here::here(dir_path, file)) |> 
  rename_at(vars(meta$sample_name), ~ meta$id)

rownames(meta) <- meta$id

all(rownames(meta) == colnames(counts_raw)[-c(1,2,3)])


# Filter only active tb?
meta_active <- dplyr::filter(meta, disease == "Tuberculosis")
counts_active <-  dplyr::select(counts_raw, c(Genes, meta_active$id))

all(rownames(meta_active) == colnames(counts_active)[-c(1)])



se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts_raw[,-c(1,2,3)])),
  colData = meta,
  metadata = list(notes = "Latent and active TB and controls")
  
)

rownames(se) <- counts_raw$Genes
se <- se[!is.na(rownames(se)),]
saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))

