library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)
library(stringr)

# Not sure if witrht to be included

GEO_accs <- "GSE155925"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE155925&format=file&file=GSE155925%5FRaw%5Fcounts%5Fmatrix%2Etxt%2Egz"
dir_path <- here::here("data", GEO_accs)
file_path <- here::here(dir_path, paste0(GEO_accs, ".txt.gz"))

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
                sample_name = title,
                sex = "Sex:ch1",
                source = "tissue:ch1"
  ) |> 
  mutate(disease = case_when(`pathogen:ch1` =="negative" ~ "healthy",
                             `pathogen:ch1` =="RSV" ~ "respiratory_syncytial_virus_infection",
                             `pathogen:ch1` =="HRV" ~ "rhinovirus_infection",
                             `pathogen:ch1` =="HMPV" ~ "metapneumovirus_infection",
                             `pathogen:ch1` =="FLU" ~ "Influenza",
                             `pathogen:ch1` =="Adenovirus" ~ "adenovirus_infection",
                             `pathogen:ch1` =="RSV,HRV" ~ "rsv and rhinovirus infection",
                             `pathogen:ch1` =="RSV,EV-D68" ~ "rsv and enterovirus D-86 infection",
                             TRUE ~ `pathogen:ch1`),
         dataset = GEO_accs, 
         pediatric = TRUE,
         age_month = as.numeric(`age (months):ch1`),
         age = round(as.numeric(age_month) / 12)
         #      processing_info = paste(
         # "growth_protocol:", growth_protocol_ch1,
         # "extract_protocol:", extract_protocol_ch1,
         # "library_prep:", extract_protocol_ch1.1, 
         # "data_processing_1:", data_processing,
         # "data_processing_2:", data_processing.1,
         # "assembly:", data_processing.2,
         # "instrument:", instrument_model,
         # sep = "\t")
  ) |>
  dplyr::select(id, sample_name, age, age_month, pediatric, sex,
                disease, 
                processing_info, source, dataset)


counts <- fread(file_path) |>
  separate_wider_delim(V1, delim = ":", names = c("gene.name", "ENS.Symbol")) |> 
  rename_at(vars(meta$sample_name), ~ meta$id)


rownames(meta) <- meta$id
all(rownames(meta) == colnames(counts)[-c(1,2)])

se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts[,-c(1,2)])),
  colData = meta,
  metadata = list(notes = "Pediatric patients with different infections and healthy controls")
)


rownames(se) <- counts$ENS.Symbol

se <- se[!is.na(rownames(se)),]
saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))
