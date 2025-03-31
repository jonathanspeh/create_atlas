library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)

GEO_accs <- "GSE261482"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE261482&format=file&file=GSE261482%5FCounts%5Fraw%5Fdata%2Ecsv%2Egz"

dir_path <- here::here("data", GEO_accs)
file_path <- here::here(dir_path, paste0(GEO_accs, ".csv.gz"))



if(!file.exists(file_path)){
  dir.create(dir_path)
  curl::curl_download(download,
                      destfile = file_path)
  }


sm <- getGEO(GEO_accs, destdir = dir_path) 

glimpse(pData(sm[[1]]))
unique(pData(sm[[1]])$`etiology:ch1`)
as_tibble(pData(sm[[1]])) |>
  count(`condition:ch1`, `etiology:ch1`)


meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = `title`,
                sex = `Sex:ch1`,
                source = `tissue:ch1`,
                sex = `Sex:ch1`
                ) |>
  mutate(#disease = ifelse(`condition:ch1` == "Control", "healthy", `condition:ch1`),
         disease = case_when(`etiology:ch1` == "PB" ~ "pneumonia_probably_bacterial",
                              `etiology:ch1` == "DB" ~ "pneumonia_bacterial",
                              `etiology:ch1` == "PV" ~ "pneumonia_probably_viral",
                              `etiology:ch1` == "DV" ~ "pneumonia_viral",
                              `etiology:ch1` == "U" ~ "pneumonia",
                              `etiology:ch1` == "C" ~ "healthy",
                              ),
         sample_num = stringr::str_match(sample_name, "(\\d+)")[,1], 
         age = NA,
         dataset = GEO_accs, 
         pediatric = TRUE
         ) |>
  dplyr::select(id, sample_name, sample_num, pediatric,
                disease, 
                processing_info, source, dataset, age, sex) 


counts <- fread(file_path, header = TRUE) |> rename("gene" = "Sample number")

ids <- meta$sample_num
names(ids) <- meta$id

counts <- counts |> 
  rename(all_of(ids))

all(colnames(counts)[-1] == meta$id)


se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts[,-1])),
  colData = meta,
  metadata = list(notes = "Pediatric patients with different causes of pneumonia"))

all(floor(assay(se)) == assay(se), na.rm = TRUE)

rownames(se) <- counts$gene

se <- se[!is.na(rownames(se)),]
#saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))
