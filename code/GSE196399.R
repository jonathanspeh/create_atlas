library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)

GEO_accs <- "GSE196399"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE196399&format=file&file=GSE196399%5Fcount%5Fmatrix%2Ecsv%2Egz"

dir_path <- here::here("data", GEO_accs)
file_path <- here::here("data", paste0(GEO_accs, ".csv.gz"))


if(!file.exists(file_path)){
  dir.create(dir_path)
  curl::curl_download(download_csv,
                    destfile = file_path)
}

sm <- getGEO(GEO_accs, destdir = dir_path) 

meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = `title`,
                source = `cell type:ch1`,
                ) |>
  mutate(disease = ifelse(grepl("SCAP", `disease state:ch1`), "pneumonia", "healthy"),
         age = NA,
         sex = NA,
         dataset = GEO_accs, 
         pediatric = TRUE
         ) |>
  dplyr::select(id, sample_name, pediatric,
                disease, 
                processing_info, source, dataset, age, sex) 



counts <- fread(file_path) |> rename("gene" = "V1")

ids <- meta$sample_name
names(ids) <- meta$id

counts <- counts |> 
  rename(all_of(ids))

all(colnames(counts)[-1] == meta$id)
count(meta, disease)

se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts[,-1])),
  colData = meta,
  metadata = list(notes = "Patients with severe CAP and healthy controls"))

all(floor(assay(se)) == assay(se), na.rm = TRUE)

rownames(se) <- counts$gene

se <- se[!is.na(rownames(se)),]
saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))
