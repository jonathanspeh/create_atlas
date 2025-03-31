library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)


# Failed - only normalised counts
GEO_accs <- "GSE196822"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE196822&format=file&file=GSE196822%5FRaw%5Fcounts%5Fmatrix%2Ecsv%2Egz"
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

unique(pData(sm[[1]])$`condition:ch1`)




meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = `sequencing id:ch1`,
                sex = `Sex:ch1`,
                source = `source_name_ch1`) |>
  mutate(disease = case_when(`condition:ch1` == "Healthy" ~ "healthy",
                             grepl("Coin", `condition:ch1`) ~ "COVID-19_and_bacterial_coinfection",
                             TRUE ~ "COVID-19"),
         severity = stringr::str_remove(`condition:ch1`, "\\s*Covid-19\\s*"),
         age = floor(as.numeric(`age:ch1`)),
         dataset = GEO_accs, 
         pediatric = age < 18) |>
  dplyr::select(id, sample_name, pediatric,
                disease, 
                processing_info, source, dataset, age, sex, severity) 

counts <- fread(file_path) 

counts <- counts |> 
  dplyr::select(gene, any_of(meta$sample_name)) #

meta <- meta |>
  dplyr::filter(sample_name %in% colnames(counts))

ids <- meta$sample_name
names(ids) <- meta$id

counts <- counts |>
  rename(any_of(ids))
  
rownames(meta) <- meta$id

all(rownames(meta) == colnames(counts)[-1])

se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts[,-1])),
  colData = meta,
  metadata = list(notes = "Patients with covid-19 and other respiratory diseases"))

all(floor(assay(se)) == assay(se), na.rm = TRUE)

rownames(se) <- counts$gene

se <- se[!is.na(rownames(se)),]
#saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))
