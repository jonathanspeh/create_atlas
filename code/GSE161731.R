library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)


# Failed - only normalised counts
GEO_accs <- "GSE161731"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE161731&format=file&file=GSE161731%5Fcounts%2Ecsv%2Egz"
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

unique(pData(sm[[1]])$`title`)




meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = title,
                sex = "gender:ch1",
                timepoint = `characteristics_ch1.5`,
                group = `characteristics_ch1.5`,
                source = `tissue:ch1`
  ) |>
  mutate(disease = case_when(`cohort:ch1` == "Bacterial" ~ "bacterial_pneumonia",
                             `cohort:ch1` == "CoV other" ~ "seasonal_coronavirus_infection",
                             TRUE ~ `cohort:ch1`), 
         age = floor(as.numeric(`age:ch1`)),
         dataset = GEO_accs, 
         pediatric = ifelse(is.na(age), FALSE, age < 18),   # has on adult with unknown age
         source = "whole blood",
         sample_name = stringr::str_remove(sample_name, "whole blood, ")
         ) |>
  dplyr::select(id, sample_name, pediatric,
                disease, 
                processing_info, source, dataset, group, age, sex) |>
  separate_longer_delim(sample_name, ", ") # Has some strange duplication going on



counts <- fread(file_path) |> rename("gene" = V1)


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
saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))






