library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)
library(stringr)

# Not sure if witrht to be included

GEO_accs <- "E-MTAB-13307"
download <- "https://www.ebi.ac.uk/biostudies/files/E-MTAB-13307/gene_count_matrix.tsv.gz"
sdrf <- "https://www.ebi.ac.uk/biostudies/files/E-MTAB-13307/E-MTAB-13307.sdrf.txt"
dir_path <- here::here("data", GEO_accs)
file_path <- here::here(dir_path, paste0(GEO_accs, ".tsv.gz"))
sdrf_path <- here::here(dir_path, paste0(GEO_accs, ".sdrf.txt"))

if(!file.exists(file_path)){
  dir.create(dir_path)
  curl::curl_download(download,
                      destfile = file_path)
  curl::curl_download(sdrf,
                      destfile = sdrf_path)
  }


counts <- fread(file_path)  |> rename("gene" = "V1")
meta <- fread(sdrf_path) 

colnames(meta) <- janitor::make_clean_names(colnames(meta))
meta$source_name <- janitor::make_clean_names(meta$source_name)
colnames(counts) <- paste0(GEO_accs, "_", janitor::make_clean_names(colnames(counts)))

meta$age

meta <- meta[rep(c(T,F), nrow(meta)/2),] |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(
    id = source_name,
    age = characteristics_age,
    sex = characteristics_sex) |>
  mutate(
    sample_name = id,
    id = paste0(GEO_accs, "_", id),
    disease = case_when(
      factor_value_pathogen == "HC" ~ "healthy",
      factor_value_pathogen == "SARS-CoV-2" ~ "COVID-19",
      factor_value_pathogen == "Measles" ~ "Measles",
      factor_value_pathogen == "Escherichia coli (ESBL)" ~ "Escherichia coli infection",
      factor_value_pathogen == "Herpes zoster myeloma" ~ "Herpes zoster myeloma",
      TRUE ~ paste(factor_value_pathogen, "infection")),
    dataset = GEO_accs,
    pediatric = characteristics_developmental_stage != "adult", 
    source = "PBMC"
  ) |>
  dplyr::select(id, sample_name, age, sex, pediatric, disease, processing_info, source, dataset)


se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts[,-1])),
  colData = meta,
  metadata = list(notes = "various diseases, adult patients")
)    


rownames(se) <- counts$`E-MTAB-13307_gene`

se <- se[!is.na(rownames(se)),]

saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))
