library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)
library(stringr)


GEO_accs <- "GSE69529"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE69529&format=file&file=GSE69529%5Fmexico%5Fprocessed%5Fdata%5Ffor%5FGEO%2Ecsv%2Egz"
download_2 <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE69529&format=file&file=GSE69529%5FGSM2241808%2DGSM2241855%5Fmexico2%5Fdata%2Ecsv%2Egz"
dir_path <- here::here("data", GEO_accs)
file_path <- here::here(dir_path, paste0(GEO_accs, ".csv.gz"))
file_path_2 <- here::here(dir_path, paste0(GEO_accs, "_2.csv.gz"))


if(!file.exists(file_path)){
  dir.create(dir_path)
  curl::curl_download(download,
                      destfile = file_path)
  curl::curl_download(download_2,
                      destfile = file_path_2)
  }



sm <- getGEO(GEO_accs, destdir = dir_path) 

meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = title,
                age_group = "age:ch1",
                pathogen = `organisms:ch1`,
                severity = `severity:ch1`, 
                case_or_control = `group:ch1`
  ) |>
  mutate(path_group = ifelse(grepl("EC", pathogen), "Escherichia coli", pathogen),
         group = paste(case_or_control, pathogen, severity, sep = "_"),
         disease = case_when(case_or_control == "Control" & pathogen == "Negative" ~ "healthy",
                             case_or_control == "Control" & pathogen != "Negative" ~ paste(path_group, "colonisation", sep = "_"),  # change path_group to pathogen if E coli strain is required
                             case_or_control == "Diarrhea" ~ paste(path_group, "infection", sep = "_")),  # change path_group to pathogen if E coli strain is required
         dataset = GEO_accs, 
         pediatric = TRUE,
         source = "whole blood",
         sex = NA,
         age = NA
         ) |>
  dplyr::select(id, sample_name, age_group, pediatric,
                disease, 
                processing_info, source, dataset, pathogen, group, case_or_control, severity, age, sex)


counts <- fread(file_path)
counts_2 <- fread(file_path_2)

counts_all <- left_join(counts, counts_2, 
                        by = join_by("Ensembl Gene ID" == "Ensembl.Gene.ID")) 

meta_filtered <- meta |>
  dplyr::filter(sample_name %in% colnames(counts_all) & !grepl("repeat", sample_name))

counts_selected <- counts_all |>
  dplyr::select(`Ensembl Gene ID`, meta_filtered$sample_name) |>
  rename_at(vars(meta_filtered$sample_name), ~ meta_filtered$id)

rownames(meta_filtered) <- meta_filtered$id


se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts_selected[,-1])),
  colData = meta_filtered,
  metadata = list(notes = "Pediatric patients with different causes of diarrhea",
                  full_counts = counts_all,
                  fill_meta = meta))

all(floor(assay(se)) == assay(se), na.rm = TRUE)

rownames(se) <- counts_selected$`Ensembl Gene ID`

se <- se[!is.na(rownames(se)),]
saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))



