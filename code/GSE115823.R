library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)
library(stringr)

#. Failed -only transcipt level annotations in GEO object -  could use rawer data


GEO_accs <- "GSE115823"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE115823&format=file&file=GSE115823%5Fraw%5Fcounts%5Fblood%5Fmuppits511%2Etxt%2Egz"
dir_path <- here::here("data", GEO_accs)
file_path <- here::here(dir_path, paste0(GEO_accs, ".txt.gz"))




if(!file.exists(file_path)){
  dir.create(dir_path)
  curl::curl_download(download,
                      destfile = file_path)
}


sm <- getGEO(GEO_accs, destdir = dir_path)

# Visist 0 --> Baseline, no infection, no exacerbation



 pData(sm[[1]]) |> dplyr:: select(geo_accession, starts_with("char")) |> 
  pivot_longer(starts_with("char"), names_to = NULL) |>
   filter(value != "") |>
  separate_wider_delim(value, delim = ": ", names = c("name", "value")) |>
   pivot_wider(names_from = name, values_from = value) |>
   select(age.in.years, Sex, analysis.visit, case.or.control.status.event, csteroid.start.relative.to.visit, starts_with("virus")) |>
   #filter(analysis.visit != "Visit 0") |>
   print(n = 550)
   
    #select(starts_with("virus")) -> virmat
  

 

meta <- pData(sm[[1]]) |> 
  as_tibble(rownames = "id")
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(sample_name = title,
                sample_type =  characteristics_ch1.1,
                group = "group:ch1",
                source = "tissue:ch1") |> 
  mutate(individual = stringr::str_remove(sample_name, "Berry_SouthAfrica_Sample"),
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








