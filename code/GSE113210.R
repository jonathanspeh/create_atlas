library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)
library(stringr)

#. Failed -only transcipt level annotations in GEO object -  could use rawer data


GEO_accs <- "GSE113210"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE113210&format=file&file=GSE113210%5FAVB%5FPBMC%5Fraw%5Fcounts%2Ecsv%2Egz"
file <- paste0(GEO_accs, ".csv.gz")
dir_path <- paste0(here::here("data", GEO_accs), "/")

dir.create(dir_path)
curl::curl_download(download,
                    destfile = here::here(dir_path, file))


sm <- getGEO(GEO_accs, destdir = dir_path)

pData(sm[[1]]) |> dplyr:: select(starts_with("char"))

last_n <- function(str, n=1){
  substr(str, nchar(str)-n+1, nchar(str))
}


pData(sm[[1]]) |>
  reframe(subject = characteristics_ch1,
          n_coinfections = last_n(characteristics_ch1.7),
          rsv_pos = last_n(characteristics_ch1.8),
          rhinovirus_pos = last_n(characteristics_ch1.9),
          adenovirus_pos = last_n(characteristics_ch1.10),
          influenzavirus_pos = last_n(characteristics_ch1.11),
          coronavirus_pos = last_n(characteristics_ch1.12),
          polyomavirus_pos = last_n(characteristics_ch1.13),
          parainfluenzavirus_pos = last_n(characteristics_ch1.14),
          bocavirus_pos = last_n(characteristics_ch1.15),
          enterovirus_pos = last_n(characteristics_ch1.16),
          metapneumovirus_pos = last_n(characteristics_ch1.17))



meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = title,
                sample_type =  `characteristics_ch1.3`,
                group = "characteristics_ch1.2",
                timepoint = `characteristics_ch1.1`) |> 
  mutate(individual = stringr::str_remove(characteristics_ch1, "patient: "),
         disease = case_when(stringr::str_remove(characteristics_ch1.8, "rsv positive at av (1=yes, 0=no, 9=not measured): ") == 1 ~"RSV",
                             
                             
                             
         )
         dataset = GEO_accs,
         age = stringr::str_remove(characteristics_ch1.1, "age at visit (years): ")
         pediatric = TRUE,
         sex = ifelse(stringr::str_remove("characteristics_ch1.4", "sex(0=female, 1=male): ") == 0, "female", "male"),
         source = "PBMC"
  ) |>
  dplyr::select(id, individual, sample_name, sample_type, age, pediatric, sex, group, 
                disease, processing_info, source, dataset, timepoint)
