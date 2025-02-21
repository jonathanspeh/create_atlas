library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)

#####USELESS - COUNTS DONT CONTAIN GENE NAMES #####


GEO_accs <- "GSE152641"
file <- "GSE152641_Inflammatix_COVID19_counts_entrez.csv.gz"
dir_path <- paste0(here::here("data", GEO_accs), "/")
dir.create(dir_path)
#curl::curl_download("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE152641&format=file&file=GSE152641%5FInflammatix%5FCOVID19%5Fcounts%5Fentrez%2Ecsv%2Egz",
#                    destfile = here::here(dir_path, file))


sm <- getGEO(GEO_accs, destdir = dir_path)

meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::select(id = geo_accession,
                sample_name = title,
                sample_type =  characteristics_ch1.2,
                age = "age:ch1",
                disease = "disease:ch1",
                sex = "Sex:ch1",
                group = "source_name_ch1",
                source = "tissue:ch1",
                processing_info) |> 
  mutate(individual = sub(".*(\\d{2})\\D*$", "\\1", sample_name),
         disease = case_when(grepl("control", disease) ~ "healthy",
                             TRUE ~ "COVID-19"),
         dataset = GEO_accs,
         pediatric = FALSE
         )


files <- list.files(dir_path, full.names = FALSE)
files <- files[grepl("GSM", files)]



counts <- fread(here::here(dir_path, file))


# se <- SummarizedExperiment(
#   assays = list(counts = as.matrix(counts_raw[,-1])),
#   colData = meta,
#   metadata = list(notes = "covid 19 (Omicron) patienst with ans without previous challenge (infection / vaccine) and healthy controls")
# )
# 
# rownames(se) <- mapIds(EnsDb.Hsapiens.v86,
#                        keys = counts_raw$gene,
#                        keytype = "SYMBOL",
#                        columns = "ENSEMBL")
# 
# se <- se[!is.na(rownames(se)),]
# saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))



