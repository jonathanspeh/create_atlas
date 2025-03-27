library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)


GEO_accs <- "GSE205244"
file <- paste0(GEO_accs, ".tar")
dir_path <- paste0(here::here("data", GEO_accs), "/")

# curl::curl_download(paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=", 
#                            GEO_accs, 
#                            "&format=file"),
#                     destfile = here::here("data", file))
# 
# untar(here::here("data", file), 
#       exdir = dir_path)

sm <- getGEO(filename = here::here(dir_path, "GSE205244_series_matrix.txt.gz"))


meta <- pData(sm) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = title,
                sample_type =  characteristics_ch1.2,
                age = "age:ch1",
                sex = "gender:ch1",
                sampling_time = "days after positive pcr results:ch1",
                group = "disease state:ch1",
                omicron_lineage = "omicron sublineage:ch1",
                source = "cell type:ch1") |> 
  mutate(individual = sub("^.*_(\\d+)_.*$", "\\1", sample_name),
         sampling_point = sub("^.*?(\\d+)\\D*$", "\\1", sample_name), 
         disease = "COVID-19",
         variant = "Omicron",
         dataset = GEO_accs,
         pediatric = as.numeric(age < 18)
         ) |>
    dplyr::select(id, individual, sample_name, sample_type, age, pediatric,sex, group, 
                  omicron_lineage, sampling_time, sampling_point, disease, variant, 
                  processing_info, source, dataset)


files <- list.files(dir_path, full.names = FALSE)
files <- files[grepl("GSM", files)]

counts <- lapply(files, read_counts)
counts_raw <- purrr::reduce(counts, full_join, by = "gene") 


meta_point1 <- dplyr::filter(meta, sampling_point == 1)
counts_point1 <-  dplyr::select(counts_raw, c(gene, meta_point1$id))

rownames(meta) <- meta$id
rownames(meta) %in% colnames(counts_raw)

se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts_point1[,-1])),
  colData = meta_point1,
  metadata = list(notes = "covid 19 (Omicron), only sampling point one is included in assay, full data in meta data",
                  full_counts = counts_raw,
                  full_meta = meta)
  
)

rownames(se) <- mapIds(EnsDb.Hsapiens.v86,
                       keys = counts_raw$gene,
                       keytype = "SYMBOL",
                       columns = "ENSEMBL")

se <- se[!is.na(rownames(se)),]
saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))




