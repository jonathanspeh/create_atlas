library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)


GEO_accs <- "GSE163151"
file <- paste0(GEO_accs, ".tar")
dir_path <- paste0(here::here("data", GEO_accs), "/")

#curl::curl_download(paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=", 
#                           GEO_accs, 
#                           "&format=file"),
#                    destfile = here::here("data", file))

untar(here::here("data", file), 
      exdir = dir_path)

sm <- getGEO(GEO_accs, destdir = dir_path)


meta <- pData(sm[[1]]) |>  
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = title,
                sample_type =  "disease state:ch1",
                disease = "pathogen:ch1",
                source = "source_name_ch1",
                severity = "disease severity:ch1"
  ) |> 
  mutate(individual = sub("^\\D+", "\\1", sample_name),
         dataset = GEO_accs, 
         age = NA, 
         sex = NA,
         pediatric = FALSE,
         #      processing_info = paste(
         # "growth_protocol:", growth_protocol_ch1,
         # "extract_protocol:", extract_protocol_ch1,
         # "library_prep:", extract_protocol_ch1.1, 
         # "data_processing_1:", data_processing,
         # "data_processing_2:", data_processing.1,
         # "assembly:", data_processing.2,
         # "instrument:", instrument_model,
         # sep = "\t")
  ) |>
  dplyr::select(id, individual, sample_name, sample_type,
                disease, severity,
                processing_info, source, dataset, age, pediatric, sex)



files <- list.files(dir_path, full.names = FALSE)
files <- files[grepl("GSM", files)]

read_counts <- function(name){
  sample <- stringr::str_remove(name, "_.*.txt.gz")
  file <- fread(paste0(dir_path, name))  #|>
    #dplyr::filter(!grepl("^__", V1)) 
  colnames(file) <- c("gene", substitute(sample))
  file
}

counts <- lapply(files, read_counts)

cl <- lapply(counts, function(x) x$gene)
common_genes <- Reduce(intersect, cl)
count_list <- lapply(counts, dplyr::filter, gene %in% common_genes)
counts_raw <- purrr::reduce(count_list, full_join, by = "gene") 


meta_wb <- dplyr::filter(meta, source == "Whole Blood")
counts_wb <-  dplyr::select(counts_raw, c("gene", meta_wb$id))

meta_wb <- meta_wb |>
  mutate(disease = 
           case_when(grepl("No path", disease) ~ "healthy",
                     grepl("CoV-2", disease) ~ "COVID-19",
                     grepl("Escherichia", disease) ~ "escherichia coli infection",
                     grepl("Group C Strep", disease) ~ "group c streptococcal infection",
                     grepl("Klebsiella pn", disease) ~ "klebsiella pneumoniae infection",
                     grepl("Klebsiella ox", disease) ~ "klebsiella oxytoca infection",
                     TRUE ~ disease))


se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts_wb[,-1])),
  colData = meta_wb,
  metadata = list(notes = "Contains only whole blood samples, full data in meta data",
                  full_counts = counts_raw,
                  full_meta = meta
                  )
)


rownames(se) <- mapIds(EnsDb.Hsapiens.v86,
                       keys = counts_raw$gene,
                       keytype = "SYMBOL",
                       columns = "ENSEMBL")

se <- se[!is.na(rownames(se)),]
saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))


