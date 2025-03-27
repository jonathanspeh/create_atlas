library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)



accs <- "E-MTAB-11671"
#dir.create(here::here("data", accs))
dir_path <- paste0(here::here("data", accs), "/")

# curl::curl_download("https://www.ebi.ac.uk/biostudies/files/E-MTAB-11671/RNA-Seq_genecounts.csv",
#                     destfile = here::here("data", accs, paste0(accs, "_counts.csv")))
# 
# curl::curl_download("https://www.ebi.ac.uk/biostudies/files/E-MTAB-11671/E-MTAB-11671.sdrf.txt",
#                     destfile = here::here("data", accs, paste0(accs, "_sdrf.txt")))

counts <- fread(paste0(dir_path, accs, "_counts.csv")) |> rename("gene" = "V1")
meta <- fread(paste0(dir_path, accs, "_sdrf.txt")) 

colnames(counts) <- paste0(accs, "_", janitor::make_clean_names(colnames(counts)))


colnames(meta) <- janitor::make_clean_names(colnames(meta))
meta$source_name <- janitor::make_clean_names(meta$source_name)


#TODO fix duplicated rows
meta_fixed <- meta[rep(c(T,F), 411),] |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(
    id = source_name,
    age_month = characteristics_age,
    sex = characteristics_sex,
    disease = characteristics_disease) |>
  mutate(
          sample_name = id,
          id = paste0(accs, "_", id),
          age = as.character(round(age_month/12)),
          # processing_info = paste(
          #   "comment_library_layout", comment_library_layout,
          #   "comment_library_selection", comment_library_selection,
          #   "comment_library_source", comment_library_source, 
          #   "comment_library_strand", comment_library_strand,
          #   "comment_library_strategy", comment_library_strategy,
          #   sep = "\t"),
          dataset = accs,
          source = "whole_blood",
          pediatric = TRUE
          ) |>
  dplyr::select(id, sample_name, age_month, age, sex, pediatric, disease, processing_info, source, dataset)

rownames(counts) <- counts$`E-MTAB-11671_gene`
rownames(meta_fixed) <- meta_fixed$id

counts <- counts |> dplyr::select(`E-MTAB-11671_gene`, meta_fixed$id)

all(rownames(meta_fixed) == colnames(counts)[-1])


se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts[,-1])),
  colData = meta_fixed,
  metadata = list(notes = "various diseases, peadiatric patients")
)    
  
rownames(se) <- rownames(counts)

all(colData(se)$id  %in% colnames(se))


se <- se[!is.na(rownames(se)),]
saveRDS(se, paste0("data/ses/", accs, "_ped_se.RDS"))



