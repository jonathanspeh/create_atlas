library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)


GEO_accs <- "GSE131411"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE131411&format=file&file=GSE131411%5Frawcounts%5FCS%5FSS%2Exlsx"
download_ctrl <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE131411&format=file&file=GSE131411%5FControls%2Exlsx"

dir_path <- here::here("data", GEO_accs)
file_path <- here::here(dir_path, paste0(GEO_accs, ".xlsx"))
ctrl_path <- here::here(dir_path, paste0(GEO_accs, "_ctrl.xlsx"))



if(!file.exists(file_path)){
  dir.create(dir_path)
  curl::curl_download(download,
                      destfile = file_path)
  curl::curl_download(download_ctrl,
                      destfile = ctrl_path)
  
  }


sm <- getGEO(GEO_accs, destdir = dir_path) 

meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = `title`,
                timepoint = `timepoint:ch1`,
                patient = `patient:ch1`,
                source = `tissue:ch1`) |>
  mutate(disease = ifelse(`clinical classification:ch1` == "SS", 
                          "septic_shock",
                          "cardiogenic_shock"),
         reanalysis = grepl("re-analysis", sample_name),
         age = NA,
         sex = NA,
         dataset = GEO_accs, 
         pediatric = FALSE, 
         sample_name = case_when(grepl("Control", sample_name) ~  stringr::str_remove_all(stringr::str_remove(sample_name, "Control"), "\\s"),
                                 TRUE ~ stringr::str_remove(sample_name, "_.*"))) |>

  dplyr::select(id, sample_name, pediatric,
                disease, 
                processing_info, source, dataset, age, sex,, reanalysis, timepoint) 


meta_2 <- pData(sm[[2]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = `title`,
                timepoint = `timepoint:ch1`,
                patient = `patient:ch1`,
                source = `tissue:ch1`) |>
  mutate(disease = case_when(`clinical classification:ch1` == "SS" ~ "septic_shock",
                             `clinical classification:ch1` == "CS" ~ "cardiogenic_shock",
                             TRUE ~ "healthy"
                             ),
         reanalysis = grepl("re-analysis", sample_name),
         age = NA,
         sex = NA,
         dataset = GEO_accs, 
         pediatric = FALSE,
         sample_name = case_when(grepl("Control", sample_name) ~  stringr::str_remove_all(stringr::str_remove(sample_name, "Control"), "\\s"),
                                 TRUE ~ stringr::str_remove(sample_name, "_.*")
                                 )) |>
  dplyr::select(id, sample_name, pediatric,
                disease, 
                processing_info, source, dataset, age, sex, reanalysis, timepoint) 


meta <- rbind(meta, meta_2)


counts <- readxl::read_xlsx(file_path) |> rename("gene" = "...1")
ctrl <- readxl::read_xlsx(ctrl_path) |> rename("gene" = "Geneid")
counts <- left_join(counts, ctrl, by = join_by(gene))
colnames(counts) <- stringr::str_remove(colnames(counts), "B")

ids <- meta$sample_name
names(ids) <- meta$id


counts <- counts |> 
  dplyr::select(gene, any_of(meta$sample_name)) |>
  rename(any_of(ids))

all(colnames(counts)[-1] == meta$id)

meta_t1 <- meta |>
  filter(timepoint == "T1",
         !reanalysis) #some samples are  re-analysis from GSE110487


counts_t1 <- counts |>
  dplyr::select(gene, all_of(meta_t1$id))
  
rownames(meta_t1) <- meta_t1$id

all(rownames(meta_t1) == colnames(counts_t1)[-1])

se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts_t1[,-1])),
  colData = meta_t1,
  metadata = list(notes = "Patients with Septic and cardiogenic shock, only Timepoint 1 in assay, some samples have been excluded as they are re-analyses from GSE110487",
                  full_meta = counts,
                  full_meta = meta))

all(floor(assay(se)) == assay(se), na.rm = TRUE)

rownames(se) <- counts_t1$gene

se <- se[!is.na(rownames(se)),]
saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))
