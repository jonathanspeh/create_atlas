library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)


GEO_accs <- "GSE152418"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE152418&format=file&file=GSE152418%5Fp20047%5FStudy1%5FRawCounts%2Etxt%2Egz"

dir_path <- here::here("data", GEO_accs)
file_path <- here::here(dir_path, paste0(GEO_accs, ".tsv.gz"))



if(!file.exists(file_path)){
  dir.create(dir_path)
  curl::curl_download(download,
                      destfile = file_path)
  }


sm <- getGEO(GEO_accs, destdir = dir_path) 

glimpse(pData(sm[[1]]))
unique(pData(sm[[1]])$`characteristics_ch1.2`)

as.data.frame(pData(sm[[1]])) |>
  count(`condition:ch1`)

meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = `title`,
                sex = `gender:ch1`,
                source = `cell type:ch1`,
                timepoint = `days_post_symptom_onset:ch1`,
                covid_severity = `severity:ch1`
                ) |>
  mutate(disease = case_when(grepl("Conval", characteristics_ch1.2) ~ "covid-convalescent",
                             grepl("COVID", characteristics_ch1.2) ~ "covid-19",
                             grepl("Healthy", characteristics_ch1.2) ~ "healthy"),
         age = NA,
         dataset = GEO_accs, 
         pediatric = FALSE, 
         source = "PBMCs") |>

  dplyr::select(id, sample_name, pediatric,
                disease, timepoint, 
                processing_info, source, dataset, age, sex) 

counts <- fread(file_path, data.table = FALSE) |> rename("gene" = "ENSEMBLID")


ids <- meta$sample_name
names(ids) <- meta$id

counts <- counts |> 
  dplyr::select(gene, all_of(meta$sample_name)) |>
  rename(all_of(ids))

all(colnames(counts)[-1] == meta$id)


meta_selected <- meta |> dplyr::filter(disease != "covid-convalescent")
counts_selected <- counts[, c("gene", meta_selected$id)]


se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts_selected[,-1])),
  colData = meta_selected,
  metadata = list(notes = "Patients with covid-19 and healthy control, full data (including covalescent) in metadata",
                  full_counts = counts, 
                  full_meta = meta))

all(floor(assay(se)) == assay(se), na.rm = TRUE)

rownames(se) <- counts$gene

se <- se[!is.na(rownames(se)),]
saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))
