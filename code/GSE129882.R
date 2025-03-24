library(GEOquery)
library(data.table)
library(SummarizedExperiment)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)

#. Failed -only transcipt level annotations in GEO object -  could use rawer data


GEO_accs <- "GSE129882"
download <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE129882&format=file&file=GSE129882%5FCountsMatrix%2Etranscript%2Ecsv%2Egz"
file <- paste0(GEO_accs, ".csv.gz")
dir_path <- paste0(here::here("data", GEO_accs), "/")

dir.create(dir_path)
curl::curl_download(download,
                    destfile = here::here(dir_path, file))


sm <- getGEO(GEO_accs, destdir = dir_path)

meta <- pData(sm[[1]]) |> 
  rowwise() |>
  mutate(processing_info = list(across(everything()))) |>
  ungroup() |>
  dplyr::rename(id = geo_accession,
                sample_name = title,
                sample_type =  `characteristics_ch1.3`,
                group = "exposure:ch1",
                age = `age:ch1`,
                timepoint = `time:ch1`) |> 
  mutate(individual = stringr::str_remove(characteristics_ch1, "patient: "),
         disease = "denguevirus infection",
         dataset = GEO_accs,
         pediatric = age < 18,
         sex = stringr::str_remove("characteristics_ch1.4", "Sex: "),
         source = "PBMC"
         ) |>
  dplyr::select(id, individual, sample_name, sample_type, age, pediatric, sex, group, 
                disease, processing_info, source, dataset, timepoint)


counts_raw <- fread(here::here(dir_path, file)) |> rename("Genes" = "V1")


meta_active <- dplyr::filter(meta, grepl("Early", timepoint))
counts_active <-  dplyr::select(counts_raw, c(Genes, meta_active$sample_name))


length(counts_active$Genes)
length(unique(counts_active$Genes))

length(rownames(se))
length(unique(rownames(se)))

filter(EnsDb.Hsapiens.v86,
       filter = ~ gene_biotype == "protein_coding")


counts_active$Genes <- mapIds(db_coding,
       keys = counts_active$Genes,
       keytype = "TXID",
       columns = "ENSEMBL")


counts_active |>
  count(Genes) |>
  dplyr::filter(n > 1) |> 
  arrange(desc(n))

counts_active |>
  dplyr::filter(Genes == "ENSG00000131051")


se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts_active[,-1])),
  colData = meta_active,
  metadata = list(notes = "Active denguevirus infection, late active and convalescent stored in metadata",
                  full_counts = counts_raw,
                  full_meta = meta))

rownames(se) <- counts_active$Genes


tx2gene <- mapIds(EnsDb.Hsapiens.v86,
                  keys = counts_active$Genes,
                  keytype = "TXID",
                  columns = "ENSEMBL")
  
tx2gene_df <- data.frame(txids = names(tx2gene),
                         genids = tx2gene)
                  
  
tximport::tximport(files = here::here(dir_path, file), 
                   type = "stringtie",
                   txOut = FALSE,
                   tx2gene = tx2gene_df
                   )


tximport::summarizeToGene(list(assay(se)))

db_coding <- filter(EnsDb.Hsapiens.v86,
       filter = ~ gene_biotype == "protein_coding")

rownames(se) <- mapIds(db_coding,
       keys = counts_active$Genes,
       keytype = "TXID",
       columns = "ENSEMBL")


se <- se[!is.na(rownames(se)),]

saveRDS(se, paste0("data/ses/", GEO_accs, "_se.RDS"))


