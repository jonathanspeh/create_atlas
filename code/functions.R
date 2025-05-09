
# Read count txt files and adjsut column names to contain sample accession
read_counts <- function(name){
  sample <- stringr::str_remove(name, "_.*.txt.gz")
  file <- data.table::fread(here::here(dir_path, name)) |>
    dplyr::filter(!grepl("^__", V1)) 
  colnames(file) <- c("gene", substitute(sample))
  file
}



# Verify that diseases are found in manually curated list of allowed inputs
verify_disease <- function(meta){
  allowed <- readxl::read_excel(here::here("data", "curated_diseases.xlsx"))$disease
  not_verified <- unique(meta[!meta$disease %in% allowed, c("disease", "dataset")])
  if(nrow(not_verified) != 0){
    warning("Not all diseases could be verified: ", paste(capture.output(print(not_verified)), 
                                                          collapse = "\n"))
  }
}


# Takes list of summarised experiments and combines them
combine_se <- function(se_list, common_only = TRUE){
  stopifnot(
    "alle SEs must have unique colnames" = !any(duplicated(unlist(lapply(se_list, colnames)))),
    "all SEs must have a counts assay" = all(unlist(lapply(se_list, assayNames)) %in% "counts"),
    "rownames of the SE cannot be NA" = all(!is.na(unlist(lapply(se_list, rownames)))),
    "all dataset need id, disease, dataset, age, sex, source and processing_info" = {
      coldata_list <- lapply(se_list, colData)
      all(
        c("id", "disease", "dataset", 
          "id", "age", "sex", "source", "processing_info") %in% 
          Reduce(intersect, lapply(coldata_list, colnames)))})
  
  ### Combine sample informations 
  coldata_list <- lapply(se_list, colData)
  coldata_list <- lapply(coldata_list, as.data.frame)  
  coldata_list <- lapply(coldata_list, function(x) {x$age = as.integer(x$age); x})
  
  ## could use this for cleaner metadata - only use columns that are in all 
  #common_cols <- Reduce(intersect, lapply(coldata_list, colnames))
  #combined_meta <- dplyr::bind_rows(lapply(coldata_list, dplyr::select, dplyr::all_of(common))) 
  combined_meta <- dplyr::bind_rows(coldata_list)
  
  #TODO - Make more elegant
  combined_meta$disease <- tolower(combined_meta$disease)
  combined_meta$disease <- stringr::str_replace_all(combined_meta$disease, " ", "_")
  combined_meta$source <- tolower(combined_meta$source)
  combined_meta$source <- stringr::str_replace_all(combined_meta$source, " ", "_")
  combined_meta$sex <- tolower(combined_meta$sex)
  combined_meta <- combined_meta |>
    dplyr::mutate(sex = dplyr::case_when(sex == "f" ~ "female",
                                         sex == "m" ~ "male",
                                         TRUE ~ sex))
  
  combined_meta$pediatric <- as.logical(combined_meta$pediatric)
  stopifnot("pediatric column cannot contain NA" = !is.na(combined_meta$pediatric))
  
  
  verify_disease(combined_meta)
  
  
  ### Combine Assays 
  assay_list <- lapply(se_list, assay, "counts")
  common_genes <- Reduce(intersect, lapply(assay_list, rownames))
  assay_list <- lapply(assay_list, dplyr::as_tibble, rownames="gene")
  if(common_only) {
    message(paste("joining on", length(common_genes), "common genes"))
    assay_list <- lapply(assay_list, dplyr::filter, gene %in% common_genes)
  }
  
  combined_assay <- purrr::reduce(assay_list, dplyr::full_join, by = "gene") 
  assay_mat <- as.matrix(combined_assay[,-1])
  rownames(assay_mat) <- combined_assay$gene
  
  stopifnot(
    "Rownames of metadata and colnames of assay must be equal" = all(rownames(combined_meta) == colnames(assay_mat)),
    "All counts must be integer" = all(round(combined_assay[,-1]) == combined_assay[,-1], na.rm = TRUE))
    
  combined_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list("counts" = assay_mat),
    colData = combined_meta)
 
  combined_se
}




