
# Read count txt files and adjsut column names to contain sample accession
read_counts <- function(name){
  sample <- stringr::str_remove(name, "_.*.txt.gz")
  file <- data.table::fread(paste0(dir_path, name)) |>
    dplyr::filter(!grepl("^__", V1)) 
  colnames(file) <- c("gene", substitute(sample))
  file
}



# Takes list of summarised experiments and combines them
# TODO - add functions and look-up table to harmonise diseases


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
  
  ### Combine Assays - c
  #TODO - make sure datasets are joined on common genes
  assay_list <- lapply(se_list, assay, "counts")
  common_genes <- Reduce(intersect, lapply(assay_list, rownames))
  assay_list <- lapply(assay_list, dplyr::as_tibble, rownames="gene")
  if(common_only) {
    message(paste("joining on", length(common_genes), "common genes"))
    assay_list <- lapply(assay_list, dplyr::filter, gene %in% common_genes)
  }
  
  combined_assay <- purrr::reduce(assay_list, dplyr::full_join, by = "gene")
  
  
  
  ### Combine sample informations 
  coldata_list <- lapply(se_list, colData)
  coldata_list <- lapply(coldata_list, dplyr::as_tibble)  
  
  ## could use this for cleaner metadata - only use columns that are in all 
  #common_cols <- Reduce(intersect, lapply(coldata_list, colnames))
  #combined_meta <- dplyr::bind_rows(lapply(coldata_list, dplyr::select, dplyr::all_of(common))) 
  combined_meta <- dplyr::bind_rows(coldata_list) 
  combined_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list("counts" = combined_assay[,-1]),
    colData = combined_meta)
  rownames(combined_se) <- combined_assay$gene
  colnames(combined_se) <- combined_meta$id
  combined_se
}




