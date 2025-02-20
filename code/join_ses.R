library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
source("code/functions.R")


files <- list.files(here::here("data", "ses"), full.names = TRUE)
files <- files[grepl("_se.RDS", files)]

counts_list <- lapply(files, readRDS)


assay_list <- lapply(counts_list, assay, "counts")
common_genes <- Reduce(intersect, lapply(assay_list, rownames))
assay_list <- lapply(assay_list, dplyr::as_tibble, rownames="gene")

assays_filtered <- lapply(assay_list, dplyr::filter, gene %in% common_genes)

combined_assay <- purrr::reduce(assays_filtered, dplyr::full_join, by = "gene")



#TODO - add messages to combine_se()?
combined_se <- combine_se(counts_list)

norm_assays <- assay(combined_se, "counts")
norm_factors <- colSums(norm_assays)

norm_assays <- t(t(norm_assays) / norm_factors) * 1e6

sample_cor <- cor(assay(combined_se, "counts"))
norm_cor = cor(norm_assays)


as_tibble(norm_cor, rownames = "x") |>
  pivot_longer(!x, names_to = "y") |>
  left_join(as.data.frame(colData(combined_se)), by = join_by("x" == "id")) |>
  #dplyr::filter(value > .99) |>
  ggplot(aes(x = x, y = y )) +
  geom_tile(aes(fill = value)) +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  scale_y_discrete(labels = NULL, breaks = NULL) +
  scale_fill_gradient2(high = "red", mid = "white", low = "blue", midpoint = .6) + 
  ggnewscale::new_scale_fill() + 
  geom_tile(aes(y = -5, fill = dataset, height = 5)) + 
  scale_fill_brewer(palette = "Set1") + 
  ggnewscale::new_scale_fill() + 
  geom_tile(aes(y = -10, fill = disease, height = 5))





pca <- prcomp(t(assay(combined_se, "counts")))

as_tibble(pca$x, rownames = "id") |>
  left_join(as.data.frame(colData(combined_se))) |>
  ggplot(aes(PC1, PC2, colour = disease, shape = dataset)) +
  geom_point(size = 1.5)
