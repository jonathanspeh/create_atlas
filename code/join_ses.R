library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
source("code/functions.R")


files <- list.files(here::here("data", "ses"), full.names = TRUE)
files <- files[grepl("_se.RDS", files)]

counts_list <- lapply(files, readRDS)

combined_se <- combine_se(counts_list)

adult_se <- combined_se[,!colData(combined_se)$pediatric]
ped_se <- combined_se[,colData(combined_se)$pediatric]

sample_cor <- cor(assay(adult_se, "counts"))

as_tibble(sample_cor, rownames = "x") |>
  tidyr::pivot_longer(!x, names_to = "y") |>
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
  geom_tile(aes(y = -11, fill = disease, height = 5)) +
  ggnewscale::new_scale_fill() + 
  geom_tile(aes(y = -17, fill = source, height = 5))





pca <- prcomp(t(assay(adult_se , "counts")))

as_tibble(pca$x, rownames = "id") |>
  left_join(as.data.frame(colData(combined_se))) |>
  ggplot(aes(PC1, PC2, colour = disease, shape = dataset)) +
  geom_point(size = 1.5)
