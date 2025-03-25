library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(gplots)
library(dendextend)
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
  #ggnewscale::new_scale_fill() + 
  #geom_tile(aes(y = -11, fill = disease, height = 5)) +
  ggnewscale::new_scale_fill() + 
  geom_tile(aes(y = -17, fill = source, height = 5))


# hirarchical clustering - adults
sample_cor <- cor(assay(adult_se, "counts"))
sample_dist <- as.dist(1-sample_cor)
sample_tree <- hclust(sample_dist)
sample_dend <- as.dendrogram(sample_tree)

clusters <- cutree(sample_dend, k = 10)
clusters.df <- data.frame(sample = names(clusters), cluster = clusters)


color.scheme <- rev(RColorBrewer::brewer.pal(10,"RdBu"))
col_annot <- colData(adult_se) |>
  as.data.frame() |>
  left_join(clusters.df, by = join_by("id" == "sample")) |>
  select(dataset, cluster) |>
  mutate(cluster = factor(cluster))

rownames(col_annot) <- rownames(colData(adult_se))

row_annot <- colData(adult_se) |>
  as.data.frame() |>
  select(sex) 
  
  
pheatmap::pheatmap(sample_cor,
                   color = color.scheme,
                   cluster_cols = sample_tree,
                   cluster_rows = sample_tree,
                   annotation_col = col_annot, 
                   annotation_row = row_annot,
                   legend = TRUE,
                   show_rownames = FALSE,
                   show_colnames = FALSE
                   )

# hirarchical clustering - kids
sample_cor <- cor(assay(ped_se, "counts"))
sample_dist <- as.dist(1-sample_cor)
sample_tree <- hclust(sample_dist)
sample_dend <- as.dendrogram(sample_tree)

clusters <- cutree(sample_dend, k = 10)
clusters.df <- data.frame(sample = names(clusters), cluster = clusters)


color.scheme <- rev(RColorBrewer::brewer.pal(10,"RdBu"))
col_annot <- colData(ped_se) |>
  as.data.frame() |>
  left_join(clusters.df, by = join_by("id" == "sample")) |>
  select(dataset, cluster) |>
  mutate(cluster = factor(cluster))

rownames(col_annot) <- rownames(colData(ped_se))

row_annot <- colData(ped_se) |>
  as.data.frame() |>
  select(sex, disease) 


pheatmap::pheatmap(sample_cor,
                   color = color.scheme,
                   cluster_cols = sample_tree,
                   cluster_rows = sample_tree,
                   annotation_col = col_annot, 
                   annotation_row = row_annot,
                   legend = TRUE,
                   show_rownames = FALSE,
                   show_colnames = FALSE
)

colData(adult_se) |>
  as.data.frame() |>
  #mutate(sex = dplyr::case_when(sex == "f" ~ "female",
  #                       sex == "m" ~ "male",
  #                       TRUE ~ sex)) |>
  count(sex)




pca <- prcomp(t(assay(adult_se , "counts")))

as_tibble(pca$x, rownames = "id") |>
  left_join(as.data.frame(colData(combined_se))) |>
  ggplot(aes(PC1, PC2, colour = disease, 
             #shape = disease
             )) +
  geom_point(size = 1.5)



colData(combined_se) |>
  as_tibble() |>
  group_by(disease, dataset) |>
  tally() |>
  group_by(disease) |>
  summarise(n = sum(n),
            n_datasets = n()) |>
  arrange(desc(n)) |>
  print(n = 50) 

  
colData(adult_se) |>
  as_tibble() |>
  count(disease) |>
  arrange(desc(n))


colData(combined_se) |>
  as_tibble() |>
  count(disease) |>
  arrange(desc(n))

