library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(gplots)
library(pheatmap)
library(dendextend)


# Prepare #####
combined_se <- readRDS(here::here("data", "ses", "combined_atlas.RDS"))
adult_se <- combined_se[,!colData(combined_se)$pediatric]
ped_se <- combined_se[,colData(combined_se)$pediatric]


# all
selected <- as.data.frame(colData(combined_se)) |>
  count(disease) |>
  dplyr::filter(n>=10) 
combined_min_10 <- combined_se[,colData(combined_se)$disease %in% selected$disease]

# adults only
selected <- as.data.frame(colData(adult_se)) |>
  count(disease) |>
  dplyr::filter(n>=10)
adult_min_10 <- adult_se[,colData(adult_se)$disease %in% selected$disease]

# kids only
selected <- as.data.frame(colData(ped_se)) |>
  count(disease) |>
  dplyr::filter(n>=10) 
ped_min_10 <- ped_se[,colData(ped_se)$disease %in% selected$disease]


# heatmaps #####
## heatmaps all #####
# All diseases with at least 10 individuals
meta <- as.data.frame(colData(combined_min_10))

sample_cor <- cor(assay(combined_min_10, "counts"))
sample_dist <- as.dist(1-sample_cor)
sample_tree <- hclust(sample_dist)
sample_dend <- as.dendrogram(sample_tree)

clusters <- cutree(sample_dend, k = 10)
clusters.df <- data.frame(sample = names(clusters), cluster = clusters)


color.scheme <- rev(RColorBrewer::brewer.pal(10,"RdBu"))
col_annot <- meta |>
  left_join(clusters.df, by = join_by("id" == "sample")) |>
  dplyr::select(dataset, cluster) |>
  mutate(cluster = factor(cluster))

rownames(col_annot) <- rownames(meta)

row_annot <- meta |>
  dplyr::select(disease) 


pheatmap::pheatmap(sample_cor,
                   color = color.scheme,
                   cluster_rows = sample_tree,
                   cluster_cols = sample_tree,
                   annotation_col = col_annot, 
                   annotation_row = row_annot,
                   #annotation_legend = FALSE,
                   legend = TRUE,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   width = 20,
                   height = 15,
                   filename = here::here("results", "heatmap_all_min10.png")
                   )

## heatmaps adults #####

meta <- as.data.frame(colData(adult_min_10))
sample_cor <- cor(assay(adult_min_10, "counts"))

sample_dist <- as.dist(1-sample_cor)
sample_tree <- hclust(sample_dist)
sample_dend <- as.dendrogram(sample_tree)


clusters <- cutree(sample_dend, k = 5)
clusters.df <- data.frame(sample = names(clusters), cluster = clusters)


color.scheme <- rev(RColorBrewer::brewer.pal(10,"RdBu"))
col_annot <- meta |>
  left_join(clusters.df, by = join_by("id" == "sample")) |>
  dplyr::select(dataset, cluster) |>
  mutate(cluster = factor(cluster))

rownames(col_annot) <- rownames(meta)

row_annot <- meta |>
  dplyr::select(disease) 


pheatmap::pheatmap(sample_cor,
                   color = color.scheme,
                   cluster_rows = sample_tree,
                   cluster_cols = sample_tree,
                   annotation_col = col_annot, 
                   annotation_row = row_annot,
                   #annotation_legend = FALSE,
                   legend = TRUE,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   width = 20,
                   height = 15,
                   filename = here::here("results", "heatmap_adult_min10.png")
)


## heatmap kids #####
meta <- as.data.frame(colData(ped_min_10))
sample_cor <- cor(assay(ped_min_10, "counts"))

sample_dist <- as.dist(1-sample_cor)
sample_tree <- hclust(sample_dist)
sample_dend <- as.dendrogram(sample_tree)


clusters <- cutree(sample_dend, k = 5)
clusters.df <- data.frame(sample = names(clusters), cluster = clusters)


color.scheme <- rev(RColorBrewer::brewer.pal(10,"RdBu"))
col_annot <- meta |>
  left_join(clusters.df, by = join_by("id" == "sample")) |>
  dplyr::select(dataset, cluster) |>
  mutate(cluster = factor(cluster))

rownames(col_annot) <- rownames(meta)

row_annot <- meta |>
  dplyr::select(disease) 


pheatmap::pheatmap(sample_cor,
                   color = color.scheme,
                   cluster_rows = sample_tree,
                   cluster_cols = sample_tree,
                   annotation_col = col_annot, 
                   annotation_row = row_annot,
                   #annotation_legend = FALSE,
                   legend = TRUE,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   width = 20,
                   height = 15,
                   filename = here::here("results", "heatmap_ped_min10.png")
)







# PCA #####
## all #####
use_se <- combined_se
pca <- prcomp(t(assay(use_se , "counts")))

pcs_joined <- as_tibble(pca$x, rownames = "id") |>
  left_join(as.data.frame(colData(use_se)))


p1 <- pcs_joined |>
  ggplot(aes(PC1, PC2, colour = disease, 
             #shape = disease
  )) +
  geom_point(size = 1.5) +
  ggtitle("Full data")


xdens <-
  cowplot::axis_canvas(p1, axis = "x") +
  geom_density(data = pcs_joined, 
               aes(PC1, fill = dataset, colour = dataset), 
               alpha = .3)

ydens <- 
  cowplot::axis_canvas(p1, axis = "y", coord_flip = TRUE) +
  geom_density(data = pcs_joined, 
               aes(PC2, fill = dataset, colour = dataset), 
               alpha = .3) +
  coord_flip()

p1 |>
  cowplot::insert_xaxis_grob(xdens, 
                             grid::unit(1, "cm"), 
                             position = "top") |>
  cowplot::insert_yaxis_grob(ydens, 
                             grid::unit(2, "cm"), 
                             position = "right") |>
  
  cowplot::ggdraw()



## adults #####
use_se <- adult_se
pca <- prcomp(t(assay(use_se , "counts")))
pcs_joined <- as_tibble(pca$x, rownames = "id") |>
  left_join(as.data.frame(colData(use_se)))


p1 <- as_tibble(pca$x, rownames = "id") |>
  left_join(as.data.frame(colData(use_se))) |>
  ggplot(aes(PC1, PC2, colour = disease, 
             #shape = disease
  )) +
  geom_point(size = 1.5) +
  ggtitle("Only adults")


as_tibble(pca$x, rownames = "id") |>
  left_join(as.data.frame(colData(use_se))) |>
  ggplot(aes(PC1, fill = dataset, colour = dataset)) +
  geom_density(alpha = .3)

as_tibble(pca$x, rownames = "id") |>
  left_join(as.data.frame(colData(use_se))) |>
  ggplot(aes(PC4, fill = disease, colour = disease)) +
  geom_density(alpha = .3)

as_tibble(pca$x, rownames = "id") |>
  left_join(as.data.frame(colData(use_se))) |>
  ggplot(aes(PC1, fill = sex, colour = sex)) +
  geom_density(alpha = .3)

as_tibble(pca$x, rownames = "id") |>
  left_join(as.data.frame(colData(use_se))) |>
  ggplot(aes(PC1, colour = dataset, fill = sex)) +
  geom_density(alpha = .3)


xdens <-
  cowplot::axis_canvas(p1, axis = "x") +
  geom_density(data = pcs_joined, 
               aes(PC1, fill = dataset, colour = dataset), 
               alpha = .3)

ydens <- 
  cowplot::axis_canvas(p1, axis = "y", coord_flip = TRUE) +
  geom_density(data = pcs_joined, 
               aes(PC2, fill = dataset, colour = dataset), 
               alpha = .3) +
  coord_flip()

p1 |>
  cowplot::insert_xaxis_grob(xdens, 
                             grid::unit(1, "cm"), 
                             position = "top") |>
  cowplot::insert_yaxis_grob(ydens, 
                             grid::unit(2, "cm"), 
                             position = "right") |>
  
  cowplot::ggdraw()




## kids #####
use_se <- ped_se
pca <- prcomp(t(assay(use_se , "counts")))

as_tibble(pca$x, rownames = "id") |>
  left_join(as.data.frame(colData(use_se))) |>
  ggplot(aes(PC1, PC2, colour = dataset, 
             #shape = dataset
  )) +
  geom_point(size = 1.5) +
  ggtitle("Only pediatrics")


use_se <- ped_min_10
pca <- prcomp(t(assay(use_se , "counts")))

pcs_joined <- as_tibble(pca$x, rownames = "id") |>
  left_join(as.data.frame(colData(use_se)))


p1 <- pcs_joined |>
  ggplot(aes(PC1, PC2, colour = disease, 
             #shape = disease
  )) +
  geom_point(size = 1.5) +
  ggtitle("Only kids")


xdens <-
  cowplot::axis_canvas(p1, axis = "x") +
  geom_density(data = pcs_joined, 
               aes(PC1, fill = dataset, colour = dataset), 
               alpha = .3)

ydens <- 
  cowplot::axis_canvas(p1, axis = "y", coord_flip = TRUE) +
  geom_density(data = pcs_joined, 
               aes(PC2, fill = dataset, colour = dataset), 
               alpha = .3) +
  coord_flip()

p1 |>
  cowplot::insert_xaxis_grob(xdens, 
                             grid::unit(1, "cm"), 
                             position = "top") |>
  cowplot::insert_yaxis_grob(ydens, 
                             grid::unit(2, "cm"), 
                             position = "right") |>
  
  cowplot::ggdraw()

