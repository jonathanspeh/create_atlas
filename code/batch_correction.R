library(sva)
library(SummarizedExperiment)
library(tidyverse)



adult_se <- readRDS(here::here("data", "ses", "selected_adult_no_hg19.RDS"))
counts <- assay(adult_se)
meta <- as.data.frame(colData(adult_se))

batch <- as.numeric(factor(meta$dataset))
group <- as.numeric(factor(meta$disease))


adjusted_counts <- ComBat_seq(counts, 
                              batch=batch, 
                              group=group, 
                              full_mod=TRUE)

adjusted_counts_no_group <- ComBat_seq(counts, 
                              batch=batch, 
                              group=NULL, 
                              full_mod=FALSE)



pca_raw_obj = prcomp(t(counts))
pca_raw_join <- as_tibble(pca_raw_obj$x, rownames = "id") |>
  left_join(meta)

pca_adjusted_obj = prcomp(t(adjusted_counts))
pca_adjusted_join <- as_tibble(pca_adjusted_obj$x, rownames = "id") |>
  left_join(meta)

pca_adjusted_no_group_obj = prcomp(t(adjusted_counts_no_group))
pca_adjusted_no_group_join <- as_tibble(pca_adjusted_no_group_obj$x, rownames = "id") |>
  left_join(meta)



pca_raw_obj_healthy = prcomp(t(counts[, filter(meta, disease == "healthy")$id]))
pca_raw_join_healthy <- as_tibble(pca_raw_obj_healthy$x, rownames = "id") |>
  left_join(meta)

pca_adjusted_obj_healthy = prcomp(t(adjusted_counts[, filter(meta, disease == "healthy")$id]))
pca_adjusted_join_healthy <- as_tibble(pca_adjusted_obj_healthy$x, rownames = "id") |>
  left_join(meta)

pca_adjusted_no_group_obj_healthy = prcomp(t(adjusted_counts_no_group[, filter(meta, disease == "healthy")$id]))
pca_adjusted_no_group_join_healthy <- as_tibble(pca_adjusted_no_group_obj_healthy$x, rownames = "id") |>
  left_join(meta)




p1 <- pca_raw_join |>
  ggplot(aes(x = PC1, y = PC2, colour = dataset)) +
  geom_point() + 
  #stat_ellipse() +
  ggtitle("raw")

p2 <- pca_raw_join |>
  ggplot(aes(x = PC1, y = PC2, colour = disease)) +
  geom_point() + 
  #stat_ellipse() +
  ggtitle("raw")


p1_adj <- pca_adjusted_join |>
  ggplot(aes(x = -PC1, y = PC2, colour = dataset)) +
  geom_point() +  
  #stat_ellipse() +
  ggtitle("adjusted")

p2_adj <- pca_adjusted_join |>
  ggplot(aes(x = -PC1, y = PC2, colour = disease)) +
  geom_point() +
  #stat_ellipse() +
  ggtitle("adjusted")


p1_adj_no_group <- pca_adjusted_no_group_join |>
  ggplot(aes(x = -PC1, y = PC2, colour = dataset)) +
  geom_point() +  
  #stat_ellipse() +
  ggtitle("adjusted - no group variable")

p2_adj_no_group <- pca_adjusted_no_group_join |>
  ggplot(aes(x = -PC1, y = PC2, colour = disease)) +
  geom_point() +
  #stat_ellipse() +
  ggtitle("adjusted - no group variable")



cowplot::plot_grid(p1, p2, ncol =1)
cowplot::plot_grid(p1_adj, p2_adj, ncol =1)
cowplot::plot_grid(p1_adj_no_group, p2_adj_no_group, ncol =1)

cowplot::plot_grid(p1, p1_adj, p1_adj_no_group, ncol =1)
cowplot::plot_grid(p2, p2_adj, p2_adj_no_group, ncol =1)


plot_pc_density <- function(group, dataset, PC="PC1") {
  dataset |>
    ggplot(aes(x = !!sym(PC), fill = !!sym(group), color = !!sym(group))) +
    geom_density(alpha = .5) + 
    ggtitle(group) +
    theme(legend.title=element_blank())
}

plot_pc_density("dataset", filter(pca_raw_join, disease == "healthy"))
plot_pc_density("dataset", pca_raw_join_healthy) +
  ggtitle("Healthy controls - no adjustment")

plot_pc_density("dataset", filter(pca_adjusted_join, disease == "healthy"))
plot_pc_density("dataset", pca_adjusted_join_healthy)

plot_pc_density("dataset", filter(pca_adjusted_no_group_join, disease == "healthy"))
plot_pc_density("dataset", pca_adjusted_no_group_join_healthy) + 
  ggtitle("Healthy controls - combat-seq")



plots_raw <- lapply(c("disease", "dataset"), 
                plot_pc_density, 
                dataset = pca_raw_join, 
                PC = "PC1")

cowplot::plot_grid(plotlist = plots_raw, ncol = 1) 


plots_adjusted <- lapply(c("disease", "dataset"), 
                    plot_pc_density, 
                    dataset = pca_adjusted_join, 
                    PC = "PC1")
cowplot::plot_grid(plotlist = plots_adjusted, ncol = 1) 


plots_adjusted_no_group <- lapply(c("disease", "dataset"), 
                         plot_pc_density, 
                         dataset = pca_adjusted_no_group_join, 
                         PC = "PC1")
cowplot::plot_grid(plotlist = plots_adjusted_no_group, ncol = 1) 



plots_raw_2 <- lapply(c("PC1", "PC2", "PC3", "PC4"), 
                    plot_pc_density, 
                    group = "dataset",
                    dataset = pca_raw_join)
cowplot::plot_grid(plotlist = plots_raw_2, ncol = 1) 


plots_adjusted_2 <- lapply(c("PC1", "PC2", "PC3", "PC4"), 
                         plot_pc_density, 
                         group = "dataset",
                         dataset = pca_adjusted_join)
cowplot::plot_grid(plotlist = plots_adjusted_2, ncol = 1) 


plots_adjusted_no_group_2 <- lapply(c("PC1", "PC2", "PC3", "PC4"), 
                           plot_pc_density, 
                           group = "dataset",
                           dataset = pca_adjusted_no_group_join)

cowplot::plot_grid(plotlist = plots_adjusted_no_group_2, ncol = 1) 



raw_cor <- cor(counts)

raw_dist <- as.dist(1-raw_cor)
raw_tree <- hclust(raw_dist)
raw_dend <- as.dendrogram(raw_tree)


clusters_raw <- cutree(raw_dend, k = 5)
clusters.df_raw <- data.frame(sample = names(clusters_raw), cluster = clusters_raw)

col_annot_raw <- meta |>
  left_join(clusters.df_raw, by = join_by("id" == "sample")) |>
  dplyr::select(dataset, cluster) |>
  mutate(cluster = factor(cluster))

rownames(col_annot_raw) <- rownames(meta)

row_annot <- meta |>
  dplyr::select(disease) 

color.scheme <- rev(RColorBrewer::brewer.pal(10,"RdBu"))
pheatmap::pheatmap(raw_cor,
                   color = color.scheme,
                   cluster_rows = raw_tree,
                   cluster_cols = raw_tree,
                   annotation_col = col_annot_raw, 
                   annotation_row = row_annot,
                   #annotation_legend = FALSE,
                   legend = FALSE,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
)





adjusted_cor <- cor(adjusted_counts)

adjusted_dist <- as.dist(1-adjusted_cor)
adjusted_tree <- hclust(adjusted_dist)
adjusted_dend <- as.dendrogram(adjusted_tree)


clusters_adjusted <- cutree(adjusted_dend, k = 5)
clusters.df_adjusted <- data.frame(sample = names(clusters_adjusted), cluster = clusters_adjusted)

col_annot_adjusted <- meta |>
  left_join(clusters.df_adjusted, by = join_by("id" == "sample")) |>
  dplyr::select(dataset, cluster) |>
  mutate(cluster = factor(cluster))

rownames(col_annot_adjusted) <- rownames(meta)

row_annot <- meta |>
  dplyr::select(disease) 

color.scheme <- rev(RColorBrewer::brewer.pal(10,"RdBu"))
pheatmap::pheatmap(adjusted_cor,
                   color = color.scheme,
                   cluster_rows = adjusted_tree,
                   cluster_cols = adjusted_tree,
                   annotation_col = col_annot_adjusted, 
                   annotation_row = row_annot,
                   #annotation_legend = FALSE,
                   legend = FALSE,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
)





adjusted_no_group_cor <- cor(adjusted_counts_no_group)

adjusted_no_group_dist <- as.dist(1-adjusted_no_group_cor)
adjusted_no_group_tree <- hclust(adjusted_no_group_dist)
adjusted_no_group_dend <- as.dendrogram(adjusted_no_group_tree)


clusters_adjusted_no_group <- cutree(adjusted_no_group_dend, k = 5)
clusters.df_adjusted_no_group <- data.frame(sample = names(clusters_adjusted_no_group),
                                            cluster = clusters_adjusted_no_group)

col_annot_adjusted_no_group <- meta |>
  left_join(clusters.df_adjusted_no_group, 
            by = join_by("id" == "sample")) |>
  dplyr::select(dataset, cluster) |>
  mutate(cluster = factor(cluster))

rownames(col_annot_adjusted_no_group) <- rownames(meta)

row_annot <- meta |>
  dplyr::select(disease) 

color.scheme <- rev(RColorBrewer::brewer.pal(10,"RdBu"))
pheatmap::pheatmap(adjusted_no_group_cor,
                   color = color.scheme,
                   cluster_rows = adjusted_no_group_tree,
                   cluster_cols = adjusted_no_group_tree,
                   annotation_col = col_annot_adjusted_no_group, 
                   annotation_row = row_annot,
                   #annotation_legend = FALSE,
                   legend = FALSE,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
)



batch
mod0 <- model.matrix(~batch, data = meta) # could also add platform as covariate
mod1 <- model.matrix(~disease + batch + platform_id, data = meta)

counts_sv <- svaseq(counts, 
                    mod = mod1, 
                    mod0 = mod0,
                    method = "irw"
                    )


adju_sv <- svaseq(adjusted_counts, 
                    mod = mod1, 
                    mod0 = mod0,
                    method = "irw"
)

plot(counts_sv$sv[,1])

library(factoextra)

fviz_eig(pca_raw_obj)
fviz_eig(pca_adjusted_obj)
fviz_eig(pca_adjusted_no_group_obj)


raw_var <- get_pca_var(pca_raw_obj)
raw_var$coord
raw_var$contrib[1:10, 1:10]
raw_var$contrib |>
  as_tibble(rownames = "gene") |>
  filter(Dim.1 > 0.1) |>
  ggplot(aes(x = reorder(gene, Dim.1), 
             y = Dim.1
             )) +
  geom_point()

    raw_var$cos2   
    
fviz_contrib(pca_raw_obj, choice = "var", axes = 1, top = 10)

#TODO - GSEA on PC1 contributions - it's kind of an ordered list

#Results for individuals
raw_ind <- get_pca_ind(pca_raw_obj)
raw_ind$coord          
raw_ind$contrib        
raw_ind$cos2

