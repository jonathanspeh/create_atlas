library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
#library(gplots)
library(pheatmap)
library(dendextend)
library(sva)

combined_se <- readRDS(here::here("data", "ses", "combined_atlas.RDS"))
disease_info <- readxl::read_xlsx(here::here("data", "curated_diseases.xlsx"))

colData(combined_se) <- colData(combined_se) |>
  as.data.frame() |>
  left_join(disease_info, by = join_by(disease)) |>
  DataFrame()

colnames(combined_se) <- colData(combined_se)$id

adult_se <- combined_se[,!colData(combined_se)$pediatric]
ped_se <- combined_se[,as.logical(colData(combined_se)$pediatric)]

colData(combined_se) |>
  as_tibble() |>
  group_by(disease, pathogen_curated, dataset) |>
  tally() |>
  group_by(disease, pathogen_curated) |>
  summarise(n = sum(n),
            n_datasets = n()) |>
  arrange(desc(n)) |>
  dplyr::filter(n >= 10,
                n_datasets >= 2,
                pathogen_curated != "unknown_pathogen") |>
  flextable::flextable() |>
  flextable::set_caption("all samples") |>
  flextable::add_footer_lines("All diseases with at least 10 samples and 2 or more datasets where the pathogen is known")


colData(adult_se) |>
  as_tibble() |>
  group_by(disease, pathogen_curated, dataset) |>
  tally() |>
  group_by(disease, pathogen_curated) |>
  summarise(n = sum(n),
            n_datasets = n()) |>
  arrange(desc(n)) |>
  dplyr::filter(n >= 10,
                n_datasets >= 2,
                pathogen_curated != "unknown_pathogen") |>
  flextable::flextable() |>
  flextable::set_caption("Adult samples") |>
  flextable::add_footer_lines("All diseases with at least 10 samples and 2 or more datasets where the pathogen is known")

colData(ped_se) |>
  as_tibble() |>
  group_by(disease, pathogen_curated, dataset) |>
  tally() |>
  group_by(disease, pathogen_curated) |>
  summarise(n = sum(n),
            n_datasets = n()) |>
  arrange(desc(n)) |>
  dplyr::filter(n >= 10,
                n_datasets >= 2,
                pathogen_curated != "unknown_pathogen"
                ) |>
  flextable::flextable() |>
  flextable::set_caption("paediatric samples") |>
  flextable::add_footer_lines("All diseases with at least 10 samples and 2 or more datasets where the pathogen is known")



adult_selected <- colData(adult_se) |>
  as_tibble() |>
  group_by(disease, pathogen_curated, dataset) |>
  tally() |>
  group_by(disease, pathogen_curated) |>
  summarise(n = sum(n),
            n_datasets = n()) |>
  arrange(desc(n)) |>
  dplyr::filter(n >= 10,
                n_datasets >= 2,
                pathogen_curated != "unknown_pathogen") 

ped_selected <- colData(ped_se) |>
  as_tibble() |>
  group_by(disease, pathogen_curated, dataset) |>
  tally() |>
  group_by(disease, pathogen_curated) |>
  summarise(n = sum(n),
            n_datasets = n()) |>
  arrange(desc(n)) |>
  dplyr::filter(n >= 10,
                n_datasets >= 2,
                pathogen_curated != "unknown_pathogen") 

adult_selected_se <- adult_se[, colData(adult_se)$disease %in% adult_selected$disease]
ped_selected_se <- ped_se[, colData(ped_se)$disease %in% ped_selected$disease]
  
# Heatmaps ####
## Adults ####
meta <- as.data.frame(colData(adult_selected_se))

sample_cor <- cor(assay(adult_selected_se, "counts"), use = "na.or.complete")
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
  dplyr::select(disease, genome) 


pheatmap::pheatmap(sample_cor,
                   color = color.scheme,
                   cluster_rows = sample_tree,
                   cluster_cols = sample_tree,
                   annotation_col = col_annot, 
                   annotation_row = row_annot,
                   #annotation_legend = FALSE,
                   legend = FALSE,
                   show_rownames = FALSE,
                   show_colnames = FALSE)


## Kids ####
meta <- as.data.frame(colData(ped_selected_se))

sample_cor <- cor(assay(ped_selected_se, "counts"), use = "na.or.complete")
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
                   legend = FALSE,
                   show_rownames = FALSE,
                   show_colnames = FALSE)




# PCAs #####
## Adults #####
pca <- prcomp(na.omit(t(assay(adult_selected_se , "counts"))))
pcs_adult <- as_tibble(pca$x, rownames = "id") |>
  left_join(as.data.frame(colData(adult_selected_se)))

p1 <- pcs_adult |>
  ggplot(aes(PC1, PC2, colour = disease, 
             #shape = disease
  )) +
  geom_point(size = 1.5) +
  ggtitle("Adults - selected diseases")


xdens <-
  cowplot::axis_canvas(p1, axis = "x") +
  geom_density(data = pcs_adult, 
               aes(PC1, fill = dataset, colour = dataset), 
               alpha = .3)

ydens <- 
  cowplot::axis_canvas(p1, axis = "y", coord_flip = TRUE) +
  geom_density(data = pcs_adult, 
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




## Kids #####
pca <- prcomp(na.omit(t(assay(ped_selected_se , "counts"))))
pcs_ped <- as_tibble(pca$x, rownames = "id") |>
  left_join(as.data.frame(colData(ped_selected_se)))


p1 <- pcs_ped |>
  ggplot(aes(PC1, PC2, colour = disease, 
             #shape = disease
  )) +
  geom_point(size = 1.5) +
  ggtitle("Kids - selected diseases")


xdens <-
  cowplot::axis_canvas(p1, axis = "x") +
  geom_density(data = pcs_ped, 
               aes(PC1, fill = dataset, colour = dataset), 
               alpha = .3)

ydens <- 
  cowplot::axis_canvas(p1, axis = "y", coord_flip = TRUE) +
  geom_density(data = pcs_ped, 
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


# PCA densities #####
## Adults #####
to_plot <- c("sex", "disease", "source", 
             #"dataset", "pathogen_curated", 
             "disease_class", 
             "platform_id", "genome")


plot_pc_density <- function(group, dataset, PC="PC1") {
  dataset |>
    ggplot(aes(x = !!sym(PC), fill = !!sym(group), color = !!sym(group))) +
    geom_density(alpha = .5) + 
    ggtitle(group) +
    theme(legend.title=element_blank())
}


plots <- lapply(to_plot, 
                plot_pc_density, 
                dataset = pcs_adult, 
                PC = "PC1")
cowplot::plot_grid(plotlist = plots, ncol = 2)


plots_dataset <- lapply(c("PC1", "PC2"), 
                        plot_pc_density,
                        dataset = pcs_adult,
                        group = "dataset")


cowplot::plot_grid(plotlist = plots, ncol = 2)
cowplot::plot_grid(plotlist = plots_dataset, ncol = 1)

pcs_adult |>
  tidyr::pivot_longer(starts_with("PC")) |>
  filter(name %in% c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")) |>
  ggplot(aes(x = value, fill = disease, colour = disease)) + 
  geom_density(alpha = .3) +
  facet_wrap(~name, ncol = 2, scales = "free")



pcs_adult |>
  filter(grepl("37", genome)) |>
  count(dataset, genome, disease) |>
  flextable::flextable()

pcs_adult |>
  ggplot(aes(x = PC1, y = PC2, colour = source, shape = genome)) +
  geom_point(size = 2) 


# #kids#####

plots <- lapply(to_plot, 
                plot_pc_density, 
                dataset = pcs_ped, 
                PC = "PC1")
cowplot::plot_grid(plotlist = plots, ncol = 2)


plots_dataset <- lapply(c("PC1", "PC2"), 
                        plot_pc_density,
                        dataset = pcs_adult,
                        group = "dataset")


cowplot::plot_grid(plotlist = plots, ncol = 2)


pcs_ped |>
  tidyr::pivot_longer(starts_with("PC")) |>
  filter(name %in% c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")) |>
  ggplot(aes(x = value, fill = dataset, colour = dataset)) + 
  geom_density(alpha = .3) +
  facet_wrap(~name, ncol = 2, scales = "free")

  
  
pcs_ped |>
  filter(grepl("37", genome)) |>
  count(dataset, genome, disease) |>
  flextable::flextable()

pcs_ped |>
  ggplot(aes(x = PC1, y = PC2, colour = source, shape = genome)) +
  geom_point(size = 2) 

# no hg19 #####
## Adults#####

adult_selected_no_hg19 <- adult_selected_se[, colData(adult_selected_se)$genome != "GRCh37"]

pca_no_hg19 <- prcomp(na.omit(t(assay(adult_selected_no_hg19 , "counts"))))
pca_no_hg19 <- as_tibble(pca_no_hg19$x, rownames = "id") |>
  left_join(as.data.frame(colData(adult_selected_no_hg19)))


plots <- lapply(c(to_plot), 
                plot_pc_density, 
                dataset = pca_no_hg19, 
                PC = "PC1")
                
cowplot::plot_grid(plotlist = plots, ncol = 2)

pca_no_hg19 |>
  tidyr::pivot_longer(starts_with("PC")) |>
  filter(name %in% c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")) |>
  ggplot(aes(x = value, fill = disease, colour = disease)) + 
  geom_density(alpha = .3) +
  facet_wrap(~name, ncol = 2, scales = "free")


meta <- as.data.frame(colData(adult_selected_no_hg19))
sample_cor <- cor(assay(adult_selected_no_hg19, "counts"), use = "na.or.complete")
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
                   legend = FALSE,
                   show_rownames = FALSE,
                   show_colnames = FALSE)



colData(adult_selected_no_hg19) |>
  as_tibble() |>
  group_by(disease, pathogen_curated, dataset) |>
  tally() |>
  group_by(disease, pathogen_curated) |>
  summarise(n = sum(n),
            n_datasets = n()) |>
  arrange(desc(n)) |>
  dplyr::filter(n >= 10,
                n_datasets >= 2,
                pathogen_curated != "unknown_pathogen") |>
  flextable::flextable() |>
  flextable::set_caption("adult samples") |>
  flextable::add_footer_lines("All diseases with at least 10 samples and 2 or more datasets where the pathogen is known and the assembly is hg38")


meta <- colData(adult_selected_no_hg19)
batch <- as.numeric(factor(meta$dataset))
group <- as.numeric(factor(meta$disease))
counts <- assay(adult_selected_no_hg19)

adjusted_counts <- ComBat_seq(counts, 
                              batch=batch, 
                              group=group, 
                              full_mod=TRUE)

adjusted_counts_no_group <- ComBat_seq(counts, 
                                       batch=batch, 
                                       group=NULL, 
                                       full_mod=FALSE)

assays(adult_selected_no_hg19) <- list(counts = assay(adult_selected_no_hg19), 
                                       combat_seq = adjusted_counts,
                                       combat_seq_null = adjusted_counts_no_group)


pca_combat <- prcomp(na.omit(t(assay(adult_selected_no_hg19 , "combat_seq"))))
pca_combat <- as_tibble(pca_combat$x, rownames = "id") |>
  left_join(as.data.frame(colData(adult_selected_no_hg19)))


plots <- lapply(c(pca_combat, pca_no_hg19), 
                plot_pc_density, 
                group = "dataset", 
                PC = "PC1")

plot_pc_density("dataset", pca_no_hg19)
pca_combat |>
  ggplot(aes(x = PC1, y = PC2, colour = dataset)) +
  geom_point()


cowplot::plot_grid(plotlist = plots, ncol = 1)

assays(adult_selected_no_hg19)$combat_seq

saveRDS(adult_selected_no_hg19, here::here("data", "ses", "selected_adult_no_hg19.RDS"))
adult_selected_no_hg19 <- readRDS(here::here("data", "ses", "selected_adult_no_hg19.RDS"))



