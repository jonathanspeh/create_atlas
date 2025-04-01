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


pca_raw_obj = prcomp(t(counts))
pca_raw_join <- as_tibble(pca_raw_obj$x, rownames = "id") |>
  left_join(meta)

pca_adjusted_obj = prcomp(t(adjusted_counts))
pca_adjusted_join <- as_tibble(pca_adjusted_obj$x, rownames = "id") |>
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

cowplot::plot_grid(p1, p2, ncol =1)
cowplot::plot_grid(p1_adj, p2_adj, ncol =1)

cowplot::plot_grid(p1, p1_adj, ncol =1)
cowplot::plot_grid(p2, p2_adj, ncol =1)




plot_pc_density <- function(group, dataset, PC="PC1") {
  dataset |>
    ggplot(aes(x = !!sym(PC), fill = !!sym(group), color = !!sym(group))) +
    geom_density(alpha = .5) + 
    ggtitle(group) +
    theme(legend.title=element_blank())
}


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




plots_raw_2 <- lapply(c("PC1", "PC2", "PC3", "PC4"), 
                    plot_pc_density, 
                    group = "dataset",
                    dataset = pca_raw_join)

cowplot::plot_grid(plotlist = plots_raw_2, ncol = 2) 


plots_adjusted_2 <- lapply(c("PC1", "PC2", "PC3", "PC4"), 
                         plot_pc_density, 
                         group = "dataset",
                         dataset = pca_adjusted_join)

cowplot::plot_grid(plotlist = plots_adjusted_2, ncol = 1) 




