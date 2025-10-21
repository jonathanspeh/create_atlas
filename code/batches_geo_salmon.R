library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
source("code/functions.R")

tximeta.se.1 <- readRDS("data/tximeta/SRP350281_tximeta.rds") # Filter to have only TP1
tximeta.se.2 <- readRDS("data/tximeta/SRP372307_tximeta.rds")
#tximeta.se.3 <- readRDS("data/tximeta/SRP460615_tximeta.rds")


geo.se.1 <- readRDS("data/ses/GSE190680_se.RDS")
geo.se.2 <- readRDS("data/ses/GSE201530_se.RDS")

##Join tximeta #####
meta.1 <- colData(geo.se.1)
meta.2 <- colData(geo.se.2)

meta.1 <- colData(tximeta.se.1) |>
  as.data.frame() |>
  right_join(as.data.frame(meta.1), by = join_by(sample_alias_ENA == id)) #|>
  rename(id = names)

meta.2 <- colData(tximeta.se.2) |>
  as.data.frame() |>
  right_join(as.data.frame(meta.2), by = join_by(sample_alias_ENA == id)) #|>
  rename(id = names)



rownames(meta.1) <- meta.1$id
rownames(meta.2) <- meta.2$id

tximeta.se.1 <- tximeta.se.1[,colnames(tximeta.se.1) %in% rownames(meta.1)]
colData(tximeta.se.1) <-DataFrame(meta.1)

tximeta.se.2 <- tximeta.se.2[,colnames(tximeta.se.2) %in% rownames(meta.2)]
colData(tximeta.se.2) <- DataFrame(meta.2)


se_list <- list(tximeta.se.1, tximeta.se.2)

tximeta.se.1$disease
tximeta.se.2$disease


tximeta_se <- combine_se(se_list)

## Join GEO #####
geo_se <- combine_se(list(geo.se.1, geo.se.2))
  
  
## Plot #####
counts_txi <- assay(tximeta_se)
meta_txi <- as.data.frame(colData(tximeta_se))

counts_txi <- counts_txi[rowSums(counts_txi) != 0,]
pca_obj_txi <- prcomp(t(counts_txi), scale. = TRUE)
pca_join_txi <- as_tibble(pca_obj_txi$x, rownames = "id") |>
  left_join(meta_txi)



p1_txi <- pca_join_txi |>
  ggplot(aes(x = PC1, y = PC2, colour = dataset)) +
  geom_point() + 
  stat_ellipse() +
  ggtitle("TXI")

p2_txi <- pca_join_txi |>
  ggplot(aes(x = PC1, y = PC2, colour = disease)) +
  geom_point() + 
  stat_ellipse() +
  ggtitle("txi")


counts_geo <- assay(geo_se)
meta_geo <- as.data.frame(colData(geo_se))

counts_geo <- counts_geo[rowSums(counts_geo) != 0,]


pca_obj_geo <- prcomp(t(counts_geo), scale. = TRUE)
pca_join_geo <- as_tibble(pca_obj_geo$x, rownames = "id") |>
  left_join(meta_geo)



p1_geo <- pca_join_geo |>
  ggplot(aes(x = PC1, y = PC2, colour = dataset)) +
  geom_point() + 
  stat_ellipse() +
  ggtitle("geo")

p2_geo <- pca_join_geo |>
  ggplot(aes(x = PC1, y = PC2, colour = disease)) +
  geom_point() + 
  stat_ellipse() +
  ggtitle("geo")


p1_geo
p1_txi

plot_pc_density <- function(group, dataset, PC="PC1") {
  dataset |>
    ggplot(aes(x = !!sym(PC), fill = !!sym(group), color = !!sym(group))) +
    geom_density(alpha = .5) + 
    ggtitle(group) +
    theme(legend.title=element_blank())
}

plot_pc_density("dataset", pca_join_txi, PC = "PC1")
plot_pc_density("dataset", pca_join_geo, PC = "PC1")


## Combat-Seq ####


#TODO










