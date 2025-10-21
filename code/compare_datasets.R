library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)
library(DESeq2)

tidyverse::

# Load Data #####
from_GEO <- readRDS("data/ses/GSE201530_se.RDS")
tximport <- readRDS("data/tximports/SRP372307_tximport.rds") #FIXME - Uses wrong tx2gene
tximeta <- readRDS("data/tximeta/SRP372307_tximeta.rds") 

fr <- read_delim("data/filereport_read_run_SRP372307_tsv.txt") |>
  dplyr::select(run_accession, sample_alias)



#from_GEO <- readRDS("data/ses/GSE190680_se.RDS")
#tximport <- readRDS("data/tximports/SRP350281_tximport.rds") 
#tximeta <- readRDS("data/tximeta/SRP350281_tximeta.rds")
#
#fr <- colData(tximeta) |>
#  as_tibble(rownames = "run_accession") |>
#  dplyr::select(run_accession, sample_alias = sample_alias_ENA)

colData(from_GEO)



GEO_counts <- assay(from_GEO) |> 
  as_tibble(rownames = "Gene") |>
  pivot_longer(!Gene, names_to = "sample", values_to = "count_geo")

salmon_counts_tximport <- tximport$counts |> 
  as_tibble(rownames = "Gene") |>
  pivot_longer(!Gene, names_to = "sample", values_to = "count_salmon_tximport") |>
  mutate(Gene = stringr::str_remove(Gene, "\\.[^.]+$")) |>
  left_join(fr, by = join_by(sample == run_accession))

salmon_counts_tximeta <- assay(tximeta, "counts") |> 
  as_tibble(rownames = "Gene") |>
  pivot_longer(!Gene, names_to = "sample", values_to = "count_salmon_tximeta") |>
  mutate(Gene = stringr::str_remove(Gene, "\\.[^.]+$")) |>
  left_join(fr, by = join_by(sample == run_accession))



joined <- inner_join(GEO_counts, salmon_counts_tximport, 
                     by = join_by(Gene, sample == sample_alias))

joined <- inner_join(joined, salmon_counts_tximeta, 
                     by = join_by(Gene, sample == sample_alias))


joined_salmon <- full_join(salmon_counts_tximeta, salmon_counts_tximport)


# Compare Distributeions #####
joined.non.zero <- joined |>
  dplyr::filter(count_geo > 0,
         count_salmon_tximport > 0,
         count_salmon_tximeta > 0)
  
joined.non.zero |>
  dplyr::filter(count_geo < max(count_salmon_tximeta)) |>
  slice_sample(prop = .01) |> 
  pivot_longer(starts_with("count")) |>
  ggplot(aes(y = value, x = name, colour = name)) +
  geom_violin() + 
  geom_jitter(alpha = .1) + 
  scale_y_log10()


joined.non.zero |>  
  dplyr::filter(count_geo < max(count_salmon_tximport)) |>
  slice_sample(prop = .1) |>
  ggplot(aes(x = count_geo,
             y = count_salmon_tximport)) +
  geom_point() +
  geom_abline(colour = "red")

joined.non.zero  |>  
  #filter(count_geo < max(count_salmon_tximport)) |>
  slice_sample(prop = .1) |>
  ggplot(aes(x = count_geo,
             y = count_salmon_tximeta)) +
  geom_point() +
  geom_abline(colour = "red")

joined.non.zero  |>  
  #filter(count_geo < max(count_salmon_tximport)) |>
  slice_sample(prop = .1) |>
  ggplot(aes(x = count_salmon_tximeta,
             y = count_salmon_tximport)) +
  geom_point() +
  geom_abline(colour = "red")


cor.test(joined.non.zero$count_geo,
         joined.non.zero$count_salmon_tximport)

cor.test(joined.non.zero$count_geo,
         joined.non.zero$count_salmon_tximeta)


joined |> 
  dplyr::filter(count_geo < max(count_salmon_tximeta)) |>
  slice_sample(prop = .01) |> 
  pivot_longer(starts_with("count")) |>
  mutate(log_count = log2(value+1)) |>
  ggplot(aes(x = log_count, colour = name, fill = name)) +
  geom_density(alpha = .5) 

joined |> 
  slice_sample(prop = .01) |> 
  pivot_longer(starts_with("count")) |>
  filter(value < 10000) |>
  ggplot(aes(x = value, colour = name, fill = name)) +
  geom_density(alpha = .5) 


joined |> 
  slice_sample(prop = .01) |> 
  pivot_longer(starts_with("count")) |>
  dplyr::filter(value < 10000) |>
  mutate(log_count = log2(value+1)) |>
  ggplot(aes(x = log_count, colour = sample)) +
  geom_density(alpha = .5) +
  facet_wrap(~name)



# Prpare DESeq2 Objects ####

dds_geo <- DESeqDataSet(from_GEO, design = ~disease)
colData(dds_geo)

colDat_tximeta <- as.data.frame(colData(tximeta))

colDat_tximeta <- 
  left_join(colDat_tximeta, as_tibble(colData(dds_geo)),
          by = join_by(sample_alias_ENA == id))



rownames(colDat_tximeta) <- colDat_tximeta$names
all(rownames(colDat_tximeta) == colnames(tximeta))


colData(tximeta) <- DataFrame(colDat_tximeta)
dds_slm <- DESeqDataSet(tximeta, design = ~disease)

rename_vec <- fr$run_accession
names(rename_vec) <- fr$sample_alias



smallestGroupSize <- 3
keep <- rowSums(counts(dds_geo) >= 10) >= smallestGroupSize
dds_geo <- dds_geo[keep,]

keep <- rowSums(counts(dds_slm) >= 10) >= smallestGroupSize
dds_slm <- dds_slm[keep,]


dds_geo <- DESeq(dds_geo)
dds_slm <- DESeq(dds_slm)

res_geo <- results(dds_geo, contrast = c("disease", "COVID-19", "healthy"))
res_slm <- results(dds_slm, contrast = c("disease", "COVID-19", "healthy"))


LFC_geo <- lfcShrink(dds_geo, 
                     coef = "disease_healthy_vs_COVID.19", 
                     type = "apeglm")


LFC_slm <- lfcShrink(dds_slm, 
                     coef = "disease_healthy_vs_COVID.19", 
                     type = "apeglm")


nrow(dds_geo)
nrow(dds_slm)

plotMA(dds_geo)
plotMA(dds_slm)

plotMA(LFC_geo)
plotMA(LFC_slm)


as_tibble(res_geo) |>
  mutate(group = case_when(
    padj >= 0.05 ~ "padj above treshold",
    abs(log2FoldChange) > 1 ~ "log2FC > 1",
    TRUE ~ "log2FC < 1")) |>
  ggplot(aes(x = log2FoldChange, 
             y = -log10(padj),
             colour = group)) +
  geom_point() +
  scale_color_discrete(palette = c("blue", "darkred", "gray"))



as_tibble(res_slm) |>
  mutate(group = case_when(
    padj >= 0.05 ~ "padj above treshold",
    abs(log2FoldChange) > 1 ~ "log2FC > 1",
    TRUE ~ "log2FC < 1")) |>
  ggplot(aes(x = log2FoldChange, 
             y = -log10(padj),
             colour = group)) +
  geom_point() +
  scale_color_discrete(palette = c("blue", "darkred", "gray"))


res_slm |>
  as_tibble(rownames = "gene") |>
  dplyr::filter(log2FoldChange > 20)



slm_sig <- as_tibble(res_slm, rownames = "gene") |>
  dplyr::filter(padj < 0.05,
                abs(log2FoldChange) > 1) |>
  mutate(gene = stringr::str_remove(gene, "\\.[^.]+$"))

geo_sig <- as_tibble(res_geo, rownames = "gene") |>
  dplyr::filter(padj < 0.05,
                abs(log2FoldChange) > 1)




length(intersect(slm_sig$gene,
          geo_sig$gene))


nrow(slm_sig)
nrow(geo_sig)


gene_list <-list(salmon = slm_sig$gene,
           geo = geo_sig$gene)


library(ggVennDiagram)

ggVennDiagram(gene_list) +
  scale_fill_gradient(low="grey90",high = "red")










