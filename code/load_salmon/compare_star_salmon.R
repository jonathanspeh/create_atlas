library(tidyverse)

rename_counts <- function(count_df) {
  colnames(count_df) <- ifelse(str_detect(colnames(counts), "SRP"),
                               str_extract(colnames(counts), "SRR[0-9]*"),
                               colnames(counts)
  )
  count_df
}



counts <- read.table("/Volumes/sc-project-cc12-pbmc_lens/infection_atlas/star_test/SRP350281/counts.txt",
           header = TRUE,
           skip = 1
           )

salmon_counts <- read.table("/Volumes/sc-project-cc12-pbmc_lens/infection_atlas/star_test/SRP350281/salmon_quant/SRR17197391_trimmed/quant.sf",
                            header = TRUE
                            )
salmon_counts |>
  separate_wider_delim(Name, delim = ".", names=c("id", "suffix")) |>
  dplyr::count(id, suffix) |>
  dplyr::filter(n>1)

#salmon_counts <- salmon_counts |>
#  separate_wider_delim(Name, delim = "|", names = c("Ensembl_transcript", "Ensembl_gene", "Genecode_gene", "Genecode_transcript",
#                                                    "HGNC_Transcript", "HGNC_Symbol", "Entrez_id", "Biotype", NA),
#                       too_many = "debug"
#                       )

salmon_counts_proc <- salmon_counts |>
  separate_wider_delim(Name, ".", names = c("Geneid", NA)) |>
  group_by(Geneid) |>
  summarise(count = sum(NumReads))



se_original <- readRDS("/Users/jonathanspeh/Desktop/create_atlas/data/ses/GSE190680_se.RDS")
run_tbl <- read.csv("/Users/jonathanspeh/Downloads/SraRunTable(1).csv") |> 
  select(Run, Sample.Name)

glimpse(run_tbl)


counts_original <- assay(se_original)
colData(se_original)$processing_info[1] |> glimpse()

counts <- rename_counts(counts)

head(counts)
counts_long <- counts |>
  select("Geneid", starts_with("SRR")) |>
  pivot_longer(!"Geneid", 
               names_to="sample", 
               values_to="count")



counts_original_long <- counts_original |> 
  as_tibble(rownames = "Geneid") |>
  pivot_longer(!"Geneid", 
               names_to="sample", 
               values_to="count")

counts_original_long <- counts_original_long |>
  left_join(run_tbl, by = join_by("sample" == "Sample.Name"))



counts_joined <- counts_long |>
  left_join(counts_original_long, by = join_by(Geneid, "sample" == "Run"))

options(scipen = 10000)

acc_list <- colnames(counts)[7:17]

counts_joined |>
  #filter(Geneid %in% keep$Geneid) |>
  filter(!is.na(count.y),
         sample %in% acc_list) |>
  ggplot(aes(x = count.x, y = count.y, colour = sample)) +
  geom_point() +
  geom_abline(slope = 1) + 
  facet_wrap(~sample, scales = "free")


cor.test(counts_joined$count.x,
         counts_joined$count.y,
         na.action = "na.omit")


counts_long |>
  group_by(Geneid) |>
  summarise(sum_count = sum(count)) |>
  filter(sum_count >=10) -> keep


nrow(counts_original)



salmon_star_joined <- counts_long |>
  filter(sample == "SRR17197391") |>
  left_join(salmon_counts_proc, by = join_by(Geneid==Geneid)) 

salmon_star_joined |>
  #filter(count.x < 50000,
  #       count.y < 50000) |>
  ggplot(aes(x = count.x, y = count.y)) +
  geom_point() +
  geom_abline(slope = 1) 


cor.test(salmon_star_joined$count.x,
         salmon_star_joined$count.y
         )



mean(round(salmon_counts_proc$count) == salmon_counts_proc$count)


library(tximport)
library(EnsDb.Hsapiens.v86)

tx2gene <- ensembldb::transcripts(EnsDb.Hsapiens.v86,
                       return.type="data.frame"
                       ) |>
  dplyr::select(tx_id, gene_id)


#tx2gene <- data.frame(salmon_counts$Name, salmon_counts$Ensembl_gene)

#tx2gene$salmon_counts.Ensembl_gene <- str_remove(tx2gene$salmon_counts.Ensembl_gene, "\\..*")


txi <- tximport("/Volumes/sc-project-cc12-pbmc_lens/infection_atlas/star_test/SRP350281/salmon_quant/SRR17197391/quant.sf",
                type = "salmon", 
                tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM"
                )

names(txi)

txi$countsFromAbundance
proc <- as_tibble(txi$counts, rownames = "Geneid") |>
  dplyr::rename(count = V1) 

counts_long |>
  filter(sample == "SRR17197391") |>
  left_join(proc, by = "Geneid") |>
  filter(count.x < 50000,
         count.y < 50000) |>
  ggplot(aes(x = count.x, y = count.y)) +
  geom_point()



salmon_star_joined |>
  filter(count.x < 50000,
         count.y < 50000) |>
  ggplot(aes(x = count.x, y = count.y)) +
  geom_point() +
  geom_abline(slope = 1) 


salmon_star_joined |> arrange(Geneid) |> tail()
txi$counts |> tail()
