## ---------------------------
## Script name: physeq_parse.R
## Author: McKenna Farmer
## Date Created: 2022-12-15
## ---------------------------
## Notes:
##   Making phyloseq parsing and saving to RDS 
##
## ---------------------------

library(qiime2R)
library(phyloseq)
library(tidyverse)
library(lubridate)

physeq <- qza_to_phyloseq(
  features="./qiime/table_dada2_all.qza",
  tree="./qiime/rooted_tree.qza",
  taxonomy="./qiime/taxonomy.qza",
  metadata = "./qiime/metadata_all.txt"
)

physeq <- subset_taxa(physeq, !Genus=="Mitochondria" & !Genus=="Chloroplast")

# rarefy
physeq2 <- rarefy_even_depth(physeq, sample.size=min(sample_sums(physeq)), rngseed=2)

# relative abundance
rel <- transform_sample_counts(physeq, function(x) x*100/sum(x))
rel_df <- psmelt(rel)
rel_df$date <- mdy(rel_df$date)

saveRDS(physeq, file=file.path("data","physeq.rds"))
saveRDS(physeq2, file=file.path("data","physeq2.rds"))
saveRDS(rel, file=file.path("data","rel.rds"))
saveRDS(rel_df, file=file.path("data","rel_df.rds"))

