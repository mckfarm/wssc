## ---------------------------
## Script name: dechloro_prep.R
## Author: McKenna Farmer
## Date Created: 2022-10-26
## ---------------------------
## Notes:
## preparing dechloromonas sequences for blast comparison
##
## ---------------------------

library(Biostrings)
library(readr)
library(dplyr)
library(qiime2R)
library(phyloseq)

# read in
## midas
midas_seqs <- readDNAStringSet("./data/QIIMEMiDAS_4.8.1.fa")
midas_taxonomy <- read_csv("./data/midas_taxonomy.csv")

## picking reference seqs from midas
phosphoritropha <- midas_taxonomy %>%
  filter(species=="Ca_Dechloromonas_phosphoritropha")

phosphoritropha_ref_seqs <- midas_seqs[names(midas_seqs) %in% phosphoritropha$id]

writeXStringSet(phosphoritropha_ref_seqs, "ref_phosphoritropha.fasta")


phosphorivorans <- midas_taxonomy %>%
  filter(species=="Ca_Dechloromonas_phosphorivorans")

phosphorivorans_ref_seqs <- midas_seqs[names(midas_seqs) %in% phosphorivorans$id]

writeXStringSet(phosphorivorans_ref_seqs, "ref_phosphorivorans.fasta")


# sample data ----------
physeq <- qza_to_phyloseq(
  features="./qiime/table_dada2_all.qza",
  tree="./qiime/rooted_tree.qza",
  taxonomy="./qiime/taxonomy.qza",
  metadata = "./qiime/metadata_all.txt"
)

asvs <- readDNAStringSet("./qiime/sequences.fasta")

physeq@refseq <- asvs

# subset to dechloromonas
physeq_dechloro <- subset_taxa(physeq, Genus=="Dechloromonas")

# extract dechloro sequences
dechloro_seqs <- physeq_dechloro@refseq

asv_dechloro_all <- dechloro_seqs@ranges@NAMES

saveRDS(asv_dechloro_all, file=file.path("data","asv_dechloro_all.rds"))
  
writeXStringSet(dechloro_seqs, "seqs_dechloro.fasta")

