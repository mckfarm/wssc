## ---------------------------
## Script name: comp_blast_parse.R
## Author: McKenna Farmer
## Date Created: 2022-12-01
## ---------------------------
## Notes:
## parsing blast results, where competibacter denitrificans reference sequences were made into a blastdb
##
## ---------------------------

# packages -------
library(readr)
library(dplyr)


# read in --------
blast_colnames <- c("qseqid", "sseqid", "pident", "length", "mismatch",
                    "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")


blast_comp <- read_table(file.path("results","blast_comp_denit.txt"), 
                         col_names = blast_colnames)

# view blast comp
asv_comp_match <- c("6dfa1d3cd280ae1124f11510300535db")
saveRDS(asv_comp_match, file=file.path("data","asv_comp_denit.rds"))

 
# check for denit relevance
rel_phos <- readRDS(file=file.path("data","rel_phos.rds"))

rel_phos_comp <- rel_phos %>% filter(Genus=="Ca_Competibacter")
