## ---------------------------
## Script name: dechloro_blast_parse.R
## Author: McKenna Farmer
## Date Created: 2022-12-01
## ---------------------------
## Notes:
## parsing blast results, where dechloromonas reference sequences were made into a blastdb
##
## ---------------------------

# packages -------
library(readr)
library(dplyr)


# read in --------
blast_colnames <- c("qseqid", "sseqid", "pident", "length", "mismatch",
                    "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")


blast_troph <- read_table(file.path("results","blast_troph.txt"), col_names = blast_colnames)
blast_vans <- read_table(file.path("results","blast_vans.txt"), col_names = blast_colnames)

blast_troph$species <- "troph"
blast_vans$species <- "vans"


# parse -------

blast_troph_match <- blast_troph %>% filter(pident >= 99) # species classification
blast_vans_match <- blast_vans %>% filter(pident >= 99) # species classification

asv_troph_match <- unique(blast_troph_match$qseqid)
asv_vans_match <- unique(blast_vans_match$qseqid)

intersect(asv_troph_match, asv_vans_match) # check for overlap

# making the call that ASV 8717744be592c8beffade9c24a1887b2 is a closer match to phosphorivorans
# due to multiple 100% identity matches

asv_troph_match <- asv_troph_match[asv_troph_match!="8717744be592c8beffade9c24a1887b2"]

# save for reference in 16S analysis files
saveRDS(asv_troph_match, file=file.path("data","asv_troph.rds"))
saveRDS(asv_vans_match, file=file.path("data","asv_vans.rds"))


