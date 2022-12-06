## ---------------------------
## Script name: 16S_analysis.R
## Author: McKenna Farmer
## Date Created: 2022-11-29
## ---------------------------
## Notes:
## running analysis and calculations on test basin only
##
## ---------------------------


# packages ------------
## data manipulation
library(qiime2R)
library(phyloseq)
library(readxl)
library(tidyverse)
library(lubridate)

# 
# library(zoo)
library(ggpubr)
# library(microbiome)

## plotting
library(ggplot2)
library(MetBrewer)
library(cowplot)
# library(ggcorrplot)
# library(gridExtra)
# library(ggdist)
# library(ggConvexHull)

source(file.path("scripts","plotting.R"))


# dechloro data read in ------
# reading in asv names that were identified as dechloromonas PAO
asv_troph <- readRDS(file=file.path("data","asv_troph.rds"))
asv_vans <- readRDS(file=file.path("data","asv_vans.rds"))
asv_dechloro_all <- readRDS(file=file.path("data","asv_dechloro_all.rds"))

asv_dechloro_pao <- c(asv_troph,asv_vans)
asv_dechloro_remove <- setdiff(asv_dechloro_all, asv_dechloro_pao)


# 16S data read in and basic parsing ------
physeq <- qza_to_phyloseq(
  features="./qiime/table_dada2_all.qza",
  tree="./qiime/rooted_tree.qza",
  taxonomy="./qiime/taxonomy.qza",
  metadata = "./qiime/metadata_all.txt"
)

physeq <- subset_taxa(physeq, !Genus=="Mitochondria" & !Genus=="Chloroplast")
physeq <- subset_samples(physeq, battery=="test") # test only
physeq2 <- rarefy_even_depth(physeq, sample.size=min(sample_sums(physeq)), rngseed=2)


# relative abundance
rel <- transform_sample_counts(physeq, function(x) x*100/sum(x))

rel_df <- psmelt(rel)
rel_df$date <- mdy(rel_df$date)


# PAO/GAO --------
# parsing
rel_phos_sum <- rel_df %>% 
  filter(Genus %in% pao_gao_select) %>% 
  filter(! OTU %in% asv_dechloro_remove) %>%
  group_by(date,Genus) %>%
  summarise(sum=sum(Abundance))

rel_phos_sum$Genus <- factor(rel_phos_sum$Genus, levels=pao_gao_list)

# stats
rel_phos_wide <- rel_phos_sum %>%
  pivot_wider(id_cols=date, names_from=Genus, values_from=sum)

# plots
ggplot(data=rel_phos_sum, aes(x=date, y=sum, fill=Genus)) +
  geom_bar(stat="identity") + 
  theme_classic() +
  scale_fill_paogao + 
  scale_x_main +
  ylim(0,15) + 
  labs(x="Date", y="Relative abundance [%]") +
  theme(legend.key.size = unit(.75,"line"),
        legend.text=element_text(size=8),
        legend.title=element_text(size=9),
        legend.position=c(0.15,0.8),
        legend.box.background = element_rect(colour = "black"))
ggsave(filename=file.path("results","paogao_barplot.tiff"), 
       height=4, width=6, units="in")


ggplot(data=rel_phos_sum, aes(x=Genus, y=sum, color=Genus)) +
  geom_boxplot(outlier.color=NA) + 
  geom_point(alpha=0.5) + 
  theme_classic() +
  coord_flip() + 
  scale_color_paogao +
  scale_x_paogao + 
  labs(x="Genus", y="Relative abundance [%]") +
  theme(legend.position="none")
ggsave(filename=file.path("results","paogao_boxplot.tiff"), 
       height=4, width=6, units="in")


# Nitrifiers -------------

rel_nit_sum <- rel_df %>% 
  filter(Genus %in% nitrifier_list) %>%
  group_by(date,Genus) %>%
  summarise(sum=sum(Abundance))

ggplot(data=rel_nit_sum, aes(x=date, y=sum, fill=Genus)) +
  geom_bar(stat="identity") + 
  theme_classic() +
  scale_fill_manual(values=met.brewer("Egypt")) + 
  scale_x_main +
  ylim(0,5) + 
  labs(x="Date", y="Relative abundance [%]") +
  theme(legend.key.size = unit(.75,"line"),
        legend.text=element_text(size=8),
        legend.title=element_text(size=9),
        legend.position=c(0.13,0.85),
        legend.box.background = element_rect(colour = "black"))

ggsave(filename=file.path("results","nitrifier_barplot.tiff"), 
       height=4, width=6, units="in")

