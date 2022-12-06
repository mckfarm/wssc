## ---------------------------
## Script name: 16S_screening.R
## Author: McKenna Farmer
## Date Created: 2022-11-29
## ---------------------------
## Notes:
##   comparing test and control batteries
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

# alpha and beta diversity ---------
# rarefy
physeq2 <- rarefy_even_depth(physeq, sample.size=min(sample_sums(physeq)), rngseed=2)

## richness
observed <- plot_richness(physeq2, x="date", measures=c("Observed"))
observed_df <- as_tibble(observed$data)
observed_df$date <- mdy(observed_df$date)


ggplot(data=observed_df, aes(x=date, y=value, color=battery, shape=battery)) +
  geom_point() +
  theme_classic() +
  labs(x="Date", y="Number of ASVs observed", title="Richness") +
  color_battery +
  shape_battery +
  scale_x_main

ggsave(filename=file.path("results","screening","richness.tiff"), 
       height=4, width=5, units="in")


## shannon diversity
shannon <- plot_richness(physeq2, x="date", measures=c("Shannon"))

shannon_df <- as_tibble(shannon$data)
shannon_df$date <- mdy(shannon_df$date)

ggplot(data=shannon_df, aes(x=date, y=value, color=battery, shape=battery)) +
  geom_point() +
  theme_classic() +
  labs(x="Date", y="Shannon diversity index", title="Shannon diversity index") +
  color_battery +
  shape_battery +
  scale_x_main

ggsave(filename=file.path("results","screening","shannon.tiff"), 
       height=4, width=5, units="in")


# beta diversity
ord <- ordinate(physeq2,  "NMDS", "unifrac", weighted=TRUE)
plot_ordination(physeq2, ord, color="battery", label="timepoint") +
   theme_bw() + 
  color_battery +
  shape_battery

ggsave(filename=file.path("results","screening","ordination.tiff"), 
       height=6, width=7, units="in")


# relative abundance
rel <- transform_sample_counts(physeq, function(x) x*100/sum(x))

rel_df <- psmelt(rel)
rel_df$date <- mdy(rel_df$date)
rel_df$battery <- factor(rel_df$battery,levels=c("control","test"))


# PAO/GAO --------


rel_phos_sum <- rel_df %>% 
  filter(Genus %in% pao_gao_select) %>% 
  filter(! OTU %in% asv_dechloro_remove) %>%
  group_by(date,Genus,battery) %>%
  summarise(sum=sum(Abundance))

rel_phos_sum$Genus <- factor(rel_phos_sum$Genus, levels=pao_gao_list)

ggplot(data=rel_phos_sum, aes(x=date, y=sum, fill=Genus)) +
  facet_wrap(~battery) +
  geom_bar(stat="identity")

# compare PAO and GAO abundances from control and test batteries -------
control_dates <- unique(rel_df$date[rel_df$battery =="control"])

rel_phos_compare <- rel_phos_sum %>% filter(date %in% control_dates)

ggplot(data=rel_phos_compare, aes(x=battery, y=sum, color=battery)) +
  facet_wrap(~Genus) + 
  theme_classic() + 
  geom_boxplot(outlier.color=NA) +
  geom_point() +
  labs(y="relative ab") +
  color_battery + 
  stat_compare_means(label="p.signif")
ggsave(filename=file.path("results","screening","paogao_controltest.tiff"), 
       height=6, width=7, units="in")



# Nitrifiers -------------

rel_nit_sum <- rel_df %>% 
  filter(Genus %in% nitrifier_list) %>%
  group_by(date,Genus,battery) %>%
  summarise(sum=sum(Abundance))

rel_nit_compare <- rel_nit_sum %>%  filter(date %in% control_dates)

ggplot(data=rel_nit_compare, aes(x=battery, y=sum)) +
  facet_wrap(~Genus) + 
  geom_boxplot(outlier.color=NA) +
  geom_point() +
  labs(y="relative ab") +
  stat_compare_means(label="p.signif")
ggsave(filename=file.path("results","screening","nit_controltest.tiff"), 
       height=2, width=4, units="in")

