## ---------------------------
## Script name: test_only.R
## Author: McKenna Farmer
## Date Created: 2022-11-29
## ---------------------------
## Notes:
## running analysis and calculations on test basin only
##
## ---------------------------


# packages ------------
## data manipulation
library(phyloseq)
library(readxl)
library(tidyverse)
library(lubridate)

## plotting
library(ggplot2)
library(MetBrewer)
library(patchwork)
library(ggpubr)

source(file.path("scripts","plotting_manuscript.R"))


# dechloro data read in ------
# reading in asv names that were identified as dechloromonas PAO
asv_troph <- readRDS(file=file.path("data","asv_troph.rds"))
asv_vans <- readRDS(file=file.path("data","asv_vans.rds"))
asv_dechloro_all <- readRDS(file=file.path("data","asv_dechloro_all.rds"))
asv_comp_denit <- readRDS(file=file.path("data","asv_comp_denit.rds"))

asv_dechloro_pao <- c(asv_troph,asv_vans)
asv_dechloro_remove <- setdiff(asv_dechloro_all, asv_dechloro_pao)


# 16S data read in and basic parsing ------
rel_df <- readRDS(file.path("data","rel_df.rds"))
rel_df <- rel_df %>% filter(battery=="test")

date_range <- c(min(rel_df$date), max(rel_df$date))


# PAO/GAO --------
# parsing

rel_phos_sum <- rel_df %>% 
  filter(Genus %in% pao_gao_list) %>% 
  filter(! OTU %in% asv_dechloro_remove) %>%
  group_by(date,Genus) %>%
  summarise(sum=sum(Abundance))

tapply(rel_phos_sum$sum, rel_phos_sum$Genus, median)
tapply(rel_phos_sum$sum, rel_phos_sum$Genus, max)

rel_phos_sum$Genus <- factor(rel_phos_sum$Genus, levels=pao_gao_list)


# plots
## all boxplot
pao_gao_box <- 
ggplot(data=rel_phos_sum, aes(x=Genus, y=sum, color=Genus)) +
  geom_boxplot(outlier.color=NA) + 
  geom_point(alpha=0.5) + 
  theme_classic() +
  scale_color_paogao +
  scale_x_paogao + 
  coord_flip() + 
  labs(x="Genus", y="Relative abundance [%]") +
  theme(legend.position="none")

## accumulibacter, tetrasphaera, and competibacter only
pao_gao_select_bars <- 
rel_phos_sum %>% filter(Genus %in% pao_gao_select) %>%
  ggplot(data=., aes(x=date, y=sum, fill=Genus)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    scale_fill_paogao_select + 
    scale_x_main +
    ylim(0, 15) + 
    labs(x="Date", y="Relative abundance [%]") +
    theme(legend.position=c(0.2, 0.85),
          legend.background = element_rect(colour = 'grey', fill = 'white', linetype='solid')) +
    guides(fill = guide_legend(override.aes = list(size = 0.5)))

pao_gao_box + pao_gao_select_bars + plot_layout(widths = c(0.7, 1)) + plot_annotation(tag_levels = "A")

ggsave(filename=file.path("results/manuscript","paogao_bygenus_test.png"), 
       height=5, width=10, units="in", dpi=300)


# Nitrifiers -------------

rel_nit_sum <- rel_df %>% 
  filter(Genus %in% nitrifier_list) %>%
  group_by(date,Genus) %>%
  summarise(sum=sum(Abundance))

rel_nitrospira <- rel_df %>% 
  filter(Genus == "Nitrospira") %>%
  group_by(date,Species) %>%
  summarise(sum=sum(Abundance)) %>%
  replace_na(list(Species="Unknown")) %>%
  filter(! Species %in% c("midas_s_21815", "midas_s_22526"))

nitrifier_box <- ggplot(rel_nit_sum, aes(x=Genus, y=sum, color=Genus)) +
  geom_boxplot() +
  geom_point(alpha=0.5) + 
  theme_classic() +
  scale_color_manual(values = met.brewer("Egypt", 3)) + 
  ylim(0, 4) + 
  coord_flip() + 
  labs(y="", x = "") +
  theme(legend.position="none")

nitrospira_box <- ggplot(rel_nitrospira, aes(x=Species, y=sum, fill=Species)) +
  geom_point(alpha=0.3) + 
  geom_boxplot() +
  theme_classic() + 
  scale_x_nitrospira + 
  scale_fill_nitrospira +
  ylim(0, 4) + 
  coord_flip() + 
  labs(y="Relative abundance [%]", x = "") +
  theme(legend.position="none")

(nitrifier_box / nitrospira_box) + plot_annotation(tag_levels = "A")

ggsave(filename=file.path("results/manuscript","nitrospira_combo.png"), 
       height=5, width=5, units="in", dpi=600)

