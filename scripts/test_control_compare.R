## ---------------------------
## Script name: test_control_compare.R
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

# read in ----
physeq2 <- readRDS(file.path("data","physeq.rds")) #rarefied

rel_df <- readRDS(file.path("data","rel_df.rds"))

rel_df$battery <- factor(rel_df$battery,levels=c("control","test"))

control_dates <- unique(rel_df$date[rel_df$battery =="control"])

# alpha and beta diversity ------


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

shannon_df <- shannon_df %>% filter(date %in% control_dates)

ggplot(data=shannon_df, aes(x=battery, y=value, color=battery)) +
  geom_boxplot(outlier.color=NA) + 
  geom_point() +
  theme_classic() +
  labs(x="Battery", y="Shannon diversity index", title="") +
  theme(legend.position="none") +
  color_battery +
  scale_x_battery

ggsave(filename=file.path("results", "shannon_boxplot.png"), 
       height=3.5, width=4, units="in", dpi=500)

# beta diversity
ord <- ordinate(physeq2,  "NMDS", "unifrac", weighted=TRUE)
plot_ordination(physeq2, ord, color="battery", label="timepoint") +
   theme_bw() + 
  color_battery +
  shape_battery +
  annotate(geom="text", x=0.035, y=-0.03, label="Stress = 0.08", size=3)

ggsave(filename=file.path("results", "ordination.png"), 
       height=4, width=5, units="in")


# functional groups
# PAO/GAO --------

# reading in asv names that were identified as dechloromonas PAO
asv_troph <- readRDS(file=file.path("data","asv_troph.rds"))
asv_vans <- readRDS(file=file.path("data","asv_vans.rds"))
asv_dechloro_all <- readRDS(file=file.path("data","asv_dechloro_all.rds"))

asv_dechloro_pao <- c(asv_troph,asv_vans)
asv_dechloro_remove <- setdiff(asv_dechloro_all, asv_dechloro_pao)

# parse
rel_phos <- rel_df %>% 
  filter(Genus %in% pao_gao_list) %>% 
  filter(! OTU %in% asv_dechloro_remove) 

rel_phos_sum <- rel_df %>% 
  filter(Genus %in% pao_gao_select) %>% 
  filter(! OTU %in% asv_dechloro_remove) %>%
  group_by(date,Genus,battery) %>%
  summarise(sum=sum(Abundance))

rel_phos_sum$Genus <- factor(rel_phos_sum$Genus, levels=pao_gao_list)

rel_phos_sum %>% group_by(Genus, battery) %>% summarise(median=median(sum))

# compare PAO and GAO abundances from control and test batteries

rel_phos_compare <- rel_phos_sum %>% filter(date %in% control_dates)

ggplot(data=rel_phos_compare, aes(x=battery, y=sum, color=battery)) +
  facet_wrap(~Genus, labeller=labels_paogao) + 
  theme_bw() + 
  geom_boxplot(outlier.color=NA) +
  geom_point() +
  labs(y="relative ab") +
  color_battery + 
  theme(legend.position="none") +
  labs(x="Battery", y="Relative abundance [%]")
ggsave(filename=file.path("results","screening","paogao_controltest.tiff"), 
       height=2.5, width=4, units="in")


# plot against each other 

rel_phos_compare %>%
  pivot_wider(names_from=battery,values_from=sum) %>%
  ggplot(data=., aes(x=test, y=control, color=Genus)) +
  facet_wrap(~Genus, labeller=labels_paogao, nrow=3) +
  geom_point() +
  theme_classic() +
  theme(legend.position="none") + 
  scale_color_paogao_select + 
  labs(y="Relative abundance in control basin [%]",
       x="Relative abundance in test basin [%]") +
  stat_cor(method="spearman", label.x=c(4,0.5,4), label.y=7)
ggsave(filename=file.path("results","paogao_controltest_corr.tiff"), 
       height=5, width=3.5, units="in")


# glycogen stats

gly_dates <- read_rds(file=file.path("data","gly_dates.rds"))
gly_dates <- sort(gly_dates)

rel_phos_sum %>% 
  filter(Genus %in% pao_gao_select) %>%
  ggplot(., aes(x=date, y=sum, color=Genus)) +
  facet_wrap(~battery, nrow=2, labeller=labels_basins) + 
  geom_point() +
  geom_line() + 
  scale_color_paogao_select +
  theme_classic() +
  theme(legend.position="top") +
  ylim(0,8.5) + 
  labs(y="Relative abundance [%]", x="Date") +
  annotate(geom = "text", x = as_date(gly_dates[1]),
           y = 8.2, label = "1") +
  annotate(geom = "text", x = as_date(gly_dates[2]),
           y = 8.2, label = "2") +
  annotate(geom = "text", x = as_date(gly_dates[3]),
           y = 8.2, label = "3") +
  annotate(geom = "text", x = as_date(gly_dates[4]),
           y = 8.2, label = "4") +
  annotate(geom = "text", x = as_date(gly_dates[5]),
           y = 8.2, label = "5") +
  annotate(geom = "text", x = as_date(gly_dates[6]),
           y = 8.2, label = "6")
  
ggsave(filename=file.path("results","screening","paogao_select_ct.tiff"), 
       height=5, width=6, units="in")



# Nitrifiers -------------

rel_nit_sum <- rel_df %>% 
  filter(Genus %in% nitrifier_list) %>%
  group_by(date,Genus,battery) %>%
  summarise(sum=sum(Abundance))

rel_nit_compare <- rel_nit_sum %>%  filter(date %in% control_dates)

ggplot(data=rel_nit_sum, aes(x=date, y=sum, fill=Genus)) +
  facet_wrap(~battery, labeller=labels_basins) + 
  geom_bar(stat="identity") +
  theme_classic() + 
  scale_fill_manual(values=met.brewer("Egypt")) + 
  labs(y="Relative abundance [%]")
ggsave(filename=file.path("results","screening","nit_controltest_bar.tiff"), 
       height=2, width=4, units="in")

ggplot(data=rel_nit_compare, aes(x=battery, y=sum)) +
  facet_wrap(~Genus) + 
  geom_boxplot(outlier.color=NA) +
  geom_point() +
  labs(y="relative ab") +
  stat_compare_means(label="p.signif")
ggsave(filename=file.path("results","screening","nit_controltest.tiff"), 
       height=2, width=4, units="in")


rel_nit_sum %>%
  filter(Genus != "Nitrotoga") %>%
  pivot_wider(names_from=battery,values_from=sum) %>%
  ggplot(data=., aes(x=test, y=control, color=Genus)) +
  facet_wrap(~Genus, nrow=3) +
  geom_point() +
  theme_classic() +
  theme(legend.position="none") + 
  scale_color_manual(values=met.brewer("Egypt")) + 
  labs(y="Relative abundance in control basin [%]",
       x="Relative abundance in test basin [%]") +
  stat_cor(method="spearman")
ggsave(filename=file.path("results","nit_controltest_corr.tiff"), 
                                    height=4, width=3.5, units="in")

rel_nit_perf <- rel_nit_sum %>%
  pivot_wider(id_cols=c(date,battery), names_from=Genus, values_from=sum) %>%
  ungroup()

rel_nit_perf$AOB_NOB <- rel_nit_perf$Nitrospira / rel_nit_perf$Nitrosomonas

rel_nit_perf$AOB_NOB[is.infinite(rel_nit_perf$AOB_NOB)] <- NA

ggplot(rel_nit_perf, aes(x=date, y=AOB_NOB)) + 
  facet_wrap(~battery, labeller=labels_basins) + 
  geom_point() +
  stat_cor(method="spearman") +
  theme_bw() +
  scale_x_main + 
  labs(x="Date", y="Nitrospira:Nitrosomonas") +
  theme(axis.title.y=element_text(face="italic"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(filename=file.path("results","nitrifier_ratio.png"), 
       height=3, width=6, units="in", dpi=600)


