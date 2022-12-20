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


library(ggpubr)

## plotting
library(ggplot2)
library(MetBrewer)
library(cowplot)
library(ggcorrplot)
# library(gridExtra)
# library(ggdist)
# library(ggConvexHull)

source(file.path("scripts","plotting.R"))

# performance data read in ------
lims <- readRDS(file.path("data","lims_data.rds"))
fit <- readRDS(file.path("data","fit_data.rds"))


# glycogen date read in ---

gly_dates <- read_rds(file=file.path("data","gly_dates.rds"))
gly_dates <- sort(gly_dates)

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

# all taxa ------
# parsing 
rel_top_list <- rel_df %>%
  group_by(Genus) %>%
  summarise(sum=sum(Abundance)) %>%
  top_n(10, sum)

rel_top <- rel_df %>%
  group_by(date, Genus) %>%
  summarise(sum=sum(Abundance)) %>%
  filter(Genus %in% rel_top_list$Genus) 

ggplot(data=rel_top, aes(x=date, y=sum, fill=Genus)) +
  geom_bar(stat="identity") + 
  theme_classic() +
  scale_x_main +
  scale_fill_manual(values=met.brewer("Signac", 10)) +
  labs(x="Date", y="Relative abundance [%]") +
  theme(legend.key.size = unit(.7,"line"),
        legend.text=element_text(size=6, face="italic"),
        legend.title=element_text(size=7),
        legend.position="right") +
  pause_lines +
  annotate(geom="text", x=ymd("2022-02-01"), y=50, label="Basin maintenance", size=2) +
  annotate("segment", x = ymd("2022-02-25"), y = 49, xend = ymd("2022-03-22"), yend = 49,
           arrow = arrow(type = "open", length = unit(0.02, "npc")), size=0.4)
ggsave(filename=file.path("results","top10_barplot.png"), 
       height=3.5, width=6, units="in", dpi=600)

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



# stats
rel_phos_wide <- rel_phos_sum %>%
  pivot_wider(id_cols=date, names_from=Genus, values_from=sum)
  
rel_phos_wide %>% filter(date <= phases$x1) %>%
  ggplot(data=., aes(x=date, y=Ca_Accumulibacter)) +
    geom_point() +
    stat_cor(method="spearman") +
    theme_classic()

rel_phos_wide %>% filter(date >= phases$x2) %>%
  ggplot(data=., aes(x=date, y=Ca_Accumulibacter)) +
   geom_point() +
   stat_cor(method="spearman") +
   theme_classic()



# plots
ggplot(data=rel_phos_sum, aes(x=date, y=sum, fill=Genus)) +
  geom_bar(stat="identity") + 
  theme_classic() +
  scale_fill_paogao + 
  scale_x_main +
  labs(x="Date", y="Relative abundance [%]") +
  theme(legend.key.size = unit(.7,"line"),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7),
        legend.position="right") +
  pause_lines +
  annotate(geom="text", x=ymd("2022-02-01"), y=15, label="Basin maintenance", size=2) +
  annotate("segment", x = ymd("2022-02-25"), y = 14.3, xend = ymd("2022-03-22"), yend = 14.3,
           arrow = arrow(type = "open", length = unit(0.02, "npc")), size=0.4)

ggsave(filename=file.path("results","paogao_barplot.png"), 
       height=3.5, width=6, units="in", dpi=600)


ggplot(data=rel_phos_sum, aes(x=Genus, y=sum, color=Genus)) +
  geom_boxplot(outlier.color=NA) + 
  geom_point(alpha=0.5) + 
  theme_classic() +
  coord_flip() + 
  scale_color_paogao +
  scale_x_paogao + 
  labs(x="Genus", y="Relative abundance [%]") +
  theme(legend.position="none")
ggsave(filename=file.path("results","paogao_boxplot.png"), 
       height=4, width=6, units="in", dpi=600)


# accumulibacter, tetrasphaera, and competibacter only
rel_phos_sum %>% filter(Genus %in% pao_gao_select) %>%
  ggplot(data=., aes(x=date, y=sum, color=Genus)) +
    geom_point(alpha=0.8) +
    geom_line(linetype="dotted") +
    theme_classic() +
    scale_color_paogao_select + 
    scale_x_main +
    ylim(0,8) +
    labs(x="Date", y="Relative abundance [%]") +
    theme(legend.position="top") +
    pause_lines

ggsave(filename=file.path("results","paogao_bygenus_test.png"), 
       height=3.5, width=5, units="in", dpi=600)

rel_phos_sum %>% 
  filter(Genus %in% pao_gao_select) %>%
  ggplot(., aes(x=date, y=sum, color=Genus)) +
  geom_point() +
  geom_line() + 
  scale_color_paogao_select +
  theme_classic() +
  theme(legend.position="top") +
  labs(y="Relative abundance [%]", x="Date")

ggsave(filename=file.path("results","screening","paogao_select_gly.tiff"), 
       height=5, width=6, units="in")

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
  labs(x="Date", y="Relative abundance [%]") +
  theme(legend.key.size = unit(.75,"line"),
        legend.text=element_text(size=8),
        legend.title=element_text(size=9),
        legend.position=c(0.13,0.88),
        legend.box.background = element_rect(colour = "black")) +
  pause_lines +
  annotate(geom="text", x=ymd("2022-02-01"), y=4.6, label="Basin\nmaintenance", size=2.5) +
  annotate("segment", x = ymd("2022-02-25"), y = 4.8, xend = ymd("2022-03-22"), yend = 4.8,
           arrow = arrow(type = "open", length = unit(0.02, "npc")), size=0.4)

ggsave(filename=file.path("results","nitrifier_barplot.png"), 
       height=3.5, width=5, units="in", dpi=600)



rel_nitrospira <- rel_df %>% 
  filter(Genus == "Nitrospira") %>%
  group_by(date,Species) %>%
  summarise(sum=sum(Abundance)) %>%
  replace_na(list(Species="Unknown")) %>%
  filter(! Species %in% c("midas_s_21815", "midas_s_22526"))

nitrospira_box <- ggplot(rel_nitrospira, aes(x=Species, y=sum, color=Species)) +
  geom_boxplot() +
  geom_point() + 
  theme_classic() + 
  scale_color_manual(values=met.brewer("Lakota",3)) +
  scale_x_discrete(labels=c("Nitrospira_defluvii"="Nitrospira defluvii", 
                            "Nitrospira_nitrosa" ="Nitrospira nitrosa", 
                            "Unknown" = "Unknown")) + 
  labs(y="Relative abundance [%]") +
  theme(legend.position="none") +
  ylim(0,4)

nitrospira_bar <- ggplot(rel_nitrospira, aes(x=date, y=sum, fill=Species)) +
  geom_bar(stat="identity") +
  theme_classic() + 
  scale_x_main + 
  scale_fill_manual(values=met.brewer("Lakota",3),
                    labels=c("Nitrospira defluvii", "Nitrospira nitrosa","Unknown")) + 
  labs(y="", x="Date") +
  theme(legend.position="none") +
  ylim(0,4)

plot_grid(nitrospira_box, nitrospira_bar, align="hv", axis="tblr", nrow=1)

ggsave(filename=file.path("results","nitrospira_combo.png"), 
       height=4, width=8, units="in", dpi=600)


# Performance comparison -------
# fit params -----
# airflow parsing
air <- fit %>% filter(LocationAbbrev == "AT 1" & ParamAbbrev == "Airflow") %>% 
  select(Date2, Value)

colnames(air) <- c("date","air") 
air$date <- as_date(air$date)

air <- air %>% filter(date >= date_range[1] & date <= date_range[2]) %>%
  filter(air >= 1000)


# svi parsing
svi <- fit %>% filter(LocationAbbrev == "AT 1" & ParamAbbrev == "SVI") %>% 
  select(Date2, Value)

colnames(svi) <- c("date","svi") 
svi$date <- as_date(svi$date)

svi <- svi %>% filter(date >= date_range[1] & date <= date_range[2])

# temp parsing
temp <- fit %>% filter(LocationAbbrev == "Raw" & ParamAbbrev =="Temp.") %>% select(Date2, Value)
colnames(temp) <- c("date","temp") 
temp$date <- as_date(temp$date)

temp <- temp %>% filter(date >= date_range[1] & date <= date_range[2]) %>%
  filter(temp >= 40)

fit_select <- left_join(air, temp) %>%
  left_join(., svi)

# making a df for correlation analysis
rel_phos_perf <- rel_phos_sum %>% 
  pivot_wider(id_cols=date, names_from=Genus, values_from=sum) %>%
  left_join(., fit_select, by="date") %>%
  ungroup()

rel_phos_corr <- rel_phos_perf %>% select(-date)

# nitrifier
rel_nit_perf <- rel_nit_sum %>%
  pivot_wider(id_cols=date, names_from=Genus, values_from=sum) %>%
  left_join(., fit_select, by="date") %>%
  ungroup()

rel_nit_perf$AOB_NOB <- rel_nit_perf$Nitrospira / rel_nit_perf$Nitrosomonas

rel_nit_perf$AOB_NOB[is.infinite(rel_nit_perf$AOB_NOB)] <- NA

rel_nit_corr <- rel_nit_perf %>% select(-date)

# filtering out columns with many zeros
rel_phos_corr <- rel_phos_corr[, colSums(rel_phos_corr != 0, na.rm=TRUE) > 10] 

phos_corr_mat <- cor(rel_phos_corr, method="spearman", use="complete.obs")
phos_corr_pmat <- cor_pmat(rel_phos_corr, method="spearman", use="complete.obs")

ggcorrplot(phos_corr_mat, p.mat=phos_corr_pmat, type="lower", lab=TRUE)
## no correlation with fit params 


rel_nit_corr <- rel_nit_corr[, colSums(rel_nit_corr != 0, na.rm=TRUE) > 10] 

nit_corr_mat <- cor(rel_nit_corr, method="spearman", use="complete.obs")
nit_corr_pmat <- cor_pmat(rel_nit_corr, method="spearman", use="complete.obs")

ggcorrplot(nit_corr_mat, p.mat=nit_corr_pmat, type="lower", lab=TRUE)




# influent effluent parsing ----
inf_eff <- lims %>% filter(`LIMS Sample` %in% c("Raw Influent + Recycle", "Secondary Effluent 1")) %>%
  select(c("Collect Date","Parameter Abbrev","LIMS Sample", "Result"))
  
colnames(inf_eff) <- c("date", "parameter", "location", "value")
inf_eff$date <- as_date(inf_eff$date)
inf_eff <- inf_eff[!duplicated(inf_eff), ] # remove duplicate rows

inf_eff_sel <- inf_eff %>% filter(date >= date_range[1] & date <= date_range[2]) %>%
  filter(location != "Secondary Effluent 1") %>%
  select(c("date","value","parameter")) %>%
  pivot_wider(names_from=parameter, values_from=value)


# making a df for correlation analysis
rel_phos_perf <- rel_phos_sum %>% 
  pivot_wider(id_cols=date, names_from=Genus, values_from=sum) %>%
  left_join(., inf_eff_sel, by="date") %>%
  ungroup()

rel_phos_corr <- rel_phos_perf %>% select(-date)

# filtering out columns with many zeros
rel_phos_corr <- rel_phos_corr[, colSums(rel_phos_corr != 0, na.rm=TRUE) > 10] 

phos_corr_mat <- cor(rel_phos_corr, method="spearman", use="complete.obs")
phos_corr_pmat <- cor_pmat(rel_phos_corr, method="spearman", use="complete.obs")

ggcorrplot(phos_corr_mat, p.mat=phos_corr_pmat, type="lower", lab=TRUE)
## no correlation with fit params 

ggplot(inf_eff_sel, aes(x=date,y=TSS)) +
  geom_point()


# nitrospira
rel_nitro_perf <- rel_nitrospira %>%
  pivot_wider(id_cols=c(date), names_from=Species, values_from=sum) %>%
  left_join(., fit_select, by="date") %>%
  ungroup()


rel_nitro_corr <- rel_nitro_perf %>% select(-date)

# filtering out columns with many zeros
rel_nitro_corr <- rel_nitro_corr[, colSums(rel_nitro_corr != 0, na.rm=TRUE) > 10] 

nitro_corr_mat <- cor(rel_nitro_corr, method="spearman", use="complete.obs")
nitro_corr_pmat <- cor_pmat(rel_nitro_corr, method="spearman", use="complete.obs")


ggcorrplot(nitro_corr_mat, p.mat=nitro_corr_pmat, type="lower", lab=TRUE)
