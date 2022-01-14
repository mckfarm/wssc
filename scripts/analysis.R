### ---------
# WSSC Initial 16s rRNA sequence analysis
# McKenna Farmer  
# January 2022


### packages and working directory ---------
setwd("~/Github/wssc/")

library(qiime2R)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(lubridate)
library(MetBrewer)

### data import -----------
physeq <- qza_to_phyloseq(
  features="qiime/rarefied_table.qza",
  tree="qiime/rooted_tree.qza",
  taxonomy="qiime/taxonomy_silva.qza",
  metadata = "qiime/dec21_metadata.txt"
)

### summary stats ----------
# cleaning - remove mitochondria and chloroplasts
physeq <- subset_taxa(physeq, !Genus=="Mitochondria" & !Genus=="Chloroplast")

# relative abundance
rel <- transform_sample_counts(physeq, function(x) x*100/sum(x))
write.csv(rel@tax_table,"all_tax.csv")

# split between test batteries - used for the overall abundance comparison
rel_at1 <- subset_samples(rel, test=="AT1")
rel_at23 <- subset_samples(rel, test=="AT2/AT3")

topN <- 10
rel_10_at1_names <- sort(taxa_sums(rel_at1), decreasing=TRUE)[1:topN]
rel_10_at1 <- prune_taxa(names(rel_10_at1_names), rel_at1)

rel_10_at23_names <- sort(taxa_sums(rel_at23), decreasing=TRUE)[1:topN]
rel_10_at23 <- prune_taxa(names(rel_10_at23_names), rel_at23)

# specific genus of interest
## nitrifiers, PAOs/GAOs
rel_nit <- subset_taxa(rel,Genus=="Nitrotoga" | 
                           Genus=="Nitrospira" | 
                           Genus=="Nitrobacter" | 
                           Genus=="Nitrosomonas" | 
                           Genus=="Nitrosospira")
rel_phos <- subset_taxa(rel,Genus=="Candidatus_Accumulibacter" | 
                            Genus=="Candidatus_Competibacter" | 
                            Genus=="Tetrasphaera")

rel_nit <- tax_glom(rel_nit,"Genus")
rel_phos <- tax_glom(rel_phos,"Genus")

# prep for plotting with ggplot - convert to dataframe format
rel_nit_df <- psmelt(rel_nit) 
rel_phos_df <- psmelt(rel_phos)
rel_10_at1_df <- psmelt(rel_10_at1)
rel_10_at23_df <- psmelt(rel_10_at23)

# date formatting
rel_nit_df$date <- dmy(rel_nit_df$date)
rel_phos_df$date <- dmy(rel_phos_df$date)
rel_10_at1_df$date <- dmy(rel_10_at1_df$date)
rel_10_at23_df$date <- dmy(rel_10_at23_df$date)

### plotting ---------------
# nitrifiers
ggplot(data=rel_nit_df,mapping=aes(x=date,y=Abundance,fill=Genus)) + 
  facet_grid(test~.) + 
  geom_bar(stat="identity") +   
  theme_bw() + 
  ylim(0,7) + 
  scale_fill_manual(values=met.brewer("Isfahan2", 2)) + 
  theme(axis.text.x = element_text(angle = 0)) +
  scale_x_date(date_labels="%m-%y") +
  ylab("Relative Abundance (%)") + 
  xlab("Date (Month-Year)") +
  labs(fill="Genus")

ggsave("relative_ab_nit.tiff",width=2500,height=1500,unit="px")

# PAOs/GAOs
ggplot(data=rel_phos_df,mapping=aes(x=date,y=Abundance,
  fill=factor(Genus,levels=c("Candidatus_Accumulibacter","Tetrasphaera","Candidatus_Competibacter")))) + 
  facet_grid(test~.) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  ylim(0,7) + 
  scale_fill_manual(values=met.brewer("Lakota", 3)) + 
  theme(axis.text.x = element_text(angle = 0)) +
  scale_x_date(date_labels="%m-%y") +
  ylab("Relative Abundance (%)") + 
  xlab("Date (Month-Year)") + 
  labs(fill="Genus")
ggsave("relative_ab_phos.tiff",width=3000,height=1500,unit="px")


top10_at1 <- ggplot(data=rel_10_at1_df,mapping=aes(x=date,y=Abundance,fill=Genus)) +
  geom_bar(stat="identity") + 
  theme_bw() + 
  scale_fill_manual(values=met.brewer("Hiroshige", 10)) + 
  theme(axis.text.x = element_text(angle = 0)) +
  scale_x_date(date_labels="%m-%y") +
  ylab("Relative Abundance (%)") + 
  xlab("") + 
  labs(fill="Genus")

top10_at23 <- ggplot(data=rel_10_at23_df,mapping=aes(x=date,y=Abundance,fill=Genus)) +
  geom_bar(stat="identity") + 
  theme_bw() + 
  scale_fill_manual(values=met.brewer("Tiepolo", 10)) + 
  theme(axis.text.x = element_text(angle = 0)) +
  scale_x_date(date_labels="%m-%y") +
  ylab("Relative Abundance (%)") + 
  xlab("Date (Month-Year)") + 
  labs(fill="Genus")

plot_grid(top10_at1,top10_at23,nrow=2)
ggsave("top10.tiff",width=3000,height=1700,unit="px")

# shannon diversity
shannon <- plot_richness(physeq, x="test", measures=c("Shannon"))
shannon_df <- as.data.frame(shannon$data)
shannon_df$date <- dmy(shannon_df$date)

# used to format date - 
# https://stackoverflow.com/questions/5024798/how-can-a-color-gradient-based-on-date-be-applied-to-a-ggplot2-scatter-plot
# https://stackoverflow.com/questions/21311489/scatter-plot-with-ggplot2-colored-by-dates
as.Date_origin <- function(x){
  as.Date(x, origin = "1970-01-01")
}

ggplot(data=shannon_df) + 
  geom_point(aes(x=test, y=value, color=as.integer(date))) +
  geom_boxplot(aes(x=test, y=value), alpha=0.1) +
  theme_classic() + 
  ylab("Shannon Diversity Index") +
  xlab("Test battery") + 
  scale_colour_gradientn(colors=met.brewer("Lakota", 3), 
                        limits=as.integer(as.Date(c("2021-06-25","2021-12-27"))),
                        labels=as.Date_origin) + 
  labs(color="Date") +
  ylim(4.3,5.3)

ggsave("shannon.tiff",width=1200,height=1000,unit="px")
