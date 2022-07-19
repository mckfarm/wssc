### ---------
# WSSC Initial 16s rRNA sequence analysis - MiDAS taxonomy assignments
# McKenna Farmer  
# January 2022

### packages and working directory ---------
setwd("~/Github/wssc/16S/dec21")

library(qiime2R)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(lubridate)
library(MetBrewer)
library(ggcorrplot)

### data import -----------
physeq <- qza_to_phyloseq(
  features="qiime_files/rarefied_table.qza",
  tree="qiime_files/rooted_tree.qza",
  taxonomy="qiime_files/taxonomy_midas.qza",
  metadata = "qiime_files/dec21_metadata.txt"
)

### summary stats ----------
# cleaning - remove mitochondria and chloroplasts
physeq <- subset_taxa(physeq, !Genus=="Mitochondria" & !Genus=="Chloroplast")

# relative abundance
rel <- transform_sample_counts(physeq, function(x) x*100/sum(x))
# write.csv(rel@tax_table,"all_tax.csv")

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
                           Genus=="Nitrosospira" |
                           Genus=="Nitrosococcus")
rel_phos <- subset_taxa(rel,Genus=="Ca_Accumulibacter" | 
                            Genus=="Tetrasphaera" |
                            Genus=="Dechloromonas" |
                            Genus=="Ca_Competibacter" |
                            Genus=="Micropruina")

rel_nit <- tax_glom(rel_nit,"Genus")
rel_phos <- tax_glom(rel_phos,"Genus")

# prep for plotting with ggplot - convert to dataframe format
rel_nit_df <- psmelt(rel_nit) 
rel_phos_df <- psmelt(rel_phos)
rel_10_at1_df <- psmelt(rel_10_at1)
rel_10_at23_df <- psmelt(rel_10_at23)
rel_df <- psmelt(rel)
write.csv(rel_df,"all_tax_table.csv")

# date formatting
rel_nit_df$date <- dmy(rel_nit_df$date)
rel_phos_df$date <- dmy(rel_phos_df$date)
rel_10_at1_df$date <- dmy(rel_10_at1_df$date)
rel_10_at23_df$date <- dmy(rel_10_at23_df$date)

### correlation matrix PAO/GAO ------

# subset columns
rel_phos_df_subset_at1 <- rel_phos_df %>% filter(test=="AT1") %>% select("Abundance","date","Genus","test")
rel_phos_df_subset_at23 <- rel_phos_df %>% filter(test=="AT2/AT3") %>% select("Abundance","date","Genus","test")

# reshape to wide format for correlation function
rel_phos_df_subset_at1 <- spread(rel_phos_df_subset_at1, key = Genus, value = Abundance) 
rel_phos_df_subset_at23 <- spread(rel_phos_df_subset_at23, key = Genus, value = Abundance) 
write.csv(rbind(rel_phos_df_subset_at1,rel_phos_df_subset_at23),"./results_midas/phos_taxa.csv")

# save correlation matrix
cor_rel_phos_at1 <- cor(select(rel_phos_df_subset_at1,-c(date:test)),method="spearman")
cor_rel_phos_at23 <- cor(select(rel_phos_df_subset_at23,-c(date:test)),method="spearman")

# save significance matrix
cor_rel_phos_at1_pmat <- cor_pmat(select(rel_phos_df_subset_at1,-c(date:test)),method="spearman")
cor_rel_phos_at23_pmat <- cor_pmat(select(rel_phos_df_subset_at23,-c(date:test)),method="spearman") 

### plotting 
# nitrifiers------------------
ggplot(data=rel_nit_df,mapping=aes(x=date,y=Abundance,fill=Genus)) + 
  facet_grid(test~.) + 
  geom_bar(stat="identity") +   
  theme_bw() + 
  scale_fill_manual(values=met.brewer("Isfahan2", 2)) + 
  theme(axis.text.x = element_text(angle = 0),legend.text = element_text(face = "italic")) +
  scale_x_date(date_labels="%m-%y") +
  ylab("Relative Abundance (%)") + 
  xlab("Date (Month-Year)") +
  labs(fill="Genus")
ggsave("./results_midas/relative_ab_nit.tiff",width=2500,height=1500,unit="px")

box_nit <- ggplot(data=rel_nit_df,mapping=aes(x=Genus,color=test,y=Abundance)) +
  geom_boxplot() +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values=met.brewer("Tiepolo", 5)) +
  scale_color_manual(values=c("darkred","navyblue")) +
  ylim(0,8) +
  labs(color="Battery",y="Relative Abundance [%]",x="") 
ggsave("./results_midas/nit_box.tiff",width=6,height=3,unit="in")

# PAOs/GAOs ----------------
ggplot(data=rel_phos_df,mapping=aes(x=date,y=Abundance,
  fill=factor(Genus,levels=c("Ca_Accumulibacter","Tetrasphaera","Dechloromonas","Ca_Competibacter","Micropruina")))) +
  facet_grid(test~.) +
  geom_bar(stat="identity") +
  theme_classic() +
  scale_fill_manual(values=met.brewer("Tiepolo", 5)) +
  theme(axis.text.x = element_text(angle = 0),legend.text = element_text(face = "italic")) +
  scale_x_date(date_labels="%b-%y") +
  ylab("Relative Abundance [%]") +
  xlab("Date") +
  labs(fill="Genus")
ggsave("./results_midas/relative_ab_phos.tiff",width=6,height=3,unit="in")

box_phos <- ggplot(data=rel_phos_df,mapping=aes(x=factor(Genus,levels=rev(c("Ca_Accumulibacter","Tetrasphaera","Dechloromonas","Ca_Competibacter","Micropruina"))),
                                    color=test,y=Abundance)) +
  geom_boxplot() +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values=met.brewer("Tiepolo", 5)) +
  scale_color_manual(values=c("darkred","navyblue")) +
  ylim(0,8) +
  labs(color="Battery",y="",x="") 
ggsave("./results_midas/phos_box.tiff",width=6,height=3,unit="in")

plot_grid(box_phos,box_nit,nrow=2,align="v")
ggsave("./results_midas/boxes.tiff",width=6,height=5,unit="in")

# top 10 abundant genus ------------
top10_at1 <- ggplot(data=rel_10_at1_df,mapping=aes(x=date,y=Abundance,fill=Genus)) +
  geom_bar(stat="identity") + 
  theme_bw() + 
  scale_fill_manual(values=met.brewer("Hiroshige", 10)) + 
  theme(axis.text.x = element_text(angle = 0),legend.text = element_text(face = "italic")) +
  scale_x_date(date_labels="%b-%y") +
  ylab("Relative Abundance [%]") + 
  xlab("") + 
  labs(fill="Genus") +
  ggtitle("Test - AT1")

top10_at23 <- ggplot(data=rel_10_at23_df,mapping=aes(x=date,y=Abundance,fill=Genus)) +
  geom_bar(stat="identity") + 
  theme_bw() + 
  scale_fill_manual(values=met.brewer("Tiepolo", 10)) + 
  theme(axis.text.x = element_text(angle = 0),legend.text = element_text(face = "italic")) +
  scale_x_date(date_labels="%b-%y") +
  ylab("Relative Abundance [%]") + 
  xlab("Date") + 
  labs(fill="Genus") +
  ggtitle("Control - AT2/AT3")

plot_grid(top10_at1,top10_at23,nrow=2)
ggsave("./results_midas/top10.tiff",width=10,height=7,unit="in")

# shannon diversity ------------
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

ggsave("./results_midas/shannon.tiff",width=1200,height=1000,unit="px")


# correlation plots ----------

ggcorrplot(cor_rel_phos_at1, p.mat=cor_rel_phos_at1_pmat, 
           type="lower", lab="true", insig="blank",colors=c("red","white","lightblue"))
ggsave("./results_midas/phos_at1.tiff",width=1200,height=1000,unit="px",scale=1.5)

ggcorrplot(cor_rel_phos_at23, p.mat=cor_rel_phos_at23_pmat, 
           type="lower", lab="true",insig="blank",colors=c("red","white","lightblue"))
ggsave("./results_midas/phos_at23.tiff",width=1200,height=1000,unit="px",scale=1.5)
