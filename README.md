# WSSC low DO nutrient removal project
Analysis workflow for 16S rRNA sequence analysis

## Programs and computing resources:  
- 16S rRNA amplicon sequence analysis using QIIME2 2021.11 performed on Northwestern Quest computing cluster
- Use an interactive job and activate the module in Quest - module load qiime2/2021.11
- Data analysis using R locally

## QIIME2 workflow:  
1) create manifest file for each sequence run
- December 21 and November 22 sequencing runs at Rush U Microbiome Sequencing Center
- create manifest files in Excel and save as csv in Quest where raw reads are located

2) import paired end reads - these will be output as demultiplexed since they are imported with the manifest file

```
qiime tools import --input-path /projects/p31629/wssc/raw_reads/manifest_dec21.csv \
--input-format PairedEndFastqManifestPhred33 \
--output-path /projects/p31629/wssc/qiime/dec21.qza \
--type SampleData[PairedEndSequencesWithQuality]

qiime tools import --input-path /projects/p31629/wssc/raw_reads/manifest_nov22.csv \
--input-format PairedEndFastqManifestPhred33 \
--output-path /projects/p31629/wssc/qiime/nov22.qza \
--type SampleData[PairedEndSequencesWithQuality]
```

3) trim primers and visualize read quality

```
qiime cutadapt trim-paired --i-demultiplexed-sequences dec21.qza \
--p-front-f GTGYCAGCMGCCGCGGTAA \
--p-front-r CCGYCAATTYMTTTRAGTTT \
--p-error-rate 0.1 \
--p-overlap 3 \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--o-trimmed-sequences dec21_trimmed.qza

qiime cutadapt trim-paired --i-demultiplexed-sequences nov22.qza \
--p-front-f GTGYCAGCMGCCGCGGTAA \
--p-front-r CCGYCAATTYMTTTRAGTTT \
--p-error-rate 0.1 \
--p-overlap 3 \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--o-trimmed-sequences nov22_trimmed.qza

# error rate and overlap are default parameters - written out in the command for informational purposes

qiime demux summarize --i-data dec21_trimmed.qza \
--o-visualization dec21_trimmed_readquality.qzv

qiime demux summarize --i-data nov22_trimmed.qza \
--o-visualization nov22_trimmed_readquality.qzv

```

4) produce ASVs (dada2.sh)
- denoise, then merge sequences from separate runs before identifying ASVs
- trimming is based on read quality statistics from the previous step - keep sequences with average read quality of >30
- this command creates three files: dada2 quality filtering table (stats), data table of read info that can be coupled to metadata (table), and a list of amplicon sequence variants that will be used for blast or other commands (rep_seqs)


5) merge sequencing runs and produce qzv files
```
qiime feature-table merge \
  --i-tables table_dada2_dec21.qza \
  --i-tables table_dada2_nov22.qza \
  --o-merged-table table_dada2_all.qza

qiime feature-table merge-seqs \
  --i-data rep_seqs_dada2_dec21.qza \
  --i-data rep_seqs_dada2_nov22.qza  \
  --o-merged-data rep_seqs_dada2_all.qza

qiime feature-table tabulate-seqs \
  --i-data rep_seqs_dada2_all.qza \
  --o-visualization rep_seqs_dada2_all.qzv

```


6) create phylogenetic tree with mafft

```
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep_seqs_dada2_all.qza \
--o-alignment aligned_rep_seqs_dada2.qza \
--o-masked-alignment masked_aligned_rep_seqs_dada2.qza \
--o-tree unrooted_tree.qza \
--o-rooted-tree rooted_tree.qza
```


7) assign taxonomy (taxa.sh)
- assigns taxa from Midas classifier


8) make alpha rarefaction curves
- measure of how diversity changes with sequencing depth

```
qiime diversity alpha-rarefaction \
--i-table table_dada2_all.qza \
--i-phylogeny rooted_tree.qza \
--o-visualization alpha_rarefaction.qzv \
--p-max-depth 11000
```


## Data analysis
Check out the data analysis script 
