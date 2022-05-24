# WSSC
Analysis scripts for WSSC collaboration - 16s rRNA sequence analysis

## Programs and computing resources:  
- 16s rRNA amplicon sequence analysis using QIIME2 2021.11 performed on Northwestern Quest computing cluster
- Use an interactive job and activate the module in Quest - module load qiime2/2021.11
- Data analysis using R locally

## QIIME2 workflow:  
1) create manifest file
- created in Excel and saved to dec21_manifest.csv

2) import paired end reads - these will be output as demultiplexed since they are imported with the manifest file

```
qiime tools import --input-path /projects/b1052/mckenna/wssc/raw_reads/dec21_manifest.csv \
--input-format PairedEndFastqManifestPhred33 \
--output-path /projects/b1052/mckenna/wssc/qiime/dec21.qza \
--type SampleData[PairedEndSequencesWithQuality]
```

3) trim primers

```
qiime tools import --input-path /projects/b1052/mckenna/wssc/raw_reads/manifest_dec21.csv \
--input-format PairedEndFastqManifestPhred33 \
--output-path /projects/b1052/mckenna/wssc/qiime/dec21.qza \

qiime cutadapt trim-paired --i-demultiplexed-sequences /projects/b1052/mckenna/wssc/qiime/dec21.qza  \
--p-front-f GTGYCAGCMGCCGCGGTAA \
--p-front-r CCGYCAATTYMTTTRAGTTT \
--p-error-rate 0.1 \
--p-overlap 3 \
--o-trimmed-sequences /projects/b1052/mckenna/wssc/qiime/dec21_trimmed.qza
```
error rate and overlap are default parameters - written out in the command for informational purposes

4) visualize read quality

```
qiime demux summarize --i-data /projects/b1052/mckenna/wssc/qiime/dec21.qza \
--o-visualization /projects/b1052/mckenna/wssc/qiime/dec21_readquality.qzv

qiime demux summarize --i-data /projects/b1052/mckenna/wssc/qiime/dec21_trimmed.qza \
--o-visualization /projects/b1052/mckenna/wssc/qiime/dec21_trimmed_readquality.qzv
```

5) [dada2.sh](https://github.com/mckfarm/wssc/blob/main/16S/scripts/dec21/dada2.sh)
- denoise and trim
- trimming is based on read quality statistics from the previous step - keep sequences with average read quality of >30
- this command creates three files: dada2 quality filtering table (stats), data table of read info that can be coupled to metadata (table), and a list of amplicon sequence variants that will be used for blast or other commands (rep_seqs)


6) process dada2 outputs
- List ASVs (to keep in case I want to blast later)
- Show dada2 denoising stats

```
qiime feature-table tabulate-seqs \
--i-data /projects/b1052/mckenna/wssc/qiime/rep_seqs_dada2.qza \
--o-visualization /projects/b1052/mckenna/wssc/qiime/rep_seqs_dada2.qzv

qiime metadata tabulate --m-input-file /projects/b1052/mckenna/wssc/qiime/stats_dada2.qza \
--o-visualization /projects/b1052/mckenna/wssc/qiime/stats_dada2.qzv
```

7) create phylogenetic tree with mafft

```
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences /projects/b1052/mckenna/wssc/qiime/rep_seqs_dada2.qza \
--o-alignment /projects/b1052/mckenna/wssc/qiime/aligned_rep_seqs_dada2.qza \
--o-masked-alignment /projects/b1052/mckenna/wssc/qiime/masked_aligned_rep_seqs_dada2.qza \
--o-tree /projects/b1052/mckenna/wssc/qiime/unrooted_tree.qza \
--o-rooted-tree /projects/b1052/mckenna/wssc/qiime/rooted_tree.qza
```

8) assign taxonomy with [taxa.sh](https://github.com/mckfarm/wssc/blob/main/16S/scripts/dec21/taxa.sh)
- assigns taxa from Midas and Silva classifers
- also produces output qzv file for viewing

9) make alpha rarefaction curves
- measure of how diversity changes with sequencing depth

```
qiime diversity alpha-rarefaction \
--i-table /projects/b1052/mckenna/wssc/qiime/table_dada2.qza \
--i-phylogeny /projects/b1052/mckenna/wssc/qiime/rooted_tree.qza \
--o-visualization /projects/b1052/mckenna/wssc/qiime/alpha_rarefaction.qzv \
--p-max-depth 11000
```

10) rarefy samples
- picking a depth of 5000 based on rough estimate of plateau in faith_pd rarefaction curve
- rarefaction is optional - you randomly select data from your sequences at the same depth for each sample, this is best practice to compare sequencing data from different sampling dates and DNA extractions but you are losing data in the process 

```
qiime diversity core-metrics-phylogenetic \
--i-table /projects/b1052/mckenna/wssc/qiime/table_dada2.qza \
--i-phylogeny /projects/b1052/mckenna/wssc/qiime/rooted_tree.qza \
--p-sampling-depth 5000 \
--m-metadata-file /projects/b1052/mckenna/wssc/qiime/dec21_metadata.txt \
--output-dir /projects/b1052/mckenna/wssc/qiime/core-metrics-results-5000
```

11) Comparing classifier outputs - this actually didnt work because the classifiers use different labels for kingdom, keeping for reference

```
qiime feature-table relative-frequency \
--i-table /projects/b1052/mckenna/wssc/qiime/table_dada2.qza \
--o-relative-frequency-table /projects/b1052/mckenna/wssc/qiime/relative_table_dada2.qza

qiime quality-control evaluate-taxonomy \
--i-expected-taxa /projects/b1052/mckenna/wssc/qiime/taxonomy_silva.qza \
--i-observed-taxa /projects/b1052/mckenna/wssc/qiime/taxonomy_midas.qza \
--i-feature-table /projects/b1052/mckenna/wssc/qiime/relative_table_dada2.qza \
--p-depth 6 \
--o-visualization /projects/b1052/mckenna/wssc/qiime/taxa_compare.qzv
```

# Data analysis
[analysis_midas.R](https://github.com/mckfarm/wssc/blob/main/16S/scripts/dec21/analysis_midas.R)
- Data analysis performed in R with phyloseq and other R functions
