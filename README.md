# WSSC
Analysis scripts for WSSC collaboration - 16s rRNA sequence analysis

## Programs and computing resources:  
- 16s rRNA amplicon sequence analysis using QIIME2 performed on Northwestern Quest computing cluster
- Data analysis using R locally

## QIIME2 workflow:  
1) create manifest file
- created in Excel and saved to manifest_dec21.csv

2.0) start an interactive job in Quest
```
srun --account=b1042 --time=01:00:00 --partition=genomicsguest --mem=1G --pty bash -l
module load qiime2/2021.11
```
2) import paired end reads - these will be demultiplexed since they are imported with the manifest file
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

5) [dada2.sh]()
- denoise and trim
- trimming is based on read quality statistics from the previous step - keep sequences with average read quality of >30
- this command creates three files: dada2 quality filtering table (stats), data table of read info that can be coupled to metadata (table), and a list of amplicon sequence variants that will be used for blast or other commands (rep_seqs)



--o-representative-sequences /projects/b1052/mckenna/wssc/qiime/rep_seqs_dada2.qza \
--o-table /projects/b1052/mckenna/wssc/qiime/table_dada2.qza --o-denoising-stats /projects/b1052/mckenna/wssc/qiime/stats_dada2.qza


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

7) phylogenetic tree with mafft
```
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences /projects/b1052/mckenna/wssc/qiime/rep_seqs_dada2.qza \
--o-alignment /projects/b1052/mckenna/wssc/qiime/aligned_rep_seqs_dada2.qza \
--o-masked-alignment /projects/b1052/mckenna/wssc/qiime/masked_aligned_rep_seqs_dada2.qza \
--o-tree /projects/b1052/mckenna/wssc/qiime/unrooted_tree.qza \
--o-rooted-tree /projects/b1052/mckenna/wssc/qiime/rooted_tree.qza
```

8) [taxa.sh]()
- assigns taxa from Midas classifier
- also produces output qzv file for viewing

9) alpha rarefaction curves
- measure of how diversity changes with sequencing depth
```
qiime diversity alpha-rarefaction \
--i-table /projects/b1052/mckenna/wssc/qiime/table_dada2.qza \
--i-phylogeny /projects/b1052/mckenna/wssc/qiime/rooted_tree.qza \
--o-visualization /projects/b1052/mckenna/wssc/qiime/alpha_rarefaction.qzv \
--p-max-depth 20000
```

10) rarefaction
- picking a depth of 8000 based on rough estimate of plateau in faith_pd rarefaction curve
```
qiime diversity core-metrics-phylogenetic \
--i-table /projects/b1052/mckenna/wssc/qiime/table_dada2.qza \
--i-phylogeny /projects/b1052/mckenna/wssc/qiime/rooted_tree.qza \
--p-sampling-depth 8000 \
--m-metadata-file /projects/b1052/mckenna/wssc/qiime/dec21_metadata.txt \
--output-dir /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/core-metrics-results-8000
```

13) [analysis_notebook.Rmd](https://github.com/mckfarm/s2ebpr_16s/blob/main/analysis_notebook.Rmd)
- Data analysis performed in R with phyloseq and other R functions
