# WSSC
Analysis scripts for WSSC collaboration - 16s rRNA sequence analysis

## Programs and computing resources:  
- 16s rRNA amplicon sequence analysis using QIIME2 performed on Northwestern Quest computing cluster
- Data analysis using R locally

## QIIME2 workflow:  
1) create manifest file

2) import paired end reads
```
singularity exec -B /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s /projects/b1052/shared/qiime2-core2020.11.simg \
qiime tools import --input-path /projects/b1052/mckenna/ \
--input-format PairedEndFastqManifestPhred33 \
--output-path /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/round1_paired.qza \
--type SampleData[PairedEndSequencesWithQuality]
```

3) visualize read quality
```
singularity exec -B /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s /projects/b1052/shared/qiime2-core2020.11.simg \
qiime demux summarize --i-data /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/round1_paired.qza \
--o-visualization /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/visuals/round1_paired.qzv
```

4) [dada2.sh](https://github.com/mckfarm/s2ebpr_16s/blob/main/dada2.sh)
- denoise and trim
- this command creates three files: dada2 quality filtering table (stats), data table of read info that can be coupled to metadata (table), and a list of amplicon sequence variants that will be used for blast or other commands (rep_seqs)
- also includes merging the three separate sequencing runs

5) mapping
- couples metadata to sequences
```
singularity exec -B /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s /projects/b1052/shared/qiime2-core2020.11.simg \
qiime feature-table summarize \
--i-table /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/table_merged.qza \
--o-visualization /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/table_merged.qzv \
--m-sample-metadata-file /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/metadata.txt
```

6) list ASVs
```
singularity exec -B /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s /projects/b1052/shared/qiime2-core2020.11.simg \
qiime feature-table tabulate-seqs \
--i-data /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/rep_seqs_merged.qza \
--o-visualization /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/rep_seqs_merged.qzv
```

7) phylogenetic tree with mafft
```
singularity exec -B /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s /projects/b1052/shared/qiime2-core2020.11.simg \
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/rep_seqs_merged.qza \
--o-alignment /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/aligned_rep_seqs_merged.qza \
--o-masked-alignment /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/masked_aligned_rep_seqs_merged.qza \
--o-tree /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/unrooted_tree.qza \
--o-rooted-tree /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/rooted_tree.qza
```

7) make [Midas classifier](https://midasfieldguide.org/guide/downloads)
- Midas is a 16s sequence database specifically for wastewater
```
singularity exec -B /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s -B /projects/b1052/shared /projects/b1052/shared/qiime2-core2020.11.simg \
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path /projects/b1052/shared/midas_3.7.fa \
--output-path /projects/b1052/shared/midas_asv.qza
singularity exec -B /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s -B /projects/b1052/shared /projects/b1052/shared/qiime2-core2020.11.simg \
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path /projects/b1052/shared/midas_3.7.txt \
--output-path /projects/b1052/shared/midas_taxonomy.qza
```

8) [taxa.sh](https://github.com/mckfarm/s2ebpr_16s/blob/main/taxa.sh)
- assigns taxa from Midas classifier
- also produces output qzv file for viewing

9) alpha rarefaction curves
- measure of how diversity changes with sequencing depth
- switched to new qiime filepath!
```
singularity exec -B /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s -B /projects/b1052/shared/qiime /projects/b1052/shared/qiime/qiime2-core2020.11.simg \
qiime diversity alpha-rarefaction \
--i-table /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/table_merged.qza \
--i-phylogeny /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/rooted_tree.qza \
--o-visualization /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/alpha_rarefaction.qzv \
--p-max-depth 15000
```

10) rarefaction
- picking a depth of 8000 based on rough estimate of plateau in faith_pd rarefaction curve
```
singularity exec -B /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s -B /projects/b1052/shared/qiime /projects/b1052/shared/qiime/qiime2-core2020.11.simg \
qiime diversity core-metrics-phylogenetic \
--i-phylogeny /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/rooted_tree.qza \
--i-table /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/table_merged.qza \
--p-sampling-depth 8000 \
--m-metadata-file /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/metadata.txt \
--output-dir /projects/b1052/Wells_b1042/McKenna/s2ebpr_16s/core-metrics-results-8000
```

13) [analysis_notebook.Rmd](https://github.com/mckfarm/s2ebpr_16s/blob/main/analysis_notebook.Rmd)
- Data analysis performed in R with phyloseq and other R functions
