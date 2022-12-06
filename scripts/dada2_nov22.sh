#!/bin/bash
#SBATCH --job-name="dada2_nov22"
#SBATCH -A b1042
#SBATCH -p genomicsguest
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH --mem=30G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mckennafarmer2023@u.northwestern.edu
#SBATCH --output=dada2_nov22.out
#SBATCH --error=dada2_nov22.err

module purge all
module load qiime2/2021.11

qiime dada2 denoise-paired --verbose \
--i-demultiplexed-seqs /projects/p31629/wssc/qiime/nov22_trimmed.qza \
--p-trunc-len-f 280 --p-trunc-len-r 200 \
--o-representative-sequences /projects/p31629/wssc/qiime/rep_seqs_dada2_nov22.qza \
--o-table /projects/p31629/wssc/qiime/table_dada2_nov22.qza \
--o-denoising-stats /projects/p31629/wssc/qiime/stats_dada2_nov22.qza

qiime metadata tabulate \
  --m-input-file /projects/p31629/wssc/qiime/stats_dada2_nov22.qza \
  --o-visualization /projects/p31629/wssc/qiime/stats_dada2_nov22.qzv