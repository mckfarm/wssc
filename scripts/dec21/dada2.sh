#!/bin/bash
#SBATCH --job-name="dada2"
#SBATCH -A b1042
#SBATCH -p genomicsguest
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --mem=30G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mckennafarmer2023@u.northwestern.edu
#SBATCH --output=dada2.out
#SBATCH --error=dada2.err

module purge all
module load qiime2/2021.11

qiime dada2 denoise-paired --verbose \
--i-demultiplexed-seqs /projects/b1052/mckenna/wssc/qiime/dec21_trimmed.qza \
--p-trunc-len-f 253 --p-trunc-len-r 200 \
--o-representative-sequences /projects/b1052/mckenna/wssc/qiime/rep_seqs_dada2.qza \
--o-table /projects/b1052/mckenna/wssc/qiime/table_dada2.qza --o-denoising-stats /projects/b1052/mckenna/wssc/qiime/stats_dada2.qza
