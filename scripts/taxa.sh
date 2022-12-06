#!/bin/bash
#SBATCH --job-name="taxa"
#SBATCH -A b1042
#SBATCH -p genomicsguest
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH --mem=30G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mckennafarmer2023@u.northwestern.edu
#SBATCH --output=taxa.out
#SBATCH --error=taxa.err

module purge all
module load qiime2/2021.11 

# midas classifier

qiime feature-classifier classify-sklearn \
    --i-classifier /projects/b1052/shared/qiime/midas_4.8.1_classifier.qza \
    --i-reads /projects/p31629/wssc/qiime/rep_seqs_dada2_all.qza \
    --o-classification /projects/p31629/wssc/qiime/taxonomy.qza

qiime metadata tabulate \
    --m-input-file /projects/p31629/wssc/qiime/taxonomy.qza \
    --o-visualization /projects/p31629/wssc/qiime/taxonomy.qzv
