#!/bin/bash
#SBATCH --job-name="taxa"
#SBATCH -A b1042
#SBATCH -p genomicsguest
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --mem=30G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mckennafarmer2023@u.northwestern.edu
#SBATCH --output=taxa.out
#SBATCH --error=taxa.err

module purge all
module load qiime2/2021.2 # use older version of qiime since SILVA scikit-learn classifier is version 0.23.1

# midas classifier
# qiime feature-classifier classify-sklearn --i-classifier /projects/b1052/shared/qiime/midas_4.8.1_classifier.qza \
# --i-reads /projects/b1052/mckenna/wssc/qiime/rep_seqs_dada2.qza --o-classification /projects/b1052/mckenna/wssc/qiime/taxonomy.qza
#
# qiime metadata tabulate --m-input-file /projects/b1052/mckenna/wssc/qiime/taxonomy.qza \
# --o-visualization /projects/b1052/mckenna/wssc/qiime/taxonomy.qzv

# silva classifier
qiime feature-classifier classify-sklearn --i-classifier /projects/b1052/shared/qiime/silva-138-99-nb-classifier.qza \
--i-reads /projects/b1052/mckenna/wssc/qiime/rep_seqs_dada2.qza --o-classification /projects/b1052/mckenna/wssc/qiime/taxonomy_silva.qza

qiime metadata tabulate --m-input-file /projects/b1052/mckenna/wssc/qiime/taxonomy_silva.qza \
--o-visualization /projects/b1052/mckenna/wssc/qiime/taxonomy_silva.qzv
