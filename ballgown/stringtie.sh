#!/bin/bash
#SBATCH --job-name=stringtie
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=120G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --partition=general
#SBATCH --qos=general

export TMPDIR=/home/CAM/$USER/tmp/

mkdir -p {athaliana_wt_Rep1,athaliana_wt_Rep2,athaliana_wt_Rep3,athaliana_mutant_Rep1,athaliana_mutant_Rep2,athaliana_mutant_Rep3}

module load stringtie

stringtie -e -B -p 8 ../mapping/wt_Rep1_sort.bam -G /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/TAIR10_GFF3_genes.gtf -o athaliana_wt_Rep1/athaliana_wt_Rep1.count -A athaliana_wt_Rep1/wt_Rep1_gene_abun.out

stringtie -e -B -p 8 ../mapping/wt_Rep2_sort.bam -G /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/TAIR10_GFF3_genes.gtf -o athaliana_wt_Rep2/athaliana_wt_Rep2.count -A athaliana_wt_Rep2/wt_Rep2_gene_abun.out

stringtie -e -B -p 8 ../mapping/wt_Rep3_sort.bam -G /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/TAIR10_GFF3_genes.gtf -o athaliana_wt_Rep3/athaliana_wt_Rep3.count -A athaliana_wt_Rep3/wt_Rep3_gene_abun.out

stringtie -e -B -p 8 ../mapping/mutant_Rep1_sort.bam -G /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/TAIR10_GFF3_genes.gtf -o athaliana_mutant_Rep1/athaliana_mutant_Rep1.count -A athaliana_mutant_Rep1/mutant_Rep1_gene_abun.out

stringtie -e -B -p 8 ../mapping/mutant_Rep2_sort.bam -G /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/TAIR10_GFF3_genes.gtf -o athaliana_mutant_Rep2/athaliana_mutant_Rep2.count -A athaliana_mutant_Rep2/mutant_Rep2_gene_abun.out

stringtie -e -B -p 8 ../mapping/mutant_Rep3_sort.bam -G /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/TAIR10_GFF3_genes.gtf -o athaliana_mutant_Rep3/athaliana_mutant_Rep3.count -A athaliana_mutant_Rep3/mutant_Rep3_gene_abun.out

