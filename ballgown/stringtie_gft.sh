#!/bin/bash
#SBATCH --job-name=stringtie
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=40G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --partition=general
#SBATCH --qos=general

export TMPDIR=/home/CAM/$USER/tmp/

#mkdir -p {athaliana_wt_Rep1,athaliana_wt_Rep2,athaliana_wt_Rep3,athaliana_EE_Rep1,athaliana_EE_Rep2,athaliana_EE_Rep3}

module load stringtie

stringtie -p 8 -l wT1 -G /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/TAIR10_GFF3_genes.gtf -o athaliana_wt_Rep1/transcripts.gtf ../mapping/wt_Rep1_sort.bam
stringtie -p 8 -l wT2 -G /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/TAIR10_GFF3_genes.gtf -o athaliana_wt_Rep2/transcripts.gtf ../mapping/wt_Rep2_sort.bam
stringtie -p 8 -l wT3 -G /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/TAIR10_GFF3_genes.gtf -o athaliana_wt_Rep3/transcripts.gtf ../mapping/wt_Rep3_sort.bam

stringtie -p 8 -l EE1 -G /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/TAIR10_GFF3_genes.gtf -o athaliana_EE_Rep1/transcripts.gtf ../mapping/EE_Rep1_sort.bam
stringtie -p 8 -l EE2 -G /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/TAIR10_GFF3_genes.gtf -o athaliana_EE_Rep2/transcripts.gtf ../mapping/EE_Rep2_sort.bam
stringtie -p 8 -l EE3 -G /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/TAIR10_GFF3_genes.gtf -o athaliana_EE_Rep3/transcripts.gtf ../mapping/EE_Rep3_sort.bam

ls -1 ath*/*.gtf >> sample_assembly_gtf_list.txt
stringtie --merge -p 8 -o stringtie_merged.gtf -G /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/TAIR10_GFF3_genes.gtf sample_assembly_gtf_list.txt

module load gffcompare/0.10.4

gffcompare -r /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/TAIR10_GFF3_genes.gtf -o gffcompare stringtie_merged.gtf


