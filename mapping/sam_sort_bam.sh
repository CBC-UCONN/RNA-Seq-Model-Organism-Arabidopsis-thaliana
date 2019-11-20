#!/bin/bash
#SBATCH --job-name=sam_sort_bam
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

module load samtools

samtools view -@ 8 -bhS wt_Rep1.sam -o wt_Rep1.bam
samtools sort -@ 8 wt_Rep1.bam -o wt_Rep1_sort.bam

samtools view -@ 8 -bhS wt_Rep2.sam -o wt_Rep2.bam
samtools sort -@ 8 wt_Rep2.bam -o wt_Rep2_sort.bam

samtools view -@ 8 -bhS wt_Rep3.sam -o wt_Rep3.bam
samtools sort -@ 8 wt_Rep3.bam -o wt_Rep3_sort.bam

samtools view -@ 8 -bhS EE_Rep1.sam -o EE_Rep1.bam
samtools sort -@ 8 EE_Rep1.bam -o EE_Rep1_sort.bam

samtools view -@ 8 -bhS EE_Rep2.sam -o EE_Rep2.bam
samtools sort -@ 8 EE_Rep2.bam -o EE_Rep2_sort.bam

samtools view -@ 8 -bhS EE_Rep3.sam -o EE_Rep3.bam
samtools sort -@ 8 EE_Rep3.bam -o EE_Rep3_sort.bam

