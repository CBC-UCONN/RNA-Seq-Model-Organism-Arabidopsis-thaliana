#!/bin/bash
#SBATCH --job-name=hisat2
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=120G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --partition=xeon
#SBATCH --qos=general
export TMPDIR=/home/CAM/$USER/tmp/

module load hisat2


hisat2 -p 8 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/athaliana10/athaliana10 -1 ../trimmed_reads/trimmed_wt_Rep1_R1.fastq -2 ../trimmed_reads/trimmed_wt_Rep1_R2.fastq -S wt_Rep1.sam

hisat2 -p 8 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/athaliana10/athaliana10 -1 ../trimmed_reads/trimmed_wt_Rep2_R1.fastq -2 ../trimmed_reads/trimmed_wt_Rep2_R2.fastq -S wt_Rep2.sam

hisat2 -p 8 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/athaliana10/athaliana10 -1 ../trimmed_reads/trimmed_wt_Rep3_R1.fastq -2 ../trimmed_reads/trimmed_wt_Rep3_R2.fastq -S wt_Rep3.sam

hisat2 -p 8 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/athaliana10/athaliana10 -1 ../trimmed_reads/trimmed_mutant_Rep1_R1.fastq -2 ../trimmed_reads/trimmed_mutant_Rep1_R2.fastq -S mutant_Rep1.sam

hisat2 -p 8 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/athaliana10/athaliana10 -1 ../trimmed_reads/trimmed_mutant_Rep2_R1.fastq -2 ../trimmed_reads/trimmed_mutant_Rep2_R2.fastq -S mutant_Rep2.sam

hisat2 -p 8 --dta -x /isg/shared/databases/alignerIndex/plant/Arabidopsis/thaliana/athaliana10/athaliana10 -1 ../trimmed_reads/trimmed_mutant_Rep3_R1.fastq -2 ../trimmed_reads/trimmed_mutant_Rep3_R2.fastq -S mutant_Rep3.sam


