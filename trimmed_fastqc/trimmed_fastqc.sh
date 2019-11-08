#!/bin/bash
#SBATCH --job-name=quality_control
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --partition=general
#SBATCH --qos=general


export TMPDIR=/home/CAM/$USER/tmp/


module load fastqc
module load anaconda2/2.4.0
module load MultiQC



fastqc -t 4 ../trimmed_reads/trimmed_wt_Rep1_R1.fastq ../trimmed_reads/trimmed_wt_Rep1_R2.fastq

fastqc -t 4 ../trimmed_reads/trimmed_wt_Rep2_R1.fastq ../trimmed_reads/trimmed_wt_Rep2_R2.fastq

fastqc -t 4 ../trimmed_reads/trimmed_wt_Rep3_R1.fastq ../trimmed_reads/trimmed_wt_Rep3_R2.fastq

fastqc -t 4 ../trimmed_reads/trimmed_mutant_Rep1_R1.fastq ../trimmed_reads/trimmed_mutant_Rep1_R2.fastq

fastqc -t 4 ../trimmed_reads/trimmed_mutant_Rep2_R1.fastq ../trimmed_reads/trimmed_mutant_Rep2_R2.fastq

fastqc -t 4 ../trimmed_reads/trimmed_mutant_Rep3_R1.fastq ../trimmed_reads/trimmed_mutant_Rep3_R2.fastq

multiqc -n trimmed_fastqc .

