#!/bin/bash
#SBATCH --job-name=data_dump
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -o data_dump_%j.out
#SBATCH -e data_dump_%j.err
#SBATCH --partition=general
#SBATCH --qos=general


export TMPDIR=/home/CAM/$USER/tmp/

module load sratoolkit

fastq-dump --split-files SRR8428909 
mv SRR8428909_1.fastq wt_Rep1_R1.fastq
mv SRR8428909_2.fastq wt_Rep1_R2.fastq

fastq-dump --split-files SRR8428908
mv SRR8428908_1.fastq wt_Rep2_R1.fastq
mv SRR8428908_2.fastq wt_Rep2_R2.fastq
 
fastq-dump --split-files SRR8428907 
mv SRR8428907_1.fastq wt_Rep3_R1.fastq
mv SRR8428907_2.fastq wt_Rep3_R2.fastq

fastq-dump --split-files SRR8428906 
mv SRR8428906_1.fastq mutant_Rep1_R1.fastq
mv SRR8428906_2.fastq mutant_Rep1_R2.fastq

fastq-dump --split-files SRR8428905 
mv SRR8428905_1.fastq mutant_Rep2_R1.fastq
mv SRR8428905_2.fastq mutant_Rep2_R2.fastq

fastq-dump --split-files SRR8428904 
mv SRR8428904_1.fastq mutant_Rep3_R1.fastq
mv SRR8428904_2.fastq mutant_Rep3_R2.fastq


