#!/bin/bash
#SBATCH --array=1-4
#SBATCH --qos=short
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

# Get the fastq files from de accesion numbers download from GEO
cat SRR_Acc_List.txt | parallel -j 7 fastq-dump --split-files {} --outdir ./FASTQ
