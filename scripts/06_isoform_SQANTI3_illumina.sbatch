#!/bin/bash
#SBATCH --job-name=isoform_sqanti3_illumina        # Job name (showed with squeue)
#SBATCH --output=isoform_sqanti3_illumina_%j.out   # Standard output and error log
#SBATCH --qos=short
#SBATCH --time=24:00:00
#SBATCH --mem=90G
#SBATCH --cpus-per-task=8

annotation=/home/vgaya/estudios_in_silico/GRCm39.gtf
genome=/home/vgaya/estudios_in_silico/GRCm39.fa
shortreads=/home/vgaya/estudios_in_silico/shortreads.fofn


python /home/vgaya/SQANTI3-5.3.0/sqanti3_qc.py ./flair_out/align.isoforms.gtf \
    $annotation $genome -d ./isoform_sqanti3_illumina -o isoform_illumina --force_id_ignore \
    --report pdf --short_reads $shortreads
