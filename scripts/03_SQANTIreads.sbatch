#!/bin/bash
#SBATCH --array=1
#SBATCH --qos=short
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=4

#conda activate SQANTI3.env

annotation=/home/vgaya/estudios_in_silico/GRCm39.gtf
genome=/home/vgaya/estudios_in_silico/GRCm39.fa
design=/home/vgaya/estudios_in_silico/ONT/design.csv
QC_dir=/home/vgaya/estudios_in_silico/ONT/SQANTI_QC
reads_dir=/home/vgaya/estudios_in_silico/ONT/SQANTI_reads

python /home/vgaya/SQANTI3-5.3.0/sqanti_reads.py --design "$design" --annotation "$annotation" --genome "$genome" \
--force_id_ignore -i "$QC_dir" --dir "$reads_dir" --report both -p NEW
