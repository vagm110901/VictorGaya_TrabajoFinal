#!/bin/bash
#SBATCH --array=1-4     # to parallelize the process with 4 array = 4 fastq files
#SBATCH --qos=short
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

#conda activate SQANTI3.env

accession_list=SRR_Acc_List.txt

accession=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$accession_list")

# Generate fastq file name
file="./FASTQ/${accession}_1.fastq"

python /home/vgaya/SQANTI3-5.3.0/sqanti3_qc.py "$file" \
/home/vgaya/estudios_in_silico/GRCm39.gtf /home/vgaya/estudios_in_silico/GRCm39.fa --fasta \
--force_id_ignore --dir /home/vgaya/estudios_in_silico/ONT/SQANTI_QC_NEW -t 20 --report both



