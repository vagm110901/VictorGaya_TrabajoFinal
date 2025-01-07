#!/bin/bash
#SBATCH --array=1-4     # to parallelize the process with 4 array = 4 fastq files
#SBATCH --qos=short
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8

#mkdir -p ./pychopper_out

accession_list=SRR_Acc_List.txt

accession=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$accession_list")

# Generate fastq file name
file="./FASTQ/${accession}_1.fastq"

# Run chopper
pychopper -q 7 "$file" > ./pychopper_out/"${accession}"_1.fastq


