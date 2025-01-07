#!/bin/bash
#SBATCH --job-name=flair        # Job name (showed with squeue)
#SBATCH --output=flair_%j.out   # Standard output and error log
#SBATCH --qos=short
#SBATCH --time=24:00:00
#SBATCH --mem=45G
#SBATCH --cpus-per-task=12

#conda activate flair

# Get the corresponding accession based on the array index
readarray -t a <<< $(ls pychopper_out/*.fastq)
reads=$(echo ${a[@]} | sed -s 's/ /,/g')

annotation=/home/vgaya/estudios_in_silico/GRCm39.gtf
genome=/home/vgaya/estudios_in_silico/GRCm39.fa
manifest=/home/vgaya/estudios_in_silico/ONT/reads_manifest.tsv

flair 1234 --reads $reads --genome $genome -f $annotation --gtf $annotation \
    --nvrna --temp_dir ./flair_out_2 \
    -o flair_out_2/align --threads 12 --reads_manifest $manifest
