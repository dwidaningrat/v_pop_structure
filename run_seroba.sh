#!/bin/bash
#SBATCH -D /x/seroba/
#SBATCH -A x.prj
#SBATCH -J seroba
#SBATCH -o seroba.out    
#SBATCH -e seroba.err
#SBATCH -p long

module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate /x/.conda/envs/myproject

input='/x/fastq/'
db='/x/seroba/databases/database/'
output='/x/seroba/output/'

mkdir -p "$output"

for file in "$input"/*_1.f*q.gz; do
    if [[ -f $file ]]; then
        read2=${file/_1/_2}
        if [[ -f "$read2" ]]; then
            prefix=$(basename "$file" | sed 's/_1\.f.*q\.gz//')
            echo "Running SeroBA for $prefix..."
            seroba runSerotyping "$db" "$file" "$read2" "$output/$prefix" 2>> seroba_errors.log
        else
            echo "Paired file for $file not found!" >> seroba_errors.log
        fi
    fi
done

seroba summary "$output"
