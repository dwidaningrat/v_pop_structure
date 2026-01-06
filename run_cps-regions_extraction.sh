#!/bin/bash
#SBATCH -D /x/x
#SBATCH -A x.prj
#SBATCH -J cps-regions
#SBATCH -o cps-regions.out
#SBATCH -e cps-regions.err
#SBATCH -c 12
#SBATCH --gpus 1
#SBATCH -p gpu_short

module load Python/3.11.3-GCCcore-12.3.0
module load Biopython/1.83-foss-2023a

python3 cps-regions_extraction.py
