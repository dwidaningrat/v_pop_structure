#!/bin/bash
#SBATCH -D /x/x
#SBATCH -A aanensen.prj
#SBATCH -J reference_cps_blast
#SBATCH -o reference_cps_blast.out
#SBATCH -e reference_cps_blast.err
#SBATCH -c 12
#SBATCH --gpus 1
#SBATCH -p gpu_long

module load Python/3.11.3-GCCcore-12.3.0
module load Biopython/1.83-foss-2023a
module load BLAST+/2.14.1-gompi-2023a

python3 cps_reference_blast.py
