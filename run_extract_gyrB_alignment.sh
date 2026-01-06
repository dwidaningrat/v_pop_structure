#!/bin/bash
#SBATCH -D /x/gyrB/
#SBATCH -A x.prj
#SBATCH -J extract-gyrB
#SBATCH -o extract-gyrB.out
#SBATCH -e extract-gyrB.err
#SBATCH -c 12
#SBATCH -p short

module load Python/3.11.3-GCCcore-12.3.0
module load Biopython/1.83-foss-2023a

python3 extract_gyrB_alignment.py
