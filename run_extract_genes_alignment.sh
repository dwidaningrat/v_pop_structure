#!/bin/bash
#SBATCH -D /x/x/
#SBATCH -A x.prj
#SBATCH -J extract-genes
#SBATCH -o extract-genes.out
#SBATCH -e extract-genes.err
#SBATCH -c 12
#SBATCH -p short

module load Python/3.11.3-GCCcore-12.3.0
module load Biopython/1.83-foss-2023a

python3 extract_genes_alignment.py
