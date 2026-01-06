#!/bin/bash
#SBATCH -D /x/x
#SBATCH -A xn.prj
#SBATCH -J mafft-x
#SBATCH -o mafft-x.out
#SBATCH -e mafft-x.err
#SBATCH -c 12
#SBATCH --gpus 1
#SBATCH -p gpu_long

module load MAFFT/7.520-GCC-12.3.0-with-extensions

mafft --auto concatenated_x.fasta > aligned_x.fasta
