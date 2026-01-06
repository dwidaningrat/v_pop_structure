#!/bin/bash
#SBATCH -D /x/x
#SBATCH -A x.prj
#SBATCH -J x_clonalframe
#SBATCH -o x_clonalframe.out
#SBATCH -e x_clonalframe.err
#SBATCH --ntasks=1
#SBATCH --mem=800G
#SBATCH --partition=himem

module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate myproject

treefile="/x/x/x.treefile"
alignment="/x/x/x_renamedtip.fasta"
prefix="x"

ClonalFrameML $treefile $alignment $prefix
