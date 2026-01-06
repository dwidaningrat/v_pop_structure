#!/bin/bash
#SBATCH -D /x/x
#SBATCH -A x.prj
#SBATCH -J iqtree-x
#SBATCH -o iqtree-x.out
#SBATCH -e iqtree-x.err
#SBATCH -c 32
#SBATCH --mem=800G
#SBATCH --partition=himem

module load IQ-TREE/2.3.5-gompi-2023a

alignment="/x/x/x.fasta"
prefix="x"

iqtree2 -s $alignment \
--prefix $prefix \
-m MFP \
-bb 1000 \
-nt 32
