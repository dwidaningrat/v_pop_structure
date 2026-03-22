#!/bin/bash
#SBATCH -D /x/x/ #working directory
#SBATCH -A x.prj 
#SBATCH -J x
#SBATCH -o x.out
#SBATCH -e x.err
#SBATCH -c x #number of cores
#SBATCH -p long #partition

module load Anaconda3/2022.05
eval "$(conda shell.bash hook)"
conda activate myproject

mkdir -p results_pirate

INPUT_DIRECTORY="/x/x/"
OUTPUT_DIRECTORY="/x/x/"

PIRATE -i $INPUT_DIRECTORY -o $OUTPUT_DIRECTORY -a -r -t 24
