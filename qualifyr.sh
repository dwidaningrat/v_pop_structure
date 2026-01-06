#!/bin/bash
#SBATCH -A x.prj
#SBATCH -J qualifyr_x
#SBATCH -o qualifyr_x.out
#SBATCH -e qualifyr_x.err
#SBATCH -p short

#modify the script and run for every batch

module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate myproject

batch="/x/x/batch_x/"
output_dir="/x/x/qualifyr/"

qualifyr report -i $batch \
-o $output_dir \
-t  qualifyr_report_batch_x
