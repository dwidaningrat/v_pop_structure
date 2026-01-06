#!/bin/bash
#SBATCH -A x.prj
#SBATCH -p short
module load prokka/1.14.5-gompi-2021a
file=$1
output=$2
mkdir annotations/$output
prokka --outdir annotations/$output --force --prefix $output $file
