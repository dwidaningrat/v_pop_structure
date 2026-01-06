#!/bin/bash
#SBATCH -D /x/gtdb-tk/
#SBATCH -A x.prj
#SBATCH -J gtdbtk-x
#SBATCH -o gtdbtk-x.out
#SBATCH -e gtdbtk-x.err
#SBATCH -n 24
#SBATCH -p long

module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate /x/.conda/envs/gtdbtk-2.1.1

INPUT_DIR="/x/assemblies/"                                           
OUTPUT_DIR="/x/gtdb-tk/output/"                                                        
TEMP_DIR="/x/temp_output/"                                                     

mkdir -p ${OUTPUT_DIR}
mkdir -p ${TEMP_DIR}

gtdbtk classify_wf --genome_dir $INPUT_DIR --out_dir $OUTPUT_DIR -x fasta --prefix gtdbtk-x --cpus 24 --pplacer_cpus 24 --tmpdir
$TEMP_DIR --debug
