#!/bin/bash
#SBATCH -D /x/sha1_checksum/
#SBATCH -A x.prj
#SBATCH -J sha1-x
#SBATCH -o sha1-x.out
#SBATCH -e sha1-x.err
#SBATCH -c 24
#SBATCH -p short

module load Python/3.9.6-GCCcore-11.2.0

python sha1_checksum.py
