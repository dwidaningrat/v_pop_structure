#!/bin/bash
#SBATCH -D /x/snp-sites/
#SBATCH -A x.prj
#SBATCH -J snp-sites
#SBATCH -o snp-sites.out 
#SBATCH -e snp-sites.err 
#SBATCH -c 24
#SBATCH -p long

module load snp-sites/2.5.1-GCCcore-10.3.0

output_aln="/x/snp-sites/output/x.aln"
input_aln="/x/xalignment.fasta"

snp-sites -m -v -p -o $output_aln $input_aln
