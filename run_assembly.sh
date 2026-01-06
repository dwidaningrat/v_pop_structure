#!/bin/bash
#SBATCH -o assembly_x.out
#SBATCH -e assembly_x.err
#SBATCH -P x.prjc
#SBATCH -N assembly_x
#SBATCH -q long
#SBATCH -l h_vmem=16G

#run this script for every batch, modify the data and output directories, run in a separate tmux session

export NXF_ANSI_LOG=false
export NXF_OPTS="-Xms16G -Xmx16G"

NEXTFLOW_WORKFLOWS_DIR='/x/nextflow_workflows'
DATA_DIR='/x/batch_x/fastq/' #modify this directory for every batch
OUTPUT_DIR='/x/batch_x/assembly/' #modify this directory for every batch

/x/software/bin/nextflow run \
${NEXTFLOW_WORKFLOWS_DIR}/assembly-dev/main.nf \
--adapter_file ${NEXTFLOW_WORKFLOWS_DIR}/assembly-dev/adapters.fas \
--qc_conditions ${NEXTFLOW_WORKFLOWS_DIR}/assembly-dev/qc_conditions_nextera_relaxed_dev.yml \
--input_dir ${DATA_DIR} \
--fastq_pattern '*{R,_}{1,2}.f*q.gz' \
--output_dir ${OUTPUT_DIR} \
--depth_cutoff 100 \
--prescreen_file_size_check 12 \
--confindr_db_path /x/software/confindr_database \
--careful \
--skip_fastqc_multiqc \
--skip_quast_summary \
--skip_quast_multiqc \
--kmer_min_copy 3 \
-w ${OUTPUT_DIR}/work \
-profile bmrc -qs 1000 -resume

# clean up on exit 0
status=$?
## take some decision ##
if [[ $status -eq 0 ]]; then
  rm -r ${OUTPUT_DIR}/work
fi
