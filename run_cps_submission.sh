#!/bin/bash

for i in {1..8}; do
    DIR="dir_$i"

    sbatch \
        -D "$(pwd)/$DIR" \
        -J "$DIR" \
        -o "$DIR.out" \
        -e "$DIR.err" \
        blast_job.sbatch
done