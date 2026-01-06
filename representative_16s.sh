#!/bin/bash
#SBATCH -D /x/x/
#SBATCH -A x.prj
#SBATCH -J x
#SBATCH -o x.out
#SBATCH -e x.err
#SBATCH -c 12
#SBATCH --gpus 1
#SBATCH -p gpu_long

input="x.fasta"
output="x.fasta"

awk '
/^>/ {
    if (seq != "") {
        if (length(seq) > maxlen[acc]) {
            maxlen[acc] = length(seq)
            seqs[acc] = header "\n" seq
        }
    }
    header = $0
    acc = substr($0, 2)  # remove ">"
    seq = ""
    next
}
{
    seq = seq $0
}
END {
    if (length(seq) > maxlen[acc]) {
        seqs[acc] = header "\n" seq
    }
    for (a in seqs) {
        print seqs[a]
    }
}
' "$input" > "$output"
