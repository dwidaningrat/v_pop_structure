#!/bin/bash
#SBATCH -A x.prj
#SBATCH -D /x/x/16S/
#SBATCH -J pairwise-16S
#SBATCH -o pairwise-16S.out
#SBATCH -e pairwise-16S.err
#SBATCH -c 12
#SBATCH --gpus 1
#SBATCH -p gpu_long

module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate /x/.conda/envs/myproject

set -euo pipefail

FASTA="16S.fasta"
OUTDIR="/x/16S/ouput"
PAIRWISE_DIR="/x/16S/pairwise_vsearch"
SUMMARY="/x/16S/pairwise_summary.tsv"
THREADS=12

mkdir -p "$OUTDIR" "$PAIRWISE_DIR"

echo "Step 1 — build accession -> sequence groups and count copies..." >&2

export FASTA OUTDIR

python3 - <<'PY'
from collections import defaultdict
from Bio import SeqIO
import os, sys, re

fasta = os.environ.get("FASTA", "16S.fasta")
outdir = os.environ.get("OUTDIR", "output")

groups = defaultdict(list)

def get_accession(header):
    # header: full header string (without leading '>' and removed '|' or whitespace
    h = header.strip()
    first = h.split('|',1)[0].split()[0]
    return first

def safe_filename(name):
    # keep alphanum, dot, underscore, dash and replace others with underscore
    return re.sub(r'[^A-Za-z0-9._-]', '_', name)

for r in SeqIO.parse(fasta, "fasta"):
    hdr = r.description if r.description else r.id
    acc = get_accession(hdr)
    groups[acc].append((r.id, str(r.seq), r.description))

multi_accs = [acc for acc, seqs in groups.items() if len(seqs) > 1]
print(f"Found {len(groups)} distinct accessions total, {len(multi_accs)} with >1 copy", file=sys.stderr)

os.makedirs(outdir, exist_ok=True)
for acc in multi_accs:
    safe = safe_filename(acc) + ".fasta"
    outpath = os.path.join(outdir, safe)
    with open(outpath, "w") as fh:
        for i,(rid, seq, desc) in enumerate(groups[acc]):
            header = desc if desc and desc != rid else f"{acc}|copy{i+1}"
            fh.write(f">{header}\n")
            for j in range(0, len(seq), 80):
                fh.write(seq[j:j+80] + "\n")

print(f"Wrote {len(multi_accs)} per-accession FASTA files to {outdir}", file=sys.stderr)
PY

echo -e "accession\tn_copies\tmin_id\tmax_id\tmean_id\tvsearch_userout" > "$SUMMARY"

echo "Step 2 — run vsearch --allpairs_global on each multi-copy FASTA and summarize..." >&2

shopt -s nullglob
for f in "${OUTDIR}"/*.fasta; do
    [ -e "$f" ] || continue
    base=$(basename "$f" .fasta)
    userout="${PAIRWISE_DIR}/${base}.vsearch.tsv"
    echo "Running vsearch on $base ..." >&2

    # run vsearch
    vsearch --allpairs_global "$f" \
           --userout "$userout" \
           --userfields query+target+id+alnlen+mism+opens \
	   --acceptall \
           --iddef 1 \
           --threads "$THREADS" \
           2> "${PAIRWISE_DIR}/${base}.vsearch.log"

    # Count pair lines (skip empty/comment lines)
    n_pairs=$(awk 'NF && $1 !~ /^#/ {count++} END{print count+0}' "$userout" || echo 0)
    n_seqs=$(grep -c "^>" "$f" || true)

    if [ "$n_pairs" -eq 0 ]; then
        min_id="NA"; max_id="NA"; mean_id="NA"
    else
        # compute min, max, mean
        min_id=$(awk 'NF && $1 !~ /^#/ { if(min==""||$3+0<min) min=$3+0; n++ } END{ if(n>0) printf "%.6f", min; else print "NA"}' "$userout")
        max_id=$(awk 'NF && $1 !~ /^#/ { if(max==""||$3+0>max) max=$3+0; n++ } END{ if(n>0) printf "%.6f", max; else print "NA"}' "$userout")
        mean_id=$(awk 'NF && $1 !~ /^#/ { sum+=$3+0; n++ } END{ if(n>0) printf "%.6f", sum/n; else print "NA"}' "$userout")
    fi

    echo -e "${base}\t${n_seqs}\t${min_id}\t${max_id}\t${mean_id}\t${userout}" >> "$SUMMARY"
done

echo "Done. Summary at $SUMMARY" >&2
echo "Pairwise outputs at $PAIRWISE_DIR/ and per-accession FASTAs at $OUTDIR/" >&2
