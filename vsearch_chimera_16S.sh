python3 - <<'PY'
from collections import defaultdict
from Bio import SeqIO
f="16S_fixed-header.fasta"
d=defaultdict(list)
for r in SeqIO.parse(f,"fasta"):
    acc = r.id.split("|")[0]
    d[acc].append(r)
for acc, seqs in d.items():
    if len(seqs) > 1:
        with open(f"{acc}_copies.fasta","w") as fh:
            for i,s in enumerate(seqs):
                s.id = f"{acc}|copy{i+1}"
                s.description=""
                SeqIO.write(s,fh,"fasta")
PY
