#!/usr/bin/env python3

from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Define genes and their positions (1-based indexing)
gene_positions = {
    "geneA": (x, x),
    "geneB": (x, x),
    "geneC": (x, x),
    "geneD": (x, x),
    "geneE": (x, x),
    "geneF": (x, x),
    "geneG": (x, x),
    "geneH": (x, x)
}

# Load core gene alignment
alignment_file = "/x/core_alignment.fasta"
alignment = AlignIO.read(alignment_file, "fasta")

# Extract sequences for each gene and save in a new alignment file
output_file = "/x/x/multiple_x_alignment.fasta"
with open(output_file, "w") as out_f:
    for record in alignment:
        extracted_seq = ""
        for gene, (start, end) in gene_positions.items():
            extracted_seq += str(record.seq[start-1:end])
        new_record = SeqRecord(Seq(extracted_seq), id=record.id, description="")
        out_f.write(f">{new_record.id}\n{new_record.seq}\n")

print(f"Extracted housekeeping gene alignments saved in {output_file}")
