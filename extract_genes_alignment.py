#!/usr/bin/env python3

from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

gene_positions = {
    "pflA": (126631, 127503),
    "ppaC": (127504, 128439),
    "sodA": (201502, 202113),
    "rpoB": (250657, 254382),
    "tuf": (394624, 395820),
    "guaA": (413872, 415482),
    "map": (486682, 487542),
    "pyk": (523168, 524724)
}

# Load core gene alignment
alignment_file = "/x/core_alignment.fasta"
alignment = AlignIO.read(alignment_file, "fasta")

# Extract sequences for each gene and save in a new alignment file
output_file = "x/multiple_x_alignment.fasta"
with open(output_file, "w") as out_f:
    for record in alignment:
        extracted_seq = ""
        for gene, (start, end) in gene_positions.items():
            extracted_seq += str(record.seq[start-1:end])
        new_record = SeqRecord(Seq(extracted_seq), id=record.id, description="")
        out_f.write(f">{new_record.id}\n{new_record.seq}\n")

print(f"Extracted housekeeping gene alignments saved in {output_file}")
