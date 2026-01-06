#!/usr/bin/env python3

from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

gene_positions = {
    "gyrB": (308548, 310500)
}

# Load the concatenated core gene alignment
alignment_file = "/x/core_alignment.fasta"
alignment = AlignIO.read(alignment_file, "fasta")

# Extract sequences for each housekeeping gene and save in a new alignment file
output_file = "/x/gyrB/gyrB.fasta"
with open(output_file, "w") as out_f:
    for record in alignment:
        extracted_seq = ""
        for gene, (start, end) in gene_positions.items():
            extracted_seq += str(record.seq[start-1:end])  # Ensure 1-based index is correctly adjusted
        new_record = SeqRecord(Seq(extracted_seq), id=record.id, description="")
        out_f.write(f">{new_record.id}\n{new_record.seq}\n")

print(f"Extracted housekeeping gene alignments saved in {output_file}")
