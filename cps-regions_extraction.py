import csv
import os
from Bio import SeqIO

#Setup
csv_file = "x.csv"
fasta_dir = "/x/x/"
gff_dir = "/x/x/"
output_fasta_dir = "x/all_fasta_outpu/"
output_gff_dir = "/x/all_gff_output/"

#Directories
os.makedirs(output_fasta_dir, exist_ok=True)
os.makedirs(output_gff_dir, exist_ok=True)

#Read and Process CSV
with open(csv_file, newline="") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        accession = row["accession"]
        contig = row["contig"]
        start = int(row["start"])
        end = int(row["end"])

        #Fasta Extraction
        fasta_path = os.path.join(fasta_dir, accession + ".fasta")
        out_fasta_path = os.path.join(output_fasta_dir, accession + ".fasta")
        sequence_found = False

        with open(fasta_path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id == contig:
                    extracted_seq = record.seq[start - 1:end]
                    new_record = record[start - 1:end]
                    new_record.id = accession
                    new_record.name = accession
                    new_record.description = f"{accession} region {contig}:{start}-{end}"
                    SeqIO.write(new_record, out_fasta_path, "fasta")
                    sequence_found = True
                    break

        if not sequence_found:
            print(f"[WARNING] Contig '{contig}' not found in {accession}.fasta")
            continue

        #GFF extraction
        gff_path = os.path.join(gff_dir, accession + ".gff")
        out_gff_path = os.path.join(output_gff_dir, accession + ".gff")

        with open(gff_path) as infile, open(out_gff_path, "w") as outfile:
            outfile.write("##gff-version 3\n")
            for line in infile:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) != 9:
                    continue  # skip malformed lines
                gff_contig, source, feature_type, gff_start, gff_end, score, strand, phase, attributes = parts
                gff_start = int(gff_start)
                gff_end = int(gff_end)

                if gff_contig == contig and gff_start >= start and gff_end <= end:
                    # Shift coordinates
                    new_start = gff_start - start + 1
                    new_end = gff_end - start + 1
                    outfile.write("\t".join([
                        accession, source, feature_type,
                        str(new_start), str(new_end),
                        score, strand, phase, attributes
                    ]) + "\n")
