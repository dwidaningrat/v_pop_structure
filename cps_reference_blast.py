import os
import subprocess
from Bio import SeqIO
from Bio.Blast import NCBIXML

#setup
GENOMES_DIR = "x/x"
QUERY_NUCLEOTIDES = "x.fasta"
OUTPUT_TSV = "x.tsv"
BLASTDB_DIR = "blastdbs"
TMP_DIR = "tmp"

os.makedirs(BLASTDB_DIR, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)

with open(OUTPUT_TSV, "w") as out_f:
    out_f.write("query_id\tgenome\tcontig\tstart\tend\tstrand\tscore\tevalue\tidentities\talign_length\tpercent_identity\tgaps\tquery_start\tquery_end\tbit_score\n")

queries = list(SeqIO.parse(QUERY_NUCLEOTIDES, "fasta"))

for genome_file in os.listdir(GENOMES_DIR):
    if not genome_file.endswith((".fna", ".fa", ".fasta")):
        continue

    genome_path = os.path.join(GENOMES_DIR, genome_file)
    genome_id = os.path.splitext(genome_file)[0]
    db_prefix = os.path.join(BLASTDB_DIR, genome_id)

    subprocess.run([
        "makeblastdb", "-in", genome_path, "-dbtype", "nucl", "-out", db_prefix
    ], check=True)

    for query in queries:
        query_file = os.path.join(TMP_DIR, f"{query.id}.fna")
        with open(query_file, "w") as f:
            SeqIO.write(query, f, "fasta")

        #Run blastn
        xml_output = os.path.join(TMP_DIR, f"{genome_id}_{query.id}_blastn.xml")
        subprocess.run([
            "blastn",
            "-query", query_file,
            "-db", db_prefix,
            "-evalue", "1e-5",
            "-outfmt", "5",
            "-max_target_seqs", "1",
            "-num_threads", "2",
            "-out", xml_output
        ], check=True)

        #Parse BLAST XML
        with open(xml_output) as result_handle:
            blast_records = list(NCBIXML.parse(result_handle))

        if not blast_records or not blast_records[0].alignments:
            continue

        alignment = blast_records[0].alignments[0]
        hsp = alignment.hsps[0]
        contig = alignment.hit_def
        start = min(hsp.sbjct_start, hsp.sbjct_end)
        end = max(hsp.sbjct_start, hsp.sbjct_end)
        strand = "+" if hsp.frame[1] > 0 else "-"
        score = hsp.score
        evalue = hsp.expect
        identities = hsp.identities
        align_length = hsp.align_length
        percent_identity = 100 * identities / align_length if align_length > 0 else 0
        gaps = hsp.gaps
        query_start = hsp.query_start
        query_end = hsp.query_end
        bit_score = hsp.bits

        #write output
        with open(OUTPUT_TSV, "a") as out_f:
            out_f.write(f"{query.id}\t{genome_id}\t{contig}\t{start}\t{end}\t{strand}\t{score}\t{evalue}\t{identities}\t{align_length}\t{percent_identity:.2f}\t{gaps}\t{query_start}\t{query_end}\t{bit_score}\n")
