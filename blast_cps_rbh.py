import os
import subprocess
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq

QUERY_PROTEINS = "x/cps_genes_cdhit_98.faa"
GENOMES_DIR = "x/"
OUTPUT_TSV = os.path.basename(os.getcwd()) + ".tsv"
BLASTDB_DIR = "blastdbs"
TMP_DIR = "tmp"

EVALUE_CUTOFF = 1e-5
MIN_QUERY_COVERAGE = 0.70
THREADS = 2

os.makedirs(BLASTDB_DIR, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)

# Load queries
queries = list(SeqIO.parse(QUERY_PROTEINS, "fasta"))

# Write output header
with open(OUTPUT_TSV, "w") as out_f:
    out_f.write(
        "query_id\tgenome\tcontig\tstart\tend\tstrand\t"
        "score\tevalue\tidentities\talign_length\t"
        "percent_identity\tgaps\tquery_start\tquery_end\tbit_score\n"
    )


QUERY_DB = os.path.join(BLASTDB_DIR, "query_db")

subprocess.run([
    "makeblastdb",
    "-in", QUERY_PROTEINS,
    "-dbtype", "prot",
    "-out", QUERY_DB
], check=True)


for genome_file in os.listdir(GENOMES_DIR):
    if not genome_file.endswith((".fa", ".fna", ".fasta")):
        continue

    genome_path = os.path.join(GENOMES_DIR, genome_file)
    genome_id = os.path.splitext(genome_file)[0]
    db_prefix = os.path.join(BLASTDB_DIR, genome_id)

    # Build genome BLAST DB
    subprocess.run([
        "makeblastdb",
        "-in", genome_path,
        "-dbtype", "nucl",
        "-out", db_prefix
    ], check=True)

    # Load genome sequences (for extraction)
    genome_seqs = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))

    for query in queries:
        query_faa = os.path.join(TMP_DIR, f"{query.id}.faa")
        SeqIO.write(query, query_faa, "fasta")

        tblastn_xml = os.path.join(
            TMP_DIR, f"{genome_id}_{query.id}_tblastn.xml"
        )

        subprocess.run([
            "tblastn",
            "-query", query_faa,
            "-db", db_prefix,
            "-evalue", str(EVALUE_CUTOFF),
            "-outfmt", "5",
            "-max_target_seqs", "5",
            "-num_threads", str(THREADS),
            "-out", tblastn_xml
        ], check=True)

        with open(tblastn_xml) as handle:
            records = list(NCBIXML.parse(handle))

        if not records or not records[0].alignments:
            continue

        alignment = records[0].alignments[0]
        hsp = alignment.hsps[0]

        query_len = len(query.seq)
        query_coverage = hsp.align_length / query_len

        if query_coverage < MIN_QUERY_COVERAGE:
            continue

        contig = alignment.hit_def
        start = min(hsp.sbjct_start, hsp.sbjct_end)
        end = max(hsp.sbjct_start, hsp.sbjct_end)
        strand = "+" if hsp.frame[1] > 0 else "-"

        hit_seq = genome_seqs[contig].seq[start - 1:end]
        if strand == "-":
            hit_seq = hit_seq.reverse_complement()

        hit_fna = os.path.join(
            TMP_DIR, f"{genome_id}_{query.id}_hit.fna"
        )

        with open(hit_fna, "w") as f:
            f.write(f">{genome_id}_{query.id}\n{hit_seq}\n")

        rbh_xml = os.path.join(
            TMP_DIR, f"{genome_id}_{query.id}_rbh.xml"
        )

        subprocess.run([
            "blastx",
            "-query", hit_fna,
            "-db", QUERY_DB,
            "-evalue", str(EVALUE_CUTOFF),
            "-outfmt", "5",
            "-max_target_seqs", "1",
            "-num_threads", str(THREADS),
            "-out", rbh_xml
        ], check=True)

        with open(rbh_xml) as handle:
            rbh_records = list(NCBIXML.parse(handle))

        if not rbh_records or not rbh_records[0].alignments:
            continue

        best_back_hit = rbh_records[0].alignments[0].hit_def

        if best_back_hit != query.id:
            continue

        percent_identity = 100 * hsp.identities / hsp.align_length

        with open(OUTPUT_TSV, "a") as out_f:
            out_f.write(
                f"{query.id}\t{genome_id}\t{contig}\t"
                f"{start}\t{end}\t{strand}\t"
                f"{hsp.score}\t{hsp.expect}\t"
                f"{hsp.identities}\t{hsp.align_length}\t"
                f"{percent_identity:.2f}\t{hsp.gaps}\t"
                f"{hsp.query_start}\t{hsp.query_end}\t"
                f"{hsp.bits}\n"
            )