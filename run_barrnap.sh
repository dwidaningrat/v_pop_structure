#!/bin/bash
#SBATCH -D /x/x/16S/
#SBATCH -A x.prj
#SBATCH -J 16S-barrnap
#SBATCH -o 16S-barrnap.out
#SBATCH -e 16S-barrnap.err
#SBATCH -c 12
#SBATCH --gpus 1
#SBATCH -p gpu_long

module load Python/3.11.3-GCCcore-12.3.0
module load Biopython/1.83-foss-2023a
module load HMMER/3.4-gompi-2023a
module load BEDTools/2.31.0-GCC-12.3.0
module load Perl/5.36.1-GCCcore-12.3.0

module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate /x/.conda/envs/myproject

#setup
INPUT_DIR="/x/x"
OUTPUT_DIR="/x/x/16S_outputs"
mkdir -p "$OUTPUT_DIR"

for f in "$INPUT_DIR"/*.fasta; do
    base=$(basename "$f" .fasta)

    echo "Processing $base..."

    barrnap "$f" > "$OUTPUT_DIR/${base}.gff"

    grep "16S_rRNA" "$OUTPUT_DIR/${base}.gff" > "$OUTPUT_DIR/${base}_16s.gff"

    awk -v id="$base" 'BEGIN{OFS="\t"} !/^#/ {
        contig=$1; start=$4 - 1; end=$5; strand=$7;
        header = id "|" contig "|" $4 "_" $5;
        print contig, start, end, header, ".", strand
    }' "$OUTPUT_DIR/${base}_16s.gff" > "$OUTPUT_DIR/${base}_16s.bed"

    bedtools getfasta -fi "$f" -bed "$OUTPUT_DIR/${base}_16s.bed" -s -name \
        -fo "$OUTPUT_DIR/${base}_16s.fasta"
done

# all FASTA files with full headers
cat "$OUTPUT_DIR"/*_16S.fasta > "$OUTPUT_DIR/16S.fasta"

# all FASTA files only accession ID as header
awk '/^>/ {
    split($0, a, "|");
    print ">" a[1]
    next
}1' "$OUTPUT_DIR/16S.fasta" > "$OUTPUT_DIR/16S_fixed-header.fasta"

echo "✅ Done:"
echo "• Full headers:      $OUTPUT_DIR/16S.fasta"
echo "• Accession only:    $OUTPUT_DIR/16S_fixed-header.fasta"
