#!/usr/bin/env bash
# Pool all trimmed barcode reads, then dereplicate the combined dataset.
# This ensures ASVs (ZOTUs) are resolved consistently across all samples.
#
# Workflow:
#   1. Pool all barcode fastq.gz into a single uncompressed FASTQ
#      (read headers are prefixed with the barcode name for traceability)
#   2. Dereplicate with usearch -fastx_uniques -sizeout
#      -minuniquesize 2 : discard singletons
#      -sizeout         : add ;size=N abundance annotations for unoise3
#
# Output:
#   pooled/all_samples.fastq       — pooled reads (intermediate)
#   pooled/all_samples_derep.fasta — dereplicated sequences, ready for unoise3

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV_MAIN"

INPUT_DIR="fastq_trimmed"
POOLED_DIR="pooled"
LOG_DIR="logs/derep"

mkdir -p "$POOLED_DIR" "$LOG_DIR"

POOLED_FASTQ="$POOLED_DIR/all_samples.fastq"
DEREP_FASTA="$POOLED_DIR/all_samples_derep.fasta"
LOG="$LOG_DIR/derep_pooled.log"

# --- Step 1: Pool all barcodes ---
echo "=== Pooling all barcodes ==="
rm -f "$POOLED_FASTQ"

for input_file in "$INPUT_DIR"/barcode*.fastq.gz; do
    barcode=$(basename "$input_file" .fastq.gz)
    echo "  Adding $barcode..."
    # Tag each read header with sample=barcode (required by usearch -otutab)
    zcat "$input_file" | awk -v bc="$barcode" '
        NR%4==1 { print "@" bc "." substr($0,2) ";sample=" bc }
        NR%4!=1 { print }
    ' >> "$POOLED_FASTQ"
done

total_reads=$(awk 'NR%4==1' "$POOLED_FASTQ" | wc -l)
echo "  Total pooled reads: $total_reads"

# --- Step 2: Dereplicate pooled reads ---
echo ""
echo "=== Dereplicating pooled reads ==="

usearch \
    -fastx_uniques "$POOLED_FASTQ" \
    -fastaout "$DEREP_FASTA" \
    -minuniquesize 2 \
    -sizeout \
    > "$LOG" 2>&1

cat "$LOG"

n_uniques=$(grep -c '^>' "$DEREP_FASTA")
echo ""
echo "Unique sequences (size >= 2): $n_uniques"
echo "Done. Output: $DEREP_FASTA"
