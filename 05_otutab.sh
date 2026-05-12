#!/usr/bin/env bash
# Generate per-sample ZOTU abundance table by mapping pooled reads to ZOTUs.
# usearch -otutab assigns reads to samples via the "sample=" tag in read headers
# (added during pooling in 03_dereplicate.sh).
#
# Usage: bash 05_otutab.sh <zotus.fasta>
#
# Input:  pooled/all_samples.fastq  (pooled reads from 03_dereplicate.sh)
#         <zotus.fasta>             (ASV sequences from 04_unoise3.sh, e.g. pooled/zotus_minsize3.fasta)
# Output: results/zotu_table.txt    (samples x ZOTUs count matrix)

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: bash 05_otutab.sh <zotus.fasta>" >&2
    echo "  e.g. bash 05_otutab.sh pooled/zotus_minsize3.fasta" >&2
    exit 1
fi

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV_MAIN"

POOLED_DIR="pooled"
RESULTS_DIR="results"
LOG_DIR="logs/otutab"

mkdir -p "$RESULTS_DIR" "$LOG_DIR"

POOLED_FASTQ="$POOLED_DIR/all_samples.fastq"
ZOTUS_FASTA="$1"
# Name the output table after the input ZOTU file
ZOTU_TABLE="$RESULTS_DIR/zotu_table_$(basename "${ZOTUS_FASTA%.fasta}").txt"
LOG="$LOG_DIR/otutab.log"

THREADS="${THREADS:-$(nproc)}"

echo "=== Generating ZOTU table (threads: $THREADS) ==="

usearch \
    -otutab "$POOLED_FASTQ" \
    -zotus "$ZOTUS_FASTA" \
    -otutabout "$ZOTU_TABLE" \
    -threads "$THREADS" \
    > "$LOG" 2>&1

cat "$LOG"

echo ""
echo "ZOTU table written: $ZOTU_TABLE"
echo "Dimensions: $(awk 'NR==1{print NF-1, "samples x", NR-1, "ZOTUs"}' "$ZOTU_TABLE") (approx)"
