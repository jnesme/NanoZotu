#!/bin/bash
### General options
#BSUB -q hpcspecial
#BSUB -J otutab
#BSUB -n 20
#BSUB -R "span[hosts=1] rusage[mem=1GB]"
#BSUB -M 1000MB
#BSUB -W 1:00
#BSUB -u jnesme@gmail.com
#BSUB -B
#BSUB -N
#BSUB -o logs/otutab_%J.out
#BSUB -e logs/otutab_%J.err

# [OPTIONAL — QC step] Generate per-sample ZOTU abundance table by mapping
# pooled reads to ZOTUs. Not required for the NanoASV branch (steps 08–10).
#
# Primary uses:
#   1. QC: the read mapping rate here should agree with the fraction of reads
#      classified by NanoASV — a large discrepancy indicates a pipeline issue.
#   2. minsize sweep: run with different ZOTU FASTAs to compare how threshold
#      choice affects abundance profiles before committing to MINSIZE_WORKING.
#
# usearch -otutab assigns reads to samples via the "sample=" tag in read headers
# (added during pooling in 03_dereplicate.sh).
#
# Defaults to pooled/zotus_minsize${MINSIZE_WORKING}.fasta from config.sh.
# Override by passing a ZOTU FASTA as argument (interactive use only):
#   bash 05_otutab.sh pooled/zotus_minsize5.fasta
#
# Input:  pooled/all_samples.fastq
#         pooled/zotus_minsize${MINSIZE_WORKING}.fasta  (or argument)
# Output: results/zotu_table_zotus_minsize${MINSIZE_WORKING}.txt
#
# Usage: bsub < 05_otutab.sh

set -euo pipefail

# PROJECT_DIR must be an absolute path — LSF copies this script to /tmp before
# execution, so relative paths and ${BASH_SOURCE[0]} are not reliable.
PROJECT_DIR="/work3/josne/github/NanoZotu"
source "$PROJECT_DIR/config.sh"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV_MAIN"

THREADS=$LSB_DJOB_NUMPROC

POOLED_DIR="$PROJECT_DIR/pooled"
RESULTS_DIR="$PROJECT_DIR/results"
LOG_DIR="$PROJECT_DIR/logs/otutab"

mkdir -p "$RESULTS_DIR" "$LOG_DIR"

POOLED_FASTQ="$POOLED_DIR/all_samples.fastq"
ZOTUS_FASTA="${1:-$POOLED_DIR/zotus_minsize${MINSIZE_WORKING}.fasta}"
ZOTU_TABLE="$RESULTS_DIR/zotu_table_$(basename "${ZOTUS_FASTA%.fasta}").txt"
LOG="$LOG_DIR/otutab.log"

echo "=== Generating ZOTU table ==="
echo "  ZOTUs:   $ZOTUS_FASTA"
echo "  Threads: $THREADS"

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
