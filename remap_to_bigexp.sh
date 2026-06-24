#!/bin/bash
### General options
#BSUB -q hpcspecial
#BSUB -J remap_bigexp
#BSUB -n 20
#BSUB -R "span[hosts=1] rusage[mem=1024]"
#BSUB -W 1:00
#BSUB -u jnesme@gmail.com
#BSUB -B
#BSUB -N
#BSUB -o logs/remap_bigexp_%J.out
#BSUB -e logs/remap_bigexp_%J.err

# Remap Prelim16S reads against BigExp ZOTUs to produce comparable
# per-sample abundance tables across datasets.
#
# BigExp's ZOTUs (higher coverage, stricter minsize) serve as the reference
# set. Mapping Prelim reads against them bypasses Prelim's UNOISE3
# fragmentation — reads land on the correct BigExp ZOTU regardless of how
# many near-identical Prelim ZOTUs exist.
#
# Runs two mappings: against BigExp minsize5 and minsize8.
#
# Input:  pooled/all_samples.fastq  (Prelim16S, with sample= tags)
#         BigExp ZOTUs (paths below)
# Output: results/remap_bigexp/
#           zotu_table_bigexp_minsize5.txt
#           zotu_table_bigexp_minsize8.txt
#           mapping_stats.txt
#
# Usage: bsub < remap_to_bigexp.sh

set -euo pipefail

source "/work3/josne/github/NanoZotu/config.sh"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV_MAIN"

THREADS=$LSB_DJOB_NUMPROC

# ─── Configuration ────────────────────────────────────────────────────────────

BIGEXP_DIR="/work3/josne/Projects/Igalbana_Microbiome_GRF/Igalbana_Microbiome_16S/NanoZotu"

POOLED_FASTQ="$PROJECT_DIR/pooled/all_samples.fastq"
OUT_DIR="$PROJECT_DIR/results/remap_bigexp"
LOG_DIR="$PROJECT_DIR/logs/remap_bigexp"

BIGEXP_MINSIZES=(5 8)

# ─── Validate ─────────────────────────────────────────────────────────────────

[ -f "$POOLED_FASTQ" ] || { echo "ERROR: not found: $POOLED_FASTQ"; exit 1; }

for ms in "${BIGEXP_MINSIZES[@]}"; do
    f="$BIGEXP_DIR/pooled/zotus_minsize${ms}.fasta"
    [ -f "$f" ] || { echo "ERROR: not found: $f"; exit 1; }
done

mkdir -p "$OUT_DIR" "$LOG_DIR"

> "$OUT_DIR/mapping_stats.txt"

# ─── Map against each BigExp ZOTU set ────────────────────────────────────────

for ms in "${BIGEXP_MINSIZES[@]}"; do
    ZOTUS_FASTA="$BIGEXP_DIR/pooled/zotus_minsize${ms}.fasta"
    OTUTAB="$OUT_DIR/zotu_table_bigexp_minsize${ms}.txt"
    LOG="$LOG_DIR/otutab_bigexp_minsize${ms}.log"

    echo "=== Mapping Prelim reads -> BigExp minsize${ms} ==="
    echo "  Reference: $ZOTUS_FASTA ($(grep -c '^>' "$ZOTUS_FASTA") ZOTUs)"
    echo "  Reads:     $POOLED_FASTQ"
    echo "  Threads:   $THREADS"

    usearch \
        -otutab "$POOLED_FASTQ" \
        -zotus "$ZOTUS_FASTA" \
        -otutabout "$OTUTAB" \
        -threads "$THREADS" \
        > "$LOG" 2>&1

    cat "$LOG"

    n_samples=$(awk 'NR==1{print NF-1}' "$OTUTAB")
    n_zotus=$(awk 'END{print NR-1}' "$OTUTAB")
    echo ""
    echo "  Output: $OTUTAB ($n_zotus ZOTUs x $n_samples samples)"

    {
        echo "BigExp minsize${ms}: $n_zotus ZOTUs x $n_samples samples"
        grep -i "mapped\|matched\|total" "$LOG" | sed 's/^/  /'
        echo ""
    } >> "$OUT_DIR/mapping_stats.txt"

    echo ""
done

echo "=== Done ==="
echo "  Results: $OUT_DIR/"
cat "$OUT_DIR/mapping_stats.txt"
