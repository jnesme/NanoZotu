#!/bin/bash
#BSUB -J savont_pooled
#BSUB -q hpcspecial
#BSUB -n 20
#BSUB -R "span[hosts=1] rusage[mem=32768]"
#BSUB -M 34000MB
#BSUB -W 2:00
#BSUB -u jnesme@gmail.com
#BSUB -N
#BSUB -o logs/savont/savont_pooled_%J.out
#BSUB -e logs/savont/savont_pooled_%J.err

# Savont pooled mode — all 21 samples concatenated before ASV inference.
#
# All trimmed reads are passed together to a single savont asv run, matching
# the UNOISE3 pooling strategy (step 03). Global ASVs are defined once across
# all samples. No per-sample abundance table is produced — this run is for
# taxonomy and ASV count comparison against UNOISE3 (14 ZOTUs, minsize=3)
# and against the per-sample savont run (test_savont.sh).
#
# Input:  fastq_trimmed/barcode*.fastq.gz  (all 21 barcodes concatenated)
# Output: results/savont_pooled/
#   final_asvs.fasta        — global ASV sequences
#   species_abundance.tsv   — taxonomy profile (after classify)
#   genus_abundance.tsv
#
# Usage: bsub < test_savont_pooled.sh

set -euo pipefail

source "/work3/josne/github/NanoZotu/config.sh"
SILVA_DB="/work3/josne/Databases/silva-138.2"
TRIMMED_DIR="$PROJECT_DIR/fastq_trimmed"
OUT_DIR="$PROJECT_DIR/results/savont_pooled"
LOG_DIR="$PROJECT_DIR/logs/savont"
THREADS=$LSB_DJOB_NUMPROC

mkdir -p "$LOG_DIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate savont

echo "=== savont ASV inference — pooled mode ==="
echo "  Threads: $THREADS"
echo "  Input:   $TRIMMED_DIR/barcode*.fastq.gz ($(ls "$TRIMMED_DIR"/barcode*.fastq.gz | wc -l) files)"
echo "  Output:  $OUT_DIR"
echo ""

# --- Pooled ASV inference (all barcodes concatenated) ---
savont asv "$TRIMMED_DIR"/barcode*.fastq.gz \
    --output-dir "$OUT_DIR" \
    --threads "$THREADS" \
    --fl-16s

echo ""
echo "=== Classifying against SILVA 138.2 ==="
savont classify \
    -i "$OUT_DIR" \
    -d "$SILVA_DB" \
    --threads "$THREADS"

echo ""
echo "=== Done ==="
echo "  ASVs:            $OUT_DIR/final_asvs.fasta"
echo "  Species profile: $OUT_DIR/species_abundance.tsv"
echo "  Genus profile:   $OUT_DIR/genus_abundance.tsv"
echo ""
echo "ASV count:"
grep -c "^>" "$OUT_DIR/final_asvs.fasta" || true
