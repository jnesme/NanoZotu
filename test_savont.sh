#!/bin/bash
#BSUB -J savont_test
#BSUB -q hpcspecial
#BSUB -n 20
#BSUB -R "span[hosts=1] rusage[mem=8192]"
#BSUB -M 8500MB
#BSUB -W 4:00
#BSUB -u jnesme@gmail.com
#BSUB -N
#BSUB -o logs/savont/savont_%J.out
#BSUB -e logs/savont/savont_%J.err

# Standalone savont test — independent of the UNOISE3 pipeline.
#
# Runs savont ASV inference on trimmed reads (step 02 output), one job per
# barcode, then merges all samples and classifies against SILVA 138.2.
# Results are compared against the 14 ZOTUs from UNOISE3 (minsize=3).
#
# Input:  fastq_trimmed/barcode*.fastq.gz
# Output: results/savont/
#   per_sample/<barcode>/  — per-sample ASVs and cluster assignments
#   merged/                — merged feature table across all samples
#   merged/species_abundance.tsv / genus_abundance.tsv — taxonomy profiles
#
# Usage: bsub < test_savont.sh

set -euo pipefail

source "/work3/josne/github/NanoZotu/config.sh"
SILVA_DB="/work3/josne/Databases/silva-138.2"
TRIMMED_DIR="$PROJECT_DIR/fastq_trimmed"
OUT_DIR="$PROJECT_DIR/results/savont"
LOG_DIR="$PROJECT_DIR/logs/savont"
THREADS=$LSB_DJOB_NUMPROC

mkdir -p "$OUT_DIR/per_sample" "$LOG_DIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate savont

echo "=== savont ASV inference ==="
echo "  Threads: $THREADS"
echo "  SILVA DB: $SILVA_DB"
echo ""

# --- Per-sample ASV inference ---
sample_dirs=()
sample_labels=()

for fq in "$TRIMMED_DIR"/barcode*.fastq.gz; do
    barcode=$(basename "$fq" .fastq.gz)
    out="$OUT_DIR/per_sample/$barcode"
    echo "--- $barcode ---"
    savont asv "$fq" \
        --output-dir "$out" \
        --threads "$THREADS" \
        --fl-16s
    sample_dirs+=("$out")
    sample_labels+=("$barcode")
    echo ""
done

# --- Merge all samples ---
echo "=== Merging $(${#sample_dirs[@]}) samples ==="
savont merge \
    -i "${sample_dirs[@]}" \
    --relabel "${sample_labels[@]}" \
    -o "$OUT_DIR/merged"

echo ""

# --- Taxonomy classification against SILVA 138.2 ---
echo "=== Classifying against SILVA 138.2 ==="
savont classify \
    -i "$OUT_DIR/merged" \
    -d "$SILVA_DB" \
    --threads "$THREADS"

echo ""
echo "=== Done ==="
echo "  ASVs:            $OUT_DIR/merged/final_asvs.fasta"
echo "  Feature table:   $OUT_DIR/merged/feature-table.tsv"
echo "  Species profile: $OUT_DIR/merged/species_abundance.tsv"
echo "  Genus profile:   $OUT_DIR/merged/genus_abundance.tsv"
echo ""
echo "ASV count:"
grep -c "^>" "$OUT_DIR/merged/final_asvs.fasta" || true
