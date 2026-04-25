#!/bin/bash
### General options
#BSUB -q hpcspecial
#BSUB -J nanoasv_rrndb
#BSUB -n 20
#BSUB -R "span[hosts=1] rusage[mem=8GB]"
#BSUB -M 8500MB
#BSUB -W 8:00
#BSUB -u jnesme@gmail.com
#BSUB -B
#BSUB -N
#BSUB -o logs/nanoasv_rrndb_%J.out
#BSUB -e logs/nanoasv_rrndb_%J.err

# NanoASV with rrnDB 5.10 as reference database (parallel to 11_nanoasv.sh which uses SILVA).
#
# rrnDB is curated for 16S copy number and has better species-level coverage than SILVA
# for common bacteria, but no Elusimicrobiota entries — those will fall through to the
# vsearch de novo clustering of unmatched reads.
#
# Prerequisite: run 12_build_rrndb_nanoasv.sh first to generate the NanoASV-compatible
# rrnDB FASTA at db/rrndb/SINGLELINE_rrnDB_5.10_species_taxid.fasta.gz
#
# --no-r-cleaning: R cleaning step is tuned to SILVA taxonomy keywords; skip with rrnDB.
#
# Input:  fastq_merged/barcodeXX.fastq.gz
#         metadata/nanoasv_metadata.csv
# Output: results/nanoasv_rrndb/
#
# Usage: bsub < 11b_nanoasv_rrndb.sh

set -euo pipefail

export NANOASV_PATH="/work3/josne/github/NanoASV"

source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate NanoASV

THREADS=$LSB_DJOB_NUMPROC

PROJECT_DIR="/work3/josne/Projects/AstaMSc_GRF_Igalbana"
MERGED_DIR="$PROJECT_DIR/fastq_merged"
INPUT_STAGE="$PROJECT_DIR/results/nanoasv_rrndb/input_stage"
OUT_DIR="$PROJECT_DIR/results/nanoasv_rrndb/output"
LOG_DIR="$PROJECT_DIR/logs"
METADATA_CSV="$PROJECT_DIR/metadata/nanoasv_metadata.csv"
DATABASE="$PROJECT_DIR/db/rrndb/SINGLELINE_rrnDB_5.10_plus_zotus.fasta"

if [[ ! -s "$DATABASE" ]]; then
    echo "ERROR: rrnDB NanoASV database not found at $DATABASE"
    echo "Run: bash 12_build_rrndb_nanoasv.sh && bash 13_augment_database.sh"
    exit 1
fi

mkdir -p "$OUT_DIR" "$LOG_DIR"

# --- Stage input ---
echo "=== Staging input directory ==="
rm -rf "$INPUT_STAGE"
mkdir -p "$INPUT_STAGE"

for gz in "$MERGED_DIR"/barcode*.fastq.gz; do
    barcode=$(basename "$gz" .fastq.gz)
    mkdir -p "$INPUT_STAGE/$barcode"
    ln -s "$gz" "$INPUT_STAGE/$barcode/$barcode.fastq.gz"
done

n_barcodes=$(ls -d "$INPUT_STAGE"/barcode*/ | wc -l)
echo "  Staged $n_barcodes barcodes"

# --- Metadata ---
if [[ ! -s "$METADATA_CSV" ]]; then
    echo "  WARNING: $METADATA_CSV not found — generating stub metadata."
    mkdir -p "$(dirname "$METADATA_CSV")"
    {
        echo ",SampleID"
        for gz in "$MERGED_DIR"/barcode*.fastq.gz; do
            barcode=$(basename "$gz" .fastq.gz)
            echo "$barcode,$barcode"
        done
    } > "$METADATA_CSV"
fi

cp "$METADATA_CSV" "$INPUT_STAGE/metadata.csv"
echo "  Metadata: $METADATA_CSV"

# --- Run NanoASV ---
echo ""
echo "=== Running NanoASV (rrnDB database) ==="
cd "$OUT_DIR"

bash "$NANOASV_PATH/workflow/run.sh" \
    --dir         "$INPUT_STAGE" \
    --out         "$OUT_DIR" \
    --num-process "$THREADS" \
    --database    "$DATABASE" \
    --metadata    "$INPUT_STAGE" \
    --minab       2 \
    --subsampling 100000 \
    --no-r-cleaning \
    --sam-qual    0

echo ""
echo "Done. Results in: $OUT_DIR"
ls -lh "$OUT_DIR/Results/" 2>/dev/null || true
