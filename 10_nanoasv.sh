#!/bin/bash
### General options
#BSUB -q hpcspecial
#BSUB -J nanoasv
#BSUB -n 20
#BSUB -R "span[hosts=1] rusage[mem=8GB]"
#BSUB -M 8500MB
#BSUB -W 8:00
#BSUB -u jnesme@gmail.com
#BSUB -B
#BSUB -N
#BSUB -o logs/nanoasv_%J.out
#BSUB -e logs/nanoasv_%J.err

# Run NanoASV using the GTDB SSU + project ZOTUs reference database.
#
# The database is GTDB-curated, consistent with the taxonomy from script 06
# (BLASTn). The 14 project ZOTUs are added so reads map to their exact
# consensus sequences at ~100% identity rather than to distant GTDB
# representatives. Novel Elusimicrobiota (Zotu1/2/12/14) are labelled
# unclassified_UBA9628 (confident to family; novel at genus/species level).
#
# --no-r-cleaning: NanoASV's R cleaning step uses SILVA taxonomy keywords
# (Mitochondria, Chloroplast, etc.); GTDB uses different strings, so the
# cleaning step is skipped to avoid discarding legitimate sequences.
#
# Prerequisite: run 08_build_gtdb_nanoasv.sh and 09_augment_nanoasv_db.sh first.
#
# Input:  fastq_merged/barcodeXX.fastq.gz
#         metadata/nanoasv_metadata.csv
#         db/gtdb/SINGLELINE_GTDB_SSU_plus_zotus.fasta
# Output: results/nanoasv/output/
#
# Usage: bsub < 10_nanoasv.sh

set -euo pipefail

export NANOASV_PATH="/work3/josne/github/NanoASV"

source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate NanoASV

# minimap2 has no threads: declaration in NanoASV's snakefile, so Snakemake treats each
# barcode job as 1-thread and runs --num-process jobs in parallel. minimap2 defaults to
# -t 4, so actual CPU use is num-process × 4. Divide allocated CPUs by 4 to stay within
# the LSF allocation and avoid node oversubscription.
THREADS=$(( LSB_DJOB_NUMPROC / 4 ))

PROJECT_DIR="/work3/josne/Projects/AstaMSc_GRF_Igalbana"
MERGED_DIR="$PROJECT_DIR/fastq_merged"
INPUT_STAGE="$PROJECT_DIR/results/nanoasv/input_stage"
OUT_DIR="$PROJECT_DIR/results/nanoasv/output"
LOG_DIR="$PROJECT_DIR/logs"
METADATA_CSV="$PROJECT_DIR/metadata/nanoasv_metadata.csv"
DATABASE="$PROJECT_DIR/db/gtdb/SINGLELINE_GTDB_SSU_plus_zotus.fasta"

if [[ ! -s "$DATABASE" ]]; then
    echo "ERROR: GTDB NanoASV database not found at $DATABASE"
    echo "Run: bash 08_build_gtdb_nanoasv.sh && bash 09_augment_nanoasv_db.sh"
    exit 1
fi

mkdir -p "$OUT_DIR" "$LOG_DIR"

# --- Stage input: NanoASV requires INPUT_DIR/{barcode}/{file}.fastq.gz ---
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
echo "=== Running NanoASV (GTDB SSU + ZOTUs) ==="
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
