#!/usr/bin/env bash
# Trim 16S primers with cutadapt, retain only reads with BOTH primers,
# size-filtered to 1300-1800 bp (after trimming).
#
# Primers:
#   FWD (27F):  AGRGTTYGATYMTGGCTCAG
#   REV (1492R): RGYTACCTTGTTACGACTT
#
# --revcomp (cutadapt >= 4.0): tries both the original read and its reverse
# complement; outputs whichever orientation the primers were found on.
# This replaces the previous two-pass approach and ensures all output reads
# are in forward orientation.
# --discard-untrimmed ensures BOTH primers must be present.

set -euo pipefail

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qiime2-amplicon-2026.1

INPUT_DIR="fastq_merged"
OUTPUT_DIR="fastq_trimmed"
LOG_DIR="logs/cutadapt"

# Total CPUs available to this job (override via env: THREADS=16 bash 02_trim_primers.sh)
THREADS="${THREADS:-$(nproc)}"
# Threads given to each cutadapt call; barcodes run in parallel to fill remaining cores
CUTADAPT_THREADS=4
BARCODE_JOBS=$(( THREADS / CUTADAPT_THREADS ))
BARCODE_JOBS=$(( BARCODE_JOBS < 1 ? 1 : BARCODE_JOBS ))

echo "Using $THREADS total threads: $BARCODE_JOBS barcodes x $CUTADAPT_THREADS threads each"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# Primer sequences
FWD="AGRGTTYGATYMTGGCTCAG"
RC_REV="AAGTCGTAACAAGGTARCY"    # reverse complement of REV (RGYTACCTTGTTACGACTT)

MIN_LEN=1300
MAX_LEN=1800

trim_barcode() {
    local input_file="$1"
    local barcode
    barcode=$(basename "$input_file" .fastq.gz)
    local log="$LOG_DIR/${barcode}.log"
    local out_file="$OUTPUT_DIR/${barcode}.fastq.gz"

    echo "=== Processing $barcode ==="

    cutadapt \
        -j "$CUTADAPT_THREADS" \
        -g "$FWD" \
        -a "$RC_REV" \
        --revcomp \
        --discard-untrimmed \
        --minimum-length "$MIN_LEN" \
        --maximum-length "$MAX_LEN" \
        --error-rate 0.15 \
        --overlap 15 \
        --output "$out_file" \
        "$input_file" \
        > "$log" 2>&1

    local n_reads
    n_reads=$(zcat "$out_file" | awk 'NR%4==1' | wc -l)
    if [[ "$n_reads" -eq 0 ]]; then
        echo "  [$barcode] WARNING: no reads passed filters — output file is empty" >&2
    else
        echo "  [$barcode] Retained reads: $n_reads"
    fi
    grep -E "Reads written|Reads with adapter" "$log" | sed "s/^/  [$barcode] /"
}

export -f trim_barcode
export OUTPUT_DIR LOG_DIR FWD RC_REV MIN_LEN MAX_LEN CUTADAPT_THREADS

# Run barcodes in parallel (GNU parallel preferred, fallback to xargs)
if command -v parallel &>/dev/null; then
    parallel -j "$BARCODE_JOBS" trim_barcode ::: "$INPUT_DIR"/barcode*.fastq.gz
else
    printf '%s\n' "$INPUT_DIR"/barcode*.fastq.gz \
        | xargs -P "$BARCODE_JOBS" -I{} bash -c 'trim_barcode "$@"' _ {}
fi

echo ""
echo "Done. Trimmed files in: $OUTPUT_DIR"
ls -lh "$OUTPUT_DIR"/*.fastq.gz
