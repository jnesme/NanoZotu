#!/bin/bash
### General options
#BSUB -q hpcspecial
#BSUB -J trim_primers
#BSUB -n 20
#BSUB -R "span[hosts=1] rusage[mem=4GB]"
#BSUB -M 4500MB
#BSUB -W 1:00
#BSUB -u jnesme@gmail.com
#BSUB -B
#BSUB -N
#BSUB -o logs/trim_primers_%J.out
#BSUB -e logs/trim_primers_%J.err

# Trim 16S primers with cutadapt, retain only reads with BOTH primers,
# size-filtered to MIN_LEN-MAX_LEN bp (after trimming).
#
# --revcomp (cutadapt >= 4.0): tries both the original read and its reverse
# complement; outputs whichever orientation the primers were found on.
# This replaces the previous two-pass approach and ensures all output reads
# are in forward orientation.
# --discard-untrimmed ensures BOTH primers must be present.
#
# Input:  fastq_merged/barcode*.fastq.gz
# Output: fastq_trimmed/barcode*.fastq.gz
#
# Usage: bsub < 02_trim_primers.sh

set -euo pipefail

# PROJECT_DIR must be an absolute path — LSF copies this script to /tmp before
# execution, so relative paths and ${BASH_SOURCE[0]} are not reliable.
PROJECT_DIR="/work3/josne/github/NanoZotu"
source "$PROJECT_DIR/config.sh"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV_MAIN"

THREADS=$LSB_DJOB_NUMPROC
# Threads given to each cutadapt call; barcodes run in parallel to fill remaining cores
CUTADAPT_THREADS=4
BARCODE_JOBS=$(( THREADS / CUTADAPT_THREADS ))
BARCODE_JOBS=$(( BARCODE_JOBS < 1 ? 1 : BARCODE_JOBS ))

INPUT_DIR="$PROJECT_DIR/fastq_merged"
OUTPUT_DIR="$PROJECT_DIR/fastq_trimmed"
LOG_DIR="$PROJECT_DIR/logs/cutadapt"

echo "Using $THREADS total threads: $BARCODE_JOBS barcodes x $CUTADAPT_THREADS threads each"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

trim_barcode() {
    local input_file="$1"
    local barcode
    barcode=$(basename "$input_file" .fastq.gz)
    local log="$LOG_DIR/${barcode}.log"
    local out_file="$OUTPUT_DIR/${barcode}.fastq.gz"

    echo "=== Processing $barcode ==="

    cutadapt \
        -j "$CUTADAPT_THREADS" \
        -g "$PRIMER_FWD" \
        -a "$PRIMER_RC_REV" \
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
export OUTPUT_DIR LOG_DIR PRIMER_FWD PRIMER_RC_REV MIN_LEN MAX_LEN CUTADAPT_THREADS

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
