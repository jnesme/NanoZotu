#!/usr/bin/env bash
# Concatenate all fastq.gz files within each barcode folder in fastq_pass/
# Output: one merged fastq.gz per barcode in fastq_merged/

set -euo pipefail

FASTQ_PASS_DIR="fastq_pass"
OUTPUT_DIR="fastq_merged"

mkdir -p "$OUTPUT_DIR"

for barcode_dir in "$FASTQ_PASS_DIR"/barcode*; do
    barcode=$(basename "$barcode_dir")
    out_file="$OUTPUT_DIR/${barcode}.fastq.gz"

    echo "Concatenating $barcode_dir -> $out_file"
    cat "$barcode_dir"/*.fastq.gz > "$out_file"
    echo "  Done: $(ls "$barcode_dir"/*.fastq.gz | wc -l) files merged"
done

echo ""
echo "All barcodes concatenated. Output in: $OUTPUT_DIR"
ls -lh "$OUTPUT_DIR"
