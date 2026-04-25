#!/usr/bin/env bash
# Denoise dereplicated reads with usearch -unoise3 across a range of -minsize
# thresholds (1 to 8) to assess sensitivity of ZOTU recovery to this parameter.
#
# Input:  pooled/all_samples_derep.fasta
# Output: pooled/zotus_minsize{2..8}.fasta  — one file per threshold
#         logs/unoise3/unoise3_minsize{2..8}.log

set -euo pipefail

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qiime2-amplicon-2026.1

POOLED_DIR="pooled"
LOG_DIR="logs/unoise3"

mkdir -p "$LOG_DIR"

INPUT_FASTA="$POOLED_DIR/all_samples_derep.fasta"

SWEEP_REPORT="$POOLED_DIR/minsize_sweep.txt"

echo "=== UNOISE3 minsize sweep (1 -> 8) ==="
echo ""
printf "%-10s %s\n" "minsize" "ZOTUs" | tee "$SWEEP_REPORT"
printf "%-10s %s\n" "-------" "-----" | tee -a "$SWEEP_REPORT"

for minsize in $(seq 1 8); do
    out_fasta="$POOLED_DIR/zotus_minsize${minsize}.fasta"
    log="$LOG_DIR/unoise3_minsize${minsize}.log"

    usearch \
        -unoise3 "$INPUT_FASTA" \
        -zotus "$out_fasta" \
        -minsize "$minsize" \
        > "$log" 2>&1

    n_zotus=$(grep -c '^>' "$out_fasta" 2>/dev/null || echo 0)
    printf "%-10s %s\n" "$minsize" "$n_zotus" | tee -a "$SWEEP_REPORT"
done

echo ""
echo "Sweep report saved: $SWEEP_REPORT"
echo ""
echo "Output files:"
ls -lh "$POOLED_DIR"/zotus_minsize*.fasta
