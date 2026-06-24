#!/usr/bin/env bash
# compare_zotus.sh — Cross-dataset ZOTU comparison (Prelim16S vs BigExp)
#
# Bidirectional BLASTn of ZOTUs between two independent datasets from the
# same biological system. Identifies shared ZOTUs as cross-validation:
# a ZOTU found independently in both datasets is unlikely to be artifact.
#
# Reports best hit per query per direction at two thresholds (99% and 100%).
# Reciprocal best hits (A→B and B→A) are flagged as strongest matches.
#
# Input:
#   Prelim16S:  pooled/zotus_minsize${MINSIZE_WORKING}.fasta  (this project)
#   BigExp:     BIGEXP_DIR/pooled/zotus_minsize5.fasta
#               BIGEXP_DIR/pooled/zotus_minsize8.fasta
#
# Output:
#   results/compare_zotus/prelim_vs_bigexp5_100pct.tsv
#   results/compare_zotus/prelim_vs_bigexp5_99pct.tsv
#   results/compare_zotus/prelim_vs_bigexp8_100pct.tsv
#   results/compare_zotus/prelim_vs_bigexp8_99pct.tsv
#   results/compare_zotus/summary.txt
#   results/compare_zotus/blast/   — raw BLAST output for re-filtering
#
# Usage: bash compare_zotus.sh
#   Run interactively from the project root. No LSF needed.
#   To tweak filters without re-running BLAST, edit MIN_ALN_LENGTH / thresholds
#   below and re-run — existing blast/ files are reused when present.
#
# Dependencies: BLAST+ (makeblastdb, blastn), awk

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

# ─── Configuration ────────────────────────────────────────────────────────────

BIGEXP_DIR="/work3/josne/Projects/Igalbana_Microbiome_GRF/Igalbana_Microbiome_16S/NanoZotu"

PRELIM_FASTA="$PROJECT_DIR/pooled/zotus_minsize${MINSIZE_WORKING}.fasta"
PRELIM_TAX="$PROJECT_DIR/results/taxonomy_zotus_minsize${MINSIZE_WORKING}/taxonomy_all.tsv"

BIGEXP_MINSIZES=(5 8)

MIN_ALN_LENGTH=1300   # reject partial alignments shorter than this (matches MIN_LEN in config.sh)

OUT_DIR="$PROJECT_DIR/results/compare_zotus"
BLAST_DIR="$OUT_DIR/blast"
TMP_DIR="$OUT_DIR/tmp"

# ─── Validate inputs ─────────────────────────────────────────────────────────

for f in "$PRELIM_FASTA" "$PRELIM_TAX"; do
    [ -f "$f" ] || { echo "ERROR: not found: $f"; exit 1; }
done

for ms in "${BIGEXP_MINSIZES[@]}"; do
    for f in "$BIGEXP_DIR/pooled/zotus_minsize${ms}.fasta" \
             "$BIGEXP_DIR/results/taxonomy_zotus_minsize${ms}/taxonomy_all.tsv"; do
        [ -f "$f" ] || { echo "ERROR: not found: $f"; exit 1; }
    done
done

mkdir -p "$OUT_DIR" "$BLAST_DIR" "$TMP_DIR"

# ─── Prefix FASTA headers ────────────────────────────────────────────────────

echo "=== Preparing prefixed FASTAs ==="

awk '/^>/{sub(/>/,">Prelim_")} {print}' "$PRELIM_FASTA" > "$TMP_DIR/prelim.fasta"
echo "  Prelim16S: $(grep -c '^>' "$TMP_DIR/prelim.fasta") ZOTUs"

for ms in "${BIGEXP_MINSIZES[@]}"; do
    awk '/^>/{sub(/>/,">BigExp_")} {print}' \
        "$BIGEXP_DIR/pooled/zotus_minsize${ms}.fasta" > "$TMP_DIR/bigexp${ms}.fasta"
    echo "  BigExp minsize${ms}: $(grep -c '^>' "$TMP_DIR/bigexp${ms}.fasta") ZOTUs"
done

# ─── Build taxonomy lookup files (zotu<TAB>taxonomy) ─────────────────────────

echo ""
echo "=== Building taxonomy lookups ==="

awk -F'\t' 'NR>1 {print "Prelim_"$1"\t"$8}' "$PRELIM_TAX" \
    | sort -t$'\t' -k1,1 > "$TMP_DIR/tax_prelim.tsv"

for ms in "${BIGEXP_MINSIZES[@]}"; do
    awk -F'\t' 'NR>1 {print "BigExp_"$1"\t"$8}' \
        "$BIGEXP_DIR/results/taxonomy_zotus_minsize${ms}/taxonomy_all.tsv" \
        | sort -t$'\t' -k1,1 > "$TMP_DIR/tax_bigexp${ms}.tsv"
done

# ─── Function: run one comparison ────────────────────────────────────────────

run_comparison() {
    local ms=$1
    local label="bigexp${ms}"

    echo ""
    echo "=== Prelim16S vs BigExp minsize${ms} (min_aln_length=${MIN_ALN_LENGTH}) ==="

    local fwd_raw="$BLAST_DIR/fwd_${label}_raw.tsv"
    local rev_raw="$BLAST_DIR/rev_${label}_raw.tsv"

    # Build BLAST DBs (always rebuild — cheap)
    makeblastdb -in "$TMP_DIR/prelim.fasta" -dbtype nucl \
        -out "$TMP_DIR/db_prelim" -logfile /dev/null
    makeblastdb -in "$TMP_DIR/${label}.fasta" -dbtype nucl \
        -out "$TMP_DIR/db_${label}" -logfile /dev/null

    # Bidirectional BLASTn — reuse saved output if present
    # perc_identity 97 captures extra hits for re-filtering at higher thresholds
    if [ -s "$fwd_raw" ] && [ -s "$rev_raw" ]; then
        echo "  Reusing saved BLAST output in $BLAST_DIR/"
    else
        echo "  Running BLASTn..."
        blastn -query "$TMP_DIR/prelim.fasta" -db "$TMP_DIR/db_${label}" \
            -perc_identity 97 -outfmt "6 qseqid sseqid pident length evalue bitscore" \
            -evalue 1e-10 \
            > "$fwd_raw" || true

        blastn -query "$TMP_DIR/${label}.fasta" -db "$TMP_DIR/db_prelim" \
            -perc_identity 97 -outfmt "6 qseqid sseqid pident length evalue bitscore" \
            -evalue 1e-10 \
            > "$rev_raw" || true
    fi

    # Filter by minimum alignment length, then best hit per query (highest pident, then bitscore)
    awk -F'\t' -v min="$MIN_ALN_LENGTH" '$4 >= min' "$fwd_raw" \
        | sort -t$'\t' -k1,1 -k3,3nr -k6,6nr \
        | awk -F'\t' '!seen[$1]++ {print}' > "$TMP_DIR/fwd_${label}_best.tsv"

    awk -F'\t' -v min="$MIN_ALN_LENGTH" '$4 >= min' "$rev_raw" \
        | sort -t$'\t' -k1,1 -k3,3nr -k6,6nr \
        | awk -F'\t' '!seen[$1]++ {print}' > "$TMP_DIR/rev_${label}_best.tsv"

    # For each threshold, produce output
    for pct in 100 99; do
        local out_file="$OUT_DIR/prelim_vs_${label}_${pct}pct.tsv"

        # Collect forward best hits at this threshold
        awk -F'\t' -v pct="$pct" '$3 >= pct' "$TMP_DIR/fwd_${label}_best.tsv" \
            > "$TMP_DIR/fwd_${label}_${pct}.tsv"

        # Collect reverse best hits at this threshold
        awk -F'\t' -v pct="$pct" '$3 >= pct' "$TMP_DIR/rev_${label}_best.tsv" \
            > "$TMP_DIR/rev_${label}_${pct}.tsv"

        # Build combined match table with reciprocal flag and taxonomy
        awk -F'\t' -v pct="$pct" \
            -v tax_prelim="$TMP_DIR/tax_prelim.tsv" \
            -v tax_bigexp="$TMP_DIR/tax_bigexp${ms}.tsv" \
        'BEGIN {
            OFS = "\t"
            # Load taxonomy lookups
            while ((getline line < tax_prelim) > 0) {
                split(line, a, "\t"); tax_p[a[1]] = a[2]
            }
            close(tax_prelim)
            while ((getline line < tax_bigexp) > 0) {
                split(line, a, "\t"); tax_b[a[1]] = a[2]
            }
            close(tax_bigexp)
        }
        FILENAME ~ /fwd_/ {
            fwd[$1] = $2; fwd_pident[$1] = $3; fwd_len[$1] = $4
            pair[$1 FS $2] = 1
            next
        }
        FILENAME ~ /rev_/ {
            rev[$1] = $2; rev_pident[$1] = $3; rev_len[$1] = $4
            pair[$2 FS $1] = 1
            next
        }
        END {
            print "prelim_zotu", "bigexp_zotu", "pident", "aln_length", "direction", "reciprocal", "prelim_taxonomy", "bigexp_taxonomy"

            # Forward hits
            for (q in fwd) {
                s = fwd[q]
                recip = (rev[s] == q) ? "yes" : "no"
                tp = (q in tax_p) ? tax_p[q] : "NA"
                tb = (s in tax_b) ? tax_b[s] : "NA"
                print q, s, fwd_pident[q], fwd_len[q], "prelim->bigexp", recip, tp, tb
            }

            # Reverse hits not already covered by forward
            for (q in rev) {
                s = rev[q]
                if (!(s in fwd) || fwd[s] != q) {
                    tp = (s in tax_p) ? tax_p[s] : "NA"
                    tb = (q in tax_b) ? tax_b[q] : "NA"
                    recip = "no"
                    print s, q, rev_pident[q], rev_len[q], "bigexp->prelim", recip, tp, tb
                }
            }
        }' "$TMP_DIR/fwd_${label}_${pct}.tsv" "$TMP_DIR/rev_${label}_${pct}.tsv" \
            > "$out_file"

        # Counts
        local n_total=$(( $(wc -l < "$out_file") - 1 ))
        local n_recip=$(awk -F'\t' 'NR>1 && $6=="yes"' "$out_file" | wc -l)
        local n_prelim_hit=$(awk -F'\t' 'NR>1 {print $1}' "$out_file" | sort -u | wc -l)
        local n_bigexp_hit=$(awk -F'\t' 'NR>1 {print $2}' "$out_file" | sort -u | wc -l)

        local n_prelim_total=$(grep -c '^>' "$TMP_DIR/prelim.fasta")
        local n_bigexp_total=$(grep -c '^>' "$TMP_DIR/${label}.fasta")

        echo "  ${pct}%: ${n_total} matches (${n_recip} reciprocal), ${n_prelim_hit}/${n_prelim_total} Prelim ZOTUs, ${n_bigexp_hit}/${n_bigexp_total} BigExp ZOTUs"

        # Append to summary
        {
            echo "Prelim16S (minsize${MINSIZE_WORKING}) vs BigExp (minsize${ms}) — ${pct}% identity, min_aln=${MIN_ALN_LENGTH}"
            echo "  Total matches:        $n_total"
            echo "  Reciprocal best hits: $n_recip"
            echo "  Prelim ZOTUs matched: $n_prelim_hit / $n_prelim_total"
            echo "  BigExp ZOTUs matched: $n_bigexp_hit / $n_bigexp_total"
            echo ""
        } >> "$OUT_DIR/summary.txt"
    done
}

# ─── Run comparisons ─────────────────────────────────────────────────────────

> "$OUT_DIR/summary.txt"   # truncate summary

for ms in "${BIGEXP_MINSIZES[@]}"; do
    run_comparison "$ms"
done

# ─── Cleanup ──────────────────────────────────────────────────────────────────

rm -rf "$TMP_DIR"

echo ""
echo "=== Done ==="
echo "  Results: $OUT_DIR/"
echo "  Summary: $OUT_DIR/summary.txt"
echo ""
cat "$OUT_DIR/summary.txt"
