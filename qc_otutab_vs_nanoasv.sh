#!/bin/bash
# QC: compare per-sample ZOTU read fractions between step 05 (usearch -otutab)
# and step 10 (NanoASV). The two methods should agree: ZOTUs should represent
# the majority of classified reads, and per-ZOTU relative abundances should be
# consistent. A large discrepancy indicates a pipeline issue (e.g. wrong ZOTU
# FASTA, mismatched minsize, NanoASV database not augmented with ZOTUs).
#
# Input:  results/zotu_table_zotus_minsize${MINSIZE_WORKING}.txt   (step 05)
#         results/nanoasv/output/Results/CSV/Taxonomy-Abundance_table.csv (step 10)
# Output: printed comparison tables
#
# Usage: bash qc_otutab_vs_nanoasv.sh

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

ZOTU_TABLE="results/zotu_table_zotus_minsize${MINSIZE_WORKING}.txt"
NANOASV_CSV="results/nanoasv/output/Results/CSV/Taxonomy-Abundance_table.csv"

if [[ ! -s "$ZOTU_TABLE" ]]; then
    echo "ERROR: $ZOTU_TABLE not found — run step 05 first"
    exit 1
fi
if [[ ! -s "$NANOASV_CSV" ]]; then
    echo "ERROR: $NANOASV_CSV not found — run step 10 first"
    exit 1
fi

MINSIZE_WORKING="$MINSIZE_WORKING" python3 - "$ZOTU_TABLE" "$NANOASV_CSV" << 'PYEOF'
import sys, csv, os

zotu_table_path = sys.argv[1]
nanoasv_path    = sys.argv[2]

# --- Load otutab ---
otutab   = {}   # {zotu: {sample: count}}
samples_otu = []
with open(zotu_table_path) as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    samples_otu = header[1:]
    for row in reader:
        zotu = row[0]
        otutab[zotu] = {samples_otu[i]: int(row[i+1]) for i in range(len(samples_otu))}

# --- Load NanoASV CSV ---
nanoasv  = {}   # {ref_id: {sample: count}}
samples_na = []
with open(nanoasv_path) as f:
    reader = csv.reader(f)
    header = [c.strip('"') for c in next(reader)]
    sample_idx = [(i, h) for i, h in enumerate(header) if h.startswith('barcode')]
    samples_na = [h for _, h in sample_idx]
    for row in reader:
        ref_id = row[0].strip('"')
        nanoasv[ref_id] = {}
        for i, sample in sample_idx:
            val = row[i].strip('"')
            nanoasv[ref_id][sample] = int(float(val)) if val not in ('', 'NA') else 0

# Use the intersection of samples present in both
samples = [s for s in samples_otu if s in samples_na]

zotus = sorted(otutab.keys())

# --- Per-sample summary ---
print("=== Per-sample ZOTU read counts ===")
print(f"{'Sample':<14} {'OTUtab_ZOTUs':>13} {'NanoASV_ZOTUs':>14} {'NanoASV_total':>14} {'ZOTU%_otutab':>13} {'ZOTU%_nanoasv':>14}")
print("-" * 84)

otu_total_all = 0
na_zotu_total_all = 0
na_total_all = 0

for sample in samples:
    otu_zotu = sum(otutab[z][sample] for z in zotus)
    na_zotu  = sum(nanoasv[z][sample] for z in zotus if z in nanoasv)
    na_total = sum(nanoasv[ref][sample] for ref in nanoasv)

    otu_pct = 100 * otu_zotu / na_total if na_total > 0 else 0
    na_pct  = 100 * na_zotu  / na_total if na_total > 0 else 0

    otu_total_all    += otu_zotu
    na_zotu_total_all += na_zotu
    na_total_all     += na_total

    print(f"{sample:<14} {otu_zotu:>13,} {na_zotu:>14,} {na_total:>14,} {otu_pct:>12.1f}% {na_pct:>13.1f}%")

print("-" * 84)
print(f"{'TOTAL':<14} {otu_total_all:>13,} {na_zotu_total_all:>14,} {na_total_all:>14,} "
      f"{100*otu_total_all/na_total_all:>12.1f}% {100*na_zotu_total_all/na_total_all:>13.1f}%")

# --- Per-ZOTU relative abundance ---
print()
print("=== Per-ZOTU relative abundance (% of total ZOTU pool, summed across samples) ===")
print(f"{'ZOTU':<10} {'OTUtab%':>9} {'NanoASV%':>9} {'delta':>7}")
print("-" * 40)

otu_pool  = sum(otutab[z][s] for z in zotus for s in samples)
na_pool   = sum(nanoasv.get(z, {}).get(s, 0) for z in zotus for s in samples)

for zotu in zotus:
    otu_count = sum(otutab[zotu][s] for s in samples)
    na_count  = sum(nanoasv.get(zotu, {}).get(s, 0) for s in samples)
    otu_pct = 100 * otu_count / otu_pool if otu_pool > 0 else 0
    na_pct  = 100 * na_count  / na_pool  if na_pool  > 0 else 0
    delta   = na_pct - otu_pct
    flag    = " !" if abs(delta) > 5 else ""
    print(f"{zotu:<10} {otu_pct:>8.2f}% {na_pct:>8.2f}% {delta:>+6.2f}%{flag}")

PYEOF
