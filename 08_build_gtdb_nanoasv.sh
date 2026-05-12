#!/bin/bash
# Convert GTDB SSU representative sequences to NanoASV-compatible singleline FASTA.
#
# NanoASV requires: >UniqueID Rank1;Rank2;...;Species  (singleline sequences)
#
# GTDB SSU header format:
#   >AccessionID d__Kingdom;p__Phylum;...;s__Species [locus_tag=...] [location=...]
# Sequences are already singleline in the GTDB distribution.
#
# This script strips the rank prefixes (d__, p__, ...) and discards bracketed
# metadata so each header becomes: >AccessionID Kingdom;Phylum;...;Species
#
# --no-r-cleaning is required when using this database in NanoASV: the R
# cleaning step filters with SILVA taxonomy keywords (Mitochondria, Chloroplast,
# etc.); GTDB uses different strings so filtering is skipped to avoid data loss.
#
# Input:  db/gtdb/bac120_ssu_reps.fna.gz   (~47k bacterial SSU reps)
#         db/gtdb/ar53_ssu_reps.fna.gz     (~3.4k archaeal SSU reps)
# Output: db/gtdb/SINGLELINE_GTDB_SSU_nanoasv.fasta  (uncompressed)
#
# Note: uncompressed intentionally — NanoASV format check uses plain grep (not zgrep).
#
# Usage: bash 08_build_gtdb_nanoasv.sh

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

BAC_GZ="db/gtdb/bac120_ssu_reps.fna.gz"
ARC_GZ="db/gtdb/ar53_ssu_reps.fna.gz"
OUT_FASTA="db/gtdb/SINGLELINE_GTDB_SSU_nanoasv.fasta"

echo "=== Converting GTDB SSU to NanoASV format ==="
echo "  Bacterial SSU: $BAC_GZ"
echo "  Archaeal SSU:  $ARC_GZ"
echo "  Output:        $OUT_FASTA"
echo ""

for f in "$BAC_GZ" "$ARC_GZ"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: $f not found"
        exit 1
    fi
done

python3 - "$BAC_GZ" "$ARC_GZ" "$OUT_FASTA" << 'PYEOF'
import sys
import gzip

bac_gz  = sys.argv[1]
arc_gz  = sys.argv[2]
out_path = sys.argv[3]

RANK_PREFIXES = ('d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')

def strip_taxonomy(raw):
    # raw: "d__Bacteria;p__Pseudomonadota;...;s__Species [locus_tag=...] ..."
    # Discard everything from the first '[' onward (bracketed GTDB metadata)
    tax_field = raw.split('[')[0].strip()
    parts = tax_field.split(';')
    stripped = [p.split('__', 1)[1] if '__' in p else p for p in parts]
    return ';'.join(stripped)

n_written = 0
n_skipped = 0

with open(out_path, 'w') as out:
    for gz_path in (bac_gz, arc_gz):
        with gzip.open(gz_path, 'rt') as f:
            header   = None
            sequence = None
            for line in f:
                line = line.rstrip()
                if line.startswith('>'):
                    if header is not None and sequence:
                        out.write(f'{header}\n{sequence}\n')
                        n_written += 1
                    parts = line[1:].split(None, 1)
                    if len(parts) < 2 or not any(parts[1].startswith(p) for p in RANK_PREFIXES):
                        header = None
                        sequence = None
                        n_skipped += 1
                        continue
                    accession = parts[0]
                    taxonomy  = strip_taxonomy(parts[1])
                    header    = f'>{accession} {taxonomy}'
                    sequence  = None
                else:
                    if header is not None:
                        sequence = line  # sequences are singleline in GTDB distribution
            if header is not None and sequence:
                out.write(f'{header}\n{sequence}\n')
                n_written += 1

print(f"  Sequences written: {n_written}")
print(f"  Sequences skipped: {n_skipped}")
PYEOF

echo ""
echo "Done."
ls -lh "$OUT_FASTA"
echo ""
echo "Use with NanoASV:"
echo "  --database $OUT_FASTA --no-r-cleaning"
echo ""
echo "Next: bash 09_augment_nanoasv_db.sh"
