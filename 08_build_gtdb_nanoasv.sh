#!/bin/bash
# Build the NanoASV base reference database: GTDB SSU + SILVA organellar sequences.
#
# NanoASV requires: >UniqueID Rank1;Rank2;...;Species  (singleline sequences)
#
# Sources:
#   GTDB SSU — bacterial and archaeal representative sequences; rank prefixes
#     (d__, p__, ...) stripped, bracketed metadata discarded.
#   SILVA 138.2 organellar subset — Chloroplast (;Chloroplast; lineage node) and
#     Mitochondria (;Mitochondria lineage node) sequences extracted from the SILVA
#     SSU file; already in NanoASV-compatible singleline format, appended as-is.
#
# This combined database is used by both step 06 (BLAST taxonomy) and step 09
# (NanoASV mapping). A ZOTU is only injected in step 09 if its exact sequence is
# absent from this database (pident < 100%).
#
# --no-r-cleaning is required with this database in NanoASV: the R cleaning step
# filters on SILVA taxonomy keywords; GTDB uses different strings. The SILVA
# organellar sequences DO use the expected keywords (Chloroplast, Mitochondria)
# so they are retained correctly.
#
# Input:  db/gtdb/bac120_ssu_reps.fna.gz   (~47k bacterial SSU reps)
#         db/gtdb/ar53_ssu_reps.fna.gz     (~3.4k archaeal SSU reps)
#         db/silva/SINGLELINE_SILVA_138.2_plus_zotus.fasta  (organellar subset extracted)
# Output: db/gtdb/SINGLELINE_GTDB_SSU_nanoasv.fasta  (uncompressed)
#
# Note: uncompressed intentionally — NanoASV format check uses plain grep (not zgrep).
#
# Usage: bash 08_build_gtdb_nanoasv.sh

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

BAC_GZ="db/gtdb/bac120_ssu_reps.fna.gz"
ARC_GZ="db/gtdb/ar53_ssu_reps.fna.gz"
SILVA_FASTA="db/silva/SINGLELINE_SILVA_138.2_plus_zotus.fasta"
OUT_FASTA="db/gtdb/SINGLELINE_GTDB_SSU_nanoasv.fasta"

echo "=== Building NanoASV base database: GTDB SSU + SILVA organellar ==="
echo "  Bacterial SSU: $BAC_GZ"
echo "  Archaeal SSU:  $ARC_GZ"
echo "  SILVA:         $SILVA_FASTA"
echo "  Output:        $OUT_FASTA"
echo ""

for f in "$BAC_GZ" "$ARC_GZ" "$SILVA_FASTA"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: $f not found"
        exit 1
    fi
done

python3 - "$BAC_GZ" "$ARC_GZ" "$SILVA_FASTA" "$OUT_FASTA" << 'PYEOF'
import sys
import gzip

bac_gz      = sys.argv[1]
arc_gz      = sys.argv[2]
silva_fasta = sys.argv[3]
out_path    = sys.argv[4]

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

print(f"  GTDB sequences written: {n_written}")
print(f"  GTDB sequences skipped: {n_skipped}")

# Append SILVA organellar sequences (Chloroplast + Mitochondria lineages).
# SILVA file is already in NanoASV-compatible singleline format; append as-is.
# Use ';Chloroplast;' (both semicolons) to match the Chloroplast lineage node
# without matching Chloroplastida (eukaryotic nuclear sequences).
n_organellar = 0
with open(silva_fasta) as f:
    include = False
    for line in f:
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            include = ';Chloroplast;' in line or ';Mitochondria' in line
            if include:
                n_organellar += 1
        if include:
            out.write(f'{line}\n')

print(f"  SILVA organellar appended: {n_organellar}")
PYEOF

echo ""
echo "Done."
ls -lh "$OUT_FASTA"
echo ""
echo "Use with NanoASV:"
echo "  --database $OUT_FASTA --no-r-cleaning"
echo ""
echo "Next: bash 09_augment_nanoasv_db.sh"
