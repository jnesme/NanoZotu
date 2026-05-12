#!/bin/bash
# Augment the GTDB NanoASV database with the 14 project ZOTUs (UNOISE3, minsize=3).
#
# Adding ZOTUs to the reference database ensures that project organisms map to their
# exact consensus sequences at ~100% identity, rather than to the nearest GTDB
# representative at lower identity. Without this step, novel organisms (e.g.
# Elusimicrobiota at 88% to GWA2-66-18) map spuriously to unrelated references —
# minimap2 does not refuse to map at low identity.
#
# Taxonomy labelling follows the 97% pident threshold used in script 06:
#   pident >= 97%  → full GTDB taxonomy (known organism)
#   pident <  97%  → confident to family; genus+species → unclassified_<family>
#                    (Zotu1/2/12/14: unclassified_UBA9628)
#
# Prerequisite: run 08_build_gtdb_nanoasv.sh first.
#
# Input:  db/gtdb/SINGLELINE_GTDB_SSU_nanoasv.fasta
#         pooled/zotus_minsize3.fasta
#         results/taxonomy_zotus_minsize3/taxonomy_all.tsv
# Output: db/gtdb/SINGLELINE_GTDB_SSU_plus_zotus.fasta
#
# Usage: bash 09_augment_nanoasv_db.sh

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

BASE_DB="db/gtdb/SINGLELINE_GTDB_SSU_nanoasv.fasta"
ZOTU_FASTA="pooled/zotus_minsize${MINSIZE_WORKING}.fasta"
TAXONOMY_TSV="results/taxonomy_zotus_minsize${MINSIZE_WORKING}/taxonomy_all.tsv"
OUT_DB="db/gtdb/SINGLELINE_GTDB_SSU_plus_zotus.fasta"

if [[ ! -s "$BASE_DB" ]]; then
    echo "ERROR: GTDB NanoASV base database not found at $BASE_DB"
    echo "Run: bash 08_build_gtdb_nanoasv.sh"
    exit 1
fi

echo "=== Augmenting GTDB NanoASV database with project ZOTUs ==="

ZOTU_BLOCK=$(python3 - "$ZOTU_FASTA" "$TAXONOMY_TSV" "$BLAST_IDENTITY_THRESHOLD" << 'PYEOF'
import sys
import csv

fasta_in = sys.argv[1]
tax_tsv  = sys.argv[2]

PIDENT_THRESHOLD = float(sys.argv[3])

tax = {}
with open(tax_tsv) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        zotu   = row['zotu']
        gtdb   = row['taxonomy']
        pident = float(row['pident'])

        # Strip GTDB rank prefixes
        parts = [p.split('__', 1)[1] if '__' in p else p for p in gtdb.split(';')]

        if pident < PIDENT_THRESHOLD:
            # Novel below genus: confident to family only
            family = parts[4] if len(parts) > 4 else 'unclassified'
            label  = f'unclassified_{family}'
            parts[5] = label
            parts[6] = label

        tax[zotu] = ';'.join(parts)

header    = None
seq_parts = []

def emit(header, seq_parts, tax):
    if header is None:
        return
    zotu = header[1:].strip()
    taxonomy = tax.get(zotu)
    if taxonomy is None:
        return
    print(f'>{zotu} {taxonomy}')
    print(''.join(seq_parts))

with open(fasta_in) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            emit(header, seq_parts, tax)
            header    = line
            seq_parts = []
        else:
            seq_parts.append(line)
emit(header, seq_parts, tax)
PYEOF
)

N_ZOTUS=$(echo "$ZOTU_BLOCK" | grep -c "^>" || true)

echo ""
echo "  ZOTU entries:"
echo "$ZOTU_BLOCK" | grep "^>" | sed 's/^/    /'
echo ""

cat "$BASE_DB" <(echo "$ZOTU_BLOCK") > "$OUT_DB"

echo "  Base DB:   $(grep -c "^>" "$BASE_DB") sequences"
echo "  ZOTUs:     $N_ZOTUS"
echo "  Output:    $(grep -c "^>" "$OUT_DB") sequences total"
echo "  File:      $OUT_DB"
echo ""
echo "Done."
echo ""
echo "Use with NanoASV:"
echo "  --database $OUT_DB --no-r-cleaning"
echo ""
echo "Next: bsub < 10_nanoasv.sh"
