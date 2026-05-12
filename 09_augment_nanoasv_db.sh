#!/bin/bash
# Augment the GTDB NanoASV database with ZOTUs whose sequences are not already
# represented in GTDB (pident < 100%). ZOTUs at 100% pident are skipped — the
# identical GTDB entry is sufficient and injection would cause minimap2 to split
# reads between two identical entries, distorting relative abundances.
#
# Injection + taxonomy labelling rules:
#   pident == 100%                         → skip (GTDB entry sufficient)
#   pident >= BLAST_SPECIES_THRESHOLD (98.7%), < 100%
#                                          → inject, full GTDB taxonomy
#   pident >= BLAST_IDENTITY_THRESHOLD (97%), < 98.7%
#                                          → inject, genus only (Genus sp.)
#   pident <  BLAST_IDENTITY_THRESHOLD (97%)
#                                          → inject, unclassified_<family>
#                                            (verify manually + tree in step 07)
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

ZOTU_BLOCK=$(python3 - "$ZOTU_FASTA" "$TAXONOMY_TSV" "$BLAST_IDENTITY_THRESHOLD" "$BLAST_SPECIES_THRESHOLD" << 'PYEOF'
import sys
import csv

fasta_in         = sys.argv[1]
tax_tsv          = sys.argv[2]
PIDENT_THRESHOLD = float(sys.argv[3])   # 97  — genus minimum
SPECIES_THRESHOLD = float(sys.argv[4])  # 98.7 — species minimum (Kim et al. 2014)

skipped = []
tax = {}
with open(tax_tsv) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        zotu   = row['zotu']
        gtdb   = row['taxonomy']
        pident = float(row['pident'])

        if pident == 100.0:
            skipped.append(zotu)
            continue  # exact sequence already in GTDB — skip to avoid read splitting

        parts = [p.split('__', 1)[1] if '__' in p else p for p in gtdb.split(';')]

        if pident < PIDENT_THRESHOLD:
            # Novel organism — family call may also be uncertain at low pident.
            # Cross-check taxonomy_unknown.tsv manually and verify with the
            # phylogenetic tree (07_elusimicrobiota_tree.sh) for any unknowns.
            family = parts[4] if len(parts) > 4 else 'unclassified'
            label  = f'unclassified_{family}'
            parts[5] = label
            parts[6] = label
        elif pident < SPECIES_THRESHOLD:
            # Same genus, species uncertain — keep genus, truncate species
            genus    = parts[5] if len(parts) > 5 else 'unclassified'
            parts[6] = f'{genus} sp.'

        tax[zotu] = ';'.join(parts)

if skipped:
    print(f'Skipped (100% pident): {", ".join(skipped)}', file=sys.stderr)

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
N_SKIPPED=$(python3 -c "
import csv
with open('$TAXONOMY_TSV') as f:
    print(sum(1 for r in csv.DictReader(f, delimiter='\t') if float(r['pident']) == 100.0))
")

echo ""
echo "  Injected ($N_ZOTUS):"
echo "$ZOTU_BLOCK" | grep "^>" | sed 's/^/    /'
echo ""

cat "$BASE_DB" <(echo "$ZOTU_BLOCK") > "$OUT_DB"

echo "  Base DB:   $(grep -c "^>" "$BASE_DB") sequences"
echo "  Skipped:   $N_SKIPPED (100% pident — GTDB entry sufficient)"
echo "  Injected:  $N_ZOTUS"
echo "  Output:    $(grep -c "^>" "$OUT_DB") sequences total"
echo "  File:      $OUT_DB"
echo ""
echo "Done."
echo ""
echo "Use with NanoASV:"
echo "  --database $OUT_DB --no-r-cleaning"
echo ""
echo "Next: bsub < 10_nanoasv.sh"
