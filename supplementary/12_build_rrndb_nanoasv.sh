#!/bin/bash
# Convert rrnDB species_taxid.fasta to NanoASV-compatible format.
#
# NanoASV requires:
#   >UniqueID Rank1;Rank2;...;Species
# with singleline (non-interleaved) sequences.
#
# rrnDB header format:  >TAXID:source:copy [description]
# rrnDB taxonomy TSV:   tax_id, species, genus, family, order, class, phylum, superkingdom, ...
#
# The unique ID per sequence is preserved (TAXID:source:copy).
# All copies of the same taxon receive the same taxonomy string so NanoASV's
# R step aggregates them correctly.
#
# Note: use --no-r-cleaning in NanoASV when using this database — the R
# cleaning step is tuned to SILVA taxonomy keyword patterns.
#
# Usage: bash 12_build_rrndb_nanoasv.sh
# (fast enough to run interactively; ~97k sequences)

set -euo pipefail

PROJECT_DIR="/work3/josne/Projects/AstaMSc_GRF_Igalbana"
RRNDB_DIR="$PROJECT_DIR/rrnDBv5.10_March2026"
OUT_DIR="$PROJECT_DIR/db/rrndb"

mkdir -p "$OUT_DIR"

# Uncompressed: NanoASV's format check uses plain grep, which cannot read gzip.
# minimap2 and zcat-based steps both handle uncompressed FASTA fine.
OUT_FASTA="$OUT_DIR/SINGLELINE_rrnDB_5.10_species_taxid.fasta"

echo "=== Converting rrnDB to NanoASV format ==="
echo "  Input FASTA:    $RRNDB_DIR/species_taxid.fasta"
echo "  Input taxonomy: $RRNDB_DIR/taxonomy.tsv"
echo "  Output:         $OUT_FASTA"
echo ""

python3 - "$RRNDB_DIR/species_taxid.fasta" "$RRNDB_DIR/taxonomy.tsv" "$OUT_FASTA" << 'PYEOF'
import sys
import csv
import re

fasta_in  = sys.argv[1]
tax_tsv   = sys.argv[2]
out_gz    = sys.argv[3]

# Build tax_id -> taxonomy string
# Rank order: superkingdom > phylum > class > order > family > genus > species
RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

tax = {}
with open(tax_tsv) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if not row.get('tax_id'):
            continue
        tid = row['tax_id'].strip()
        parts = [(row.get(r) or '').strip() for r in RANKS]
        # Replace empty ranks with 'unclassified' so taxonomy string is always 7 levels
        parts = [p if p else 'unclassified' for p in parts]
        tax[tid] = ';'.join(parts)

print(f"  Loaded taxonomy for {len(tax)} taxa")

# Convert FASTA: singleline, new headers
n_written = 0
n_no_tax  = 0

with open(fasta_in) as f_in, open(out_gz, 'w') as f_out:
    header   = None
    seq_parts = []

    def write_record(f_out, header, seq_parts, tax, stats):
        if header is None:
            return
        m = re.match(r'>(\d+):', header)
        if not m:
            stats['no_tax'] += 1
            return
        tid = m.group(1)
        uid = header[1:].split()[0].replace(':', '_')  # TAXID_source_copy — colons are Newick branch-length separator, FastTree truncates at them
        taxonomy = tax.get(tid)
        if taxonomy is None:
            stats['no_tax'] += 1
            return
        seq = ''.join(seq_parts)
        f_out.write(f'>{uid} {taxonomy}\n{seq}\n')
        stats['written'] += 1

    stats = {'written': 0, 'no_tax': 0}

    for line in f_in:
        line = line.rstrip()
        if line.startswith('>'):
            write_record(f_out, header, seq_parts, tax, stats)
            header    = line
            seq_parts = []
        else:
            seq_parts.append(line)
    write_record(f_out, header, seq_parts, tax, stats)

print(f"  Sequences written:          {stats['written']}")
print(f"  Sequences skipped (no tax): {stats['no_tax']}")
PYEOF

echo ""
echo "Done."
ls -lh "$OUT_FASTA"
echo ""
echo "Use with NanoASV:"
echo "  --database $OUT_FASTA --no-r-cleaning"
echo ""
echo "Note: uncompressed intentionally — NanoASV format check uses plain grep (not zgrep)."
