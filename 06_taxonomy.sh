#!/usr/bin/env bash
# Taxonomy assignment for ZOTUs using BLASTn against the NanoASV base database.
#
# The reference database (GTDB SSU + SILVA organellar) is built once by
# 08_build_gtdb_nanoasv.sh. Using the same database for both BLAST taxonomy
# and NanoASV mapping ensures consistency: a ZOTU at pident == 100% is already
# represented in the mapping database and does not need injection (step 09).
#
# Assignment rules:
#   pident >= BLAST_IDENTITY_THRESHOLD AND evalue <= BLAST_EVALUE_THRESHOLD
#     → assigned (taxonomy_assigned.tsv)
#   Hit to organellar sequence (Chloroplast / Mitochondria lineage)
#     → always assigned regardless of pident — any organellar hit is informative
#   Otherwise → unknown (taxonomy_unknown.tsv)
#
# Usage: bash 06_taxonomy.sh <zotus.fasta>
#   e.g. bash 06_taxonomy.sh pooled/zotus_minsize3.fasta
#
# Prerequisite: bash 08_build_gtdb_nanoasv.sh  (builds the reference database)
#
# Input:
#   <zotus.fasta>                                   — ZOTUs from 04_unoise3.sh
#   db/gtdb/SINGLELINE_GTDB_SSU_nanoasv.fasta       — NanoASV base database
#
# Output:
#   results/taxonomy_<zotus>/
#     taxonomy_assigned.tsv  — ZOTUs with assigned taxonomy
#     taxonomy_unknown.tsv   — ZOTUs with no hit or below threshold
#     taxonomy_all.tsv       — full table

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: bash 06_taxonomy.sh <zotus.fasta>" >&2
    echo "  e.g. bash 06_taxonomy.sh pooled/zotus_minsize3.fasta" >&2
    exit 1
fi

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV_MAIN"

ZOTUS_FASTA="$1"
ZOTU_BASE=$(basename "${ZOTUS_FASTA%.fasta}")

DB_DIR="db/gtdb"
BASE_DB="$DB_DIR/SINGLELINE_GTDB_SSU_nanoasv.fasta"
BLAST_DB="$DB_DIR/gtdb_ssu"
RESULTS_DIR="results/taxonomy_${ZOTU_BASE}"
LOG_DIR="logs/taxonomy"

mkdir -p "$RESULTS_DIR" "$LOG_DIR"

THREADS="${THREADS:-$(nproc)}"

if [[ ! -f "$BASE_DB" ]]; then
    echo "ERROR: Base database not found: $BASE_DB" >&2
    echo "Run: bash 08_build_gtdb_nanoasv.sh" >&2
    exit 1
fi

# --- Step 1: Build BLAST database from base DB (once) ---
if [[ ! -f "${BLAST_DB}.nhr" ]]; then
    echo ""
    echo "=== Building BLAST database ==="
    makeblastdb \
        -in "$BASE_DB" \
        -dbtype nucl \
        -out "$BLAST_DB" \
        -title "GTDB_SSU_plus_organellar" \
        > "$LOG_DIR/makeblastdb.log" 2>&1
    echo "  BLAST database built: $BLAST_DB"
fi

# --- Step 2: BLASTn ---
echo ""
echo "=== Running BLASTn (threads: $THREADS) ==="

BLAST_OUT="$RESULTS_DIR/blast_results.tsv"

blastn \
    -query "$ZOTUS_FASTA" \
    -db "$BLAST_DB" \
    -out "$BLAST_OUT" \
    -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs" \
    -perc_identity 80 \
    -max_target_seqs 500 \
    -num_threads "$THREADS" \
    2> "$LOG_DIR/${ZOTU_BASE}_blastn.log"
# Note: -max_target_seqs 1 does NOT guarantee the best hit (Shah et al. 2019).
# We retrieve 500 hits per query and select the true best hit by bitscore in the parser.

echo "  BLASTn done."

# --- Step 3: Parse results and assign taxonomy ---
python3 - "$BLAST_OUT" "$BASE_DB" "$ZOTUS_FASTA" \
          "$RESULTS_DIR" "$BLAST_IDENTITY_THRESHOLD" "$BLAST_EVALUE_THRESHOLD" << 'PYEOF'
import sys
import csv

blast_file  = sys.argv[1]
ref_fasta   = sys.argv[2]
query_fasta = sys.argv[3]
out_dir     = sys.argv[4]
id_thresh   = float(sys.argv[5])
ev_thresh   = float(sys.argv[6])

# Load all ZOTU IDs from the query FASTA
zotu_ids = []
with open(query_fasta) as f:
    for line in f:
        if line.startswith('>'):
            zotu_ids.append(line[1:].strip().split()[0])

# Load taxonomy from reference FASTA headers: >accession taxonomy_string
tax = {}
with open(ref_fasta) as f:
    for line in f:
        if line.startswith('>'):
            parts = line[1:].strip().split(None, 1)
            if len(parts) == 2:
                tax[parts[0]] = parts[1]

# Load BLAST results and select true best hit per query by bitscore.
# -max_target_seqs 1 does NOT guarantee the best hit (Shah et al. 2019,
# Bioinformatics) — it returns the first hit found during internal processing,
# which depends on database order. We retrieve 500 hits and select explicitly.
# Fields: qseqid sseqid pident length evalue bitscore qcovs
blast_hits = {}
with open(blast_file) as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        qseqid   = row[0]
        bitscore = float(row[5])
        if qseqid not in blast_hits or bitscore > blast_hits[qseqid]['bitscore']:
            blast_hits[qseqid] = {
                'ref':      row[1],
                'pident':   float(row[2]),
                'length':   int(row[3]),
                'evalue':   float(row[4]),
                'bitscore': bitscore,
                'qcovs':    float(row[6]),
            }

assigned = []
unknown  = []

for zotu in zotu_ids:
    hit = blast_hits.get(zotu)
    if hit is None:
        unknown.append({
            'zotu': zotu, 'ref': 'no_hit', 'pident': None,
            'length': None, 'evalue': None, 'bitscore': None,
            'qcovs': None, 'taxonomy': 'Unassigned'
        })
        continue

    lineage = tax.get(hit['ref'], 'Not in taxonomy file')
    record = {'zotu': zotu, 'taxonomy': lineage, **hit}

    is_organellar = ';Chloroplast;' in lineage or ';Mitochondria' in lineage
    if (hit['pident'] >= id_thresh and hit['evalue'] <= ev_thresh) or is_organellar:
        assigned.append(record)
    else:
        unknown.append(record)

# Write outputs
header = ['zotu', 'ref', 'pident', 'length', 'evalue', 'bitscore', 'qcovs', 'taxonomy']

def write_tsv(records, path):
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter='\t')
        writer.writeheader()
        writer.writerows(records)

write_tsv(assigned,            f"{out_dir}/taxonomy_assigned.tsv")
write_tsv(unknown,             f"{out_dir}/taxonomy_unknown.tsv")
write_tsv(assigned + unknown,  f"{out_dir}/taxonomy_all.tsv")

print(f"  Assigned (pident >= {id_thresh}%, evalue <= {ev_thresh}): {len(assigned)}")
print(f"  Unknown  (no hit or below threshold): {len(unknown)}")

# Print assigned taxonomy summary
if assigned:
    print("\n  Assigned ZOTUs:")
    for r in assigned:
        print(f"    {r['zotu']}: {r['pident']}% identity — {r['taxonomy']}")
if unknown:
    print("\n  Unknown ZOTUs:")
    for r in unknown:
        print(f"    {r['zotu']}: {r.get('ref','no_hit')} (pident={r['pident']}, evalue={r['evalue']})")
PYEOF

echo ""
echo "Done. Results in: $RESULTS_DIR"
