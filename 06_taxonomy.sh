#!/usr/bin/env bash
# Taxonomy assignment for ZOTUs using BLASTn against GTDB SSU reference.
#
# Mixed approach:
#   - High identity + low e-value hits → direct GTDB taxonomy assignment
#   - No hit / low identity            → retained as de novo novel sequences
#
# Usage: bash 06_taxonomy.sh <zotus.fasta>
#   e.g. bash 06_taxonomy.sh pooled/zotus_minsize3.fasta
#
# Input:
#   <zotus.fasta>                                           — ZOTUs from 04_unoise3.sh
#   db/gtdb/bac120_ssu_reps.fna.gz + ar53_ssu_reps.fna.gz  — GTDB SSU sequences
#   db/gtdb/bac120_taxonomy.tsv.gz + ar53_taxonomy.tsv.gz   — GTDB taxonomy
#
# Output:
#   results/taxonomy_<zotus>/
#     taxonomy_assigned.tsv  — ZOTUs with high-confidence taxonomy
#     taxonomy_unknown.tsv   — ZOTUs with no/low-confidence hit
#     taxonomy_all.tsv       — full table including all ZOTUs

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
BLAST_DB="$DB_DIR/gtdb_ssu"
RESULTS_DIR="results/taxonomy_${ZOTU_BASE}"
LOG_DIR="logs/taxonomy"

mkdir -p "$RESULTS_DIR" "$LOG_DIR"

BAC_SSU="$DB_DIR/bac120_ssu_reps.fna.gz"
AR_SSU="$DB_DIR/ar53_ssu_reps.fna.gz"
BAC_TAX="$DB_DIR/bac120_taxonomy.tsv.gz"
AR_TAX="$DB_DIR/ar53_taxonomy.tsv.gz"

COMBINED_SSU="$DB_DIR/gtdb_ssu_combined.fna"
COMBINED_TAX="$DB_DIR/gtdb_taxonomy_combined.tsv"

THREADS="${THREADS:-$(nproc)}"

# --- Step 1: Build combined reference (once) ---
if [[ ! -f "$COMBINED_SSU" ]]; then
    echo "=== Decompressing and combining GTDB SSU sequences ==="
    zcat "$BAC_SSU" "$AR_SSU" > "$COMBINED_SSU"
    echo "  $(grep -c '^>' "$COMBINED_SSU") sequences written to $COMBINED_SSU"
fi

if [[ ! -f "$COMBINED_TAX" ]]; then
    echo "=== Building combined taxonomy table ==="
    zcat "$BAC_TAX" > "$COMBINED_TAX"
    zcat "$AR_TAX" | tail -n +2 >> "$COMBINED_TAX"
    echo "  $(wc -l < "$COMBINED_TAX") entries written to $COMBINED_TAX"
fi

# --- Step 2: Build BLAST database (once) ---
if [[ ! -f "${BLAST_DB}.nhr" ]]; then
    echo ""
    echo "=== Building BLAST database ==="
    makeblastdb \
        -in "$COMBINED_SSU" \
        -dbtype nucl \
        -out "$BLAST_DB" \
        -title "GTDB_SSU" \
        > "$LOG_DIR/makeblastdb.log" 2>&1
    echo "  BLAST database built: $BLAST_DB"
fi

# --- Step 3: BLASTn ---
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

# --- Step 4: Parse results and assign taxonomy ---
python3 - "$BLAST_OUT" "$COMBINED_TAX" "$ZOTUS_FASTA" \
          "$RESULTS_DIR" "$BLAST_IDENTITY_THRESHOLD" "$BLAST_EVALUE_THRESHOLD" << 'PYEOF'
import sys
import csv

blast_file  = sys.argv[1]
tax_file    = sys.argv[2]
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

# Load taxonomy: accession -> lineage
tax = {}
with open(tax_file) as f:
    reader = csv.reader(f, delimiter='\t')
    next(reader)  # skip header
    for row in reader:
        if len(row) >= 2:
            tax[row[0].strip()] = row[1].strip()

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

    if hit['pident'] >= id_thresh and hit['evalue'] <= ev_thresh:
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
