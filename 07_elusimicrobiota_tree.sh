#!/bin/bash
### General options
#BSUB -q hpcspecial
#BSUB -J elusimicrobiota_tree
#BSUB -n 20
#BSUB -R "span[hosts=1] rusage[mem=8GB]"
#BSUB -M 8500MB
#BSUB -W 4:00
#BSUB -u jnesme@gmail.com
#BSUB -B
#BSUB -N
#BSUB -o logs/elusimicrobiota_tree_%J.out
#BSUB -e logs/elusimicrobiota_tree_%J.err

# Phylogenetic tree of Elusimicrobiota to place novel ZOTUs within the phylum.
#
# Sequences:
#   - All GTDB SSU representative sequences from phylum Elusimicrobiota (271 seqs)
#   - 4 novel ZOTUs (Zotu1, Zotu2, Zotu12, Zotu14)
#   - Archaeal outgroup (from GTDB ar53 reps)
#
# Usage: bsub < 10_elusimicrobiota_tree.sh

set -euo pipefail

# PROJECT_DIR must be an absolute path — LSF copies this script to /tmp before
# execution, so relative paths and ${BASH_SOURCE[0]} are not reliable.
PROJECT_DIR="/work3/josne/github/NanoZotu"
source "$PROJECT_DIR/config.sh"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV_MAIN"

THREADS=$LSB_DJOB_NUMPROC

GTDB_TAX="$PROJECT_DIR/db/gtdb/gtdb_taxonomy_combined.tsv"
GTDB_SSU="$PROJECT_DIR/db/gtdb/bac120_ssu_reps.fna.gz"
GTDB_AR_SSU="$PROJECT_DIR/db/gtdb/ar53_ssu_reps.fna.gz"
ZOTUS_FASTA="$PROJECT_DIR/pooled/zotus_minsize${MINSIZE_WORKING}.fasta"

TREE_DIR="$PROJECT_DIR/results/elusimicrobiota_tree"
LOG_DIR="$PROJECT_DIR/logs/elusimicrobiota_tree"

mkdir -p "$TREE_DIR" "$LOG_DIR"

SEQS_FASTA="$TREE_DIR/sequences_for_tree.fasta"
ALIGNMENT="$TREE_DIR/alignment.fasta"
TREE="$TREE_DIR/tree.nwk"

# --- Step 1: Extract Elusimicrobiota accessions from GTDB taxonomy ---
echo "=== Extracting Elusimicrobiota sequences from GTDB SSU reps ==="

python3 - "$GTDB_TAX" "$GTDB_SSU" "$GTDB_AR_SSU" "$ZOTUS_FASTA" "$SEQS_FASTA" << 'PYEOF'
import sys
import csv
import gzip
import re

tax_file    = sys.argv[1]
gtdb_ssu    = sys.argv[2]
gtdb_ar_ssu = sys.argv[3]
zotus_fasta = sys.argv[4]
out_fasta   = sys.argv[5]

# Elusimicrobiota ZOTUs to extract from the ZOTU file
ELUSIMICROBIOTA_ZOTUS = {'Zotu1', 'Zotu2', 'Zotu12', 'Zotu14'}

# Get Elusimicrobiota accessions from GTDB taxonomy
elusimicro_accessions = {}   # accession -> species label
with open(tax_file) as f:
    reader = csv.reader(f, delimiter='\t')
    next(reader)  # skip header
    for row in reader:
        if len(row) < 2:
            continue
        accession = row[0].strip()
        lineage   = row[1].strip()
        if 'elusimicro' in lineage.lower():
            species = re.search(r's__(.+)$', lineage)
            species = species.group(1).replace(' ', '_') if species else 'unknown'
            elusimicro_accessions[accession] = species

print(f"  Elusimicrobiota accessions in taxonomy: {len(elusimicro_accessions)}")

# Extract matching sequences from GTDB bacterial SSU reps
found = {}
current_acc = None
current_seq = []

def flush(found, current_acc, current_seq, elusimicro_accessions):
    if current_acc and current_acc in elusimicro_accessions:
        found[current_acc] = ''.join(current_seq)

with gzip.open(gtdb_ssu, 'rt') as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            flush(found, current_acc, current_seq, elusimicro_accessions)
            current_acc = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    flush(found, current_acc, current_seq, elusimicro_accessions)

print(f"  Sequences extracted from GTDB SSU reps: {len(found)}")

# Write Elusimicrobiota sequences
with open(out_fasta, 'w') as out:
    for acc, seq in sorted(found.items()):
        species = elusimicro_accessions[acc]
        out.write(f">GTDB_{acc}_{species}\n{seq}\n")

    # Add the 4 Elusimicrobiota ZOTUs
    with open(zotus_fasta) as f:
        in_target = False
        seq_lines = []
        header    = None
        for line in f:
            if line.startswith('>'):
                if header and header.lstrip('>').split()[0] in ELUSIMICROBIOTA_ZOTUS:
                    out.write(f"{header}{(''.join(seq_lines))}\n")
                header    = line
                seq_lines = []
            else:
                seq_lines.append(line.rstrip())
        if header and header.lstrip('>').split()[0] in ELUSIMICROBIOTA_ZOTUS:
            out.write(f"{header}{''.join(seq_lines)}\n")

    # Add archaeal outgroup (first sequence from GTDB ar53 reps)
    with gzip.open(gtdb_ar_ssu, 'rt') as f:
        header = None
        seq    = []
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if header:
                    break
                header = line
            else:
                seq.append(line)
        out.write(f">OUTGROUP_Archaea\n{''.join(seq)}\n")

print(f"  Added 4 Elusimicrobiota ZOTUs + archaeal outgroup")
print(f"  Total sequences: {len(found) + 4 + 1}")
print(f"  Written to: {out_fasta}")
PYEOF

# --- Step 2: Align with MAFFT (G-INS-i) ---
echo ""
echo "=== Aligning with MAFFT (G-INS-i) ==="

mafft \
    --globalpair \
    --maxiterate 16 \
    --thread "$THREADS" \
    "$SEQS_FASTA" \
    > "$ALIGNMENT" \
    2> "$LOG_DIR/mafft.log"

n_aligned=$(grep -c '^>' "$ALIGNMENT")
echo "  Aligned $n_aligned sequences."

# --- Step 3: Build tree with FastTree (GTR) ---
echo ""
echo "=== Building tree with FastTree (GTR) ==="

FastTree \
    -nt \
    -gtr \
    -log "$LOG_DIR/fasttree.log" \
    "$ALIGNMENT" \
    > "$TREE" \
    2>> "$LOG_DIR/fasttree.log"

echo "  Tree written to: $TREE"
echo ""
echo "Done. Results in: $TREE_DIR"
ls -lh "$TREE_DIR"
