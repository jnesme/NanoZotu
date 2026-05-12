#!/usr/bin/env bash
# config.sh — pipeline parameters
#
# Edit this file to adapt the pipeline to a new project.
# All scripts source this file automatically.
#
# NOTE: LSF job headers (#BSUB lines) in 07_elusimicrobiota_tree.sh and
# 10_nanoasv.sh cannot use shell variables and must be edited directly.
# LSF scripts also require PROJECT_DIR set to an absolute path — that is
# the one hardcoded value that remains in those scripts.

# =============================================================================
# System
# =============================================================================

# Conda environment for steps 01–07 (read QC, UNOISE3, taxonomy, tree)
CONDA_ENV_MAIN="qiime2-amplicon-2026.1"

# Conda environment for step 10 (NanoASV)
CONDA_ENV_NANOASV="NanoASV"

# Path to the NanoASV installation
NANOASV_PATH="/work3/josne/github/NanoASV"

# =============================================================================
# Primers  (step 02 — cutadapt)
# =============================================================================

PRIMER_FWD="AGRGTTYGATYMTGGCTCAG"    # 27F forward primer
PRIMER_REV="RGYTACCTTGTTACGACTT"     # 1492R reverse primer (reference only)
PRIMER_RC_REV="AAGTCGTAACAAGGTARCY"  # reverse complement of PRIMER_REV (cutadapt -a)

# =============================================================================
# Read QC  (step 02 — size filter on trimmed insert)
# =============================================================================

MIN_LEN=1300
MAX_LEN=1800

# =============================================================================
# UNOISE3 denoising  (step 04)
# =============================================================================

MINSIZE_MIN=1           # lower bound of minsize sweep
MINSIZE_MAX=8           # upper bound of minsize sweep
MINSIZE_WORKING=3       # working threshold — used by steps 05, 06, 07, 09, 10

# =============================================================================
# Taxonomy  (steps 06 and 09)
# =============================================================================

BLAST_IDENTITY_THRESHOLD=97    # minimum % identity for genus-level assignment
BLAST_SPECIES_THRESHOLD=98.7   # minimum % identity for confident species-level assignment (Kim et al. 2014)
BLAST_EVALUE_THRESHOLD="1e-10" # maximum e-value for confident assignment

# =============================================================================
# NanoASV  (step 10)
# =============================================================================

NANOASV_MINAB=2           # minimum read abundance per taxon
NANOASV_SUBSAMPLING=100000  # reads per sample cap
NANOASV_SAM_QUAL=0        # minimum mapping quality (0 = accept all)
