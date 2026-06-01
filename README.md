# NanoZotu — Full-length 16S ZOTU pipeline for Oxford Nanopore

MSc thesis project (Asta). De novo ZOTU analysis of full-length 16S rRNA gene sequences (V1–V9, ~1500 bp) from a cultured algae microbiome, sequenced on Oxford Nanopore R10.4 with super accuracy basecalling. Produces GTDB-curated taxonomy and a Phyloseq R object via NanoASV.

---

## Sequencing overview

| Parameter | Value |
|---|---|
| Platform | Oxford Nanopore R10.4 |
| Basecalling | Super accuracy (SUP) |
| Target | Full-length 16S rRNA (~1500 bp, V1-V9) |
| Forward primer (27F) | `AGRGTTYGATYMTGGCTCAG` |
| Reverse primer (1492R) | `RGYTACCTTGTTACGACTT` |
| Samples | 21 barcodes |
| Reads per sample | 20,000 – 60,000 |

---

## Pipeline overview

```
01 → 02 → 03 → 04 → 05 → 06
                         ↓
                    08 → 09 → 10
```

| Script | Purpose |
|---|---|
| `01_concatenate_barcodes.sh` | Merge per-run fastq.gz files per barcode |
| `02_trim_primers.sh` | Primer trimming and size filtering |
| `03_dereplicate.sh` | Pool samples and dereplicate |
| `04_unoise3.sh` | UNOISE3 denoising across minsize thresholds |
| `05_otutab.sh` | Per-sample ZOTU abundance table |
| `06_taxonomy.sh` | BLASTn taxonomy against GTDB SSU + SILVA organellar |
| `08_build_gtdb_nanoasv.sh` | Build NanoASV base database: GTDB + SILVA organellar (one-time) |
| `09_augment_nanoasv_db.sh` | Inject novel ZOTUs into NanoASV database (one-time) |
| `10_nanoasv.sh` | NanoASV run → Phyloseq output |

Quick start:
```bash
# Sequential: merge chunks → trim → pool → denoise
bash 01_concatenate_barcodes.sh
bsub < 02_trim_primers.sh    # wait for completion
bsub < 03_dereplicate.sh     # wait for completion
bash 04_unoise3.sh

# Parallel: after step 04, three independent tasks can run simultaneously
bsub < 05_otutab.sh                                        # cluster: otutab
bash 06_taxonomy.sh pooled/zotus_minsize3.fasta            # interactive: BLAST (fast)
bash 08_build_gtdb_nanoasv.sh                              # interactive: DB build (run once)

# Converge: 09 requires both 06 (taxonomy_all.tsv) and 08 (GTDB NanoASV fasta)
bash 09_augment_nanoasv_db.sh

# Final cluster job (submit after 09 completes)
bsub < 10_nanoasv.sh                # NanoASV → Phyloseq
```

---

## Setup

### 1. Clone the repository

```bash
git clone git@github.com:jnesme/NanoZotu.git
cd NanoZotu
```

### 2. Install conda environments

**Main environment** (steps 01–07) — follow the QIIME 2 amplicon distribution
install instructions for your platform, then install additional tools:

```bash
# After creating qiime2-amplicon-2026.1 per QIIME 2 docs:
conda activate qiime2-amplicon-2026.1
conda install -c bioconda blast mafft fasttree
```

**usearch v12** — download the binary from https://drive5.com/usearch/ and place
it on your `PATH` as `usearch`.

**NanoASV environment** — follow the install instructions in the NanoASV repository:

```bash
git clone https://github.com/ImagoXV/NanoASV.git /path/to/NanoASV
# Follow NanoASV README to create its conda environment
```

### 3. Download GTDB SSU reference files

Download the following four files from the GTDB releases page
(https://gtdb.ecogenomics.org/downloads) for your chosen GTDB release
and place them in `db/gtdb/`:

```
db/gtdb/
├── bac120_ssu_reps.fna.gz       # bacterial SSU representative sequences
├── ar53_ssu_reps.fna.gz         # archaeal SSU representative sequences
├── bac120_taxonomy.tsv.gz       # bacterial taxonomy table
└── ar53_taxonomy.tsv.gz         # archaeal taxonomy table
```

```bash
mkdir -p db/gtdb
cd db/gtdb

# Replace RXX with the GTDB release number (e.g. R226)
wget https://data.gtdb.ecogenomics.org/releases/latest/genomic_files_reps/bac120_ssu_reps.fna.gz
wget https://data.gtdb.ecogenomics.org/releases/latest/genomic_files_reps/ar53_ssu_reps.fna.gz
wget https://data.gtdb.ecogenomics.org/releases/latest/bac120_taxonomy.tsv.gz
wget https://data.gtdb.ecogenomics.org/releases/latest/ar53_taxonomy.tsv.gz
cd ../..
```

### 4. Edit config.sh

Open `config.sh` and update the system paths for your environment:

```bash
CONDA_ENV_MAIN="qiime2-amplicon-2026.1"   # name of your QIIME 2 env
CONDA_ENV_NANOASV="NanoASV"               # name of your NanoASV env
NANOASV_PATH="/path/to/NanoASV"           # absolute path to NanoASV clone
```

Update primers, size filter, and thresholds if your amplicon target differs
from full-length 16S (V1–V9, 27F/1492R).

### 5. Edit the LSF script

In `10_nanoasv.sh`, set `PROJECT_DIR` to the absolute path of your project
clone and update the `#BSUB` headers for your HPC environment (queue name,
email, walltime):

```bash
# Change this line in 10_nanoasv.sh:
PROJECT_DIR="/absolute/path/to/NanoZotu"

# Update these headers for your HPC:
#BSUB -q your_queue
#BSUB -u your@email.com
```

### 6. Prepare raw data

Place basecaller output in `fastq_pass/` with one subdirectory per barcode:

```
fastq_pass/
├── barcode08/
│   ├── chunk_0.fastq.gz
│   └── chunk_1.fastq.gz
├── barcode16/
│   └── ...
```

The pipeline starts at step 01 which merges the per-barcode chunks into
`fastq_merged/barcode*.fastq.gz`.

---

## Configuration

All project-specific parameters live in **`config.sh`** in the project root. Edit this file before running any script on a new dataset.

```bash
# System
CONDA_ENV_MAIN="qiime2-amplicon-2026.1"   # conda env for steps 01–07
CONDA_ENV_NANOASV="NanoASV"               # conda env for step 10
NANOASV_PATH="/work3/josne/github/NanoASV"

# Primers (step 02)
PRIMER_FWD="AGRGTTYGATYMTGGCTCAG"         # 27F
PRIMER_REV="RGYTACCTTGTTACGACTT"          # 1492R (reference)
PRIMER_RC_REV="AAGTCGTAACAAGGTARCY"       # reverse complement of 1492R (cutadapt -a)

# Read QC (step 02)
MIN_LEN=1300
MAX_LEN=1800

# UNOISE3 (step 04)
MINSIZE_MIN=1
MINSIZE_MAX=8
MINSIZE_WORKING=3    # working threshold — propagates to steps 05, 06, 07, 09, 10

# Taxonomy (steps 06 and 09)
BLAST_IDENTITY_THRESHOLD=97
BLAST_EVALUE_THRESHOLD="1e-10"

# NanoASV (step 10)
NANOASV_MINAB=2
NANOASV_SUBSAMPLING=100000
NANOASV_SAM_QUAL=0
```

Scripts 01–09 source `config.sh` automatically relative to their own location. The LSF script (`10_nanoasv.sh`) additionally requires `PROJECT_DIR` to be set as an absolute path — this is the **only line that must be edited directly in that script**, since LSF copies it to `/tmp` before execution and relative paths are not reliable:

```bash
# In 10_nanoasv.sh — change this line:
PROJECT_DIR="/work3/josne/Projects/AstaMSc_GRF_Igalbana"
```

The `#BSUB` headers (queue, email, resources) in script 10 also require manual editing for your HPC environment.

---

## Dependencies

### Conda environments

| Environment | Used by |
|---|---|
| `qiime2-amplicon-2026.1` | Scripts 01–07 |
| `NanoASV` | Script 10 (activated inside bsub) |

Activate before running interactive scripts:
```bash
source /work3/josne/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2026.1
```

### Tools

- **cutadapt** >= 4.0 (required for `--revcomp`)
- **usearch** v12
- **blastn** (NCBI BLAST+), **mafft**, **FastTree**
- **minimap2**, **samtools**, **vsearch**, **chopper**, **porechop** (NanoASV env)

---

## Directory structure

```
.
├── fastq_pass/               # Raw basecaller output (one folder per barcode)
├── fastq_merged/             # One merged fastq.gz per barcode (step 01)
├── fastq_trimmed/            # Primer-trimmed reads (step 02)
├── pooled/                   # Pooled reads and ZOTU files (steps 03–04)
│   ├── all_samples.fastq
│   ├── all_samples_derep.fasta
│   ├── minsize_sweep.txt
│   └── zotus_minsize{1..8}.fasta
├── db/
│   ├── gtdb/                 # GTDB SSU reps + BLAST index + NanoASV FASTA
│   └── rrndb/                # rrnDB NanoASV FASTA (supplementary)
├── metadata/
│   └── nanoasv_metadata.csv  # 21 samples — update Measure_1/Measure_2 before NanoASV
├── results/
│   ├── taxonomy_zotus_minsize3/    # BLAST taxonomy results (step 06)
│   └── nanoasv/output/            # NanoASV outputs including Phyloseq (step 10)
├── supplementary/            # rrnDB scripts (future copy-number correction)
├── logs/
├── config.sh                 # ← edit this for each new project
├── 01_concatenate_barcodes.sh
├── 02_trim_primers.sh
├── 03_dereplicate.sh
├── 04_unoise3.sh
├── 05_otutab.sh
├── 06_taxonomy.sh
├── 08_build_gtdb_nanoasv.sh
├── 09_augment_nanoasv_db.sh
└── 10_nanoasv.sh
```

---

## Pipeline

### Step 01 — Concatenate barcode chunks

Merges the multiple fastq.gz chunk files produced by the basecaller into a single file per barcode.

```bash
bash 01_concatenate_barcodes.sh
```

**Input:** `fastq_pass/barcode*/` — multiple `*.fastq.gz` chunks per barcode  
**Output:** `fastq_merged/barcode{08..95}.fastq.gz`

---

### Step 02 — Primer trimming

Trims 16S primers using cutadapt. Retains only reads where both primers are found, with insert size between 1300 and 1800 bp. Uses `--revcomp` to handle reads in either orientation in a single pass — all output reads are in forward orientation.

```bash
bash 02_trim_primers.sh

# Override thread count (default: all available CPUs)
THREADS=16 bash 02_trim_primers.sh
```

**Input:** `fastq_merged/barcode*.fastq.gz`  
**Output:** `fastq_trimmed/barcode*.fastq.gz`  
**Logs:** `logs/cutadapt/barcode*.log`

Key parameters:
- `--revcomp` — single-pass handling of both strand orientations (requires cutadapt >= 4.0)
- `--discard-untrimmed` — both primers must be found
- `--error-rate 0.15` — tolerant of ONT read errors and IUPAC degenerate bases
- `--overlap 15` — minimum primer overlap to avoid false matches
- `-m 1300 -M 1800` — size filter on trimmed insert

---

### Step 03 — Pool and dereplicate

Pools all barcode reads into a single FASTQ (tagging each read header with `sample=barcodeXX` for downstream sample assignment), then dereplicates the pooled dataset. Pooling before dereplication is required so that ZOTUs are defined consistently across all samples.

```bash
bsub < 03_dereplicate.sh
```

**Input:** `fastq_trimmed/barcode*.fastq.gz`  
**Output:**
- `pooled/all_samples.fastq` — pooled reads, kept for ZOTU table mapping in step 05
- `pooled/all_samples_derep.fasta` — dereplicated sequences with `size=N` abundance annotations

**Logs:** `logs/derep/derep_pooled.log`

Key parameters:
- `-minuniquesize 2` — discards singletons
- `-sizeout` — adds `;size=N` to headers, required by UNOISE3

> Note: ~99.9% of unique sequences per sample are singletons — expected for ONT full-length 16S due to error spreading across 1500 bp reads.

---

### Step 04 — UNOISE3 denoising (minsize sweep)

Runs UNOISE3 across a range of `-minsize` thresholds (1 to 8) and reports ZOTU counts for each. One ZOTU FASTA is saved per threshold for inspection.

```bash
bash 04_unoise3.sh
```

**Input:** `pooled/all_samples_derep.fasta`  
**Output:**
- `pooled/zotus_minsize{1..8}.fasta` — denoised ZOTUs per threshold
- `pooled/minsize_sweep.txt` — summary table of ZOTU counts

**Logs:** `logs/unoise3/unoise3_minsize{1..8}.log`

Example output for this dataset:

```
minsize    ZOTUs
-------    -----
1          272
2          272
3          14
4          11
5          9
6          9
7          7
8          7
```

Key observations:
- `minsize=1` and `minsize=2` yield identical results: UNOISE3's error model implicitly rejects all singletons regardless of the minsize pre-filter
- The cliff between minsize=2 (272) and minsize=3 (14) reflects systematic ONT error pairs — sequences seen exactly twice due to recurring error patterns, not true biology
- `minsize=3` (14 ZOTUs) is the working threshold for this dataset: biologically plausible for this low-diversity system
- Some ZOTUs represent intragenomic 16S copy variants from the same organism
- Several ZOTUs match full-length reference sequences at 100% identity, validating pipeline accuracy

---

### Step 05 — ZOTU abundance table

Maps pooled reads back to ZOTUs and generates a per-sample abundance table. The ZOTU FASTA is passed as an argument, allowing easy comparison across minsize thresholds.

```bash
bash 05_otutab.sh pooled/zotus_minsize3.fasta

# Override thread count (default: all available CPUs)
THREADS=16 bash 05_otutab.sh pooled/zotus_minsize3.fasta
```

**Input:**
- `pooled/all_samples.fastq` — pooled reads (fixed)
- `pooled/zotus_minsize3.fasta` — ZOTU sequences (user-specified)

**Output:** `results/zotu_table_zotus_minsize3.txt` — samples × ZOTUs count matrix  
**Logs:** `logs/otutab/otutab.log`

> Per-sample assignment is based on the `sample=barcodeXX` tag added to read headers during pooling in step 03.

---

---

### Step 06 — Taxonomy assignment

BLASTn of ZOTUs against the NanoASV base database (GTDB SSU + SILVA organellar).
Best hit is selected by bitscore across 500 candidate hits — not by the first BLAST result,
which is not guaranteed to be the best (Shah et al. 2019). Organellar hits (Chloroplast,
Mitochondria) are always assigned regardless of pident — any organellar match is informative.

```bash
bash 06_taxonomy.sh pooled/zotus_minsize3.fasta
```

**Prerequisite:** `bash 08_build_gtdb_nanoasv.sh` (builds the reference database)  
**Input:** `<zotus.fasta>`, `db/gtdb/SINGLELINE_GTDB_SSU_nanoasv.fasta`  
**Output:** `results/taxonomy_<zotus>/taxonomy_assigned.tsv`, `taxonomy_unknown.tsv`, `taxonomy_all.tsv`

| Parameter | Value | Rationale |
|---|---|---|
| pident threshold | ≥ 97% | Standard species-level 16S identity cutoff (bacteria) |
| evalue threshold | ≤ 1e-10 | Eliminates spurious low-complexity matches |
| max_target_seqs | 500 | Ensures true best hit is found before bitscore selection |

Results (minsize=3, 14 ZOTUs):

**Bacterial (pident ≥ 97%, 10 ZOTUs):**

| ZOTU | Species | pident |
|---|---|---|
| Zotu3 | Sulfitobacter pontiacus | 100% |
| Zotu4 | Roseivirga pacifica | 100% |
| Zotu5, 7, 8, 11 | Phaeobacter piscinae | 100% — intragenomic 16S copies |
| Zotu6 | Lentilitoribacter sp. | 99.6% |
| Zotu9 | Alteromonas marina_A | 100% |
| Zotu10 | Alteromonas abrolhosensis | 99.8% |
| Zotu13 | Alteromonas sp. | 99.9% |

**Host chloroplast (4 ZOTUs):** Zotu1, 2, 12, 14 are I. galbana chloroplast 16S (confirmed
100% identity to NC_049168.1 via NCBI BLAST). The 4 ZOTUs correspond to the 4 16S rRNA gene
copies in the chloroplast genome. These reads are quantified separately — they are not part
of the bacterial microbiome.

> Note: GTDB contains no organellar sequences. Without the SILVA organellar pre-filter, these
> ZOTUs receive a misleading bacterial hit at ~88% identity. The unified GTDB + SILVA organellar
> database ensures correct assignment in any culture system.

---

---

### NanoASV branch — overview (steps 08–10)

Steps 08–10 produce a Phyloseq R object from the same reads, using NanoASV's Snakemake
workflow. The approach combines the de novo ZOTUs from step 04 with reference-based read
classification:

1. **Step 08** builds the NanoASV base database: GTDB SSU representative sequences (~93k
   bacterial + archaeal) plus the SILVA organellar subset (Chloroplast + Mitochondria lineages,
   ~5.8k sequences). This is the reference that minimap2 maps reads against.

2. **Step 09** injects novel project ZOTUs — those whose exact sequence is absent from the base
   database (pident < 100%) — labelled with their step 06 taxonomy. This is critical: without it,
   minimap2 maps reads from novel organisms to whatever distant reference produces the best local
   alignment score, giving wrong taxonomy. With the ZOTUs in the database, those reads map to
   their own ZOTU consensus at ~100% identity, outcompeting any spurious hit.

3. **Step 10** runs NanoASV against the augmented database and exports a Phyloseq Rdata object.

Steps 08 and 09 are one-time setup per project. Step 09 depends on step 06 having been run.

---

### Step 08 — Build NanoASV base database

Builds the NanoASV base database: GTDB SSU representative sequences (~93k bacterial + archaeal)
plus SILVA organellar sequences (~5.8k Chloroplast and Mitochondria lineages). Run once before
the first NanoASV job; re-run when GTDB or SILVA is updated.

```bash
bash 08_build_gtdb_nanoasv.sh
```

**Input:** `db/gtdb/bac120_ssu_reps.fna.gz`, `db/gtdb/ar53_ssu_reps.fna.gz`,
`db/silva/SINGLELINE_SILVA_138.2_plus_zotus.fasta`  
**Output:** `db/gtdb/SINGLELINE_GTDB_SSU_nanoasv.fasta` (uncompressed)

GTDB SSU headers carry rank-prefixed taxonomy (`d__Bacteria;p__Pseudomonadota;...`) and
bracketed metadata (`[locus_tag=...]`). The script strips rank prefixes and discards bracketed
fields. SILVA organellar sequences are already in NanoASV-compatible format and are appended
as-is. Output is uncompressed because NanoASV's format validation uses plain `grep`.

---

### Step 09 — Augment NanoASV database with project ZOTUs

Appends the project ZOTUs to the GTDB NanoASV database, labelling each with the
GTDB taxonomy string produced by step 06. **This step must run after step 06.**

The BLAST result from step 06 (`taxonomy_all.tsv`) is the taxonomy source: each ZOTU
is written into the NanoASV database with exactly the same GTDB label it received from
BLAST. This ensures that NanoASV assigns consistent taxonomy whether reads are classified
via the broad GTDB reference or via a ZOTU entry.

```bash
bash 09_augment_nanoasv_db.sh
```

**Input:**
- `db/gtdb/SINGLELINE_GTDB_SSU_nanoasv.fasta` — base database from step 08
- `pooled/zotus_minsize3.fasta` — ZOTU sequences from step 04
- `results/taxonomy_zotus_minsize3/taxonomy_all.tsv` — GTDB taxonomy from step 06

**Output:** `db/gtdb/SINGLELINE_GTDB_SSU_plus_zotus.fasta`

Taxonomy labelling rules (mirroring the 97% threshold from step 06):
- pident = 100%: skip — exact sequence is already in the base database
- pident ≥ 98.7%: full taxonomy string (species-level, Kim et al. 2014)
- pident ≥ 97%: genus-level only → `Genus sp.`
- pident < 97%: family-level only → `unclassified_<family>`

Organellar ZOTUs are handled automatically: if their sequence is present in the
SILVA organellar subset at 100%, they are skipped; if absent or diverged, they
are injected with their SILVA organellar taxonomy.

---

### Step 10 — NanoASV (Phyloseq output)

Runs the NanoASV Snakemake workflow: chopper quality filter → porechop adapter trim →
subsampling → minimap2 mapping against GTDB+ZOTUs → vsearch clustering of unmatched reads →
MAFFT+FastTree phylogeny → Phyloseq R export.

```bash
bsub < 10_nanoasv.sh
```

**Input:** `fastq_merged/barcodeXX.fastq.gz`, `metadata/nanoasv_metadata.csv`, `db/gtdb/SINGLELINE_GTDB_SSU_plus_zotus.fasta`  
**Output:** `results/nanoasv/output/Results/` — Phyloseq Rdata, per-sample taxonomy CSVs, abundance tables, phylogenetic tree

| Parameter | Value | Rationale |
|---|---|---|
| `--subsampling` | 100000 reads | Cap per-sample depth for comparability |
| `--minab` | 2 | Below default 5; appropriate for low-depth samples |
| `--sam-qual` | 0 | Accept all mapping quality levels |
| `--no-r-cleaning` | enabled | GTDB taxonomy keywords differ from SILVA; skip to avoid false filtering |

**Before running:** update `metadata/nanoasv_metadata.csv` — the `Measure_1` and
`Measure_2` columns are placeholders. Replace with real experimental variables
(timepoint, treatment, etc.) before interpreting Phyloseq output.

---

## Design notes

### Why --revcomp instead of two passes
ONT reads are sequenced from both strands at roughly equal frequency. A naive cutadapt run with `-g FWD -a RC_REV` would silently discard all reverse-complement reads (~50% of valid reads). `--revcomp` handles both orientations in one pass and normalises all output to forward orientation.

### Why pool before dereplicating
If dereplication and UNOISE3 are run per-sample, each sample produces an independent ZOTU set with no shared identifiers — a per-sample abundance table cannot be constructed. Pooling first ensures ZOTUs are defined globally across all samples.

### Why not DADA2
DADA2 uses per-base quality scores to build a substitution error model. ONT errors are dominated by indels (especially in homopolymers), not substitutions, making the DADA2 error model a poor fit. UNOISE3's abundance-ratio error model is platform-agnostic and orders of magnitude faster.

### Why not EMU
EMU is a closed-reference tool: reads are distributed across reference taxa via expectation-maximisation. It detects more taxa (including rare ones) but produces no sequences — ZOTUs cannot be extracted, BLASTed, or mapped to a genome of interest. For this project, where genome mapping is a downstream goal, UNOISE3 is the appropriate choice. EMU also cannot detect organisms absent from its reference database — novel lineages at low identity would be silently misassigned.

### Why GTDB over SILVA for taxonomy
SILVA species-level taxonomy is not curated: many entries carry names like "uncultured bacterium"
or have uncertain assignments. GTDB taxonomy is phylogenomically curated using genome-based methods,
providing reliable species resolution for cultivated organisms and a principled framework for
naming novel lineages. Using GTDB throughout (BLAST in step 06, NanoASV database in steps 08–10)
ensures consistent taxonomy across the entire pipeline.

### Why add ZOTUs to the NanoASV database
NanoASV uses minimap2, which does not refuse to map reads at low identity — it assigns every read
to whatever reference produces the best local alignment score. A novel organism at 88% identity
to any database entry will map spuriously to some unrelated reference rather than being flagged
as unknown. Adding the 14 project ZOTUs as explicit database entries ensures reads map to
their correct consensus sequences at ~100% identity, completely outcompeting any alternative hit.

### Sequencing depth and rare taxa
At 20–60k reads per sample, rare community members fall below the detection threshold for UNOISE3.
Riisgaard-Jensen et al. (2026) showed that ONT requires 4.1–5.6× more reads than PacBio for
V1-V8 amplicons to fully resolve communities. The 14 ZOTUs recovered at minsize=3 represent
the dominant community; rare taxa are not captured at this depth.

---

## Experimental: savont comparison

`test_savont.sh` and `test_savont_pooled.sh` benchmark [savont](https://github.com/bluenote-1577/savont) as a potential replacement for the UNOISE3 denoising step (step 04). savont uses SNPmer-based clustering designed for long-read indel error profiles, as opposed to UNOISE3's substitution-focused error model.

### Results

| Mode | ASVs | Notes |
|---|---|---|
| Per-sample (`test_savont.sh`) | 427 | Artifact — independent per-sample runs produce overlapping but inconsistent ASV sets that inflate counts on merge ([savont issue #2](https://github.com/bluenote-1577/savont/issues/2)) |
| Pooled (`test_savont_pooled.sh`) | 29 | Biologically plausible; consistent with UNOISE3's 14 ZOTUs + rare tail |

The pooled run also provided the first quantitative estimate of chloroplast read fraction: **51% of pooled reads map to I. galbana chloroplast** at this sequencing depth.

### Proposed pipeline integration (not yet implemented)

The natural integration preserves the sample-tagging strategy already in place:

1. **Step 03 unchanged** — `pooled/all_samples.fastq` with `sample=barcodeXX` read headers
2. **Replace step 04** — `savont asv pooled/all_samples.fastq --output-dir results/savont_pooled --fl-16s`
3. **Step 05 unchanged** — `usearch -otutab pooled/all_samples.fastq -zotus results/savont_pooled/final_asvs.fasta`

The `sample=` tags survive through savont (it does not strip read headers), so usearch `-otutab` can still demultiplex per-sample abundances from the mapping. Steps 06–10 are unaffected.

Next step: re-run `test_savont_pooled.sh` using `pooled/all_samples.fastq` (tagged reads from step 03) as input instead of raw trimmed files, then run `05_otutab.sh` against `final_asvs.fasta` and compare the per-sample table to the UNOISE3 result.

---

## Supplementary scripts

`supplementary/` contains scripts for a potential future analysis: per-species 16S copy number
correction using rrnDB. rrnDB tracks copy numbers per species, allowing read counts to be
normalised to genome equivalents. This has not been implemented and would require integrating
rrnDB taxonomy with the GTDB-based pipeline.

- `12_build_rrndb_nanoasv.sh` — converts rrnDB 5.10 to NanoASV singleline FASTA
- `11b_nanoasv_rrndb.sh` — NanoASV run using rrnDB as reference

---

## Key reference

Riisgaard-Jensen et al. 2026 — *Nanopore sequencing reaches amplicon sequence variant (ASV) resolution*  
bioRxiv: https://www.biorxiv.org/content/10.64898/2026.02.26.708165v1
