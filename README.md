# NanoZotu — Full-length 16S ZOTU pipeline for Oxford Nanopore

MSc thesis project (Asta). De novo ZOTU analysis of full-length 16S rRNA gene sequences (V1–V9, ~1500 bp) from a cultured algae microbiome, sequenced on Oxford Nanopore R10.4 with super accuracy basecalling. Produces GTDB-curated taxonomy, phylogenetic placement of a novel Elusimicrobiota lineage, and a Phyloseq R object via NanoASV.

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
01 → 02 → 03 → 04 → 05 → 06 → 07
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
| `06_taxonomy.sh` | BLASTn taxonomy against GTDB SSU |
| `07_elusimicrobiota_tree.sh` | Phylogenetic placement of novel Elusimicrobiota ZOTUs |
| `08_build_gtdb_nanoasv.sh` | Convert GTDB SSU to NanoASV format (one-time) |
| `09_augment_nanoasv_db.sh` | Append project ZOTUs to GTDB NanoASV database (one-time) |
| `10_nanoasv.sh` | NanoASV run → Phyloseq output |

Quick start:
```bash
# Step 1: build GTDB NanoASV base database (needs only db/gtdb/*.fna.gz — run once)
bash 08_build_gtdb_nanoasv.sh

# Step 2: de novo ZOTU pipeline
bash 01_concatenate_barcodes.sh
bash 02_trim_primers.sh
bash 03_dereplicate.sh
bash 04_unoise3.sh
bash 05_otutab.sh pooled/zotus_minsize3.fasta
bash 06_taxonomy.sh pooled/zotus_minsize3.fasta   # produces taxonomy_all.tsv

# Step 3: inject ZOTUs + their BLAST taxonomy into the NanoASV database (depends on step 06)
bash 09_augment_nanoasv_db.sh

# Step 4: optional phylogenetic placement of novel ZOTUs
bsub < 07_elusimicrobiota_tree.sh

# Step 5: NanoASV → Phyloseq (cluster job)
bsub < 10_nanoasv.sh
```

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
│   ├── taxonomy_zotus_minsize3/    # GTDB BLAST results (step 06)
│   ├── elusimicrobiota_tree/       # Elusimicrobiota phylogeny (step 07)
│   └── nanoasv/output/            # NanoASV outputs including Phyloseq (step 10)
├── supplementary/            # rrnDB scripts (future copy-number correction)
├── logs/
├── 01_concatenate_barcodes.sh
├── 02_trim_primers.sh
├── 03_dereplicate.sh
├── 04_unoise3.sh
├── 05_otutab.sh
├── 06_taxonomy.sh
├── 07_elusimicrobiota_tree.sh
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
bash 03_dereplicate.sh
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

### Step 06 — GTDB taxonomy assignment

BLASTn of ZOTUs against the GTDB SSU representative sequences (bac120 + ar53 combined).
Best hit is selected by bitscore across 500 candidate hits — not by the first BLAST result,
which is not guaranteed to be the best (Shah et al. 2019).

```bash
bash 06_taxonomy.sh pooled/zotus_minsize3.fasta
```

**Input:** `<zotus.fasta>`, `db/gtdb/` (BLAST index auto-built on first run)  
**Output:** `results/taxonomy_<zotus>/taxonomy_assigned.tsv`, `taxonomy_unknown.tsv`, `taxonomy_all.tsv`

| Parameter | Value | Rationale |
|---|---|---|
| pident threshold | ≥ 97% | Standard species-level 16S identity cutoff |
| evalue threshold | ≤ 1e-10 | Eliminates spurious low-complexity matches |
| max_target_seqs | 500 | Ensures true best hit is found before bitscore selection |

Results (minsize=3, 14 ZOTUs):

**Assigned (pident ≥ 97%, 10 ZOTUs):**

| ZOTU | Species | pident |
|---|---|---|
| Zotu3 | Sulfitobacter pontiacus | 100% |
| Zotu4 | Roseivirga pacifica | 100% |
| Zotu5, 7, 8, 11 | Phaeobacter piscinae | 100% — intragenomic 16S copies |
| Zotu6 | Lentilitoribacter sp. | 99.6% |
| Zotu9 | Alteromonas marina_A | 100% |
| Zotu10 | Alteromonas abrolhosensis | 99.8% |
| Zotu13 | Alteromonas sp. | 99.9% |

**Novel (pident ~88%, 4 ZOTUs):** Zotu1, 2, 12, 14 all hit GWA2-66-18 sp965283875
(phylum Elusimicrobiota, order UBA1565, family UBA9628). At 88% identity these represent
a novel genus/species — see design decisions below.

---

### Step 07 — Elusimicrobiota phylogenetic tree

Builds a maximum-likelihood phylogenetic tree of all 271 GTDB SSU representatives from
phylum Elusimicrobiota, the 4 novel project ZOTUs, and an archaeal outgroup (MAFFT + FastTree).

```bash
bsub < 07_elusimicrobiota_tree.sh
```

**Input:** GTDB SSU FASTA, `pooled/zotus_minsize3.fasta`  
**Output:** `results/elusimicrobiota_tree/` (FASTA, alignment, tree)

The 4 ZOTUs form a tight monophyletic clade with GWA2-66-18 (the closest GTDB reference),
confirming order-level placement within UBA1565. This is independent validation of the BLAST
result and demonstrates the novelty is real divergence, not a pipeline artifact.

---

---

### NanoASV branch — overview (steps 08–10)

Steps 08–10 produce a Phyloseq R object from the same reads, using NanoASV's Snakemake
workflow. The approach combines the de novo ZOTUs from step 04 with reference-based read
classification:

1. **Step 08** converts the GTDB SSU representative sequences into the singleline FASTA
   format required by NanoASV. This is the reference database that minimap2 will map reads against.

2. **Step 09** takes the ZOTUs from step 04 and injects them — labelled with their step 06
   GTDB taxonomy — into the database built in step 08. This is the critical step: without it,
   minimap2 maps reads from novel organisms (e.g. Elusimicrobiota at 88% identity) to whatever
   distant reference produces the best local alignment score, giving wrong taxonomy. With the
   ZOTUs in the database, those reads map to their own ZOTU at ~100% identity, outcompeting
   any spurious hit.

3. **Step 10** runs NanoASV against the augmented database and exports a Phyloseq Rdata object.

Steps 08 and 09 are one-time setup per project. Step 09 depends on step 06 having been run.

---

### Step 08 — Build GTDB NanoASV database

Converts GTDB SSU representative sequences (~93k bacterial + archaeal) into the singleline
FASTA format required by NanoASV. Run once before the first NanoASV job.

```bash
bash 08_build_gtdb_nanoasv.sh
```

**Input:** `db/gtdb/bac120_ssu_reps.fna.gz`, `db/gtdb/ar53_ssu_reps.fna.gz`  
**Output:** `db/gtdb/SINGLELINE_GTDB_SSU_nanoasv.fasta` (~129 MB, uncompressed)

GTDB SSU headers carry rank-prefixed taxonomy (`d__Bacteria;p__Pseudomonadota;...`) and
bracketed metadata (`[locus_tag=...]`). The script strips rank prefixes and discards bracketed
fields. Output is uncompressed because NanoASV's format validation uses plain `grep`.

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
- pident ≥ 97%: full GTDB taxonomy string (known organism)
- pident < 97%: confident to family only; genus+species → `unclassified_<family>`
  (Zotu1/2/12/14 at ~88% → `unclassified_UBA9628`)

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
EMU is a closed-reference tool: reads are distributed across reference taxa via expectation-maximisation. It detects more taxa (including rare ones) but produces no sequences — ZOTUs cannot be extracted, BLASTed, or mapped to a genome of interest. For this project, where genome mapping is a downstream goal, UNOISE3 is the appropriate choice. More critically, EMU cannot detect organisms absent from its reference database — the novel Elusimicrobiota (88% identity to the nearest GTDB entry) would be undetectable by any closed-reference pipeline.

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

### Elusimicrobiota taxonomy: why `unclassified_UBA9628`
The four Elusimicrobiota ZOTUs hit GWA2-66-18 sp965283875 at ~88% 16S identity.
Standard thresholds: ≥98.7% = same species, ≥94.5% = same genus, ≥86.5% = same family.
At 88%, we are confident the organism belongs to family UBA9628 but genus and species are
genuinely novel. Labelling these ZOTUs as sp965283875 would be factually wrong — GWA2-66-18
is the closest known relative, not the organism itself. `unclassified_UBA9628` correctly
represents the resolution supported by the data.

### Sequencing depth and rare taxa
At 20–60k reads per sample, rare community members fall below the detection threshold for UNOISE3.
Riisgaard-Jensen et al. (2026) showed that ONT requires 4.1–5.6× more reads than PacBio for
V1-V8 amplicons to fully resolve communities. The 14 ZOTUs recovered at minsize=3 represent
the dominant community; rare taxa are not captured at this depth.

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
