# CLAUDE.md — Project Context for Claude Code

## Project
MSc thesis project (Asta). Full-length 16S rRNA amplicon sequencing of the microbiome associated with a cultured algae. Low-to-moderate diversity environment.

## Sequencing
- Platform: Oxford Nanopore, R10.4 flow cells
- Basecalling: Super accuracy (SUP) model
- Target: Full-length 16S (~1500 bp, V1-V9)
- Primers: FWD 27F (AGRGTTYGATYMTGGCTCAG), REV 1492R (RGYTACCTTGTTACGACTT)
- 21 barcoded samples (barcode08 to barcode95, non-contiguous)

## Environment
- Conda: qiime2-amplicon-2026.1 (steps 01–06), NanoASV (step 10)
- Key tools: cutadapt (>=4.0 for --revcomp), usearch v12

## Configuration
All project-specific parameters are in config.sh (project root). Scripts source it automatically.
Key parameters: PRIMER_FWD/REV/RC_REV, MIN_LEN/MAX_LEN, MINSIZE_MIN/MAX/WORKING,
BLAST_IDENTITY_THRESHOLD, BLAST_EVALUE_THRESHOLD, NANOASV_MINAB/SUBSAMPLING/SAM_QUAL,
CONDA_ENV_MAIN, CONDA_ENV_NANOASV, NANOASV_PATH.
PROJECT_DIR is defined in config.sh (single source of truth). LSF scripts source config.sh via a hardcoded absolute path (unavoidable — LSF copies scripts to /tmp).
#BSUB headers (queue, email) must also be edited directly in each LSF script.

## Pipeline overview
01 → 02 → 03 → 04 → 05 → 06
                         ↓
                    08 → 09 → 10

- 01–06: core UNOISE3 pipeline (read processing → taxonomy)
- 08–10: NanoASV branch (GTDB SSU + SILVA organellar + ZOTUs → Phyloseq output)
  - 08: build NanoASV base DB (GTDB + SILVA organellar, run once, interactive)
  - 09: inject novel ZOTUs into base DB (run once, interactive)
  - 10: NanoASV cluster job → Phyloseq Rdata

supplementary/ contains rrnDB scripts (future copy-number correction).

All scripts are run from the project root: /work3/josne/Projects/AstaMSc_GRF_Igalbana/

## Key decisions and rationale

### Primer trimming (02)
- cutadapt --revcomp used instead of two-pass approach: handles both read orientations in one pass, outputs all reads in forward orientation
- --error-rate 0.15, --overlap 15: tolerant settings for ONT reads with IUPAC degenerate primers
- Size filter 1300-1800 bp applied post-trimming (insert size)

### Dereplication (03)
- All barcodes pooled before dereplication (not per-sample): required for consistent ZOTU IDs across samples
- Read headers tagged with sample=barcode (e.g. sample=barcode08): required by usearch -otutab for per-sample abundance table
- 99.9% of unique sequences per sample are singletons — expected for ONT full-length 16S

### UNOISE3 denoising (04)
- minsize sweep 1-8 revealed a hard cliff at minsize=3 (272 → 14 ZOTUs)
- minsize=1 and minsize=2 give identical results (272): UNOISE3 error model implicitly rejects all singletons regardless of minsize setting
- minsize=3 (14 ZOTUs) is the working threshold: biologically plausible for this low-diversity system
- minsize=8 (from reference paper) is too aggressive for this sequencing depth
- Some ZOTUs represent intragenomic 16S copy variants from the same organism
- Some ZOTUs match reference database at 100% identity — validates pipeline accuracy
- 4 ZOTUs (Zotu1/2/12/14) are Isochrysis galbana chloroplast 16S (see taxonomy results) — excluded from microbiome analysis

### ZOTU table (05)
- Script takes ZOTU FASTA as argument to allow comparison across minsize thresholds
- usearch -otutab ignores quality scores — mapping only, quality filtering already done upstream
- Output named after input ZOTU file (e.g. zotu_table_zotus_minsize3.txt)

## Known limitations and conclusions
- Low sequencing depth (20-60k reads/sample) for full-length 16S ONT: Riisgaard-Jensen et al. 2026 showed ONT requires 4.1-5.6x more reads than PacBio for V1-V8 amplicons to fully resolve communities
- Rare taxa below detection threshold — EMU (closed-reference) identified ~62 species vs 14 ZOTUs at minsize=3
- Intragenomic 16S variants inflate ZOTU count relative to true species count
- UNOISE3 error model designed for Illumina substitution errors, not ONT indel-dominated errors
- 27F/1492R primers amplify host chloroplast 16S in algae cultures (Zotu1/2/12/14 = I. galbana chloroplast). GTDB-only BLAST cannot detect this — always cross-check unknown ZOTUs against NCBI nt. A chloroplast blocking/filtering step should be added to the pipeline.
- ONT full-length 16S is NOT suitable for ASV-level resolution of complex, diverse communities at realistic sequencing depths. The limiting factor is depth, not chemistry accuracy. R10.4 SUP produces reads accurate enough for ASV resolution (100% identity matches confirmed), but the depth required scales prohibitively with community diversity. ONT full-length 16S is appropriate for low-diversity, controlled systems (cultured algae, bioreactors, clinical isolates) where dominant taxa can be recovered at available depths.

## Step 06 — Mixed taxonomy assignment (implemented)

Inspired by NanoASV (https://github.com/ImagoXV/NanoASV), but adapted for de novo ZOTUs rather than read-level reference mapping.

Strategy:
1. BLASTn of ZOTUs against GTDB SSU reference (best hit by bitscore — not -max_target_seqs 1 which does not guarantee best hit per Shah et al. 2019)
2. pident >= 97% + evalue <= 1e-10 → direct taxonomy assignment (known taxa)
3. Below threshold → retain as de novo novel sequences for manual inspection or genome mapping
4. No vsearch clustering of unknowns — dataset is small enough to handle individually

GTDB chosen over SILVA: SILVA species-level taxonomy is not curated; GTDB provides phylogenomically curated taxonomy.
GTDB SSU database: db/gtdb/ (bac120_ssu_reps.fna.gz, ar53_ssu_reps.fna.gz, taxonomy TSVs)

Rationale:
- Avoids closed-reference limitation of NanoASV and EMU for novel/uncharacterised taxa
- BLASTn preferred over Minimap2 for 14 denoised ZOTUs: exhaustive alignment, direct % identity output, no need for heuristic speed
- Genome-mappable sequences are preserved throughout
- Phyloseq output from NanoASV approach worth adopting for downstream R analysis

### Taxonomy results (minsize=3, 14 ZOTUs)
Assigned (pident >= 97%): 10 ZOTUs
- Zotu3:  Sulfitobacter pontiacus (100%)
- Zotu4:  Roseivirga pacifica (100%)
- Zotu5, 7, 8, 11: Phaeobacter piscinae (100%) — intragenomic 16S copy variants
- Zotu6:  Lentilitoribacter sp. (99.6%)
- Zotu9:  Alteromonas marina_A (100%)
- Zotu10: Alteromonas abrolhosensis (99.8%)
- Zotu13: Alteromonas sp. (99.9%)

Host chloroplast (confirmed): 4 ZOTUs (Zotu1, 2, 12, 14)
- NCBI nt BLAST: 100% identity to Isochrysis galbana chloroplast genome (NC_049168.1)
- GTDB BLAST had returned Elusimicrobiota GWA2-66-18 at 88% — a false lead; GTDB contains no
  organellar sequences, so it reported the best prokaryotic hit instead
- 27F/1492R primers amplify plastid 16S (conserved primer sites due to cyanobacterial ancestry)
- I. galbana has a secondary plastid (red-algal origin via secondary endosymbiosis); its 16S
  has diverged from modern cyanobacteria for >1 Ga, explaining why no cyanobacteria hit better
  than 88% in GTDB — at that identity level BLAST hits are not phylogenetically informative
- These ZOTUs are EXCLUDED from all microbiome analyses
- Lesson: GTDB-only BLAST is blind to organellar sequences; any "unknown" ZOTU in an algae
  culture context must be cross-checked against NCBI nt before drawing biological conclusions

## Key reference
Riisgaard-Jensen et al. 2026 — "Nanopore sequencing reaches amplicon sequence variant (ASV) resolution"
bioRxiv: https://www.biorxiv.org/content/10.64898/2026.02.26.708165v1
