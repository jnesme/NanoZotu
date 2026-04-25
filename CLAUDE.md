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
- Conda: qiime2-amplicon-2026.1 (must be activated in all scripts)
- Key tools: cutadapt (>=4.0 for --revcomp), usearch v12

## Pipeline overview
01 → 02 → 03 → 04 → 05 → 06 → 07
                              ↓
                         08 → 09 → 10

- 01–06: core UNOISE3 pipeline (read processing → taxonomy)
- 07: Elusimicrobiota phylogenetic placement tree
- 08–10: NanoASV branch (GTDB SSU + ZOTUs → Phyloseq output)
  - 08: convert GTDB SSU to NanoASV format (run once, interactive)
  - 09: augment with project ZOTUs (run once, interactive)
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

### ZOTU table (05)
- Script takes ZOTU FASTA as argument to allow comparison across minsize thresholds
- usearch -otutab ignores quality scores — mapping only, quality filtering already done upstream
- Output named after input ZOTU file (e.g. zotu_table_zotus_minsize3.txt)

## Known limitations and conclusions
- Low sequencing depth (20-60k reads/sample) for full-length 16S ONT: Riisgaard-Jensen et al. 2026 showed ONT requires 4.1-5.6x more reads than PacBio for V1-V8 amplicons to fully resolve communities
- Rare taxa below detection threshold — EMU (closed-reference) identified ~62 species vs 14 ZOTUs at minsize=3
- Intragenomic 16S variants inflate ZOTU count relative to true species count
- UNOISE3 error model designed for Illumina substitution errors, not ONT indel-dominated errors
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

Unknown (pident ~88%): 4 ZOTUs (Zotu1, 2, 12, 14)
- All hit the same reference: GWA2-66-18 sp965283875 (phylum Elusimicrobiota)
- 88% identity = novel genus/family level — not in GTDB at species resolution
- Elusimicrobiota are known endosymbionts of flagellates, commonly associated with algae — biologically coherent
- Likely intragenomic variants of the same novel organism
- Candidate for genome mapping
- Completely absent from rrnDB (only 3 Elusimicrobiota entries total: Elusimicrobium minutum,
  Endomicrobium proavitum, Candidatus Endomicrobiellum trichonymphae — none from order UBA1565)
  → EMU cannot detect this organism regardless of sequencing depth; confirmed closed-reference limitation

### Phylogenetic placement of Elusimicrobiota ZOTUs (script 07)
- Dedicated Elusimicrobiota tree built with all 271 GTDB SSU reps from the phylum + 4 ZOTUs + archaeal outgroup
- ZOTUs form a tight monophyletic clade with GB_GCA_965283875.1 (GWA2-66-18 sp965283875)
- Clade leaf labels: GB_GCA_965283875.1, Zotu1, Zotu2, Zotu12, Zotu14
- Phylogenetic placement confirms BLAST assignment: the 88% identity is real divergence, not artifact
- The organism is a genuinely novel lineage within Elusimicrobiota with only one close relative in GTDB

## Key reference
Riisgaard-Jensen et al. 2026 — "Nanopore sequencing reaches amplicon sequence variant (ASV) resolution"
bioRxiv: https://www.biorxiv.org/content/10.64898/2026.02.26.708165v1
