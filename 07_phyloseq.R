#!/usr/bin/env Rscript
# Build a phyloseq object from UNOISE3 pipeline outputs (steps 05 + 06).
#
# Replaces the NanoASV branch (steps 08-10) for downstream R analysis.
# Chloroplast ZOTUs are retained — filter via tax_table(ps)[,"Order"] == "Chloroplast"
# in downstream analysis. No sample metadata is included (add later via sample_data<-).
#
# Dependencies (Bioconductor):
#   BiocManager::install(c("phyloseq", "Biostrings"))
#
# Input:
#   pooled/zotus_minsize{N}.fasta
#   results/zotu_table_zotus_minsize{N}.txt
#   results/taxonomy_zotus_minsize{N}/taxonomy_all.tsv
#
# Output:
#   results/phyloseq_minsize{N}.rds
#   results/phyloseq_minsize{N}_otutable.tsv
#   results/phyloseq_minsize{N}_taxonomy.tsv
#
# Usage: Rscript 07_phyloseq.R
#   Run from the project root directory.

suppressPackageStartupMessages({
    library(phyloseq)
    library(Biostrings)
})

# =============================================================================
# Configuration — edit these two lines to adapt to a different run
# =============================================================================

PROJECT_DIR    <- "/work3/josne/Projects/Igalbana_Microbiome_GRF/Igalbana_Microbiome_16S/NanoZotu"
MINSIZE_WORKING <- 5

# =============================================================================
# Paths
# =============================================================================

fasta_path  <- file.path(PROJECT_DIR, "pooled",
                         paste0("zotus_minsize", MINSIZE_WORKING, ".fasta"))
otutab_path <- file.path(PROJECT_DIR, "results",
                         paste0("zotu_table_zotus_minsize", MINSIZE_WORKING, ".txt"))
tax_path    <- file.path(PROJECT_DIR, "results",
                         paste0("taxonomy_zotus_minsize", MINSIZE_WORKING),
                         "taxonomy_all.tsv")
out_rds     <- file.path(PROJECT_DIR, "results",
                         paste0("phyloseq_minsize", MINSIZE_WORKING, ".rds"))
out_otu_tsv <- file.path(PROJECT_DIR, "results",
                         paste0("phyloseq_minsize", MINSIZE_WORKING, "_otutable.tsv"))
out_tax_tsv <- file.path(PROJECT_DIR, "results",
                         paste0("phyloseq_minsize", MINSIZE_WORKING, "_taxonomy.tsv"))

for (f in c(fasta_path, otutab_path, tax_path)) {
    if (!file.exists(f)) stop("Input not found: ", f)
}

# =============================================================================
# 1. OTU table
# =============================================================================

cat("Loading OTU table:", otutab_path, "\n")
otu_raw <- read.table(otutab_path, header = TRUE, sep = "\t",
                      comment.char = "", check.names = FALSE, row.names = 1)
otu_mat <- as.matrix(otu_raw)
storage.mode(otu_mat) <- "integer"
cat("  ZOTUs:", nrow(otu_mat), " Samples:", ncol(otu_mat), "\n")

# =============================================================================
# 2. Taxonomy
# =============================================================================

cat("Loading taxonomy:", tax_path, "\n")
tax_raw <- read.table(tax_path, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "")

ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

tax_split <- do.call(rbind, lapply(tax_raw$taxonomy, function(x) {
    parts <- strsplit(x, ";", fixed = TRUE)[[1]]
    length(parts) <- 7   # pad / truncate to exactly 7 ranks
    parts
}))
colnames(tax_split)  <- ranks
rownames(tax_split)  <- tax_raw$zotu

# Flag chloroplast for easy downstream filtering
is_chloroplast <- tax_split[, "Order"] == "Chloroplast"
is_chloroplast[is.na(is_chloroplast)] <- FALSE
tax_split <- cbind(tax_split,
                   is_chloroplast = ifelse(is_chloroplast, "TRUE", "FALSE"))

# Fill ZOTUs present in OTU table but missing from taxonomy (no BLAST hit)
missing <- setdiff(rownames(otu_mat), rownames(tax_split))
if (length(missing) > 0) {
    cat("  WARNING:", length(missing), "ZOTUs in OTU table have no BLAST hit — filled Unassigned\n")
    fill <- matrix("Unassigned", nrow = length(missing), ncol = ncol(tax_split),
                   dimnames = list(missing, colnames(tax_split)))
    fill[, "is_chloroplast"] <- "FALSE"
    tax_split <- rbind(tax_split, fill)
}

# Subset taxonomy to ZOTUs in OTU table
extra <- setdiff(rownames(tax_split), rownames(otu_mat))
if (length(extra) > 0) {
    cat("  NOTE:", length(extra), "ZOTUs in taxonomy not in OTU table (excluded)\n")
    tax_split <- tax_split[rownames(otu_mat), , drop = FALSE]
}

cat("  Chloroplast ZOTUs:", sum(tax_split[, "Order"] == "Chloroplast", na.rm = TRUE), "\n")

# =============================================================================
# 3. ZOTU sequences
# =============================================================================

cat("Loading FASTA:", fasta_path, "\n")
seqs <- readDNAStringSet(fasta_path)
names(seqs) <- sub(" .*", "", names(seqs))   # strip to >ZotuN

# Align sequence set to OTU table ZOTUs
seqs <- seqs[names(seqs) %in% rownames(otu_mat)]
missing_seqs <- setdiff(rownames(otu_mat), names(seqs))
if (length(missing_seqs) > 0)
    cat("  WARNING:", length(missing_seqs), "ZOTUs in OTU table have no FASTA sequence\n")

cat("  Sequences loaded:", length(seqs), "\n")

# =============================================================================
# 4. Build phyloseq object
# =============================================================================

cat("Building phyloseq object...\n")

ps <- phyloseq(
    otu_table(otu_mat, taxa_are_rows = TRUE),
    tax_table(tax_split[rownames(otu_mat), , drop = FALSE])
)

if (length(seqs) > 0)
    ps <- merge_phyloseq(ps, refseq(seqs))

cat("  ntaxa:", ntaxa(ps), "\n")
cat("  nsamples:", nsamples(ps), "\n")

# =============================================================================
# 5. Save
# =============================================================================

saveRDS(ps, out_rds)
cat("Phyloseq object saved:", out_rds, "\n")

write.table(as.data.frame(otu_table(ps)), out_otu_tsv,
            sep = "\t", quote = FALSE, col.names = NA)
write.table(as.data.frame(tax_table(ps)), out_tax_tsv,
            sep = "\t", quote = FALSE, col.names = NA)
cat("TSV exports written.\n")

cat("\nDone. To load:\n")
cat("  ps <- readRDS(\"", out_rds, "\")\n", sep = "")
cat("  ntaxa(ps); nsamples(ps); head(tax_table(ps))\n")
