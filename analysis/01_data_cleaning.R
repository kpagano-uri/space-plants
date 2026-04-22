# =============================================================================
# 01_data_cleaning.R
# Space Plants Project — GLDS-120 (CARA Experiment)
#
# Research Question:
#   Do different Arabidopsis thaliana genotypes (Col-0, WS, and phyD) show
#   different transcriptional responses to spaceflight microgravity, and does
#   light environment (light vs. dark) modulate the magnitude or direction
#   of these differences?
#
# Purpose of this script:
#   Load the raw NASA OSDR files, inspect their structure, and reshape them
#   into a clean "tidy" format ready for analysis. All column names and
#   contrast labels are discovered from the data itself — nothing is hardcoded.
#
# Input:
#   data/GLDS-120_rna_seq_differential_expression_rRNArm_GLbulkRNAseq.csv
#   data/GLDS-120_rna_seq_contrasts_GLbulkRNAseq.csv
#   data/GLDS-120_rna_seq_SampleTable_GLbulkRNAseq.csv
#
# Output:
#   data/cleaned_de_long.csv       — tidy long-format differential expression
#   data/cleaned_sample_table.csv  — cleaned sample metadata
#   data/cleaned_summary.csv       — per-group summary statistics
#
# Usage:
#   Run from the project root in RStudio with the .Rproj open, or:
#   Rscript analysis/01_data_cleaning.R
# =============================================================================

library(tidyverse)

# ── 0. Locate data files relative to project root ────────────────────────────

if (requireNamespace("here", quietly = TRUE)) {
  data_dir <- here::here("data")
} else {
  data_dir <- file.path(getwd(), "data")
}

message("Looking for data in: ", data_dir)

data_path <- function(filename) {
  path <- file.path(data_dir, filename)
  if (!file.exists(path)) {
    stop(
      "File not found: ", path,
      "\n  Run data/download_data.sh first to download the data."
    )
  }
  path
}

# ── 1. Load raw files ─────────────────────────────────────────────────────────

message("\n--- Loading raw data files ---")

message("Loading differential expression table (this may take a moment)...")
de_raw <- read_csv(
  data_path("GLDS-120_rna_seq_differential_expression_rRNArm_GLbulkRNAseq.csv"),
  show_col_types = FALSE
)
message("  Loaded: ", nrow(de_raw), " rows x ", ncol(de_raw), " columns")

message("Loading contrasts table...")
contrasts_raw <- read_csv(
  data_path("GLDS-120_rna_seq_contrasts_GLbulkRNAseq.csv"),
  show_col_types = FALSE
)

message("Loading sample table...")
samples_raw <- read_csv(
  data_path("GLDS-120_rna_seq_SampleTable_GLbulkRNAseq.csv"),
  show_col_types = FALSE
)

# ── 2. Identify column types ──────────────────────────────────────────────────

message("\n--- Identifying column types ---")

all_cols    <- colnames(de_raw)
gene_col    <- "TAIR"
annot_cols  <- c("TAIR", "SYMBOL", "GENENAME", "REFSEQ", "ENTREZID",
                 "STRING_id", "GOSLIM_IDS")

log2fc_cols <- all_cols[str_detect(all_cols, "^Log2fc_")]
pval_cols   <- all_cols[str_detect(all_cols, "^P\\.value_")]
padj_cols   <- all_cols[str_detect(all_cols, "^Adj\\.p\\.value_")]

message("Log2FC columns found: ",      length(log2fc_cols))
message("P-value columns found: ",     length(pval_cols))
message("Adj p-value columns found: ", length(padj_cols))

# ── 3. Identify the 6 focal contrasts ────────────────────────────────────────
# The file contains ALL pairwise comparisons. We need exactly 6:
#   Space Flight vs Ground Control, within each genotype x light condition.
#
# Three filtering steps:
#   Step A — keep only contrasts involving both Space Flight and Ground Control
#   Step B — keep only within-genotype contrasts (same genotype on both sides)
#   Step C — keep only contrasts where Space Flight is on the LEFT (numerator)
#             so that: positive log2fc = up-regulated in spaceflight
#                      negative log2fc = down-regulated in spaceflight

message("\n--- Identifying focal contrasts (FLT vs GC only) ---")

contrast_id_from_col <- function(col) str_remove(col, "^Log2fc_")

# Helper: split contrast into left and right groups around ")v("
split_contrast <- function(col_name) {
  contrast <- str_remove(col_name, "^Log2fc_\\(")
  parts    <- str_split(contrast, fixed(")v("))[[1]]
  if (length(parts) != 2) return(list(left = NA, right = NA))
  list(
    left  = parts[1],
    right = str_remove(parts[2], "\\)$")
  )
}

# Step A: both Space Flight and Ground Control must appear
focal <- log2fc_cols[
  str_detect(log2fc_cols, "Space Flight") &
    str_detect(log2fc_cols, "Ground Control")
]

# Step B: same genotype on both sides
is_within_genotype <- function(col_name) {
  parts      <- split_contrast(col_name)
  geno_left  <- str_trim(str_extract(parts$left,  "^[^&]+"))
  geno_right <- str_trim(str_extract(parts$right, "^[^&]+"))
  !is.na(geno_left) && geno_left == geno_right
}
focal <- focal[sapply(focal, is_within_genotype)]

# Step C: Space Flight must be on the LEFT (numerator) so positive log2fc
# means up-regulated in spaceflight. This also eliminates the duplicate
# mirror contrasts that otherwise produce perfectly symmetric 50/50 splits.
is_spaceflight_left <- function(col_name) {
  parts <- split_contrast(col_name)
  str_detect(parts$left, "Space Flight")
}
focal_log2fc_cols <- focal[sapply(focal, is_spaceflight_left)]

message("Focal contrasts identified: ", length(focal_log2fc_cols))
message("(Should be exactly 6 — one per genotype x light condition)")
for (col in focal_log2fc_cols) message("  ", contrast_id_from_col(col))

# Build matching p-value and adj p-value column names.
# contrast ID has outer parens, e.g. "(FLT group)v(GC group)"
# p-value columns are named P.value_(FLT group)v(GC group) — no double parens
focal_contrast_ids       <- contrast_id_from_col(focal_log2fc_cols)
focal_contrast_ids_clean <- str_remove(focal_contrast_ids, "^\\(") %>%
  str_remove("\\)$")

focal_pval_cols <- paste0("P.value_(", focal_contrast_ids_clean, ")")
focal_padj_cols <- paste0("Adj.p.value_(", focal_contrast_ids_clean, ")")

focal_pval_cols <- focal_pval_cols[focal_pval_cols %in% all_cols]
focal_padj_cols <- focal_padj_cols[focal_padj_cols %in% all_cols]

message("\nMatching p-value columns:     ", length(focal_pval_cols))
message("Matching adj p-value columns: ", length(focal_padj_cols))

# ── 4. Reshape to long (tidy) format ─────────────────────────────────────────

message("\n--- Reshaping to long (tidy) format ---")

de_log2fc <- de_raw %>%
  select(all_of(c(annot_cols, focal_log2fc_cols))) %>%
  pivot_longer(
    cols      = all_of(focal_log2fc_cols),
    names_to  = "col_log2fc",
    values_to = "log2fc"
  ) %>%
  mutate(contrast = contrast_id_from_col(col_log2fc)) %>%
  select(-col_log2fc)

if (length(focal_pval_cols) > 0) {
  de_pval <- de_raw %>%
    select(all_of(c(gene_col, focal_pval_cols))) %>%
    pivot_longer(
      cols      = all_of(focal_pval_cols),
      names_to  = "col_pval",
      values_to = "pvalue"
    ) %>%
    mutate(contrast = str_remove(col_pval, "^P\\.value_\\(") %>%
             str_remove("\\)$") %>%
             { paste0("(", ., ")") }) %>%
    select(all_of(gene_col), contrast, pvalue)
} else {
  de_pval <- NULL
}

if (length(focal_padj_cols) > 0) {
  de_padj <- de_raw %>%
    select(all_of(c(gene_col, focal_padj_cols))) %>%
    pivot_longer(
      cols      = all_of(focal_padj_cols),
      names_to  = "col_padj",
      values_to = "padj"
    ) %>%
    mutate(contrast = str_remove(col_padj, "^Adj\\.p\\.value_\\(") %>%
             str_remove("\\)$") %>%
             { paste0("(", ., ")") }) %>%
    select(all_of(gene_col), contrast, padj)
} else {
  de_padj <- NULL
}

de_long <- de_log2fc
if (!is.null(de_pval)) de_long <- left_join(de_long, de_pval, by = c(gene_col, "contrast"))
if (!is.null(de_padj)) de_long <- left_join(de_long, de_padj, by = c(gene_col, "contrast"))

message("Long-format rows: ", nrow(de_long),
        " (expected ~", 24725 * 6, " = 24725 genes x 6 contrasts)")

# ── 5. Parse genotype and light condition ─────────────────────────────────────
# Genotype names in contrast strings:
#   "Col-0 PhyD"            → phyD
#   "Wassilewskija ecotype" → WS
#   "Col-0"                 → Col-0
#
# Since Space Flight is now always on the LEFT side:
#   positive log2fc = gene is MORE expressed in spaceflight (up-regulated)
#   negative log2fc = gene is LESS expressed in spaceflight (down-regulated)

message("\n--- Parsing genotype and light condition ---")

de_long <- de_long %>%
  mutate(
    # Left side of contrast = Space Flight group; extract genotype from it
    contrast_left = str_extract(contrast, "^[^&]+"),
    
    genotype = case_when(
      str_detect(contrast_left, "Col-0 PhyD")           ~ "phyD",
      str_detect(contrast_left, "Wassilewskija ecotype") ~ "WS",
      str_detect(contrast_left, "Col-0")                ~ "Col-0",
      TRUE                                               ~ "Unknown"
    ),
    
    light_condition = case_when(
      str_detect(contrast, "Light Treatment") ~ "Light",
      str_detect(contrast, "Dark Treatment")  ~ "Dark",
      TRUE                                   ~ "Unknown"
    ),
    
    # positive log2fc = up in spaceflight (Space Flight is numerator)
    direction = case_when(
      is.na(log2fc) ~ NA_character_,
      log2fc > 0    ~ "Up in spaceflight",
      log2fc < 0    ~ "Down in spaceflight",
      TRUE          ~ "No change"
    ),
    
    significant = !is.na(padj) & padj < 0.05
  ) %>%
  select(-contrast_left)

message("\nGenotypes detected:")
print(table(de_long$genotype))
message("\nLight conditions detected:")
print(table(de_long$light_condition))

# Sanity check: each TAIR gene should appear exactly 6 times (once per contrast)
n_per_gene <- de_long %>% count(TAIR) %>% pull(n)
message("\nRows per gene — should all be 6:")
print(table(n_per_gene))

# ── 6. Clean sample table ─────────────────────────────────────────────────────

message("\n--- Cleaning sample table ---")
samples_clean <- samples_raw %>%
  rename_with(~ str_to_lower(str_replace_all(., "\\s+", "_")))
print(samples_clean)

# ── 7. Summary ────────────────────────────────────────────────────────────────

message("\n--- Summary of cleaned data ---")

summary_table <- de_long %>%
  group_by(genotype, light_condition) %>%
  summarise(
    total_genes     = n_distinct(TAIR),
    n_significant   = sum(significant, na.rm = TRUE),
    pct_significant = round(mean(significant, na.rm = TRUE) * 100, 1),
    n_up            = sum(significant & direction == "Up in spaceflight",   na.rm = TRUE),
    n_down          = sum(significant & direction == "Down in spaceflight",  na.rm = TRUE),
    median_log2fc   = round(median(log2fc, na.rm = TRUE), 3),
    .groups         = "drop"
  )

message("\nSummary (significant DE genes, padj < 0.05):")
print(summary_table)

# ── 8. Save outputs ───────────────────────────────────────────────────────────

message("\n--- Saving cleaned files ---")

write_csv(de_long,       file.path(data_dir, "cleaned_de_long.csv"))
write_csv(samples_clean, file.path(data_dir, "cleaned_sample_table.csv"))
write_csv(summary_table, file.path(data_dir, "cleaned_summary.csv"))

message("Saved: data/cleaned_de_long.csv")
message("Saved: data/cleaned_sample_table.csv")
message("Saved: data/cleaned_summary.csv")

message("\n========================================")
message(" Cleaning complete! Ready for analysis.")
message(" Next: run analysis/02_exploratory.R")
message("========================================")
