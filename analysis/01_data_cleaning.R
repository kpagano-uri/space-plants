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
message("  Loaded: ", nrow(contrasts_raw), " rows x ", ncol(contrasts_raw), " columns")

message("Loading sample table...")
samples_raw <- read_csv(
  data_path("GLDS-120_rna_seq_SampleTable_GLbulkRNAseq.csv"),
  show_col_types = FALSE
)
message("  Loaded: ", nrow(samples_raw), " rows x ", ncol(samples_raw), " columns")

# ── 2. Identify column types ──────────────────────────────────────────────────
# The DE file is wide: gene annotation columns first, then per-sample count
# columns, then DE statistic columns (Log2fc_, P.value_, Adj.p.value_) for
# each pairwise contrast. We detect each type automatically by prefix.

message("\n--- Identifying column types ---")

all_cols <- colnames(de_raw)

gene_col   <- "TAIR"
annot_cols <- c("TAIR", "SYMBOL", "GENENAME", "REFSEQ", "ENTREZID",
                "STRING_id", "GOSLIM_IDS")

log2fc_cols <- all_cols[str_detect(all_cols, "^Log2fc_")]
pval_cols   <- all_cols[str_detect(all_cols, "^P\\.value_")]
padj_cols   <- all_cols[str_detect(all_cols, "^Adj\\.p\\.value_")]

message("Log2FC columns found: ",      length(log2fc_cols))
message("P-value columns found: ",     length(pval_cols))
message("Adj p-value columns found: ", length(padj_cols))
message("\nExample Log2FC column:     ", log2fc_cols[1])
message("Example P-value column:     ", pval_cols[1])
message("Example Adj p-value column: ", padj_cols[1])

# ── 3. Identify the 6 focal contrasts ────────────────────────────────────────
# The file contains ALL pairwise comparisons among groups. We only want the
# 6 biologically relevant ones for our research question:
#   Spaceflight (FLT) vs Ground Control (GC), within each genotype x light.
#
# These contrasts contain BOTH "Space Flight" AND "Ground Control", and
# the same genotype on both sides.

message("\n--- Identifying focal contrasts (FLT vs GC only) ---")

# Extract the contrast ID from a Log2fc column by stripping the prefix.
# The contrast ID retains its outer parentheses, e.g.:
#   "(Col-0 & Ground Control & Dark Treatment)v(Col-0 & Space Flight & Dark Treatment)"
contrast_id_from_col <- function(col) str_remove(col, "^Log2fc_")

# Step 1: keep only contrasts comparing Space Flight to Ground Control
focal_log2fc_cols <- log2fc_cols[
  str_detect(log2fc_cols, "Space Flight") &
    str_detect(log2fc_cols, "Ground Control")
]

# Step 2: keep only within-genotype contrasts (same genotype on both sides)
is_within_genotype <- function(col_name) {
  contrast <- str_remove(col_name, "^Log2fc_\\(")
  parts    <- str_split(contrast, fixed(")v("))[[1]]
  if (length(parts) != 2) return(FALSE)
  left       <- parts[1]
  right      <- str_remove(parts[2], "\\)$")
  geno_left  <- str_trim(str_extract(left,  "^[^&]+"))
  geno_right <- str_trim(str_extract(right, "^[^&]+"))
  geno_left == geno_right
}

focal_log2fc_cols <- focal_log2fc_cols[sapply(focal_log2fc_cols, is_within_genotype)]

message("Focal contrasts identified: ", length(focal_log2fc_cols))
message("Contrast names:")
for (col in focal_log2fc_cols) message("  ", contrast_id_from_col(col))

# Build matching p-value and adj p-value column names.
# The contrast ID includes outer parens, e.g. "(Col-0 & ...)v(Col-0 & ...)"
# The p-value columns are named: P.value_(Col-0 & ...)v(Col-0 & ...)
# So we strip the outer parens from the contrast ID before concatenating.

focal_contrast_ids       <- contrast_id_from_col(focal_log2fc_cols)
focal_contrast_ids_clean <- str_remove(focal_contrast_ids, "^\\(") %>%
  str_remove("\\)$")

focal_pval_cols <- paste0("P.value_(", focal_contrast_ids_clean, ")")
focal_padj_cols <- paste0("Adj.p.value_(", focal_contrast_ids_clean, ")")

# Verify columns exist in the data
focal_pval_cols <- focal_pval_cols[focal_pval_cols %in% all_cols]
focal_padj_cols <- focal_padj_cols[focal_padj_cols %in% all_cols]

message("\nMatching p-value columns found: ",    length(focal_pval_cols))
message("Matching adj p-value columns found: ", length(focal_padj_cols))

# ── 4. Reshape focal contrasts to long format ─────────────────────────────────

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
    mutate(contrast = str_remove(col_pval, "^P\\.value_\\(") %>% str_remove("\\)$")) %>%
    mutate(contrast = paste0("(", contrast, ")")) %>%
    select(all_of(gene_col), contrast, pvalue)
} else {
  message("[WARN] No matching p-value columns found — pvalue will be NA")
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
    mutate(contrast = str_remove(col_padj, "^Adj\\.p\\.value_\\(") %>% str_remove("\\)$")) %>%
    mutate(contrast = paste0("(", contrast, ")")) %>%
    select(all_of(gene_col), contrast, padj)
} else {
  message("[WARN] No matching adj p-value columns found — padj will be NA")
  de_padj <- NULL
}

de_long <- de_log2fc
if (!is.null(de_pval)) de_long <- left_join(de_long, de_pval, by = c(gene_col, "contrast"))
if (!is.null(de_padj)) de_long <- left_join(de_long, de_padj, by = c(gene_col, "contrast"))

message("Long-format table: ", nrow(de_long), " rows x ", ncol(de_long), " columns")

# ── 5. Parse genotype and light condition ─────────────────────────────────────
# Actual genotype names found in the contrast strings:
#   "Col-0"                 → Col-0 (wild-type Columbia ecotype)
#   "Col-0 PhyD"            → phyD  (phytochrome D loss-of-function mutant)
#   "Wassilewskija ecotype" → WS    (wild-type Wassilewskija ecotype)

message("\n--- Parsing genotype and light condition ---")

de_long <- de_long %>%
  mutate(
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
    median_log2fc   = round(median(log2fc, na.rm = TRUE), 3),
    .groups         = "drop"
  )

message("\nSignificant DE genes (padj < 0.05) per group:")
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
