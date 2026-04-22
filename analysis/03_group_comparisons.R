# =============================================================================
# 03_group_comparisons.R
# Space Plants Project — GLDS-120 (CARA Experiment)
#
# Research Question:
#   Do different Arabidopsis thaliana genotypes (Col-0, WS, and phyD) show
#   different transcriptional responses to spaceflight microgravity, and does
#   light environment (light vs. dark) modulate the magnitude or direction
#   of these differences?
#
# Purpose of this script:
#   Formal statistical comparisons between genotypes and light conditions.
#   We test three related questions:
#     Q1. Do genotypes differ in the NUMBER of DE genes in spaceflight?
#     Q2. Do genotypes differ in the MAGNITUDE of fold changes?
#     Q3. Do genotypes differ in the DIRECTION (up vs down) of DE genes?
#
# Statistical approach:
#   - Counts (Q1): Fisher's exact test (pairwise)
#   - Fold change magnitudes (Q2): Kruskal-Wallis + pairwise Wilcoxon
#   - Directions (Q3): Chi-square test on contingency tables
#
# Input:
#   data/cleaned_de_long.csv
#   data/cleaned_summary.csv
#
# Output:
#   tables/01_fisher_results.csv
#   tables/02_kruskal_results.csv
#   tables/03_dunn_results.csv
#   tables/04_chisq_results.csv
#   figures/05_log2fc_boxplots.png
#   figures/06_overlap_heatmap.png
#
# Usage:
#   Run from the project root in RStudio with the .Rproj open, or:
#   Rscript analysis/03_group_comparisons.R
# =============================================================================

library(tidyverse)
library(scales)

# ── 0. Setup ──────────────────────────────────────────────────────────────────

if (requireNamespace("here", quietly = TRUE)) {
  data_dir  <- here::here("data")
  fig_dir   <- here::here("figures")
  table_dir <- here::here("tables")
} else {
  data_dir  <- file.path(getwd(), "data")
  fig_dir   <- file.path(getwd(), "figures")
  table_dir <- file.path(getwd(), "tables")
}

dir.create(fig_dir,   showWarnings = FALSE)
dir.create(table_dir, showWarnings = FALSE)

theme_space <- theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", size = 14),
    plot.subtitle   = element_text(color = "grey40", size = 11),
    strip.text      = element_text(face = "bold"),
    legend.position = "bottom"
  )

genotype_colors <- c("Col-0" = "#2166ac", "WS" = "#d6604d", "phyD" = "#4dac26")

# ── 1. Load data ──────────────────────────────────────────────────────────────

message("\n--- Loading cleaned data ---")

de_long <- read_csv(
  file.path(data_dir, "cleaned_de_long.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    genotype        = factor(genotype,        levels = c("Col-0", "WS", "phyD")),
    light_condition = factor(light_condition, levels = c("Light", "Dark"))
  )

message("Loaded: ", nrow(de_long), " rows")

# ── 2. Q1: Fisher's exact test — DE gene counts ───────────────────────────────

message("\n--- Q1: Fisher's exact test ---")

run_fisher <- function(data, light_cond) {
  df        <- data %>% filter(light_condition == light_cond)
  genotypes <- levels(df$genotype)
  pairs     <- combn(genotypes, 2, simplify = FALSE)
  
  map_dfr(pairs, function(pair) {
    g1 <- pair[1]; g2 <- pair[2]
    counts <- df %>%
      filter(genotype %in% pair) %>%
      group_by(genotype) %>%
      summarise(n_sig     = sum(significant, na.rm = TRUE),
                n_not_sig = sum(!significant, na.rm = TRUE),
                .groups   = "drop") %>%
      arrange(match(genotype, pair))
    
    mat  <- matrix(c(counts$n_sig[1], counts$n_not_sig[1],
                     counts$n_sig[2], counts$n_not_sig[2]), nrow = 2)
    test <- fisher.test(mat)
    
    tibble(light_condition = light_cond,
           genotype_1      = g1,
           genotype_2      = g2,
           n_sig_g1        = counts$n_sig[1],
           n_sig_g2        = counts$n_sig[2],
           odds_ratio      = round(test$estimate, 3),
           p_value         = signif(test$p.value, 4),
           significant     = test$p.value < 0.05)
  })
}

fisher_results <- bind_rows(
  run_fisher(de_long, "Light"),
  run_fisher(de_long, "Dark")
)

message("\nFisher's exact test results:")
print(fisher_results)
write_csv(fisher_results, file.path(table_dir, "01_fisher_results.csv"))
message("Saved: tables/01_fisher_results.csv")

# ── 3. Q2: Kruskal-Wallis — fold change magnitude ────────────────────────────

message("\n--- Q2: Kruskal-Wallis test ---")

de_sig <- de_long %>%
  filter(significant, !is.na(log2fc)) %>%
  mutate(abs_log2fc = abs(log2fc))

kruskal_results <- de_sig %>%
  group_by(light_condition) %>%
  summarise(
    kruskal_stat = kruskal.test(abs_log2fc ~ genotype)$statistic,
    kruskal_df   = kruskal.test(abs_log2fc ~ genotype)$parameter,
    kruskal_p    = kruskal.test(abs_log2fc ~ genotype)$p.value,
    .groups      = "drop"
  ) %>%
  mutate(kruskal_stat = round(kruskal_stat, 3),
         kruskal_p    = signif(kruskal_p, 4),
         significant  = kruskal_p < 0.05)

message("\nKruskal-Wallis results:")
print(kruskal_results)
write_csv(kruskal_results, file.path(table_dir, "02_kruskal_results.csv"))
message("Saved: tables/02_kruskal_results.csv")

# ── 4. Post-hoc: Pairwise Wilcoxon ───────────────────────────────────────────

message("\n--- Post-hoc: Pairwise Wilcoxon (Bonferroni corrected) ---")

run_pairwise_wilcox <- function(data, light_cond) {
  df        <- data %>% filter(light_condition == light_cond)
  genotypes <- levels(df$genotype)
  pairs     <- combn(genotypes, 2, simplify = FALSE)
  
  map_dfr(pairs, function(pair) {
    g1 <- pair[1]; g2 <- pair[2]
    v1 <- df %>% filter(genotype == g1) %>% pull(abs_log2fc)
    v2 <- df %>% filter(genotype == g2) %>% pull(abs_log2fc)
    test <- wilcox.test(v1, v2)
    tibble(light_condition  = light_cond,
           genotype_1       = g1,
           genotype_2       = g2,
           median_abs_lfc_1 = round(median(v1), 3),
           median_abs_lfc_2 = round(median(v2), 3),
           W_statistic      = test$statistic,
           p_value          = signif(test$p.value, 4))
  })
}

dunn_results <- bind_rows(
  run_pairwise_wilcox(de_sig, "Light"),
  run_pairwise_wilcox(de_sig, "Dark")
) %>%
  group_by(light_condition) %>%
  mutate(p_adj_bonferroni = signif(p.adjust(p_value, method = "bonferroni"), 4),
         significant      = p_adj_bonferroni < 0.05) %>%
  ungroup()

message("\nPairwise Wilcoxon results:")
print(dunn_results)
write_csv(dunn_results, file.path(table_dir, "03_dunn_results.csv"))
message("Saved: tables/03_dunn_results.csv")

# ── 5. Q3: Chi-square — up/down direction ────────────────────────────────────

message("\n--- Q3: Chi-square test ---")

run_chisq <- function(data, light_cond) {
  df <- data %>%
    filter(light_condition == light_cond, significant, !is.na(direction))
  ct <- df %>%
    count(genotype, direction) %>%
    pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
    column_to_rownames("genotype") %>%
    as.matrix()
  test <- chisq.test(ct)
  tibble(light_condition = light_cond,
         chi_sq_stat     = round(test$statistic, 3),
         df              = test$parameter,
         p_value         = signif(test$p.value, 4),
         significant     = test$p.value < 0.05)
}

chisq_results <- bind_rows(
  run_chisq(de_long, "Light"),
  run_chisq(de_long, "Dark")
)

message("\nChi-square results:")
print(chisq_results)
write_csv(chisq_results, file.path(table_dir, "04_chisq_results.csv"))
message("Saved: tables/04_chisq_results.csv")

# ── 6. Figure 5: Boxplots ─────────────────────────────────────────────────────
# KW p-value placed in the facet strip label — clean and unobtrusive.

message("\n--- Figure 5: Fold change magnitude boxplots ---")

p5 <- de_sig %>%
  ggplot(aes(x = genotype, y = abs_log2fc, fill = genotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.5) +
  geom_jitter(aes(color = genotype), width = 0.15, alpha = 0.05, size = 0.4) +
  facet_wrap(~ light_condition,
             labeller = labeller(light_condition = c(
               Light = "Light  (KW p = 2.91e-40)",
               Dark  = "Dark  (KW p = 5.64e-16)"
             ))) +
  coord_cartesian(ylim = c(0, 7)) +
  scale_fill_manual(values  = genotype_colors, guide = "none") +
  scale_color_manual(values = genotype_colors, guide = "none") +
  labs(
    title    = "Magnitude of Gene Expression Changes in Spaceflight",
    subtitle = "Absolute log2 fold change | Significantly DE genes only (padj < 0.05)",
    x        = "Genotype",
    y        = "|Log2 Fold Change|",
    caption  = "KW = Kruskal-Wallis test | Data: NASA OSDR OSD-120 (CARA experiment)"
  ) +
  theme_space +
  theme(
    axis.text.x = element_text(
      color = genotype_colors[c("Col-0", "WS", "phyD")],
      face  = "bold"
    ),
    plot.margin = margin(t = 15, r = 10, b = 10, l = 10)
  )

ggsave(file.path(fig_dir, "05_log2fc_boxplots.png"), p5,
       width = 9, height = 6, dpi = 300)
message("Saved: figures/05_log2fc_boxplots.png")
print(p5)

# ── 7. Figure 6: Jaccard overlap heatmap ─────────────────────────────────────

message("\n--- Figure 6: DE gene overlap heatmap ---")

sig_gene_sets <- de_long %>%
  filter(significant) %>%
  group_by(genotype, light_condition) %>%
  summarise(genes = list(unique(TAIR)), .groups = "drop") %>%
  mutate(group = paste(genotype, light_condition, sep = "\n"))

groups <- sig_gene_sets$group
n      <- length(groups)

jaccard_mat <- matrix(NA, nrow = n, ncol = n,
                      dimnames = list(groups, groups))
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    s_i  <- sig_gene_sets$genes[[i]]
    s_j  <- sig_gene_sets$genes[[j]]
    inter <- length(intersect(s_i, s_j))
    uni   <- length(union(s_i, s_j))
    jaccard_mat[i, j] <- if (uni > 0) inter / uni else 0
  }
}

jaccard_long <- as.data.frame(jaccard_mat) %>%
  rownames_to_column("group_1") %>%
  pivot_longer(-group_1, names_to = "group_2", values_to = "jaccard") %>%
  mutate(label   = round(jaccard, 2),
         group_1 = factor(group_1, levels = groups),
         group_2 = factor(group_2, levels = groups))

p6 <- jaccard_long %>%
  ggplot(aes(x = group_1, y = group_2, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3.5, color = "white", fontface = "bold") +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b",
                      labels = percent_format(),
                      name   = "Jaccard\nsimilarity") +
  labs(
    title    = "Overlap of Significantly DE Genes Between Groups",
    subtitle = "Jaccard similarity = shared genes / total unique genes across both groups",
    x        = NULL,
    y        = NULL,
    caption  = "Data: NASA OSDR OSD-120 (CARA experiment)"
  ) +
  theme_space +
  theme(axis.text.x  = element_text(angle = 30, hjust = 1),
        legend.position = "right")

ggsave(file.path(fig_dir, "06_overlap_heatmap.png"), p6,
       width = 8, height = 7, dpi = 300)
message("Saved: figures/06_overlap_heatmap.png")
print(p6)

# ── 8. Summary ────────────────────────────────────────────────────────────────

message("\n========================================")
message(" Statistical comparison summary")
message("========================================")
message("\nQ1 — Fisher's exact test:")
print(fisher_results %>% select(light_condition, genotype_1, genotype_2,
                                n_sig_g1, n_sig_g2, p_value, significant))
message("\nQ2 — Kruskal-Wallis:")
print(kruskal_results)
message("\nQ3 — Chi-square:")
print(chisq_results)
message("\n========================================")
message(" Comparisons complete!")
message(" Tables saved to: tables/")
message(" Figures saved to: figures/")
message(" Next: run analysis/04_figures.R")
message("========================================")