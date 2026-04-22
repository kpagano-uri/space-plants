# =============================================================================
# 04_figures.R
# Space Plants Project — GLDS-120 (CARA Experiment)
#
# Research Question:
#   Do different Arabidopsis thaliana genotypes (Col-0, WS, and phyD) show
#   different transcriptional responses to spaceflight microgravity, and does
#   light environment (light vs. dark) modulate the magnitude or direction
#   of these differences?
#
# Purpose of this script:
#   Combine the best figures from exploratory and comparison scripts into
#   two clean multi-panel figures suitable for the paper.
#
#   Figure A (figures/FigureA_overview.png):
#     Panel A — DE gene counts (bar chart)
#     Panel B — Log2FC distributions (density)
#     Panel C — Volcano grid
#     Panel D — Up/down direction proportions
#
#   Figure B (figures/FigureB_statistics.png):
#     Panel A — Fold change magnitude boxplots
#     Panel B — Jaccard overlap heatmap
#
# Input:
#   data/cleaned_de_long.csv
#   data/cleaned_summary.csv
#   tables/02_kruskal_results.csv
#
# Output:
#   figures/FigureA_overview.png
#   figures/FigureB_statistics.png
#
# Usage:
#   Run from the project root in RStudio with the .Rproj open, or:
#   Rscript analysis/04_figures.R
# =============================================================================

library(tidyverse)
library(scales)
library(patchwork)

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

dir.create(fig_dir, showWarnings = FALSE)

theme_space <- theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(color = "grey40", size = 9),
    strip.text      = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.text     = element_text(size = 9),
    legend.title    = element_text(size = 9)
  )

genotype_colors <- c("Col-0" = "#2166ac", "WS" = "#d6604d", "phyD" = "#4dac26")

# ── 1. Load data ──────────────────────────────────────────────────────────────

message("Loading data...")

de_long <- read_csv(
  file.path(data_dir, "cleaned_de_long.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    genotype        = factor(genotype,        levels = c("Col-0", "WS", "phyD")),
    light_condition = factor(light_condition, levels = c("Light", "Dark"))
  )

summary_table <- read_csv(
  file.path(data_dir, "cleaned_summary.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    genotype        = factor(genotype,        levels = c("Col-0", "WS", "phyD")),
    light_condition = factor(light_condition, levels = c("Light", "Dark"))
  )

# ── 2. Build panels ───────────────────────────────────────────────────────────

message("Building panels...")

# ── Panel A: DE gene counts ───────────────────────────────────────────────────

pA <- summary_table %>%
  ggplot(aes(x = genotype, y = n_significant,
             fill = genotype, alpha = light_condition)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "white") +
  geom_text(aes(label = comma(n_significant)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3, fontface = "bold") +
  scale_fill_manual(values = genotype_colors, guide = "none") +
  scale_alpha_manual(values = c("Light" = 1, "Dark" = 0.35),
                     labels = c("Light" = "Light", "Dark" = "Dark")) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.18))) +
  labs(title = "A — DE Gene Counts",
       x = "Genotype", y = "Number of DE genes",
       alpha = "Light condition") +
  theme_space +
  theme(axis.text.x = element_text(
    color = genotype_colors[c("Col-0", "WS", "phyD")], face = "bold"))

# ── Panel B: Log2FC density distributions ────────────────────────────────────

de_plot <- de_long %>%
  filter(!is.na(log2fc)) %>%
  filter(log2fc > quantile(log2fc, 0.005) &
           log2fc < quantile(log2fc, 0.995))

pB <- de_plot %>%
  ggplot(aes(x = log2fc, color = genotype, fill = genotype)) +
  geom_density(alpha = 0.2, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ light_condition, ncol = 1,
             labeller = labeller(light_condition = c(
               Light = "Light", Dark = "Dark"))) +
  scale_color_manual(values = genotype_colors) +
  scale_fill_manual(values  = genotype_colors) +
  coord_cartesian(xlim = c(-3, 3)) +
  labs(title = "B — Log2FC Distributions",
       x = "Log2FC (Spaceflight / Ground Control)",
       y = "Density", color = "Genotype", fill = "Genotype") +
  theme_space

# ── Panel C: Volcano grid ─────────────────────────────────────────────────────

de_volcano <- de_long %>%
  filter(!is.na(log2fc), !is.na(padj)) %>%
  mutate(
    neg_log10_padj = pmin(-log10(padj), 50),
    de_class = case_when(
      padj < 0.05 & log2fc >  1 ~ "Up (sig.)",
      padj < 0.05 & log2fc < -1 ~ "Down (sig.)",
      padj < 0.05               ~ "Sig. (small FC)",
      TRUE                      ~ "Not significant"
    ),
    de_class = factor(de_class, levels = c("Up (sig.)", "Down (sig.)",
                                           "Sig. (small FC)", "Not significant"))
  )

de_class_colors <- c("Up (sig.)"       = "#d6604d",
                     "Down (sig.)"     = "#2166ac",
                     "Sig. (small FC)" = "#f4a582",
                     "Not significant" = "grey85")

panel_counts <- de_long %>%
  filter(significant) %>%
  count(genotype, light_condition, name = "n_sig") %>%
  mutate(label = paste0("n=", comma(n_sig)),
         log2fc = -5.5, neg_log10_padj = Inf)

pC <- de_volcano %>%
  ggplot(aes(x = log2fc, y = neg_log10_padj, color = de_class)) +
  geom_point(alpha = 0.35, size = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_text(data = panel_counts,
            aes(x = log2fc, y = neg_log10_padj, label = label),
            inherit.aes = FALSE,
            vjust = 1.8, hjust = 0, size = 2.5,
            color = "grey30", fontface = "italic") +
  facet_grid(light_condition ~ genotype, scales = "free_y") +
  scale_color_manual(values = de_class_colors) +
  coord_cartesian(xlim = c(-6, 6)) +
  labs(title = "C — Volcano Plots",
       x = "Log2 Fold Change", y = "-log10(padj)",
       color = "Gene class") +
  theme_space +
  theme(legend.position = "right")

# ── Panel D: Up/down proportions ──────────────────────────────────────────────

de_direction <- de_long %>%
  filter(significant) %>%
  count(genotype, light_condition, direction) %>%
  group_by(genotype, light_condition) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

pD <- de_direction %>%
  ggplot(aes(x = genotype, y = pct, fill = direction)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = paste0(round(pct), "%")),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white", fontface = "bold") +
  facet_wrap(~ light_condition,
             labeller = labeller(light_condition = c(
               Light = "Light", Dark = "Dark"))) +
  scale_fill_manual(values = c("Up in spaceflight"   = "#d6604d",
                               "Down in spaceflight" = "#2166ac")) +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  labs(title = "D — Up/Down Proportions",
       x = "Genotype", y = "% of significant DE genes",
       fill = "Direction") +
  theme_space

# ── Panel E: Fold change magnitude boxplots ───────────────────────────────────
# KW p-value in strip label — clean placement, no overlap with boxes

de_sig <- de_long %>%
  filter(significant, !is.na(log2fc)) %>%
  mutate(abs_log2fc = abs(log2fc))

pE <- de_sig %>%
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
  labs(title    = "A — Fold Change Magnitude",
       x        = "Genotype",
       y        = "|Log2 Fold Change|",
       caption  = "KW = Kruskal-Wallis test") +
  theme_space +
  theme(
    axis.text.x = element_text(
      color = genotype_colors[c("Col-0", "WS", "phyD")],
      face  = "bold"
    ),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
  )

# ── Panel F: Jaccard overlap heatmap ─────────────────────────────────────────

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
    s_i   <- sig_gene_sets$genes[[i]]
    s_j   <- sig_gene_sets$genes[[j]]
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

pF <- jaccard_long %>%
  ggplot(aes(x = group_1, y = group_2, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3, color = "white", fontface = "bold") +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b",
                      labels = percent_format(),
                      name   = "Jaccard\nsimilarity") +
  labs(title = "B — DE Gene Overlap (Jaccard Similarity)",
       x = NULL, y = NULL) +
  theme_space +
  theme(axis.text.x     = element_text(angle = 30, hjust = 1, size = 9),
        axis.text.y     = element_text(size = 9),
        legend.position = "right")

# ── 3. Stitch Figure A ────────────────────────────────────────────────────────

message("Stitching Figure A...")

figA <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title    = "Transcriptional Responses of Arabidopsis thaliana to Spaceflight",
    subtitle = "CARA Experiment (NASA OSDR OSD-120) | Spaceflight vs. Ground Control",
    caption  = "Data: NASA OSDR OSD-120 | FDR-adjusted p-value < 0.05",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey40", size = 11),
      plot.caption  = element_text(color = "grey50", size = 8)
    )
  )

ggsave(file.path(fig_dir, "FigureA_overview.png"), figA,
       width = 14, height = 12, dpi = 300)
message("Saved: figures/FigureA_overview.png")

# ── 4. Stitch Figure B ────────────────────────────────────────────────────────

message("Stitching Figure B...")

figB <- (pE | pF) +
  plot_annotation(
    title    = "Statistical Comparison of Spaceflight Responses Across Genotypes",
    subtitle = "CARA Experiment (NASA OSDR OSD-120)",
    caption  = "Data: NASA OSDR OSD-120 | Significantly DE genes (padj < 0.05) only",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey40", size = 11),
      plot.caption  = element_text(color = "grey50", size = 8)
    )
  )

ggsave(file.path(fig_dir, "FigureB_statistics.png"), figB,
       width = 14, height = 7, dpi = 300)
message("Saved: figures/FigureB_statistics.png")

message("\n========================================")
message(" Figure assembly complete!")
message(" figures/FigureA_overview.png")
message(" figures/FigureB_statistics.png")
message("========================================")