# =============================================================================
# 02_exploratory.R
# Space Plants Project — GLDS-120 (CARA Experiment)
#
# Research Question:
#   Do different Arabidopsis thaliana genotypes (Col-0, WS, and phyD) show
#   different transcriptional responses to spaceflight microgravity, and does
#   light environment (light vs. dark) modulate the magnitude or direction
#   of these differences?
#
# Purpose of this script:
#   Exploratory data analysis (EDA) of the cleaned differential expression
#   data. The goal is to understand distributions, spot patterns, and identify
#   any data quality issues BEFORE running formal statistical tests.
#   All figures are saved to the figures/ directory.
#
# Input:
#   data/cleaned_de_long.csv
#   data/cleaned_summary.csv
#
# Output (saved to figures/):
#   01_de_gene_counts.png       — bar chart of significant DE genes per group
#   02_log2fc_distributions.png — density plots of log2FC per group
#   03_volcano_grid.png         — volcano plots for all 6 contrasts
#   04_direction_proportions.png — up vs down proportions per group
#
# Usage:
#   Run from the project root in RStudio with the .Rproj open, or:
#   Rscript analysis/02_exploratory.R
# =============================================================================

library(tidyverse)
library(scales)

# ── 0. Setup: paths and plot theme ───────────────────────────────────────────

if (requireNamespace("here", quietly = TRUE)) {
  data_dir <- here::here("data")
  fig_dir  <- here::here("figures")
} else {
  data_dir <- file.path(getwd(), "data")
  fig_dir  <- file.path(getwd(), "figures")
}

dir.create(fig_dir, showWarnings = FALSE)

# Consistent ggplot2 theme for all figures
theme_space <- theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 11),
    strip.text    = element_text(face = "bold"),
    legend.position = "bottom"
  )

# Consistent color palette for the 3 genotypes
genotype_colors <- c(
  "Col-0" = "dodgerblue",
  "WS"    = "firebrick1",
  "phyD"  = "seagreen2"
)

message("Figures will be saved to: ", fig_dir)

# ── 1. Load cleaned data ──────────────────────────────────────────────────────

message("\n--- Loading cleaned data ---")

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

message("Loaded: ", nrow(de_long), " rows")

# ── 2. Figure 1: Bar chart of significant DE genes per group ──────────────────
# Improvement: single clean legend using linetype for light condition
# instead of the redundant dual-legend from the previous version.

message("\n--- Figure 1: DE gene counts ---")

p1 <- summary_table %>%
  ggplot(aes(x = genotype, y = n_significant,
             fill = genotype, alpha = light_condition)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "white") +
  geom_text(
    aes(label = comma(n_significant)),
    position = position_dodge(width = 0.7),
    vjust = -0.4, size = 3.5, fontface = "bold"
  ) +
  scale_fill_manual(values = genotype_colors, guide = "none") +
  scale_alpha_manual(
    values = c("Light" = 1, "Dark" = 0.35),
    labels = c("Light" = "Light", "Dark" = "Dark")
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.15))
  ) +
  # Color the genotype labels on the x-axis to match the bars
  scale_x_discrete() +
  labs(
    title   = "Number of Significantly Differentially Expressed Genes in Spaceflight",
    subtitle = "Spaceflight vs. Ground Control | FDR-adjusted p-value < 0.05",
    x       = "Genotype",
    y       = "Number of DE genes",
    alpha   = "Light condition",
    caption = "Data: NASA OSDR OSD-120 (CARA experiment)"
  ) +
  # Color each x-axis label to match its genotype
  theme_space +
  theme(
    axis.text.x = element_text(
      color = genotype_colors[levels(summary_table$genotype)],
      face  = "bold"
    )
  )

ggsave(file.path(fig_dir, "01_de_gene_counts.png"), p1,
       width = 8, height = 5, dpi = 300)
message("Saved: figures/01_de_gene_counts.png")
print(p1)

# ── 3. Figure 2: Log2FC distribution per group ────────────────────────────────
# No changes needed — looks good. Just re-saving cleanly.

message("\n--- Figure 2: Log2FC distributions ---")

de_plot <- de_long %>%
  filter(!is.na(log2fc)) %>%
  filter(log2fc > quantile(log2fc, 0.005) &
           log2fc < quantile(log2fc, 0.995))

p2 <- de_plot %>%
  ggplot(aes(x = log2fc, color = genotype, fill = genotype)) +
  geom_density(alpha = 0.2, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ light_condition, ncol = 1,
             labeller = labeller(light_condition = c(
               Light = "Light Treatment",
               Dark  = "Dark Treatment"
             ))) +
  scale_color_manual(values = genotype_colors) +
  scale_fill_manual(values  = genotype_colors) +
  coord_cartesian(xlim = c(-3, 3)) +
  labs(
    title    = "Distribution of Log2 Fold Changes in Spaceflight vs. Ground Control",
    subtitle = "All genes shown | Dashed line = no change",
    x        = "Log2 Fold Change (Spaceflight / Ground Control)",
    y        = "Density",
    color    = "Genotype",
    fill     = "Genotype",
    caption  = "Data: NASA OSDR OSD-120 (CARA experiment)"
  ) +
  theme_space

ggsave(file.path(fig_dir, "02_log2fc_distributions.png"), p2,
       width = 8, height = 7, dpi = 300)
message("Saved: figures/02_log2fc_distributions.png")
print(p2)

# ── 4. Figure 3: Volcano plots ────────────────────────────────────────────────
# Improvement: free y-axis scales per row so the Dark panels are not
# squashed flat by the much larger y-range in the Light panels.
# Each row now uses its own y-axis scale, making the Dark results readable.

message("\n--- Figure 3: Volcano plots ---")

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

de_class_colors <- c(
  "Up (sig.)"       = "#d6604d",
  "Down (sig.)"     = "#2166ac",
  "Sig. (small FC)" = "#f4a582",
  "Not significant" = "grey80"
)

# Annotation: add DE gene counts per panel as a label in the top corner
panel_counts <- de_long %>%
  filter(significant) %>%
  count(genotype, light_condition, name = "n_sig") %>%
  mutate(
    label = paste0("n = ", comma(n_sig)),
    # Position label in top-left of each panel
    log2fc         = -5.5,
    neg_log10_padj = Inf
  )

p3 <- de_volcano %>%
  ggplot(aes(x = log2fc, y = neg_log10_padj, color = de_class)) +
  geom_point(alpha = 0.4, size = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1),    linetype = "dashed", color = "grey40") +
  # Count label per panel
  geom_text(
    data    = panel_counts,
    aes(x = log2fc, y = neg_log10_padj, label = label),
    inherit.aes = FALSE,
    vjust = 1.5, hjust = 0, size = 3, color = "grey30", fontface = "italic"
  ) +
  # Free y scales: Light and Dark rows get independent y-axes so Dark is readable
  facet_grid(light_condition ~ genotype, scales = "free_y") +
  scale_color_manual(values = de_class_colors) +
  coord_cartesian(xlim = c(-6, 6)) +
  labs(
    title    = "Volcano Plots: Spaceflight vs. Ground Control",
    subtitle = "Each panel = one genotype x light condition | Dashed lines: padj = 0.05, |log2FC| = 1",
    x        = "Log2 Fold Change (Spaceflight / Ground Control)",
    y        = "-log10(Adjusted p-value)",
    color    = "Gene class",
    caption  = "Data: NASA OSDR OSD-120 (CARA experiment)\nNote: y-axis scales differ between Light and Dark rows"
  ) +
  theme_space +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "03_volcano_grid.png"), p3,
       width = 11, height = 7, dpi = 300)
message("Saved: figures/03_volcano_grid.png")
print(p3)

# ── 5. Figure 4: Up vs down proportions ──────────────────────────────────────
# No changes needed — looks good. Just re-saving cleanly.

message("\n--- Figure 4: Up vs Down proportions ---")

de_direction <- de_long %>%
  filter(significant) %>%
  count(genotype, light_condition, direction) %>%
  group_by(genotype, light_condition) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p4 <- de_direction %>%
  ggplot(aes(x = genotype, y = pct, fill = direction)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(
    aes(label = paste0(round(pct), "%")),
    position = position_stack(vjust = 0.5),
    size = 3.5, color = "white", fontface = "bold"
  ) +
  facet_wrap(~ light_condition,
             labeller = labeller(light_condition = c(
               Light = "Light Treatment",
               Dark  = "Dark Treatment"
             ))) +
  scale_fill_manual(
    values = c("Up in spaceflight"   = "firebrick1",
               "Down in spaceflight" = "dodgerblue")
  ) +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  labs(
    title    = "Direction of Significant Gene Expression Changes in Spaceflight",
    subtitle = "Among significantly DE genes (padj < 0.05) only",
    x        = "Genotype",
    y        = "Percentage of significant DE genes",
    fill     = "Direction",
    caption  = "Data: NASA OSDR OSD-120 (CARA experiment)"
  ) +
  theme_space

ggsave(file.path(fig_dir, "04_direction_proportions.png"), p4,
       width = 9, height = 5, dpi = 300)
message("Saved: figures/04_direction_proportions.png")
print(p4)

# ── 6. Top genes summary ──────────────────────────────────────────────────────

message("\n--- Top 10 most significant genes across all contrasts ---")

top_genes <- de_long %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  select(TAIR, SYMBOL, GENENAME, genotype, light_condition, log2fc, padj) %>%
  distinct() %>%
  head(10)

print(top_genes)

message("\n========================================")
message(" Exploratory analysis complete!")
message(" Figures saved to: figures/")
message(" Next: run analysis/03_group_comparisons.R")
message("========================================")
