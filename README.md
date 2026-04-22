# Plants in Space
## Transcriptional Responses of *Arabidopsis thaliana* to Spaceflight

**Author:** Kathryn Pagano  
**Course:** BIO539  
**Date:** April 2026  

---

## Research Question

Do different *Arabidopsis thaliana* genotypes (Col-0, WS, and phyD) show different transcriptional responses to spaceflight microgravity, and does light environment (light vs. dark) modulate the magnitude or direction of these differences?

---

## Background

The CARA (Characterizing Arabidopsis Root Attraction) experiment flew three genotypes of *Arabidopsis thaliana* aboard the International Space Station (ISS) and compared their gene expression to ground controls. The three genotypes are:

- **Col-0** — wild-type Columbia ecotype
- **WS** — wild-type Wassilewskija ecotype  
- **phyD** — phytochrome D loss-of-function mutant in the Col-0 background

Plants were grown under two light conditions (light and dark), producing six experimental groups. This analysis uses pre-processed RNA-seq differential expression data from NASA's Open Science Data Repository (OSDR) to compare how each genotype responds to spaceflight at the transcriptional level.

---

## Data Source

**NASA Open Science Data Repository (OSDR)**  
Study: OSD-120 (GLDS-120) — CARA Experiment  
URL: https://osdr.nasa.gov/bio/repo/data/studies/OSD-120  

Data are downloaded automatically by the download script — no manual downloads required. The ribosomal RNA-removed differential expression file (`rRNArm`) is used for all analyses as it provides cleaner signal.

> Raw data files (>100 MB) are excluded from this repository via `.gitignore`. They are downloaded automatically by running `data/download_data.sh`.

---

## How to Reproduce This Analysis

### 1. Clone the repository
```bash
git clone https://github.com/kpagano-uri/space-plants.git
cd space-plants
```

### 2. Install R dependencies
Open R or RStudio and run:
```r
install.packages(c("tidyverse", "here", "scales", "patchwork"))
```

### 3. Download the data
Requires Python 3 (checks automatically). Run from the project root:
```bash
bash data/download_data.sh
```
This will download ~500 MB of data from NASA OSDR into the `data/` folder.

### 4. Run the analysis scripts in order
Open `space-plants.Rproj` in RStudio, then source each script:

```r
source("analysis/01_data_cleaning.R")    # Load, reshape, and clean data
source("analysis/02_exploratory.R")      # Exploratory figures
source("analysis/03_group_comparisons.R") # Statistical comparisons
source("analysis/04_figures.R")          # Final multi-panel figures
```

Output figures will appear in `figures/` and statistical tables in `tables/`.

---

## Repository Structure

```
space-plants/
│
├── README.md                          # This file
├── space-plants.Rproj                 # RStudio project file
├── .gitignore                         # Excludes large data files and figures
│
├── data/
│   └── download_data.sh               # Downloads all data from NASA OSDR
│
├── analysis/
│   ├── 01_data_cleaning.R             # Load, filter, reshape, and clean data
│   ├── 02_exploratory.R               # Exploratory plots and distributions
│   ├── 03_group_comparisons.R         # Statistical tests across groups
│   └── 04_figures.R                   # Final combined multi-panel figures
│
├── figures/                           # Output figures (auto-generated, gitignored)
│   ├── FigureA_overview.png           # 4-panel overview figure
│   ├── FigureB_statistics.png         # 2-panel statistics figure
│   └── ...                            # Individual panel figures
│
├── tables/                            # Statistical result tables
│   ├── 01_fisher_results.csv          # Fisher's exact test — DE gene counts
│   ├── 02_kruskal_results.csv         # Kruskal-Wallis — fold change magnitude
│   ├── 03_dunn_results.csv            # Pairwise Wilcoxon post-hoc tests
│   └── 04_chisq_results.csv           # Chi-square — up/down direction
│
├── paper/
│   └── paper.md                       # ~1000 word analysis paper
│
└── ai_statement.md                    # Statement on AI tool usage
```

---

## Key Results

| Group | DE Genes (padj < 0.05) | % Significant | Dominant Direction |
|-------|------------------------|---------------|--------------------|
| Col-0 Light | 777 | 3.1% | Down in spaceflight (56%) |
| Col-0 Dark  | 69  | 0.3% | Up in spaceflight (75%) |
| WS Light    | 3,488 | 14.1% | Up in spaceflight (60%) |
| WS Dark     | 56  | 0.2% | Up in spaceflight (84%) |
| phyD Light  | 1,030 | 4.2% | Up in spaceflight (64%) |
| phyD Dark   | 425 | 1.7% | Up in spaceflight (62%) |

All pairwise genotype comparisons were statistically significant (Fisher's exact test, p < 0.05) except Col-0 vs WS in the dark condition. Fold change magnitudes differed significantly across genotypes in both light conditions (Kruskal-Wallis, p < 2.91e-40). Up/down proportions differed significantly across genotypes in both conditions (chi-square, p < 1.14e-19).

---

## Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| R | ≥ 4.2 | All analysis |
| tidyverse | ≥ 2.0 | Data wrangling and plotting |
| here | ≥ 1.0 | Reproducible file paths |
| scales | ≥ 1.3 | Plot formatting |
| patchwork | ≥ 1.2 | Multi-panel figure assembly |
| Python | ≥ 3.8 | Data download script |

---

## Citation

Paul, A.L., Zupanska, A.K., Schultz, E.R., and Ferl, R.J. (2013). Organ-specific remodeling of the *Arabidopsis* transcriptome in response to spaceflight. *BMC Plant Biology*, 13, 112.

Data accessed from NASA OSDR: https://osdr.nasa.gov/bio/repo/data/studies/OSD-120

---

## AI Statement

See `ai_statement.md` for a full description of how AI tools were used in this project.
