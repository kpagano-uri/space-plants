# Statement on AI Tool Usage

## Project: Plants in Space: Transcriptional Responses of *Arabidopsis thaliana* to Spaceflight

**Author:** Kathryn E. Pagano  
**Course:** BIO539/DSP539 
**Date:** April 2026

---

## Overview

Claude (Anthropic; Sonnet 4.6 Adaptive) was used extensively throughout this project. Below is a detailed description of how AI assistance was used at each stage, and what decisions and interpretations were made independently by me.

---

## How AI Was Used

### 1. GitHub Repository Setup
Claude guided me through setting up Git on Windows, creating the GitHub repository, building the folder structure, and troubleshooting issues including working directory errors, large file size limits, and git history rewrites to remove accidentally committed large files.

### 2. Data Download Script (`data/download_data.sh`)
Claude wrote the bash/Python download script that queries the NASA OSDR API dynamically to discover and download the relevant CSV files. I verified that the script ran correctly and downloaded the expected files.

### 3. Data Cleaning Script (`analysis/01_data_cleaning.R`)
Claude wrote the initial version of the cleaning script. I ran it iteratively and reported errors back to Claude, including issues with contrast name parsing, such as genotype labels like "Wassilewskija ecotype" and "Col-0 PhyD", duplicate mirror contrasts producing artificial 50/50 up/down splits, and mismatched parentheses in column name construction. Claude revised the script based on my error reports and console output. I verified each fix by examining the printed output.

### 4. Exploratory Analysis (`analysis/02_exploratory.R`)
Claude wrote this script. I reviewed all four output figures and identified improvements including a redundant dual legend on the bar chart and flat Dark panels on the volcano grid. Claude revised the script to address both issues (colored x-axis labels, free y-axis scales per row). All figure interpretation was done by me.

### 5. Statistical Comparisons (`analysis/03_group_comparisons.R`)
Claude wrote this script with appropriate statistical tests that I had selected myself (Fisher's exact, Kruskal-Wallis, pairwise Wilcoxon with Bonferroni correction, chi-square, Jaccard similarity). I identified a problem with the Kruskal-Wallis p-value annotation floating awkwardly over the boxplot and requested it be moved to the facet strip label. Claude revised accordingly.

### 6. Figure Assembly (`analysis/04_figures.R`)
Claude wrote the patchwork figure assembly script combining individual panels into Figure 1 and Figure 2. I identified a y-axis clipping issue on the boxplot caused by long strip labels. Claude revised said issue.

### 7. README
Claude outlined the README document by suggesting document structure.

### 8. Troubleshooting
Claude helped diagnose and fix numerous technical issues throughout the project including: PowerShell vs bash syntax differences, Python not being installed, NASA API key mismatches, git push failures due to large files, working directory errors in RStudio, and git history rewrites.

---

## What I Did Independently


- Project design and dataset selection
- Ran every script and verified outputs at each step
- Identified bugs and errors by reading console output and reporting them
- Interpreted the biological meaning of results
- Selected appropriate statistical tests
- Managed the GitHub repository and resolved all git issues
- All written content found in the README, paper, and AI statement

---

## Reflection

AI assistance significantly accelerated the technical setup and coding portions of this project. Without it, setting up the full reproducible workflow, including the API-based download script, tidy data reshaping, and multi-panel figure assembly, would have taken considerably longer. However, the iterative debugging process required active engagement on my part. Every error required me to read the output, understand what went wrong, and communicate it clearly. The scientific interpretation of results was entirely my own work.
