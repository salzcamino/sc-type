# ScType Uncertainty Quantification

Comprehensive uncertainty scoring and confidence metrics for cell type annotations.

## Overview

After running ScType, it's critical to assess the confidence of cell type assignments. The uncertainty quantification functions provide:

- **Top N candidates** per cluster (not just the best match)
- **Confidence scores** (0-1 normalized)
- **Confidence levels** (High/Medium/Low)
- **Score differences** between top candidates
- **Visualization** of uncertainty metrics

This enables you to:
- Identify ambiguous annotations that need manual review
- Prioritize high-confidence cell types for downstream analysis
- Detect potentially novel or transitional cell states
- Make informed decisions about annotation quality

## Installation

```r
# No additional packages required beyond ScType dependencies
# But recommended for visualization:
install.packages(c("ggplot2", "dplyr", "patchwork"))
```

## Quick Start

### For Seurat Objects

```r
library(Seurat)

# Load uncertainty functions
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_uncertainty.R")

# Add uncertainty scores to your Seurat object
seurat_obj <- add_sctype_uncertainty(
    seurat_obj,
    known_tissue_type = "Immune system",
    database_file = "ScTypeDB_full.xlsx",
    top_n = 3  # Report top 3 candidates
)

# Visualize uncertainty
plots <- visualize_sctype_uncertainty(
    seurat_obj,
    plot_types = c("candidates", "umap", "distribution", "heatmap"),
    save_plots = TRUE
)

# View results
print(plots$candidates)  # Top candidates per cluster
print(plots$umap$score)  # UMAP colored by confidence
print(plots$distribution$levels)  # Confidence level distribution
print(plots$heatmap)  # Uncertainty heatmap
```

### For SingleCellExperiment Objects

```r
library(SingleCellExperiment)

# Load uncertainty functions
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_uncertainty_sce.R")

# Add uncertainty scores
sce <- add_sctype_uncertainty_sce(
    sce,
    known_tissue_type = "Brain",
    assay_name = "logcounts",
    cluster_col = "cluster",
    top_n = 3
)

# Visualize uncertainty
plots <- visualize_sctype_uncertainty_sce(
    sce,
    plot_types = c("candidates", "umap", "distribution", "heatmap"),
    save_plots = TRUE
)
```

## New Metadata Columns

After running `add_sctype_uncertainty()` or `add_sctype_uncertainty_sce()`, the following columns are added:

| Column | Description | Type | Range |
|--------|-------------|------|-------|
| `sctype_top1` | Primary cell type assignment | Character | - |
| `sctype_score1` | ScType score for top candidate | Numeric | 0 to ∞ |
| `sctype_top2` | Second most likely cell type | Character | - |
| `sctype_score2` | ScType score for 2nd candidate | Numeric | 0 to ∞ |
| `sctype_top3` | Third most likely cell type | Character | - |
| `sctype_score3` | ScType score for 3rd candidate | Numeric | 0 to ∞ |
| `sctype_confidence` | Normalized confidence score | Numeric | 0 to 1 |
| `sctype_confidence_level` | Confidence category | Factor | High/Medium/Low |
| `sctype_score_diff` | Difference: score1 - score2 | Numeric | 0 to ∞ |

### Accessing the Data

```r
# For Seurat
head(seurat_obj@meta.data[, c("sctype_top1", "sctype_confidence", "sctype_confidence_level")])

# View cluster with multiple candidates
seurat_obj@meta.data %>%
    select(seurat_clusters, sctype_top1, sctype_top2, sctype_top3,
           sctype_score1, sctype_score2, sctype_confidence_level) %>%
    distinct()

# For SingleCellExperiment
head(colData(sce)[, c("sctype_top1", "sctype_confidence", "sctype_confidence_level")])
```

## Confidence Metrics Explained

### 1. Raw Scores (`sctype_score1`, `sctype_score2`, `sctype_score3`)

- Direct output from ScType algorithm
- Higher = better match
- Formula: `sum(positive_markers)/sqrt(n_pos) - sum(negative_markers)/sqrt(n_neg)`
- Absolute values depend on cluster size and marker expression

### 2. Score Difference (`sctype_score_diff`)

- Difference between top candidate and second candidate
- Formula: `score1 - score2`
- **Large difference** (>50% of score1) = high confidence, clear winner
- **Small difference** (<20% of score1) = low confidence, ambiguous
- **Interpretation**:
  - If score_diff is large, top annotation is distinctly better
  - If score_diff is small, cluster may be transitional or mixed

### 3. Confidence Score (`sctype_confidence`)

- Normalized 0-1 metric combining absolute and relative confidence
- Formula: `(abs_confidence + rel_confidence) / 2`
  - `abs_confidence = min(score1 / threshold, 1)`
  - `rel_confidence = min(score_diff / (score1 * 0.5), 1)`
- **0.0-0.4** = Low confidence
- **0.4-0.7** = Medium confidence
- **0.7-1.0** = High confidence

### 4. Confidence Level (`sctype_confidence_level`)

- Categorical classification based on confidence score
- **High**: Confidence ≥ 0.7 (trustworthy annotation)
- **Medium**: 0.4 ≤ Confidence < 0.7 (acceptable, but review recommended)
- **Low**: Confidence < 0.4 (manual review required)
- Automatically set to "Low" if score1 < ncells/4 (ScType's Unknown threshold)

## Visualization Types

### 1. Top Candidates Plot

Shows top 1-3 cell type candidates per cluster with their scores.

**Features:**
- Bar plot with clusters on x-axis
- Bars colored by rank (1st=red, 2nd=blue, 3rd=green)
- Cell type labels on bars
- Quickly identify clusters with multiple candidates

**Interpretation:**
- Clusters with similar-height bars for rank 1 and 2 = ambiguous
- Clusters with only one tall bar = confident assignment
- Missing 2nd/3rd bars = no other reasonable candidates

```r
# Generate only this plot
plots <- visualize_sctype_uncertainty(seurat_obj, plot_types = "candidates")
print(plots$candidates)
```

### 2. UMAP Confidence Plots

Two UMAP plots showing confidence spatially.

**Plot A - Confidence Score:**
- Continuous color scale: blue (low) → yellow (mid) → red (high)
- Identifies spatial regions of uncertain annotations
- Helps detect gradients or transitions

**Plot B - Confidence Level:**
- Categorical colors: green (High), yellow (Medium), red (Low)
- Clear visualization of quality zones
- Highlights problematic regions

**Interpretation:**
- Blue/red regions on Plot A = uncertain annotations
- Red regions on Plot B = need manual review
- Gradients may indicate transitional states

```r
plots <- visualize_sctype_uncertainty(seurat_obj, plot_types = "umap")
print(plots$umap$score)
print(plots$umap$level)
```

### 3. Confidence Distribution

Two plots showing overall annotation confidence.

**Plot A - Confidence Levels:**
- Bar chart showing number of cells in each confidence category
- Assess overall annotation quality
- Ideal: most cells in "High" category

**Plot B - Score Differences:**
- Bar chart per cluster showing score_diff
- Colored by confidence level
- Tall bars = confident, short bars = ambiguous

**Interpretation:**
- High proportion of "Low" confidence = challenging dataset
- Clusters with short bars in Plot B = review these first

```r
plots <- visualize_sctype_uncertainty(seurat_obj, plot_types = "distribution")
print(plots$distribution$levels)
print(plots$distribution$score_diff)
```

### 4. Uncertainty Heatmap

Compact visualization showing cell type and confidence per cluster.

**Features:**
- Each tile = one cluster
- Tile color = confidence (red=low, yellow=mid, green=high)
- Cell type name displayed in tile (top text)
- ScType score displayed in tile (bottom text)
- At-a-glance summary of all clusters

**Interpretation:**
- Red tiles = uncertain, need review
- Green tiles = confident, trustworthy
- Quickly scan for problematic clusters

```r
plots <- visualize_sctype_uncertainty(seurat_obj, plot_types = "heatmap")
print(plots$heatmap)
```

## Complete Workflow Example

### Seurat Workflow

```r
library(Seurat)
library(dplyr)
library(ggplot2)

# Standard Seurat preprocessing
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# Add ScType uncertainty scores
source("R/sctype_uncertainty.R")
seurat_obj <- add_sctype_uncertainty(
    seurat_obj,
    known_tissue_type = "Immune system",
    top_n = 3
)

# Generate all visualizations
plots <- visualize_sctype_uncertainty(
    seurat_obj,
    save_plots = TRUE,
    output_dir = "immune_uncertainty"
)

# Examine results
## View top candidates
print(plots$candidates)

## Check confidence on UMAP
print(plots$umap$score)
print(plots$umap$level)

## Overall quality
print(plots$distribution$levels)

## Cluster summary
print(plots$heatmap)

# Identify clusters needing review
low_conf_clusters <- seurat_obj@meta.data %>%
    filter(sctype_confidence_level == "Low") %>%
    pull(seurat_clusters) %>%
    unique()

message("Clusters with low confidence: ", paste(low_conf_clusters, collapse = ", "))

# Examine specific low-confidence cluster
cluster_5_cells <- seurat_obj@meta.data %>%
    filter(seurat_clusters == 5) %>%
    select(sctype_top1, sctype_top2, sctype_top3,
           sctype_score1, sctype_score2, sctype_score3,
           sctype_confidence) %>%
    head(1)

print(cluster_5_cells)

# Decision making based on confidence
## High confidence: use sctype_top1 directly
## Medium confidence: validate with marker visualization
## Low confidence: consider leaving as "Unknown" or manual curation

# Filter to high-confidence annotations only
high_conf_cells <- seurat_obj@meta.data %>%
    filter(sctype_confidence_level == "High")

message(paste0(nrow(high_conf_cells), " cells (",
               round(100*nrow(high_conf_cells)/ncol(seurat_obj), 1),
               "%) have high-confidence annotations"))
```

### SingleCellExperiment Workflow

```r
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)

# Standard SCE preprocessing
sce <- logNormCounts(sce)
sce <- runPCA(sce)
g <- buildSNNGraph(sce, use.dimred = "PCA")
colData(sce)$cluster <- igraph::cluster_louvain(g)$membership
sce <- runUMAP(sce)

# Add uncertainty scores
source("R/sctype_uncertainty_sce.R")
sce <- add_sctype_uncertainty_sce(
    sce,
    known_tissue_type = "Brain",
    assay_name = "logcounts",
    cluster_col = "cluster",
    top_n = 3
)

# Visualize
plots <- visualize_sctype_uncertainty_sce(
    sce,
    save_plots = TRUE,
    output_dir = "brain_uncertainty"
)

# Examine confidence
table(colData(sce)$sctype_confidence_level)

# View cluster-level results
metadata(sce)[["sctype_uncertainty_clusters"]]
```

## Advanced Usage

### Custom Top N

```r
# Get top 5 candidates instead of 3
seurat_obj <- add_sctype_uncertainty(
    seurat_obj,
    known_tissue_type = "Immune system",
    top_n = 5
)

# Columns added: sctype_top4, sctype_score4, sctype_top5, sctype_score5
```

### Custom Annotation Prefix

```r
# Use different prefix for multiple analyses
seurat_obj <- add_sctype_uncertainty(
    seurat_obj,
    known_tissue_type = "Immune system",
    annotation_prefix = "immune"
)

# Columns: immune_top1, immune_score1, immune_confidence, etc.
```

### Integrating with Standard ScType

```r
# Run standard ScType first
source("R/sctype_wrapper.R")
seurat_obj <- run_sctype(seurat_obj, known_tissue_type = "Immune system")

# Then add uncertainty scores
source("R/sctype_uncertainty.R")
seurat_obj <- add_sctype_uncertainty(seurat_obj, known_tissue_type = "Immune system")

# Compare
table(Standard = seurat_obj$sctype_classification,
      WithUncertainty = seurat_obj$sctype_top1)
```

### Filtering by Confidence

```r
# Keep only high-confidence cells
Idents(seurat_obj) <- "sctype_confidence_level"
high_conf_subset <- subset(seurat_obj, idents = "High")

# Downstream analysis on high-confidence cells only
```

### Exporting Results

```r
# Create summary table
summary_table <- seurat_obj@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(
        n_cells = n(),
        primary_type = first(sctype_top1),
        secondary_type = first(sctype_top2),
        primary_score = first(sctype_score1),
        secondary_score = first(sctype_score2),
        confidence = first(sctype_confidence),
        confidence_level = first(sctype_confidence_level)
    )

write.csv(summary_table, "sctype_uncertainty_summary.csv", row.names = FALSE)

# Export full cell-level data
cell_level <- seurat_obj@meta.data %>%
    select(starts_with("sctype_"))

write.csv(cell_level, "sctype_uncertainty_cells.csv")
```

## Interpretation Guide

### Scenario 1: High Confidence, Single Candidate

```
Cluster 0:
- sctype_top1: "CD4+ T cells" (score: 120)
- sctype_top2: "CD8+ T cells" (score: 45)
- sctype_score_diff: 75
- sctype_confidence: 0.85
- sctype_confidence_level: "High"
```

**Interpretation**: Clear, confident annotation. CD4+ T cells is distinctly the best match. Use this annotation directly.

### Scenario 2: Medium Confidence, Two Close Candidates

```
Cluster 3:
- sctype_top1: "Monocytes" (score: 80)
- sctype_top2: "Macrophages" (score: 65)
- sctype_score_diff: 15
- sctype_confidence: 0.55
- sctype_confidence_level: "Medium"
```

**Interpretation**: Ambiguous. Cluster may be transitional monocytes→macrophages, or contain mixed population. **Action**: Visualize marker genes to distinguish.

### Scenario 3: Low Confidence, Multiple Candidates

```
Cluster 7:
- sctype_top1: "Unknown" (score: 12)
- sctype_top2: "Endothelial" (score: 10)
- sctype_top3: "Fibroblasts" (score: 9)
- sctype_score_diff: 2
- sctype_confidence: 0.15
- sctype_confidence_level: "Low"
```

**Interpretation**: Very uncertain. Score below threshold (ncells/4). All candidates similar. **Action**: Manual curation required. Check for novel cell type, doublets, or low-quality cells.

### Scenario 4: Confident Unknown

```
Cluster 9:
- sctype_top1: "Unknown" (score: 8)
- sctype_top2: "Hepatocytes" (score: 3)
- sctype_score_diff: 5
- sctype_confidence: 0.25
- sctype_confidence_level: "Low"
```

**Interpretation**: None of the known cell types match well. Clear rejection rather than ambiguity. **Action**: Potential novel cell type. Perform differential expression analysis to identify marker genes.

## Troubleshooting

### All clusters show "Unknown"

**Cause**: Wrong tissue type or markers don't match species/dataset

**Solution**:
```r
# Try auto-detection
seurat_obj <- add_sctype_uncertainty(seurat_obj, known_tissue_type = NULL)

# Or check available tissues
library(openxlsx)
db <- read.xlsx("ScTypeDB_full.xlsx")
unique(db$tissueType)
```

### Confidence scores all very low

**Possible causes:**
1. Markers in database don't match your data (species, gene naming)
2. Dataset needs better preprocessing (normalization, scaling)
3. Clusters are truly ambiguous/transitional
4. Resolution too high (over-clustering)

**Solutions:**
```r
# Check gene name overlap
gs_list <- gene_sets_prepare("ScTypeDB_full.xlsx", "Immune system")
markers <- unique(unlist(gs_list$gs_positive))
overlap <- mean(markers %in% rownames(seurat_obj))
message(paste0(round(100*overlap, 1), "% of markers found in dataset"))

# Try different resolution
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)  # Lower resolution
```

### UMAP plots not generated

**Error**: `UMAP not found. Skipping UMAP plots.`

**Solution**:
```r
# For Seurat
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# For SCE
sce <- runUMAP(sce, dimred = "PCA")
```

## Best Practices

1. **Always visualize uncertainty** - Don't trust annotations blindly
2. **Review low-confidence clusters** - These need manual curation
3. **Consider top 2-3 candidates** - May reveal biological insights
4. **Use confidence for filtering** - High-confidence subset for sensitive analyses
5. **Document thresholds** - Record what confidence levels you accept
6. **Combine with marker visualization** - Validate medium-confidence annotations
7. **Check for patterns** - Spatial clusters of low confidence may indicate biology
8. **Report confidence in publications** - Transparency about annotation quality

## Citation

If you use ScType uncertainty quantification in your research, please cite:

Ianevski, A., Giri, A.K. & Aittokallio, T. Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data. Nat Commun 13, 1246 (2022). https://doi.org/10.1038/s41467-022-28803-w

---

*Last Updated: 2025-11-15*
