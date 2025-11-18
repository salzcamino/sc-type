# ScType Installation Guide

This guide explains how to install and use the ScType package for automated cell type annotation in single-cell RNA-seq data.

---

## Installation Methods

### Method 1: Install from GitHub (Recommended)

The easiest way to install ScType is directly from GitHub using `devtools`:

```r
# Install devtools if you don't have it
install.packages("devtools")

# Install ScType from GitHub
devtools::install_github("IanevskiAleksandr/sc-type")

# Load the package
library(ScType)
```

**Benefits:**
- One-time installation
- All functions automatically available
- No need to source remote files
- Secure (no HTTP code execution)
- Automatic dependency management

---

### Method 2: Install from Local Directory

If you've cloned the repository locally:

```r
# Install devtools if needed
install.packages("devtools")

# Install from local directory
devtools::install("/path/to/sc-type")

# Load the package
library(ScType)
```

**Use this method if:**
- You want to modify the code
- You're developing new features
- You don't have internet access

---

### Method 3: Build and Install Package Archive

For advanced users who want to create a distributable package:

```bash
# In the sc-type directory
R CMD build .
R CMD INSTALL ScType_2.0.0.tar.gz
```

Then in R:
```r
library(ScType)
```

---

## Dependencies

ScType will automatically install all required dependencies when you use Method 1 or 2. If installing manually, you'll need:

### Required Packages
```r
install.packages(c(
  "dplyr",
  "HGNChelper",
  "openxlsx",
  "ggplot2",
  "patchwork"
))
```

### For Seurat Workflows
```r
install.packages("Seurat")
```

### For SingleCellExperiment Workflows
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "SingleCellExperiment",
  "SummarizedExperiment",
  "scater"
))
```

### Optional Packages
For advanced features like pathway enrichment and visualization:

```r
# For pathway enrichment integration
install.packages("enrichR")
BiocManager::install(c("fgsea", "msigdbr", "clusterProfiler", "org.Hs.eg.db"))

# For heatmap visualizations
BiocManager::install(c("ComplexHeatmap", "circlize"))
```

---

## Quick Start

### Basic Usage (Seurat)

```r
# Load package and dependencies
library(ScType)
library(Seurat)

# Assuming you have a clustered Seurat object
# Run ScType annotation
seurat_obj <- run_sctype(
  seurat_object = seurat_obj,
  known_tissue_type = "Immune system",  # or NULL for auto-detection
  plot = TRUE
)

# View results
table(seurat_obj@meta.data$sctype_classification)
DimPlot(seurat_obj, group.by = "sctype_classification", label = TRUE)
```

### Basic Usage (SingleCellExperiment)

```r
# Load package and dependencies
library(ScType)
library(SingleCellExperiment)
library(scater)

# Run ScType annotation
sce <- run_sctype_sce(
  sce_object = sce,
  known_tissue_type = "Immune system",
  cluster_col = "cluster",
  plot = TRUE
)

# View results
table(colData(sce)$sctype_classification)
plotUMAP(sce, colour_by = "sctype_classification")
```

### Hierarchical Annotation (Broad + Fine)

```r
# Seurat
seurat_obj <- run_sctype_hierarchical(
  seurat_object = seurat_obj,
  known_tissue_type = "Immune system",
  plot = TRUE
)

# View both annotation levels
table(seurat_obj@meta.data$sctype_broad)   # Broad categories
table(seurat_obj@meta.data$sctype_fine)    # Fine subtypes

# SingleCellExperiment
sce <- run_sctype_hierarchical_sce(
  sce_object = sce,
  known_tissue_type = "Immune system",
  plot = TRUE
)
```

---

## Advanced Features

### 1. Marker Visualization

```r
# Load visualization functions
library(ScType)

# Generate all marker visualizations
plots <- visualize_sctype_markers(
  seurat_object = seurat_obj,
  annotation_col = "sctype_classification",
  tissue_type = "Immune system",
  top_n = 5,
  save_plots = TRUE,
  output_dir = "marker_plots"
)

# View individual plots
print(plots$dotplot)
print(plots$heatmap)
```

### 2. Uncertainty Scoring

```r
# Add confidence metrics
seurat_obj <- add_sctype_uncertainty(
  seurat_obj,
  known_tissue_type = "Immune system",
  top_n = 3  # Get top 3 candidates per cluster
)

# Visualize uncertainty
uncertainty_plots <- visualize_sctype_uncertainty(
  seurat_obj,
  save_plots = TRUE
)

# View results
print(uncertainty_plots$candidates)      # Top candidates per cluster
print(uncertainty_plots$umap$score)      # UMAP by confidence
print(uncertainty_plots$distribution$levels)  # Confidence distribution

# Identify low-confidence clusters for manual review
low_conf_cells <- seurat_obj@meta.data %>%
  filter(sctype_confidence_level == "Low")
```

### 3. Pathway Enrichment Integration

```r
# Add pathway-weighted annotations
# Requires: enrichR, fgsea, and/or clusterProfiler
seurat_obj <- add_pathway_weighted_scores(
  seurat_obj,
  known_tissue_type = "Immune system",
  enrichment_tools = c("enrichr", "fgsea", "go"),
  top_n_genes = 200
)

# Visualize pathway support
pathway_plots <- visualize_pathway_support(
  seurat_obj,
  save_plots = TRUE
)

# Compare ScType vs pathway scores
print(pathway_plots$sctype_vs_pathway)
print(pathway_plots$combined_umap)
```

---

## Custom Marker Databases

### Using Custom Markers

If you have your own cell type marker database:

```r
seurat_obj <- run_sctype(
  seurat_object = seurat_obj,
  known_tissue_type = "Immune system",
  custom_marker_file = "/path/to/your/custom_markers.xlsx"
)
```

### Database Format

Your Excel file should have these columns:

| Column | Description | Example |
|--------|-------------|---------|
| `tissueType` | Tissue/organ | "Immune system" |
| `cellName` | Cell type name | "CD8+ T cells" |
| `geneSymbolmore1` | Positive markers (comma-separated) | "CD8A,CD8B,CD3D" |
| `geneSymbolmore2` | Negative markers (comma-separated) | "CD4,CD19" |

### Available Databases

ScType includes several pre-built databases:

```r
# Standard database (most comprehensive)
db_full <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"

# Short database (essential markers only)
db_short <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx"

# Enhanced database (122 cell types with curated markers)
db_enhanced <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_enhanced.xlsx"

# Hierarchical database (broad + fine annotations)
db_hierarchical <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_hierarchical.xlsx"
```

---

## Supported Tissue Types

The default database includes markers for:

- Immune system
- Brain
- Liver
- Pancreas
- Kidney
- Eye
- Lung
- Adrenal
- Heart
- Intestine
- Muscle
- Placenta
- Spleen
- Stomach
- Thymus

If your tissue type is not listed, you can either:
1. Use `known_tissue_type = NULL` for auto-detection
2. Create a custom marker database

---

## Integration with Existing Workflows

### Standard Seurat Pipeline

```r
library(Seurat)
library(ScType)

# Standard Seurat preprocessing
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Add ScType annotation HERE
pbmc <- run_sctype(pbmc, known_tissue_type = "Immune system", plot = TRUE)

# Continue with downstream analysis
DimPlot(pbmc, group.by = "sctype_classification")
```

### Standard SingleCellExperiment Pipeline

```r
library(SingleCellExperiment)
library(scater)
library(scran)
library(ScType)

# Standard SCE preprocessing
sce <- logNormCounts(sce)
sce <- runPCA(sce, ncomponents = 50)

# Clustering
g <- buildSNNGraph(sce, use.dimred = "PCA")
clusters <- igraph::cluster_louvain(g)$membership
colData(sce)$cluster <- clusters

# Dimensionality reduction
sce <- runUMAP(sce, dimred = "PCA")

# Add ScType annotation HERE
sce <- run_sctype_sce(sce, known_tissue_type = "Immune system",
                       cluster_col = "cluster", plot = TRUE)

# Continue with downstream analysis
plotUMAP(sce, colour_by = "sctype_classification")
```

---

## Troubleshooting

### Issue: "Package 'ScType' not found"

**Solution:**
```r
# Make sure devtools is installed
install.packages("devtools")

# Try installing again
devtools::install_github("IanevskiAleksandr/sc-type")
```

### Issue: "Function not found" errors

**Solution:** Make sure you've loaded the package:
```r
library(ScType)
```

### Issue: All clusters assigned "Unknown"

**Possible causes and solutions:**

1. **Wrong tissue type**
   ```r
   # Use auto-detection
   result <- auto_detect_tissue_type(
     path_to_db_file = "ScTypeDB_full.xlsx",
     seuratObject = seurat_obj
   )
   print(result)
   ```

2. **Unscaled data with scaled=TRUE**
   ```r
   # Make sure data is scaled
   seurat_obj <- ScaleData(seurat_obj)
   # Then run ScType
   seurat_obj <- run_sctype(seurat_obj, ...)
   ```

3. **Species mismatch** (human markers on mouse data)
   - Use species-appropriate marker database
   - Create custom markers for your species

### Issue: Dependency installation failures

**Solution for Bioconductor packages:**
```r
# Update Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")

# Then install required packages
BiocManager::install(c("SingleCellExperiment", "scater"))
```

### Issue: Plot generation fails

**Solution:**
```r
# Make sure UMAP/tSNE coordinates exist
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
# Or for SCE
sce <- runUMAP(sce, dimred = "PCA")

# Then run ScType with plot = TRUE
```

---

## Migration from Old Source-Based Loading

### Old Method (NOT RECOMMENDED - Security Risk)

```r
# ❌ OLD WAY - Don't use this anymore
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
```

### New Method (RECOMMENDED - Secure)

```r
# ✅ NEW WAY - Install once, use forever
devtools::install_github("IanevskiAleksandr/sc-type")
library(ScType)

# All functions are now available automatically!
# - gene_sets_prepare()
# - sctype_score()
# - auto_detect_tissue_type()
# - run_sctype()
# - run_sctype_hierarchical()
# - visualize_sctype_markers()
# - add_sctype_uncertainty()
# - add_pathway_weighted_scores()
# And more!
```

**Why the change?**
- ✅ **Security**: No remote code execution vulnerability
- ✅ **Speed**: Functions loaded once, not every time
- ✅ **Reliability**: Works offline after installation
- ✅ **Version control**: Know exactly which version you're using
- ✅ **Dependency management**: R handles all dependencies automatically

---

## Getting Help

- **Documentation**: See README.md and function help pages (`?run_sctype`)
- **Issues**: https://github.com/IanevskiAleksandr/sc-type/issues
- **Publication**: https://doi.org/10.1038/s41467-022-28803-w
- **Web Portal**: http://sctype.app
- **Contact**: aleksandr.ianevski@helsinki.fi

---

## Citation

If you use ScType in your research, please cite:

```
Ianevski, A., Giri, A.K. & Aittokallio, T.
Fully-automated and ultra-fast cell-type identification using specific marker
combinations from single-cell transcriptomic data.
Nat Commun 13, 1246 (2022).
https://doi.org/10.1038/s41467-022-28803-w
```

---

**Version**: 2.0.0
**Last Updated**: 2025-11-18
**License**: GNU General Public License v3.0
