# ScType Marker Gene Visualization

Comprehensive visualization functions to explore the marker genes used for cell type annotation.

## Overview

After running ScType annotation, you can visualize the expression of the marker genes (both positive and negative) that were used to determine each cell type assignment. This helps you:

- **Validate annotations** by confirming marker expression patterns
- **Quality control** to identify potential misannotations
- **Explore biology** by examining marker co-expression
- **Publication-ready figures** with multiple visualization types

## Available Visualization Types

1. **Violin Plots** - Expression distribution across cell types
2. **UMAP Plots** - Spatial expression patterns on UMAP
3. **Dotplots** - Combined view of all markers across all cell types
4. **Heatmaps** - Hierarchically clustered expression patterns

## Installation

### Required Packages

```r
# CRAN packages
install.packages(c("ggplot2", "dplyr", "openxlsx", "patchwork"))

# Bioconductor packages (for heatmaps)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("ComplexHeatmap", "circlize"))

# For Seurat
install.packages("Seurat")

# For SingleCellExperiment
BiocManager::install(c("SingleCellExperiment", "scater"))
```

## Quick Start

### For Seurat Objects

```r
# Load visualization functions
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_visualize.R")

# Generate all visualization types
plots <- visualize_sctype_markers(
    seurat_object = seurat_obj,
    annotation_col = "sctype_classification",
    database_file = "ScTypeDB_full.xlsx",
    tissue_type = "Immune system",
    top_n = 5,
    plot_types = c("violin", "umap", "dotplot", "heatmap"),
    save_plots = TRUE,
    output_dir = "sctype_plots"
)

# Access individual plots
print(plots$dotplot)
print(plots$heatmap)
plots$violin[["CD4+ T cells"]]  # Violin plots for specific cell type
plots$umap[["B cells"]]         # UMAP plots for specific cell type
```

### For SingleCellExperiment Objects

```r
# Load visualization functions
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_visualize_sce.R")

# Generate all visualization types
plots <- visualize_sctype_markers_sce(
    sce_object = sce,
    annotation_col = "sctype_classification",
    database_file = "ScTypeDB_full.xlsx",
    tissue_type = "Immune system",
    assay_name = "logcounts",
    top_n = 5,
    plot_types = c("violin", "umap", "dotplot", "heatmap"),
    save_plots = TRUE,
    output_dir = "sctype_plots"
)
```

## Detailed Usage

### Function: `visualize_sctype_markers()`

**Parameters:**

- `seurat_object` / `sce_object`: Annotated object with ScType results
- `annotation_col`: Column name with cell type annotations (default: "sctype_classification")
- `database_file`: Path to marker database used for annotation (default: ScTypeDB_full.xlsx from GitHub)
- `tissue_type`: Tissue type used for annotation (e.g., "Immune system", "Brain")
- `assay` / `assay_name`: Assay to use for expression data (default: "RNA" for Seurat, "logcounts" for SCE)
- `top_n`: Number of top markers to show per cell type (default: 5)
- `plot_types`: Vector of plot types to generate (default: all four types)
- `save_plots`: Save plots to files (default: FALSE)
- `output_dir`: Directory for saved plots (default: "sctype_plots")

**Returns:** List of ggplot objects organized by plot type

### 1. Violin Plots

Shows expression distribution of each marker gene across all cell types.

```r
# Generate only violin plots
plots <- visualize_sctype_markers(
    seurat_obj,
    tissue_type = "Immune system",
    plot_types = "violin"
)

# View violin plots for a specific cell type
print(plots$violin[["CD4+ T cells"]])

# The plots show:
# - Positive markers labeled as "(Positive)"
# - Negative markers labeled as "(Negative)"
# - Expression across all annotated cell types
```

**Interpretation:**
- **Positive markers** should be highly expressed in the target cell type
- **Negative markers** should be lowly expressed in the target cell type
- Clear separation indicates good marker quality

### 2. UMAP Plots

Spatial visualization of marker expression on UMAP coordinates.

```r
# Generate only UMAP plots
plots <- visualize_sctype_markers(
    seurat_obj,
    tissue_type = "Immune system",
    plot_types = "umap"
)

# View UMAP plots for a specific cell type
print(plots$umap[["B cells"]])
```

**Interpretation:**
- **Positive markers** should show high expression (red) in clusters of that cell type
- **Negative markers** should show low expression (blue) in those clusters
- Spatial patterns reveal marker specificity

### 3. Dotplot

Combined view of all markers across all cell types in a single plot.

```r
# Generate dotplot
plots <- visualize_sctype_markers(
    seurat_obj,
    tissue_type = "Immune system",
    plot_types = "dotplot"
)

print(plots$dotplot)
```

**Dotplot features:**
- **Dot size**: Percentage of cells expressing the gene
- **Dot color**: Average expression level
- **Rows**: Cell types
- **Columns**: Marker genes

**Interpretation:**
- Large, red dots indicate high expression in many cells
- Look for strong expression of positive markers in corresponding cell types
- Confirm low expression of negative markers

### 4. Heatmap

Hierarchically clustered heatmap of average marker expression.

```r
# Generate heatmap
plots <- visualize_sctype_markers(
    seurat_obj,
    tissue_type = "Immune system",
    plot_types = "heatmap"
)

print(plots$heatmap)
```

**Heatmap features:**
- **Rows**: Marker genes (hierarchically clustered)
- **Columns**: Cell types (hierarchically clustered)
- **Colors**: Scaled expression (blue = low, white = medium, red = high)
- Clustering reveals co-expression patterns

**Interpretation:**
- Cell types with similar expression patterns cluster together
- Marker genes specific to a cell type show clear red blocks
- Helps identify potential misannotations (unexpected patterns)

## Quick Visualization Functions

For rapid exploration of a single cell type:

### Seurat

```r
# Quick visualization for one cell type
source("R/sctype_visualize.R")

plots <- quick_marker_viz(
    seurat_obj,
    cell_type = "CD4+ T cells",
    tissue_type = "Immune system",
    plot_type = "both"  # "violin", "umap", or "both"
)

print(plots$violin)
print(plots$umap)
```

### SingleCellExperiment

```r
# Quick visualization for one cell type
source("R/sctype_visualize_sce.R")

plots <- quick_marker_viz_sce(
    sce,
    cell_type = "CD4+ T cells",
    tissue_type = "Immune system",
    plot_type = "both"
)
```

## Complete Example

### Seurat Workflow

```r
library(Seurat)
library(ggplot2)

# 1. Run ScType annotation
source("R/sctype_wrapper.R")
seurat_obj <- run_sctype(
    seurat_obj,
    known_tissue_type = "Immune system",
    custom_marker_file = "ScTypeDB_full.xlsx",
    plot = TRUE
)

# 2. Visualize all markers
source("R/sctype_visualize.R")
all_plots <- visualize_sctype_markers(
    seurat_obj,
    annotation_col = "sctype_classification",
    tissue_type = "Immune system",
    top_n = 5,
    save_plots = TRUE,
    output_dir = "immune_markers"
)

# 3. Examine overall patterns
print(all_plots$dotplot)
print(all_plots$heatmap)

# 4. Dive into specific cell types
print(all_plots$violin[["CD4+ T cells"]])
print(all_plots$umap[["CD8+ T cells"]])
print(all_plots$violin[["B cells"]])

# 5. Quick check for any cell type
quick_plots <- quick_marker_viz(
    seurat_obj,
    cell_type = "NK cells",
    tissue_type = "Immune system"
)
```

### SingleCellExperiment Workflow

```r
library(SingleCellExperiment)
library(scater)
library(ggplot2)

# 1. Run ScType annotation
source("R/sctype_wrapper_sce.R")
sce <- run_sctype_sce(
    sce,
    known_tissue_type = "Immune system",
    custom_marker_file = "ScTypeDB_full.xlsx",
    plot = TRUE
)

# 2. Visualize all markers
source("R/sctype_visualize_sce.R")
all_plots <- visualize_sctype_markers_sce(
    sce,
    annotation_col = "sctype_classification",
    tissue_type = "Immune system",
    assay_name = "logcounts",
    top_n = 5,
    save_plots = TRUE,
    output_dir = "immune_markers"
)

# 3. Examine results
print(all_plots$dotplot)
print(all_plots$heatmap)
```

## Hierarchical Annotations

For hierarchical annotations (broad + fine), visualize both levels:

```r
# Run hierarchical annotation
source("R/sctype_hierarchical.R")
seurat_obj <- run_sctype_hierarchical(
    seurat_obj,
    known_tissue_type = "Immune system",
    custom_marker_file = "ScTypeDB_hierarchical.xlsx"
)

# Visualize broad level
plots_broad <- visualize_sctype_markers(
    seurat_obj,
    annotation_col = "sctype_broad",
    database_file = "ScTypeDB_hierarchical.xlsx",
    tissue_type = "Immune system",
    top_n = 3
)

# Visualize fine level
plots_fine <- visualize_sctype_markers(
    seurat_obj,
    annotation_col = "sctype_fine",
    database_file = "ScTypeDB_hierarchical.xlsx",
    tissue_type = "Immune system",
    top_n = 5
)

# Compare dotplots
library(patchwork)
plots_broad$dotplot / plots_fine$dotplot
```

## Customization

### Custom Number of Markers

```r
# Show more markers per cell type
plots <- visualize_sctype_markers(
    seurat_obj,
    tissue_type = "Brain",
    top_n = 10  # Show top 10 markers instead of 5
)
```

### Custom Databases

```r
# Use custom marker database
plots <- visualize_sctype_markers(
    seurat_obj,
    database_file = "path/to/custom_markers.xlsx",
    tissue_type = "My Custom Tissue",
    top_n = 8
)
```

### Selective Plot Generation

```r
# Generate only specific plot types
plots <- visualize_sctype_markers(
    seurat_obj,
    tissue_type = "Immune system",
    plot_types = c("dotplot", "heatmap")  # Skip violin and UMAP
)
```

### Save High-Resolution Plots

```r
# Save plots with custom settings
plots <- visualize_sctype_markers(
    seurat_obj,
    tissue_type = "Immune system",
    save_plots = TRUE,
    output_dir = "publication_figures"
)

# Plots are saved as PNG files:
# - violin_[CellType].png (12" x 8", 300 dpi)
# - umap_[CellType].png (14" x 10", 300 dpi)
# - dotplot_all_markers.png (16" x 10", 300 dpi)
# - heatmap_all_markers.png (14" x 10", 300 dpi)
```

## Troubleshooting

### "UMAP not found"

**Error:** `UMAP not found in reducedDims. Skipping UMAP plots.`

**Solution:**
```r
# For Seurat
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# For SingleCellExperiment
sce <- runUMAP(sce, dimred = "PCA")
```

### "No markers found"

**Error:** `No markers found for tissue type: [tissue]`

**Solution:** Check tissue type spelling. Available tissues:
```r
library(openxlsx)
db <- read.xlsx("ScTypeDB_full.xlsx")
unique(db$tissueType)
```

### "Package not installed" warnings

**Solution:** Install missing packages:
```r
# For heatmaps
BiocManager::install("ComplexHeatmap")

# For plot arrangement
install.packages("patchwork")

# For SCE violin/UMAP plots
BiocManager::install("scater")
```

### Empty plots

**Issue:** Plots are generated but show no data

**Possible causes:**
1. Markers not present in dataset (different gene naming)
2. All cells annotated as "Unknown"
3. Wrong assay specified

**Solutions:**
```r
# Check available genes
head(rownames(seurat_obj))  # Seurat
head(rownames(sce))          # SCE

# Check annotations
table(seurat_obj$sctype_classification)
table(colData(sce)$sctype_classification)

# Try different assay
plots <- visualize_sctype_markers(seurat_obj, assay = "SCT")
plots <- visualize_sctype_markers_sce(sce, assay_name = "counts")
```

## Advanced Usage

### Combining Multiple Visualizations

```r
library(patchwork)

# Generate plots for two cell types
plots_cd4 <- quick_marker_viz(seurat_obj, "CD4+ T cells", tissue_type = "Immune system")
plots_cd8 <- quick_marker_viz(seurat_obj, "CD8+ T cells", tissue_type = "Immune system")

# Combine violin plots side-by-side
plots_cd4$violin | plots_cd8$violin

# Combine UMAP plots
plots_cd4$umap / plots_cd8$umap
```

### Custom Marker Selection

If you want to visualize specific markers (not from database):

```r
# Manual marker visualization using Seurat functions
library(Seurat)

custom_markers <- c("CD4", "CD8A", "CD3D", "CD19", "CD14")

VlnPlot(seurat_obj, features = custom_markers, group.by = "sctype_classification")
FeaturePlot(seurat_obj, features = custom_markers)
DotPlot(seurat_obj, features = custom_markers, group.by = "sctype_classification")
```

### Exporting Plot Data

```r
# Get plot data for custom analysis
plots <- visualize_sctype_markers(seurat_obj, tissue_type = "Immune system")

# Extract data from dotplot
dotplot_data <- plots$dotplot$data

# Save as table
write.csv(dotplot_data, "marker_expression_summary.csv")
```

## Best Practices

1. **Always validate annotations** - Use these visualizations to confirm ScType results
2. **Check both positive and negative markers** - Negative markers are equally important
3. **Compare hierarchical levels** - Visualize both broad and fine annotations
4. **Save plots early** - Use `save_plots = TRUE` to preserve results
5. **Customize top_n** - Adjust based on marker availability (3-10 is typical)
6. **Use dotplots for overview** - Best for comparing all cell types at once
7. **Use UMAP for spatial patterns** - Reveals unexpected co-localization
8. **Use heatmaps for clustering** - Identifies related cell types

## Citation

If you use these visualization functions in your research, please cite:

Ianevski, A., Giri, A.K. & Aittokallio, T. Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data. Nat Commun 13, 1246 (2022). https://doi.org/10.1038/s41467-022-28803-w

---

*Last Updated: 2025-11-15*
