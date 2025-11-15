# ScType for SingleCellExperiment Objects

This guide explains how to use ScType with **SingleCellExperiment (SCE)** objects from the Bioconductor ecosystem.

## Overview

ScType now provides full support for SingleCellExperiment objects through parallel implementations of all major functions:

- `run_sctype_sce()` - Basic cell type annotation
- `run_sctype_hierarchical_sce()` - Two-level hierarchical annotation (broad + fine)
- `get_cluster_hierarchy_sce()` - Query annotations for specific clusters
- `print_hierarchy_table_sce()` - Display formatted annotation table
- `sctype_source_sce()` - Load ScType functions and return database URLs

## Quick Start

### Installation

```r
# Install required Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SingleCellExperiment", "scater", "scran"))

# Install CRAN packages
install.packages(c("dplyr", "HGNChelper", "openxlsx", "ggplot2"))
```

### Basic Usage

```r
library(SingleCellExperiment)
library(scater)

# Load your SingleCellExperiment object
# sce <- readRDS("my_sce_object.rds")

# Ensure you have clustering performed
# Assuming clusters are in colData(sce)$cluster

# Run ScType annotation
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper_sce.R")

sce <- run_sctype_sce(
    sce_object = sce,
    known_tissue_type = "Immune system",
    assay_name = "logcounts",
    cluster_col = "cluster",
    plot = TRUE
)

# View results
table(colData(sce)$sctype_classification)
```

### Hierarchical Annotation

```r
# Load hierarchical annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_hierarchical_sce.R")

# Run hierarchical annotation (broad + fine levels)
sce <- run_sctype_hierarchical_sce(
    sce_object = sce,
    known_tissue_type = "Immune system",
    assay_name = "logcounts",
    cluster_col = "cluster",
    plot = TRUE,
    broad_name = "sctype_broad",
    fine_name = "sctype_fine"
)

# View both annotation levels
print(table(colData(sce)$sctype_broad))   # Broad categories
print(table(colData(sce)$sctype_fine))    # Fine subtypes

# Display hierarchical table
print_hierarchy_table_sce(sce)

# Query specific cluster
cluster_0_info <- get_cluster_hierarchy_sce(sce, cluster_id = 0)
print(cluster_0_info)
```

## Complete Workflow Example

```r
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)

# 1. Load or create SingleCellExperiment object
# For this example, let's assume you have raw counts
# sce <- SingleCellExperiment(assays = list(counts = raw_counts))

# 2. Quality control (example - adjust thresholds for your data)
sce <- addPerCellQC(sce)
sce <- sce[, sce$total > 1000 & sce$detected > 500]

# 3. Normalization
sce <- logNormCounts(sce)

# 4. Feature selection
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, n = 2000)

# 5. Dimensionality reduction
sce <- runPCA(sce, subset_row = hvg, ncomponents = 50)
sce <- runUMAP(sce, dimred = "PCA")

# 6. Clustering
g <- buildSNNGraph(sce, use.dimred = "PCA", k = 20)
clusters <- igraph::cluster_louvain(g)$membership
colData(sce)$cluster <- clusters

# 7. Run hierarchical ScType annotation
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_hierarchical_sce.R")

sce <- run_sctype_hierarchical_sce(
    sce_object = sce,
    known_tissue_type = NULL,  # Auto-detect tissue type
    assay_name = "logcounts",
    scaled = TRUE,
    cluster_col = "cluster",
    plot = TRUE
)

# 8. Visualize results
library(patchwork)

# Clusters
p1 <- plotUMAP(sce, colour_by = "cluster", text_by = "cluster") +
      ggtitle("Clusters")

# Broad cell types
p2 <- plotUMAP(sce, colour_by = "sctype_broad") +
      ggtitle("Broad Cell Categories")

# Fine cell types
p3 <- plotUMAP(sce, colour_by = "sctype_fine") +
      ggtitle("Fine Cell Subtypes")

# Combined plot
(p1 | p2) / p3

# 9. View annotation summary
print_hierarchy_table_sce(sce)
```

## Function Reference

### run_sctype_sce()

Basic cell type annotation for SingleCellExperiment objects.

**Parameters:**
- `sce_object`: SingleCellExperiment object with clustering
- `known_tissue_type`: Tissue type (e.g., "Immune system", "Brain"). NULL for auto-detection.
- `assay_name`: Assay to use (default: "logcounts")
- `scaled`: Whether data is scaled (default: TRUE)
- `cluster_col`: Column in colData with clusters (default: "cluster")
- `custom_marker_file`: Path to marker database (default: ScTypeDB_full.xlsx from GitHub)
- `plot`: Generate UMAP plot (default: FALSE)
- `name`: Column name for results (default: "sctype_classification")

**Returns:** Modified SCE object with new colData column

**Example:**
```r
sce <- run_sctype_sce(sce, known_tissue_type = "Immune system", plot = TRUE)
```

### run_sctype_hierarchical_sce()

Two-level hierarchical cell type annotation (broad categories + fine subtypes).

**Parameters:**
- `sce_object`: SingleCellExperiment object with clustering
- `known_tissue_type`: Tissue type. NULL for auto-detection.
- `assay_name`: Assay to use (default: "logcounts")
- `scaled`: Whether data is scaled (default: TRUE)
- `cluster_col`: Column in colData with clusters (default: "cluster")
- `custom_marker_file`: Path to hierarchical database (default: ScTypeDB_hierarchical.xlsx)
- `plot`: Generate side-by-side UMAP plots (default: FALSE)
- `broad_name`: Column name for broad categories (default: "sctype_broad")
- `fine_name`: Column name for fine subtypes (default: "sctype_fine")

**Returns:** Modified SCE object with two new colData columns

**Example:**
```r
sce <- run_sctype_hierarchical_sce(sce, known_tissue_type = "Brain", plot = TRUE)
```

### get_cluster_hierarchy_sce()

Query both broad and fine annotations for a specific cluster.

**Parameters:**
- `sce_object`: Annotated SingleCellExperiment object
- `cluster_id`: Cluster ID to query
- `cluster_col`: Column with clusters (default: "cluster")
- `broad_name`: Column with broad annotations (default: "sctype_broad")
- `fine_name`: Column with fine annotations (default: "sctype_fine")

**Returns:** List with cluster information and both annotation levels

**Example:**
```r
info <- get_cluster_hierarchy_sce(sce, cluster_id = 0)
# Returns: list(cluster = 0, broad_category = "T cells",
#               fine_subtype = "CD4+ T cells", n_cells = 500, is_refined = TRUE)
```

### print_hierarchy_table_sce()

Display formatted table of all hierarchical annotations.

**Parameters:**
- `sce_object`: Annotated SingleCellExperiment object
- `cluster_col`: Column with clusters (default: "cluster")
- `broad_name`: Column with broad annotations (default: "sctype_broad")
- `fine_name`: Column with fine annotations (default: "sctype_fine")

**Returns:** Invisible data.frame with the table (also prints to console)

**Example:**
```r
print_hierarchy_table_sce(sce)
```

## Key Differences from Seurat

| Aspect | Seurat | SingleCellExperiment |
|--------|--------|---------------------|
| **Data storage** | `seurat_obj[["RNA"]]@scale.data` or `$scale.data` | `assay(sce, "logcounts")` |
| **Metadata** | `seurat_obj@meta.data` | `colData(sce)` |
| **Reduced dimensions** | `Reductions(seurat_obj)$umap@cell.embeddings` | `reducedDim(sce, "UMAP")` |
| **Default assay** | "RNA" (automatic) | "logcounts" or "counts" (must specify) |
| **Cluster column** | `seurat_clusters` (created automatically) | Must specify in `cluster_col` parameter |
| **Plotting** | `DimPlot()` from Seurat | `plotUMAP()` from scater |

## Tissue Types Supported

ScType databases include markers for:
- **Immune system** - 44 cell types (enhanced database)
- **Brain** - 17 cell types (neurons, glia, etc.)
- **Liver** - 8 cell types
- **Pancreas** - 9 cell types
- **Kidney** - 10 cell types
- **Lung** - 11 cell types
- **Heart** - 8 cell types
- **Intestine** - 6 cell types
- **Muscle** - 3 cell types
- **Skin** - 4 cell types
- **Adipose** - 2 cell types

Plus additional tissues in the full database.

## Available Databases

```r
# Load database URLs
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper_sce.R")
db_urls <- sctype_source_sce()

# Available databases:
# - full: Original ScType database (comprehensive)
# - short: Abbreviated database (quick annotation)
# - enhanced: High-confidence markers from multiple sources (122 cell types)
# - hierarchical: Two-level annotation database (broad + fine)
```

## Troubleshooting

### "Cluster column not found"

**Error:** `Cluster column 'cluster' not found in colData`

**Solution:** Perform clustering first or specify the correct column name:
```r
# Option 1: Perform clustering
g <- buildSNNGraph(sce, use.dimred = "PCA")
colData(sce)$cluster <- igraph::cluster_louvain(g)$membership

# Option 2: Specify correct column
sce <- run_sctype_sce(sce, cluster_col = "seurat_clusters")
```

### "Assay not found"

**Error:** `Assay 'logcounts' not found`

**Solution:** Normalize data or specify correct assay:
```r
# Option 1: Normalize
sce <- logNormCounts(sce)

# Option 2: Use raw counts (less accurate)
sce <- run_sctype_sce(sce, assay_name = "counts", scaled = FALSE)
```

### "UMAP not found" (when plot = TRUE)

**Warning:** `UMAP not found in reducedDims. Skipping plot.`

**Solution:** Run UMAP before annotation:
```r
sce <- runUMAP(sce, dimred = "PCA")
```

### All clusters assigned as "Unknown"

**Possible causes:**
1. Wrong tissue type selected
2. Using wrong assay (raw counts instead of logcounts)
3. Species mismatch (human markers on mouse data)

**Solutions:**
```r
# Auto-detect tissue type
sce <- run_sctype_sce(sce, known_tissue_type = NULL)

# Ensure using scaled data
sce <- run_sctype_sce(sce, assay_name = "logcounts", scaled = TRUE)

# For mouse data, use mouse-compatible databases
# (Gene symbols should be uppercase in both data and database)
```

## Advanced Usage

### Custom Marker Database

```r
# Create custom Excel file with columns:
# - tissueType
# - broadCategory (for hierarchical)
# - cellName
# - geneSymbolmore1 (positive markers, comma-separated)
# - geneSymbolmore2 (negative markers, comma-separated)

sce <- run_sctype_hierarchical_sce(
    sce,
    custom_marker_file = "path/to/custom_markers.xlsx",
    known_tissue_type = "My Custom Tissue"
)
```

### Using with Specific Assays

```r
# Use counts (raw) instead of logcounts
sce <- run_sctype_sce(sce, assay_name = "counts", scaled = FALSE)

# Use normalized but not log-transformed data
sce <- run_sctype_sce(sce, assay_name = "normcounts", scaled = FALSE)
```

### Integration with Other Bioconductor Tools

```r
library(SingleR)
library(celldex)

# Compare ScType with SingleR
ref <- celldex::HumanPrimaryCellAtlasData()
singler_results <- SingleR(test = sce, ref = ref, labels = ref$label.main)

colData(sce)$SingleR <- singler_results$labels

# Compare annotations
table(ScType = colData(sce)$sctype_classification,
      SingleR = colData(sce)$SingleR)
```

## Citation

If you use ScType with SingleCellExperiment in your research, please cite:

Ianevski, A., Giri, A.K. & Aittokallio, T. Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data. Nat Commun 13, 1246 (2022). https://doi.org/10.1038/s41467-022-28803-w

## Support

- **GitHub Issues**: https://github.com/IanevskiAleksandr/sc-type/issues
- **Email**: aleksandr.ianevski@helsinki.fi
- **Documentation**: See CLAUDE.md for detailed technical documentation

---

*Last Updated: 2025-11-15*
