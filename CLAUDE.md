# CLAUDE.md - ScType Repository Guide for AI Assistants

## Project Overview

**ScType** is a computational method for fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data (scRNA-seq).

- **Publication**: [Nature Communications (2022)](https://doi.org/10.1038/s41467-022-28803-w)
- **Primary Author**: Aleksandr Ianevski (aleksandr.ianevski@helsinki.fi)
- **License**: GNU General Public License v3.0
- **Web Portal**: http://sctype.app
- **Repository**: https://github.com/IanevskiAleksandr/sc-type

### Purpose

ScType enables automated selection and annotation of cell types from scRNA-seq data by:
1. Using positive and negative marker gene sets for each cell type
2. Computing enrichment scores for each cell based on marker expression
3. Assigning cell types to clusters based on aggregated scores
4. Supporting multiple tissue types and custom marker databases

---

## Repository Structure

```
sc-type/
├── R/                           # Core R functions
│   ├── gene_sets_prepare.R      # Gene set preparation from marker database
│   ├── sctype_score_.R          # ScType scoring algorithm
│   ├── auto_detect_tissue_type.R # Tissue type auto-detection
│   ├── sctype_wrapper.R         # Convenience wrapper for Seurat objects
│   ├── sctype_wrapper_sce.R     # Convenience wrapper for SingleCellExperiment
│   ├── sctype_hierarchical.R    # Hierarchical annotation (Seurat)
│   └── sctype_hierarchical_sce.R # Hierarchical annotation (SingleCellExperiment)
├── ScTypeDB_full.xlsx           # Complete cell marker database
├── ScTypeDB_short.xlsx          # Abbreviated marker database
├── ScTypeDB_enhanced.xlsx       # Enhanced marker database (122 cell types)
├── ScTypeDB_hierarchical.xlsx   # Hierarchical marker database (broad + fine)
├── exampleData.RDS              # Example scRNA-seq dataset (PBMC 3k)
├── filtered_gene_bc_matrices/   # Example 10X Genomics data
├── livercellatlass/             # Web interface for liver cell atlas
│   ├── css/                     # Stylesheets
│   ├── js/                      # JavaScript files
│   └── *.html                   # HTML pages
├── index.html                   # Main web interface
├── index2.html                  # Alternative web interface
├── images/                      # Image assets
├── fonts/                       # Font files
├── fig*.png                     # Documentation figures
├── README.md                    # User documentation
└── LICENSE                      # GNU GPL v3.0 license
```

---

## Core Components

### 1. Gene Sets Preparation (`R/gene_sets_prepare.R`)

**Purpose**: Prepares positive and negative gene sets from Excel marker databases.

**Function**: `gene_sets_prepare(path_to_db_file, cell_type)`

**Parameters**:
- `path_to_db_file`: Path to Excel file with cell type markers
- `cell_type`: Tissue type (e.g., "Immune system", "Liver", "Pancreas", "Brain")

**Process**:
1. Reads Excel database filtered by tissue type
2. Validates and corrects gene symbols using `HGNChelper::checkGeneSymbols()`
3. Returns list with `gs_positive` (marker genes) and `gs_negative` (genes not expressed)

**Database Format** (Excel columns):
- `tissueType`: Tissue/organ classification
- `cellName`: Cell type name
- `geneSymbolmore1`: Positive marker genes (comma-separated)
- `geneSymbolmore2`: Negative marker genes (comma-separated)

### 2. ScType Scoring (`R/sctype_score_.R`)

**Purpose**: Calculates enrichment scores for each cell type in each cell.

**Function**: `sctype_score(scRNAseqData, scaled, gs, gs2, gene_names_to_uppercase, ...)`

**Parameters**:
- `scRNAseqData`: Expression matrix (genes × cells)
- `scaled`: Boolean indicating if data is z-scaled (default: TRUE)
- `gs`: List of positive marker gene sets
- `gs2`: List of negative marker gene sets (NULL if not applicable)
- `gene_names_to_uppercase`: Convert gene names to uppercase (default: TRUE)

**Algorithm**:
1. Calculate marker sensitivity scores (rescaled based on gene frequency)
2. Optionally z-scale the expression matrix
3. Weight expression by marker sensitivity
4. Compute combined scores: `score = sum(positive_markers)/sqrt(n) - sum(negative_markers)/sqrt(n)`
5. Return matrix of cell types × cells with enrichment scores

### 3. Tissue Type Auto-Detection (`R/auto_detect_tissue_type.R`)

**Purpose**: Automatically detects the most likely tissue type for a dataset.

**Function**: `auto_detect_tissue_type(path_to_db_file, seuratObject, scaled, assay, ...)`

**Process**:
1. Iterates through all tissue types in database
2. Runs ScType scoring for each tissue
3. Aggregates scores by cluster
4. Returns tissues ranked by mean score
5. Generates barplot visualization

**Usage**: When tissue type is unknown, use this before cell type annotation.

### 4. Wrapper Function (`R/sctype_wrapper.R`)

**Purpose**: Simplified interface for ScType analysis on Seurat objects.

**Main Functions**:
- `sctype_source()`: Loads all ScType functions and returns database URL
- `run_sctype()`: Complete workflow from Seurat object to annotated results

**Function**: `run_sctype(seurat_object, known_tissue_type, assay, scaled, custom_marker_file, plot, name)`

**Parameters**:
- `seurat_object`: Seurat object (v4 or v5 compatible)
- `known_tissue_type`: Tissue type (optional, auto-detected if NULL)
- `assay`: Assay name (default: "RNA")
- `scaled`: Use scaled data (default: TRUE)
- `custom_marker_file`: Custom marker database (optional)
- `plot`: Generate UMAP plot (default: FALSE)
- `name`: Metadata column name (default: "sctype_classification")

**Returns**: Modified Seurat object with new metadata column containing cell type annotations.

### 5. Hierarchical Annotation (`R/sctype_hierarchical.R`)

**Purpose**: Provides two-level hierarchical cell type annotation with both broad categories and fine-grained subtypes.

**Main Function**: `run_sctype_hierarchical(seurat_object, known_tissue_type, assay, scaled, custom_marker_file, plot, broad_name, fine_name)`

**Parameters**:
- `seurat_object`: Seurat object (v4 or v5 compatible)
- `known_tissue_type`: Tissue type (optional, auto-detected if NULL)
- `assay`: Assay name (default: "RNA")
- `scaled`: Use scaled data (default: TRUE)
- `custom_marker_file`: Path to hierarchical marker database (default: ScTypeDB_hierarchical.xlsx)
- `plot`: Generate side-by-side UMAP plots (default: FALSE)
- `broad_name`: Metadata column name for broad categories (default: "sctype_broad")
- `fine_name`: Metadata column name for fine subtypes (default: "sctype_fine")

**Process**:
1. **Step 1 - Broad Annotation**: Aggregates all markers within each broad category (e.g., "T cells", "B cells", "Neurons") and runs ScType scoring at the broad level
2. **Step 2 - Fine Annotation**: Runs ScType scoring at the fine level using individual cell type markers (e.g., "CD4+ T cells", "CD8+ T cells")
3. **Confidence-Based Assignment**: If fine-level confidence is low (score < ncells/4), falls back to broad category assignment
4. **Dual Metadata Columns**: Adds both `sctype_broad` and `sctype_fine` columns to the Seurat object

**Returns**: Modified Seurat object with two new metadata columns for hierarchical annotations.

**Example Hierarchical Annotations**:
- Broad: "T cells" → Fine: "CD4+ T cells"
- Broad: "Neurons" → Fine: "Excitatory neurons"
- Broad: "Endothelial" → Fine: "Arterial endothelial cells"
- Broad: "Unknown" → Fine: "Unknown" (low confidence at both levels)

**Helper Functions**:
- `get_cluster_hierarchy(seurat_object, cluster_id, broad_name, fine_name)`: Query both annotations for a specific cluster
- `print_hierarchy_table(seurat_object, broad_name, fine_name)`: Display formatted table of all hierarchical annotations

**Database Structure**:
The hierarchical database (ScTypeDB_hierarchical.xlsx) includes an additional `broadCategory` column:
- `tissueType`: Tissue/organ classification
- `broadCategory`: Broad cell type category (NEW)
- `cellName`: Fine cell type name
- `geneSymbolmore1`: Positive marker genes
- `geneSymbolmore2`: Negative marker genes

---

## Supported Tissue Types

The ScTypeDB database includes markers for:
- Immune system
- Liver
- Pancreas
- Kidney
- Eye
- Brain
- Lung
- Adrenal
- Heart
- Intestine
- Muscle
- Placenta
- Spleen
- Stomach
- Thymus

---

## Development Workflows

### Standard Usage Pattern

```r
# 1. Load dependencies
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)

# 2. Source ScType functions
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# 3. Prepare gene sets
gs_list <- gene_sets_prepare("path/to/database.xlsx", "Immune system")

# 4. Extract scaled data from Seurat object
seurat_v5 <- isFALSE('counts' %in% names(attributes(seurat_obj[["RNA"]])))
scRNAseqData <- if (seurat_v5) as.matrix(seurat_obj[["RNA"]]$scale.data)
                else as.matrix(seurat_obj[["RNA"]]@scale.data)

# 5. Calculate scores
es.max <- sctype_score(scRNAseqData, scaled = TRUE,
                       gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# 6. Aggregate by cluster and assign cell types
cL_results <- do.call("rbind", lapply(unique(seurat_obj@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[, rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters==cl, ])]),
                     decreasing = TRUE)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl,
                    ncells = sum(seurat_obj@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores <- cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# 7. Filter low-confidence assignments
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
```

### Simplified Wrapper Approach

```r
# Load wrapper
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R")

# Run complete analysis
seurat_obj <- run_sctype(seurat_obj,
                         known_tissue_type = "Immune system",
                         custom_marker_file = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",
                         name = "sctype_classification",
                         plot = TRUE)
```

### Hierarchical Annotation Workflow

**Purpose**: Obtain both broad and fine-grained cell type annotations in separate metadata columns.

```r
# Load hierarchical annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_hierarchical.R")

# Run hierarchical annotation
seurat_obj <- run_sctype_hierarchical(
    seurat_object = seurat_obj,
    known_tissue_type = "Immune system",  # or NULL for auto-detection
    custom_marker_file = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_hierarchical.xlsx",
    assay = "RNA",
    scaled = TRUE,
    plot = TRUE,  # Shows side-by-side UMAP plots
    broad_name = "sctype_broad",  # Column name for broad categories
    fine_name = "sctype_fine"     # Column name for fine subtypes
)

# View results
print(table(seurat_obj@meta.data$sctype_broad))   # Broad categories
print(table(seurat_obj@meta.data$sctype_fine))    # Fine subtypes

# Display hierarchical annotation table
print_hierarchy_table(seurat_obj)

# Query specific cluster
cluster_info <- get_cluster_hierarchy(seurat_obj, cluster_id = 0)
print(cluster_info)
# Output: list(cluster = 0, broad_category = "T cells",
#              fine_subtype = "CD4+ T cells", n_cells = 500, is_refined = TRUE)

# Visualize both levels
library(patchwork)
p1 <- DimPlot(seurat_obj, group.by = "sctype_broad", label = TRUE, repel = TRUE) +
      ggtitle("Broad Cell Categories")
p2 <- DimPlot(seurat_obj, group.by = "sctype_fine", label = TRUE, repel = TRUE) +
      ggtitle("Fine Cell Subtypes")
p1 / p2
```

**Key Features**:
- **Automatic fallback**: If fine-level annotation has low confidence, uses broad category instead
- **Dual visualization**: See both annotation levels side-by-side
- **Flexible naming**: Customize metadata column names
- **Query functions**: Easily retrieve annotations for specific clusters
- **Comprehensive table**: View all cluster assignments in formatted table

**Use Cases**:
1. **Initial exploration**: Use broad categories to understand major cell populations
2. **Detailed analysis**: Examine fine subtypes for specific cell lineages of interest
3. **Quality control**: Identify clusters with confident fine annotations (where broad ≠ fine)
4. **Publication**: Report both levels to provide context and specificity

### SingleCellExperiment Support

**Purpose**: Use ScType with SingleCellExperiment objects instead of Seurat.

ScType now provides full support for **SingleCellExperiment (SCE)** objects through parallel implementations of all major functions. This allows users who prefer Bioconductor workflows to use ScType seamlessly.

#### Basic SCE Annotation

```r
# Load required packages
library(SingleCellExperiment)
library(scater)

# Load wrapper
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper_sce.R")

# Run ScType annotation
sce <- run_sctype_sce(
    sce_object = sce,
    known_tissue_type = "Immune system",  # or NULL for auto-detection
    assay_name = "logcounts",             # Use "counts" for raw data
    scaled = TRUE,                         # TRUE if using logcounts
    cluster_col = "cluster",               # Column in colData with clusters
    custom_marker_file = NULL,            # NULL uses default database
    plot = TRUE,                           # Requires UMAP in reducedDims
    name = "sctype_classification"         # Column name for results
)

# View results
table(colData(sce)$sctype_classification)
```

#### Hierarchical SCE Annotation

```r
# Load hierarchical function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_hierarchical_sce.R")

# Run hierarchical annotation
sce <- run_sctype_hierarchical_sce(
    sce_object = sce,
    known_tissue_type = "Immune system",
    assay_name = "logcounts",
    scaled = TRUE,
    cluster_col = "cluster",
    custom_marker_file = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_hierarchical.xlsx",
    plot = TRUE,
    broad_name = "sctype_broad",
    fine_name = "sctype_fine"
)

# View both annotation levels
print(table(colData(sce)$sctype_broad))
print(table(colData(sce)$sctype_fine))

# Display hierarchical table
print_hierarchy_table_sce(sce)

# Query specific cluster
cluster_info <- get_cluster_hierarchy_sce(sce, cluster_id = 0)
```

#### Key Differences from Seurat

| Aspect | Seurat | SingleCellExperiment |
|--------|--------|---------------------|
| **Data access** | `seurat_obj[["RNA"]]@scale.data` | `assay(sce, "logcounts")` |
| **Metadata** | `seurat_obj@meta.data` | `colData(sce)` |
| **Reduced dims** | `Reductions(seurat_obj)` | `reducedDim(sce, "UMAP")` |
| **Default assay** | "RNA" | "logcounts" or "counts" |
| **Cluster column** | `seurat_clusters` (automatic) | Must specify `cluster_col` |

#### SCE Workflow Example

```r
# Complete SCE workflow with hierarchical annotation
library(SingleCellExperiment)
library(scater)
library(scran)

# Assuming you have a SingleCellExperiment object 'sce'

# 1. Normalize if not already done
sce <- logNormCounts(sce)

# 2. Perform clustering if not already done
sce <- runPCA(sce)
sce <- runUMAP(sce)
clusters <- clusterCells(sce, use.dimred = "PCA")
colData(sce)$cluster <- clusters

# 3. Run hierarchical ScType annotation
source("R/sctype_hierarchical_sce.R")
sce <- run_sctype_hierarchical_sce(sce, known_tissue_type = "Immune system", plot = TRUE)

# 4. Visualize results
library(ggplot2)
plotUMAP(sce, colour_by = "sctype_broad")
plotUMAP(sce, colour_by = "sctype_fine")
```

#### Function Correspondence

| Seurat Version | SingleCellExperiment Version |
|---------------|------------------------------|
| `run_sctype()` | `run_sctype_sce()` |
| `run_sctype_hierarchical()` | `run_sctype_hierarchical_sce()` |
| `get_cluster_hierarchy()` | `get_cluster_hierarchy_sce()` |
| `print_hierarchy_table()` | `print_hierarchy_table_sce()` |
| `sctype_source()` | `sctype_source_sce()` |

---

## Key Conventions and Design Patterns

### 1. Seurat Version Compatibility

ScType supports both Seurat v4 and v5 with different data access patterns:

**Detection Pattern**:
```r
# Check if Seurat v5
seurat_v5 <- isFALSE('counts' %in% names(attributes(seurat_obj[["RNA"]])))

# Version-specific data extraction
if (seurat_v5) {
    data <- seurat_obj[["RNA"]]$scale.data  # v5 uses $
} else {
    data <- seurat_obj[["RNA"]]@scale.data  # v4 uses @
}
```

### 2. Gene Symbol Validation

Always validate gene symbols using `HGNChelper::checkGeneSymbols()`:
- Converts to uppercase
- Corrects outdated symbols
- Removes NA values
- Ensures unique symbols

### 3. Score Calculation Logic

**Positive contribution**: `sum(expression_of_positive_markers) / sqrt(n_positive_markers)`
**Negative contribution**: `sum(expression_of_negative_markers) / sqrt(n_negative_markers) * -1`
**Final score**: positive + negative

The square root normalization prevents bias toward cell types with many markers.

### 4. Low-Confidence Filtering

Default threshold: `score < ncells_in_cluster / 4`
- If cluster score is too low relative to cluster size, assign "Unknown"
- This prevents false positives in ambiguous clusters

### 5. Marker Sensitivity Weighting

Genes are weighted by inverse frequency across all cell types:
- Rare markers (specific to few cell types) get higher weights
- Common markers (present in many cell types) get lower weights
- Rescaled to [0, 1] range

---

## Testing and Validation

### Example Datasets

1. **PBMC 3k** (Seurat tutorial dataset):
   - Located in `filtered_gene_bc_matrices/hg19/`
   - Available at: https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

2. **Example RDS** (`exampleData.RDS`):
   - Pre-processed scaled scRNA-seq matrix
   - Can be loaded directly for testing

### Validation Approach

When modifying core functions:
1. Load example data
2. Run ScType with standard parameters
3. Verify expected cell types are identified (e.g., T cells, B cells, Monocytes in PBMC)
4. Check that scores are reasonable (typically between -10 and 50)
5. Ensure low-confidence clusters are marked as "Unknown"

---

## Common Pitfalls and Solutions

### 1. Matrix Format Issues

**Problem**: ScType expects genes as rows, cells as columns
**Solution**: Transpose if needed: `t(matrix)`

### 2. Gene Name Mismatches

**Problem**: Mouse genes vs human genes, or lowercase vs uppercase
**Solution**:
- Set `gene_names_to_uppercase = TRUE` in `sctype_score()`
- Use species-specific marker databases
- Validate with `HGNChelper`

### 3. Unscaled Data

**Problem**: Using raw counts instead of scaled data
**Solution**:
- If using raw counts, set `scaled = FALSE` in `sctype_score()`
- ScType will z-scale internally
- For best results, use pre-scaled data from Seurat's `ScaleData()`

### 4. Missing Negative Markers

**Problem**: `gs2` parameter errors when negative markers absent
**Solution**: Set `gs2 = NULL` explicitly when no negative markers available

### 5. Empty Clusters

**Problem**: Clusters with no cells cause errors
**Solution**: Filter `seurat_clusters` to only include non-empty clusters before aggregation

---

## Custom Marker Databases

### Creating Custom Databases

To use custom markers, create an Excel file with these columns:

| Column | Description | Example |
|--------|-------------|---------|
| tissueType | Tissue/organ | "Immune system" |
| cellName | Cell type name | "CD8+ T cells" |
| geneSymbolmore1 | Positive markers (comma-separated) | "CD8A,CD8B,CD3D,CD3E" |
| geneSymbolmore2 | Negative markers (comma-separated) | "CD4,CD19,CD14" |
| shortName | Abbreviated name (optional) | "CD8 T" |

**Best Practices**:
- Use well-validated, specific markers
- Include negative markers to improve specificity
- 3-10 positive markers per cell type is optimal
- Avoid markers expressed in many cell types
- Use official gene symbols (HGNC for human, MGI for mouse)

### Extending Existing Database

To add new cell types:
1. Open `ScTypeDB_full.xlsx`
2. Add rows with appropriate tissue type and markers
3. Ensure gene symbols are comma-separated without spaces
4. Save and use updated file path in `gene_sets_prepare()`

---

## Web Interface Components

### File Locations

- **Main portal**: `index.html`, `index2.html`
- **Liver atlas**: `livercellatlass/` directory
- **Styling**: `css/style.css`, `livercellatlass/css/`
- **Scripts**: `livercellatlass/js/`

### Purpose

The web interfaces provide:
- Interactive cell type annotation without R
- Visualization of results
- Database browsing
- Quality control metrics (liver atlas)

**Note**: Web interfaces are separate from core R functionality. R functions can be used independently.

---

## Integration with Analysis Pipelines

### Seurat Pipeline Integration

ScType fits into Seurat workflow after clustering:

```r
# Standard Seurat workflow
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Insert ScType here
seurat_obj <- run_sctype(seurat_obj, known_tissue_type = "Immune system")

# Continue with downstream analysis
DimPlot(seurat_obj, group.by = "sctype_classification")
```

### SingleCellExperiment Pipeline Integration

ScType integrates seamlessly with Bioconductor's SingleCellExperiment workflow:

```r
# Standard SingleCellExperiment/scater workflow
library(SingleCellExperiment)
library(scater)
library(scran)

# Normalize and preprocess
sce <- logNormCounts(sce)
sce <- runPCA(sce, ncomponents = 50)

# Clustering
g <- buildSNNGraph(sce, use.dimred = "PCA")
clusters <- igraph::cluster_louvain(g)$membership
colData(sce)$cluster <- clusters

# Dimensionality reduction for visualization
sce <- runUMAP(sce, dimred = "PCA")
sce <- runTSNE(sce, dimred = "PCA")

# Insert ScType here
source("R/sctype_wrapper_sce.R")
sce <- run_sctype_sce(sce, known_tissue_type = "Immune system",
                       cluster_col = "cluster", plot = TRUE)

# Continue with downstream analysis
plotUMAP(sce, colour_by = "sctype_classification")
```

### Other Single-Cell Tools (Scanpy, etc.)

ScType core functions can work with any tool that provides:
1. Scaled expression matrix (genes × cells)
2. Cluster assignments

For Python/Scanpy users, extract the expression matrix and cluster assignments, then use `sctype_score()` directly in R:

```r
# Example: Using ScType with Scanpy-processed data
# (Assume data exported from Python as CSV files)
expression_matrix <- read.csv("scaled_expression.csv", row.names = 1)
clusters <- read.csv("clusters.csv")$cluster

# Run ScType scoring
source("R/gene_sets_prepare.R")
source("R/sctype_score_.R")
gs_list <- gene_sets_prepare("ScTypeDB_full.xlsx", "Immune system")
es.max <- sctype_score(as.matrix(expression_matrix), scaled = TRUE,
                       gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# Aggregate by cluster and assign cell types
# (Same as standard workflow)
```

---

## Performance Considerations

### Speed Optimization

- ScType is designed to be ultra-fast
- Matrix operations are vectorized
- For very large datasets (>100k cells), consider:
  - Subsampling for initial exploration
  - Running on clusters separately
  - Using scaled data (faster than unscaled)

### Memory Usage

- Expression matrices are memory-intensive
- Seurat objects can be large with multiple assays
- For large datasets:
  - Extract only necessary assay data
  - Clear unnecessary objects with `rm()`
  - Use `gc()` to free memory

---

## Contribution Guidelines

### Code Style

Based on existing code:
- Use `<-` for assignment (not `=`)
- Function names use snake_case
- Comments use `#` with space after
- Include function documentation at file top
- Use meaningful variable names

### License Compliance

- All code is GNU GPL v3.0
- Include license header in new R files:
```r
# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/master/LICENSE)
# Written by [Your Name] <email>, [Date]
```

### Version Control

- This is a Git repository
- Main branch contains stable releases
- Include clear commit messages
- Reference issues in commits when applicable

---

## Troubleshooting Guide

### Error: "object not found"

**Cause**: Missing dependencies or incorrect sourcing
**Solution**: Load all required packages and source functions in correct order

### Error: "subscript out of bounds"

**Cause**: Gene names don't match between data and markers
**Solution**: Check gene name formatting (case, symbols)

### Warning: "scRNAseqData doesn't seem to be a matrix"

**Cause**: Input is data frame or other format
**Solution**: Convert to matrix: `as.matrix(data)`

### Low scores or all "Unknown"

**Causes**:
- Wrong tissue type selected
- Unscaled data used with `scaled = TRUE`
- Markers don't match species (human vs mouse)

**Solutions**:
- Use `auto_detect_tissue_type()` to find correct tissue
- Verify scaling parameter matches data
- Use species-appropriate marker database

### No visualization appearing

**Cause**: Seurat's `DimPlot()` requires UMAP/tSNE coordinates
**Solution**: Run `RunUMAP()` or `RunTSNE()` before plotting

---

## AI Assistant Best Practices

### When Helping Users

1. **Verify Seurat version**: Check if user has v4 or v5 before suggesting code
2. **Confirm tissue type**: Ask what tissue/organ the data comes from
3. **Check data format**: Ensure data is properly scaled and formatted
4. **Suggest wrapper first**: For beginners, recommend `run_sctype()` wrapper
5. **Validate markers**: If using custom markers, verify format matches requirements

### When Debugging

1. Check that all dependencies are installed: `dplyr`, `Seurat`, `HGNChelper`, `openxlsx`
2. Verify gene symbols match between data and markers
3. Confirm matrix orientation (genes as rows)
4. Check that clusters exist and are properly assigned
5. Ensure Seurat object has required slots filled

### When Extending Functionality

1. Maintain compatibility with both Seurat v4 and v5
2. Follow existing coding conventions
3. Include parameter validation
4. Add informative error messages
5. Test with example data before suggesting to users
6. Document new features clearly

---

## Additional Resources

- **Publication**: https://doi.org/10.1038/s41467-022-28803-w
- **Web Portal**: http://sctype.app
- **Seurat Tutorials**: https://satijalab.org/seurat/
- **GitHub Issues**: https://github.com/IanevskiAleksandr/sc-type/issues
- **Contact**: aleksandr.ianevski@helsinki.fi

---

## Version Information

### Core Dependencies
- **Required**: `dplyr`, `HGNChelper`, `openxlsx`
- **For Seurat**: `Seurat`, `SeuratObject` (v4 or v5)
- **For SingleCellExperiment**: `SingleCellExperiment`, `SummarizedExperiment`, `scater` (optional, for visualization)
- **Optional**: `ggraph`, `igraph`, `tidyverse`, `data.tree` (for bubble plots), `patchwork` (for multi-panel plots)

### Compatibility
- **Seurat**: v4 and v5 fully supported
- **SingleCellExperiment**: Bioconductor 3.14+ (tested with 3.18)
- **R Version**: R ≥ 4.0.0 (tested with R 4.1+)
- **Tested Package Versions**: HGNChelper_0.8.1, SeuratObject_4.0.2, Seurat_4.0.3, dplyr_1.0.6, SingleCellExperiment_1.18+

---

*Last Updated: 2025-11-15*
*This document is intended for AI assistants to better understand and work with the ScType codebase.*
