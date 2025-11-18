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
│   ├── sctype_hierarchical_sce.R # Hierarchical annotation (SingleCellExperiment)
│   ├── sctype_visualize.R       # Marker visualization for Seurat
│   ├── sctype_visualize_sce.R   # Marker visualization for SingleCellExperiment
│   ├── sctype_uncertainty.R     # Uncertainty scoring for Seurat
│   ├── sctype_uncertainty_sce.R # Uncertainty scoring for SingleCellExperiment
│   └── sctype_pathway_enrichment.R # Pathway enrichment integration
├── python/                      # Python/scanpy integration
│   ├── sctype_python.py         # Python wrapper using rpy2
│   └── requirements.txt         # Python package dependencies
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
├── CLAUDE.md                    # AI assistant guide (this file)
├── SINGLECELLEXPERIMENT_README.md # SCE-specific documentation
├── VISUALIZATION_README.md      # Marker visualization guide
├── UNCERTAINTY_README.md        # Uncertainty scoring guide
├── PATHWAY_ENRICHMENT_README.md # Pathway enrichment integration guide
├── PYTHON_README.md             # Python/scanpy integration guide
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

### 6. Marker Visualization (`R/sctype_visualize.R` and `R/sctype_visualize_sce.R`)

**Purpose**: Comprehensive visualization of marker genes used for cell type annotation.

After running ScType annotation, these functions generate multiple visualization types to validate and explore the marker genes (positive and negative) that determined each cell type assignment.

**Main Functions**:
- `visualize_sctype_markers()` (Seurat) / `visualize_sctype_markers_sce()` (SCE): Generate all visualization types
- `quick_marker_viz()` (Seurat) / `quick_marker_viz_sce()` (SCE): Quick visualization for specific cell types

**Function**: `visualize_sctype_markers(seurat_object, annotation_col, database_file, tissue_type, assay, top_n, plot_types, save_plots, output_dir)`

**Parameters**:
- `seurat_object` / `sce_object`: Annotated object with ScType results
- `annotation_col`: Column with cell type annotations (default: "sctype_classification")
- `database_file`: Path to marker database (default: GitHub URL)
- `tissue_type`: Tissue type used for annotation (required)
- `assay` / `assay_name`: Assay for expression data (default: "RNA" / "logcounts")
- `top_n`: Number of top markers per cell type (default: 5)
- `plot_types`: Vector of "violin", "umap", "dotplot", "heatmap" (default: all)
- `save_plots`: Save to files (default: FALSE)
- `output_dir`: Directory for saved plots (default: "sctype_plots")

**Visualization Types**:

1. **Violin Plots**
   - Expression distribution of each marker across all cell types
   - Separate plots for each annotated cell type
   - Positive and negative markers labeled
   - Confirms marker expression patterns

2. **UMAP Plots**
   - Spatial visualization of marker expression on UMAP
   - Feature plots for each marker gene
   - Reveals marker specificity and co-localization
   - Requires UMAP coordinates in object

3. **Dotplot**
   - Combined view of all markers across all cell types
   - Dot size = percentage of cells expressing
   - Dot color = average expression level
   - Best for overall comparison

4. **Heatmap**
   - Hierarchically clustered heatmap of average expression
   - Rows = marker genes, Columns = cell types
   - Scaled expression (blue-white-red)
   - Reveals co-expression patterns and relationships

**Usage Example** (Seurat):
```r
source("R/sctype_visualize.R")

# Generate all visualizations
plots <- visualize_sctype_markers(
    seurat_object = seurat_obj,
    annotation_col = "sctype_classification",
    tissue_type = "Immune system",
    top_n = 5,
    save_plots = TRUE,
    output_dir = "marker_plots"
)

# Access individual plots
print(plots$dotplot)
print(plots$heatmap)
plots$violin[["CD4+ T cells"]]
plots$umap[["B cells"]]

# Quick visualization for one cell type
quick_plots <- quick_marker_viz(
    seurat_obj,
    cell_type = "NK cells",
    tissue_type = "Immune system",
    plot_type = "both"
)
```

**Usage Example** (SingleCellExperiment):
```r
source("R/sctype_visualize_sce.R")

plots <- visualize_sctype_markers_sce(
    sce_object = sce,
    annotation_col = "sctype_classification",
    tissue_type = "Brain",
    assay_name = "logcounts",
    top_n = 5,
    save_plots = TRUE
)
```

**Key Features**:
- Validates annotation quality by showing marker expression
- Identifies potential misannotations (unexpected patterns)
- Publication-ready figures
- Works with both hierarchical and standard annotations
- Saves high-resolution plots (300 dpi)

**Dependencies**:
- Required: `ggplot2`, `dplyr`, `openxlsx`, `patchwork`
- Optional: `ComplexHeatmap`, `circlize` (for heatmaps)
- Seurat: `Seurat` package
- SCE: `SingleCellExperiment`, `scater` packages

**Output Files** (when `save_plots = TRUE`):
- `violin_[CellType].png` - Violin plots for each cell type (12" × 8", 300 dpi)
- `umap_[CellType].png` - UMAP plots for each cell type (14" × 10", 300 dpi)
- `dotplot_all_markers.png` - Combined dotplot (16" × 10", 300 dpi)
- `heatmap_all_markers.png` - Hierarchical heatmap (14" × 10", 300 dpi)

**Returns**: List of ggplot objects (violin, umap, dotplot) and ComplexHeatmap object (heatmap)

**See Also**: VISUALIZATION_README.md for detailed usage guide and examples

### 7. Uncertainty Scoring (`R/sctype_uncertainty.R` and `R/sctype_uncertainty_sce.R`)

**Purpose**: Quantify confidence and uncertainty in cell type annotations.

These functions calculate comprehensive confidence metrics for ScType annotations, providing top N candidate cell types per cluster, confidence scores, and visualizations of annotation uncertainty.

**Main Functions**:
- `add_sctype_uncertainty()` (Seurat) / `add_sctype_uncertainty_sce()` (SCE): Add uncertainty metrics to object
- `visualize_sctype_uncertainty()` (Seurat) / `visualize_sctype_uncertainty_sce()` (SCE): Visualize uncertainty

**Function**: `add_sctype_uncertainty(seurat_object, known_tissue_type, database_file, assay, scaled, cluster_col, top_n, annotation_prefix)`

**Parameters**:
- `seurat_object` / `sce_object`: Object with clustering
- `known_tissue_type`: Tissue type (required)
- `database_file`: Path to marker database (default: GitHub URL)
- `assay` / `assay_name`: Assay for expression data (default: "RNA" / "logcounts")
- `scaled`: Use scaled data (default: TRUE)
- `cluster_col`: Column with clusters (default: "seurat_clusters" / "cluster")
- `top_n`: Number of top candidates to report (default: 3)
- `annotation_prefix`: Prefix for new columns (default: "sctype")

**New Metadata Columns Added**:
- `sctype_top1`, `sctype_top2`, `sctype_top3`: Top N cell type candidates
- `sctype_score1`, `sctype_score2`, `sctype_score3`: Raw ScType scores for each candidate
- `sctype_confidence`: Normalized confidence score (0-1)
- `sctype_confidence_level`: Categorical confidence (High/Medium/Low)
- `sctype_score_diff`: Score difference between top 2 candidates

**Confidence Metrics**:

1. **Raw Scores**: Direct ScType algorithm output (higher = better match)
2. **Score Difference**: `score1 - score2` (larger = more confident, clear winner)
3. **Confidence Score**: 0-1 normalized metric combining absolute and relative confidence
   - 0.0-0.4 = Low confidence
   - 0.4-0.7 = Medium confidence
   - 0.7-1.0 = High confidence
4. **Confidence Level**: Categorical classification (High/Medium/Low/Unknown)

**Visualization Types** (from `visualize_sctype_uncertainty()`):

1. **Top Candidates Plot**
   - Bar plot showing top 1-3 cell types per cluster with scores
   - Bars colored by rank (1st=red, 2nd=blue, 3rd=green)
   - Cell type labels on bars
   - Identifies ambiguous clusters

2. **UMAP Confidence Plots**
   - Plot A: Continuous confidence score (blue=low → yellow → red=high)
   - Plot B: Categorical confidence level (green/yellow/red)
   - Reveals spatial patterns in annotation quality

3. **Confidence Distribution**
   - Plot A: Bar chart of confidence levels across all cells
   - Plot B: Score differences per cluster
   - Assesses overall annotation quality

4. **Uncertainty Heatmap**
   - Compact tile plot: one tile per cluster
   - Tile color = confidence (red=low, yellow=mid, green=high)
   - Cell type and score displayed in each tile
   - At-a-glance summary

**Usage Example** (Seurat):
```r
source("R/sctype_uncertainty.R")

# Add uncertainty scores
seurat_obj <- add_sctype_uncertainty(
    seurat_obj,
    known_tissue_type = "Immune system",
    top_n = 3
)

# Visualize uncertainty
plots <- visualize_sctype_uncertainty(
    seurat_obj,
    plot_types = c("candidates", "umap", "distribution", "heatmap"),
    save_plots = TRUE,
    output_dir = "uncertainty_plots"
)

# Access plots
print(plots$candidates)      # Top candidates per cluster
print(plots$umap$score)      # UMAP by confidence score
print(plots$distribution$levels)  # Confidence level distribution
print(plots$heatmap)         # Uncertainty heatmap

# Identify low-confidence clusters
low_conf <- seurat_obj@meta.data %>%
    filter(sctype_confidence_level == "Low") %>%
    pull(seurat_clusters) %>%
    unique()
```

**Usage Example** (SingleCellExperiment):
```r
source("R/sctype_uncertainty_sce.R")

sce <- add_sctype_uncertainty_sce(
    sce,
    known_tissue_type = "Brain",
    assay_name = "logcounts",
    top_n = 3
)

plots <- visualize_sctype_uncertainty_sce(sce, save_plots = TRUE)

# View cluster-level results
metadata(sce)[["sctype_uncertainty_clusters"]]
```

**Key Features**:
- Identifies ambiguous annotations requiring manual review
- Provides top N candidates (not just best match)
- Quantifies annotation confidence objectively
- Helps prioritize high-confidence cells for downstream analysis
- Detects potentially novel or transitional cell states
- Publication-ready uncertainty visualizations

**Use Cases**:
1. **Quality Control**: Filter to high-confidence annotations for sensitive analyses
2. **Manual Curation**: Prioritize low-confidence clusters for expert review
3. **Novel Cell Types**: Identify clusters with low scores for all known types
4. **Transitional States**: Detect clusters with similar scores for multiple types
5. **Publication**: Report confidence metrics to demonstrate annotation quality

**Output Files** (when `save_plots = TRUE`):
- `top_candidates.png` - Top 3 candidates per cluster (14" × 10", 300 dpi)
- `umap_confidence.png` - UMAP plots colored by confidence (16" × 8", 300 dpi)
- `confidence_distribution.png` - Distribution plots (14" × 10", 300 dpi)
- `uncertainty_heatmap.png` - Compact heatmap summary (12" × 8", 300 dpi)

**Returns**: List of ggplot objects for all visualization types

**See Also**: UNCERTAINTY_README.md for detailed usage guide and interpretation

### 8. Pathway Enrichment Integration (`R/sctype_pathway_enrichment.R`)

**Purpose**: Integrate pathway enrichment and gene ontology analyses to validate and weight ScType annotations.

This advanced feature enhances ScType by running differential expression, performing pathway enrichment using multiple tools (EnrichR, fgsea, clusterProfiler/GO), and combining pathway support with ScType marker scores for improved confidence metrics.

**Main Functions**:
- `add_pathway_weighted_scores()`: Add pathway enrichment-weighted annotations
- `visualize_pathway_support()`: Visualize pathway support for annotations

**Function**: `add_pathway_weighted_scores(seurat_object, known_tissue_type, database_file, assay, cluster_col, enrichment_tools, top_n_genes, min_pct, logfc_threshold, annotation_prefix)`

**Parameters**:
- `seurat_object`: Seurat object with clustering
- `known_tissue_type`: Tissue type (required)
- `database_file`: Path to marker database (default: GitHub URL)
- `assay`: Assay to use (default: "RNA")
- `cluster_col`: Cluster column (default: "seurat_clusters")
- `enrichment_tools`: Vector of "enrichr", "fgsea", "go" (default: all)
- `top_n_genes`: Top DE genes for enrichment (default: 200)
- `min_pct`: Min % cells expressing (default: 0.25)
- `logfc_threshold`: LogFC threshold (default: 0.25)
- `annotation_prefix`: Column prefix (default: "sctype")

**Workflow**:

1. **Run ScType Uncertainty**: First adds standard uncertainty scores (calls `add_sctype_uncertainty()`)
2. **Differential Expression**: Identifies cluster markers using `FindAllMarkers()`
3. **Pathway Enrichment**: Runs multiple enrichment tools on marker genes:
   - **EnrichR**: Web-based enrichment against CellMarker, Azimuth, PanglaoDB, GO, KEGG, Reactome
   - **fgsea**: Fast GSEA using MSigDB Hallmark and cell type signatures (C8)
   - **clusterProfiler**: GO Biological Process and KEGG pathway enrichment
4. **Pathway-to-Cell Type Matching**: Compares enriched pathways to expected cell type functions using knowledge base
5. **Combined Scoring**: Calculates weighted confidence (60% ScType + 40% pathway)

**New Metadata Columns Added**:
- `sctype_pathway_score`: Pathway enrichment support (0-1)
- `sctype_pathway_support`: Pathway support level (High/Medium/Low/None)
- `sctype_combined_confidence`: Weighted average of ScType and pathway scores
- `sctype_top_pathways`: Top enriched pathways for cluster

**Pathway Support Scoring**:
- Matches enriched pathways to expected cell type functions
- Example: T cells should enrich for "T cell activation", "TCR signaling", "lymphocyte"
- Score = (matching pathways) / (total pathways), normalized 0-1
- **High** (≥0.7): Strong pathway support
- **Medium** (0.4-0.7): Moderate support
- **Low** (<0.4): Weak or no support

**Usage Example**:
```r
source("R/sctype_pathway_enrichment.R")

# Add pathway-weighted scores
seurat_obj <- add_pathway_weighted_scores(
    seurat_obj,
    known_tissue_type = "Immune system",
    enrichment_tools = c("enrichr", "fgsea", "go"),
    top_n_genes = 200
)

# Visualize pathway support
plots <- visualize_pathway_support(
    seurat_obj,
    save_plots = TRUE
)

print(plots$sctype_vs_pathway)  # Scatter plot comparing scores
print(plots$combined_umap)       # UMAP by combined confidence
print(plots$support_levels)      # Bar chart of support levels

# Access enrichment results
pathway_results <- attr(seurat_obj, "sctype_pathway_results")
cluster_3_enrichment <- pathway_results[["3"]]
```

**Key Features**:
- Validates marker-based annotations with functional evidence
- Identifies annotations with weak pathway support (potential issues)
- Detects missed annotations (high pathway score, low ScType)
- Combines multiple enrichment tools for robust results
- Knowledge base of cell type → pathway mappings
- Publication-ready visualizations

**Interpreting Combined Scores**:

*High ScType + High Pathway*: Excellent annotation (both marker and function support)
*High ScType + Low Pathway*: Review - markers present but unexpected pathways (may indicate doublets or low-quality cells)
*Low ScType + High Pathway*: Missed annotation - pathway clearly indicates cell type but markers not in database
*Low ScType + Low Pathway*: Unknown cell type, low quality, or novel population

**Required Packages** (install at least one enrichment tool):
- **EnrichR**: `install.packages("enrichR")` (web-based, easiest)
- **fgsea**: `BiocManager::install(c("fgsea", "msigdbr"))` (offline, fast)
- **GO/KEGG**: `BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))` (comprehensive)

**Visualization Types**:
1. **ScType vs Pathway Scatter**: Compares ScType confidence with pathway support scores
2. **Combined Confidence UMAP**: Spatial distribution of weighted confidence
3. **Pathway Support Levels**: Bar chart showing High/Medium/Low/None distribution

**Output Files** (when `save_plots = TRUE`):
- `sctype_vs_pathway.png` - Confidence comparison scatter plot (12" × 8", 300 dpi)
- `combined_confidence_umap.png` - UMAP colored by combined score (10" × 8", 300 dpi)
- `pathway_support_levels.png` - Support level distribution (10" × 6", 300 dpi)

**Limitations**:
- Requires differential expression (computationally intensive)
- Pathway databases are primarily human/mouse
- Knowledge base coverage varies by cell type
- EnrichR requires internet connection
- Enrichment quality depends on sequencing depth

**See Also**: PATHWAY_ENRICHMENT_README.md for detailed usage, interpretation guide, and troubleshooting

### 9. Python/Scanpy Integration (`python/sctype_python.py`)

**Purpose**: Python wrapper for ScType enabling seamless integration with scanpy workflows and AnnData objects.

This module provides a bridge between Python/scanpy and R-based ScType functions using `rpy2`, allowing users to run ScType cell type annotation directly on scanpy `AnnData` objects without leaving Python.

**Main Classes and Functions**:
- `ScType` class: Main interface with methods for annotation, hierarchical annotation, and uncertainty scoring
- `run_sctype()`: Convenience function for quick annotation

**Class**: `ScType(github_repo="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master")`

**Methods**:
- `annotate(adata, tissue_type, ...)`: Basic cell type annotation
- `annotate_hierarchical(adata, tissue_type, ...)`: Hierarchical annotation (broad + fine)
- `add_uncertainty(adata, tissue_type, ...)`: Add uncertainty scores
- `visualize_markers(adata, tissue_type, ...)`: Generate marker visualization plots

**Function**: `run_sctype(adata, tissue_type, database_file, layer, scaled, cluster_key, annotation_key, plot)`

**Parameters**:
- `adata` (AnnData): Annotated data object with clustering (scanpy format)
- `tissue_type` (str): Tissue type (e.g., "Immune system", "Brain", "Liver")
- `database_file` (str, optional): Path or URL to custom marker database
- `layer` (str, optional): Layer to use for expression data (None = .X)
- `scaled` (bool): Whether data is scaled (default: True)
- `cluster_key` (str): Key in adata.obs with cluster assignments (default: 'leiden')
- `annotation_key` (str): Key to store annotations in adata.obs (default: 'sctype_classification')
- `plot` (bool): Generate UMAP plot with annotations (default: False)

**Returns**: Modified AnnData object with annotations in `.obs[annotation_key]`

**Usage Example**:
```python
import scanpy as sc
import sys
sys.path.append('path/to/sc-type/python')
from sctype_python import run_sctype

# Load your data
adata = sc.read_h5ad("pbmc3k.h5ad")

# Standard scanpy preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata, resolution=0.8)
sc.tl.umap(adata)

# Run ScType annotation
adata = run_sctype(
    adata,
    tissue_type="Immune system",
    scaled=True,
    cluster_key='leiden',
    annotation_key='sctype_classification',
    plot=True
)

# View results
print(adata.obs['sctype_classification'].value_counts())
sc.pl.umap(adata, color=['leiden', 'sctype_classification'])
```

**Class Interface Example**:
```python
from sctype_python import ScType

# Initialize ScType
sctype = ScType()

# Run annotation
adata = sctype.annotate(
    adata,
    tissue_type="Brain",
    cluster_key='louvain',
    annotation_key='cell_type',
    plot=True
)

# Run hierarchical annotation
adata = sctype.annotate_hierarchical(
    adata,
    tissue_type="Immune system",
    broad_key='sctype_broad',
    fine_key='sctype_fine'
)

# Add uncertainty scores
adata = sctype.add_uncertainty(
    adata,
    tissue_type="Immune system",
    top_n=3
)

# Visualize marker genes
plots = sctype.visualize_markers(
    adata,
    tissue_type="Immune system",
    annotation_col='cell_type',
    top_n=5,
    plot_types=['violin', 'umap', 'dotplot', 'heatmap'],
    save_plots=True
)

# Access results
cell_types = adata.obs['cell_type']
confidence = adata.obs['sctype_confidence']
```

**Key Features**:
- **Native scanpy Integration**: Works directly with AnnData objects
- **Automatic Conversion**: Seamlessly converts between Python and R data formats
- **rpy2 Bridge**: Leverages R ScType functions without manual data export
- **Full Functionality**: Supports basic annotation, hierarchical annotation, and uncertainty scoring
- **Comprehensive Visualization**: Violin plots, UMAP feature plots, dotplots, and heatmaps of marker expression
- **Publication-Ready Figures**: High-resolution plots (300 dpi) with customizable sizes

**Data Conversion**:
The wrapper automatically handles conversion between Python and R formats:
- **AnnData → R matrix**: Extracts expression matrix (genes × cells)
- **Sparse → Dense**: Converts sparse matrices to dense for R compatibility
- **Python lists → R vectors**: Uses rpy2 converters (pandas2ri, numpy2ri)
- **R results → Python**: Converts annotations back to pandas categorical

**Implementation Details**:
```python
# Core conversion pattern
def _adata_to_matrix(self, adata, layer=None):
    # Extract expression data
    if layer is None:
        mat = adata.X
    else:
        mat = adata.layers[layer]

    # Convert sparse to dense if needed
    if scipy.sparse.issparse(mat):
        mat = mat.toarray()

    # Transpose to genes × cells (R format)
    mat = mat.T

    # Get gene and cell names
    genes = adata.var_names.tolist()
    cells = adata.obs_names.tolist()

    return mat, genes, cells

# Run R ScType scoring
ro.r('es.max <- sctype_score(scRNAseqData, scaled, gs_list$gs_positive, gs_list$gs_negative)')

# Convert results back to Python
annotations = np.array([cluster_assignments.get(cl, "Unknown") for cl in clusters])
adata.obs[annotation_key] = pd.Categorical(annotations)
```

**Prerequisites**:
1. **R Installation**: R must be installed and accessible
2. **R Packages**: `dplyr`, `HGNChelper`, `openxlsx`
3. **Python Packages**: `numpy`, `pandas`, `scanpy`, `anndata`, `rpy2`, `matplotlib`, `openpyxl`
4. **Optional**: `seaborn` (for enhanced visualization)

**Installation**:
```bash
# Install Python dependencies
cd sc-type/python
pip install -r requirements.txt

# Or install individually
pip install numpy pandas scanpy anndata rpy2 matplotlib openpyxl seaborn
```

**Troubleshooting**:

**Error: Cannot find R**
```python
# Set R_HOME environment variable before importing rpy2
import os
os.environ['R_HOME'] = '/usr/lib/R'  # Linux
# os.environ['R_HOME'] = '/Library/Frameworks/R.framework/Resources'  # macOS
```

**Memory Issues**
```python
# Downsample or process in batches
adata_subset = sc.pp.subsample(adata, fraction=0.1, copy=True)
adata_subset = run_sctype(adata_subset, tissue_type="Immune system")

# Or increase R memory limit
import rpy2.robjects as ro
ro.r('memory.limit(size=16000)')  # 16 GB
```

**All Cells Annotated as "Unknown"**
- Check tissue type spelling (case-sensitive)
- Verify data is scaled if `scaled=True`
- Ensure cluster_key exists in adata.obs
- Check that expression data format is correct (genes as rows after conversion)

**Limitations**:
1. **R Dependency**: Requires R installation and rpy2
2. **Advanced Features**: Pathway enrichment integration not yet available in Python wrapper
3. **Performance**: Slightly slower than native R due to conversion overhead
4. **Hierarchical Annotation**: Currently falls back to standard annotation (R version has full support)

**Comparison: Python vs R**:

| Feature | Python (scanpy) | R (Seurat/SCE) |
|---------|----------------|----------------|
| **Basic annotation** | ✓ Supported | ✓ Supported |
| **Hierarchical annotation** | ⚠ Fallback to standard | ✓ Full support |
| **Uncertainty scoring** | ⚠ Fallback to standard | ✓ Full support |
| **Marker visualization** | ✓ Supported (4 plot types) | ✓ Full support |
| **Pathway enrichment** | ✗ Not yet implemented | ✓ Full support |
| **Data object** | AnnData | Seurat/SCE |
| **Plotting** | matplotlib/scanpy | ggplot2/Seurat |

**Marker Visualization**:

The Python wrapper includes comprehensive marker visualization capabilities:

```python
# After annotation
plots = sctype.visualize_markers(
    adata,
    tissue_type="Immune system",
    annotation_col="sctype_classification",
    top_n=5,
    plot_types=['violin', 'umap', 'dotplot', 'heatmap'],
    save_plots=True,
    output_dir="marker_plots"
)

# Access individual plots
plots['dotplot'].show()
plots['heatmap'].show()
plots['violin']['CD4+ T cells'].show()
```

**Visualization Types**:
1. **Violin Plots**: Expression distribution across cell types for each marker
2. **UMAP Feature Plots**: Spatial marker expression patterns
3. **Dotplot**: Combined view of all markers (uses scanpy.pl.dotplot)
4. **Heatmap**: Mean expression heatmap across cell types

**Output**: Dictionary with plot objects and optionally saved high-resolution PNG files (300 dpi)

**See Also**: PYTHON_README.md for complete installation guide, detailed examples, and troubleshooting

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

### Python/Scanpy Pipeline Integration

ScType now provides native Python integration through the `python/sctype_python.py` wrapper, enabling seamless use with scanpy workflows:

```python
import scanpy as sc
import sys
sys.path.append('path/to/sc-type/python')
from sctype_python import run_sctype

# Standard scanpy workflow
adata = sc.read_h5ad("data.h5ad")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata, resolution=0.8)
sc.tl.umap(adata)

# Insert ScType here
adata = run_sctype(adata, tissue_type="Immune system", cluster_key='leiden', plot=True)

# Continue with downstream analysis
sc.pl.umap(adata, color='sctype_classification')
```

**Alternative: Manual R Integration**

For other Python tools or custom workflows, ScType core functions can work with any tool that provides:
1. Scaled expression matrix (genes × cells)
2. Cluster assignments

Extract the expression matrix and cluster assignments, then use `sctype_score()` directly in R:

```r
# Example: Using ScType with Python-processed data
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

**See Also**: PYTHON_README.md for complete Python/scanpy integration guide

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

### Core R Dependencies
- **Required**: `dplyr`, `HGNChelper`, `openxlsx`
- **For Seurat**: `Seurat`, `SeuratObject` (v4 or v5)
- **For SingleCellExperiment**: `SingleCellExperiment`, `SummarizedExperiment`, `scater` (for SCE operations and visualization)
- **For Visualization**: `ggplot2`, `patchwork` (required); `ComplexHeatmap`, `circlize` (optional, for heatmaps)
- **For Pathway Enrichment**: `enrichR` (optional), `fgsea`, `msigdbr` (optional), `clusterProfiler`, `org.Hs.eg.db` (optional)
- **Other Optional**: `ggraph`, `igraph`, `tidyverse`, `data.tree` (for bubble plots and networks)

### Python Dependencies
- **Required**: `numpy`, `pandas`, `scanpy`, `anndata`, `rpy2`, `matplotlib`, `openpyxl`
- **Optional**: `seaborn` (enhanced visualization), `scipy`, `scikit-learn`
- **R Installation**: R ≥ 4.0.0 with core R packages (required for rpy2 bridge)

### Compatibility
- **R Version**: R ≥ 4.0.0 (tested with R 4.1+)
- **Python Version**: Python ≥ 3.7 (tested with Python 3.8+)
- **Seurat**: v4 and v5 fully supported
- **SingleCellExperiment**: Bioconductor 3.14+ (tested with 3.18)
- **scanpy**: Version 1.8+ (tested with 1.9+)
- **rpy2**: Version 3.4+ (tested with 3.5+)
- **Tested Package Versions**:
  - R: HGNChelper_0.8.1, SeuratObject_4.0.2, Seurat_4.0.3, dplyr_1.0.6, SingleCellExperiment_1.18+
  - Python: scanpy_1.9.1, anndata_0.8.0, rpy2_3.5.1, numpy_1.21+, pandas_1.3+

---

*Last Updated: 2025-11-15*
*This document is intended for AI assistants to better understand and work with the ScType codebase.*
