# ScType for Python/Scanpy

Python wrapper for ScType cell type annotation, enabling seamless integration with scanpy workflows.

## Overview

This Python module provides a bridge between Python/scanpy and R-based ScType functions using `rpy2`. It allows you to run ScType cell type annotation directly on scanpy `AnnData` objects without leaving Python.

**Key Features**:
- Native scanpy `AnnData` support
- Automatic conversion between Python and R formats
- Compatible with standard scanpy workflows
- Access to all ScType databases and marker sets
- UMAP visualization integration

## Installation

### Prerequisites

1. **R installation** (required for rpy2)
   ```bash
   # Check if R is installed
   R --version

   # If not installed:
   # Ubuntu/Debian
   sudo apt-get install r-base r-base-dev

   # macOS
   brew install r

   # Windows: Download from https://cran.r-project.org/
   ```

2. **R packages** (required for ScType)
   ```r
   # In R console:
   install.packages(c("dplyr", "HGNChelper", "openxlsx"))
   ```

### Python Installation

```bash
# Install ScType Python wrapper and dependencies
cd sc-type/python
pip install -r requirements.txt

# Or install packages individually:
pip install numpy pandas scanpy anndata rpy2 matplotlib
```

### Verify Installation

```python
import rpy2.robjects as ro
import scanpy as sc
print("R version:", ro.r('R.version.string')[0])
print("scanpy version:", sc.__version__)
```

## Quick Start

### Complete Workflow Example

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

# Visualize
sc.pl.umap(adata, color=['leiden', 'sctype_classification'])

# Save annotated data
adata.write("pbmc3k_sctype_annotated.h5ad")
```

### Using the Class Interface

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

# Access results
cell_types = adata.obs['cell_type']
print(cell_types.value_counts())
```

## Detailed Usage

### Function: `run_sctype()`

Main entry point for ScType annotation.

**Parameters**:
- `adata` (AnnData): Annotated data object with clustering
- `tissue_type` (str): Tissue type (e.g., "Immune system", "Brain", "Liver")
- `database_file` (str, optional): Path to custom marker database
- `layer` (str, optional): Layer to use. If None, uses .X
- `scaled` (bool): Whether data is scaled (default: True)
- `cluster_key` (str): Key in adata.obs with clusters (default: 'leiden')
- `annotation_key` (str): Key to store annotations (default: 'sctype_classification')
- `plot` (bool): Generate UMAP plot (default: False)

**Returns**:
- `adata` (AnnData): Modified AnnData with annotations in .obs[annotation_key]

**Example**:
```python
adata = run_sctype(
    adata,
    tissue_type="Immune system",
    scaled=True,
    cluster_key='leiden',
    plot=True
)
```

### Available Tissue Types

Same as R version:
- "Immune system"
- "Brain"
- "Liver"
- "Pancreas"
- "Kidney"
- "Lung"
- "Heart"
- "Intestine"
- "Muscle"
- "Eye"
- "Adrenal"
- "Placenta"
- "Spleen"
- "Stomach"
- "Thymus"

### Custom Marker Databases

```python
# Use custom database
adata = run_sctype(
    adata,
    tissue_type="My Custom Tissue",
    database_file="/path/to/custom_markers.xlsx"
)

# Use enhanced database
adata = run_sctype(
    adata,
    tissue_type="Immune system",
    database_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_enhanced.xlsx"
)
```

### Working with Different Clustering Methods

```python
# Leiden clustering
sc.tl.leiden(adata, resolution=0.8)
adata = run_sctype(adata, tissue_type="Brain", cluster_key='leiden')

# Louvain clustering
sc.tl.louvain(adata, resolution=0.8)
adata = run_sctype(adata, tissue_type="Brain", cluster_key='louvain')

# Custom clusters
adata.obs['my_clusters'] = custom_cluster_labels
adata = run_sctype(adata, tissue_type="Brain", cluster_key='my_clusters')
```

### Using Different Data Layers

```python
# Use scaled data (default)
adata = run_sctype(adata, tissue_type="Immune system", scaled=True)

# Use raw counts
adata = run_sctype(adata, tissue_type="Immune system", layer=None, scaled=False)

# Use specific layer
adata = run_sctype(adata, tissue_type="Immune system", layer='log1p_transformed')
```

## Integration with Scanpy Workflows

### After Cell Type Annotation

```python
# Differential expression by cell type
sc.tl.rank_genes_groups(adata, 'sctype_classification', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

# Cell type-specific analyses
cd4_tcells = adata[adata.obs['sctype_classification'] == 'CD4+ T cells']
sc.tl.paga(cd4_tcells)
sc.pl.paga(cd4_tcells)

# Export for Seurat (if needed)
adata.write("data_with_sctype.h5ad")
# Can be read in R with: adata <- ReadH5AD("data_with_sctype.h5ad")
```

### Visualization

```python
import matplotlib.pyplot as plt

# UMAP with annotations
sc.pl.umap(adata, color='sctype_classification', legend_loc='on data')

# Multiple visualizations
sc.pl.umap(adata, color=['leiden', 'sctype_classification', 'n_genes', 'n_counts'],
          ncols=2)

# Dotplot of marker genes
sc.pl.dotplot(adata, marker_genes, groupby='sctype_classification')

# Heatmap
sc.pl.heatmap(adata, marker_genes, groupby='sctype_classification')

# Save figures
plt.savefig('sctype_umap.png', dpi=300, bbox_inches='tight')
```

## Complete Example: PBMC Dataset

```python
import scanpy as sc
import numpy as np
import pandas as pd
from sctype_python import run_sctype

# Download example data (PBMC 3k)
adata = sc.datasets.pbmc3k()

# Preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Normalize and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]

# Scale
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

# PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)

# Neighbors and clustering
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata, resolution=0.8)

# UMAP
sc.tl.umap(adata)

# ScType annotation
adata = run_sctype(
    adata,
    tissue_type="Immune system",
    cluster_key='leiden',
    annotation_key='cell_type',
    plot=False
)

# Visualizations
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(adata, color='leiden', ax=axes[0], show=False, title='Leiden Clusters')
sc.pl.umap(adata, color='cell_type', ax=axes[1], show=False, title='ScType Annotations')
plt.tight_layout()
plt.savefig('pbmc_sctype_comparison.png', dpi=300)

# Summary statistics
print("\nCell type distribution:")
print(adata.obs['cell_type'].value_counts())

# Cell type per cluster
ct_per_cluster = pd.crosstab(adata.obs['leiden'], adata.obs['cell_type'])
print("\nCell types per cluster:")
print(ct_per_cluster)

# Marker genes for each cell type
marker_genes = {
    'CD4+ T cells': ['CD4', 'IL7R', 'CCR7'],
    'CD8+ T cells': ['CD8A', 'CD8B'],
    'B cells': ['CD19', 'MS4A1', 'CD79A'],
    'NK cells': ['GNLY', 'NKG7', 'KLRB1'],
    'Monocytes': ['CD14', 'LYZ', 'CD68']
}

# Dotplot
all_markers = [g for genes in marker_genes.values() for g in genes]
sc.pl.dotplot(adata, all_markers, groupby='cell_type', figsize=(10, 6))
plt.savefig('sctype_marker_dotplot.png', dpi=300, bbox_inches='tight')

# Save results
adata.write("pbmc3k_sctype_annotated.h5ad")
```

## Troubleshooting

### rpy2 Installation Issues

**Error**: `Cannot find R`

**Solution**:
```bash
# Set R_HOME environment variable
export R_HOME=/usr/lib/R  # Linux
export R_HOME=/Library/Frameworks/R.framework/Resources  # macOS

# Or in Python before importing rpy2:
import os
os.environ['R_HOME'] = '/usr/lib/R'
```

**Error**: `rpy2 version mismatch`

**Solution**:
```bash
pip install --upgrade rpy2
```

### Memory Issues

**Error**: `MemoryError` or R crashes

**Solution**:
```python
# Process data in batches or downsample
adata_subset = sc.pp.subsample(adata, fraction=0.1, copy=True)
adata_subset = run_sctype(adata_subset, tissue_type="Immune system")

# Or increase R memory limit
import rpy2.robjects as ro
ro.r('memory.limit(size=16000)')  # 16 GB
```

### Cluster Key Not Found

**Error**: `Cluster key 'leiden' not found in adata.obs`

**Solution**:
```python
# Run clustering first
sc.tl.leiden(adata, resolution=0.8)

# Or specify correct key
print(adata.obs.columns)  # See available columns
adata = run_sctype(adata, tissue_type="Brain", cluster_key='your_cluster_column')
```

### All Cells Annotated as "Unknown"

**Possible causes**:
1. Wrong tissue type
2. Unscaled data when `scaled=True`
3. Wrong cluster key

**Solutions**:
```python
# Check tissue type spelling
adata = run_sctype(adata, tissue_type="Immune system")  # Correct capitalization

# Check if data is scaled
if 'X_scaled' not in adata.layers:
    sc.pp.scale(adata)

# Verify data format
print("Data shape:", adata.X.shape)
print("Scaled:", adata.X.mean(), adata.X.std())  # Should be ~0 and ~1
```

## Comparison: Python vs R

### Python (scanpy)
```python
import scanpy as sc
from sctype_python import run_sctype

adata = sc.read_h5ad("data.h5ad")
sc.tl.leiden(adata)
adata = run_sctype(adata, tissue_type="Immune system")
sc.pl.umap(adata, color='sctype_classification')
```

### R (Seurat)
```r
library(Seurat)
source("R/sctype_wrapper.R")

seurat_obj <- readRDS("data.rds")
seurat_obj <- FindClusters(seurat_obj)
seurat_obj <- run_sctype(seurat_obj, known_tissue_type = "Immune system")
DimPlot(seurat_obj, group.by = "sctype_classification")
```

## Limitations

1. **R dependency**: Requires R installation and rpy2
2. **Advanced features**: Hierarchical annotation, uncertainty scoring, pathway enrichment currently only in R
3. **Performance**: Slightly slower than native R due to conversion overhead
4. **Visualization**: Limited to scanpy plotting (use R for ScType-specific visualizations)

## Advanced: Direct R Access

For advanced users who need features not yet wrapped:

```python
import rpy2.robjects as ro

# Source R functions directly
ro.r('source("R/sctype_uncertainty.R")')

# Run R commands
ro.r('''
# Your R code here
''')

# Or convert AnnData to Seurat in R
# (requires reticulate package in R)
```

## Citation

If you use ScType in your research, please cite:

Ianevski, A., Giri, A.K. & Aittokallio, T. Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data. Nat Commun 13, 1246 (2022). https://doi.org/10.1038/s41467-022-28803-w

## Support

- **GitHub Issues**: https://github.com/IanevskiAleksandr/sc-type/issues
- **Documentation**: See CLAUDE.md and other README files
- **Contact**: aleksandr.ianevski@helsinki.fi

---

*Last Updated: 2025-11-15*
