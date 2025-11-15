# ScType Pathway Enrichment Integration

Integrate pathway enrichment and gene ontology analyses to validate and weight ScType cell type annotations.

## Overview

This feature enhances ScType annotations by:
- Running differential expression analysis per cluster
- Performing pathway enrichment using multiple tools (EnrichR, fgsea, GO/KEGG)
- Matching enriched pathways to expected cell type functions
- Calculating pathway support scores
- Combining ScType and pathway scores for improved confidence metrics

**Key Benefit**: Annotations supported by both marker expression AND functional pathways are more trustworthy.

## Installation

### Required Packages

```r
# Core packages
install.packages(c("Seurat", "dplyr", "ggplot2"))

# Enrichment tools (install at least one)
install.packages("BiocManager")

# EnrichR
install.packages("enrichR")

# fgsea
BiocManager::install("fgsea")
BiocManager::install("msigdbr")

# GO/KEGG enrichment
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")  # Human
# BiocManager::install("org.Mm.eg.db")  # Mouse
```

## Quick Start

### For Seurat Objects

```r
library(Seurat)

# Load pathway enrichment functions
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_pathway_enrichment.R")

# Add pathway-weighted scores
seurat_obj <- add_pathway_weighted_scores(
    seurat_obj,
    known_tissue_type = "Immune system",
    enrichment_tools = c("enrichr", "fgsea", "go"),  # Use all available
    top_n_genes = 200  # Top 200 DE genes per cluster
)

# Visualize pathway support
plots <- visualize_pathway_support(
    seurat_obj,
    save_plots = TRUE,
    output_dir = "pathway_support_plots"
)

# View results
print(plots$sctype_vs_pathway)
print(plots$combined_umap)
print(plots$support_levels)
```

## New Metadata Columns

After running `add_pathway_weighted_scores()`, these columns are added:

| Column | Description | Type | Range |
|--------|-------------|------|-------|
| `sctype_pathway_score` | Pathway enrichment support score | Numeric | 0 to 1 |
| `sctype_pathway_support` | Level of pathway support | Factor | High/Medium/Low/None |
| `sctype_combined_confidence` | ScType + pathway combined score | Numeric | 0 to 1 |
| `sctype_top_pathways` | Top enriched pathways for cluster | Character | - |

**Plus all columns from** `add_sctype_uncertainty()`:
- `sctype_top1`, `sctype_top2`, `sctype_top3` (top candidates)
- `sctype_score1`, `sctype_score2`, `sctype_score3` (ScType scores)
- `sctype_confidence` (ScType confidence)
- `sctype_confidence_level` (High/Medium/Low)

## How It Works

### Step 1: Differential Expression

For each cluster, identifies marker genes using `FindAllMarkers()`:
- Only positive markers (upregulated in cluster)
- Minimum 25% of cells expressing (customizable)
- Log fold-change threshold 0.25 (customizable)
- Returns top N genes per cluster (default: 200)

### Step 2: Pathway Enrichment

Runs multiple enrichment tools on cluster marker genes:

**EnrichR**:
- Web-based enrichment against multiple databases
- Includes: CellMarker, Azimuth, PanglaoDB, GO, KEGG, Reactome
- Returns terms with adjusted p-value < 0.05

**fgsea** (Fast Gene Set Enrichment Analysis):
- Ranked-based enrichment using MSigDB gene sets
- Uses Hallmark pathways and cell type signatures (C8)
- Efficient for large gene lists

**clusterProfiler** (GO/KEGG):
- Gene Ontology Biological Process enrichment
- KEGG pathway analysis
- Converts gene symbols to Entrez IDs
- Returns enriched terms with q-value < 0.05

### Step 3: Pathway-to-Cell Type Matching

Uses knowledge base mapping cell types to expected pathways:

**Examples**:
- **T cells** → "T cell", "lymphocyte", "TCR", "CD3", "adaptive immune"
- **CD4+ T cells** → "CD4", "T helper", "Th1", "Th2"
- **Neurons** → "neuron", "synapse", "axon", "neurotransmitter"
- **Macrophages** → "macrophage", "phagocytosis", "inflammation"

**Scoring**:
- Counts how many enriched pathways match expected cell type functions
- Score = (matching pathways) / (total pathways)
- Normalized to 0-1 range

### Step 4: Combined Confidence

Calculates weighted average of ScType and pathway scores:

```
Combined Confidence = 0.6 × ScType Confidence + 0.4 × Pathway Score
```

**Rationale**:
- ScType marker-based scoring is primary (60% weight)
- Pathway support provides validation (40% weight)
- Adjustable weights if needed

## Pathway Support Levels

Categorized based on pathway score:

- **High** (≥0.7): Strong pathway support - enriched pathways match cell type well
- **Medium** (0.4-0.7): Moderate support - some pathway evidence
- **Low** (<0.4): Weak support - enriched pathways don't match expectations
- **None**: No enrichment or no pathway analysis

## Interpreting Results

### Scenario 1: High ScType + High Pathway Support

```
Cluster 0:
- sctype_top1: "CD4+ T cells"
- sctype_confidence: 0.85
- sctype_pathway_score: 0.80
- sctype_pathway_support: "High"
- sctype_combined_confidence: 0.83
- sctype_top_pathways: "T cell activation; CD4+ T helper; TCR signaling"
```

**Interpretation**: **Excellent annotation**. Both marker expression and functional pathways strongly support CD4+ T cells. High confidence for downstream use.

### Scenario 2: High ScType + Low Pathway Support

```
Cluster 3:
- sctype_top1: "Monocytes"
- sctype_confidence: 0.75
- sctype_pathway_score: 0.25
- sctype_pathway_support: "Low"
- sctype_combined_confidence: 0.55
- sctype_top_pathways: "ribosomal proteins; oxidative phosphorylation"
```

**Interpretation**: **Questionable annotation**. Markers suggest monocytes, but enriched pathways are generic metabolic processes, not myeloid-specific. **Action**: Review markers, check for low-quality cells or doublets.

### Scenario 3: Low ScType + High Pathway Support

```
Cluster 5:
- sctype_top1: "Unknown"
- sctype_confidence: 0.35
- sctype_pathway_score: 0.70
- sctype_pathway_support: "High"
- sctype_combined_confidence: 0.49
- sctype_top_pathways: "B cell activation; immunoglobulin; BCR signaling"
```

**Interpretation**: **Missed annotation**. ScType didn't find good marker match (possibly missing from database), but pathways clearly indicate B cells. **Action**: Manual annotation as B cells, update marker database.

### Scenario 4: Low ScType + Low Pathway Support

```
Cluster 7:
- sctype_top1: "Unknown"
- sctype_confidence: 0.20
- sctype_pathway_score: 0.15
- sctype_pathway_support: "Low"
- sctype_combined_confidence: 0.18
- sctype_top_pathways: "No enrichment"
```

**Interpretation**: **Unidentifiable cluster**. No marker support, no pathway support. **Possible causes**: Novel cell type, low-quality cells, doublets, contamination. **Action**: QC check, differential expression analysis, literature search.

## Visualization

### 1. ScType vs Pathway Confidence Plot

Scatter plot comparing ScType confidence and pathway support scores.

**Features**:
- X-axis: ScType confidence
- Y-axis: Pathway support score
- Color: Cell type
- Diagonal line: Equal confidence

**Interpretation**:
- **Above diagonal**: Pathway support > ScType (pathway validates weak ScType)
- **On diagonal**: Agreement between methods
- **Below diagonal**: ScType strong but weak pathway support (investigate)

```r
plots <- visualize_pathway_support(seurat_obj)
print(plots$sctype_vs_pathway)
```

### 2. Combined Confidence UMAP

UMAP colored by combined confidence score.

**Features**:
- Color scale: blue (low) → yellow (mid) → red (high)
- Shows spatial distribution of confidence
- Identifies high/low confidence regions

```r
print(plots$combined_umap)
```

### 3. Pathway Support Levels

Bar chart showing distribution of pathway support categories.

**Interpretation**:
- Ideal: Most cells in "High" category
- If many "Low"/"None": Dataset may lack characteristic pathways or need deeper sequencing

```r
print(plots$support_levels)
```

## Complete Workflow Example

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

# Add pathway-weighted ScType annotations
source("R/sctype_pathway_enrichment.R")

seurat_obj <- add_pathway_weighted_scores(
    seurat_obj,
    known_tissue_type = "Immune system",
    enrichment_tools = c("enrichr", "fgsea", "go"),
    top_n_genes = 200,
    min_pct = 0.25,
    logfc_threshold = 0.25
)

# Visualize pathway support
plots <- visualize_pathway_support(
    seurat_obj,
    save_plots = TRUE,
    output_dir = "immune_pathway_support"
)

# Examine results
## Overall agreement
print(plots$sctype_vs_pathway)

## Combined confidence spatial pattern
print(plots$combined_umap)

## Pathway support distribution
print(plots$support_levels)

# Filter to high combined confidence
high_conf_cells <- seurat_obj@meta.data %>%
    filter(sctype_combined_confidence >= 0.7)

message(paste0(nrow(high_conf_cells), " cells (",
               round(100*nrow(high_conf_cells)/ncol(seurat_obj), 1),
               "%) have high combined confidence"))

# Identify clusters with pathway-ScType disagreement
disagreement <- seurat_obj@meta.data %>%
    mutate(diff = abs(sctype_confidence - sctype_pathway_score)) %>%
    group_by(seurat_clusters) %>%
    summarise(
        mean_diff = mean(diff),
        celltype = first(sctype_top1),
        sctype_conf = first(sctype_confidence),
        pathway_score = first(sctype_pathway_score),
        top_pathways = first(sctype_top_pathways)
    ) %>%
    arrange(desc(mean_diff))

print("Clusters with largest ScType-Pathway disagreement:")
print(head(disagreement))

# View pathway enrichment results
pathway_results <- attr(seurat_obj, "sctype_pathway_results")

# Examine specific cluster
cluster_3_pathways <- pathway_results[["3"]]

# EnrichR results
if (!is.null(cluster_3_pathways$enrichr)) {
    print("EnrichR results for cluster 3:")
    print(cluster_3_pathways$enrichr$CellMarker_Augmented_2021)
}

# fgsea results
if (!is.null(cluster_3_pathways$fgsea)) {
    print("fgsea results for cluster 3:")
    print(head(cluster_3_pathways$fgsea))
}

# GO results
if (!is.null(cluster_3_pathways$go)) {
    print("GO results for cluster 3:")
    print(cluster_3_pathways$go)
}
```

## Advanced Usage

### Custom Enrichment Tools

```r
# Use only EnrichR (fastest, web-based)
seurat_obj <- add_pathway_weighted_scores(
    seurat_obj,
    known_tissue_type = "Brain",
    enrichment_tools = "enrichr"
)

# Use only fgsea (no internet required)
seurat_obj <- add_pathway_weighted_scores(
    seurat_obj,
    known_tissue_type = "Liver",
    enrichment_tools = "fgsea"
)

# Use GO and KEGG only
seurat_obj <- add_pathway_weighted_scores(
    seurat_obj,
    known_tissue_type = "Heart",
    enrichment_tools = "go"
)
```

### Custom DE Parameters

```r
# More stringent DE
seurat_obj <- add_pathway_weighted_scores(
    seurat_obj,
    known_tissue_type = "Immune system",
    top_n_genes = 100,      # Fewer genes
    min_pct = 0.5,          # More stringent expression
    logfc_threshold = 0.5   # Higher fold-change
)

# More permissive DE
seurat_obj <- add_pathway_weighted_scores(
    seurat_obj,
    known_tissue_type = "Brain",
    top_n_genes = 500,      # More genes
    min_pct = 0.1,          # Lower expression threshold
    logfc_threshold = 0.1   # Lower fold-change
)
```

### Adjusting Combined Confidence Weights

Currently hardcoded as 60% ScType, 40% pathway. To modify:

```r
# After running add_pathway_weighted_scores()
# Recalculate with custom weights (e.g., 70% ScType, 30% pathway)
seurat_obj@meta.data$sctype_combined_confidence_custom <-
    0.7 * seurat_obj@meta.data$sctype_confidence +
    0.3 * seurat_obj@meta.data$sctype_pathway_score
```

## Troubleshooting

### Error: "No enrichment tools available"

**Cause**: None of the enrichment packages are installed

**Solution**:
```r
# Install at least one:
install.packages("enrichR")  # Easiest
# OR
BiocManager::install("fgsea")
BiocManager::install("msigdbr")
# OR
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
```

### Warning: "No enrichment for cluster X"

**Cause**: Cluster has too few marker genes or genes don't match pathway databases

**Result**: Pathway score = 0, pathway support = "None"

**Actions**:
- Check if cluster is real (not low-quality cells)
- Lower `min_pct` and `logfc_threshold` for more genes
- Increase `top_n_genes`

### All pathway scores are low

**Possible causes**:
1. **Wrong species**: Using human databases on mouse data
2. **Gene naming**: Gene symbols don't match database
3. **Non-specific clusters**: Clusters lack characteristic pathways
4. **Shallow sequencing**: Not enough genes detected

**Solutions**:
```r
# Check gene overlap with pathway databases
marker_genes <- FindAllMarkers(seurat_obj, only.pos = TRUE)
unique_genes <- unique(marker_genes$gene)

# For EnrichR check
library(enrichR)
dbs <- listEnrichrDbs()
# Genes should be in standard symbol format

# For GO check
library(org.Hs.eg.db)
gene_ids <- bitr(unique_genes, fromType = "SYMBOL",
                toType = "ENTREZID", OrgDb = org.Hs.eg.db)
message(paste0(nrow(gene_ids), " out of ", length(unique_genes),
               " genes have Entrez IDs"))
```

### EnrichR API errors

**Symptoms**: "EnrichR error: timeout" or connection issues

**Solutions**:
```r
# Use offline tools instead
seurat_obj <- add_pathway_weighted_scores(
    seurat_obj,
    known_tissue_type = "Immune system",
    enrichment_tools = c("fgsea", "go")  # Skip EnrichR
)
```

## Best Practices

1. **Use multiple tools**: Different tools have different strengths
2. **Check pathway support**: Don't blindly trust combined scores
3. **Investigate disagreements**: ScType-pathway conflicts reveal biology
4. **Consider tissue context**: Some tissues have less well-defined pathways
5. **Review low pathway support**: May indicate novel cell types
6. **Document thresholds**: Record which confidence cutoffs you use
7. **Report both scores**: ScType and pathway independently in publications

## Limitations

- **Knowledge base coverage**: Not all cell types have well-defined pathway signatures
- **Database bias**: Enrichment depends on well-studied pathways
- **Computational time**: Running enrichment on many clusters is slow
- **Internet required**: EnrichR needs web access
- **Species-specific**: Pathway databases are primarily human/mouse

## Citation

If you use pathway-weighted ScType in your research, please cite:

**ScType**:
Ianevski, A., Giri, A.K. & Aittokallio, T. Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data. Nat Commun 13, 1246 (2022). https://doi.org/10.1038/s41467-022-28803-w

**EnrichR**:
Chen et al. (2013) BMC Bioinformatics. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool.

**fgsea**:
Korotkevich et al. (2016) bioRxiv. Fast gene set enrichment analysis.

**clusterProfiler**:
Yu et al. (2012) OMICS. clusterProfiler: an R package for comparing biological themes among gene clusters.

---

*Last Updated: 2025-11-15*
