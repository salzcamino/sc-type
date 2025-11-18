# Refactoring TODO: Library Calls

## Status: ✅ COMPLETE

### Fixed (9/9 files)
- ✅ `R/sctype_hierarchical.R` - Fixed patchwork library() call
- ✅ `R/sctype_doublet_detection.R` - Fixed scDblFinder and SingleCellExperiment library() calls
- ✅ `R/sctype_visualize.R` - Fixed ggplot2, dplyr, ComplexHeatmap, circlize library() calls
- ✅ `R/sctype_uncertainty.R` - Fixed ggplot2, dplyr, Seurat, tidyr library() calls
- ✅ `R/sctype_pathway_enrichment.R` - Fixed Seurat, dplyr, enrichR, fgsea, msigdbr, clusterProfiler, org.Hs.eg.db, ggplot2 library() calls
- ✅ `R/sctype_visualize_sce.R` - Fixed ggplot2, dplyr, SingleCellExperiment, scater, ComplexHeatmap, circlize library() calls
- ✅ `R/sctype_uncertainty_sce.R` - Fixed ggplot2, dplyr, SingleCellExperiment, scater library() calls
- ✅ `R/sctype_hierarchical_sce.R` - Fixed ggplot2, patchwork library() calls
- ✅ `R/sctype_wrapper_sce.R` - Fixed ggplot2 library() call

### Summary

All library() calls inside functions have been refactored to use requireNamespace() checks with :: notation. This refactoring brings the package into compliance with R package best practices and CRAN guidelines.

**Total refactored:** ~45 library() calls across 7 files (in addition to the 2 files previously completed)

## Recommended Refactoring Pattern

```r
# BEFORE
library(ggplot2)
plot <- ggplot(data) + geom_point()

# AFTER
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
}
plot <- ggplot2::ggplot(data) + ggplot2::geom_point()
```

## Why This Matters

1. **CRAN Compliance**: `library()` inside functions can cause namespace conflicts
2. **Best Practice**: Packages should use `::` notation or `@importFrom` in NAMESPACE
3. **User Experience**: Clear error messages when optional dependencies missing
4. **Testing**: Easier to test with specific package versions

## Impact

- **Priority**: Medium
- **Effort**: ~2-3 hours for complete refactoring
- **Risk**: Low (mostly affects optional visualization features)

## Note

Core functionality (gene_sets_prepare, sctype_score, auto_detect) does NOT use library() calls and follows best practices. The issues are limited to visualization and advanced features that are optional.

## Automated Refactoring

A semi-automated script could be created to handle this, but manual review is recommended to ensure correct :: operator usage and proper error messages.
