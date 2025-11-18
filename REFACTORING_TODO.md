# Refactoring TODO: Library Calls

## Status: Partial Completion

### Fixed (2/7 files)
- ✅ `R/sctype_hierarchical.R` - Fixed patchwork library() call
- ✅ `R/sctype_doublet_detection.R` - Fixed scDblFinder and SingleCellExperiment library() calls

### Remaining (5/7 files with 40+ library() calls)

These files still use `library()` inside functions, which is not best practice for R packages:

1. **R/sctype_visualize.R** (4 calls)
   - Lines 49-50: ggplot2, dplyr
   - Lines 345-346: ComplexHeatmap, circlize

2. **R/sctype_uncertainty.R** (11 calls)
   - Lines 256-257: ggplot2, dplyr
   - Lines 346-347: ggplot2, dplyr
   - Lines 413-414: ggplot2, Seurat
   - Lines 448-449: ggplot2, dplyr
   - Lines 489-491: ggplot2, dplyr, tidyr

3. **R/sctype_pathway_enrichment.R** (9 calls)
   - Lines 55-56: Seurat, dplyr
   - Line 205: enrichR
   - Lines 258-259: fgsea, msigdbr
   - Lines 302-303: clusterProfiler, org.Hs.eg.db
   - Lines 591-592: ggplot2, dplyr

4. **R/sctype_visualize_sce.R** (8 calls)
   - Lines 50-52: ggplot2, dplyr, SingleCellExperiment
   - Lines 214, 271, 330: scater
   - Lines 361-362: ComplexHeatmap, circlize

5. **R/sctype_uncertainty_sce.R** (10 calls)
   - Lines 254-256: ggplot2, dplyr, SingleCellExperiment
   - Lines 345-346: ggplot2, dplyr
   - Lines 412-413: ggplot2, SingleCellExperiment
   - Line 428: scater
   - Lines 454-456: ggplot2, dplyr, SingleCellExperiment
   - Lines 496-497: ggplot2, dplyr

6. **R/sctype_hierarchical_sce.R** (2 calls)
   - Lines 241-242: ggplot2, patchwork

7. **R/sctype_wrapper_sce.R** (1 call)
   - Line 149: ggplot2

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
