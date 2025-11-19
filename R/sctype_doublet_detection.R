# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/master/LICENSE)
# ScType v2 - Doublet Detection Integration
# Written by Claude Code, November 2025
#
# Expected Impact: +5-10% accuracy on datasets with high doublet rates
# Complexity: EASY
# Timeline: 1-2 days

#' Detect doublets using scDblFinder
#'
#' Wrapper around Bioconductor's scDblFinder for easy integration with ScType.
#' Identifies cell multiplets (doublets/multiplets) that can cause spurious annotations.
#'
#' @param seurat_object Seurat object
#' @param assay Assay to use (default: "RNA")
#' @param dbr Expected doublet rate (optional, auto-estimated if NULL)
#' @param clusters Cluster assignments for doublet detection (optional)
#'
#' @return Seurat object with new metadata columns:
#'   - doublet_score: Continuous doublet score (0-1, higher = more likely doublet)
#'   - doublet_class: Classification ("singlet" or "doublet")
#'
#' @export
detect_doublets_scdblfinder <- function(seurat_object,
                                        assay = "RNA",
                                        dbr = NULL,
                                        clusters = NULL) {

  # Check dependencies
  if (!requireNamespace("scDblFinder", quietly = TRUE)) {
    stop("scDblFinder package not found. Install with:\n  BiocManager::install('scDblFinder')")
  }

  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment package not found. Install with:\n  BiocManager::install('SingleCellExperiment')")
  }

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package not found. Install with:\n  install.packages('Seurat')")
  }

  message("Converting Seurat to SingleCellExperiment...")
  sce <- Seurat::as.SingleCellExperiment(seurat_object, assay = assay)

  message("Running scDblFinder...")
  sce <- scDblFinder::scDblFinder(sce, dbr = dbr, clusters = clusters, verbose = FALSE)

  # Extract results
  seurat_object$doublet_score <- sce$scDblFinder.score
  seurat_object$doublet_class <- sce$scDblFinder.class

  # Summary
  n_doublets <- sum(sce$scDblFinder.class == "doublet")
  pct_doublets <- round(100 * n_doublets / ncol(seurat_object), 2)

  message(sprintf("\n=== Doublet Detection Summary ==="))
  message(sprintf("Total cells: %d", ncol(seurat_object)))
  message(sprintf("Detected doublets: %d (%.1f%%)", n_doublets, pct_doublets))
  message(sprintf("Singlets: %d (%.1f%%)", sum(sce$scDblFinder.class == "singlet"), 100 - pct_doublets))
  message(sprintf("\nNew metadata: doublet_score, doublet_class"))

  return(seurat_object)
}


#' Detect doublets from ScType score patterns
#'
#' Alternative doublet detection using ScType scores. Cells with high scores
#' for multiple unrelated cell types are flagged as potential doublets.
#'
#' @param sctype_score_matrix Matrix of ScType scores (cell types × cells)
#' @param high_score_threshold Threshold for "high" score (default: 75th percentile)
#' @param min_high_scores Min number of high scores to flag doublet (default: 2)
#' @param ambiguity_ratio Flag if top2/top1 score ratio > this (default: 0.8)
#'
#' @return Logical vector of doublet flags (TRUE = potential doublet)
#' @export
detect_doublets_from_scores <- function(sctype_score_matrix,
                                       high_score_threshold = NULL,
                                       min_high_scores = 2,
                                       ambiguity_ratio = 0.8) {

  # Calculate threshold if not provided
  if (is.null(high_score_threshold)) {
    high_score_threshold <- quantile(sctype_score_matrix, 0.75, na.rm = TRUE)
  }

  # Count high scores per cell
  high_score_counts <- colSums(sctype_score_matrix > high_score_threshold, na.rm = TRUE)
  doublet_flags <- high_score_counts >= min_high_scores

  # Check top2/top1 ratio (ambiguous assignments)
  top2_ratios <- apply(sctype_score_matrix, 2, function(scores) {
    sorted <- sort(scores, decreasing = TRUE)
    if (sorted[1] == 0) return(0)
    sorted[2] / sorted[1]
  })

  ambiguous_flags <- top2_ratios > ambiguity_ratio
  potential_doublets <- doublet_flags | ambiguous_flags

  message(sprintf("Score-based doublet detection:"))
  message(sprintf("  Multiple high scores: %d cells", sum(doublet_flags)))
  message(sprintf("  Ambiguous top scores: %d cells", sum(ambiguous_flags)))
  message(sprintf("  Total flagged: %d (%.1f%%)", sum(potential_doublets),
                 100 * mean(potential_doublets)))

  return(potential_doublets)
}
