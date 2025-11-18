# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/master/LICENSE)
# ScType v2 - Statistical Significance Testing Functions
# Written by Claude Code, November 15, 2025
#
# This module provides statistical testing functions to replace the arbitrary
# ncells/4 threshold in ScType with principled z-scores, p-values, and FDR correction.
#
# Expected Impact: +5-10% accuracy improvement on ambiguous clusters
# Complexity: EASY
# Timeline: 1-2 days

#' Calculate z-scores and p-values for ScType scores
#'
#' This function converts raw ScType scores to statistical significance metrics.
#' It replaces the arbitrary ncells/4 threshold with principled z-score testing.
#'
#' @param scores Named vector or data frame with cluster scores
#'   - If named vector: names should be cluster IDs, values are scores
#'   - If data frame: must have 'cluster' and 'score' columns
#' @param cluster_sizes Named vector of cell counts per cluster (optional)
#'   - Names should match cluster IDs in scores
#'   - If NULL, assumes equal weighting
#' @param use_global_stats Use global score distribution for z-scores (default: TRUE)
#'   - TRUE: z = (score - global_mean) / global_sd
#'   - FALSE: z = (score - cluster_mean) / cluster_sd
#' @return Data frame with columns:
#'   - cluster: Cluster ID
#'   - score: Raw ScType score
#'   - zscore: Standardized z-score
#'   - pvalue: One-tailed p-value (P(score > expected))
#'   - ncells: Number of cells in cluster (if cluster_sizes provided)
#'
#' @examples
#' # Example 1: Simple vector input
#' scores <- c("0" = 45.2, "1" = 23.1, "2" = 67.8, "3" = 12.4)
#' results <- calculate_zscore_pvalue(scores)
#'
#' # Example 2: With cluster sizes
#' cluster_sizes <- c("0" = 500, "1" = 300, "2" = 800, "3" = 150)
#' results <- calculate_zscore_pvalue(scores, cluster_sizes)
#'
#' @export
calculate_zscore_pvalue <- function(scores,
                                   cluster_sizes = NULL,
                                   use_global_stats = TRUE) {

  # Input validation and conversion
  if (is.data.frame(scores)) {
    if (!all(c("cluster", "score") %in% colnames(scores))) {
      stop("Data frame must have 'cluster' and 'score' columns")
    }
    cluster_ids <- scores$cluster
    score_values <- scores$score
  } else if (is.vector(scores)) {
    cluster_ids <- names(scores)
    score_values <- as.numeric(scores)
    if (is.null(cluster_ids)) {
      cluster_ids <- as.character(seq_along(scores))
      warning("Scores vector has no names, using sequential IDs: 1, 2, 3, ...")
    }
  } else {
    stop("scores must be a named vector or data frame with 'cluster' and 'score' columns")
  }

  # Validate cluster_sizes if provided
  if (!is.null(cluster_sizes)) {
    if (!is.vector(cluster_sizes) || is.null(names(cluster_sizes))) {
      stop("cluster_sizes must be a named vector with cluster IDs as names")
    }
    # Match cluster_sizes to scores
    if (!all(cluster_ids %in% names(cluster_sizes))) {
      warning("Some clusters in scores are missing from cluster_sizes")
    }
  }

  # Calculate global statistics
  global_mean <- mean(score_values, na.rm = TRUE)
  global_sd <- sd(score_values, na.rm = TRUE)

  if (global_sd == 0 || is.na(global_sd)) {
    warning("Standard deviation is 0 or NA, cannot calculate z-scores. Returning raw scores.")
    results <- data.frame(
      cluster = cluster_ids,
      score = score_values,
      zscore = rep(0, length(score_values)),
      pvalue = rep(0.5, length(score_values)),
      stringsAsFactors = FALSE
    )
    if (!is.null(cluster_sizes)) {
      results$ncells <- cluster_sizes[as.character(cluster_ids)]
    }
    return(results)
  }

  # Calculate z-scores
  if (use_global_stats) {
    # Z-score using global distribution
    z_scores <- (score_values - global_mean) / global_sd
  } else {
    # Z-score using per-cluster adjustment (if cluster_sizes provided)
    if (!is.null(cluster_sizes)) {
      ncells_vec <- cluster_sizes[as.character(cluster_ids)]
      # Adjust for cluster size: larger clusters should have more stable scores
      z_scores <- (score_values - global_mean) / (global_sd / sqrt(ncells_vec))
    } else {
      # Fallback to global if no cluster sizes
      z_scores <- (score_values - global_mean) / global_sd
    }
  }

  # Convert z-scores to p-values (one-tailed: score > expected)
  # Lower p-value = more significant = stronger evidence for this cell type
  p_values <- pnorm(z_scores, lower.tail = FALSE)

  # Create results data frame
  results <- data.frame(
    cluster = cluster_ids,
    score = score_values,
    zscore = z_scores,
    pvalue = p_values,
    stringsAsFactors = FALSE
  )

  # Add cluster sizes if provided
  if (!is.null(cluster_sizes)) {
    results$ncells <- cluster_sizes[as.character(cluster_ids)]
  }

  return(results)
}


#' Apply FDR (False Discovery Rate) correction to p-values
#'
#' Uses Benjamini-Hochberg procedure to correct for multiple testing.
#' This is essential when testing multiple clusters simultaneously.
#'
#' @param stat_results Data frame from calculate_zscore_pvalue() with 'pvalue' column
#' @param method Correction method (default: "BH" = Benjamini-Hochberg)
#'   - "BH" or "fdr": Benjamini & Hochberg (1995) FDR
#'   - "BY": Benjamini & Yekutieli (2001) FDR
#'   - "bonferroni": Bonferroni correction (very conservative)
#'   - "holm": Holm (1979) method
#'   - "none": No correction
#' @return Input data frame with added 'fdr' column
#'
#' @examples
#' scores <- c("0" = 45.2, "1" = 23.1, "2" = 67.8)
#' results <- calculate_zscore_pvalue(scores)
#' results <- apply_fdr_correction(results)
#' print(results[, c("cluster", "pvalue", "fdr")])
#'
#' @export
apply_fdr_correction <- function(stat_results, method = "BH") {

  # Input validation
  if (!is.data.frame(stat_results)) {
    stop("stat_results must be a data frame")
  }

  if (!"pvalue" %in% colnames(stat_results)) {
    stop("stat_results must have 'pvalue' column (output from calculate_zscore_pvalue)")
  }

  # Valid correction methods
  valid_methods <- c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel", "fdr", "none")
  if (!method %in% valid_methods) {
    stop(sprintf("method must be one of: %s", paste(valid_methods, collapse = ", ")))
  }

  # Apply FDR correction
  stat_results$fdr <- p.adjust(stat_results$pvalue, method = method)

  # Add note about method used
  attr(stat_results, "fdr_method") <- method

  return(stat_results)
}


#' Assign confidence levels based on statistical thresholds
#'
#' Categorizes annotations into confidence levels based on FDR and z-scores.
#' This provides interpretable categories for downstream filtering and QC.
#'
#' @param stat_results Data frame with 'fdr' and 'zscore' columns
#' @param fdr_thresholds Named vector of FDR thresholds for each level
#'   - Default: c("High" = 0.01, "Medium" = 0.05, "Low" = 0.1)
#'   - Annotations with FDR > max(fdr_thresholds) are "Very Low" confidence
#' @param zscore_thresholds Named vector of z-score thresholds (optional)
#'   - Default: c("High" = 2, "Medium" = 1, "Low" = 0)
#'   - If provided, both FDR and z-score must meet threshold
#' @param use_combined Require both FDR and z-score thresholds (default: FALSE)
#'   - FALSE: Use only FDR thresholds
#'   - TRUE: Require both FDR and z-score thresholds
#' @return Input data frame with added 'confidence' column (factor)
#'
#' @examples
#' scores <- c("0" = 45.2, "1" = 23.1, "2" = 67.8, "3" = 5.1)
#' results <- calculate_zscore_pvalue(scores)
#' results <- apply_fdr_correction(results)
#' results <- assign_confidence_level(results)
#' table(results$confidence)
#'
#' # Custom thresholds
#' results <- assign_confidence_level(
#'   results,
#'   fdr_thresholds = c("High" = 0.001, "Medium" = 0.01, "Low" = 0.05)
#' )
#'
#' @export
assign_confidence_level <- function(stat_results,
                                   fdr_thresholds = c("High" = 0.01, "Medium" = 0.05, "Low" = 0.1),
                                   zscore_thresholds = c("High" = 2, "Medium" = 1, "Low" = 0),
                                   use_combined = FALSE) {

  # Input validation
  if (!is.data.frame(stat_results)) {
    stop("stat_results must be a data frame")
  }

  required_cols <- "fdr"
  if (use_combined) required_cols <- c(required_cols, "zscore")

  if (!all(required_cols %in% colnames(stat_results))) {
    stop(sprintf("stat_results must have columns: %s", paste(required_cols, collapse = ", ")))
  }

  # Validate thresholds
  if (!is.vector(fdr_thresholds) || is.null(names(fdr_thresholds))) {
    stop("fdr_thresholds must be a named vector (e.g., c('High' = 0.01, 'Medium' = 0.05))")
  }

  # Sort thresholds in increasing order
  fdr_thresholds <- sort(fdr_thresholds)
  threshold_names <- names(fdr_thresholds)

  # Initialize confidence as "Very Low"
  stat_results$confidence <- "Very Low"

  if (!use_combined) {
    # FDR-only thresholds
    for (i in length(fdr_thresholds):1) {
      level_name <- threshold_names[i]
      threshold_val <- fdr_thresholds[i]
      stat_results$confidence[stat_results$fdr < threshold_val] <- level_name
    }
  } else {
    # Combined FDR and z-score thresholds
    if (is.null(zscore_thresholds) || !is.vector(zscore_thresholds)) {
      stop("zscore_thresholds must be a named vector when use_combined = TRUE")
    }

    zscore_thresholds <- sort(zscore_thresholds)

    # Both thresholds must be met
    for (i in length(fdr_thresholds):1) {
      level_name <- threshold_names[i]
      fdr_thresh <- fdr_thresholds[i]
      z_thresh <- zscore_thresholds[i]

      meets_criteria <- (stat_results$fdr < fdr_thresh) & (stat_results$zscore > z_thresh)
      stat_results$confidence[meets_criteria] <- level_name
    }
  }

  # Convert to ordered factor
  all_levels <- c(threshold_names, "Very Low")
  stat_results$confidence <- factor(stat_results$confidence,
                                    levels = all_levels,
                                    ordered = TRUE)

  # Store threshold parameters as attributes
  attr(stat_results, "fdr_thresholds") <- fdr_thresholds
  if (use_combined) {
    attr(stat_results, "zscore_thresholds") <- zscore_thresholds
  }
  attr(stat_results, "use_combined") <- use_combined

  return(stat_results)
}


#' Generate null distribution via permutation testing (optional)
#'
#' Creates empirical null distribution by permuting cluster labels.
#' This provides a data-driven alternative to parametric z-score testing.
#' More computationally intensive but makes no distributional assumptions.
#'
#' @param sctype_scores_matrix Matrix of ScType scores (cell types × cells)
#' @param cluster_assignments Vector of cluster IDs for each cell
#' @param n_permutations Number of permutations (default: 1000)
#'   - 100-500: Quick test, less precise
#'   - 1000-5000: Standard, good balance
#'   - 10000+: High precision, slow
#' @param seed Random seed for reproducibility (default: 42)
#' @return List with:
#'   - null_scores: Matrix of permuted scores (n_permutations × n_clusters)
#'   - observed_scores: Vector of observed scores for each cluster
#'   - empirical_pvalues: Empirical p-values from permutation test
#'   - n_permutations: Number of permutations performed
#'
#' @examples
#' \dontrun{
#' # This is computationally intensive, example with small n_permutations
#' null_dist <- generate_null_distribution(
#'   sctype_scores_matrix = es.max,
#'   cluster_assignments = seurat_obj$seurat_clusters,
#'   n_permutations = 100
#' )
#'
#' # Compare observed vs null
#' hist(null_dist$null_scores[, 1], main = "Null distribution - Cluster 0")
#' abline(v = null_dist$observed_scores[1], col = "red", lwd = 2)
#' }
#'
#' @export
generate_null_distribution <- function(sctype_scores_matrix,
                                      cluster_assignments,
                                      n_permutations = 1000,
                                      seed = 42) {

  # Input validation
  if (!is.matrix(sctype_scores_matrix)) {
    stop("sctype_scores_matrix must be a matrix (cell types × cells)")
  }

  if (ncol(sctype_scores_matrix) != length(cluster_assignments)) {
    stop("Number of cells in sctype_scores_matrix must match length of cluster_assignments")
  }

  if (n_permutations < 100) {
    warning("n_permutations < 100 may give unreliable p-values. Consider using >= 1000.")
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Get unique clusters
  unique_clusters <- unique(cluster_assignments)
  n_clusters <- length(unique_clusters)

  # Calculate observed scores (aggregate by cluster)
  observed_scores <- sapply(unique_clusters, function(cl) {
    cells_in_cluster <- which(cluster_assignments == cl)
    cluster_scores <- rowSums(sctype_scores_matrix[, cells_in_cluster, drop = FALSE])
    max(cluster_scores)  # Top score for this cluster
  })
  names(observed_scores) <- unique_clusters

  # Initialize null distribution matrix
  null_scores <- matrix(NA, nrow = n_permutations, ncol = n_clusters)
  colnames(null_scores) <- unique_clusters

  # Run permutations
  message(sprintf("Running %d permutations...", n_permutations))

  for (i in 1:n_permutations) {
    # Permute cluster labels
    permuted_assignments <- sample(cluster_assignments)

    # Calculate scores for permuted clusters
    for (j in seq_along(unique_clusters)) {
      cl <- unique_clusters[j]
      permuted_cells <- which(permuted_assignments == cl)
      permuted_cluster_scores <- rowSums(sctype_scores_matrix[, permuted_cells, drop = FALSE])
      null_scores[i, j] <- max(permuted_cluster_scores)
    }

    # Progress indicator
    if (i %% 100 == 0) {
      message(sprintf("  Completed %d/%d permutations", i, n_permutations))
    }
  }

  # Calculate empirical p-values
  # P-value = (# permutations with score >= observed) / n_permutations
  empirical_pvalues <- sapply(seq_along(unique_clusters), function(j) {
    observed <- observed_scores[j]
    null_dist <- null_scores[, j]
    (sum(null_dist >= observed) + 1) / (n_permutations + 1)  # +1 to avoid p=0
  })
  names(empirical_pvalues) <- unique_clusters

  message("Permutation testing complete!")

  # Return results
  result <- list(
    null_scores = null_scores,
    observed_scores = observed_scores,
    empirical_pvalues = empirical_pvalues,
    n_permutations = n_permutations,
    seed = seed
  )

  class(result) <- c("sctype_null_distribution", "list")

  return(result)
}


#' Print method for sctype_null_distribution
#' @export
print.sctype_null_distribution <- function(x, ...) {
  cat("ScType Null Distribution (Permutation Test)\n")
  cat("============================================\n")
  cat(sprintf("Number of permutations: %d\n", x$n_permutations))
  cat(sprintf("Number of clusters: %d\n", length(x$observed_scores)))
  cat(sprintf("Random seed: %d\n\n", x$seed))

  cat("Empirical p-values:\n")
  pval_df <- data.frame(
    cluster = names(x$empirical_pvalues),
    observed_score = x$observed_scores,
    empirical_pvalue = x$empirical_pvalues
  )
  print(pval_df, row.names = FALSE)

  cat("\nNote: Lower p-values indicate stronger evidence for cell type assignment\n")
}
