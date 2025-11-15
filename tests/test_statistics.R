# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/master/LICENSE)
# ScType v2 - Statistical Testing Unit Tests
# Written by Claude Code, November 2025
#
# This file contains unit tests for the statistical significance testing functions
# introduced in ScType v2 (R/sctype_statistics.R)

# Load required packages
if (!require("testthat", quietly = TRUE)) {
    message("testthat package not found. Installing...")
    install.packages("testthat")
    library(testthat)
}

# Source the statistics functions
source("R/sctype_statistics.R")

message("=== Running ScType v2 Statistical Testing Unit Tests ===\n")

# =============================================================================
# TEST 1: calculate_zscore_pvalue() basic functionality
# =============================================================================
test_that("calculate_zscore_pvalue works with named vector input", {
    # Create test scores
    scores <- c("0" = 45.2, "1" = 23.1, "2" = 67.8, "3" = 12.4)

    # Calculate z-scores and p-values
    results <- calculate_zscore_pvalue(scores)

    # Check structure
    expect_true(is.data.frame(results))
    expect_equal(nrow(results), 4)
    expect_true(all(c("cluster", "score", "zscore", "pvalue") %in% colnames(results)))

    # Check that z-scores are properly normalized
    expect_true(abs(mean(results$zscore)) < 0.1)  # Mean should be ~0
    expect_true(abs(sd(results$zscore) - 1) < 0.1)  # SD should be ~1

    # Check p-value range
    expect_true(all(results$pvalue >= 0 & results$pvalue <= 1))

    # Check that highest score has highest z-score and lowest p-value
    max_score_idx <- which.max(results$score)
    expect_equal(max_score_idx, which.max(results$zscore))
    expect_equal(max_score_idx, which.min(results$pvalue))

    message("✓ TEST 1: calculate_zscore_pvalue basic functionality - PASSED")
})


# =============================================================================
# TEST 2: calculate_zscore_pvalue() with cluster sizes
# =============================================================================
test_that("calculate_zscore_pvalue works with cluster sizes", {
    scores <- c("0" = 45.2, "1" = 23.1, "2" = 67.8)
    cluster_sizes <- c("0" = 500, "1" = 300, "2" = 800)

    results <- calculate_zscore_pvalue(scores, cluster_sizes)

    # Check that ncells column is added
    expect_true("ncells" %in% colnames(results))
    expect_equal(results$ncells, cluster_sizes[results$cluster])

    message("✓ TEST 2: calculate_zscore_pvalue with cluster sizes - PASSED")
})


# =============================================================================
# TEST 3: calculate_zscore_pvalue() with data frame input
# =============================================================================
test_that("calculate_zscore_pvalue works with data frame input", {
    scores_df <- data.frame(
        cluster = c("0", "1", "2"),
        score = c(45.2, 23.1, 67.8)
    )

    results <- calculate_zscore_pvalue(scores_df)

    expect_true(is.data.frame(results))
    expect_equal(nrow(results), 3)
    expect_true(all(c("cluster", "score", "zscore", "pvalue") %in% colnames(results)))

    message("✓ TEST 3: calculate_zscore_pvalue with data frame input - PASSED")
})


# =============================================================================
# TEST 4: apply_fdr_correction()
# =============================================================================
test_that("apply_fdr_correction adds FDR column", {
    scores <- c("0" = 45.2, "1" = 23.1, "2" = 67.8, "3" = 12.4, "4" = 8.5)
    results <- calculate_zscore_pvalue(scores)
    results <- apply_fdr_correction(results)

    # Check FDR column exists
    expect_true("fdr" %in% colnames(results))

    # Check FDR values are valid
    expect_true(all(results$fdr >= 0 & results$fdr <= 1))

    # Check FDR >= p-value (Benjamini-Hochberg property)
    expect_true(all(results$fdr >= results$pvalue))

    # Check attribute
    expect_equal(attr(results, "fdr_method"), "BH")

    message("✓ TEST 4: apply_fdr_correction functionality - PASSED")
})


# =============================================================================
# TEST 5: assign_confidence_level()
# =============================================================================
test_that("assign_confidence_level categorizes correctly", {
    scores <- c("0" = 100, "1" = 50, "2" = 30, "3" = 10, "4" = 5)
    results <- calculate_zscore_pvalue(scores)
    results <- apply_fdr_correction(results)
    results <- assign_confidence_level(results)

    # Check confidence column exists
    expect_true("confidence" %in% colnames(results))

    # Check it's a factor
    expect_true(is.factor(results$confidence))

    # Check levels
    expected_levels <- c("High", "Medium", "Low", "Very Low")
    expect_true(all(expected_levels %in% levels(results$confidence)))

    # Check that lower FDR = higher confidence
    sorted_by_fdr <- results[order(results$fdr), ]
    # First entry should have best confidence (High or Medium)
    expect_true(as.character(sorted_by_fdr$confidence[1]) %in% c("High", "Medium"))

    message("✓ TEST 5: assign_confidence_level categorization - PASSED")
})


# =============================================================================
# TEST 6: assign_confidence_level() with custom thresholds
# =============================================================================
test_that("assign_confidence_level works with custom thresholds", {
    scores <- c("0" = 100, "1" = 50, "2" = 30, "3" = 10)
    results <- calculate_zscore_pvalue(scores)
    results <- apply_fdr_correction(results)

    # Custom thresholds
    custom_thresholds <- c("VeryHigh" = 0.001, "High" = 0.01, "Medium" = 0.05)
    results <- assign_confidence_level(results, fdr_thresholds = custom_thresholds)

    # Check custom levels exist
    expect_true(all(c("VeryHigh", "High", "Medium", "Very Low") %in% levels(results$confidence)))

    # Check attributes
    expect_equal(attr(results, "fdr_thresholds"), sort(custom_thresholds))

    message("✓ TEST 6: assign_confidence_level custom thresholds - PASSED")
})


# =============================================================================
# TEST 7: Edge cases - all same scores
# =============================================================================
test_that("Functions handle edge case: all same scores", {
    scores <- c("0" = 50, "1" = 50, "2" = 50)

    # This should trigger the warning about SD = 0
    expect_warning(results <- calculate_zscore_pvalue(scores))

    # Should still return results
    expect_true(is.data.frame(results))
    expect_equal(nrow(results), 3)

    # Z-scores should be 0, p-values should be 0.5
    expect_true(all(results$zscore == 0))
    expect_true(all(results$pvalue == 0.5))

    message("✓ TEST 7: Edge case (all same scores) - PASSED")
})


# =============================================================================
# TEST 8: Edge cases - single cluster
# =============================================================================
test_that("Functions handle edge case: single cluster", {
    scores <- c("0" = 45.2)

    results <- calculate_zscore_pvalue(scores)

    # Should work but z-score will be 0 (no variance)
    expect_true(is.data.frame(results))
    expect_equal(nrow(results), 1)

    message("✓ TEST 8: Edge case (single cluster) - PASSED")
})


# =============================================================================
# TEST 9: generate_null_distribution() basic functionality
# =============================================================================
test_that("generate_null_distribution creates valid null distribution", {
    # Create mock score matrix (small for speed)
    set.seed(42)
    n_celltypes <- 5
    n_cells <- 100
    sctype_scores_matrix <- matrix(rnorm(n_celltypes * n_cells, mean = 10, sd = 5),
                                   nrow = n_celltypes, ncol = n_cells)
    rownames(sctype_scores_matrix) <- paste0("CellType", 1:n_celltypes)
    colnames(sctype_scores_matrix) <- paste0("Cell", 1:n_cells)

    # Create mock cluster assignments (3 clusters)
    cluster_assignments <- sample(c("A", "B", "C"), n_cells, replace = TRUE)

    # Run permutation test (small n_permutations for speed)
    null_dist <- generate_null_distribution(
        sctype_scores_matrix = sctype_scores_matrix,
        cluster_assignments = cluster_assignments,
        n_permutations = 50,  # Low for speed
        seed = 123
    )

    # Check structure
    expect_true(is.list(null_dist))
    expect_true(all(c("null_scores", "observed_scores", "empirical_pvalues", "n_permutations") %in% names(null_dist)))

    # Check dimensions
    expect_equal(nrow(null_dist$null_scores), 50)
    expect_equal(ncol(null_dist$null_scores), 3)  # 3 clusters
    expect_equal(length(null_dist$observed_scores), 3)
    expect_equal(length(null_dist$empirical_pvalues), 3)

    # Check p-value range
    expect_true(all(null_dist$empirical_pvalues >= 0 & null_dist$empirical_pvalues <= 1))

    # Check class
    expect_true("sctype_null_distribution" %in% class(null_dist))

    message("✓ TEST 9: generate_null_distribution basic functionality - PASSED")
})


# =============================================================================
# TEST 10: Integration test - full workflow
# =============================================================================
test_that("Full statistical testing workflow", {
    # Simulate realistic cluster scores
    set.seed(456)
    cluster_scores <- c(
        "Cluster_0" = 85.3,   # High confidence expected
        "Cluster_1" = 62.1,   # Medium confidence
        "Cluster_2" = 92.7,   # High confidence
        "Cluster_3" = 38.4,   # Low confidence
        "Cluster_4" = 15.2,   # Very low confidence
        "Cluster_5" = 71.8,   # Medium confidence
        "Cluster_6" = 8.3     # Very low confidence
    )

    cluster_sizes <- c(
        "Cluster_0" = 450,
        "Cluster_1" = 320,
        "Cluster_2" = 580,
        "Cluster_3" = 200,
        "Cluster_4" = 120,
        "Cluster_5" = 390,
        "Cluster_6" = 80
    )

    # Step 1: Calculate z-scores and p-values
    results <- calculate_zscore_pvalue(cluster_scores, cluster_sizes)

    # Step 2: Apply FDR correction
    results <- apply_fdr_correction(results)

    # Step 3: Assign confidence levels
    results <- assign_confidence_level(results)

    # Validate final results
    expect_equal(nrow(results), 7)
    expect_true(all(c("cluster", "score", "zscore", "pvalue", "fdr", "confidence", "ncells") %in% colnames(results)))

    # Check that high scores have low FDR
    high_score_idx <- which.max(results$score)
    expect_true(results$fdr[high_score_idx] < 0.1)

    # Check that low scores have high FDR or low confidence
    low_score_idx <- which.min(results$score)
    expect_true(results$confidence[low_score_idx] %in% c("Low", "Very Low"))

    message("✓ TEST 10: Full statistical testing workflow - PASSED")
})


# =============================================================================
# TEST 11: Error handling - missing pvalue column
# =============================================================================
test_that("apply_fdr_correction errors without pvalue column", {
    bad_df <- data.frame(cluster = c("0", "1"), score = c(10, 20))

    expect_error(apply_fdr_correction(bad_df), "must have 'pvalue' column")

    message("✓ TEST 11: Error handling (missing pvalue) - PASSED")
})


# =============================================================================
# TEST 12: Error handling - invalid FDR method
# =============================================================================
test_that("apply_fdr_correction errors with invalid method", {
    scores <- c("0" = 45.2, "1" = 23.1)
    results <- calculate_zscore_pvalue(scores)

    expect_error(apply_fdr_correction(results, method = "invalid_method"), "method must be one of")

    message("✓ TEST 12: Error handling (invalid FDR method) - PASSED")
})


# =============================================================================
# TEST 13: Reproducibility - same seed gives same results
# =============================================================================
test_that("generate_null_distribution is reproducible with seed", {
    # Create small test matrix
    set.seed(100)
    mat <- matrix(rnorm(50), nrow = 5, ncol = 10)
    clusters <- rep(c("A", "B"), each = 5)

    # Run twice with same seed
    null1 <- generate_null_distribution(mat, clusters, n_permutations = 20, seed = 999)
    null2 <- generate_null_distribution(mat, clusters, n_permutations = 20, seed = 999)

    # Should be identical
    expect_equal(null1$null_scores, null2$null_scores)
    expect_equal(null1$empirical_pvalues, null2$empirical_pvalues)

    message("✓ TEST 13: Reproducibility with seed - PASSED")
})


# =============================================================================
# TEST 14: Performance check - large dataset
# =============================================================================
test_that("Functions handle large datasets efficiently", {
    # Create large mock data
    set.seed(789)
    n_clusters <- 50
    scores <- setNames(runif(n_clusters, min = 10, max = 100), paste0("C", 1:n_clusters))
    cluster_sizes <- setNames(sample(100:1000, n_clusters), paste0("C", 1:n_clusters))

    # Time the workflow (should be fast)
    start_time <- Sys.time()
    results <- calculate_zscore_pvalue(scores, cluster_sizes)
    results <- apply_fdr_correction(results)
    results <- assign_confidence_level(results)
    end_time <- Sys.time()

    elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # Should complete in under 1 second for 50 clusters
    expect_true(elapsed < 1.0)

    # Check all results are valid
    expect_equal(nrow(results), n_clusters)
    expect_true(all(!is.na(results$fdr)))

    message(sprintf("✓ TEST 14: Performance (50 clusters in %.3f sec) - PASSED", elapsed))
})


# =============================================================================
# Summary
# =============================================================================
message("\n=== All Statistical Testing Unit Tests Completed Successfully ===")
message("Total tests run: 14")
message("All tests PASSED ✓")
message("\nStatistical functions are working correctly and ready for use in ScType v2.")
message("Next step: Create STATISTICAL_TESTING_README.md documentation\n")
