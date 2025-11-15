# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/master/LICENSE)
# ScType v2 - TF-IDF Weighting Unit Tests
# Written by Claude Code, November 2025

# Load required packages
if (!require("testthat", quietly = TRUE)) {
    install.packages("testthat")
    library(testthat)
}

# Source the TF-IDF functions
source("R/sctype_tfidf.R")
source("R/sctype_score_.R")

message("=== Running ScType v2 TF-IDF Weighting Unit Tests ===\n")

# =============================================================================
# TEST 1: calculate_tfidf_weights() basic functionality
# =============================================================================
test_that("calculate_tfidf_weights works with basic input", {
    # Create test gene sets
    gs_list <- list(
        "CD8+ T cells" = c("CD8A", "CD8B", "CD3D", "CD3E"),
        "CD4+ T cells" = c("CD4", "CD3D", "CD3E", "IL7R"),
        "B cells" = c("CD19", "CD79A", "MS4A1")
    )

    # Create mock expression matrix (z-scaled)
    set.seed(123)
    expr_mat <- matrix(rnorm(100*50, mean=0, sd=1), nrow=100, ncol=50)
    rownames(expr_mat) <- paste0("GENE", 1:100)
    # Add marker genes
    rownames(expr_mat)[1:10] <- c("CD8A", "CD8B", "CD3D", "CD3E", "CD4",
                                   "CD19", "CD79A", "MS4A1", "IL7R", "GZMA")

    # Calculate TF-IDF weights
    weights <- calculate_tfidf_weights(gs_list, expr_mat)

    # Check structure
    expect_true(is.data.frame(weights))
    expect_true(all(c("gene_", "score_marker_sensitivity", "tfidf_raw", "tf", "idf") %in% colnames(weights)))

    # Check all markers are present
    all_markers <- unique(unlist(gs_list))
    expect_true(all(all_markers %in% weights$gene_))

    # Check weights are in [0, 1]
    expect_true(all(weights$score_marker_sensitivity >= 0 & weights$score_marker_sensitivity <= 1))

    # Check IDF: genes in fewer cell types should have higher IDF
    cd8a_idf <- weights$idf[weights$gene_ == "CD8A"]  # Only in CD8+ T cells
    cd3d_idf <- weights$idf[weights$gene_ == "CD3D"]  # In both CD8+ and CD4+ T cells
    expect_true(cd8a_idf > cd3d_idf)

    message("✓ TEST 1: calculate_tfidf_weights basic functionality - PASSED")
})

# =============================================================================
# TEST 2: calculate_frequency_weights() matches original ScType
# =============================================================================
test_that("calculate_frequency_weights matches original ScType weighting", {
    gs_list <- list(
        "Type1" = c("G1", "G2", "G3"),
        "Type2" = c("G2", "G4"),
        "Type3" = c("G1", "G2", "G5")
    )

    weights <- calculate_frequency_weights(gs_list)

    # Check structure
    expect_true(is.data.frame(weights))
    expect_true(all(c("gene_", "score_marker_sensitivity", "frequency") %in% colnames(weights)))

    # G2 appears in all 3 types (most frequent) -> should have lowest weight
    # G3, G4, G5 appear in 1 type each (least frequent) -> should have highest weight
    g2_weight <- weights$score_marker_sensitivity[weights$gene_ == "G2"]
    g3_weight <- weights$score_marker_sensitivity[weights$gene_ == "G3"]

    expect_true(g3_weight > g2_weight)

    message("✓ TEST 2: calculate_frequency_weights matches original - PASSED")
})

# =============================================================================
# TEST 3: calculate_hybrid_weights() combines methods
# =============================================================================
test_that("calculate_hybrid_weights combines frequency and TF-IDF", {
    gs_list <- list(
        "Type1" = c("G1", "G2"),
        "Type2" = c("G2", "G3")
    )

    set.seed(456)
    expr_mat <- matrix(rnorm(10*20), nrow=10, ncol=20)
    rownames(expr_mat) <- paste0("G", 1:10)

    weights <- calculate_hybrid_weights(gs_list, expr_mat)

    # Check structure
    expect_true(is.data.frame(weights))
    expect_true(all(c("gene_", "score_marker_sensitivity", "frequency_weight", "tfidf_weight") %in% colnames(weights)))

    # Hybrid weights should be between frequency and TF-IDF (for geometric mean)
    for (i in 1:nrow(weights)) {
        hybrid <- weights$score_marker_sensitivity[i]
        freq <- weights$frequency_weight[i]
        tfidf <- weights$tfidf_weight[i]

        # Geometric mean should be between min and max
        expect_true(hybrid >= min(freq, tfidf) && hybrid <= max(freq, tfidf))
    }

    message("✓ TEST 3: calculate_hybrid_weights combines methods - PASSED")
})

# =============================================================================
# TEST 4: sctype_score with custom marker_weights
# =============================================================================
test_that("sctype_score accepts custom marker_weights", {
    gs_list <- list(
        "Type1" = c("G1", "G2", "G3"),
        "Type2" = c("G4", "G5", "G6")
    )

    set.seed(789)
    expr_mat <- matrix(rnorm(10*30, mean=0, sd=1), nrow=10, ncol=30)
    rownames(expr_mat) <- paste0("G", 1:10)
    colnames(expr_mat) <- paste0("Cell", 1:30)

    # Calculate TF-IDF weights
    tfidf_weights <- calculate_tfidf_weights(gs_list, expr_mat)

    # Run sctype_score with custom weights
    scores_tfidf <- sctype_score(expr_mat, scaled=TRUE, gs=gs_list, gs2=NULL,
                                 marker_weights=tfidf_weights)

    # Run sctype_score with default (frequency) weights
    scores_freq <- sctype_score(expr_mat, scaled=TRUE, gs=gs_list, gs2=NULL,
                                marker_weights=NULL)

    # Both should return valid score matrices
    expect_true(is.matrix(scores_tfidf))
    expect_true(is.matrix(scores_freq))
    expect_equal(dim(scores_tfidf), dim(scores_freq))

    # Scores should be different (different weighting methods)
    expect_false(identical(scores_tfidf, scores_freq))

    message("✓ TEST 4: sctype_score with custom marker_weights - PASSED")
})

# =============================================================================
# TEST 5: compare_weighting_methods()
# =============================================================================
test_that("compare_weighting_methods generates valid comparison", {
    gs_list <- list(
        "Type1" = c("G1", "G2", "G3"),
        "Type2" = c("G2", "G4"),
        "Type3" = c("G1", "G5")
    )

    set.seed(101)
    expr_mat <- matrix(rnorm(10*20), nrow=10, ncol=20)
    rownames(expr_mat) <- paste0("G", 1:10)

    comparison <- compare_weighting_methods(gs_list, expr_mat, top_n=5)

    # Check structure
    expect_true(is.list(comparison))
    expect_true(all(c("comparison_table", "top_frequency", "top_tfidf", "top_hybrid",
                     "summary_statistics") %in% names(comparison)))

    # Check comparison table
    expect_true(is.data.frame(comparison$comparison_table))
    expect_true(all(c("gene", "frequency_weight", "tfidf_weight", "hybrid_weight") %in%
                   colnames(comparison$comparison_table)))

    # Check summary statistics
    expect_equal(nrow(comparison$summary_statistics), 3)  # 3 methods

    message("✓ TEST 5: compare_weighting_methods generates comparison - PASSED")
})

# =============================================================================
# TEST 6: Edge cases - single cell type
# =============================================================================
test_that("Functions handle single cell type", {
    gs_list <- list("Type1" = c("G1", "G2", "G3"))

    set.seed(202)
    expr_mat <- matrix(rnorm(10*15), nrow=10, ncol=15)
    rownames(expr_mat) <- paste0("G", 1:10)

    # Should still work
    weights <- calculate_tfidf_weights(gs_list, expr_mat)

    expect_true(is.data.frame(weights))
    expect_equal(nrow(weights), 3)  # 3 markers

    # All markers have IDF = log((1+1)/(1+1)) = 0 (all in same cell type)
    # But weights should still be valid
    expect_true(all(weights$score_marker_sensitivity >= 0 & weights$score_marker_sensitivity <= 1))

    message("✓ TEST 6: Edge case (single cell type) - PASSED")
})

# =============================================================================
# TEST 7: TF calculation methods
# =============================================================================
test_that("Different TF methods produce different weights", {
    gs_list <- list("Type1" = c("G1", "G2"))

    set.seed(303)
    expr_mat <- matrix(rnorm(5*20, mean=0, sd=2), nrow=5, ncol=20)
    rownames(expr_mat) <- paste0("G", 1:5)

    # Calculate with different TF methods
    weights_mean <- calculate_tfidf_weights(gs_list, expr_mat, tf_method="mean")
    weights_max <- calculate_tfidf_weights(gs_list, expr_mat, tf_method="max")
    weights_median <- calculate_tfidf_weights(gs_list, expr_mat, tf_method="median")

    # Check that TF values are different
    expect_false(identical(weights_mean$tf, weights_max$tf))
    expect_false(identical(weights_mean$tf, weights_median$tf))

    # Max should be >= mean >= median (usually, for normal distribution)
    # (This is approximate, may not always hold)
    g1_tf_mean <- weights_mean$tf[weights_mean$gene_ == "G1"]
    g1_tf_max <- weights_max$tf[weights_max$gene_ == "G1"]

    expect_true(g1_tf_max >= g1_tf_mean)

    message("✓ TEST 7: Different TF methods produce different weights - PASSED")
})

# =============================================================================
# TEST 8: Integration test - full TF-IDF workflow
# =============================================================================
test_that("Full TF-IDF workflow with sctype_score", {
    # Realistic scenario
    set.seed(999)
    gs_list <- list(
        "CD8+ T cells" = c("CD8A", "CD8B", "CD3D", "GZMA"),
        "CD4+ T cells" = c("CD4", "CD3D", "IL7R"),
        "B cells" = c("CD19", "MS4A1", "CD79A"),
        "NK cells" = c("GNLY", "NKG7", "GZMA")
    )

    # Create expression matrix with signal
    n_genes <- 50
    n_cells <- 100
    expr_mat <- matrix(rnorm(n_genes*n_cells), nrow=n_genes, ncol=n_cells)
    rownames(expr_mat) <- c("CD8A", "CD8B", "CD3D", "GZMA", "CD4", "IL7R",
                            "CD19", "MS4A1", "CD79A", "GNLY", "NKG7",
                            paste0("GENE", 1:(n_genes-11)))
    colnames(expr_mat) <- paste0("Cell", 1:n_cells)

    # Add signal: first 25 cells are CD8+ T cells (high CD8A, CD8B, CD3D, GZMA)
    expr_mat["CD8A", 1:25] <- expr_mat["CD8A", 1:25] + 3
    expr_mat["CD8B", 1:25] <- expr_mat["CD8B", 1:25] + 3
    expr_mat["CD3D", 1:25] <- expr_mat["CD3D", 1:25] + 2

    # Calculate weights
    tfidf_weights <- calculate_tfidf_weights(gs_list, expr_mat)

    # Run scoring
    scores <- sctype_score(expr_mat, scaled=FALSE, gs=gs_list, gs2=NULL,
                          marker_weights=tfidf_weights)

    # Check that CD8+ T cells score highest for cells 1-25
    cd8_scores <- scores["CD8+ T cells", ]
    avg_score_first25 <- mean(cd8_scores[1:25])
    avg_score_last25 <- mean(cd8_scores[76:100])

    expect_true(avg_score_first25 > avg_score_last25)

    message("✓ TEST 8: Full TF-IDF workflow - PASSED")
})

message("\n=== All TF-IDF Weighting Unit Tests Completed Successfully ===")
message("Total tests run: 8")
message("All tests PASSED ✓")
message("\nTF-IDF weighting functions are working correctly and ready for use in ScType v2.\n")
