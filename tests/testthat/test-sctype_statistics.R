# Tests for ScType v2 statistical functions

test_that("calculate_zscore_pvalue handles vector input correctly", {
  scores <- c("0" = 45.2, "1" = 23.1, "2" = 67.8, "3" = 12.4)

  result <- calculate_zscore_pvalue(scores)

  expect_s3_class(result, "data.frame")
  expect_named(result, c("cluster", "score", "zscore", "pvalue"))
  expect_equal(nrow(result), 4)
  expect_equal(result$cluster, c("0", "1", "2", "3"))

  # Check z-scores are finite
  expect_true(all(is.finite(result$zscore)))

  # Check p-values are in [0, 1]
  expect_true(all(result$pvalue >= 0 & result$pvalue <= 1))
})

test_that("calculate_zscore_pvalue handles cluster sizes", {
  scores <- c("0" = 45.2, "1" = 23.1, "2" = 67.8)
  cluster_sizes <- c("0" = 500, "1" = 300, "2" = 800)

  result <- calculate_zscore_pvalue(scores, cluster_sizes)

  expect_true("ncells" %in% colnames(result))
  expect_equal(result$ncells, c(500, 300, 800))
})

test_that("calculate_zscore_pvalue handles data frame input", {
  scores_df <- data.frame(
    cluster = c("0", "1", "2"),
    score = c(45.2, 23.1, 67.8)
  )

  result <- calculate_zscore_pvalue(scores_df)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
})

test_that("apply_fdr_correction works with different methods", {
  scores <- c("0" = 45.2, "1" = 23.1, "2" = 67.8, "3" = 12.4)
  stat_results <- calculate_zscore_pvalue(scores)

  # Test BH method
  result_bh <- apply_fdr_correction(stat_results, method = "BH")
  expect_true("fdr" %in% colnames(result_bh))
  expect_equal(attr(result_bh, "fdr_method"), "BH")

  # Test bonferroni method
  result_bonf <- apply_fdr_correction(stat_results, method = "bonferroni")
  expect_true("fdr" %in% colnames(result_bonf))
  expect_equal(attr(result_bonf, "fdr_method"), "bonferroni")

  # FDR values should be >= p-values for BH
  expect_true(all(result_bh$fdr >= result_bh$pvalue))
})

test_that("assign_confidence_level categorizes correctly", {
  scores <- c("0" = 45.2, "1" = 23.1, "2" = 67.8, "3" = 12.4, "4" = 5.0)
  stat_results <- calculate_zscore_pvalue(scores)
  stat_results <- apply_fdr_correction(stat_results)

  result <- assign_confidence_level(stat_results)

  expect_true("confidence" %in% colnames(result))
  expect_s3_class(result$confidence, "factor")
  expect_true(is.ordered(result$confidence))

  # Check levels are correct
  expect_true(all(levels(result$confidence) %in% c("High", "Medium", "Low", "Very Low")))
})

test_that("assign_confidence_level handles custom thresholds", {
  scores <- c("0" = 45.2, "1" = 23.1, "2" = 67.8)
  stat_results <- calculate_zscore_pvalue(scores)
  stat_results <- apply_fdr_correction(stat_results)

  custom_thresholds <- c("VeryHigh" = 0.001, "High" = 0.01, "Medium" = 0.05)

  result <- assign_confidence_level(stat_results, fdr_thresholds = custom_thresholds)

  expect_true("confidence" %in% colnames(result))
  expect_true(all(levels(result$confidence) %in% c("VeryHigh", "High", "Medium", "Very Low")))
})

test_that("generate_null_distribution creates valid output", {
  skip_if_not_installed("Matrix")

  # Small test case
  set.seed(42)
  n_celltypes <- 5
  n_cells <- 100
  sctype_scores <- matrix(rnorm(n_celltypes * n_cells),
                         nrow = n_celltypes, ncol = n_cells)
  rownames(sctype_scores) <- paste0("CellType", 1:n_celltypes)
  cluster_assignments <- sample(1:3, n_cells, replace = TRUE)

  # Run with small number of permutations for speed
  result <- generate_null_distribution(
    sctype_scores,
    cluster_assignments,
    n_permutations = 50,
    seed = 42
  )

  expect_type(result, "list")
  expect_named(result, c("null_scores", "observed_scores", "empirical_pvalues",
                        "n_permutations", "seed"))

  # Check dimensions
  expect_equal(dim(result$null_scores), c(50, 3))  # 50 perms, 3 clusters
  expect_equal(length(result$observed_scores), 3)  # 3 clusters
  expect_equal(length(result$empirical_pvalues), 3)  # 3 clusters

  # Check p-values are valid
  expect_true(all(result$empirical_pvalues > 0 & result$empirical_pvalues <= 1))
})
