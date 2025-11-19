# Tests for TF-IDF weighting functions

test_that("calculate_tfidf_weights returns correct format", {
  # Mock gene sets
  gs_list <- list(
    "CD8+ T cells" = c("CD8A", "CD8B", "CD3D", "CD3E"),
    "CD4+ T cells" = c("CD4", "CD3D", "CD3E"),
    "B cells" = c("CD19", "CD79A", "MS4A1")
  )

  # Mock expression matrix
  set.seed(42)
  expr_mat <- matrix(rnorm(1000 * 100), nrow = 1000, ncol = 100)
  all_genes <- c("CD8A", "CD8B", "CD3D", "CD3E", "CD4", "CD19", "CD79A", "MS4A1",
                paste0("GENE", 1:992))
  rownames(expr_mat) <- all_genes
  colnames(expr_mat) <- paste0("CELL", 1:100)

  # Calculate weights
  weights <- calculate_tfidf_weights(gs_list, expr_mat)

  expect_s3_class(weights, "data.frame")
  expect_named(weights, c("gene_", "score_marker_sensitivity", "tfidf_raw", "tf", "idf"))
  expect_true(all(weights$score_marker_sensitivity >= 0 & weights$score_marker_sensitivity <= 1))
  expect_equal(attr(weights, "weighting_method"), "tfidf")
})

test_that("calculate_tfidf_weights handles different TF methods", {
  gs_list <- list("CellType1" = c("GENE1", "GENE2"))

  set.seed(42)
  expr_mat <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
  rownames(expr_mat) <- paste0("GENE", 1:100)

  # Test different TF methods
  for (method in c("mean", "max", "median")) {
    weights <- calculate_tfidf_weights(gs_list, expr_mat, tf_method = method)
    expect_s3_class(weights, "data.frame")
    expect_equal(attr(weights, "tf_method"), method)
  }
})

test_that("calculate_frequency_weights matches original ScType", {
  gs_list <- list(
    "CellType1" = c("GENE1", "GENE2", "GENE3"),
    "CellType2" = c("GENE2", "GENE3", "GENE4"),  # GENE2, GENE3 shared
    "CellType3" = c("GENE5", "GENE6")
  )

  weights <- calculate_frequency_weights(gs_list)

  expect_s3_class(weights, "data.frame")
  expect_named(weights, c("gene_", "score_marker_sensitivity", "frequency"))

  # Rare markers should have higher weights
  gene5_weight <- weights$score_marker_sensitivity[weights$gene_ == "GENE5"]
  gene2_weight <- weights$score_marker_sensitivity[weights$gene_ == "GENE2"]
  expect_true(gene5_weight > gene2_weight)  # GENE5 appears in 1 type, GENE2 in 2
})

test_that("calculate_hybrid_weights combines methods correctly", {
  gs_list <- list(
    "CellType1" = c("GENE1", "GENE2"),
    "CellType2" = c("GENE3", "GENE4")
  )

  set.seed(42)
  expr_mat <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
  rownames(expr_mat) <- paste0("GENE", 1:100)

  # Test different combination methods
  for (method in c("geometric_mean", "arithmetic_mean", "max", "min")) {
    weights <- calculate_hybrid_weights(gs_list, expr_mat, combine_method = method)
    expect_s3_class(weights, "data.frame")
    expect_true("frequency_weight" %in% colnames(weights))
    expect_true("tfidf_weight" %in% colnames(weights))
    expect_equal(attr(weights, "combine_method"), method)
  }
})

test_that("compare_weighting_methods returns comprehensive comparison", {
  gs_list <- list(
    "CellType1" = c("GENE1", "GENE2", "GENE3"),
    "CellType2" = c("GENE4", "GENE5")
  )

  set.seed(42)
  expr_mat <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
  rownames(expr_mat) <- paste0("GENE", 1:100)

  comparison <- compare_weighting_methods(gs_list, expr_mat, top_n = 5)

  expect_type(comparison, "list")
  expect_named(comparison, c("comparison_table", "top_frequency", "top_tfidf",
                             "top_hybrid", "largest_rank_differences",
                             "summary_statistics"))

  expect_s3_class(comparison$comparison_table, "data.frame")
  expect_s3_class(comparison$summary_statistics, "data.frame")
})

test_that(".rescale_values helper works correctly", {
  # Test with scales package available
  x <- c(1, 2, 3, 4, 5)
  result <- ScType:::.rescale_values(x, to = c(0, 1))

  expect_length(result, 5)
  expect_true(min(result) >= 0)
  expect_true(max(result) <= 1)
  expect_equal(min(result), 0)
  expect_equal(max(result), 1)

  # Test zero range
  x_constant <- c(5, 5, 5, 5)
  result_constant <- ScType:::.rescale_values(x_constant, to = c(0, 1))
  expect_true(all(result_constant == 0.5))
})
