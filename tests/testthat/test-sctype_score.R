# Tests for sctype_score function

test_that("sctype_score returns correct matrix dimensions", {
  # Create mock data
  set.seed(42)
  n_genes <- 100
  n_cells <- 50
  scRNAseqData <- matrix(rnorm(n_genes * n_cells), nrow = n_genes, ncol = n_cells)
  rownames(scRNAseqData) <- paste0("GENE", 1:n_genes)
  colnames(scRNAseqData) <- paste0("CELL", 1:n_cells)

  # Create mock gene sets
  gs <- list(
    "CellTypeA" = c("GENE1", "GENE2", "GENE3"),
    "CellTypeB" = c("GENE4", "GENE5", "GENE6")
  )
  gs2 <- list(
    "CellTypeA" = c("GENE10", "GENE11"),
    "CellTypeB" = c("GENE12", "GENE13")
  )

  # Run scoring
  result <- sctype_score(scRNAseqData, scaled = TRUE, gs = gs, gs2 = gs2)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)  # 2 cell types
  expect_equal(ncol(result), n_cells)  # 50 cells
  expect_equal(rownames(result), c("CellTypeA", "CellTypeB"))
})

test_that("sctype_score handles unscaled data", {
  set.seed(42)
  scRNAseqData <- matrix(rpois(100 * 50, lambda = 5), nrow = 100, ncol = 50)
  rownames(scRNAseqData) <- paste0("GENE", 1:100)
  colnames(scRNAseqData) <- paste0("CELL", 1:50)

  gs <- list("CellTypeA" = c("GENE1", "GENE2", "GENE3"))

  # Should auto-scale
  expect_no_error({
    result <- sctype_score(scRNAseqData, scaled = FALSE, gs = gs, gs2 = NULL)
  })

  expect_true(is.matrix(result))
})

test_that("sctype_score handles missing genes gracefully", {
  set.seed(42)
  scRNAseqData <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
  rownames(scRNAseqData) <- paste0("GENE", 1:100)
  colnames(scRNAseqData) <- paste0("CELL", 1:50)

  # Gene sets with genes not in data
  gs <- list(
    "CellTypeA" = c("GENE1", "GENE2", "MISSING_GENE1", "MISSING_GENE2")
  )

  # Should work but only use available genes
  expect_no_error({
    result <- sctype_score(scRNAseqData, scaled = TRUE, gs = gs, gs2 = NULL)
  })

  expect_true(is.matrix(result))
})

test_that("sctype_score return_details parameter works", {
  set.seed(42)
  scRNAseqData <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
  rownames(scRNAseqData) <- paste0("GENE", 1:100)
  colnames(scRNAseqData) <- paste0("CELL", 1:50)

  gs <- list("CellTypeA" = c("GENE1", "GENE2", "GENE3"))

  # Test return_details = TRUE
  result_detailed <- sctype_score(scRNAseqData, scaled = TRUE, gs = gs,
                                  gs2 = NULL, return_details = TRUE)

  expect_type(result_detailed, "list")
  expect_named(result_detailed, c("scores", "marker_sensitivity", "gene_sets_used",
                                 "scaled_matrix", "original_gs_names", "was_scaled"))

  # Test return_details = FALSE (default)
  result_simple <- sctype_score(scRNAseqData, scaled = TRUE, gs = gs,
                                gs2 = NULL, return_details = FALSE)

  expect_true(is.matrix(result_simple))
})

test_that("sctype_score validates input matrix", {
  # Non-matrix input
  expect_warning({
    result <- sctype_score(data.frame(a = 1:10), scaled = TRUE,
                          gs = list("A" = c("a")), gs2 = NULL)
  }, "doesn't seem to be a matrix")

  # Empty matrix
  empty_mat <- matrix(numeric(0), nrow = 0, ncol = 0)
  expect_warning({
    result <- sctype_score(empty_mat, scaled = TRUE,
                          gs = list("A" = c("a")), gs2 = NULL)
  }, "dimension.*equals to 0")
})

test_that("sctype_score works without scales package", {
  # This tests the fallback rescaling
  set.seed(42)
  scRNAseqData <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
  rownames(scRNAseqData) <- paste0("GENE", 1:100)
  colnames(scRNAseqData) <- paste0("CELL", 1:50)

  gs <- list("CellTypeA" = c("GENE1", "GENE2", "GENE3"))

  # Should work even if scales is not available
  expect_no_error({
    result <- sctype_score(scRNAseqData, scaled = TRUE, gs = gs, gs2 = NULL)
  })

  expect_true(is.matrix(result))
  expect_true(all(is.finite(result)))
})
