# Test suite for sctype_score function
# Tests the core scoring algorithm

context("ScType scoring function")

test_that("sctype_score rejects non-matrix input", {
  # Create test data
  gs_list <- list("Cell_Type_1" = c("GENE1", "GENE2"))

  # Test with data frame
  df <- data.frame(cell1 = c(1, 2), cell2 = c(3, 4))
  expect_error(
    sctype_score(df, scaled = TRUE, gs = gs_list),
    "scRNAseqData must be a matrix"
  )

  # Test with vector
  vec <- c(1, 2, 3, 4)
  expect_error(
    sctype_score(vec, scaled = TRUE, gs = gs_list),
    "scRNAseqData must be a matrix"
  )
})

test_that("sctype_score rejects empty matrix", {
  gs_list <- list("Cell_Type_1" = c("GENE1", "GENE2"))

  # Empty matrix
  empty_mat <- matrix(nrow = 0, ncol = 0)
  expect_error(
    sctype_score(empty_mat, scaled = TRUE, gs = gs_list),
    "zero dimensions"
  )

  # Matrix with no rows
  no_rows <- matrix(nrow = 0, ncol = 5)
  expect_error(
    sctype_score(no_rows, scaled = TRUE, gs = gs_list),
    "rows.*columns"
  )

  # Matrix with no columns
  no_cols <- matrix(nrow = 5, ncol = 0)
  expect_error(
    sctype_score(no_cols, scaled = TRUE, gs = gs_list),
    "rows.*columns"
  )
})

test_that("sctype_score works with valid input", {
  # Create simple test matrix
  test_matrix <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3,
    ncol = 2,
    dimnames = list(c("GENE1", "GENE2", "GENE3"), c("CELL1", "CELL2"))
  )

  # Create gene sets
  gs_positive <- list("Cell_Type_1" = c("GENE1", "GENE2"))
  gs_negative <- list("Cell_Type_1" = c("GENE3"))

  # Run scoring
  result <- sctype_score(
    scRNAseqData = test_matrix,
    scaled = TRUE,
    gs = gs_positive,
    gs2 = gs_negative
  )

  # Check output structure
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)  # One cell type
  expect_equal(ncol(result), 2)  # Two cells
  expect_equal(rownames(result), "Cell_Type_1")
})

test_that("sctype_score handles gene name case conversion", {
  # Test matrix with lowercase genes
  test_matrix <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    ncol = 2,
    dimnames = list(c("gene1", "gene2"), c("CELL1", "CELL2"))
  )

  # Gene sets with uppercase
  gs_positive <- list("Cell_Type_1" = c("GENE1", "GENE2"))

  # Should work with gene_names_to_uppercase = TRUE (default)
  result <- sctype_score(
    scRNAseqData = test_matrix,
    scaled = TRUE,
    gs = gs_positive,
    gene_names_to_uppercase = TRUE
  )

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
})

test_that("sctype_score handles multiple cell types", {
  # Create test matrix
  test_matrix <- matrix(
    runif(20),
    nrow = 5,
    ncol = 4,
    dimnames = list(paste0("GENE", 1:5), paste0("CELL", 1:4))
  )

  # Multiple cell types
  gs_positive <- list(
    "Cell_Type_1" = c("GENE1", "GENE2"),
    "Cell_Type_2" = c("GENE3", "GENE4"),
    "Cell_Type_3" = c("GENE5")
  )

  result <- sctype_score(
    scRNAseqData = test_matrix,
    scaled = TRUE,
    gs = gs_positive
  )

  expect_equal(nrow(result), 3)  # Three cell types
  expect_equal(ncol(result), 4)  # Four cells
  expect_equal(rownames(result), names(gs_positive))
})
