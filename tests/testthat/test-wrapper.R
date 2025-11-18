# Test suite for wrapper functions
# Tests high-level interface functions

context("Wrapper functions")

test_that("run_sctype validates Seurat object input", {
  skip_if_not_installed("Seurat")

  # Test with NULL
  expect_error(
    run_sctype(NULL),
    "missing"
  )

  # Test with non-Seurat object
  expect_error(
    run_sctype(list(a = 1, b = 2)),
    "Seurat object"
  )
})

test_that("sctype_source returns database URL", {
  db_url <- sctype_source()

  expect_type(db_url, "character")
  expect_match(db_url, "https://")
  expect_match(db_url, "ScTypeDB")
})
