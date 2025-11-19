# Tests for gene_sets_prepare function

test_that("gene_sets_prepare handles valid input correctly", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("HGNChelper")

  # Use the short database for testing
  db_path <- system.file("ScTypeDB_short.xlsx", package = "ScType")
  if (db_path == "") {
    db_path <- "../../ScTypeDB_short.xlsx"
  }

  skip_if_not(file.exists(db_path), "Database file not found")

  # Test with Immune system tissue type
  result <- gene_sets_prepare(db_path, "Immune system")

  expect_type(result, "list")
  expect_named(result, c("gs_positive", "gs_negative"))
  expect_type(result$gs_positive, "list")
  expect_type(result$gs_negative, "list")
  expect_true(length(result$gs_positive) > 0)
})

test_that("gene_sets_prepare validates gene symbols", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("HGNChelper")

  db_path <- system.file("ScTypeDB_short.xlsx", package = "ScType")
  if (db_path == "") {
    db_path <- "../../ScTypeDB_short.xlsx"
  }

  skip_if_not(file.exists(db_path), "Database file not found")

  result <- gene_sets_prepare(db_path, "Immune system")

  # Check that gene symbols are uppercase
  all_genes <- unique(unlist(result$gs_positive))
  expect_true(all(toupper(all_genes) == all_genes))

  # Check no empty strings
  expect_false(any(all_genes == ""))
})

test_that("gene_sets_prepare handles empty markers gracefully", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("HGNChelper")

  db_path <- system.file("ScTypeDB_short.xlsx", package = "ScType")
  if (db_path == "") {
    db_path <- "../../ScTypeDB_short.xlsx"
  }

  skip_if_not(file.exists(db_path), "Database file not found")

  # Should handle tissue types with minimal markers
  expect_no_error({
    result <- gene_sets_prepare(db_path, "Immune system")
  })
})
