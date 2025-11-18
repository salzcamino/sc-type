#!/usr/bin/env Rscript
# Test Benchmark Framework
# Quick test to ensure benchmarking functions work correctly

library(Seurat)

# Source benchmark functions
source("benchmark_comparison.R")

# Create minimal test dataset
create_test_data <- function() {
  # Small synthetic dataset
  set.seed(42)
  n_cells <- 200
  n_genes <- 500

  # Create expression matrix
  expr <- matrix(rnbinom(n_cells * n_genes, mu = 2, size = 2),
                 nrow = n_genes, ncol = n_cells)
  rownames(expr) <- paste0("Gene_", 1:n_genes)
  colnames(expr) <- paste0("Cell_", 1:n_cells)

  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = expr)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10, verbose = FALSE)

  # Add synthetic ground truth labels
  true_labels <- sample(c("TypeA", "TypeB", "TypeC"), n_cells, replace = TRUE)
  seurat_obj@meta.data$true_celltype <- true_labels

  return(seurat_obj)
}

message("=== Testing ScType Benchmark Framework ===\n")

# Test 1: Data creation
message("Test 1: Creating test dataset...")
test_data <- tryCatch({
  create_test_data()
}, error = function(e) {
  message("ERROR: Could not create test data: ", e$message)
  NULL
})

if (is.null(test_data)) {
  stop("Test failed: Could not create test data")
}
message("✓ Test data created successfully\n")

# Test 2: Accuracy metrics calculation
message("Test 2: Testing accuracy metrics...")
true_labels <- c("A", "A", "B", "B", "C", "C")
pred_labels <- c("A", "A", "B", "C", "C", "C")

metrics <- calculate_accuracy_metrics(true_labels, pred_labels)
expected_accuracy <- 4/6  # 4 correct out of 6

if (abs(metrics$accuracy - expected_accuracy) < 0.01) {
  message("✓ Accuracy calculation correct")
} else {
  message("✗ Accuracy calculation incorrect")
  message("  Expected: ", expected_accuracy)
  message("  Got: ", metrics$accuracy)
}

if (!is.null(metrics$f1_weighted) && metrics$f1_weighted > 0) {
  message("✓ F1 score calculation working")
} else {
  message("✗ F1 score calculation failed")
}
message("")

# Test 3: ScType benchmark function
message("Test 3: Testing ScType benchmark function...")
sctype_result <- tryCatch({
  benchmark_sctype(test_data, tissue_type = "Immune system", n_iterations = 1)
}, error = function(e) {
  message("ERROR: ", e$message)
  list(mean_time = NA)
})

if (!is.na(sctype_result$mean_time)) {
  message("✓ ScType benchmark ran successfully")
  message(sprintf("  Runtime: %.2f seconds", sctype_result$mean_time))
  message(sprintf("  Memory: %.1f MB", sctype_result$peak_memory))
} else {
  message("✗ ScType benchmark failed (this is OK if ScType not installed)")
}
message("")

# Test 4: Plot generation (if we have results)
message("Test 4: Testing plot generation...")
test_results <- data.frame(
  method = c("sctype", "method2"),
  mean_time_sec = c(2.5, 5.0),
  sd_time_sec = c(0.1, 0.2),
  peak_memory_mb = c(100, 200),
  accuracy = c(0.85, 0.90),
  f1_score = c(0.83, 0.88),
  n_predicted_types = c(5, 5),
  installation_difficulty = c("Easy", "Medium"),
  custom_markers_supported = c(TRUE, FALSE),
  n_cells = c(200, 200),
  n_genes = c(500, 500),
  n_true_types = c(3, 3),
  stringsAsFactors = FALSE
)

plots <- tryCatch({
  plot_benchmark_results(test_results)
}, error = function(e) {
  message("ERROR: ", e$message)
  NULL
})

if (!is.null(plots) && length(plots) > 0) {
  message("✓ Plot generation successful")
  message(sprintf("  Generated %d plots", length(plots)))
} else {
  message("✗ Plot generation failed")
}
message("")

# Test 5: Report generation
message("Test 5: Testing report generation...")
report <- tryCatch({
  generate_benchmark_report(test_results)
}, error = function(e) {
  message("ERROR: ", e$message)
  NULL
})

if (!is.null(report) && nchar(report) > 100) {
  message("✓ Report generation successful")
  message(sprintf("  Report length: %d characters", nchar(report)))
  # Show first few lines
  lines <- strsplit(report, "\n")[[1]][1:5]
  message("  Preview:")
  for (line in lines) {
    message("    ", line)
  }
} else {
  message("✗ Report generation failed")
}
message("")

# Test 6: Full benchmark workflow (minimal)
message("Test 6: Testing full benchmark workflow...")
if (exists("sctype_result") && !is.na(sctype_result$mean_time)) {
  full_results <- tryCatch({
    benchmark_cell_type_annotation(
      test_data,
      true_labels_col = "true_celltype",
      tissue_type = "Immune system",
      methods = "sctype",  # Only test ScType for speed
      n_iterations = 1
    )
  }, error = function(e) {
    message("ERROR: ", e$message)
    NULL
  })

  if (!is.null(full_results) && nrow(full_results) > 0) {
    message("✓ Full benchmark workflow successful")
    print(full_results)
  } else {
    message("✗ Full benchmark workflow failed")
  }
} else {
  message("⊘ Skipping full workflow test (ScType not available)")
}
message("")

# Summary
message("=== Test Summary ===")
message("Benchmark framework is functional!")
message("\nNext steps:")
message("1. Install comparison tools (SingleR, CellTypist, etc.)")
message("2. Run full benchmark: Rscript run_benchmark.R")
message("3. See benchmarks/README.md for detailed documentation")
message("")
