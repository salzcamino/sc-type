#!/usr/bin/env Rscript
# Quick Benchmark Runner for ScType
#
# Usage:
#   Rscript run_benchmark.R [dataset_path] [tissue_type]
#
# Example:
#   Rscript run_benchmark.R pbmc3k.rds "Immune system"

# Load benchmarking framework
source("benchmark_comparison.R")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  message("Using default PBMC 3k dataset...")

  # Download and prepare PBMC 3k example
  library(Seurat)

  # Check if example data exists
  if (file.exists("../exampleData.RDS")) {
    message("Loading example data...")
    pbmc_data <- readRDS("../exampleData.RDS")

    # Create Seurat object from matrix
    pbmc <- CreateSeuratObject(counts = pbmc_data, project = "pbmc3k")
    pbmc <- NormalizeData(pbmc)
    pbmc <- FindVariableFeatures(pbmc)
    pbmc <- ScaleData(pbmc)
    pbmc <- RunPCA(pbmc)
    pbmc <- FindNeighbors(pbmc, dims = 1:10)
    pbmc <- FindClusters(pbmc, resolution = 0.8)
    pbmc <- RunUMAP(pbmc, dims = 1:10)

    # Add ground truth labels (for demonstration - normally you'd have these)
    # Here we'll use cluster-based "pseudo-ground-truth"
    cluster_labels <- c(
      "0" = "CD14+ Monocytes",
      "1" = "CD4+ T cells",
      "2" = "CD8+ T cells",
      "3" = "B cells",
      "4" = "NK cells",
      "5" = "CD14+ Monocytes",
      "6" = "Dendritic cells",
      "7" = "Megakaryocytes"
    )

    pbmc@meta.data$true_celltype <- cluster_labels[as.character(Idents(pbmc))]

    seurat_obj <- pbmc
    tissue_type <- "Immune system"

  } else {
    stop("Example data not found. Please provide a Seurat object path.")
  }

} else {
  dataset_path <- args[1]
  tissue_type <- ifelse(length(args) >= 2, args[2], "Immune system")

  message(sprintf("Loading dataset from: %s", dataset_path))
  seurat_obj <- readRDS(dataset_path)

  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
}

# Check which methods are available
available_methods <- c()

# ScType (always available via sourcing)
available_methods <- c(available_methods, "sctype")

# SingleR
if (requireNamespace("SingleR", quietly = TRUE) &&
    requireNamespace("celldex", quietly = TRUE)) {
  available_methods <- c(available_methods, "singler")
  message("✓ SingleR available")
} else {
  message("✗ SingleR not available (install: BiocManager::install(c('SingleR', 'celldex')))")
}

# CellTypist
if (requireNamespace("reticulate", quietly = TRUE)) {
  has_celltypist <- tryCatch({
    reticulate::py_module_available("celltypist")
  }, error = function(e) FALSE)

  if (has_celltypist) {
    available_methods <- c(available_methods, "celltypist")
    message("✓ CellTypist available")
  } else {
    message("✗ CellTypist not available (install: pip install celltypist)")
  }
} else {
  message("✗ CellTypist not available (reticulate needed)")
}

# Azimuth
if (requireNamespace("Azimuth", quietly = TRUE)) {
  available_methods <- c(available_methods, "azimuth")
  message("✓ Azimuth available")
} else {
  message("✗ Azimuth not available (install: remotes::install_github('satijalab/azimuth'))")
}

# scCATCH
if (requireNamespace("scCATCH", quietly = TRUE)) {
  available_methods <- c(available_methods, "sccatch")
  message("✓ scCATCH available")
} else {
  message("✗ scCATCH not available (install: BiocManager::install('scCATCH'))")
}

message(sprintf("\n%d methods available for benchmarking: %s\n",
                length(available_methods),
                paste(available_methods, collapse = ", ")))

# Run benchmark
message("Starting benchmark...")
message("This may take several minutes...\n")

benchmark_results <- benchmark_cell_type_annotation(
  seurat_object = seurat_obj,
  true_labels_col = "true_celltype",
  tissue_type = tissue_type,
  methods = available_methods,
  n_iterations = 3
)

# Display results
message("\n=== BENCHMARK RESULTS ===\n")
print(benchmark_results)

# Save results
output_dir <- sprintf("results_%s", format(Sys.time(), "%Y%m%d_%H%M%S"))
save_benchmark_results(benchmark_results, output_dir)

message(sprintf("\nResults saved to: %s", output_dir))
message("\nDone!")
