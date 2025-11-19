# ScType Benchmarking Framework
# Comparison with major cell type annotation tools
#
# This script compares ScType against leading alternatives:
# - SingleR (reference-based)
# - CellTypist (ML-based, via reticulate)
# - Azimuth (Seurat reference mapping)
# - scCATCH (tissue-specific markers)
# - Garnett (marker-based supervised)
#
# Metrics evaluated:
# 1. Speed (wall clock time)
# 2. Memory usage (peak RAM)
# 3. Accuracy (against known labels)
# 4. Ease of installation/use
# 5. Support for custom markers

#' @title Benchmark Cell Type Annotation Tools
#' @description Comprehensive benchmarking of cell type annotation methods
#' @param seurat_object Seurat object with known cell types in metadata
#' @param true_labels_col Column name containing ground truth labels
#' @param tissue_type Tissue type for ScType (e.g., "Immune system")
#' @param methods Vector of methods to benchmark (default: all available)
#' @param n_iterations Number of iterations for timing (default: 3)
#' @return Data frame with benchmark results
benchmark_cell_type_annotation <- function(seurat_object,
                                          true_labels_col = "cell_type",
                                          tissue_type = "Immune system",
                                          methods = c("sctype", "singler", "celltypist",
                                                     "azimuth", "sccatch", "garnett"),
                                          n_iterations = 3) {

  # Check dependencies
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat required for benchmarking")
  }

  # Validate input
  if (!true_labels_col %in% colnames(seurat_object@meta.data)) {
    stop(paste0("Column '", true_labels_col, "' not found in metadata"))
  }

  true_labels <- seurat_object@meta.data[[true_labels_col]]
  n_cells <- ncol(seurat_object)
  n_genes <- nrow(seurat_object)

  message(sprintf("Benchmarking on dataset: %d cells, %d genes", n_cells, n_genes))
  message(sprintf("True cell types: %d unique types", length(unique(true_labels))))

  results <- list()

  # === ScType ===
  if ("sctype" %in% methods) {
    message("\n=== Benchmarking ScType ===")
    results$sctype <- benchmark_sctype(seurat_object, tissue_type, n_iterations)
  }

  # === SingleR ===
  if ("singler" %in% methods) {
    message("\n=== Benchmarking SingleR ===")
    results$singler <- benchmark_singler(seurat_object, n_iterations)
  }

  # === CellTypist ===
  if ("celltypist" %in% methods) {
    message("\n=== Benchmarking CellTypist ===")
    results$celltypist <- benchmark_celltypist(seurat_object, n_iterations)
  }

  # === Azimuth ===
  if ("azimuth" %in% methods) {
    message("\n=== Benchmarking Azimuth ===")
    results$azimuth <- benchmark_azimuth(seurat_object, n_iterations)
  }

  # === scCATCH ===
  if ("sccatch" %in% methods) {
    message("\n=== Benchmarking scCATCH ===")
    results$sccatch <- benchmark_sccatch(seurat_object, tissue_type, n_iterations)
  }

  # === Garnett ===
  if ("garnett" %in% methods) {
    message("\n=== Benchmarking Garnett ===")
    results$garnett <- benchmark_garnett(seurat_object, n_iterations)
  }

  # Compile results
  benchmark_df <- do.call(rbind, lapply(names(results), function(method) {
    res <- results[[method]]
    data.frame(
      method = method,
      mean_time_sec = res$mean_time,
      sd_time_sec = res$sd_time,
      peak_memory_mb = res$peak_memory,
      accuracy = res$accuracy,
      f1_score = res$f1_score,
      n_predicted_types = res$n_predicted_types,
      installation_difficulty = res$installation_difficulty,
      custom_markers_supported = res$custom_markers,
      stringsAsFactors = FALSE
    )
  }))

  # Add dataset info
  benchmark_df$n_cells <- n_cells
  benchmark_df$n_genes <- n_genes
  benchmark_df$n_true_types <- length(unique(true_labels))

  return(benchmark_df)
}


#' Benchmark ScType
#' @keywords internal
benchmark_sctype <- function(seurat_object, tissue_type, n_iterations) {

  if (!requireNamespace("ScType", quietly = TRUE)) {
    # Try sourcing from GitHub if not installed
    tryCatch({
      source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R")
    }, error = function(e) {
      stop("ScType not available")
    })
  }

  times <- numeric(n_iterations)
  predictions <- NULL
  peak_mem <- 0

  for (i in 1:n_iterations) {
    gc()  # Clean memory
    mem_before <- as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE))

    start_time <- Sys.time()
    seurat_annotated <- tryCatch({
      run_sctype(seurat_object, known_tissue_type = tissue_type,
                name = "sctype_pred", plot = FALSE)
    }, error = function(e) {
      message("ScType error: ", e$message)
      return(NULL)
    })
    end_time <- Sys.time()

    mem_after <- as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE))
    mem_used <- (mem_before - mem_after) / 1024  # Convert to MB
    peak_mem <- max(peak_mem, mem_used)

    times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))

    if (!is.null(seurat_annotated)) {
      predictions <- seurat_annotated@meta.data$sctype_pred
    }
  }

  # Calculate metrics
  if (is.null(predictions)) {
    return(list(mean_time = NA, sd_time = NA, peak_memory = NA,
                accuracy = NA, f1_score = NA, n_predicted_types = NA,
                installation_difficulty = "Easy (R package)",
                custom_markers = TRUE))
  }

  true_labels <- seurat_object@meta.data[[1]]  # Adjust as needed
  metrics <- calculate_accuracy_metrics(true_labels, predictions)

  list(
    mean_time = mean(times),
    sd_time = sd(times),
    peak_memory = peak_mem,
    accuracy = metrics$accuracy,
    f1_score = metrics$f1_weighted,
    n_predicted_types = length(unique(predictions)),
    installation_difficulty = "Easy (R package)",
    custom_markers = TRUE
  )
}


#' Benchmark SingleR
#' @keywords internal
benchmark_singler <- function(seurat_object, n_iterations) {

  if (!requireNamespace("SingleR", quietly = TRUE)) {
    message("SingleR not installed. Install with: BiocManager::install('SingleR')")
    return(list(mean_time = NA, sd_time = NA, peak_memory = NA,
                accuracy = NA, f1_score = NA, n_predicted_types = NA,
                installation_difficulty = "Medium (Bioconductor)",
                custom_markers = FALSE))
  }

  if (!requireNamespace("celldex", quietly = TRUE)) {
    message("celldex not installed (needed for references)")
    return(list(mean_time = NA, sd_time = NA, peak_memory = NA,
                accuracy = NA, f1_score = NA, n_predicted_types = NA,
                installation_difficulty = "Medium (Bioconductor)",
                custom_markers = FALSE))
  }

  times <- numeric(n_iterations)
  predictions <- NULL
  peak_mem <- 0

  # Get reference dataset
  ref <- tryCatch({
    celldex::HumanPrimaryCellAtlasData()
  }, error = function(e) {
    message("Could not load reference data")
    return(NULL)
  })

  if (is.null(ref)) {
    return(list(mean_time = NA, sd_time = NA, peak_memory = NA,
                accuracy = NA, f1_score = NA, n_predicted_types = NA,
                installation_difficulty = "Medium (Bioconductor)",
                custom_markers = FALSE))
  }

  # Convert Seurat to SingleCellExperiment
  sce <- Seurat::as.SingleCellExperiment(seurat_object)

  for (i in 1:n_iterations) {
    gc()
    mem_before <- as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE))

    start_time <- Sys.time()
    pred <- tryCatch({
      SingleR::SingleR(test = sce, ref = ref, labels = ref$label.main)
    }, error = function(e) {
      message("SingleR error: ", e$message)
      return(NULL)
    })
    end_time <- Sys.time()

    mem_after <- as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE))
    mem_used <- (mem_before - mem_after) / 1024
    peak_mem <- max(peak_mem, mem_used)

    times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))

    if (!is.null(pred)) {
      predictions <- pred$labels
    }
  }

  if (is.null(predictions)) {
    return(list(mean_time = NA, sd_time = NA, peak_memory = NA,
                accuracy = NA, f1_score = NA, n_predicted_types = NA,
                installation_difficulty = "Medium (Bioconductor)",
                custom_markers = FALSE))
  }

  true_labels <- seurat_object@meta.data[[1]]
  metrics <- calculate_accuracy_metrics(true_labels, predictions)

  list(
    mean_time = mean(times),
    sd_time = sd(times),
    peak_memory = peak_mem,
    accuracy = metrics$accuracy,
    f1_score = metrics$f1_weighted,
    n_predicted_types = length(unique(predictions)),
    installation_difficulty = "Medium (Bioconductor)",
    custom_markers = FALSE
  )
}


#' Benchmark CellTypist (Python-based)
#' @keywords internal
benchmark_celltypist <- function(seurat_object, n_iterations) {

  # CellTypist requires Python via reticulate
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    message("reticulate not installed (needed for CellTypist)")
    return(list(mean_time = NA, sd_time = NA, peak_memory = NA,
                accuracy = NA, f1_score = NA, n_predicted_types = NA,
                installation_difficulty = "Hard (Python + R)",
                custom_markers = TRUE))
  }

  # Check if celltypist is available in Python
  has_celltypist <- tryCatch({
    reticulate::py_module_available("celltypist")
  }, error = function(e) FALSE)

  if (!has_celltypist) {
    message("CellTypist not available in Python. Install with: pip install celltypist")
    return(list(mean_time = NA, sd_time = NA, peak_memory = NA,
                accuracy = NA, f1_score = NA, n_predicted_types = NA,
                installation_difficulty = "Hard (Python + R)",
                custom_markers = TRUE))
  }

  times <- numeric(n_iterations)
  predictions <- NULL
  peak_mem <- 0

  # Convert to AnnData format for Python
  adata <- tryCatch({
    SeuratDisk::SaveH5Seurat(seurat_object, filename = "/tmp/temp_seurat.h5seurat", overwrite = TRUE)
    SeuratDisk::Convert("/tmp/temp_seurat.h5seurat", dest = "h5ad", overwrite = TRUE)

    ct <- reticulate::import("celltypist")
    sc <- reticulate::import("scanpy")
    sc$read_h5ad("/tmp/temp_seurat.h5ad")
  }, error = function(e) {
    message("Error converting to AnnData: ", e$message)
    return(NULL)
  })

  if (is.null(adata)) {
    return(list(mean_time = NA, sd_time = NA, peak_memory = NA,
                accuracy = NA, f1_score = NA, n_predicted_types = NA,
                installation_difficulty = "Hard (Python + R)",
                custom_markers = TRUE))
  }

  for (i in 1:n_iterations) {
    gc()
    mem_before <- as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE))

    start_time <- Sys.time()
    result <- tryCatch({
      ct <- reticulate::import("celltypist")
      model <- ct$models$Model$load(model = "Immune_All_Low.pkl")
      predictions_obj <- ct$annotate(adata, model = model, majority_voting = TRUE)
      predictions_obj$predicted_labels$majority_voting
    }, error = function(e) {
      message("CellTypist error: ", e$message)
      return(NULL)
    })
    end_time <- Sys.time()

    mem_after <- as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE))
    mem_used <- (mem_before - mem_after) / 1024
    peak_mem <- max(peak_mem, mem_used)

    times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
    predictions <- result
  }

  if (is.null(predictions)) {
    return(list(mean_time = NA, sd_time = NA, peak_memory = NA,
                accuracy = NA, f1_score = NA, n_predicted_types = NA,
                installation_difficulty = "Hard (Python + R)",
                custom_markers = TRUE))
  }

  true_labels <- seurat_object@meta.data[[1]]
  metrics <- calculate_accuracy_metrics(true_labels, predictions)

  list(
    mean_time = mean(times),
    sd_time = sd(times),
    peak_memory = peak_mem,
    accuracy = metrics$accuracy,
    f1_score = metrics$f1_weighted,
    n_predicted_types = length(unique(predictions)),
    installation_difficulty = "Hard (Python + R)",
    custom_markers = TRUE
  )
}


#' Benchmark Azimuth
#' @keywords internal
benchmark_azimuth <- function(seurat_object, n_iterations) {

  if (!requireNamespace("Azimuth", quietly = TRUE)) {
    message("Azimuth not installed. Install with: remotes::install_github('satijalab/azimuth')")
    return(list(mean_time = NA, sd_time = NA, peak_memory = NA,
                accuracy = NA, f1_score = NA, n_predicted_types = NA,
                installation_difficulty = "Medium (GitHub)",
                custom_markers = FALSE))
  }

  times <- numeric(n_iterations)
  predictions <- NULL
  peak_mem <- 0

  # Azimuth uses web API or local reference
  for (i in 1:n_iterations) {
    gc()
    mem_before <- as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE))

    start_time <- Sys.time()
    result <- tryCatch({
      # Using PBMC reference as example
      Azimuth::RunAzimuth(seurat_object, reference = "pbmcref")
    }, error = function(e) {
      message("Azimuth error: ", e$message)
      return(NULL)
    })
    end_time <- Sys.time()

    mem_after <- as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE))
    mem_used <- (mem_before - mem_after) / 1024
    peak_mem <- max(peak_mem, mem_used)

    times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))

    if (!is.null(result)) {
      predictions <- result@meta.data$predicted.celltype.l2
    }
  }

  if (is.null(predictions)) {
    return(list(mean_time = NA, sd_time = NA, peak_memory = NA,
                accuracy = NA, f1_score = NA, n_predicted_types = NA,
                installation_difficulty = "Medium (GitHub)",
                custom_markers = FALSE))
  }

  true_labels <- seurat_object@meta.data[[1]]
  metrics <- calculate_accuracy_metrics(true_labels, predictions)

  list(
    mean_time = mean(times),
    sd_time = sd(times),
    peak_memory = peak_mem,
    accuracy = metrics$accuracy,
    f1_score = metrics$f1_weighted,
    n_predicted_types = length(unique(predictions)),
    installation_difficulty = "Medium (GitHub)",
    custom_markers = FALSE
  )
}


#' Benchmark scCATCH
#' @keywords internal
benchmark_sccatch <- function(seurat_object, tissue_type, n_iterations) {

  if (!requireNamespace("scCATCH", quietly = TRUE)) {
    message("scCATCH not installed. Install with: BiocManager::install('scCATCH')")
    return(list(mean_time = NA, sd_time = NA, peak_memory = NA,
                accuracy = NA, f1_score = NA, n_predicted_types = NA,
                installation_difficulty = "Medium (Bioconductor)",
                custom_markers = TRUE))
  }

  times <- numeric(n_iterations)
  predictions <- NULL
  peak_mem <- 0

  for (i in 1:n_iterations) {
    gc()
    mem_before <- as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE))

    start_time <- Sys.time()
    result <- tryCatch({
      obj <- scCATCH::createscCATCH(data = Seurat::GetAssayData(seurat_object, slot = "data"),
                                    cluster = as.character(Idents(seurat_object)))
      obj <- scCATCH::findmarkergene(object = obj, species = "Human",
                                     tissue = tissue_type)
      obj <- scCATCH::findcelltype(object = obj)
      obj
    }, error = function(e) {
      message("scCATCH error: ", e$message)
      return(NULL)
    })
    end_time <- Sys.time()

    mem_after <- as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE))
    mem_used <- (mem_before - mem_after) / 1024
    peak_mem <- max(peak_mem, mem_used)

    times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))

    if (!is.null(result)) {
      # Map cluster annotations to cells
      cluster_celltype <- result@celltype
      cell_clusters <- as.character(Idents(seurat_object))
      predictions <- cluster_celltype$cell_type[match(cell_clusters, cluster_celltype$cluster)]
    }
  }

  if (is.null(predictions)) {
    return(list(mean_time = NA, sd_time = NA, peak_memory = NA,
                accuracy = NA, f1_score = NA, n_predicted_types = NA,
                installation_difficulty = "Medium (Bioconductor)",
                custom_markers = TRUE))
  }

  true_labels <- seurat_object@meta.data[[1]]
  metrics <- calculate_accuracy_metrics(true_labels, predictions)

  list(
    mean_time = mean(times),
    sd_time = sd(times),
    peak_memory = peak_mem,
    accuracy = metrics$accuracy,
    f1_score = metrics$f1_weighted,
    n_predicted_types = length(unique(predictions)),
    installation_difficulty = "Medium (Bioconductor)",
    custom_markers = TRUE
  )
}


#' Benchmark Garnett
#' @keywords internal
benchmark_garnett <- function(seurat_object, n_iterations) {

  if (!requireNamespace("garnett", quietly = TRUE)) {
    message("Garnett not installed. Install with: BiocManager::install('garnett')")
    return(list(mean_time = NA, sd_time = NA, peak_memory = NA,
                accuracy = NA, f1_score = NA, n_predicted_types = NA,
                installation_difficulty = "Hard (Bioconductor + marker file)",
                custom_markers = TRUE))
  }

  # Garnett requires a marker file in specific format - complex to benchmark
  message("Garnett requires custom marker file configuration - skipping automated benchmark")

  return(list(mean_time = NA, sd_time = NA, peak_memory = NA,
              accuracy = NA, f1_score = NA, n_predicted_types = NA,
              installation_difficulty = "Hard (Bioconductor + marker file)",
              custom_markers = TRUE))
}


#' Calculate accuracy metrics
#' @keywords internal
calculate_accuracy_metrics <- function(true_labels, predicted_labels) {

  # Remove any NA or Unknown predictions
  valid_idx <- !is.na(predicted_labels) & predicted_labels != "Unknown" &
               !is.na(true_labels) & true_labels != "Unknown"

  true_labels <- true_labels[valid_idx]
  predicted_labels <- predicted_labels[valid_idx]

  if (length(true_labels) == 0) {
    return(list(accuracy = 0, precision = 0, recall = 0,
                f1_weighted = 0, confusion = NULL))
  }

  # Simple accuracy (exact match)
  accuracy <- sum(true_labels == predicted_labels) / length(true_labels)

  # Weighted F1 score (handles multi-class)
  unique_labels <- unique(c(true_labels, predicted_labels))
  f1_scores <- numeric(length(unique_labels))
  weights <- numeric(length(unique_labels))

  for (i in seq_along(unique_labels)) {
    label <- unique_labels[i]
    tp <- sum(true_labels == label & predicted_labels == label)
    fp <- sum(true_labels != label & predicted_labels == label)
    fn <- sum(true_labels == label & predicted_labels != label)

    precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
    recall <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
    f1_scores[i] <- ifelse(precision + recall > 0,
                           2 * precision * recall / (precision + recall), 0)
    weights[i] <- sum(true_labels == label)
  }

  f1_weighted <- sum(f1_scores * weights) / sum(weights)

  # Confusion matrix
  confusion <- table(True = true_labels, Predicted = predicted_labels)

  list(
    accuracy = accuracy,
    f1_weighted = f1_weighted,
    confusion = confusion
  )
}


#' Plot benchmark results
#' @param benchmark_df Data frame from benchmark_cell_type_annotation
#' @export
plot_benchmark_results <- function(benchmark_df) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required for plotting")
  }

  library(ggplot2)

  # Remove NA rows
  benchmark_df <- benchmark_df[!is.na(benchmark_df$mean_time_sec), ]

  plots <- list()

  # 1. Speed comparison
  p1 <- ggplot(benchmark_df, aes(x = reorder(method, mean_time_sec),
                                 y = mean_time_sec, fill = method)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = mean_time_sec - sd_time_sec,
                      ymax = mean_time_sec + sd_time_sec), width = 0.2) +
    coord_flip() +
    labs(title = "Runtime Comparison",
         x = "Method", y = "Time (seconds)") +
    theme_minimal() +
    theme(legend.position = "none")

  plots$runtime <- p1

  # 2. Accuracy comparison
  p2 <- ggplot(benchmark_df, aes(x = reorder(method, accuracy),
                                 y = accuracy, fill = method)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    ylim(0, 1) +
    labs(title = "Accuracy Comparison",
         x = "Method", y = "Accuracy") +
    theme_minimal() +
    theme(legend.position = "none")

  plots$accuracy <- p2

  # 3. F1 Score comparison
  p3 <- ggplot(benchmark_df, aes(x = reorder(method, f1_score),
                                 y = f1_score, fill = method)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    ylim(0, 1) +
    labs(title = "F1 Score Comparison",
         x = "Method", y = "Weighted F1 Score") +
    theme_minimal() +
    theme(legend.position = "none")

  plots$f1_score <- p3

  # 4. Memory usage comparison
  p4 <- ggplot(benchmark_df, aes(x = reorder(method, peak_memory_mb),
                                 y = peak_memory_mb, fill = method)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "Peak Memory Usage",
         x = "Method", y = "Memory (MB)") +
    theme_minimal() +
    theme(legend.position = "none")

  plots$memory <- p4

  # 5. Speed vs Accuracy tradeoff
  p5 <- ggplot(benchmark_df, aes(x = mean_time_sec, y = accuracy,
                                 color = method, size = 3)) +
    geom_point() +
    geom_text(aes(label = method), vjust = -1, hjust = 0.5, size = 3) +
    labs(title = "Speed vs Accuracy Tradeoff",
         x = "Runtime (seconds)", y = "Accuracy") +
    theme_minimal() +
    theme(legend.position = "none")

  plots$tradeoff <- p5

  return(plots)
}


#' Save benchmark results
#' @param benchmark_df Data frame from benchmark_cell_type_annotation
#' @param output_dir Output directory for results
#' @export
save_benchmark_results <- function(benchmark_df, output_dir = "benchmark_results") {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save CSV
  write.csv(benchmark_df, file.path(output_dir, "benchmark_results.csv"),
            row.names = FALSE)

  # Save plots
  plots <- plot_benchmark_results(benchmark_df)

  for (plot_name in names(plots)) {
    ggplot2::ggsave(
      filename = file.path(output_dir, paste0(plot_name, ".png")),
      plot = plots[[plot_name]],
      width = 10, height = 6, dpi = 300
    )
  }

  # Generate markdown report
  report <- generate_benchmark_report(benchmark_df)
  writeLines(report, file.path(output_dir, "benchmark_report.md"))

  message(sprintf("Results saved to: %s", output_dir))
}


#' Generate markdown benchmark report
#' @keywords internal
generate_benchmark_report <- function(benchmark_df) {

  lines <- c(
    "# Cell Type Annotation Tools Benchmark",
    "",
    sprintf("**Date:** %s", Sys.Date()),
    sprintf("**Dataset:** %d cells, %d genes, %d cell types",
            benchmark_df$n_cells[1], benchmark_df$n_genes[1],
            benchmark_df$n_true_types[1]),
    "",
    "## Summary",
    "",
    "### Performance Metrics",
    "",
    "| Method | Runtime (s) | Memory (MB) | Accuracy | F1 Score | Installation | Custom Markers |",
    "|--------|-------------|-------------|----------|----------|--------------|----------------|"
  )

  for (i in 1:nrow(benchmark_df)) {
    row <- benchmark_df[i, ]
    line <- sprintf(
      "| %s | %.2f ± %.2f | %.1f | %.3f | %.3f | %s | %s |",
      row$method,
      row$mean_time_sec,
      row$sd_time_sec,
      row$peak_memory_mb,
      row$accuracy,
      row$f1_score,
      row$installation_difficulty,
      ifelse(row$custom_markers_supported, "Yes", "No")
    )
    lines <- c(lines, line)
  }

  lines <- c(lines,
    "",
    "## Key Findings",
    "",
    "### Fastest Method",
    sprintf("**%s** - %.2f seconds",
            benchmark_df$method[which.min(benchmark_df$mean_time_sec)],
            min(benchmark_df$mean_time_sec, na.rm = TRUE)),
    "",
    "### Most Accurate Method",
    sprintf("**%s** - %.1f%% accuracy",
            benchmark_df$method[which.max(benchmark_df$accuracy)],
            max(benchmark_df$accuracy, na.rm = TRUE) * 100),
    "",
    "### Best F1 Score",
    sprintf("**%s** - %.3f weighted F1",
            benchmark_df$method[which.max(benchmark_df$f1_score)],
            max(benchmark_df$f1_score, na.rm = TRUE)),
    "",
    "### Lowest Memory Usage",
    sprintf("**%s** - %.1f MB peak memory",
            benchmark_df$method[which.min(benchmark_df$peak_memory_mb)],
            min(benchmark_df$peak_memory_mb, na.rm = TRUE)),
    "",
    "## Conclusions",
    "",
    "- **For speed**: Choose the fastest method if time is critical",
    "- **For accuracy**: Choose the most accurate method for best results",
    "- **For ease of use**: Consider installation difficulty and documentation",
    "- **For flexibility**: Choose methods supporting custom markers",
    "",
    "## Plots",
    "",
    "See generated PNG files in this directory for visual comparisons.",
    ""
  )

  paste(lines, collapse = "\n")
}
