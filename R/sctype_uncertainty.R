# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/master/LICENSE)
# ScType uncertainty scoring and confidence metrics for Seurat objects
# Written by Claude (AI assistant), 2025-11-15

#' Add uncertainty scores to ScType annotations (Seurat)
#'
#' Calculates confidence/uncertainty metrics for cell type annotations and adds
#' them to the Seurat object metadata. Provides top N candidate cell types per
#' cluster with their scores, confidence metrics, and uncertainty quantification.
#'
#' @param seurat_object Seurat object with clustering
#' @param known_tissue_type Tissue type (e.g., "Immune system", "Brain")
#' @param database_file Path to marker database (default: ScTypeDB_full.xlsx from GitHub)
#' @param assay Assay to use (default: "RNA")
#' @param scaled Use scaled data (default: TRUE)
#' @param cluster_col Column with cluster assignments (default: "seurat_clusters")
#' @param top_n Number of top candidate cell types to report (default: 3)
#' @param annotation_prefix Prefix for new metadata columns (default: "sctype")
#'
#' @return Modified Seurat object with uncertainty metrics in metadata
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- add_sctype_uncertainty(seurat_obj, known_tissue_type = "Immune system")
#' # New columns added:
#' # - sctype_top1, sctype_top2, sctype_top3 (cell type names)
#' # - sctype_score1, sctype_score2, sctype_score3 (raw scores)
#' # - sctype_confidence (0-1 normalized confidence)
#' # - sctype_confidence_level (High/Medium/Low)
#' # - sctype_score_diff (top1 - top2)
#' }
add_sctype_uncertainty <- function(seurat_object,
                                   known_tissue_type = NULL,
                                   database_file = NULL,
                                   assay = "RNA",
                                   scaled = TRUE,
                                   cluster_col = "seurat_clusters",
                                   top_n = 3,
                                   annotation_prefix = "sctype") {

  # Load required packages
  lapply(c("dplyr", "Seurat", "HGNChelper", "openxlsx"), library, character.only = TRUE)

  # Source ScType functions
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

  # Set default database
  if (is.null(database_file)) {
    database_file <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
  }

  # Auto-detect tissue type if not provided
  if (is.null(known_tissue_type)) {
    stop("Please specify known_tissue_type parameter (e.g., 'Immune system', 'Brain', etc.)")
  }

  message(paste0("Calculating uncertainty scores for tissue type: ", known_tissue_type))

  # Check cluster column
  if (!cluster_col %in% colnames(seurat_object@meta.data)) {
    stop(paste0("Cluster column '", cluster_col, "' not found in metadata."))
  }

  # Prepare gene sets
  gs_list <- gene_sets_prepare(database_file, known_tissue_type)

  # Extract scaled data
  DefaultAssay(seurat_object) <- assay
  data_type <- if (scaled) "scale.data" else "counts"

  # Check Seurat version and extract data
  seurat_v5 <- isFALSE('counts' %in% names(attributes(seurat_object[[assay]])))
  if (seurat_v5) {
    scRNAseqData <- as.matrix(slot(seurat_object[[assay]], data_type))
  } else {
    scRNAseqData <- as.matrix(seurat_object[[assay]]@scale.data)
  }

  # Run ScType scoring
  message("Computing ScType scores...")
  es.max <- sctype_score(scRNAseqData = scRNAseqData,
                        scaled = scaled,
                        gs = gs_list$gs_positive,
                        gs2 = gs_list$gs_negative)

  # Get cluster assignments
  clusters <- seurat_object@meta.data[[cluster_col]]

  # Calculate top N candidates and scores per cluster
  message(paste0("Calculating top ", top_n, " candidates per cluster..."))

  cluster_results <- do.call("rbind", lapply(unique(clusters), function(cl) {
    cells_in_cluster <- which(clusters == cl)
    es_cl <- sort(rowSums(es.max[, cells_in_cluster, drop = FALSE]), decreasing = TRUE)

    # Get top N
    top_indices <- 1:min(top_n, length(es_cl))
    top_types <- names(es_cl)[top_indices]
    top_scores <- as.numeric(es_cl[top_indices])

    # Pad with NA if fewer than top_n cell types
    if (length(top_types) < top_n) {
      top_types <- c(top_types, rep(NA, top_n - length(top_types)))
      top_scores <- c(top_scores, rep(NA, top_n - length(top_scores)))
    }

    # Calculate uncertainty metrics
    score_diff <- ifelse(length(top_scores) >= 2 && !is.na(top_scores[2]),
                         top_scores[1] - top_scores[2], NA)

    # Normalized confidence (0-1): based on score difference and top score
    # Higher score_diff = more confident
    # Use sigmoid-like transformation
    ncells <- length(cells_in_cluster)
    confidence_threshold <- ncells / 4  # Same as ScType's Unknown threshold

    if (!is.na(top_scores[1]) && top_scores[1] > 0) {
      # Confidence based on: (1) absolute score, (2) score difference
      abs_confidence <- min(top_scores[1] / confidence_threshold, 1)  # 0-1
      rel_confidence <- ifelse(!is.na(score_diff) && score_diff > 0,
                              min(score_diff / (top_scores[1] * 0.5), 1), 0)
      confidence_score <- (abs_confidence + rel_confidence) / 2
    } else {
      confidence_score <- 0
    }

    # Confidence level
    confidence_level <- ifelse(is.na(confidence_score), "Unknown",
                              ifelse(confidence_score >= 0.7, "High",
                                    ifelse(confidence_score >= 0.4, "Medium", "Low")))

    # Apply Unknown threshold
    if (!is.na(top_scores[1]) && top_scores[1] < confidence_threshold) {
      top_types[1] <- "Unknown"
      confidence_level <- "Low"
    }

    data.frame(
      cluster = cl,
      ncells = ncells,
      top1 = top_types[1],
      top2 = if(top_n >= 2) top_types[2] else NA,
      top3 = if(top_n >= 3) top_types[3] else NA,
      score1 = top_scores[1],
      score2 = if(top_n >= 2) top_scores[2] else NA,
      score3 = if(top_n >= 3) top_scores[3] else NA,
      score_diff = score_diff,
      confidence = confidence_score,
      confidence_level = confidence_level,
      stringsAsFactors = FALSE
    )
  }))

  # Assign to individual cells
  message("Adding uncertainty metrics to metadata...")

  # Initialize new columns
  n_cells <- ncol(seurat_object)
  seurat_object@meta.data[[paste0(annotation_prefix, "_top1")]] <- rep(NA, n_cells)
  seurat_object@meta.data[[paste0(annotation_prefix, "_score1")]] <- rep(NA, n_cells)
  seurat_object@meta.data[[paste0(annotation_prefix, "_confidence")]] <- rep(NA, n_cells)
  seurat_object@meta.data[[paste0(annotation_prefix, "_confidence_level")]] <- rep(NA, n_cells)
  seurat_object@meta.data[[paste0(annotation_prefix, "_score_diff")]] <- rep(NA, n_cells)

  if (top_n >= 2) {
    seurat_object@meta.data[[paste0(annotation_prefix, "_top2")]] <- rep(NA, n_cells)
    seurat_object@meta.data[[paste0(annotation_prefix, "_score2")]] <- rep(NA, n_cells)
  }
  if (top_n >= 3) {
    seurat_object@meta.data[[paste0(annotation_prefix, "_top3")]] <- rep(NA, n_cells)
    seurat_object@meta.data[[paste0(annotation_prefix, "_score3")]] <- rep(NA, n_cells)
  }

  # Fill in values
  for (i in 1:nrow(cluster_results)) {
    cl <- cluster_results$cluster[i]
    cells_in_cluster <- which(clusters == cl)

    seurat_object@meta.data[[paste0(annotation_prefix, "_top1")]][cells_in_cluster] <-
      cluster_results$top1[i]
    seurat_object@meta.data[[paste0(annotation_prefix, "_score1")]][cells_in_cluster] <-
      cluster_results$score1[i]
    seurat_object@meta.data[[paste0(annotation_prefix, "_confidence")]][cells_in_cluster] <-
      cluster_results$confidence[i]
    seurat_object@meta.data[[paste0(annotation_prefix, "_confidence_level")]][cells_in_cluster] <-
      cluster_results$confidence_level[i]
    seurat_object@meta.data[[paste0(annotation_prefix, "_score_diff")]][cells_in_cluster] <-
      cluster_results$score_diff[i]

    if (top_n >= 2) {
      seurat_object@meta.data[[paste0(annotation_prefix, "_top2")]][cells_in_cluster] <-
        cluster_results$top2[i]
      seurat_object@meta.data[[paste0(annotation_prefix, "_score2")]][cells_in_cluster] <-
        cluster_results$score2[i]
    }
    if (top_n >= 3) {
      seurat_object@meta.data[[paste0(annotation_prefix, "_top3")]][cells_in_cluster] <-
        cluster_results$top3[i]
      seurat_object@meta.data[[paste0(annotation_prefix, "_score3")]][cells_in_cluster] <-
        cluster_results$score3[i]
    }
  }

  # Store cluster-level results as an attribute
  attr(seurat_object, "sctype_uncertainty_clusters") <- cluster_results

  message("Done! Added the following columns to metadata:")
  message(paste0("  - ", annotation_prefix, "_top1: Primary cell type assignment"))
  message(paste0("  - ", annotation_prefix, "_score1: Score for primary assignment"))
  if (top_n >= 2) {
    message(paste0("  - ", annotation_prefix, "_top2: Second candidate cell type"))
    message(paste0("  - ", annotation_prefix, "_score2: Score for second candidate"))
  }
  if (top_n >= 3) {
    message(paste0("  - ", annotation_prefix, "_top3: Third candidate cell type"))
    message(paste0("  - ", annotation_prefix, "_score3: Score for third candidate"))
  }
  message(paste0("  - ", annotation_prefix, "_confidence: Confidence score (0-1)"))
  message(paste0("  - ", annotation_prefix, "_confidence_level: High/Medium/Low"))
  message(paste0("  - ", annotation_prefix, "_score_diff: Score difference (top1 - top2)"))

  # Print summary
  message("\nConfidence level distribution:")
  print(table(seurat_object@meta.data[[paste0(annotation_prefix, "_confidence_level")]]))

  return(seurat_object)
}


#' Visualize ScType uncertainty scores
#'
#' Creates comprehensive visualizations of uncertainty metrics including:
#' 1. Bar plots of top N candidates per cluster with scores
#' 2. UMAP colored by confidence score
#' 3. Distribution of confidence levels
#' 4. Uncertainty heatmap
#'
#' @param seurat_object Seurat object with uncertainty scores (from add_sctype_uncertainty)
#' @param annotation_prefix Prefix used in add_sctype_uncertainty (default: "sctype")
#' @param cluster_col Cluster column (default: "seurat_clusters")
#' @param plot_types Vector of plot types: "candidates", "umap", "distribution", "heatmap" (default: all)
#' @param save_plots Save plots to files (default: FALSE)
#' @param output_dir Directory for saved plots (default: "uncertainty_plots")
#'
#' @return List of ggplot objects
#' @export
visualize_sctype_uncertainty <- function(seurat_object,
                                        annotation_prefix = "sctype",
                                        cluster_col = "seurat_clusters",
                                        plot_types = c("candidates", "umap", "distribution", "heatmap"),
                                        save_plots = FALSE,
                                        output_dir = "uncertainty_plots") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for visualization. Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' required. Install with: install.packages('dplyr')")
  }

  # Check if uncertainty columns exist
  required_cols <- c(paste0(annotation_prefix, "_top1"),
                    paste0(annotation_prefix, "_confidence"))
  if (!all(required_cols %in% colnames(seurat_object@meta.data))) {
    stop("Uncertainty scores not found. Run add_sctype_uncertainty() first.")
  }

  # Get cluster-level results
  if (is.null(attr(seurat_object, "sctype_uncertainty_clusters"))) {
    stop("Cluster-level uncertainty data not found. Run add_sctype_uncertainty() first.")
  }
  cluster_results <- attr(seurat_object, "sctype_uncertainty_clusters")

  # Create output directory
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  plot_list <- list()

  # 1. TOP CANDIDATES BAR PLOT
  if ("candidates" %in% plot_types) {
    message("Generating top candidates plot...")
    p <- plot_top_candidates(cluster_results, annotation_prefix)
    plot_list$candidates <- p

    if (save_plots) {
      ggsave(file.path(output_dir, "top_candidates.png"), p,
             width = 14, height = 10, dpi = 300)
    }
  }

  # 2. UMAP CONFIDENCE
  if ("umap" %in% plot_types) {
    message("Generating UMAP confidence plots...")
    umap_plots <- plot_umap_confidence(seurat_object, annotation_prefix)
    plot_list$umap <- umap_plots

    if (save_plots) {
      if (requireNamespace("patchwork", quietly = TRUE)) {
        combined <- patchwork::wrap_plots(umap_plots, ncol = 2)
        ggsave(file.path(output_dir, "umap_confidence.png"), combined,
               width = 16, height = 8, dpi = 300)
      }
    }
  }

  # 3. CONFIDENCE DISTRIBUTION
  if ("distribution" %in% plot_types) {
    message("Generating confidence distribution plots...")
    dist_plots <- plot_confidence_distribution(seurat_object, cluster_results,
                                               annotation_prefix, cluster_col)
    plot_list$distribution <- dist_plots

    if (save_plots) {
      if (requireNamespace("patchwork", quietly = TRUE)) {
        combined <- patchwork::wrap_plots(dist_plots, ncol = 2)
        ggsave(file.path(output_dir, "confidence_distribution.png"), combined,
               width = 14, height = 10, dpi = 300)
      }
    }
  }

  # 4. UNCERTAINTY HEATMAP
  if ("heatmap" %in% plot_types) {
    message("Generating uncertainty heatmap...")
    p <- plot_uncertainty_heatmap(cluster_results)
    plot_list$heatmap <- p

    if (save_plots) {
      ggsave(file.path(output_dir, "uncertainty_heatmap.png"), p,
             width = 12, height = 8, dpi = 300)
    }
  }

  message("Visualization complete!")
  if (save_plots) {
    message(paste0("Plots saved to: ", output_dir))
  }

  invisible(plot_list)
}


#' Plot top candidate cell types with scores
#' @keywords internal
plot_top_candidates <- function(cluster_results, annotation_prefix) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for visualization. Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' required. Install with: install.packages('dplyr')")
  }

  # Prepare data for plotting
  plot_data <- list()

  for (i in 1:nrow(cluster_results)) {
    cluster_id <- cluster_results$cluster[i]

    # Top 1
    if (!is.na(cluster_results$top1[i])) {
      plot_data[[length(plot_data) + 1]] <- data.frame(
        cluster = cluster_id,
        rank = "1st",
        celltype = cluster_results$top1[i],
        score = cluster_results$score1[i],
        confidence = cluster_results$confidence[i]
      )
    }

    # Top 2
    if ("top2" %in% colnames(cluster_results) && !is.na(cluster_results$top2[i])) {
      plot_data[[length(plot_data) + 1]] <- data.frame(
        cluster = cluster_id,
        rank = "2nd",
        celltype = cluster_results$top2[i],
        score = cluster_results$score2[i],
        confidence = cluster_results$confidence[i]
      )
    }

    # Top 3
    if ("top3" %in% colnames(cluster_results) && !is.na(cluster_results$top3[i])) {
      plot_data[[length(plot_data) + 1]] <- data.frame(
        cluster = cluster_id,
        rank = "3rd",
        celltype = cluster_results$top3[i],
        score = cluster_results$score3[i],
        confidence = cluster_results$confidence[i]
      )
    }
  }

  plot_df <- do.call("rbind", plot_data)
  plot_df$rank <- factor(plot_df$rank, levels = c("1st", "2nd", "3rd"))

  # Create plot
  p <- ggplot(plot_df, aes(x = factor(cluster), y = score, fill = rank)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.3) +
    geom_text(aes(label = celltype), position = position_dodge(width = 0.9),
              angle = 45, hjust = -0.1, vjust = 0.5, size = 3) +
    scale_fill_manual(values = c("1st" = "#E41A1C", "2nd" = "#377EB8", "3rd" = "#4DAF4A")) +
    labs(title = "Top Cell Type Candidates per Cluster",
         subtitle = "Showing top 3 candidates with ScType scores",
         x = "Cluster", y = "ScType Score", fill = "Rank") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "right",
          plot.title = element_text(face = "bold", size = 14))

  return(p)
}


#' Plot UMAP colored by confidence
#' @keywords internal
plot_umap_confidence <- function(seurat_object, annotation_prefix) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for visualization. Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' required. Install with: install.packages('Seurat')")
  }

  plots <- list()

  # Check if UMAP exists
  if (!"umap" %in% names(seurat_object@reductions)) {
    warning("UMAP not found. Skipping UMAP plots.")
    return(plots)
  }

  # Plot 1: Confidence score
  p1 <- FeaturePlot(seurat_object,
                   features = paste0(annotation_prefix, "_confidence"),
                   reduction = "umap") +
    scale_color_gradient2(low = "blue", mid = "yellow", high = "red",
                         midpoint = 0.5) +
    ggtitle("Annotation Confidence Score")

  # Plot 2: Confidence level
  p2 <- DimPlot(seurat_object,
               group.by = paste0(annotation_prefix, "_confidence_level"),
               reduction = "umap") +
    ggtitle("Annotation Confidence Level")

  plots <- list(score = p1, level = p2)

  return(plots)
}


#' Plot confidence distribution
#' @keywords internal
plot_confidence_distribution <- function(seurat_object, cluster_results,
                                        annotation_prefix, cluster_col) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for visualization. Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' required. Install with: install.packages('dplyr')")
  }

  plots <- list()

  # Plot 1: Confidence level counts
  conf_data <- data.frame(
    level = seurat_object@meta.data[[paste0(annotation_prefix, "_confidence_level")]]
  )
  conf_data$level <- factor(conf_data$level, levels = c("High", "Medium", "Low", "Unknown"))

  p1 <- ggplot(conf_data, aes(x = level, fill = level)) +
    geom_bar() +
    scale_fill_manual(values = c("High" = "#2ECC40", "Medium" = "#FFDC00",
                                 "Low" = "#FF4136", "Unknown" = "#AAAAAA")) +
    labs(title = "Distribution of Confidence Levels",
         x = "Confidence Level", y = "Number of Cells") +
    theme_minimal() +
    theme(legend.position = "none")

  # Plot 2: Score difference distribution
  p2 <- ggplot(cluster_results, aes(x = factor(cluster), y = score_diff,
                                   fill = confidence_level)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("High" = "#2ECC40", "Medium" = "#FFDC00",
                                 "Low" = "#FF4136", "Unknown" = "#AAAAAA")) +
    labs(title = "Score Difference per Cluster (Top1 - Top2)",
         subtitle = "Larger difference = higher confidence",
         x = "Cluster", y = "Score Difference", fill = "Confidence") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  plots <- list(levels = p1, score_diff = p2)

  return(plots)
}


#' Plot uncertainty heatmap
#' @keywords internal
plot_uncertainty_heatmap <- function(cluster_results) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for visualization. Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' required. Install with: install.packages('dplyr')")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' required. Install with: install.packages('tidyr')")
  }

  # Prepare data
  heatmap_data <- cluster_results %>%
    select(cluster, top1, score1, confidence) %>%
    mutate(cluster = factor(cluster))

  # Create plot
  p <- ggplot(heatmap_data, aes(x = cluster, y = 1, fill = confidence)) +
    geom_tile(color = "white", size = 1) +
    geom_text(aes(label = top1), size = 3, fontface = "bold") +
    geom_text(aes(label = sprintf("%.2f", score1)), size = 2.5,
              vjust = 2.5, color = "gray30") +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green",
                        midpoint = 0.5, limits = c(0, 1)) +
    labs(title = "Cluster Annotation Confidence",
         subtitle = "Cell type (top) and score (bottom) shown in each tile",
         x = "Cluster", y = "", fill = "Confidence") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank())

  return(p)
}
