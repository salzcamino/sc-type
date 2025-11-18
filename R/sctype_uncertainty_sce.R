# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/master/LICENSE)
# ScType uncertainty scoring and confidence metrics for SingleCellExperiment objects
# Written by Claude (AI assistant), 2025-11-15

#' Add uncertainty scores to ScType annotations (SCE)
#'
#' Calculates confidence/uncertainty metrics for cell type annotations and adds
#' them to the SingleCellExperiment colData. Provides top N candidate cell types per
#' cluster with their scores, confidence metrics, and uncertainty quantification.
#'
#' @param sce_object SingleCellExperiment object with clustering
#' @param known_tissue_type Tissue type (e.g., "Immune system", "Brain")
#' @param database_file Path to marker database (default: ScTypeDB_full.xlsx from GitHub)
#' @param assay_name Assay to use (default: "logcounts")
#' @param scaled Use scaled data (default: TRUE)
#' @param cluster_col Column in colData with cluster assignments (default: "cluster")
#' @param top_n Number of top candidate cell types to report (default: 3)
#' @param annotation_prefix Prefix for new colData columns (default: "sctype")
#'
#' @return Modified SingleCellExperiment object with uncertainty metrics in colData
#' @export
#'
#' @examples
#' \dontrun{
#' sce <- add_sctype_uncertainty_sce(sce, known_tissue_type = "Immune system")
#' # New columns added to colData:
#' # - sctype_top1, sctype_top2, sctype_top3 (cell type names)
#' # - sctype_score1, sctype_score2, sctype_score3 (raw scores)
#' # - sctype_confidence (0-1 normalized confidence)
#' # - sctype_confidence_level (High/Medium/Low)
#' # - sctype_score_diff (top1 - top2)
#' }
add_sctype_uncertainty_sce <- function(sce_object,
                                       known_tissue_type = NULL,
                                       database_file = NULL,
                                       assay_name = "logcounts",
                                       scaled = TRUE,
                                       cluster_col = "cluster",
                                       top_n = 3,
                                       annotation_prefix = "sctype") {

  # Load required packages
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required but not installed.")
  }

  lapply(c("dplyr", "HGNChelper", "openxlsx", "SingleCellExperiment"), library, character.only = TRUE)

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
  if (!cluster_col %in% colnames(colData(sce_object))) {
    stop(paste0("Cluster column '", cluster_col, "' not found in colData."))
  }

  # Check assay
  if (!assay_name %in% SummarizedExperiment::assayNames(sce_object)) {
    stop(paste0("Assay '", assay_name, "' not found. Available assays: ",
                paste(SummarizedExperiment::assayNames(sce_object), collapse = ", ")))
  }

  # Prepare gene sets
  gs_list <- gene_sets_prepare(database_file, known_tissue_type)

  # Extract expression data
  scRNAseqData <- as.matrix(SummarizedExperiment::assay(sce_object, assay_name))

  # Run ScType scoring
  message("Computing ScType scores...")
  es.max <- sctype_score(scRNAseqData = scRNAseqData,
                        scaled = scaled,
                        gs = gs_list$gs_positive,
                        gs2 = gs_list$gs_negative)

  # Get cluster assignments
  clusters <- colData(sce_object)[[cluster_col]]

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

    # Normalized confidence (0-1)
    ncells <- length(cells_in_cluster)
    confidence_threshold <- ncells / 4

    if (!is.na(top_scores[1]) && top_scores[1] > 0) {
      abs_confidence <- min(top_scores[1] / confidence_threshold, 1)
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
  message("Adding uncertainty metrics to colData...")

  # Initialize new columns
  n_cells <- ncol(sce_object)
  colData(sce_object)[[paste0(annotation_prefix, "_top1")]] <- rep(NA, n_cells)
  colData(sce_object)[[paste0(annotation_prefix, "_score1")]] <- rep(NA, n_cells)
  colData(sce_object)[[paste0(annotation_prefix, "_confidence")]] <- rep(NA, n_cells)
  colData(sce_object)[[paste0(annotation_prefix, "_confidence_level")]] <- rep(NA, n_cells)
  colData(sce_object)[[paste0(annotation_prefix, "_score_diff")]] <- rep(NA, n_cells)

  if (top_n >= 2) {
    colData(sce_object)[[paste0(annotation_prefix, "_top2")]] <- rep(NA, n_cells)
    colData(sce_object)[[paste0(annotation_prefix, "_score2")]] <- rep(NA, n_cells)
  }
  if (top_n >= 3) {
    colData(sce_object)[[paste0(annotation_prefix, "_top3")]] <- rep(NA, n_cells)
    colData(sce_object)[[paste0(annotation_prefix, "_score3")]] <- rep(NA, n_cells)
  }

  # Fill in values
  for (i in 1:nrow(cluster_results)) {
    cl <- cluster_results$cluster[i]
    cells_in_cluster <- which(clusters == cl)

    colData(sce_object)[[paste0(annotation_prefix, "_top1")]][cells_in_cluster] <-
      cluster_results$top1[i]
    colData(sce_object)[[paste0(annotation_prefix, "_score1")]][cells_in_cluster] <-
      cluster_results$score1[i]
    colData(sce_object)[[paste0(annotation_prefix, "_confidence")]][cells_in_cluster] <-
      cluster_results$confidence[i]
    colData(sce_object)[[paste0(annotation_prefix, "_confidence_level")]][cells_in_cluster] <-
      cluster_results$confidence_level[i]
    colData(sce_object)[[paste0(annotation_prefix, "_score_diff")]][cells_in_cluster] <-
      cluster_results$score_diff[i]

    if (top_n >= 2) {
      colData(sce_object)[[paste0(annotation_prefix, "_top2")]][cells_in_cluster] <-
        cluster_results$top2[i]
      colData(sce_object)[[paste0(annotation_prefix, "_score2")]][cells_in_cluster] <-
        cluster_results$score2[i]
    }
    if (top_n >= 3) {
      colData(sce_object)[[paste0(annotation_prefix, "_top3")]][cells_in_cluster] <-
        cluster_results$top3[i]
      colData(sce_object)[[paste0(annotation_prefix, "_score3")]][cells_in_cluster] <-
        cluster_results$score3[i]
    }
  }

  # Store cluster-level results as metadata
  metadata(sce_object)[["sctype_uncertainty_clusters"]] <- cluster_results

  message("Done! Added the following columns to colData:")
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
  print(table(colData(sce_object)[[paste0(annotation_prefix, "_confidence_level")]]))

  return(sce_object)
}


#' Visualize ScType uncertainty scores (SCE)
#'
#' Creates comprehensive visualizations of uncertainty metrics including:
#' 1. Bar plots of top N candidates per cluster with scores
#' 2. UMAP colored by confidence score
#' 3. Distribution of confidence levels
#' 4. Uncertainty heatmap
#'
#' @param sce_object SingleCellExperiment object with uncertainty scores
#' @param annotation_prefix Prefix used in add_sctype_uncertainty_sce (default: "sctype")
#' @param cluster_col Cluster column (default: "cluster")
#' @param plot_types Vector of plot types: "candidates", "umap", "distribution", "heatmap" (default: all)
#' @param save_plots Save plots to files (default: FALSE)
#' @param output_dir Directory for saved plots (default: "uncertainty_plots")
#'
#' @return List of ggplot objects
#' @export
visualize_sctype_uncertainty_sce <- function(sce_object,
                                            annotation_prefix = "sctype",
                                            cluster_col = "cluster",
                                            plot_types = c("candidates", "umap", "distribution", "heatmap"),
                                            save_plots = FALSE,
                                            output_dir = "uncertainty_plots") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' required. Install with: install.packages('dplyr')")
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' required. Install with: BiocManager::install('SingleCellExperiment')")
  }

  # Check if uncertainty columns exist
  required_cols <- c(paste0(annotation_prefix, "_top1"),
                    paste0(annotation_prefix, "_confidence"))
  if (!all(required_cols %in% colnames(SingleCellExperiment::colData(sce_object)))) {
    stop("Uncertainty scores not found. Run add_sctype_uncertainty_sce() first.")
  }

  # Get cluster-level results
  if (is.null(S4Vectors::metadata(sce_object)[["sctype_uncertainty_clusters"]])) {
    stop("Cluster-level uncertainty data not found. Run add_sctype_uncertainty_sce() first.")
  }
  cluster_results <- S4Vectors::metadata(sce_object)[["sctype_uncertainty_clusters"]]

  # Create output directory
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  plot_list <- list()

  # 1. TOP CANDIDATES BAR PLOT
  if ("candidates" %in% plot_types) {
    message("Generating top candidates plot...")
    p <- plot_top_candidates_sce(cluster_results, annotation_prefix)
    plot_list$candidates <- p

    if (save_plots) {
      ggsave(file.path(output_dir, "top_candidates.png"), p,
             width = 14, height = 10, dpi = 300)
    }
  }

  # 2. UMAP CONFIDENCE
  if ("umap" %in% plot_types) {
    message("Generating UMAP confidence plots...")
    umap_plots <- plot_umap_confidence_sce(sce_object, annotation_prefix)
    plot_list$umap <- umap_plots

    if (save_plots && length(umap_plots) > 0) {
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
    dist_plots <- plot_confidence_distribution_sce(sce_object, cluster_results,
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
    p <- plot_uncertainty_heatmap_sce(cluster_results)
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


#' Plot top candidate cell types with scores (SCE)
#' @keywords internal
plot_top_candidates_sce <- function(cluster_results, annotation_prefix) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
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
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = factor(cluster), y = score, fill = rank)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = celltype), position = ggplot2::position_dodge(width = 0.9),
              angle = 45, hjust = -0.1, vjust = 0.5, size = 3) +
    ggplot2::scale_fill_manual(values = c("1st" = "#E41A1C", "2nd" = "#377EB8", "3rd" = "#4DAF4A")) +
    ggplot2::labs(title = "Top Cell Type Candidates per Cluster",
         subtitle = "Showing top 3 candidates with ScType scores",
         x = "Cluster", y = "ScType Score", fill = "Rank") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
          legend.position = "right",
          plot.title = ggplot2::element_text(face = "bold", size = 14))

  return(p)
}


#' Plot UMAP colored by confidence (SCE)
#' @keywords internal
plot_umap_confidence_sce <- function(sce_object, annotation_prefix) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' required. Install with: BiocManager::install('SingleCellExperiment')")
  }

  plots <- list()

  # Check if UMAP exists
  if (!"UMAP" %in% SingleCellExperiment::reducedDimNames(sce_object)) {
    warning("UMAP not found in reducedDims. Skipping UMAP plots.")
    return(plots)
  }

  if (!requireNamespace("scater", quietly = TRUE)) {
    warning("scater package needed for UMAP plots. Skipping.")
    return(plots)
  }

  # Plot 1: Confidence score
  p1 <- scater::plotReducedDim(sce_object, dimred = "UMAP",
                      colour_by = paste0(annotation_prefix, "_confidence"),
                      by_exprs_values = NULL) +
    ggplot2::scale_color_gradient2(low = "blue", mid = "yellow", high = "red",
                         midpoint = 0.5) +
    ggplot2::ggtitle("Annotation Confidence Score")

  # Plot 2: Confidence level
  p2 <- scater::plotReducedDim(sce_object, dimred = "UMAP",
                      colour_by = paste0(annotation_prefix, "_confidence_level"),
                      by_exprs_values = NULL) +
    ggplot2::ggtitle("Annotation Confidence Level")

  plots <- list(score = p1, level = p2)

  return(plots)
}


#' Plot confidence distribution (SCE)
#' @keywords internal
plot_confidence_distribution_sce <- function(sce_object, cluster_results,
                                            annotation_prefix, cluster_col) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' required. Install with: BiocManager::install('SingleCellExperiment')")
  }

  plots <- list()

  # Plot 1: Confidence level counts
  conf_data <- data.frame(
    level = SingleCellExperiment::colData(sce_object)[[paste0(annotation_prefix, "_confidence_level")]]
  )
  conf_data$level <- factor(conf_data$level, levels = c("High", "Medium", "Low", "Unknown"))

  p1 <- ggplot2::ggplot(conf_data, ggplot2::aes(x = level, fill = level)) +
    ggplot2::geom_bar() +
    ggplot2::scale_fill_manual(values = c("High" = "#2ECC40", "Medium" = "#FFDC00",
                                 "Low" = "#FF4136", "Unknown" = "#AAAAAA")) +
    ggplot2::labs(title = "Distribution of Confidence Levels",
         x = "Confidence Level", y = "Number of Cells") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  # Plot 2: Score difference distribution
  p2 <- ggplot2::ggplot(cluster_results, ggplot2::aes(x = factor(cluster), y = score_diff,
                                   fill = confidence_level)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = c("High" = "#2ECC40", "Medium" = "#FFDC00",
                                 "Low" = "#FF4136", "Unknown" = "#AAAAAA")) +
    ggplot2::labs(title = "Score Difference per Cluster (Top1 - Top2)",
         subtitle = "Larger difference = higher confidence",
         x = "Cluster", y = "Score Difference", fill = "Confidence") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  plots <- list(levels = p1, score_diff = p2)

  return(plots)
}


#' Plot uncertainty heatmap (SCE)
#' @keywords internal
plot_uncertainty_heatmap_sce <- function(cluster_results) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' required. Install with: install.packages('dplyr')")
  }

  # Prepare data
  heatmap_data <- cluster_results %>%
    dplyr::select(cluster, top1, score1, confidence) %>%
    dplyr::mutate(cluster = factor(cluster))

  # Create plot
  p <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = cluster, y = 1, fill = confidence)) +
    ggplot2::geom_tile(color = "white", size = 1) +
    ggplot2::geom_text(ggplot2::aes(label = top1), size = 3, fontface = "bold") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", score1)), size = 2.5,
              vjust = 2.5, color = "gray30") +
    ggplot2::scale_fill_gradient2(low = "red", mid = "yellow", high = "green",
                        midpoint = 0.5, limits = c(0, 1)) +
    ggplot2::labs(title = "Cluster Annotation Confidence",
         subtitle = "Cell type (top) and score (bottom) shown in each tile",
         x = "Cluster", y = "", fill = "Confidence") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank())

  return(p)
}
