# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/master/LICENSE)
# ScType marker gene visualization functions for SingleCellExperiment objects
# Written by Claude (AI assistant), 2025-11-15

#' Visualize marker genes used for cell type annotation (SCE)
#'
#' Creates comprehensive visualizations of the marker genes that were used to
#' determine each cell type annotation. Includes violin plots, UMAPs, dotplots,
#' and heatmaps for both positive and negative markers.
#'
#' @param sce_object SingleCellExperiment object with ScType annotations
#' @param annotation_col Column name in colData with cell type annotations (default: "sctype_classification")
#' @param database_file Path to marker database used for annotation
#' @param tissue_type Tissue type that was used for annotation
#' @param assay_name Assay to use (default: "logcounts")
#' @param top_n Number of top markers to show per cell type (default: 5)
#' @param plot_types Vector of plot types to generate: "violin", "umap", "dotplot", "heatmap" (default: all)
#' @param save_plots Save plots to files (default: FALSE)
#' @param output_dir Directory to save plots if save_plots = TRUE (default: "sctype_plots")
#'
#' @return List of ggplot objects (invisible)
#' @export
#'
#' @examples
#' \dontrun{
#' plots <- visualize_sctype_markers_sce(sce,
#'                                        annotation_col = "sctype_classification",
#'                                        database_file = "ScTypeDB_full.xlsx",
#'                                        tissue_type = "Immune system")
#' }
visualize_sctype_markers_sce <- function(sce_object,
                                         annotation_col = "sctype_classification",
                                         database_file = NULL,
                                         tissue_type = NULL,
                                         assay_name = "logcounts",
                                         top_n = 5,
                                         plot_types = c("violin", "umap", "dotplot", "heatmap"),
                                         save_plots = FALSE,
                                         output_dir = "sctype_plots") {

  # Check required packages
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' required. Install with: BiocManager::install('SingleCellExperiment')")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' required. Install with: install.packages('dplyr')")
  }
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' required. Install with: install.packages('openxlsx')")
  }

  # Check if annotation column exists
  if (!annotation_col %in% colnames(SingleCellExperiment::colData(sce_object))) {
    stop(paste0("Annotation column '", annotation_col, "' not found in colData."))
  }

  # Set default database file
  if (is.null(database_file)) {
    database_file <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
  }

  # Auto-detect tissue type if not provided
  if (is.null(tissue_type)) {
    stop("Please specify tissue_type parameter (e.g., 'Immune system', 'Brain', etc.)")
  }

  # Read marker database
  message("Reading marker database...")
  cell_markers <- openxlsx::read.xlsx(database_file)
  cell_markers <- cell_markers[cell_markers$tissueType == tissue_type, ]

  if (nrow(cell_markers) == 0) {
    stop(paste0("No markers found for tissue type: ", tissue_type))
  }

  # Get unique cell types from annotations (excluding "Unknown")
  annotated_types <- unique(SingleCellExperiment::colData(sce_object)[[annotation_col]])
  annotated_types <- annotated_types[annotated_types != "Unknown"]

  if (length(annotated_types) == 0) {
    stop("No cell types found (all are 'Unknown'). Please run ScType annotation first.")
  }

  # Check if assay exists
  if (!assay_name %in% SummarizedExperiment::assayNames(sce_object)) {
    stop(paste0("Assay '", assay_name, "' not found. Available assays: ",
                paste(SummarizedExperiment::assayNames(sce_object), collapse = ", ")))
  }

  # Extract markers for each annotated cell type
  message("Extracting markers for annotated cell types...")
  markers_list <- list()

  for (cell_type in annotated_types) {
    marker_row <- cell_markers[cell_markers$cellName == cell_type, ]

    if (nrow(marker_row) > 0) {
      # Parse positive markers
      pos_markers <- unlist(strsplit(marker_row$geneSymbolmore1[1], ","))
      pos_markers <- trimws(pos_markers)
      pos_markers <- pos_markers[pos_markers != "" & !is.na(pos_markers)]

      # Parse negative markers
      neg_markers <- character(0)
      if ("geneSymbolmore2" %in% colnames(marker_row)) {
        neg_markers <- unlist(strsplit(marker_row$geneSymbolmore2[1], ","))
        neg_markers <- trimws(neg_markers)
        neg_markers <- neg_markers[neg_markers != "" & !is.na(neg_markers)]
      }

      # Filter markers that exist in the dataset
      all_genes <- rownames(sce_object)
      pos_markers <- pos_markers[pos_markers %in% all_genes]
      neg_markers <- neg_markers[neg_markers %in% all_genes]

      # Limit to top_n
      if (length(pos_markers) > top_n) {
        pos_markers <- pos_markers[1:top_n]
      }
      if (length(neg_markers) > top_n) {
        neg_markers <- neg_markers[1:top_n]
      }

      markers_list[[cell_type]] <- list(
        positive = pos_markers,
        negative = neg_markers
      )
    }
  }

  # Create output directory if saving plots
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Generate plots
  plot_list <- list()

  # 1. VIOLIN PLOTS
  if ("violin" %in% plot_types) {
    message("Generating violin plots...")
    violin_plots <- generate_violin_plots_sce(sce_object, markers_list, annotation_col, assay_name)
    plot_list$violin <- violin_plots

    if (save_plots) {
      for (i in seq_along(violin_plots)) {
        cell_type <- names(violin_plots)[i]
        filename <- file.path(output_dir, paste0("violin_", gsub(" ", "_", cell_type), ".png"))
        ggsave(filename, violin_plots[[i]], width = 12, height = 8, dpi = 300)
      }
    }
  }

  # 2. UMAP PLOTS
  if ("umap" %in% plot_types) {
    message("Generating UMAP plots...")
    umap_plots <- generate_umap_plots_sce(sce_object, markers_list, annotation_col, assay_name)
    plot_list$umap <- umap_plots

    if (save_plots) {
      for (i in seq_along(umap_plots)) {
        cell_type <- names(umap_plots)[i]
        filename <- file.path(output_dir, paste0("umap_", gsub(" ", "_", cell_type), ".png"))
        ggsave(filename, umap_plots[[i]], width = 14, height = 10, dpi = 300)
      }
    }
  }

  # 3. DOTPLOT
  if ("dotplot" %in% plot_types) {
    message("Generating dotplot...")
    dotplot <- generate_dotplot_sce(sce_object, markers_list, annotation_col, assay_name)
    plot_list$dotplot <- dotplot

    if (save_plots) {
      filename <- file.path(output_dir, "dotplot_all_markers.png")
      ggsave(filename, dotplot, width = 16, height = 10, dpi = 300)
    }
  }

  # 4. HEATMAP
  if ("heatmap" %in% plot_types) {
    message("Generating heatmap...")
    heatmap_plot <- generate_heatmap_sce(sce_object, markers_list, annotation_col, assay_name)
    plot_list$heatmap <- heatmap_plot

    if (save_plots && !is.null(heatmap_plot)) {
      filename <- file.path(output_dir, "heatmap_all_markers.png")
      png(filename, width = 14, height = 10, units = "in", res = 300)
      print(heatmap_plot)
      dev.off()
    }
  }

  message("Visualization complete!")
  if (save_plots) {
    message(paste0("Plots saved to: ", output_dir))
  }

  invisible(plot_list)
}


#' Generate violin plots for marker genes (SCE)
#' @keywords internal
generate_violin_plots_sce <- function(sce_object, markers_list, annotation_col, assay_name) {
  if (!requireNamespace("scater", quietly = TRUE)) {
    warning("scater package not installed. Skipping violin plots.")
    return(list())
  }

  violin_plots <- list()

  for (cell_type in names(markers_list)) {
    markers <- markers_list[[cell_type]]
    all_markers <- c(markers$positive, markers$negative)

    if (length(all_markers) == 0) next

    # Create marker type labels
    marker_types <- c(
      rep("Positive", length(markers$positive)),
      rep("Negative", length(markers$negative))
    )

    plots <- list()
    for (i in seq_along(all_markers)) {
      gene <- all_markers[i]
      marker_type <- marker_types[i]

      if (gene %in% rownames(sce_object)) {
        p <- scater::plotExpression(sce_object, features = gene,
                           x = annotation_col,
                           exprs_values = assay_name,
                           colour_by = annotation_col) +
          ggplot2::ggtitle(paste0(gene, " (", marker_type, ")")) +
          ggplot2::theme(legend.position = "none",
                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
          ggplot2::ylab("Expression")
        plots[[gene]] <- p
      }
    }

    if (length(plots) > 0) {
      # Combine plots using patchwork
      if (requireNamespace("patchwork", quietly = TRUE)) {
        combined <- patchwork::wrap_plots(plots, ncol = 3) +
          patchwork::plot_annotation(title = paste0("Marker Expression: ", cell_type))
        violin_plots[[cell_type]] <- combined
      } else {
        violin_plots[[cell_type]] <- plots
      }
    }
  }

  return(violin_plots)
}


#' Generate UMAP plots for marker genes (SCE)
#' @keywords internal
generate_umap_plots_sce <- function(sce_object, markers_list, annotation_col, assay_name) {
  if (!requireNamespace("scater", quietly = TRUE)) {
    warning("scater package not installed. Skipping UMAP plots.")
    return(list())
  }

  umap_plots <- list()

  # Check if UMAP exists
  if (!"UMAP" %in% SingleCellExperiment::reducedDimNames(sce_object)) {
    warning("UMAP not found in reducedDims. Skipping UMAP plots. Run runUMAP() first.")
    return(list())
  }

  for (cell_type in names(markers_list)) {
    markers <- markers_list[[cell_type]]
    all_markers <- c(markers$positive, markers$negative)

    if (length(all_markers) == 0) next

    # Create marker type labels
    marker_types <- c(
      rep("Positive", length(markers$positive)),
      rep("Negative", length(markers$negative))
    )

    plots <- list()
    for (i in seq_along(all_markers)) {
      gene <- all_markers[i]
      marker_type <- marker_types[i]

      if (gene %in% rownames(sce_object)) {
        p <- scater::plotReducedDim(sce_object, dimred = "UMAP",
                           colour_by = gene,
                           by_exprs_values = assay_name) +
          ggplot2::ggtitle(paste0(gene, " (", marker_type, ")"))
        plots[[gene]] <- p
      }
    }

    if (length(plots) > 0) {
      # Combine plots using patchwork
      if (requireNamespace("patchwork", quietly = TRUE)) {
        combined <- patchwork::wrap_plots(plots, ncol = 3) +
          patchwork::plot_annotation(title = paste0("Marker Expression (UMAP): ", cell_type))
        umap_plots[[cell_type]] <- combined
      } else {
        umap_plots[[cell_type]] <- plots
      }
    }
  }

  return(umap_plots)
}


#' Generate dotplot for all marker genes (SCE)
#' @keywords internal
generate_dotplot_sce <- function(sce_object, markers_list, annotation_col, assay_name) {
  if (!requireNamespace("scater", quietly = TRUE)) {
    warning("scater package not installed. Skipping dotplot.")
    return(NULL)
  }

  # Collect all unique markers
  all_markers <- unique(unlist(lapply(markers_list, function(x) c(x$positive, x$negative))))
  all_markers <- all_markers[all_markers %in% rownames(sce_object)]

  if (length(all_markers) == 0) {
    warning("No valid markers found for dotplot.")
    return(NULL)
  }

  # Create dotplot
  p <- scater::plotDots(sce_object, features = all_markers,
                group = annotation_col,
                exprs_values = assay_name) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::ggtitle("Marker Gene Expression Across Cell Types") +
    ggplot2::xlab("Genes") + ggplot2::ylab("Cell Types")

  return(p)
}


#' Generate heatmap for all marker genes (SCE)
#' @keywords internal
generate_heatmap_sce <- function(sce_object, markers_list, annotation_col, assay_name) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    warning("ComplexHeatmap package not installed. Skipping heatmap. Install with: BiocManager::install('ComplexHeatmap')")
    return(NULL)
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    warning("circlize package not installed. Install with: install.packages('circlize')")
    return(NULL)
  }

  # Collect all markers
  all_markers <- unique(unlist(lapply(markers_list, function(x) c(x$positive, x$negative))))
  all_markers <- all_markers[all_markers %in% rownames(sce_object)]

  if (length(all_markers) == 0) {
    warning("No valid markers found for heatmap.")
    return(NULL)
  }

  # Get expression data
  expr_data <- as.matrix(SummarizedExperiment::assay(sce_object, assay_name)[all_markers, , drop = FALSE])

  # Average expression by cell type
  cell_types <- SingleCellExperiment::colData(sce_object)[[annotation_col]]
  unique_types <- unique(cell_types[cell_types != "Unknown"])

  avg_expr <- matrix(0, nrow = length(all_markers), ncol = length(unique_types))
  rownames(avg_expr) <- all_markers
  colnames(avg_expr) <- unique_types

  for (ct in unique_types) {
    cells <- which(cell_types == ct)
    if (length(cells) > 0) {
      avg_expr[, ct] <- rowMeans(expr_data[, cells, drop = FALSE])
    }
  }

  # Scale for visualization
  avg_expr_scaled <- t(scale(t(avg_expr)))

  # Create marker annotations
  marker_anno <- data.frame(
    Gene = all_markers,
    Type = "Unknown",
    stringsAsFactors = FALSE
  )

  for (cell_type in names(markers_list)) {
    pos_markers <- markers_list[[cell_type]]$positive
    neg_markers <- markers_list[[cell_type]]$negative

    marker_anno$Type[marker_anno$Gene %in% pos_markers] <-
      paste0(cell_type, " (Pos)")
    marker_anno$Type[marker_anno$Gene %in% neg_markers] <-
      paste0(cell_type, " (Neg)")
  }

  # Create heatmap
  col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

  ht <- ComplexHeatmap::Heatmap(avg_expr_scaled,
                name = "Scaled\nExpression",
                col = col_fun,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                row_names_gp = grid::gpar(fontsize = 8),
                column_names_gp = grid::gpar(fontsize = 10),
                heatmap_legend_param = list(title = "Scaled\nExpression"))

  return(ht)
}


#' Quick visualization of markers for a specific cell type (SCE)
#'
#' Convenience function to quickly visualize markers for one cell type
#'
#' @param sce_object SingleCellExperiment object with ScType annotations
#' @param cell_type Cell type to visualize
#' @param annotation_col Column name with cell type annotations
#' @param database_file Path to marker database
#' @param tissue_type Tissue type
#' @param plot_type Type of plot: "violin", "umap", or "both" (default: "both")
#'
#' @return ggplot object or list of plots
#' @export
quick_marker_viz_sce <- function(sce_object,
                                 cell_type,
                                 annotation_col = "sctype_classification",
                                 database_file = NULL,
                                 tissue_type = NULL,
                                 plot_type = "both") {

  plots <- visualize_sctype_markers_sce(
    sce_object,
    annotation_col = annotation_col,
    database_file = database_file,
    tissue_type = tissue_type,
    plot_types = c("violin", "umap"),
    save_plots = FALSE
  )

  result <- list()
  if (plot_type %in% c("violin", "both") && !is.null(plots$violin[[cell_type]])) {
    result$violin <- plots$violin[[cell_type]]
  }
  if (plot_type %in% c("umap", "both") && !is.null(plots$umap[[cell_type]])) {
    result$umap <- plots$umap[[cell_type]]
  }

  if (length(result) == 1) {
    return(result[[1]])
  } else {
    return(result)
  }
}
