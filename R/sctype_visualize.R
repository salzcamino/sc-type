# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/master/LICENSE)
# ScType marker gene visualization functions for Seurat objects
# Written by Claude (AI assistant), 2025-11-15

#' Visualize marker genes used for cell type annotation (Seurat)
#'
#' Creates comprehensive visualizations of the marker genes that were used to
#' determine each cell type annotation. Includes violin plots, UMAPs, dotplots,
#' and heatmaps for both positive and negative markers.
#'
#' @param seurat_object Seurat object with ScType annotations
#' @param annotation_col Column name with cell type annotations (default: "sctype_classification")
#' @param database_file Path to marker database used for annotation
#' @param tissue_type Tissue type that was used for annotation
#' @param assay Assay to use (default: "RNA")
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
#' plots <- visualize_sctype_markers(seurat_obj,
#'                                    annotation_col = "sctype_classification",
#'                                    database_file = "ScTypeDB_full.xlsx",
#'                                    tissue_type = "Immune system")
#' }
visualize_sctype_markers <- function(seurat_object,
                                     annotation_col = "sctype_classification",
                                     database_file = NULL,
                                     tissue_type = NULL,
                                     assay = "RNA",
                                     top_n = 5,
                                     plot_types = c("violin", "umap", "dotplot", "heatmap"),
                                     save_plots = FALSE,
                                     output_dir = "sctype_plots") {

  # Check required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' required. Install with: install.packages('dplyr')")
  }
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' required. Install with: install.packages('openxlsx')")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' required. Install with: install.packages('Seurat')")
  }

  # Check if annotation column exists
  if (!annotation_col %in% colnames(seurat_object@meta.data)) {
    stop(paste0("Annotation column '", annotation_col, "' not found in metadata."))
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
  annotated_types <- unique(seurat_object@meta.data[[annotation_col]])
  annotated_types <- annotated_types[annotated_types != "Unknown"]

  if (length(annotated_types) == 0) {
    stop("No cell types found (all are 'Unknown'). Please run ScType annotation first.")
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
      all_genes <- rownames(seurat_object)
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
    violin_plots <- generate_violin_plots(seurat_object, markers_list, annotation_col, assay)
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
    umap_plots <- generate_umap_plots(seurat_object, markers_list, annotation_col, assay)
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
    dotplot <- generate_dotplot(seurat_object, markers_list, annotation_col, assay)
    plot_list$dotplot <- dotplot

    if (save_plots) {
      filename <- file.path(output_dir, "dotplot_all_markers.png")
      ggsave(filename, dotplot, width = 16, height = 10, dpi = 300)
    }
  }

  # 4. HEATMAP
  if ("heatmap" %in% plot_types) {
    message("Generating heatmap...")
    heatmap_plot <- generate_heatmap(seurat_object, markers_list, annotation_col, assay)
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


#' Generate violin plots for marker genes
#' @keywords internal
generate_violin_plots <- function(seurat_object, markers_list, annotation_col, assay) {
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

      if (gene %in% rownames(seurat_object)) {
        p <- VlnPlot(seurat_object, features = gene, group.by = annotation_col,
                     pt.size = 0, assay = assay) +
          ggtitle(paste0(gene, " (", marker_type, ")")) +
          theme(legend.position = "none",
                axis.text.x = element_text(angle = 45, hjust = 1))
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


#' Generate UMAP plots for marker genes
#' @keywords internal
generate_umap_plots <- function(seurat_object, markers_list, annotation_col, assay) {
  umap_plots <- list()

  # Check if UMAP exists
  if (!"umap" %in% names(seurat_object@reductions)) {
    warning("UMAP not found. Skipping UMAP plots. Run RunUMAP() first.")
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

      if (gene %in% rownames(seurat_object)) {
        p <- FeaturePlot(seurat_object, features = gene, reduction = "umap") +
          ggtitle(paste0(gene, " (", marker_type, ")"))
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


#' Generate dotplot for all marker genes
#' @keywords internal
generate_dotplot <- function(seurat_object, markers_list, annotation_col, assay) {
  # Collect all markers
  all_pos_markers <- list()
  all_neg_markers <- list()

  for (cell_type in names(markers_list)) {
    markers <- markers_list[[cell_type]]
    if (length(markers$positive) > 0) {
      all_pos_markers[[cell_type]] <- markers$positive
    }
    if (length(markers$negative) > 0) {
      all_neg_markers[[cell_type]] <- markers$negative
    }
  }

  # Create combined marker list
  unique_pos <- unique(unlist(all_pos_markers))
  unique_neg <- unique(unlist(all_neg_markers))
  all_unique <- unique(c(unique_pos, unique_neg))

  # Filter to genes present in dataset
  all_unique <- all_unique[all_unique %in% rownames(seurat_object)]

  if (length(all_unique) == 0) {
    warning("No valid markers found for dotplot.")
    return(NULL)
  }

  # Create dotplot
  p <- DotPlot(seurat_object, features = all_unique, group.by = annotation_col, assay = assay) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggtitle("Marker Gene Expression Across Cell Types") +
    xlab("Genes") + ylab("Cell Types")

  return(p)
}


#' Generate heatmap for all marker genes
#' @keywords internal
generate_heatmap <- function(seurat_object, markers_list, annotation_col, assay) {
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
  all_markers <- all_markers[all_markers %in% rownames(seurat_object)]

  if (length(all_markers) == 0) {
    warning("No valid markers found for heatmap.")
    return(NULL)
  }

  # Get expression data
  Seurat::DefaultAssay(seurat_object) <- assay
  expr_data <- Seurat::GetAssayData(seurat_object, slot = "data")
  expr_data <- as.matrix(expr_data[all_markers, , drop = FALSE])

  # Average expression by cell type
  cell_types <- seurat_object@meta.data[[annotation_col]]
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
  col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

  ht <- Heatmap(avg_expr_scaled,
                name = "Scaled\nExpression",
                col = col_fun,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 10),
                heatmap_legend_param = list(title = "Scaled\nExpression"))

  return(ht)
}


#' Quick visualization of markers for a specific cell type
#'
#' Convenience function to quickly visualize markers for one cell type
#'
#' @param seurat_object Seurat object with ScType annotations
#' @param cell_type Cell type to visualize
#' @param annotation_col Column name with cell type annotations
#' @param database_file Path to marker database
#' @param tissue_type Tissue type
#' @param plot_type Type of plot: "violin", "umap", or "both" (default: "both")
#'
#' @return ggplot object or list of plots
#' @export
quick_marker_viz <- function(seurat_object,
                             cell_type,
                             annotation_col = "sctype_classification",
                             database_file = NULL,
                             tissue_type = NULL,
                             plot_type = "both") {

  plots <- visualize_sctype_markers(
    seurat_object,
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
