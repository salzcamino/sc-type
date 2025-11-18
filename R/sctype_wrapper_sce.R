# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/master/LICENSE)
# ScType wrapper for SingleCellExperiment objects
# Written by Claude (AI assistant), 2025-11-15

#' Run ScType annotation on SingleCellExperiment object
#'
#' This is a simplified wrapper function that performs cell type annotation
#' on a SingleCellExperiment object using the ScType algorithm.
#'
#' @param sce_object SingleCellExperiment object with clustering already performed
#' @param known_tissue_type Tissue type (e.g., "Immune system", "Brain"). If NULL, auto-detected.
#' @param assay_name Assay name to use (default: "logcounts" for scaled data, "counts" for raw)
#' @param scaled Whether the assay data is already scaled (default: TRUE)
#' @param cluster_col Column name in colData containing cluster assignments (default: "cluster")
#' @param custom_marker_file Path to marker database (default: uses ScTypeDB_full.xlsx from GitHub)
#' @param plot Generate UMAP plot (default: FALSE)
#' @param name Column name for cell type annotations (default: "sctype_classification")
#'
#' @return Modified SingleCellExperiment object with new colData column
#' @export
#'
#' @examples
#' \dontrun{
#' library(SingleCellExperiment)
#' sce <- run_sctype_sce(sce, known_tissue_type = "Immune system")
#' table(colData(sce)$sctype_classification)
#' }
run_sctype_sce <- function(sce_object,
                          known_tissue_type = NULL,
                          assay_name = "logcounts",
                          scaled = TRUE,
                          cluster_col = "cluster",
                          custom_marker_file = NULL,
                          plot = FALSE,
                          name = "sctype_classification") {

  # Load required packages
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required but not installed.")
  }

  lapply(c("dplyr", "openxlsx", "HGNChelper"), library, character.only = TRUE)

  # ScType functions are available from package namespace
  # Functions available: gene_sets_prepare, sctype_score

  # Set default marker file
  if (is.null(custom_marker_file)) {
    custom_marker_file <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
  }

  # Check if cluster column exists
  if (!cluster_col %in% colnames(SingleCellExperiment::colData(sce_object))) {
    stop(paste0("Cluster column '", cluster_col, "' not found in colData. ",
                "Please perform clustering first or specify correct cluster_col parameter."))
  }

  # Extract expression matrix
  if (!assay_name %in% SummarizedExperiment::assayNames(sce_object)) {
    stop(paste0("Assay '", assay_name, "' not found. Available assays: ",
                paste(SummarizedExperiment::assayNames(sce_object), collapse = ", ")))
  }

  scRNAseqData <- as.matrix(SummarizedExperiment::assay(sce_object, assay_name))

  # Auto-detect tissue type if not provided
  if (is.null(known_tissue_type)) {
    message("Tissue type not specified. Auto-detecting...")

    # Read database to get available tissue types
    cell_markers <- openxlsx::read.xlsx(custom_marker_file)
    available_tissues <- unique(cell_markers$tissueType)

    # Score each tissue type
    tissue_scores <- sapply(available_tissues, function(tissue) {
      gs_list <- gene_sets_prepare(custom_marker_file, tissue)
      es.max <- sctype_score(scRNAseqData = scRNAseqData,
                            scaled = scaled,
                            gs = gs_list$gs_positive,
                            gs2 = gs_list$gs_negative)
      mean(es.max, na.rm = TRUE)
    })

    known_tissue_type <- names(tissue_scores)[which.max(tissue_scores)]
    message(paste0("Auto-detected tissue type: ", known_tissue_type))
  }

  # Prepare gene sets
  gs_list <- gene_sets_prepare(custom_marker_file, known_tissue_type)

  # Run ScType scoring
  es.max <- sctype_score(scRNAseqData = scRNAseqData,
                        scaled = scaled,
                        gs = gs_list$gs_positive,
                        gs2 = gs_list$gs_negative)

  # Get cluster assignments
  clusters <- SingleCellExperiment::colData(sce_object)[[cluster_col]]

  # Aggregate scores by cluster
  cL_results <- do.call("rbind", lapply(unique(clusters), function(cl) {
    cells_in_cluster <- which(clusters == cl)
    es.max.cl <- sort(rowSums(es.max[, cells_in_cluster, drop = FALSE]), decreasing = TRUE)
    head(data.frame(
      cluster = cl,
      type = names(es.max.cl),
      scores = es.max.cl,
      ncells = length(cells_in_cluster),
      stringsAsFactors = FALSE
    ), 10)
  }))

  # Get top cell type per cluster
  sctype_scores <- cL_results %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = scores) %>%
    ungroup()

  # Filter low-confidence assignments
  sctype_scores$type[sctype_scores$scores < sctype_scores$ncells/4] <- "Unknown"

  # Assign cell types to individual cells
  cell_type_annotation <- rep("Unknown", ncol(sce_object))
  for (i in 1:nrow(sctype_scores)) {
    cl <- sctype_scores$cluster[i]
    cells_in_cluster <- which(clusters == cl)
    cell_type_annotation[cells_in_cluster] <- sctype_scores$type[i]
  }

  # Add to colData
  SingleCellExperiment::colData(sce_object)[[name]] <- cell_type_annotation

  message(paste0("Added '", name, "' to colData with ", length(unique(cell_type_annotation)),
                " unique cell types."))

  # Print summary
  cat("\nCell type distribution:\n")
  print(table(cell_type_annotation))

  # Plotting (optional)
  if (plot) {
    if (!"UMAP" %in% SingleCellExperiment::reducedDimNames(sce_object)) {
      warning("UMAP not found in reducedDims. Skipping plot. Run scater::runUMAP() first.")
    } else {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("ggplot2 package required for plotting. Skipping plot.")
      } else {
        library(ggplot2)

        # Get UMAP coordinates
        umap_coords <- SingleCellExperiment::reducedDim(sce_object, "UMAP")
        plot_df <- data.frame(
          UMAP1 = umap_coords[, 1],
          UMAP2 = umap_coords[, 2],
          CellType = SingleCellExperiment::colData(sce_object)[[name]]
        )

        p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = CellType)) +
          geom_point(size = 0.5) +
          theme_minimal() +
          ggtitle("ScType Cell Type Annotations") +
          theme(legend.position = "right")

        print(p)
      }
    }
  }

  return(sce_object)
}


#' Source ScType functions and return database URL
#'
#' Convenience function to load all ScType functions
#'
#' @return Character vector with database URLs
#' @export
sctype_source_sce <- function() {
  # ScType functions are now loaded from the package namespace
  # No need to source remote files when using as an installed package
  # All core functions (gene_sets_prepare, sctype_score) are available
  # automatically when the package is loaded

  db_urls <- c(
    full = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",
    short = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx",
    enhanced = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_enhanced.xlsx",
    hierarchical = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_hierarchical.xlsx"
  )

  message("ScType functions loaded successfully!")
  message("\nAvailable databases:")
  for (db_name in names(db_urls)) {
    message(paste0("  ", db_name, ": ", db_urls[db_name]))
  }

  return(db_urls)
}
