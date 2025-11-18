# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/master/LICENSE)
# Hierarchical cell type annotation for SingleCellExperiment objects
# Written by Claude (AI assistant), 2025-11-15

#' Run hierarchical ScType annotation on SingleCellExperiment object
#'
#' This function performs two-level hierarchical cell type annotation on a
#' SingleCellExperiment object. It first annotates broad cell categories, then
#' refines to fine-grained subtypes where confidence is sufficient.
#'
#' @param sce_object SingleCellExperiment object with clustering already performed
#' @param known_tissue_type Tissue type (e.g., "Immune system", "Brain"). If NULL, auto-detected.
#' @param assay_name Assay name to use (default: "logcounts" for scaled data, "counts" for raw)
#' @param scaled Whether the assay data is already scaled (default: TRUE)
#' @param cluster_col Column name in colData containing cluster assignments (default: "cluster")
#' @param custom_marker_file Path to hierarchical marker database (default: uses GitHub URL)
#' @param plot Generate UMAP plots (default: FALSE). Requires UMAP in reducedDims.
#' @param broad_name Column name for broad category annotations (default: "sctype_broad")
#' @param fine_name Column name for fine subtype annotations (default: "sctype_fine")
#'
#' @return Modified SingleCellExperiment object with two new colData columns
#' @export
#'
#' @examples
#' \dontrun{
#' library(SingleCellExperiment)
#' sce <- run_sctype_hierarchical_sce(sce, known_tissue_type = "Immune system")
#' table(colData(sce)$sctype_broad)
#' table(colData(sce)$sctype_fine)
#' }
run_sctype_hierarchical_sce <- function(sce_object,
                                        known_tissue_type = NULL,
                                        assay_name = "logcounts",
                                        scaled = TRUE,
                                        cluster_col = "cluster",
                                        custom_marker_file = NULL,
                                        plot = FALSE,
                                        broad_name = "sctype_broad",
                                        fine_name = "sctype_fine") {

  # Load required packages
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required but not installed.")
  }

  # Source required ScType functions
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

  # Set default marker file
  if (is.null(custom_marker_file)) {
    custom_marker_file <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_hierarchical.xlsx"
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

    # Load required packages for auto-detection
    lapply(c("dplyr", "openxlsx"), library, character.only = TRUE)

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

  # Load marker database
  lapply(c("dplyr", "openxlsx"), library, character.only = TRUE)
  cell_markers <- openxlsx::read.xlsx(custom_marker_file)
  cell_markers <- cell_markers[cell_markers$tissueType == known_tissue_type, ]

  if (nrow(cell_markers) == 0) {
    stop(paste0("No markers found for tissue type: ", known_tissue_type))
  }

  # Get cluster assignments
  clusters <- SingleCellExperiment::colData(sce_object)[[cluster_col]]

  # ============================================================
  # STEP 1: BROAD CATEGORY ANNOTATION
  # ============================================================

  message("Step 1/2: Annotating broad cell categories...")

  # Create broad category marker sets by aggregating markers
  broad_categories <- unique(cell_markers$broadCategory)
  broad_markers_pos <- list()
  broad_markers_neg <- list()

  for (broad_cat in broad_categories) {
    subset_markers <- cell_markers[cell_markers$broadCategory == broad_cat, ]

    # Aggregate positive markers
    all_pos <- unlist(strsplit(paste(subset_markers$geneSymbolmore1, collapse = ","), ","))
    all_pos <- unique(trimws(all_pos))
    all_pos <- all_pos[all_pos != "" & !is.na(all_pos)]
    broad_markers_pos[[broad_cat]] <- all_pos

    # Aggregate negative markers
    all_neg <- unlist(strsplit(paste(subset_markers$geneSymbolmore2, collapse = ","), ","))
    all_neg <- unique(trimws(all_neg))
    all_neg <- all_neg[all_neg != "" & !is.na(all_neg)]
    broad_markers_neg[[broad_cat]] <- all_neg
  }

  # Run ScType scoring for broad categories
  es_broad <- sctype_score(scRNAseqData = scRNAseqData,
                          scaled = scaled,
                          gs = broad_markers_pos,
                          gs2 = if(length(unlist(broad_markers_neg)) > 0) broad_markers_neg else NULL)

  # Aggregate scores by cluster for broad categories
  broad_results <- do.call("rbind", lapply(unique(clusters), function(cl) {
    cells_in_cluster <- which(clusters == cl)
    es_broad_cl <- sort(rowSums(es_broad[, cells_in_cluster, drop = FALSE]), decreasing = TRUE)
    data.frame(
      cluster = cl,
      type = names(es_broad_cl),
      scores = es_broad_cl,
      ncells = length(cells_in_cluster),
      stringsAsFactors = FALSE
    )
  }))

  # Get top broad category per cluster
  broad_assignments <- broad_results %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = scores) %>%
    ungroup()

  # Filter low-confidence broad assignments
  broad_assignments$type[broad_assignments$scores < broad_assignments$ncells/4] <- "Unknown"

  # ============================================================
  # STEP 2: FINE-GRAINED ANNOTATION
  # ============================================================

  message("Step 2/2: Refining to fine-grained cell subtypes...")

  # Prepare fine-grained marker sets
  gs_list <- gene_sets_prepare(custom_marker_file, known_tissue_type)

  # Run ScType scoring for fine cell types
  es_fine <- sctype_score(scRNAseqData = scRNAseqData,
                         scaled = scaled,
                         gs = gs_list$gs_positive,
                         gs2 = gs_list$gs_negative)

  # Aggregate scores by cluster for fine types
  fine_results <- do.call("rbind", lapply(unique(clusters), function(cl) {
    cells_in_cluster <- which(clusters == cl)
    es_fine_cl <- sort(rowSums(es_fine[, cells_in_cluster, drop = FALSE]), decreasing = TRUE)
    head(data.frame(
      cluster = cl,
      type = names(es_fine_cl),
      scores = es_fine_cl,
      ncells = length(cells_in_cluster),
      stringsAsFactors = FALSE
    ), 10)
  }))

  # Get top fine type per cluster
  fine_assignments <- fine_results %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = scores) %>%
    ungroup()

  # Confidence-based assignment: use broad if fine confidence is low
  fine_assignments$type_confident <- fine_assignments$type
  fine_assignments$type_confident[fine_assignments$scores < fine_assignments$ncells/4] <-
    broad_assignments$type[match(
      fine_assignments$cluster[fine_assignments$scores < fine_assignments$ncells/4],
      broad_assignments$cluster
    )]

  # ============================================================
  # ASSIGN TO CELLS
  # ============================================================

  # Create assignment vectors
  broad_annotation <- rep("Unknown", ncol(sce_object))
  fine_annotation <- rep("Unknown", ncol(sce_object))

  for (i in 1:nrow(broad_assignments)) {
    cl <- broad_assignments$cluster[i]
    cells_in_cluster <- which(clusters == cl)
    broad_annotation[cells_in_cluster] <- broad_assignments$type[i]
  }

  for (i in 1:nrow(fine_assignments)) {
    cl <- fine_assignments$cluster[i]
    cells_in_cluster <- which(clusters == cl)
    fine_annotation[cells_in_cluster] <- fine_assignments$type_confident[i]
  }

  # Add to colData
  SingleCellExperiment::colData(sce_object)[[broad_name]] <- broad_annotation
  SingleCellExperiment::colData(sce_object)[[fine_name]] <- fine_annotation

  message(paste0("Added '", broad_name, "' and '", fine_name, "' to colData."))

  # ============================================================
  # PLOTTING (OPTIONAL)
  # ============================================================

  if (plot) {
    if (!"UMAP" %in% SingleCellExperiment::reducedDimNames(sce_object)) {
      warning("UMAP not found in reducedDims. Skipping plot. Run scater::runUMAP() first.")
    } else {
      if (!requireNamespace("ggplot2", quietly = TRUE) ||
          !requireNamespace("patchwork", quietly = TRUE)) {
        warning("ggplot2 and patchwork packages required for plotting. Skipping plot.")
      } else {
        # Get UMAP coordinates
        umap_coords <- SingleCellExperiment::reducedDim(sce_object, "UMAP")
        plot_df <- data.frame(
          UMAP1 = umap_coords[, 1],
          UMAP2 = umap_coords[, 2],
          broad = SingleCellExperiment::colData(sce_object)[[broad_name]],
          fine = SingleCellExperiment::colData(sce_object)[[fine_name]]
        )

        p1 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = UMAP1, y = UMAP2, color = broad)) +
          ggplot2::geom_point(size = 0.5) +
          ggplot2::theme_minimal() +
          ggplot2::ggtitle("Broad Cell Categories") +
          ggplot2::theme(legend.position = "right")

        p2 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = UMAP1, y = UMAP2, color = fine)) +
          ggplot2::geom_point(size = 0.5) +
          ggplot2::theme_minimal() +
          ggplot2::ggtitle("Fine Cell Subtypes") +
          ggplot2::theme(legend.position = "right")

        print(p1 / p2)
      }
    }
  }

  return(sce_object)
}


#' Get hierarchical annotations for a specific cluster (SCE version)
#'
#' @param sce_object SingleCellExperiment object with hierarchical annotations
#' @param cluster_id Cluster ID to query
#' @param cluster_col Column name in colData containing cluster assignments
#' @param broad_name Column name for broad annotations
#' @param fine_name Column name for fine annotations
#'
#' @return List with cluster info and both annotation levels
#' @export
get_cluster_hierarchy_sce <- function(sce_object,
                                      cluster_id,
                                      cluster_col = "cluster",
                                      broad_name = "sctype_broad",
                                      fine_name = "sctype_fine") {

  clusters <- SingleCellExperiment::colData(sce_object)[[cluster_col]]
  cells_in_cluster <- which(clusters == cluster_id)

  if (length(cells_in_cluster) == 0) {
    stop(paste0("Cluster ", cluster_id, " not found"))
  }

  broad_cat <- unique(SingleCellExperiment::colData(sce_object)[[broad_name]][cells_in_cluster])
  fine_type <- unique(SingleCellExperiment::colData(sce_object)[[fine_name]][cells_in_cluster])

  # Check if annotation was refined (fine differs from broad)
  is_refined <- !identical(broad_cat, fine_type)

  list(
    cluster = cluster_id,
    broad_category = broad_cat,
    fine_subtype = fine_type,
    n_cells = length(cells_in_cluster),
    is_refined = is_refined
  )
}


#' Print hierarchical annotation table (SCE version)
#'
#' @param sce_object SingleCellExperiment object with hierarchical annotations
#' @param cluster_col Column name in colData containing cluster assignments
#' @param broad_name Column name for broad annotations
#' @param fine_name Column name for fine annotations
#'
#' @return Invisible data.frame with the table
#' @export
print_hierarchy_table_sce <- function(sce_object,
                                      cluster_col = "cluster",
                                      broad_name = "sctype_broad",
                                      fine_name = "sctype_fine") {

  clusters <- SingleCellExperiment::colData(sce_object)[[cluster_col]]
  unique_clusters <- sort(unique(clusters))

  results <- lapply(unique_clusters, function(cl) {
    info <- get_cluster_hierarchy_sce(sce_object, cl, cluster_col, broad_name, fine_name)
    data.frame(
      Cluster = info$cluster,
      Broad = info$broad_category,
      Fine = info$fine_subtype,
      NCells = info$n_cells,
      Refined = ifelse(info$is_refined, "Yes", "No"),
      stringsAsFactors = FALSE
    )
  })

  table_df <- do.call("rbind", results)

  cat("\n========================================\n")
  cat("  Hierarchical Cell Type Annotations\n")
  cat("========================================\n\n")
  print(table_df, row.names = FALSE)
  cat("\n")

  invisible(table_df)
}
