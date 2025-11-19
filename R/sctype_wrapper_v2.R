# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/master/LICENSE)
# ScType v2 Wrapper - Statistical Testing Integration
# Written by Claude Code, November 2025
#
# This is the v2 wrapper that integrates statistical significance testing
# to replace the arbitrary ncells/4 threshold with FDR-corrected p-values.

#' @title Load ScType v2 source files
#' @name sctype_source_v2
#' @description Loads ScType v2 functions including statistical testing
#' @details Loads all necessary functions for automated cell type annotation with statistical validation
#' @param none
#' @return URL to original ScType database
#' @export
#' @examples
#' db_ <- sctype_source_v2()
#'
sctype_source_v2 <- function(){
    # Load v2 statistical testing functions
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_statistics.R")
    # Load tissue auto detect
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
    # Load gene set preparation function
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
    # Load cell type annotation function (with return_details parameter)
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
    # Load ScType database
    db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
    return(db_)
}


#' @title Run ScType v2 analysis with statistical testing
#' @name run_sctype_v2
#' @description Run automated cell type annotation with statistical significance testing
#' @details This function replaces the arbitrary ncells/4 threshold with FDR-corrected p-values.
#'          Compatible with Seurat v4 and v5.
#'
#' @param seurat_object A Seurat object with clustering completed
#' @param known_tissue_type Tissue type (e.g., "Immune system", "Brain"). Auto-detected if NULL (default: NULL)
#' @param assay Assay name to use (default: "RNA")
#' @param scaled Whether data is scaled (default: TRUE)
#' @param custom_marker_file Path to custom marker database (default: NULL, uses ScTypeDB_full.xlsx)
#' @param plot Generate UMAP plot with annotations (default: FALSE)
#' @param name Metadata column name for cell type annotations (default: "sctype_v2_classification")
#' @param fdr_threshold FDR significance threshold (default: 0.05)
#'   - 0.01: Very stringent, high confidence only
#'   - 0.05: Standard significance level (recommended)
#'   - 0.1: More lenient, includes medium confidence
#' @param use_permutation Use permutation testing for empirical p-values (default: FALSE)
#'   - FALSE: Fast parametric z-score testing
#'   - TRUE: Slower but more robust permutation testing (requires n_permutations parameter)
#' @param n_permutations Number of permutations if use_permutation = TRUE (default: 1000)
#' @param cluster_col Metadata column with cluster assignments (default: "seurat_clusters")
#' @param top_n Number of top candidate cell types to report per cluster (default: 3)
#' @param verbose Print detailed messages (default: TRUE)
#'
#' @return Modified Seurat object with new metadata columns:
#'   - [name]: Cell type annotations (Unknown if FDR >= threshold)
#'   - [name]_score: Raw ScType scores for top candidate
#'   - [name]_zscore: Z-scores for statistical testing
#'   - [name]_pvalue: Uncorrected p-values
#'   - [name]_fdr: FDR-corrected p-values
#'   - [name]_confidence: Confidence level (High/Medium/Low/Very Low)
#'   - [name]_top1, _top2, ..., _topN: Top N cell type candidates
#'   - [name]_score1, _score2, ..., _scoreN: Scores for each top N candidate
#'
#' @import dplyr
#' @import Seurat
#' @export
#'
#' @examples
#' # Basic usage (auto-detect tissue)
#' seurat_obj <- run_sctype_v2(seurat_obj)
#'
#' # Specify tissue type and threshold
#' seurat_obj <- run_sctype_v2(seurat_obj,
#'                             known_tissue_type = "Immune system",
#'                             fdr_threshold = 0.01)
#'
#' # Get top 5 candidates instead of default 3
#' seurat_obj <- run_sctype_v2(seurat_obj,
#'                             known_tissue_type = "Immune system",
#'                             top_n = 5)
#'
#' # View top candidates for cluster 0
#' cluster_0_cells <- seurat_obj@meta.data[seurat_obj$seurat_clusters == 0, ]
#' head(cluster_0_cells[, c("sctype_v2_top1", "sctype_v2_score1",
#'                          "sctype_v2_top2", "sctype_v2_score2",
#'                          "sctype_v2_top3", "sctype_v2_score3")], 1)
#'
#' # Use permutation testing (slower but robust)
#' seurat_obj <- run_sctype_v2(seurat_obj,
#'                             use_permutation = TRUE,
#'                             n_permutations = 1000)
#'
run_sctype_v2 <- function(seurat_object,
                          known_tissue_type = NULL,
                          assay = "RNA",
                          scaled = TRUE,
                          custom_marker_file = NULL,
                          plot = FALSE,
                          name = "sctype_v2_classification",
                          fdr_threshold = 0.05,
                          use_permutation = FALSE,
                          n_permutations = 1000,
                          cluster_col = "seurat_clusters",
                          top_n = 3,
                          verbose = TRUE) {

    # Load functions
    db_ <- sctype_source_v2()

    # Input validation
    if (is.null(seurat_object)) {
        stop("Argument 'seurat_object' is missing")
    }
    if (!inherits(seurat_object, "Seurat")) {
        stop("Argument 'seurat_object' must be a Seurat object")
    }
    if (!cluster_col %in% colnames(seurat_object@meta.data)) {
        stop(sprintf("Cluster column '%s' not found in metadata. Available columns: %s",
                     cluster_col, paste(colnames(seurat_object@meta.data), collapse = ", ")))
    }

    # Set default custom marker file
    if (is.null(custom_marker_file)) {
        custom_marker_file <- db_
    }

    # Auto-detect tissue type if not provided
    if (is.null(known_tissue_type)) {
        if (verbose) message("Auto-detecting tissue type...")
        tissue_type <- auto_detect_tissue_type(path_to_db_file = custom_marker_file,
                                              seuratObject = seurat_object,
                                              scaled = scaled, assay = assay)
        rownames(tissue_type) <- NULL
        tissue_type <- tissue_type$tissue[1]
        if (verbose) message(sprintf("  Detected tissue type: %s", tissue_type))
    } else {
        tissue_type <- known_tissue_type
        if (verbose) message(sprintf("Using specified tissue type: %s", tissue_type))
    }

    # Prepare gene sets
    if (verbose) message("Preparing gene sets...")
    gs_list <- gene_sets_prepare(custom_marker_file, tissue_type)

    # Determine Seurat version and extract data
    data_type <- if (scaled) "scale.data" else "counts"
    package_type <- data_type %in% names(attributes(seurat_object[[assay]]))

    # Calculate ScType scores with details
    if (verbose) message("Calculating ScType scores...")

    if (package_type) {
        if (verbose) message("  Using Seurat v4 object")
        scRNAseqData <- slot(seurat_object[[assay]], data_type)
    } else {
        if (verbose) message("  Using Seurat v5 object")
        if (data_type == "scale.data") {
            scRNAseqData <- seurat_object[[assay]]$scale.data
        } else {
            scRNAseqData <- seurat_object[[assay]]$counts
        }
    }

    # Get scores WITH details for statistical testing
    score_results <- sctype_score(scRNAseqData = as.matrix(scRNAseqData),
                                  scaled = scaled,
                                  gs = gs_list$gs_positive,
                                  gs2 = gs_list$gs_negative,
                                  return_details = TRUE)

    es.max <- score_results$scores

    # Extract cluster assignments
    cluster_assignments <- seurat_object@meta.data[[cluster_col]]

    # Calculate aggregate scores per cluster
    if (verbose) message("Aggregating scores by cluster...")
    cluster_results <- do.call("rbind", lapply(unique(cluster_assignments), function(cl) {
        # Get cells in this cluster
        cells_in_cluster <- rownames(seurat_object@meta.data[cluster_assignments == cl, ])
        ncells <- length(cells_in_cluster)

        # Sum scores across cells in cluster
        cluster_scores <- rowSums(es.max[, cells_in_cluster, drop = FALSE])

        # Sort scores to get top N
        sorted_scores <- sort(cluster_scores, decreasing = TRUE)
        n_available <- min(top_n, length(sorted_scores))

        # Create base data frame
        result <- data.frame(
            cluster = cl,
            top_celltype = names(sorted_scores)[1],
            top_score = sorted_scores[1],
            ncells = ncells,
            stringsAsFactors = FALSE
        )

        # Add top N cell types and scores
        for (i in 1:top_n) {
            if (i <= n_available) {
                result[[paste0("top", i, "_celltype")]] <- names(sorted_scores)[i]
                result[[paste0("top", i, "_score")]] <- sorted_scores[i]
            } else {
                result[[paste0("top", i, "_celltype")]] <- "None"
                result[[paste0("top", i, "_score")]] <- 0
            }
        }

        # Calculate score differences (useful for ambiguity detection)
        if (n_available >= 2) {
            result$score_difference <- sorted_scores[1] - sorted_scores[2]
        } else {
            result$score_difference <- sorted_scores[1]
        }

        return(result)
    }))

    # Calculate statistical significance
    if (use_permutation) {
        # Permutation-based testing
        if (verbose) message(sprintf("Running permutation testing (%d permutations)...", n_permutations))
        null_dist <- generate_null_distribution(
            sctype_scores_matrix = es.max,
            cluster_assignments = cluster_assignments,
            n_permutations = n_permutations
        )

        # Add empirical p-values to cluster_results
        cluster_results$pvalue <- null_dist$empirical_pvalues[as.character(cluster_results$cluster)]

        # Calculate z-scores from observed vs null
        cluster_results$zscore <- sapply(seq_len(nrow(cluster_results)), function(i) {
            obs <- cluster_results$top_score[i]
            null_scores <- null_dist$null_scores[, as.character(cluster_results$cluster[i])]
            (obs - mean(null_scores)) / sd(null_scores)
        })

    } else {
        # Parametric z-score testing (faster)
        if (verbose) message("Calculating z-scores and p-values...")

        # Create input for calculate_zscore_pvalue
        scores_vec <- setNames(cluster_results$top_score, cluster_results$cluster)
        cluster_sizes <- setNames(cluster_results$ncells, cluster_results$cluster)

        stat_results <- calculate_zscore_pvalue(scores_vec, cluster_sizes, use_global_stats = TRUE)

        # Merge statistical results
        cluster_results$zscore <- stat_results$zscore
        cluster_results$pvalue <- stat_results$pvalue
    }

    # Apply FDR correction
    if (verbose) message("Applying FDR correction...")
    cluster_results <- apply_fdr_correction(cluster_results, method = "BH")

    # Assign confidence levels
    if (verbose) message("Assigning confidence levels...")
    cluster_results <- assign_confidence_level(
        cluster_results,
        fdr_thresholds = c("High" = 0.01, "Medium" = 0.05, "Low" = 0.1),
        use_combined = FALSE
    )

    # Assign cell types based on FDR threshold
    cluster_results$assigned_celltype <- ifelse(
        cluster_results$fdr < fdr_threshold,
        cluster_results$top_celltype,
        "Unknown"
    )

    # Add results to Seurat object metadata
    if (verbose) message("Adding results to Seurat object...")
    seurat_object_res <- seurat_object
    seurat_object_res@meta.data[[name]] <- ""
    seurat_object_res@meta.data[[paste0(name, "_score")]] <- NA
    seurat_object_res@meta.data[[paste0(name, "_zscore")]] <- NA
    seurat_object_res@meta.data[[paste0(name, "_pvalue")]] <- NA
    seurat_object_res@meta.data[[paste0(name, "_fdr")]] <- NA
    seurat_object_res@meta.data[[paste0(name, "_confidence")]] <- ""

    # Initialize top N columns
    for (i in 1:top_n) {
        seurat_object_res@meta.data[[paste0(name, "_top", i)]] <- ""
        seurat_object_res@meta.data[[paste0(name, "_score", i)]] <- NA
    }

    for (j in unique(cluster_results$cluster)) {
        cl_result <- cluster_results[cluster_results$cluster == j, ]
        cluster_cells <- seurat_object_res@meta.data[[cluster_col]] == j

        seurat_object_res@meta.data[cluster_cells, name] <- cl_result$assigned_celltype
        seurat_object_res@meta.data[cluster_cells, paste0(name, "_score")] <- cl_result$top_score
        seurat_object_res@meta.data[cluster_cells, paste0(name, "_zscore")] <- cl_result$zscore
        seurat_object_res@meta.data[cluster_cells, paste0(name, "_pvalue")] <- cl_result$pvalue
        seurat_object_res@meta.data[cluster_cells, paste0(name, "_fdr")] <- cl_result$fdr
        seurat_object_res@meta.data[cluster_cells, paste0(name, "_confidence")] <- as.character(cl_result$confidence)

        # Add top N candidates and scores
        for (i in 1:top_n) {
            seurat_object_res@meta.data[cluster_cells, paste0(name, "_top", i)] <- cl_result[[paste0("top", i, "_celltype")]]
            seurat_object_res@meta.data[cluster_cells, paste0(name, "_score", i)] <- cl_result[[paste0("top", i, "_score")]]
        }
    }

    # Print summary
    if (verbose) {
        message("\n===== ScType v2 Annotation Summary =====")
        message(sprintf("Tissue type: %s", tissue_type))
        message(sprintf("FDR threshold: %.3f", fdr_threshold))
        message(sprintf("Total clusters: %d", nrow(cluster_results)))
        message(sprintf("  High confidence (FDR < 0.01): %d clusters", sum(cluster_results$fdr < 0.01)))
        message(sprintf("  Medium confidence (FDR 0.01-0.05): %d clusters", sum(cluster_results$fdr >= 0.01 & cluster_results$fdr < 0.05)))
        message(sprintf("  Low confidence (FDR 0.05-0.1): %d clusters", sum(cluster_results$fdr >= 0.05 & cluster_results$fdr < 0.1)))
        message(sprintf("  Unknown (FDR >= %.2f): %d clusters", fdr_threshold, sum(cluster_results$fdr >= fdr_threshold)))
        message("\nNew metadata columns added:")
        message(sprintf("  - %s (cell type annotations)", name))
        message(sprintf("  - %s_score (raw ScType scores)", name))
        message(sprintf("  - %s_zscore (z-scores)", name))
        message(sprintf("  - %s_pvalue (uncorrected p-values)", name))
        message(sprintf("  - %s_fdr (FDR-corrected p-values)", name))
        message(sprintf("  - %s_confidence (High/Medium/Low/Very Low)", name))
        message(sprintf("  - %s_top1, _top2, ..., _top%d (top %d cell type candidates)", name, top_n, top_n))
        message(sprintf("  - %s_score1, _score2, ..., _score%d (scores for top %d)", name, top_n, top_n))
        message("========================================\n")
    }

    # Generate plot if requested
    if (plot) {
        if (verbose) message("Generating UMAP plot...")
        plot_ <- DimPlot(seurat_object_res, reduction = "umap", group.by = name, label = TRUE)
        print(plot_)
    }

    # Store detailed results as attribute
    attr(seurat_object_res, "sctype_v2_statistics") <- cluster_results
    attr(seurat_object_res, "sctype_v2_params") <- list(
        tissue_type = tissue_type,
        fdr_threshold = fdr_threshold,
        use_permutation = use_permutation,
        n_permutations = if (use_permutation) n_permutations else NA,
        timestamp = Sys.time()
    )

    return(seurat_object_res)
}


#' @title Get ScType v2 detailed statistics
#' @name get_sctype_v2_stats
#' @description Extract detailed cluster-level statistics from ScType v2 results
#' @param seurat_object Seurat object annotated with run_sctype_v2()
#' @return Data frame with cluster-level statistics
#' @export
#'
#' @examples
#' stats <- get_sctype_v2_stats(seurat_obj)
#' print(stats)
#'
get_sctype_v2_stats <- function(seurat_object) {
    if (!inherits(seurat_object, "Seurat")) {
        stop("Input must be a Seurat object")
    }

    stats <- attr(seurat_object, "sctype_v2_statistics")
    if (is.null(stats)) {
        stop("No ScType v2 statistics found. Run run_sctype_v2() first.")
    }

    return(stats)
}


#' @title Compare ScType v1 vs v2 annotations
#' @name compare_sctype_versions
#' @description Compare annotations from ScType v1 (ncells/4 threshold) and v2 (FDR)
#' @param seurat_object Seurat object with both v1 and v2 annotations
#' @param v1_col Metadata column with v1 annotations (default: "sctype_classification")
#' @param v2_col Metadata column with v2 annotations (default: "sctype_v2_classification")
#' @return List with comparison statistics and confusion matrix
#' @export
#'
compare_sctype_versions <- function(seurat_object,
                                   v1_col = "sctype_classification",
                                   v2_col = "sctype_v2_classification") {

    if (!all(c(v1_col, v2_col) %in% colnames(seurat_object@meta.data))) {
        stop(sprintf("Columns '%s' and '%s' must exist in metadata", v1_col, v2_col))
    }

    v1_annot <- seurat_object@meta.data[[v1_col]]
    v2_annot <- seurat_object@meta.data[[v2_col]]

    # Agreement statistics
    agreement <- mean(v1_annot == v2_annot, na.rm = TRUE)
    both_unknown <- mean(v1_annot == "Unknown" & v2_annot == "Unknown", na.rm = TRUE)
    v1_only_unknown <- mean(v1_annot == "Unknown" & v2_annot != "Unknown", na.rm = TRUE)
    v2_only_unknown <- mean(v1_annot != "Unknown" & v2_annot == "Unknown", na.rm = TRUE)

    # Confusion matrix
    conf_mat <- table(v1 = v1_annot, v2 = v2_annot)

    # Summary
    result <- list(
        agreement_rate = agreement,
        both_unknown_rate = both_unknown,
        v1_only_unknown_rate = v1_only_unknown,
        v2_only_unknown_rate = v2_only_unknown,
        confusion_matrix = conf_mat
    )

    class(result) <- c("sctype_version_comparison", "list")
    return(result)
}


#' Print method for sctype_version_comparison
#' @export
print.sctype_version_comparison <- function(x, ...) {
    cat("ScType v1 vs v2 Comparison\n")
    cat("==========================\n\n")
    cat(sprintf("Overall agreement: %.1f%%\n", x$agreement_rate * 100))
    cat(sprintf("Both versions 'Unknown': %.1f%%\n", x$both_unknown_rate * 100))
    cat(sprintf("v1 'Unknown', v2 annotated: %.1f%%\n", x$v1_only_unknown_rate * 100))
    cat(sprintf("v1 annotated, v2 'Unknown': %.1f%%\n\n", x$v2_only_unknown_rate * 100))

    cat("Confusion Matrix (top 10x10):\n")
    mat <- x$confusion_matrix
    if (nrow(mat) > 10 || ncol(mat) > 10) {
        print(mat[1:min(10, nrow(mat)), 1:min(10, ncol(mat))])
        cat(sprintf("\n(Showing top 10x10 of %dx%d matrix)\n", nrow(mat), ncol(mat)))
    } else {
        print(mat)
    }
}


#' Get top N candidates for a specific cluster
#'
#' @description Extract top N cell type candidates and their scores for a given cluster
#'
#' @param seurat_object Seurat object annotated with run_sctype_v2()
#' @param cluster_id Cluster ID to query
#' @param annotation_prefix Annotation prefix (default: "sctype_v2")
#' @param cluster_col Cluster column name (default: "seurat_clusters")
#'
#' @return Data frame with top N candidates and scores
#' @export
#'
#' @examples
#' # Get top candidates for cluster 0
#' candidates <- get_top_candidates(seurat_obj, cluster_id = 0)
#' print(candidates)
#'
get_top_candidates <- function(seurat_object,
                              cluster_id,
                              annotation_prefix = "sctype_v2",
                              cluster_col = "seurat_clusters") {

    # Get one cell from this cluster to extract values
    cluster_cells <- which(seurat_object@meta.data[[cluster_col]] == cluster_id)

    if (length(cluster_cells) == 0) {
        stop(sprintf("Cluster %s not found in column '%s'", cluster_id, cluster_col))
    }

    cell_data <- seurat_object@meta.data[cluster_cells[1], ]

    # Find all top N columns
    top_cols <- grep(paste0("^", annotation_prefix, "_top[0-9]+$"), colnames(cell_data), value = TRUE)
    score_cols <- grep(paste0("^", annotation_prefix, "_score[0-9]+$"), colnames(cell_data), value = TRUE)

    if (length(top_cols) == 0) {
        stop(sprintf("No top N columns found with prefix '%s'. Did you run run_sctype_v2()?", annotation_prefix))
    }

    # Extract rankings
    n_candidates <- length(top_cols)
    result <- data.frame(
        rank = 1:n_candidates,
        cell_type = character(n_candidates),
        score = numeric(n_candidates),
        stringsAsFactors = FALSE
    )

    for (i in 1:n_candidates) {
        result$cell_type[i] <- as.character(cell_data[[top_cols[i]]])
        result$score[i] <- as.numeric(cell_data[[score_cols[i]]])
    }

    # Add metadata
    result$cluster <- cluster_id
    result$n_cells <- length(cluster_cells)
    result$fdr <- as.numeric(cell_data[[paste0(annotation_prefix, "_fdr")]])
    result$confidence <- as.character(cell_data[[paste0(annotation_prefix, "_confidence")]])

    # Reorder columns
    result <- result[, c("cluster", "rank", "cell_type", "score", "n_cells", "fdr", "confidence")]

    # Filter out "None" entries
    result <- result[result$cell_type != "None", ]

    return(result)
}


#' Print summary of top candidates for all clusters
#'
#' @description Display top candidates for all clusters in a formatted table
#'
#' @param seurat_object Seurat object annotated with run_sctype_v2()
#' @param annotation_prefix Annotation prefix (default: "sctype_v2")
#' @param cluster_col Cluster column name (default: "seurat_clusters")
#' @param top_n_display Number of candidates to display per cluster (default: 3)
#'
#' @return Invisibly returns data frame with all candidates
#' @export
#'
#' @examples
#' # Print summary for all clusters
#' print_top_candidates_summary(seurat_obj)
#'
print_top_candidates_summary <- function(seurat_object,
                                        annotation_prefix = "sctype_v2",
                                        cluster_col = "seurat_clusters",
                                        top_n_display = 3) {

    all_clusters <- unique(seurat_object@meta.data[[cluster_col]])
    all_results <- list()

    cat("ScType v2 Top Candidates Summary\n")
    cat("=================================\n\n")

    for (cl in sort(all_clusters)) {
        candidates <- get_top_candidates(seurat_object, cl, annotation_prefix, cluster_col)
        all_results[[as.character(cl)]] <- candidates

        # Display top N
        display_candidates <- head(candidates, top_n_display)

        cat(sprintf("Cluster %s (n=%d cells, FDR=%.4f, %s confidence):\n",
                   cl, candidates$n_cells[1], candidates$fdr[1], candidates$confidence[1]))

        for (i in 1:nrow(display_candidates)) {
            cat(sprintf("  %d. %-30s (score: %.1f)\n",
                       display_candidates$rank[i],
                       display_candidates$cell_type[i],
                       display_candidates$score[i]))
        }
        cat("\n")
    }

    invisible(do.call("rbind", all_results))
}
