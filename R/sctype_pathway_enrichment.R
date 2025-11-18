# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/master/LICENSE)
# ScType pathway enrichment integration for Seurat objects
# Written by Claude (AI assistant), 2025-11-15

#' Add pathway enrichment-weighted scores to ScType annotations (Seurat)
#'
#' Integrates pathway enrichment and gene ontology analyses to validate and weight
#' ScType annotations. Runs differential expression per cluster, performs enrichment
#' using multiple tools (EnrichR, fgsea, GO), and combines pathway support with
#' ScType scores for improved confidence metrics.
#'
#' @param seurat_object Seurat object with clustering
#' @param known_tissue_type Tissue type (e.g., "Immune system", "Brain")
#' @param database_file Path to marker database (default: ScTypeDB_full.xlsx)
#' @param assay Assay to use (default: "RNA")
#' @param cluster_col Column with cluster assignments (default: "seurat_clusters")
#' @param enrichment_tools Vector of tools: "enrichr", "fgsea", "go" (default: all)
#' @param top_n_genes Number of top DE genes per cluster for enrichment (default: 200)
#' @param min_pct Minimum percentage of cells expressing gene (default: 0.25)
#' @param logfc_threshold Log fold-change threshold for DE (default: 0.25)
#' @param annotation_prefix Prefix for new metadata columns (default: "sctype")
#'
#' @return Modified Seurat object with pathway-weighted annotations
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- add_pathway_weighted_scores(
#'     seurat_obj,
#'     known_tissue_type = "Immune system",
#'     enrichment_tools = c("enrichr", "go")
#' )
#' # New columns: sctype_pathway_score, sctype_pathway_support,
#' #              sctype_combined_confidence, sctype_top_pathways
#' }
add_pathway_weighted_scores <- function(seurat_object,
                                       known_tissue_type = NULL,
                                       database_file = NULL,
                                       assay = "RNA",
                                       cluster_col = "seurat_clusters",
                                       enrichment_tools = c("enrichr", "fgsea", "go"),
                                       top_n_genes = 200,
                                       min_pct = 0.25,
                                       logfc_threshold = 0.25,
                                       annotation_prefix = "sctype") {

  # Load required packages
  required_pkgs <- c("Seurat", "dplyr")
  lapply(required_pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' is required but not installed."))
    }
  })

  library(Seurat)
  library(dplyr)

  # Check for enrichment tool packages
  if ("enrichr" %in% enrichment_tools && !requireNamespace("enrichR", quietly = TRUE)) {
    message("enrichR not installed. Removing from enrichment_tools.")
    enrichment_tools <- setdiff(enrichment_tools, "enrichr")
  }
  if ("fgsea" %in% enrichment_tools && !requireNamespace("fgsea", quietly = TRUE)) {
    message("fgsea not installed. Removing from enrichment_tools.")
    enrichment_tools <- setdiff(enrichment_tools, "fgsea")
  }
  if ("go" %in% enrichment_tools) {
    if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
        !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      message("clusterProfiler or org.Hs.eg.db not installed. Removing 'go' from tools.")
      enrichment_tools <- setdiff(enrichment_tools, "go")
    }
  }

  if (length(enrichment_tools) == 0) {
    stop("No enrichment tools available. Please install: enrichR, fgsea, and/or clusterProfiler")
  }

  # First, run standard ScType uncertainty scoring
  message("Step 1/4: Running ScType uncertainty scoring...")
  # ScType uncertainty functions are available from package namespace
  # Function available: add_sctype_uncertainty

  seurat_object <- add_sctype_uncertainty(
    seurat_object,
    known_tissue_type = known_tissue_type,
    database_file = database_file,
    assay = assay,
    cluster_col = cluster_col,
    top_n = 3,
    annotation_prefix = annotation_prefix
  )

  # Step 2: Find cluster markers
  message("Step 2/4: Identifying cluster markers via differential expression...")
  DefaultAssay(seurat_object) <- assay
  Idents(seurat_object) <- cluster_col

  cluster_markers <- FindAllMarkers(
    seurat_object,
    only.pos = TRUE,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    verbose = FALSE
  )

  # Step 3: Run enrichment analyses
  message("Step 3/4: Running pathway enrichment analyses...")
  pathway_results <- run_enrichment_per_cluster(
    cluster_markers,
    enrichment_tools = enrichment_tools,
    top_n_genes = top_n_genes
  )

  # Step 4: Match pathways to cell types and calculate support scores
  message("Step 4/4: Calculating pathway support for cell type annotations...")
  pathway_scores <- calculate_pathway_support(
    seurat_object,
    pathway_results,
    annotation_prefix = annotation_prefix
  )

  # Add pathway scores to metadata
  seurat_object <- add_pathway_scores_to_metadata(
    seurat_object,
    pathway_scores,
    annotation_prefix = annotation_prefix
  )

  message("Done! Added pathway-weighted annotations to metadata.")
  message("New columns:")
  message(paste0("  - ", annotation_prefix, "_pathway_score: Pathway enrichment support (0-1)"))
  message(paste0("  - ", annotation_prefix, "_pathway_support: Level of pathway support (High/Medium/Low)"))
  message(paste0("  - ", annotation_prefix, "_combined_confidence: ScType + pathway combined score"))
  message(paste0("  - ", annotation_prefix, "_top_pathways: Top enriched pathways"))

  # Store enrichment results
  attr(seurat_object, "sctype_pathway_results") <- pathway_results

  return(seurat_object)
}


#' Run enrichment analysis per cluster
#' @keywords internal
run_enrichment_per_cluster <- function(cluster_markers,
                                      enrichment_tools,
                                      top_n_genes = 200) {

  clusters <- unique(cluster_markers$cluster)
  enrichment_results <- list()

  for (cl in clusters) {
    message(paste0("  Enriching cluster ", cl, "..."))

    # Get top genes for this cluster
    cluster_genes <- cluster_markers %>%
      filter(cluster == cl) %>%
      arrange(desc(avg_log2FC)) %>%
      head(top_n_genes) %>%
      pull(gene)

    if (length(cluster_genes) == 0) {
      enrichment_results[[as.character(cl)]] <- NULL
      next
    }

    cluster_enrichment <- list()

    # EnrichR
    if ("enrichr" %in% enrichment_tools) {
      cluster_enrichment$enrichr <- run_enrichr(cluster_genes)
    }

    # fgsea
    if ("fgsea" %in% enrichment_tools) {
      # Get ranked genes
      ranked_genes <- cluster_markers %>%
        filter(cluster == cl) %>%
        arrange(desc(avg_log2FC)) %>%
        pull(avg_log2FC, gene)

      cluster_enrichment$fgsea <- run_fgsea(ranked_genes)
    }

    # GO enrichment
    if ("go" %in% enrichment_tools) {
      cluster_enrichment$go <- run_go_enrichment(cluster_genes)
    }

    enrichment_results[[as.character(cl)]] <- cluster_enrichment
  }

  return(enrichment_results)
}


#' Run EnrichR enrichment
#' @keywords internal
run_enrichr <- function(genes) {
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    return(NULL)
  }

  tryCatch({
    library(enrichR)

    # Use cell type-relevant databases
    dbs <- c(
      "CellMarker_Augmented_2021",
      "Azimuth_Cell_Types_2021",
      "PanglaoDB_Augmented_2021",
      "GO_Biological_Process_2021",
      "GO_Cellular_Component_2021",
      "KEGG_2021_Human",
      "Reactome_2022"
    )

    # Check which databases are available
    available_dbs <- listEnrichrDbs()$libraryName
    dbs <- dbs[dbs %in% available_dbs]

    if (length(dbs) == 0) {
      dbs <- c("GO_Biological_Process_2021", "KEGG_2021_Human")
    }

    enriched <- enrichr(genes, dbs)

    # Extract top results
    top_results <- lapply(enriched, function(df) {
      if (nrow(df) > 0) {
        df %>%
          filter(Adjusted.P.value < 0.05) %>%
          arrange(Adjusted.P.value) %>%
          head(10) %>%
          select(Term, Adjusted.P.value, Genes)
      } else {
        NULL
      }
    })

    return(top_results)
  }, error = function(e) {
    message(paste0("EnrichR error: ", e$message))
    return(NULL)
  })
}


#' Run fgsea enrichment
#' @keywords internal
run_fgsea <- function(ranked_genes) {
  if (!requireNamespace("fgsea", quietly = TRUE) ||
      !requireNamespace("msigdbr", quietly = TRUE)) {
    return(NULL)
  }

  tryCatch({
    library(fgsea)
    library(msigdbr)

    # Get gene sets
    # Use cell type signatures and hallmark pathways
    m_df <- msigdbr(species = "Homo sapiens", category = "H")  # Hallmark
    m_c8 <- msigdbr(species = "Homo sapiens", category = "C8")  # Cell type signatures

    pathways_h <- split(m_df$gene_symbol, m_df$gs_name)
    pathways_c8 <- split(m_c8$gene_symbol, m_c8$gs_name)
    pathways <- c(pathways_h, pathways_c8)

    # Run fgsea
    fgsea_res <- fgsea(
      pathways = pathways,
      stats = ranked_genes,
      minSize = 10,
      maxSize = 500
    )

    # Get significant results
    sig_results <- fgsea_res %>%
      filter(padj < 0.05) %>%
      arrange(padj) %>%
      head(20) %>%
      select(pathway, padj, NES, leadingEdge)

    return(sig_results)
  }, error = function(e) {
    message(paste0("fgsea error: ", e$message))
    return(NULL)
  })
}


#' Run GO enrichment
#' @keywords internal
run_go_enrichment <- function(genes) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    return(NULL)
  }

  tryCatch({
    library(clusterProfiler)
    library(org.Hs.eg.db)

    # Convert symbols to Entrez IDs
    gene_entrez <- bitr(
      genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )

    if (nrow(gene_entrez) == 0) {
      return(NULL)
    }

    # GO enrichment
    go_results <- enrichGO(
      gene = gene_entrez$ENTREZID,
      OrgDb = org.Hs.eg.db,
      ont = "BP",  # Biological Process
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )

    if (nrow(go_results) == 0) {
      return(NULL)
    }

    # Extract results
    go_df <- as.data.frame(go_results) %>%
      arrange(p.adjust) %>%
      head(10) %>%
      select(Description, p.adjust, GeneRatio, geneID)

    return(go_df)
  }, error = function(e) {
    message(paste0("GO enrichment error: ", e$message))
    return(NULL)
  })
}


#' Calculate pathway support for cell type annotations
#' @keywords internal
calculate_pathway_support <- function(seurat_object,
                                     pathway_results,
                                     annotation_prefix = "sctype") {

  # Get cluster assignments and ScType annotations
  clusters <- seurat_object@meta.data[[paste0(annotation_prefix, "_top1")]]
  cluster_ids <- seurat_object@meta.data[[attr(seurat_object, "sctype_uncertainty_clusters")$cluster[1] %>% class()]]

  # Cell type to pathway mapping (knowledge base)
  celltype_pathway_map <- build_celltype_pathway_map()

  # Calculate support scores
  pathway_scores <- data.frame(
    cluster = character(),
    pathway_score = numeric(),
    pathway_support = character(),
    top_pathways = character(),
    stringsAsFactors = FALSE
  )

  for (cl in names(pathway_results)) {
    enrichment <- pathway_results[[cl]]

    if (is.null(enrichment) || length(enrichment) == 0) {
      pathway_scores <- rbind(pathway_scores, data.frame(
        cluster = cl,
        pathway_score = 0,
        pathway_support = "None",
        top_pathways = "No enrichment",
        stringsAsFactors = FALSE
      ))
      next
    }

    # Get predicted cell type for this cluster
    cluster_celltype <- unique(seurat_object@meta.data[
      seurat_object@meta.data[[paste0("seurat_clusters")]] == cl,
      paste0(annotation_prefix, "_top1")
    ])[1]

    # Extract enriched terms
    all_terms <- extract_enriched_terms(enrichment)

    # Calculate pathway support score
    support_score <- calculate_support_score(
      all_terms,
      cluster_celltype,
      celltype_pathway_map
    )

    # Categorize support
    support_level <- ifelse(support_score >= 0.7, "High",
                           ifelse(support_score >= 0.4, "Medium", "Low"))

    # Get top pathways
    top_paths <- if (length(all_terms) > 0) {
      paste(head(all_terms, 3), collapse = "; ")
    } else {
      "No enrichment"
    }

    pathway_scores <- rbind(pathway_scores, data.frame(
      cluster = cl,
      pathway_score = support_score,
      pathway_support = support_level,
      top_pathways = top_paths,
      stringsAsFactors = FALSE
    ))
  }

  return(pathway_scores)
}


#' Build cell type to pathway mapping
#' @keywords internal
build_celltype_pathway_map <- function() {
  # Knowledge base mapping cell types to expected pathways/functions
  # This is a simplified version - could be expanded

  map <- list(
    # Immune cells
    "T cells" = c("T cell", "lymphocyte", "adaptive immune", "TCR", "CD3"),
    "CD4+ T cells" = c("CD4", "T helper", "Th1", "Th2", "Th17"),
    "CD8+ T cells" = c("CD8", "cytotoxic T", "CTL", "granzyme"),
    "B cells" = c("B cell", "antibody", "immunoglobulin", "BCR", "CD19"),
    "NK cells" = c("natural killer", "NK cell", "cytotoxicity", "NKG"),
    "Monocytes" = c("monocyte", "CD14", "CD16", "myeloid"),
    "Macrophages" = c("macrophage", "phagocytosis", "inflammation"),
    "Dendritic" = c("dendritic", "antigen presenting", "MHC"),

    # Neurons
    "Neurons" = c("neuron", "synapse", "axon", "dendrite", "neurotransmitter"),
    "Excitatory" = c("glutamate", "excitatory", "AMPA", "NMDA"),
    "Inhibitory" = c("GABA", "inhibitory", "interneuron"),
    "Dopaminergic" = c("dopamine", "TH", "dopaminergic"),

    # Glia
    "Astrocytes" = c("astrocyte", "glial", "GFAP"),
    "Oligodendrocytes" = c("oligodendrocyte", "myelin", "MBP", "MOG"),
    "Microglia" = c("microglia", "immune", "CX3CR1"),

    # Other
    "Endothelial" = c("endothelial", "vascular", "angiogenesis", "PECAM"),
    "Fibroblasts" = c("fibroblast", "collagen", "ECM", "connective"),
    "Hepatocytes" = c("hepatocyte", "liver", "metabolism", "cytochrome"),
    "Cardiomyocytes" = c("cardiomyocyte", "heart", "cardiac", "myosin")
  )

  return(map)
}


#' Extract enriched terms from all enrichment results
#' @keywords internal
extract_enriched_terms <- function(enrichment) {
  all_terms <- c()

  # From EnrichR
  if (!is.null(enrichment$enrichr)) {
    for (db in names(enrichment$enrichr)) {
      if (!is.null(enrichment$enrichr[[db]])) {
        all_terms <- c(all_terms, enrichment$enrichr[[db]]$Term)
      }
    }
  }

  # From fgsea
  if (!is.null(enrichment$fgsea) && nrow(enrichment$fgsea) > 0) {
    all_terms <- c(all_terms, enrichment$fgsea$pathway)
  }

  # From GO
  if (!is.null(enrichment$go) && nrow(enrichment$go) > 0) {
    all_terms <- c(all_terms, enrichment$go$Description)
  }

  return(unique(all_terms))
}


#' Calculate pathway support score
#' @keywords internal
calculate_support_score <- function(enriched_terms,
                                   celltype,
                                   celltype_pathway_map) {

  if (length(enriched_terms) == 0 || is.na(celltype) || celltype == "Unknown") {
    return(0)
  }

  # Find matching cell type in map (fuzzy matching)
  matching_keywords <- NULL
  for (ct in names(celltype_pathway_map)) {
    if (grepl(ct, celltype, ignore.case = TRUE) ||
        grepl(celltype, ct, ignore.case = TRUE)) {
      matching_keywords <- celltype_pathway_map[[ct]]
      break
    }
  }

  if (is.null(matching_keywords)) {
    # No specific mapping, return neutral score
    return(0.5)
  }

  # Count matches
  matches <- 0
  for (term in enriched_terms) {
    for (keyword in matching_keywords) {
      if (grepl(keyword, term, ignore.case = TRUE)) {
        matches <- matches + 1
        break  # Count each term once
      }
    }
  }

  # Calculate score (proportion of enriched terms supporting cell type)
  score <- min(matches / max(length(enriched_terms), 5), 1.0)

  return(score)
}


#' Add pathway scores to metadata
#' @keywords internal
add_pathway_scores_to_metadata <- function(seurat_object,
                                          pathway_scores,
                                          annotation_prefix = "sctype") {

  n_cells <- ncol(seurat_object)

  # Initialize columns
  seurat_object@meta.data[[paste0(annotation_prefix, "_pathway_score")]] <- rep(NA, n_cells)
  seurat_object@meta.data[[paste0(annotation_prefix, "_pathway_support")]] <- rep(NA, n_cells)
  seurat_object@meta.data[[paste0(annotation_prefix, "_combined_confidence")]] <- rep(NA, n_cells)
  seurat_object@meta.data[[paste0(annotation_prefix, "_top_pathways")]] <- rep(NA, n_cells)

  # Get cluster assignments
  clusters <- seurat_object@meta.data$seurat_clusters

  # Fill in values
  for (i in 1:nrow(pathway_scores)) {
    cl <- pathway_scores$cluster[i]
    cells_in_cluster <- which(as.character(clusters) == cl)

    seurat_object@meta.data[[paste0(annotation_prefix, "_pathway_score")]][cells_in_cluster] <-
      pathway_scores$pathway_score[i]
    seurat_object@meta.data[[paste0(annotation_prefix, "_pathway_support")]][cells_in_cluster] <-
      pathway_scores$pathway_support[i]
    seurat_object@meta.data[[paste0(annotation_prefix, "_top_pathways")]][cells_in_cluster] <-
      pathway_scores$top_pathways[i]

    # Calculate combined confidence (weighted average)
    sctype_conf <- seurat_object@meta.data[[paste0(annotation_prefix, "_confidence")]][cells_in_cluster[1]]
    pathway_score <- pathway_scores$pathway_score[i]

    # Weight: 60% ScType, 40% pathway
    combined_conf <- 0.6 * sctype_conf + 0.4 * pathway_score

    seurat_object@meta.data[[paste0(annotation_prefix, "_combined_confidence")]][cells_in_cluster] <-
      combined_conf
  }

  return(seurat_object)
}


#' Visualize pathway support for annotations
#'
#' Creates visualizations showing pathway enrichment support for cell type annotations
#'
#' @param seurat_object Seurat object with pathway-weighted scores
#' @param annotation_prefix Prefix for annotation columns (default: "sctype")
#' @param save_plots Save plots to files (default: FALSE)
#' @param output_dir Output directory (default: "pathway_plots")
#'
#' @return List of ggplot objects
#' @export
visualize_pathway_support <- function(seurat_object,
                                     annotation_prefix = "sctype",
                                     save_plots = FALSE,
                                     output_dir = "pathway_plots") {

  library(ggplot2)
  library(dplyr)

  # Check if pathway scores exist
  if (!paste0(annotation_prefix, "_pathway_score") %in% colnames(seurat_object@meta.data)) {
    stop("Pathway scores not found. Run add_pathway_weighted_scores() first.")
  }

  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  plot_list <- list()

  # Plot 1: ScType vs Pathway confidence
  p1 <- ggplot(seurat_object@meta.data,
               aes(x = !!sym(paste0(annotation_prefix, "_confidence")),
                   y = !!sym(paste0(annotation_prefix, "_pathway_score")),
                   color = !!sym(paste0(annotation_prefix, "_top1")))) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    labs(title = "ScType vs Pathway Enrichment Confidence",
         subtitle = "Points above diagonal have stronger pathway support",
         x = "ScType Confidence", y = "Pathway Support Score",
         color = "Cell Type") +
    theme_minimal() +
    theme(legend.position = "right")

  plot_list$sctype_vs_pathway <- p1

  # Plot 2: Combined confidence UMAP
  if ("umap" %in% names(seurat_object@reductions)) {
    p2 <- FeaturePlot(seurat_object,
                     features = paste0(annotation_prefix, "_combined_confidence"),
                     reduction = "umap") +
      scale_color_gradient2(low = "blue", mid = "yellow", high = "red",
                           midpoint = 0.5) +
      ggtitle("Combined Confidence (ScType + Pathway)")

    plot_list$combined_umap <- p2
  }

  # Plot 3: Pathway support levels
  support_data <- seurat_object@meta.data %>%
    select(!!sym(paste0(annotation_prefix, "_pathway_support"))) %>%
    rename(support = 1)

  support_data$support <- factor(support_data$support,
                                 levels = c("High", "Medium", "Low", "None"))

  p3 <- ggplot(support_data, aes(x = support, fill = support)) +
    geom_bar() +
    scale_fill_manual(values = c("High" = "#2ECC40", "Medium" = "#FFDC00",
                                 "Low" = "#FF4136", "None" = "#AAAAAA")) +
    labs(title = "Pathway Support Levels",
         x = "Pathway Support", y = "Number of Cells") +
    theme_minimal() +
    theme(legend.position = "none")

  plot_list$support_levels <- p3

  if (save_plots) {
    ggsave(file.path(output_dir, "sctype_vs_pathway.png"), p1,
           width = 12, height = 8, dpi = 300)
    if (!is.null(plot_list$combined_umap)) {
      ggsave(file.path(output_dir, "combined_confidence_umap.png"), p2,
             width = 10, height = 8, dpi = 300)
    }
    ggsave(file.path(output_dir, "pathway_support_levels.png"), p3,
           width = 10, height = 6, dpi = 300)
  }

  invisible(plot_list)
}
