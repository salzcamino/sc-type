# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Hierarchical cell type annotation for ScType
# Provides both broad and fine-grained cell type annotations

#' @title Run hierarchical sctype annotation
#' @name run_sctype_hierarchical
#' @description Run automated cell type annotation at both broad and fine levels
#' @details Provides two-level hierarchical annotation: broad categories and fine subtypes
#'
#' @param seurat_object A Seurat object
#' @param known_tissue_type The tissue type of the input data (optional)
#' @param assay Assay name (default: "RNA")
#' @param scaled Use scaled data (default: TRUE)
#' @param custom_marker_file Path to hierarchical marker database (default: ScTypeDB_hierarchical.xlsx)
#' @param plot Whether to plot the results (default: FALSE)
#' @param broad_name Metadata column name for broad categories (default: "sctype_broad")
#' @param fine_name Metadata column name for fine subtypes (default: "sctype_fine")
#'
#' @return Seurat object with two new metadata columns: broad and fine annotations
#'
#' @export

run_sctype_hierarchical <- function(seurat_object,
                                    known_tissue_type = NULL,
                                    assay = "RNA",
                                    scaled = TRUE,
                                    custom_marker_file = NULL,
                                    plot = FALSE,
                                    broad_name = "sctype_broad",
                                    fine_name = "sctype_fine") {

    # Load required functions
    if(!exists("gene_sets_prepare")) {
        source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
    }
    if(!exists("sctype_score")) {
        source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
    }
    if(!exists("auto_detect_tissue_type")) {
        source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
    }

    # Check for missing arguments
    if (is.null(seurat_object)) {
        stop("Argument 'seurat_object' is missing")
    }
    if (!inherits(seurat_object, "Seurat")) {
        stop("Argument 'seurat_object' must be a Seurat object")
    }

    # Set default database
    if (is.null(custom_marker_file)) {
        custom_marker_file = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_hierarchical.xlsx"
    }

    # Auto-detect tissue type if not provided
    if (is.null(known_tissue_type)) {
        print("Guessing tissue type...")
        tissue_result = auto_detect_tissue_type(path_to_db_file = custom_marker_file,
                                                seuratObject = seurat_object,
                                                scaled = scaled,
                                                assay = assay)
        known_tissue_type = tissue_result$tissue[1]
        print(paste("Detected tissue type:", known_tissue_type))
    }

    # Read the hierarchical database
    db_all = openxlsx::read.xlsx(custom_marker_file)
    db_tissue = db_all[db_all$tissueType == known_tissue_type,]

    if(nrow(db_tissue) == 0) {
        stop(paste("No markers found for tissue type:", known_tissue_type))
    }

    # Get unique broad categories
    broad_categories = unique(db_tissue$broadCategory)

    # Extract data based on Seurat version
    data_type <- if (scaled) "scale.data" else "counts"
    package_type <- data_type %in% names(attributes(seurat_object[[assay]]))

    if(package_type){
        print("Using Seurat v4 object")
        scRNAseqData <- as.matrix(slot(seurat_object[[assay]], data_type))
    } else {
        print("Using Seurat v5 object")
        if (data_type == "scale.data") {
            scRNAseqData <- as.matrix(seurat_object[[assay]]$scale.data)
        } else {
            scRNAseqData <- as.matrix(seurat_object[[assay]]$counts)
        }
    }

    #========================================================================#
    # Step 1: Annotate at BROAD level
    #========================================================================#

    print("Step 1/2: Annotating broad cell categories...")

    # Create broad category marker lists
    broad_markers_pos = list()
    broad_markers_neg = list()

    for(broad_cat in broad_categories) {
        # Get all cell types in this broad category
        cells_in_category = db_tissue[db_tissue$broadCategory == broad_cat,]

        # Combine all positive markers for this broad category
        pos_markers = unique(unlist(strsplit(paste(cells_in_category$geneSymbolmore1, collapse=","), ",")))
        pos_markers = pos_markers[pos_markers != "" & !is.na(pos_markers)]

        # Combine all negative markers for this broad category
        neg_markers = unique(unlist(strsplit(paste(cells_in_category$geneSymbolmore2, collapse=","), ",")))
        neg_markers = neg_markers[neg_markers != "" & !is.na(neg_markers)]

        broad_markers_pos[[broad_cat]] = pos_markers
        broad_markers_neg[[broad_cat]] = neg_markers
    }

    # Run ScType at broad level
    es_broad = sctype_score(scRNAseqData = scRNAseqData,
                           scaled = scaled,
                           gs = broad_markers_pos,
                           gs2 = broad_markers_neg)

    # Aggregate broad scores by cluster
    cL_broad = do.call("rbind", lapply(unique(seurat_object@meta.data$seurat_clusters), function(cl){
        es.max.cl = sort(rowSums(es_broad[, rownames(seurat_object@meta.data[seurat_object@meta.data$seurat_clusters==cl, ])]),
                        decreasing = TRUE)
        head(data.frame(cluster = cl,
                       type = names(es.max.cl),
                       scores = es.max.cl,
                       ncells = sum(seurat_object@meta.data$seurat_clusters==cl)), 10)
    }))

    sctype_broad = cL_broad %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

    # Set low-confident broad assignments to "Unknown"
    sctype_broad$type[as.numeric(as.character(sctype_broad$scores)) < sctype_broad$ncells/4] = "Unknown"

    #========================================================================#
    # Step 2: Annotate at FINE level within each broad category
    #========================================================================#

    print("Step 2/2: Annotating fine cell subtypes...")

    # Prepare gene sets for fine-level annotation
    gs_list_fine = gene_sets_prepare(custom_marker_file, known_tissue_type)

    # Run ScType at fine level
    es_fine = sctype_score(scRNAseqData = scRNAseqData,
                          scaled = scaled,
                          gs = gs_list_fine$gs_positive,
                          gs2 = gs_list_fine$gs_negative)

    # Aggregate fine scores by cluster
    cL_fine = do.call("rbind", lapply(unique(seurat_object@meta.data$seurat_clusters), function(cl){
        es.max.cl = sort(rowSums(es_fine[, rownames(seurat_object@meta.data[seurat_object@meta.data$seurat_clusters==cl, ])]),
                        decreasing = TRUE)
        head(data.frame(cluster = cl,
                       type = names(es.max.cl),
                       scores = es.max.cl,
                       ncells = sum(seurat_object@meta.data$seurat_clusters==cl)), 10)
    }))

    sctype_fine = cL_fine %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

    # Set low-confident fine assignments to match broad category or "Unknown"
    for(i in 1:nrow(sctype_fine)) {
        cl = sctype_fine$cluster[i]
        fine_score = as.numeric(as.character(sctype_fine$scores[i]))
        threshold = sctype_fine$ncells[i] / 4

        if(fine_score < threshold) {
            # Low confidence at fine level - use broad category instead
            broad_assignment = sctype_broad$type[sctype_broad$cluster == cl]
            sctype_fine$type[i] = as.character(broad_assignment)
        }
    }

    #========================================================================#
    # Step 3: Add annotations to Seurat object
    #========================================================================#

    seurat_object_res = seurat_object
    seurat_object_res@meta.data[[broad_name]] = ""
    seurat_object_res@meta.data[[fine_name]] = ""

    # Add broad annotations
    for(j in unique(sctype_broad$cluster)){
        cl_type = sctype_broad[sctype_broad$cluster==j,]
        seurat_object_res@meta.data[seurat_object_res@meta.data$seurat_clusters == j, broad_name] = as.character(cl_type$type[1])
    }

    # Add fine annotations
    for(j in unique(sctype_fine$cluster)){
        cl_type = sctype_fine[sctype_fine$cluster==j,]
        seurat_object_res@meta.data[seurat_object_res@meta.data$seurat_clusters == j, fine_name] = as.character(cl_type$type[1])
    }

    #========================================================================#
    # Step 4: Plot results (if requested)
    #========================================================================#

    if(plot){
        if (!requireNamespace("patchwork", quietly = TRUE)) {
            stop("Package 'patchwork' is required for plotting. Install with: install.packages('patchwork')")
        }
        p1 = DimPlot(seurat_object_res, reduction = "umap", group.by = broad_name, label = TRUE, repel = TRUE) +
            ggtitle("Broad Cell Categories")
        p2 = DimPlot(seurat_object_res, reduction = "umap", group.by = fine_name, label = TRUE, repel = TRUE) +
            ggtitle("Fine Cell Subtypes")
        print(p1 / p2)
    }

    #========================================================================#
    # Print summary
    #========================================================================#

    print("\n=== Hierarchical Annotation Summary ===")
    print(paste("Broad categories added to:", broad_name))
    print(paste("Fine subtypes added to:", fine_name))

    print("\nBroad category distribution:")
    print(table(seurat_object_res@meta.data[[broad_name]]))

    print("\nFine subtype distribution:")
    print(table(seurat_object_res@meta.data[[fine_name]]))

    return(seurat_object_res)
}


#' @title Get hierarchical annotation for single cluster
#' @name get_cluster_hierarchy
#' @description Get both broad and fine annotations for a specific cluster
#'
#' @param seurat_object Seurat object with hierarchical annotations
#' @param cluster_id Cluster ID to query
#' @param broad_name Metadata column name for broad categories
#' @param fine_name Metadata column name for fine subtypes
#'
#' @return List with broad and fine annotations
#' @export

get_cluster_hierarchy <- function(seurat_object,
                                  cluster_id,
                                  broad_name = "sctype_broad",
                                  fine_name = "sctype_fine") {

    cluster_cells = seurat_object@meta.data[seurat_object@meta.data$seurat_clusters == cluster_id,]

    if(nrow(cluster_cells) == 0) {
        stop(paste("Cluster", cluster_id, "not found"))
    }

    broad_type = unique(cluster_cells[[broad_name]])[1]
    fine_type = unique(cluster_cells[[fine_name]])[1]

    result = list(
        cluster = cluster_id,
        broad_category = broad_type,
        fine_subtype = fine_type,
        n_cells = nrow(cluster_cells),
        is_refined = broad_type != fine_type
    )

    return(result)
}


#' @title Print hierarchical annotation table
#' @name print_hierarchy_table
#' @description Print formatted table of hierarchical annotations
#'
#' @param seurat_object Seurat object with hierarchical annotations
#' @param broad_name Metadata column name for broad categories
#' @param fine_name Metadata column name for fine subtypes
#'
#' @export

print_hierarchy_table <- function(seurat_object,
                                  broad_name = "sctype_broad",
                                  fine_name = "sctype_fine") {

    hierarchy_table = seurat_object@meta.data %>%
        group_by(seurat_clusters, !!sym(broad_name), !!sym(fine_name)) %>%
        summarise(n_cells = n(), .groups = "drop") %>%
        arrange(seurat_clusters)

    colnames(hierarchy_table) = c("Cluster", "Broad Category", "Fine Subtype", "N Cells")

    print(hierarchy_table, n = Inf)

    return(hierarchy_table)
}
