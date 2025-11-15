# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
# Modified by Claude Code, November 2025 - Added return_details and marker_weights parameters for ScType v2
#
# Functions on this page:
# sctype_score: calculate ScType scores and assign cell types
#
# @params: scRNAseqData - input scRNA-seq matrix (rownames - genes, column names - cells),
# @params: scale - indicates whether the matrix is scaled (TRUE by default)
# @params: gs - list of gene sets positively expressed in the cell type
# @params: gs2 - list of gene sets that should not be expressed in the cell type (NULL if not applicable)
# @params: return_details - return detailed scoring information for statistical testing (FALSE by default for backward compatibility)
# @params: marker_weights - custom marker weights data frame (NULL uses frequency-based weighting)

sctype_score <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, return_details = FALSE, marker_weights = NULL, ...){
  
  # check input matrix
  if(!is.matrix(scRNAseqData)){
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if(sum(dim(scRNAseqData))==0){
       warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }
  
  # marker sensitivity (custom weights or frequency-based)
  if (!is.null(marker_weights)) {
    # Use custom weights (e.g., TF-IDF)
    if (!is.data.frame(marker_weights)) {
      stop("marker_weights must be a data frame with 'gene_' and 'score_marker_sensitivity' columns")
    }
    if (!all(c("gene_", "score_marker_sensitivity") %in% colnames(marker_weights))) {
      stop("marker_weights must have 'gene_' and 'score_marker_sensitivity' columns")
    }
    marker_sensitivity = marker_weights
  } else {
    # Default: frequency-based weighting
    marker_stat = sort(table(unlist(gs)), decreasing = T);
    marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                        gene_ = names(marker_stat), stringsAsFactors = !1)
  }

  # convert gene names to Uppercase
  if(gene_names_to_uppercase){
    rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
  }
  
  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(seq_along(gs), function(d_){
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  if(!is.null(gs2)){
    gs2 = lapply(seq_along(gs2), function(d_){
      GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  }
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
  
  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
  
  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }
  
  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
  
  # combine scores
  es = do.call("rbind", lapply(names(gs), function(gss_){
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j];
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z)));

      # Handle negative markers
      if(!is.null(gs2) && length(gs2[[gss_]]) > 0){
        gz_2 = Z[gs2[[gss_]], j] * -1
        sum_t2 = sum(gz_2) / sqrt(length(gz_2));
        if(is.na(sum_t2)){
          sum_t2 = 0;
        }
      } else {
        sum_t2 = 0;
      }

      sum_t1 + sum_t2
    })
  })) 
  
  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows

  # Return detailed information if requested (for ScType v2 statistical testing)
  if(return_details){
    result <- list(
      scores = es.max,
      marker_sensitivity = cell_markers_genes_score,
      gene_sets_used = list(positive = gs, negative = gs2),
      scaled_matrix = Z,
      original_gs_names = list(positive = names_gs_cp, negative = names_gs_2_cp),
      was_scaled = scaled
    )
    return(result)
  } else {
    # Default behavior: return only scores (backward compatible)
    return(es.max)
  }
}
