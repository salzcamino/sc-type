"""
ScType Python wrapper for scanpy/AnnData integration

This module provides Python/scanpy-compatible interfaces to ScType R functions
using rpy2 for seamless integration with Python single-cell analysis workflows.

Author: Claude (AI assistant)
Date: 2025-11-15
License: GNU General Public License v3.0
"""

import warnings
import numpy as np
import pandas as pd
from typing import Optional, List, Union, Dict, Tuple
import os

try:
    import scanpy as sc
    from anndata import AnnData
except ImportError:
    raise ImportError("scanpy and anndata are required. Install with: pip install scanpy anndata")

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib import gridspec
except ImportError:
    raise ImportError("matplotlib is required. Install with: pip install matplotlib")

try:
    import seaborn as sns
except ImportError:
    warnings.warn("seaborn not found. Some visualizations may be limited. Install with: pip install seaborn")
    sns = None

try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri, numpy2ri
    from rpy2.robjects.packages import importr
    from rpy2.robjects.conversion import localconverter
except ImportError:
    raise ImportError("rpy2 is required for R integration. Install with: pip install rpy2")


class ScType:
    """
    ScType cell type annotation for Python/scanpy workflows.

    This class provides a Python interface to ScType R functions, allowing
    seamless integration with scanpy AnnData objects.

    Examples
    --------
    >>> import scanpy as sc
    >>> from sctype_python import ScType
    >>>
    >>> # Load your data
    >>> adata = sc.read_h5ad("pbmc3k.h5ad")
    >>>
    >>> # Initialize ScType
    >>> sctype = ScType()
    >>>
    >>> # Run annotation
    >>> adata = sctype.annotate(adata, tissue_type="Immune system")
    >>>
    >>> # View results
    >>> sc.pl.umap(adata, color='sctype_classification')
    """

    def __init__(self, github_repo: str = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master"):
        """
        Initialize ScType Python wrapper.

        Parameters
        ----------
        github_repo : str
            Base URL for ScType R scripts on GitHub
        """
        self.github_repo = github_repo
        self._setup_r_environment()

    def _setup_r_environment(self):
        """Setup R environment and activate converters."""
        # Activate converters
        pandas2ri.activate()
        numpy2ri.activate()

        # Import R packages
        try:
            self.r_base = importr('base')
            self.r_utils = importr('utils')
        except Exception as e:
            raise RuntimeError(f"Failed to import R base packages: {e}")

    def _source_r_script(self, script_name: str):
        """Source an R script from GitHub."""
        url = f"{self.github_repo}/R/{script_name}"
        ro.r(f'source("{url}")')

    def _adata_to_matrix(self, adata: AnnData, layer: Optional[str] = None,
                        use_raw: bool = False) -> Tuple[np.ndarray, List[str], List[str]]:
        """
        Convert AnnData to matrix format for R.

        Parameters
        ----------
        adata : AnnData
            AnnData object
        layer : str, optional
            Layer to use (e.g., 'scaled'). If None, uses .X
        use_raw : bool
            Whether to use .raw.X

        Returns
        -------
        matrix : np.ndarray
            Expression matrix (genes × cells)
        genes : list
            Gene names
        cells : list
            Cell names
        """
        if use_raw and adata.raw is not None:
            mat = adata.raw.X
            genes = adata.raw.var_names.tolist()
        elif layer is not None:
            mat = adata.layers[layer]
            genes = adata.var_names.tolist()
        else:
            mat = adata.X
            genes = adata.var_names.tolist()

        # Convert to dense if sparse
        if hasattr(mat, 'toarray'):
            mat = mat.toarray()

        # Transpose to genes × cells (R format)
        mat = mat.T

        cells = adata.obs_names.tolist()

        return mat, genes, cells

    def annotate(self,
                adata: AnnData,
                tissue_type: str,
                database_file: Optional[str] = None,
                layer: Optional[str] = None,
                scaled: bool = True,
                cluster_key: str = 'leiden',
                annotation_key: str = 'sctype_classification',
                plot: bool = False) -> AnnData:
        """
        Annotate cell types using ScType.

        Parameters
        ----------
        adata : AnnData
            Annotated data object with clustering
        tissue_type : str
            Tissue type (e.g., "Immune system", "Brain", "Liver")
        database_file : str, optional
            Path to custom marker database. If None, uses default ScTypeDB_full.xlsx
        layer : str, optional
            Layer to use for expression data. If None, uses .X
        scaled : bool
            Whether the data is scaled (default: True)
        cluster_key : str
            Key in adata.obs with cluster assignments (default: 'leiden')
        annotation_key : str
            Key to store annotations in adata.obs (default: 'sctype_classification')
        plot : bool
            Whether to generate UMAP plot (default: False)

        Returns
        -------
        adata : AnnData
            Modified AnnData with annotations in .obs[annotation_key]
        """
        # Check if cluster key exists
        if cluster_key not in adata.obs.columns:
            raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs. "
                           f"Run clustering first (e.g., sc.tl.leiden(adata))")

        # Source R scripts
        print("Loading ScType R functions...")
        self._source_r_script("gene_sets_prepare.R")
        self._source_r_script("sctype_score_.R")

        # Set database file
        if database_file is None:
            database_file = f"{self.github_repo}/ScTypeDB_full.xlsx"

        # Prepare gene sets
        print(f"Preparing gene sets for tissue: {tissue_type}")
        ro.r(f'''
        gs_list <- gene_sets_prepare("{database_file}", "{tissue_type}")
        ''')

        # Convert AnnData to matrix
        print("Converting data to R format...")
        mat, genes, cells = self._adata_to_matrix(adata, layer=layer, use_raw=False)

        # Create R matrix
        with localconverter(ro.default_converter + numpy2ri.converter):
            r_mat = ro.r.matrix(mat, nrow=len(genes), ncol=len(cells))
            ro.r.assign("scRNAseqData", r_mat)
            ro.r.assign("gene_names", ro.StrVector(genes))
            ro.r.assign("cell_names", ro.StrVector(cells))
            ro.r('rownames(scRNAseqData) <- gene_names')
            ro.r('colnames(scRNAseqData) <- cell_names')

        # Run ScType scoring
        print("Running ScType scoring...")
        ro.r(f'''
        es.max <- sctype_score(scRNAseqData = scRNAseqData,
                              scaled = {str(scaled).upper()},
                              gs = gs_list$gs_positive,
                              gs2 = gs_list$gs_negative)
        ''')

        # Get clusters
        clusters = adata.obs[cluster_key].astype(str).values
        unique_clusters = np.unique(clusters)

        # Aggregate scores by cluster
        print("Aggregating scores by cluster...")
        cluster_assignments = {}

        for cl in unique_clusters:
            # Get cells in cluster
            cells_in_cluster = np.where(clusters == cl)[0]
            cell_indices = ro.IntVector([i + 1 for i in cells_in_cluster])  # R is 1-indexed

            # Aggregate scores
            ro.r.assign("cluster_cells", cell_indices)
            ro.r('es_cluster <- rowSums(es.max[, cluster_cells, drop=FALSE])')
            ro.r('es_cluster_sorted <- sort(es_cluster, decreasing=TRUE)')

            # Get top cell type
            top_types = list(ro.r('names(es_cluster_sorted)'))
            top_scores = list(ro.r('as.numeric(es_cluster_sorted)'))

            if len(top_types) > 0:
                top_type = top_types[0]
                top_score = top_scores[0]

                # Apply confidence threshold (ncells/4)
                ncells = len(cells_in_cluster)
                if top_score < ncells / 4:
                    top_type = "Unknown"

                cluster_assignments[cl] = top_type
            else:
                cluster_assignments[cl] = "Unknown"

        # Assign to cells
        annotations = np.array([cluster_assignments.get(cl, "Unknown") for cl in clusters])
        adata.obs[annotation_key] = pd.Categorical(annotations)

        print(f"Added annotations to adata.obs['{annotation_key}']")
        print(f"\nCell type distribution:")
        print(adata.obs[annotation_key].value_counts())

        # Plot if requested
        if plot:
            if 'X_umap' in adata.obsm:
                sc.pl.umap(adata, color=annotation_key, title="ScType Annotations")
            else:
                warnings.warn("UMAP not found. Run sc.tl.umap(adata) first for visualization.")

        return adata

    def annotate_hierarchical(self,
                              adata: AnnData,
                              tissue_type: str,
                              database_file: Optional[str] = None,
                              layer: Optional[str] = None,
                              scaled: bool = True,
                              cluster_key: str = 'leiden',
                              broad_key: str = 'sctype_broad',
                              fine_key: str = 'sctype_fine',
                              plot: bool = False) -> AnnData:
        """
        Hierarchical cell type annotation (broad + fine levels).

        Provides two-level annotation: broad categories (e.g., "T cells")
        and fine subtypes (e.g., "CD8+ T cells").

        Parameters
        ----------
        adata : AnnData
            Annotated data object with clustering
        tissue_type : str
            Tissue type
        database_file : str, optional
            Path to hierarchical database. If None, uses ScTypeDB_hierarchical.xlsx
        layer : str, optional
            Layer to use
        scaled : bool
            Whether data is scaled
        cluster_key : str
            Cluster column in adata.obs
        broad_key : str
            Key for broad annotations (default: 'sctype_broad')
        fine_key : str
            Key for fine annotations (default: 'sctype_fine')
        plot : bool
            Generate plots

        Returns
        -------
        adata : AnnData
            Modified AnnData with both broad and fine annotations
        """
        # Check if cluster key exists
        if cluster_key not in adata.obs.columns:
            raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")

        # Source R scripts
        print("Loading ScType R functions for hierarchical annotation...")
        self._source_r_script("gene_sets_prepare.R")
        self._source_r_script("sctype_score_.R")

        # Set database file
        if database_file is None:
            database_file = f"{self.github_repo}/ScTypeDB_hierarchical.xlsx"

        print(f"Reading hierarchical database for tissue: {tissue_type}")

        # Load hierarchical database
        ro.r(f'db_hier <- openxlsx::read.xlsx("{database_file}")')
        ro.r(f'db_tissue <- db_hier[db_hier$tissueType == "{tissue_type}", ]')

        # Check if database has broadCategory column
        has_broad = list(ro.r('("broadCategory" %in% colnames(db_tissue))'))[0]

        if not has_broad:
            warnings.warn("Database does not have 'broadCategory' column. Using standard annotation.")
            return self.annotate(adata, tissue_type, database_file, layer, scaled,
                               cluster_key, fine_key, plot)

        # Step 1: Broad annotation
        print("Step 1: Performing broad category annotation...")

        # Get unique broad categories
        ro.r('broad_categories <- unique(db_tissue$broadCategory)')
        ro.r('broad_categories <- broad_categories[!is.na(broad_categories)]')

        # Aggregate markers by broad category
        ro.r('''
        broad_markers <- lapply(broad_categories, function(bc) {
            subset_rows <- db_tissue[db_tissue$broadCategory == bc, ]
            pos_markers <- unique(unlist(strsplit(paste(subset_rows$geneSymbolmore1, collapse=","), ",")))
            neg_markers <- unique(unlist(strsplit(paste(subset_rows$geneSymbolmore2, collapse=","), ",")))
            pos_markers <- pos_markers[pos_markers != "" & !is.na(pos_markers)]
            neg_markers <- neg_markers[neg_markers != "" & !is.na(neg_markers)]
            list(positive = pos_markers, negative = neg_markers)
        })
        names(broad_markers) <- broad_categories
        ''')

        # Create broad gene sets
        ro.r('gs_broad_pos <- lapply(broad_markers, function(x) x$positive)')
        ro.r('gs_broad_neg <- lapply(broad_markers, function(x) x$negative)')

        # Convert data and run broad scoring
        mat, genes, cells = self._adata_to_matrix(adata, layer=layer)

        with localconverter(ro.default_converter + numpy2ri.converter):
            r_mat = ro.r.matrix(mat, nrow=len(genes), ncol=len(cells))
            ro.r.assign("scRNAseqData", r_mat)
            ro.r.assign("gene_names", ro.StrVector(genes))
            ro.r.assign("cell_names", ro.StrVector(cells))
            ro.r('rownames(scRNAseqData) <- gene_names')
            ro.r('colnames(scRNAseqData) <- cell_names')

        ro.r(f'''
        es_broad <- sctype_score(scRNAseqData = scRNAseqData,
                                 scaled = {str(scaled).upper()},
                                 gs = gs_broad_pos,
                                 gs2 = gs_broad_neg)
        ''')

        # Aggregate broad scores by cluster
        clusters = adata.obs[cluster_key].astype(str).values
        unique_clusters = np.unique(clusters)

        broad_assignments = {}
        for cl in unique_clusters:
            cells_in_cluster = np.where(clusters == cl)[0]
            cell_indices = ro.IntVector([i + 1 for i in cells_in_cluster])

            ro.r.assign("cluster_cells", cell_indices)
            ro.r('es_cluster_broad <- rowSums(es_broad[, cluster_cells, drop=FALSE])')
            ro.r('top_broad <- names(sort(es_cluster_broad, decreasing=TRUE))[1]')

            broad_assignments[cl] = str(list(ro.r('top_broad'))[0])

        # Step 2: Fine annotation
        print("Step 2: Performing fine-grained annotation...")

        # Prepare fine gene sets (all cell types)
        ro.r(f'gs_fine_list <- gene_sets_prepare("{database_file}", "{tissue_type}")')

        ro.r(f'''
        es_fine <- sctype_score(scRNAseqData = scRNAseqData,
                               scaled = {str(scaled).upper()},
                               gs = gs_fine_list$gs_positive,
                               gs2 = gs_fine_list$gs_negative)
        ''')

        # Aggregate fine scores by cluster and apply confidence threshold
        fine_assignments = {}
        for cl in unique_clusters:
            cells_in_cluster = np.where(clusters == cl)[0]
            cell_indices = ro.IntVector([i + 1 for i in cells_in_cluster])
            ncells = len(cells_in_cluster)

            ro.r.assign("cluster_cells", cell_indices)
            ro.r('es_cluster_fine <- rowSums(es_fine[, cluster_cells, drop=FALSE])')
            ro.r('es_sorted <- sort(es_cluster_fine, decreasing=TRUE)')

            top_type = str(list(ro.r('names(es_sorted)[1]'))[0])
            top_score = float(list(ro.r('as.numeric(es_sorted)[1]'))[0])

            # Apply confidence threshold
            if top_score < ncells / 4:
                # Low confidence, use broad category
                fine_assignments[cl] = broad_assignments[cl]
            else:
                fine_assignments[cl] = top_type

        # Assign to cells
        broad_annot = np.array([broad_assignments.get(cl, "Unknown") for cl in clusters])
        fine_annot = np.array([fine_assignments.get(cl, "Unknown") for cl in clusters])

        adata.obs[broad_key] = pd.Categorical(broad_annot)
        adata.obs[fine_key] = pd.Categorical(fine_annot)

        print(f"\nHierarchical annotation complete!")
        print(f"Broad categories added to adata.obs['{broad_key}']")
        print(f"Fine subtypes added to adata.obs['{fine_key}']")

        print(f"\nBroad category distribution:")
        print(adata.obs[broad_key].value_counts())
        print(f"\nFine subtype distribution:")
        print(adata.obs[fine_key].value_counts())

        # Plot if requested
        if plot:
            if 'X_umap' in adata.obsm:
                import matplotlib.pyplot as plt
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

                # Broad categories
                sc.pl.umap(adata, color=broad_key, ax=ax1, show=False,
                          title="Broad Cell Categories")

                # Fine subtypes
                sc.pl.umap(adata, color=fine_key, ax=ax2, show=False,
                          title="Fine Cell Subtypes")

                plt.tight_layout()
                plt.show()
            else:
                warnings.warn("UMAP not found. Run sc.tl.umap(adata) for visualization.")

        return adata

    def add_uncertainty(self,
                       adata: AnnData,
                       tissue_type: str,
                       database_file: Optional[str] = None,
                       layer: Optional[str] = None,
                       scaled: bool = True,
                       cluster_key: str = 'leiden',
                       top_n: int = 3,
                       prefix: str = 'sctype') -> AnnData:
        """
        Add uncertainty scores and top N candidates.

        Calculates confidence metrics for annotations including:
        - Top N cell type candidates per cluster
        - Raw ScType scores for each candidate
        - Normalized confidence scores
        - Confidence levels (High/Medium/Low)

        Parameters
        ----------
        adata : AnnData
            AnnData object with clustering
        tissue_type : str
            Tissue type
        database_file : str, optional
            Marker database path
        layer : str, optional
            Layer to use
        scaled : bool
            Whether data is scaled
        cluster_key : str
            Cluster column
        top_n : int
            Number of top candidates to report (default: 3)
        prefix : str
            Prefix for new columns (default: 'sctype')

        Returns
        -------
        adata : AnnData
            Modified AnnData with uncertainty metrics in .obs
        """
        # Check if cluster key exists
        if cluster_key not in adata.obs.columns:
            raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")

        # Source R scripts
        print("Loading ScType R functions...")
        self._source_r_script("gene_sets_prepare.R")
        self._source_r_script("sctype_score_.R")

        # Set database file
        if database_file is None:
            database_file = f"{self.github_repo}/ScTypeDB_full.xlsx"

        # Prepare gene sets
        print(f"Preparing gene sets for tissue: {tissue_type}")
        ro.r(f'gs_list <- gene_sets_prepare("{database_file}", "{tissue_type}")')

        # Convert data
        print("Converting data to R format...")
        mat, genes, cells = self._adata_to_matrix(adata, layer=layer)

        with localconverter(ro.default_converter + numpy2ri.converter):
            r_mat = ro.r.matrix(mat, nrow=len(genes), ncol=len(cells))
            ro.r.assign("scRNAseqData", r_mat)
            ro.r.assign("gene_names", ro.StrVector(genes))
            ro.r.assign("cell_names", ro.StrVector(cells))
            ro.r('rownames(scRNAseqData) <- gene_names')
            ro.r('colnames(scRNAseqData) <- cell_names')

        # Run ScType scoring
        print("Running ScType scoring...")
        ro.r(f'''
        es.max <- sctype_score(scRNAseqData = scRNAseqData,
                              scaled = {str(scaled).upper()},
                              gs = gs_list$gs_positive,
                              gs2 = gs_list$gs_negative)
        ''')

        # Get clusters
        clusters = adata.obs[cluster_key].astype(str).values
        unique_clusters = np.unique(clusters)

        # Calculate top N candidates and scores for each cluster
        print(f"Calculating top {top_n} candidates per cluster...")

        cluster_results = {}
        for cl in unique_clusters:
            cells_in_cluster = np.where(clusters == cl)[0]
            cell_indices = ro.IntVector([i + 1 for i in cells_in_cluster])
            ncells = len(cells_in_cluster)

            ro.r.assign("cluster_cells", cell_indices)
            ro.r('es_cluster <- rowSums(es.max[, cluster_cells, drop=FALSE])')
            ro.r(f'es_sorted <- sort(es_cluster, decreasing=TRUE)[1:{min(top_n, 10)}]')

            # Get top N types and scores
            top_types = list(ro.r('names(es_sorted)'))
            top_scores = list(ro.r('as.numeric(es_sorted)'))

            # Pad with "None" if fewer than top_n candidates
            while len(top_types) < top_n:
                top_types.append("None")
                top_scores.append(0.0)

            # Calculate confidence score (normalized)
            if top_scores[0] > 0:
                # Confidence based on:
                # 1. Absolute score strength
                # 2. Relative score (compared to threshold)
                score_strength = min(top_scores[0] / (ncells / 2), 1.0)  # Cap at 1.0

                # Difference between top 2 scores (higher diff = more confident)
                if len(top_scores) > 1 and top_scores[1] > 0:
                    score_diff = (top_scores[0] - top_scores[1]) / top_scores[0]
                else:
                    score_diff = 1.0

                confidence_score = (score_strength + score_diff) / 2.0
            else:
                confidence_score = 0.0

            # Assign confidence level
            if confidence_score >= 0.7:
                confidence_level = "High"
            elif confidence_score >= 0.4:
                confidence_level = "Medium"
            elif confidence_score >= 0.1:
                confidence_level = "Low"
            else:
                confidence_level = "Very Low"

            # Assign top cell type (or Unknown if low confidence)
            if top_scores[0] < ncells / 4:
                assigned_type = "Unknown"
            else:
                assigned_type = top_types[0]

            cluster_results[cl] = {
                'assigned_type': assigned_type,
                'top_types': top_types[:top_n],
                'top_scores': top_scores[:top_n],
                'confidence_score': confidence_score,
                'confidence_level': confidence_level,
                'ncells': ncells
            }

        # Add results to AnnData object
        print("Adding uncertainty metrics to AnnData...")

        # Initialize columns
        adata.obs[prefix] = ""
        adata.obs[f"{prefix}_confidence"] = 0.0
        adata.obs[f"{prefix}_confidence_level"] = "Unknown"

        for i in range(top_n):
            adata.obs[f"{prefix}_top{i+1}"] = ""
            adata.obs[f"{prefix}_score{i+1}"] = 0.0

        # Fill in values
        for cl in unique_clusters:
            cluster_mask = clusters == cl
            result = cluster_results[cl]

            adata.obs.loc[cluster_mask, prefix] = result['assigned_type']
            adata.obs.loc[cluster_mask, f"{prefix}_confidence"] = result['confidence_score']
            adata.obs.loc[cluster_mask, f"{prefix}_confidence_level"] = result['confidence_level']

            for i in range(top_n):
                adata.obs.loc[cluster_mask, f"{prefix}_top{i+1}"] = result['top_types'][i]
                adata.obs.loc[cluster_mask, f"{prefix}_score{i+1}"] = result['top_scores'][i]

        # Convert categorical
        adata.obs[prefix] = pd.Categorical(adata.obs[prefix])
        adata.obs[f"{prefix}_confidence_level"] = pd.Categorical(
            adata.obs[f"{prefix}_confidence_level"],
            categories=["High", "Medium", "Low", "Very Low", "Unknown"],
            ordered=True
        )

        print(f"\nUncertainty scoring complete!")
        print(f"Added annotations to adata.obs['{prefix}']")
        print(f"\nCell type distribution:")
        print(adata.obs[prefix].value_counts())
        print(f"\nConfidence level distribution:")
        print(adata.obs[f"{prefix}_confidence_level"].value_counts())

        print(f"\nNew columns added:")
        print(f"  - {prefix}: Cell type annotations")
        print(f"  - {prefix}_confidence: Confidence score (0-1)")
        print(f"  - {prefix}_confidence_level: Confidence level")
        for i in range(top_n):
            print(f"  - {prefix}_top{i+1}, {prefix}_score{i+1}: Top {i+1} candidate and score")

        return adata

    def _get_markers_from_database(self,
                                   database_file: str,
                                   tissue_type: str) -> pd.DataFrame:
        """
        Extract marker genes from ScType database.

        Parameters
        ----------
        database_file : str
            Path or URL to marker database
        tissue_type : str
            Tissue type to filter

        Returns
        -------
        markers_df : pd.DataFrame
            DataFrame with columns: cellName, geneSymbolmore1, geneSymbolmore2
        """
        try:
            # Try to load with openpyxl (for .xlsx files)
            import openpyxl
            if database_file.startswith('http'):
                import urllib.request
                tmp_file = '/tmp/sctype_db.xlsx'
                urllib.request.urlretrieve(database_file, tmp_file)
                database_file = tmp_file

            df = pd.read_excel(database_file, engine='openpyxl')
        except ImportError:
            # Fallback: use R to read the database
            ro.r(f'db <- readxl::read_excel("{database_file}")')
            with localconverter(ro.default_converter + pandas2ri.converter):
                df = ro.conversion.rpy2py(ro.r('db'))

        # Filter by tissue type
        df = df[df['tissueType'] == tissue_type].copy()

        return df

    def visualize_markers(self,
                         adata: AnnData,
                         tissue_type: str,
                         annotation_col: str = 'sctype_classification',
                         database_file: Optional[str] = None,
                         layer: Optional[str] = None,
                         top_n: int = 5,
                         plot_types: List[str] = ['violin', 'umap', 'dotplot', 'heatmap'],
                         save_plots: bool = False,
                         output_dir: str = 'sctype_plots',
                         figsize_violin: Tuple[int, int] = (12, 8),
                         figsize_umap: Tuple[int, int] = (14, 10),
                         figsize_dotplot: Tuple[int, int] = (12, 8),
                         figsize_heatmap: Tuple[int, int] = (10, 8)) -> Dict:
        """
        Visualize marker genes used for cell type annotation.

        Generates violin plots, UMAP feature plots, dotplots, and heatmaps
        showing the expression of positive and negative markers for each
        annotated cell type.

        Parameters
        ----------
        adata : AnnData
            AnnData object with ScType annotations
        tissue_type : str
            Tissue type used for annotation
        annotation_col : str
            Column in adata.obs with cell type annotations (default: 'sctype_classification')
        database_file : str, optional
            Path to marker database. If None, uses default
        layer : str, optional
            Layer to use for expression values. If None, uses .X
        top_n : int
            Number of top markers to show per cell type (default: 5)
        plot_types : list
            Types of plots to generate: 'violin', 'umap', 'dotplot', 'heatmap'
        save_plots : bool
            Whether to save plots to files (default: False)
        output_dir : str
            Directory for saved plots (default: 'sctype_plots')
        figsize_violin : tuple
            Figure size for violin plots (default: (12, 8))
        figsize_umap : tuple
            Figure size for UMAP plots (default: (14, 10))
        figsize_dotplot : tuple
            Figure size for dotplot (default: (12, 8))
        figsize_heatmap : tuple
            Figure size for heatmap (default: (10, 8))

        Returns
        -------
        plots : dict
            Dictionary with plot objects: {'violin': {...}, 'umap': {...}, 'dotplot': fig, 'heatmap': fig}

        Examples
        --------
        >>> # After running annotation
        >>> plots = sctype.visualize_markers(
        ...     adata,
        ...     tissue_type="Immune system",
        ...     annotation_col="sctype_classification",
        ...     top_n=5,
        ...     save_plots=True
        ... )
        >>>
        >>> # Access individual plots
        >>> plots['dotplot'].show()
        """
        # Check annotation column exists
        if annotation_col not in adata.obs.columns:
            raise ValueError(f"Annotation column '{annotation_col}' not found in adata.obs")

        # Set database file
        if database_file is None:
            database_file = f"{self.github_repo}/ScTypeDB_full.xlsx"

        # Get markers from database
        print("Loading marker database...")
        markers_df = self._get_markers_from_database(database_file, tissue_type)

        # Get unique cell types (excluding Unknown)
        cell_types = [ct for ct in adata.obs[annotation_col].unique() if ct != "Unknown"]

        if len(cell_types) == 0:
            raise ValueError("No annotated cell types found (all Unknown)")

        # Create output directory if saving
        if save_plots and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        plots = {}

        # Generate requested plot types
        if 'violin' in plot_types:
            print("Generating violin plots...")
            plots['violin'] = self._plot_violin(
                adata, markers_df, cell_types, annotation_col, layer,
                top_n, figsize_violin, save_plots, output_dir
            )

        if 'umap' in plot_types:
            if 'X_umap' not in adata.obsm:
                warnings.warn("UMAP coordinates not found. Skipping UMAP plots.")
            else:
                print("Generating UMAP feature plots...")
                plots['umap'] = self._plot_umap_features(
                    adata, markers_df, cell_types, annotation_col, layer,
                    top_n, figsize_umap, save_plots, output_dir
                )

        if 'dotplot' in plot_types:
            print("Generating dotplot...")
            plots['dotplot'] = self._plot_dotplot(
                adata, markers_df, cell_types, annotation_col, layer,
                top_n, figsize_dotplot, save_plots, output_dir
            )

        if 'heatmap' in plot_types:
            print("Generating heatmap...")
            plots['heatmap'] = self._plot_heatmap(
                adata, markers_df, cell_types, annotation_col, layer,
                top_n, figsize_heatmap, save_plots, output_dir
            )

        print("Visualization complete!")
        return plots

    def _plot_violin(self, adata, markers_df, cell_types, annotation_col, layer, top_n, figsize, save_plots, output_dir):
        """Generate violin plots for each cell type."""
        violin_plots = {}

        for ct in cell_types:
            # Get markers for this cell type
            ct_row = markers_df[markers_df['cellName'] == ct]

            if len(ct_row) == 0:
                continue

            pos_markers = ct_row.iloc[0]['geneSymbolmore1']
            neg_markers = ct_row.iloc[0]['geneSymbolmore2'] if 'geneSymbolmore2' in ct_row.columns else ""

            # Parse markers
            pos_markers = [m.strip() for m in str(pos_markers).split(',') if m.strip() and m.strip() != 'nan']
            neg_markers = [m.strip() for m in str(neg_markers).split(',') if m.strip() and m.strip() != 'nan']

            # Filter to genes present in dataset
            pos_markers = [g for g in pos_markers[:top_n] if g in adata.var_names]
            neg_markers = [g for g in neg_markers[:top_n] if g in adata.var_names]

            all_markers = pos_markers + neg_markers

            if len(all_markers) == 0:
                continue

            # Create violin plot
            fig, axes = plt.subplots(len(all_markers), 1, figsize=(figsize[0], figsize[1]))
            if len(all_markers) == 1:
                axes = [axes]

            for i, gene in enumerate(all_markers):
                # Get expression values
                if layer is not None:
                    expr = adata[:, gene].layers[layer].toarray().flatten() if hasattr(adata[:, gene].layers[layer], 'toarray') else adata[:, gene].layers[layer].flatten()
                else:
                    expr = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, 'toarray') else adata[:, gene].X.flatten()

                # Create dataframe for plotting
                plot_df = pd.DataFrame({
                    'expression': expr,
                    'cell_type': adata.obs[annotation_col].values
                })

                # Plot violin
                if sns is not None:
                    sns.violinplot(data=plot_df, x='cell_type', y='expression', ax=axes[i])
                else:
                    axes[i].violinplot([plot_df[plot_df['cell_type'] == ct_name]['expression'].values
                                       for ct_name in plot_df['cell_type'].unique()],
                                      positions=range(len(plot_df['cell_type'].unique())))

                marker_type = "(+)" if gene in pos_markers else "(-)"
                axes[i].set_title(f"{gene} {marker_type}")
                axes[i].set_ylabel("Expression")
                if i < len(all_markers) - 1:
                    axes[i].set_xticklabels([])
                    axes[i].set_xlabel("")
                else:
                    axes[i].tick_params(axis='x', rotation=45)

            plt.suptitle(f"Marker Expression: {ct}", fontsize=16, y=0.995)
            plt.tight_layout()

            if save_plots:
                filename = f"{output_dir}/violin_{ct.replace(' ', '_').replace('/', '_')}.png"
                plt.savefig(filename, dpi=300, bbox_inches='tight')
                print(f"  Saved: {filename}")

            violin_plots[ct] = fig

        return violin_plots

    def _plot_umap_features(self, adata, markers_df, cell_types, annotation_col, layer, top_n, figsize, save_plots, output_dir):
        """Generate UMAP feature plots for each cell type."""
        umap_plots = {}

        for ct in cell_types:
            # Get markers
            ct_row = markers_df[markers_df['cellName'] == ct]
            if len(ct_row) == 0:
                continue

            pos_markers = ct_row.iloc[0]['geneSymbolmore1']
            pos_markers = [m.strip() for m in str(pos_markers).split(',') if m.strip() and m.strip() != 'nan']
            pos_markers = [g for g in pos_markers[:top_n] if g in adata.var_names]

            if len(pos_markers) == 0:
                continue

            # Create grid of UMAP plots
            ncols = min(3, len(pos_markers))
            nrows = int(np.ceil(len(pos_markers) / ncols))

            fig = plt.figure(figsize=(figsize[0], figsize[1]))
            gs = gridspec.GridSpec(nrows, ncols, figure=fig)

            for i, gene in enumerate(pos_markers):
                ax = fig.add_subplot(gs[i])

                # Get expression
                if layer is not None:
                    expr = adata[:, gene].layers[layer].toarray().flatten() if hasattr(adata[:, gene].layers[layer], 'toarray') else adata[:, gene].layers[layer].flatten()
                else:
                    expr = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, 'toarray') else adata[:, gene].X.flatten()

                # Plot
                scatter = ax.scatter(
                    adata.obsm['X_umap'][:, 0],
                    adata.obsm['X_umap'][:, 1],
                    c=expr,
                    cmap='viridis',
                    s=5,
                    alpha=0.8
                )
                ax.set_title(gene, fontsize=12)
                ax.set_xlabel('UMAP1')
                ax.set_ylabel('UMAP2')
                plt.colorbar(scatter, ax=ax, label='Expression')

            plt.suptitle(f"Marker Expression (UMAP): {ct}", fontsize=16)
            plt.tight_layout()

            if save_plots:
                filename = f"{output_dir}/umap_{ct.replace(' ', '_').replace('/', '_')}.png"
                plt.savefig(filename, dpi=300, bbox_inches='tight')
                print(f"  Saved: {filename}")

            umap_plots[ct] = fig

        return umap_plots

    def _plot_dotplot(self, adata, markers_df, cell_types, annotation_col, layer, top_n, figsize, save_plots, output_dir):
        """Generate dotplot showing all markers across all cell types."""
        # Collect all markers
        all_markers = []
        for ct in cell_types:
            ct_row = markers_df[markers_df['cellName'] == ct]
            if len(ct_row) == 0:
                continue

            pos_markers = ct_row.iloc[0]['geneSymbolmore1']
            pos_markers = [m.strip() for m in str(pos_markers).split(',') if m.strip() and m.strip() != 'nan']
            pos_markers = [g for g in pos_markers[:top_n] if g in adata.var_names]
            all_markers.extend(pos_markers)

        all_markers = list(dict.fromkeys(all_markers))  # Remove duplicates, preserve order

        if len(all_markers) == 0:
            warnings.warn("No markers found in dataset")
            return None

        # Use scanpy's dotplot if available
        try:
            fig = sc.pl.dotplot(
                adata,
                var_names=all_markers,
                groupby=annotation_col,
                dendrogram=False,
                return_fig=True,
                figsize=figsize
            )

            if save_plots:
                filename = f"{output_dir}/dotplot_all_markers.png"
                fig.savefig(filename, dpi=300, bbox_inches='tight')
                print(f"  Saved: {filename}")

            return fig.fig if hasattr(fig, 'fig') else fig

        except Exception as e:
            warnings.warn(f"Could not generate dotplot: {e}")
            return None

    def _plot_heatmap(self, adata, markers_df, cell_types, annotation_col, layer, top_n, figsize, save_plots, output_dir):
        """Generate heatmap of average marker expression across cell types."""
        # Collect all markers
        marker_dict = {}
        for ct in cell_types:
            ct_row = markers_df[markers_df['cellName'] == ct]
            if len(ct_row) == 0:
                continue

            pos_markers = ct_row.iloc[0]['geneSymbolmore1']
            pos_markers = [m.strip() for m in str(pos_markers).split(',') if m.strip() and m.strip() != 'nan']
            pos_markers = [g for g in pos_markers[:top_n] if g in adata.var_names]

            for marker in pos_markers:
                if marker not in marker_dict:
                    marker_dict[marker] = []
                marker_dict[marker].append(ct)

        all_markers = list(marker_dict.keys())

        if len(all_markers) == 0:
            warnings.warn("No markers found in dataset")
            return None

        # Calculate mean expression per cell type
        mean_expr = pd.DataFrame(index=all_markers, columns=cell_types)

        for ct in cell_types:
            ct_cells = adata.obs[annotation_col] == ct
            for gene in all_markers:
                if layer is not None:
                    expr = adata[ct_cells, gene].layers[layer]
                else:
                    expr = adata[ct_cells, gene].X

                if hasattr(expr, 'toarray'):
                    expr = expr.toarray()

                mean_expr.loc[gene, ct] = np.mean(expr)

        mean_expr = mean_expr.astype(float)

        # Plot heatmap
        fig, ax = plt.subplots(figsize=figsize)

        if sns is not None:
            sns.heatmap(
                mean_expr,
                cmap='RdYlBu_r',
                center=0,
                cbar_kws={'label': 'Mean Expression'},
                ax=ax,
                yticklabels=True,
                xticklabels=True
            )
        else:
            im = ax.imshow(mean_expr.values, cmap='RdYlBu_r', aspect='auto')
            ax.set_xticks(range(len(cell_types)))
            ax.set_yticks(range(len(all_markers)))
            ax.set_xticklabels(cell_types, rotation=45, ha='right')
            ax.set_yticklabels(all_markers)
            plt.colorbar(im, ax=ax, label='Mean Expression')

        ax.set_title("Marker Expression Heatmap", fontsize=14)
        ax.set_xlabel("Cell Type", fontsize=12)
        ax.set_ylabel("Marker Gene", fontsize=12)

        plt.tight_layout()

        if save_plots:
            filename = f"{output_dir}/heatmap_all_markers.png"
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"  Saved: {filename}")

        return fig


def run_sctype(adata: AnnData,
              tissue_type: str,
              database_file: Optional[str] = None,
              layer: Optional[str] = None,
              scaled: bool = True,
              cluster_key: str = 'leiden',
              annotation_key: str = 'sctype_classification',
              plot: bool = False) -> AnnData:
    """
    Convenience function to run ScType annotation on scanpy AnnData.

    This is the main entry point for Python/scanpy users.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with clustering
    tissue_type : str
        Tissue type (e.g., "Immune system", "Brain", "Liver")
    database_file : str, optional
        Path to custom marker database
    layer : str, optional
        Layer to use. If None, uses .X
    scaled : bool
        Whether data is scaled (default: True)
    cluster_key : str
        Key in adata.obs with cluster assignments (default: 'leiden')
    annotation_key : str
        Key to store annotations (default: 'sctype_classification')
    plot : bool
        Generate UMAP plot (default: False)

    Returns
    -------
    adata : AnnData
        Modified AnnData with annotations in .obs[annotation_key]

    Examples
    --------
    >>> import scanpy as sc
    >>> from sctype_python import run_sctype
    >>>
    >>> # Standard scanpy workflow
    >>> adata = sc.read_h5ad("data.h5ad")
    >>> sc.pp.normalize_total(adata, target_sum=1e4)
    >>> sc.pp.log1p(adata)
    >>> sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    >>> sc.pp.scale(adata, max_value=10)
    >>> sc.tl.pca(adata, n_comps=50)
    >>> sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    >>> sc.tl.leiden(adata, resolution=0.8)
    >>> sc.tl.umap(adata)
    >>>
    >>> # Run ScType
    >>> adata = run_sctype(adata, tissue_type="Immune system", plot=True)
    >>>
    >>> # View results
    >>> print(adata.obs['sctype_classification'].value_counts())
    >>> sc.pl.umap(adata, color='sctype_classification')
    """
    sctype = ScType()
    return sctype.annotate(adata, tissue_type, database_file, layer, scaled,
                          cluster_key, annotation_key, plot)


# Convenience functions
def annotate_cells(adata: AnnData, tissue_type: str, **kwargs) -> AnnData:
    """Alias for run_sctype()."""
    return run_sctype(adata, tissue_type, **kwargs)


if __name__ == "__main__":
    print("ScType Python wrapper for scanpy")
    print("Usage:")
    print("  from sctype_python import run_sctype")
    print("  adata = run_sctype(adata, tissue_type='Immune system')")
