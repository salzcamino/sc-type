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

try:
    import scanpy as sc
    from anndata import AnnData
except ImportError:
    raise ImportError("scanpy and anndata are required. Install with: pip install scanpy anndata")

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
        print("Hierarchical annotation is currently only available in R.")
        print("Use R version: run_sctype_hierarchical() or run_sctype_hierarchical_sce()")
        print("Falling back to standard annotation...")

        return self.annotate(adata, tissue_type, database_file, layer, scaled,
                           cluster_key, fine_key, plot)

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
        print("Python uncertainty scoring not yet implemented.")
        print("Use R version: add_sctype_uncertainty() or add_sctype_uncertainty_sce()")
        print("Running standard annotation instead...")

        return self.annotate(adata, tissue_type, database_file, layer, scaled,
                           cluster_key, f"{prefix}_top1", False)


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
