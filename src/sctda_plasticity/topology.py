"""
Topological data analysis: persistent homology, Mapper, and statistical tests.
"""

import logging
from typing import Optional

import kmapper as km
import numpy as np
from persim import bottleneck, wasserstein
from ripser import ripser
from sklearn.decomposition import PCA

logger = logging.getLogger(__name__)


def compute_persistent_homology(
    X: np.ndarray,
    maxdim: int = 1,
    n_pcs: Optional[int] = 30,
    thresh: float = np.inf,
) -> dict:
    """Compute persistent homology on a point cloud.

    Parameters
    ----------
    X : np.ndarray
        Point cloud (cells × features). If from AnnData, use adata.obsm['X_pca'].
    maxdim : int
        Maximum homology dimension (1 = detect loops, 2 = detect voids)
    n_pcs : int, optional
        If set, use only the first n_pcs columns of X
    thresh : float
        Maximum filtration value

    Returns
    -------
    dict with keys:
        'dgms': list of persistence diagrams per dimension
        'result': full ripser output
        'max_persistence': dict of max persistence per dimension
    """
    if n_pcs is not None and X.shape[1] > n_pcs:
        X = X[:, :n_pcs]

    result = ripser(X, maxdim=maxdim, thresh=thresh)

    max_pers = {}
    for dim in range(maxdim + 1):
        dgm = result["dgms"][dim]
        # Remove infinite bars for H₀
        finite = dgm[np.isfinite(dgm[:, 1])] if len(dgm) > 0 else dgm
        if len(finite) > 0:
            persistence = finite[:, 1] - finite[:, 0]
            max_pers[dim] = float(np.max(persistence))
        else:
            max_pers[dim] = 0.0

    logger.info(
        f"PH computed: {X.shape[0]} points, {X.shape[1]} dims. "
        f"Max persistence: {max_pers}"
    )

    return {
        "dgms": result["dgms"],
        "result": result,
        "max_persistence": max_pers,
    }


def permutation_test_h1(
    adata_X: np.ndarray,
    observed_max_h1: float,
    n_permutations: int = 500,
    n_pcs: int = 30,
    seed: int = 42,
) -> dict:
    """Permutation test for H₁ significance.

    Permutes gene labels independently per cell, re-runs PCA, computes PH.

    Parameters
    ----------
    adata_X : np.ndarray
        Scaled expression matrix (cells × genes) BEFORE PCA
    observed_max_h1 : float
        Observed maximum H₁ persistence from real data
    n_permutations : int
        Number of permutations
    n_pcs : int
        Number of PCs for PH computation
    seed : int
        Random seed

    Returns
    -------
    dict with 'p_value', 'null_distribution', 'observed'
    """
    rng = np.random.RandomState(seed)
    null_dist = np.zeros(n_permutations)

    for i in range(n_permutations):
        if (i + 1) % 50 == 0:
            logger.info(f"Permutation {i + 1}/{n_permutations}")

        # Permute gene labels per cell
        X_perm = adata_X.copy()
        for j in range(X_perm.shape[0]):
            rng.shuffle(X_perm[j, :])

        # PCA on permuted
        pca = PCA(n_components=n_pcs, random_state=seed)
        X_pca_perm = pca.fit_transform(X_perm)

        # PH
        result = ripser(X_pca_perm, maxdim=1)
        dgm1 = result["dgms"][1]

        if len(dgm1) > 0:
            finite = dgm1[np.isfinite(dgm1[:, 1])]
            if len(finite) > 0:
                null_dist[i] = np.max(finite[:, 1] - finite[:, 0])

    p_value = np.mean(null_dist >= observed_max_h1)

    logger.info(
        f"Permutation test: observed={observed_max_h1:.4f}, "
        f"p={p_value:.4f} ({n_permutations} perms)"
    )

    return {
        "p_value": p_value,
        "null_distribution": null_dist,
        "observed": observed_max_h1,
    }


def compare_persistence_diagrams(
    dgm1: np.ndarray,
    dgm2: np.ndarray,
    metric: str = "wasserstein",
) -> float:
    """Compute distance between two persistence diagrams.

    Parameters
    ----------
    dgm1, dgm2 : np.ndarray
        Persistence diagrams (n_features × 2)
    metric : str
        'wasserstein' or 'bottleneck'

    Returns
    -------
    float distance
    """
    if metric == "wasserstein":
        return wasserstein(dgm1, dgm2)
    elif metric == "bottleneck":
        return bottleneck(dgm1, dgm2)
    else:
        raise ValueError(f"Unknown metric: {metric}")


def run_mapper(
    X: np.ndarray,
    filter_values: np.ndarray,
    n_cubes: int = 15,
    overlap: float = 0.3,
    min_cluster_size: int = 3,
) -> dict:
    """Run Mapper algorithm.

    Parameters
    ----------
    X : np.ndarray
        Point cloud (cells × features)
    filter_values : np.ndarray
        Filter function values (1D array, one per cell)
    n_cubes : int
        Number of intervals for the cover
    overlap : float
        Overlap fraction between intervals
    min_cluster_size : int
        Minimum cluster size

    Returns
    -------
    dict with 'graph' (Mapper output) and 'mapper' (KeplerMapper instance)
    """
    from sklearn.cluster import DBSCAN

    mapper = km.KeplerMapper(verbose=0)

    # Project using provided filter
    lens = filter_values.reshape(-1, 1)

    # Estimate eps from pairwise distances (10th percentile of NN distances)
    from sklearn.neighbors import NearestNeighbors

    nn = NearestNeighbors(n_neighbors=min(10, X.shape[0] - 1))
    nn.fit(X)
    distances, _ = nn.kneighbors(X)
    eps = float(np.percentile(distances[:, -1], 50))

    graph = mapper.map(
        lens,
        X,
        cover=km.Cover(n_cubes=n_cubes, perc_overlap=overlap),
        clusterer=DBSCAN(eps=eps, min_samples=min_cluster_size),
    )

    n_nodes = len(graph["nodes"])
    n_edges = sum(len(v) for v in graph["links"].values())
    logger.info(f"Mapper: {n_nodes} nodes, {n_edges} edges")

    return {"graph": graph, "mapper": mapper}


def gene_connectivity(
    mapper_graph: dict,
    expression_matrix: np.ndarray,
    gene_names: np.ndarray,
    n_top: int = 100,
) -> dict:
    """Compute gene connectivity on a Mapper graph.

    A gene has high connectivity if cells that share similar global transcriptional
    profiles (same Mapper node) also share high expression of that gene.

    Parameters
    ----------
    mapper_graph : dict
        Mapper graph output
    expression_matrix : np.ndarray
        Expression matrix (cells × genes)
    gene_names : np.ndarray
        Gene names
    n_top : int
        Number of top genes to return

    Returns
    -------
    dict with 'scores' (DataFrame) and 'top_genes' (list)
    """
    import pandas as pd

    nodes = mapper_graph["nodes"]
    n_genes = expression_matrix.shape[1]
    connectivity = np.zeros(n_genes)

    for node_id, cell_indices in nodes.items():
        if len(cell_indices) < 2:
            continue
        node_expr = expression_matrix[cell_indices, :]
        node_mean = np.mean(node_expr, axis=0)
        global_mean = np.mean(expression_matrix, axis=0)
        connectivity += np.abs(node_mean - global_mean) * len(cell_indices)

    connectivity /= expression_matrix.shape[0]

    scores = pd.DataFrame({
        "gene": gene_names,
        "connectivity": connectivity,
    }).sort_values("connectivity", ascending=False)

    top_genes = scores.head(n_top)["gene"].tolist()

    logger.info(f"Gene connectivity: top gene = {top_genes[0]} ({connectivity.max():.4f})")

    return {"scores": scores, "top_genes": top_genes}
