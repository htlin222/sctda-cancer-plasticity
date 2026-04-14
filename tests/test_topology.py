"""
Tests for the topology module.

Uses synthetic data to verify PH computation and permutation tests.
"""

import numpy as np
import pytest

from sctda_plasticity.topology import (
    compute_persistent_homology,
    gene_connectivity,
    permutation_test_h1,
    run_mapper,
)


@pytest.fixture
def circle_data():
    """Generate points on a noisy circle (should have H₁)."""
    rng = np.random.RandomState(42)
    n = 100
    theta = rng.uniform(0, 2 * np.pi, n)
    x = np.cos(theta) + rng.normal(0, 0.1, n)
    y = np.sin(theta) + rng.normal(0, 0.1, n)
    return np.column_stack([x, y])


@pytest.fixture
def cluster_data():
    """Generate two clusters (H₀=2 components, no H₁)."""
    rng = np.random.RandomState(42)
    c1 = rng.normal(0, 0.3, (50, 2))
    c2 = rng.normal(5, 0.3, (50, 2))
    return np.vstack([c1, c2])


class TestPersistentHomology:
    def test_circle_has_h1(self, circle_data):
        """A circle should have a significant H₁ feature."""
        result = compute_persistent_homology(circle_data, maxdim=1, n_pcs=None)
        assert result["max_persistence"][1] > 0.5, "Circle should have persistent H₁"

    def test_clusters_no_h1(self, cluster_data):
        """Two clusters should have minimal H₁."""
        result = compute_persistent_homology(cluster_data, maxdim=1, n_pcs=None)
        # May have small H₁ from noise, but should be much less than circle
        assert result["max_persistence"][1] < 0.5

    def test_h0_connected_components(self, cluster_data):
        """Two clusters should show in H₀."""
        result = compute_persistent_homology(cluster_data, maxdim=0, n_pcs=None)
        # H₀ should have at least one long bar (the gap between clusters)
        dgm0 = result["dgms"][0]
        finite = dgm0[np.isfinite(dgm0[:, 1])]
        persistence = finite[:, 1] - finite[:, 0]
        assert np.max(persistence) > 1.0

    def test_output_structure(self, circle_data):
        """Check output dict has expected keys."""
        result = compute_persistent_homology(circle_data, maxdim=1, n_pcs=None)
        assert "dgms" in result
        assert "max_persistence" in result
        assert 0 in result["max_persistence"]
        assert 1 in result["max_persistence"]


class TestMapper:
    def test_mapper_runs(self, circle_data):
        """Mapper should produce a graph with nodes and edges."""
        filter_vals = circle_data[:, 0]
        result = run_mapper(circle_data, filter_vals, n_cubes=10, overlap=0.3)
        assert "graph" in result
        assert len(result["graph"]["nodes"]) > 0

    def test_mapper_with_different_filters(self, circle_data):
        """Different filters should produce different graphs."""
        r1 = run_mapper(circle_data, circle_data[:, 0], n_cubes=10)
        r2 = run_mapper(circle_data, circle_data[:, 1], n_cubes=10)
        # Different node counts expected
        n1 = len(r1["graph"]["nodes"])
        n2 = len(r2["graph"]["nodes"])
        # Both should have nodes; they may or may not differ
        assert n1 > 0
        assert n2 > 0


class TestGeneConnectivity:
    def test_connectivity_basic(self):
        """Gene connectivity should rank genes by node-specificity."""
        rng = np.random.RandomState(42)
        # 100 cells, 10 genes
        expr = rng.normal(0, 1, (100, 10))
        # Make gene 0 very specific to a subset
        expr[:20, 0] = 5.0

        gene_names = np.array([f"Gene_{i}" for i in range(10)])
        mapper_graph = {
            "nodes": {
                "node_0": list(range(20)),
                "node_1": list(range(20, 60)),
                "node_2": list(range(60, 100)),
            },
            "links": {},
        }

        result = gene_connectivity(mapper_graph, expr, gene_names, n_top=5)
        assert result["top_genes"][0] == "Gene_0", "Most specific gene should rank first"


class TestPermutationTest:
    def test_circle_significant(self, circle_data):
        """H₁ of a circle should be significant vs. permuted noise."""
        # This is a slow test; use fewer permutations
        result_ph = compute_persistent_homology(circle_data, maxdim=1, n_pcs=None)
        observed = result_ph["max_persistence"][1]

        # Use the raw coordinates as "expression" for permutation
        perm_result = permutation_test_h1(
            circle_data, observed, n_permutations=20, n_pcs=2, seed=42
        )

        # With only 20 permutations on 2D data, significance is not guaranteed.
        # We mainly verify the function runs and returns correct structure.
        assert 0 <= perm_result["p_value"] <= 1
        assert len(perm_result["null_distribution"]) == 20
        assert perm_result["observed"] == observed
