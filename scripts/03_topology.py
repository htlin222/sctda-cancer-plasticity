#!/usr/bin/env python
"""
Step 4: Topological data analysis — persistent homology + Mapper.

Per-timepoint PH, permutation tests for H1 significance,
cross-timepoint comparison, and Mapper with multiple filter functions.

Run: python scripts/03_topology.py --seed 42
"""

import argparse
import logging
import sys

sys.path.insert(0, "src")

import pickle

import numpy as np
import pandas as pd
import scanpy as sc

from sctda_plasticity.config import (
    DATA_PROCESSED,
    FIGURES_DIR,
    MAPPER_N_CUBES,
    MAPPER_OVERLAP,
    N_PERMUTATIONS,
    PH_MAX_DIM,
    PH_N_PCS,
    RESULTS_DIR,
    SEED,
    TABLES_DIR,
)
from sctda_plasticity.topology import (
    compare_persistence_diagrams,
    compute_persistent_homology,
    permutation_test_h1,
    run_mapper,
)
from sctda_plasticity.visualize import (
    plot_permutation_test,
    plot_persistence_barcodes,
    plot_persistence_comparison,
    plot_wasserstein_heatmap,
    save_figure,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def main(seed: int = SEED):
    logger.info(f"=== Topological Analysis Pipeline (seed={seed}) ===")

    # ── Load data ───────────────────────────────────────────
    h5ad_path = DATA_PROCESSED / "pc9_erlotinib_analyzed.h5ad"
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Run standard analysis first: {h5ad_path}")

    adata = sc.read_h5ad(h5ad_path)
    logger.info(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # Create output directories
    ph_dir = RESULTS_DIR / "ph"
    mapper_dir = RESULTS_DIR / "mapper"
    ph_dir.mkdir(parents=True, exist_ok=True)
    mapper_dir.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # ── 4a. Per-timepoint persistent homology ───────────────
    logger.info("\n--- Per-timepoint Persistent Homology ---")

    if "timepoint" in adata.obs.columns:
        timepoints = sorted(adata.obs["timepoint"].unique())
    else:
        timepoints = ["all"]
    dgms_by_tp = {}
    ph_results = {}

    for tp in timepoints:
        if tp == "all":
            cells = adata
        else:
            cells = adata[adata.obs["timepoint"] == tp]

        logger.info(f"\nTimepoint {tp}: {cells.n_obs} cells")
        X_pca = cells.obsm["X_pca"]

        ph = compute_persistent_homology(X_pca, maxdim=PH_MAX_DIM, n_pcs=PH_N_PCS)
        dgms_by_tp[tp] = ph["dgms"]
        ph_results[tp] = {
            "n_cells": cells.n_obs,
            "max_persistence": ph["max_persistence"],
            "n_h0_features": len(ph["dgms"][0]),
            "n_h1_features": len(ph["dgms"][1]) if PH_MAX_DIM >= 1 else 0,
        }

        # Save diagrams
        for dim in range(PH_MAX_DIM + 1):
            np.save(ph_dir / f"dgm{dim}_{tp}.npy", ph["dgms"][dim])

        logger.info(f"  Max persistence: {ph['max_persistence']}")
        logger.info(f"  H1 features: {ph_results[tp]['n_h1_features']}")

    # Save summary
    ph_summary = pd.DataFrame(ph_results).T
    ph_summary.to_csv(TABLES_DIR / "ph_summary_by_timepoint.csv")

    # ── Plot persistence barcodes ───────────────────────────
    fig = plot_persistence_barcodes(
        dgms_by_tp, max_dim=PH_MAX_DIM, title="Persistence Barcodes by Timepoint"
    )
    save_figure(fig, "fig2_a_persistence_barcodes")

    # ── Plot H1 persistence diagrams overlay ────────────────
    h1_dgms = {}
    for tp, dgms in dgms_by_tp.items():
        if PH_MAX_DIM >= 1 and len(dgms[1]) > 0:
            h1_dgms[tp] = dgms[1]
    if h1_dgms:
        fig = plot_persistence_comparison(h1_dgms, dim=1)
        save_figure(fig, "fig2_b_h1_persistence_diagrams")

    # ── 4b. Permutation test for H1 ────────────────────────
    logger.info("\n--- Permutation Test for H₁ ---")

    perm_results = {}
    for tp in timepoints:
        if tp == "all":
            cells = adata
        else:
            cells = adata[adata.obs["timepoint"] == tp]

        max_h1 = ph_results[tp]["max_persistence"].get(1, 0.0)
        if max_h1 == 0:
            logger.info(f"  {tp}: no H₁ features, skipping permutation test")
            continue

        logger.info(f"\n  Permutation test for {tp} (observed max H₁ = {max_h1:.4f})")

        # Get scaled expression for permutation (need pre-PCA data)
        X_scaled = cells.X
        if hasattr(X_scaled, "toarray"):
            X_scaled = X_scaled.toarray()
        X_scaled = np.array(X_scaled, dtype=np.float64)

        perm = permutation_test_h1(
            X_scaled,
            observed_max_h1=max_h1,
            n_permutations=N_PERMUTATIONS,
            n_pcs=PH_N_PCS,
            seed=seed,
        )
        perm_results[tp] = perm

        # Plot
        fig = plot_permutation_test(
            perm["null_distribution"], perm["observed"], perm["p_value"]
        )
        save_figure(fig, f"fig2_c_permtest_{tp}")

        # Save null distribution
        np.save(ph_dir / f"null_dist_h1_{tp}.npy", perm["null_distribution"])

    # Save permutation results
    if perm_results:
        perm_summary = pd.DataFrame({
            tp: {"observed_max_h1": r["observed"], "p_value": r["p_value"]}
            for tp, r in perm_results.items()
        }).T
        perm_summary.to_csv(TABLES_DIR / "permutation_test_h1.csv")

    # ── 4c. Cross-timepoint comparison ──────────────────────
    logger.info("\n--- Cross-Timepoint Comparison ---")

    if len(timepoints) > 1 and all(tp in h1_dgms for tp in timepoints if tp != "all"):
        tp_list = [tp for tp in timepoints if tp in h1_dgms]
        n = len(tp_list)
        wass_dist = np.zeros((n, n))

        for i in range(n):
            for j in range(i + 1, n):
                d = compare_persistence_diagrams(
                    h1_dgms[tp_list[i]], h1_dgms[tp_list[j]], metric="wasserstein"
                )
                wass_dist[i, j] = d
                wass_dist[j, i] = d

        fig = plot_wasserstein_heatmap(wass_dist, tp_list)
        save_figure(fig, "fig2_d_wasserstein_heatmap")

        np.save(ph_dir / "wasserstein_distances_h1.npy", wass_dist)
        pd.DataFrame(wass_dist, index=tp_list, columns=tp_list).to_csv(
            TABLES_DIR / "wasserstein_distances_h1.csv"
        )

    # ── 4d. Mapper with multiple filter functions ───────────
    logger.info("\n--- Mapper Graphs ---")

    X_pca = adata.obsm["X_pca"][:, :PH_N_PCS]

    filters = {"pc1": adata.obsm["X_pca"][:, 0]}

    if "emt_score" in adata.obs.columns:
        filters["emt_score"] = adata.obs["emt_score"].values

    if "X_diffmap" in adata.obsm:
        filters["diffusion_1"] = adata.obsm["X_diffmap"][:, 1]  # skip DC0 (trivial)

    mapper_results = {}
    for filter_name, filter_values in filters.items():
        logger.info(f"\n  Mapper with filter={filter_name}")

        result = run_mapper(
            X_pca,
            filter_values,
            n_cubes=MAPPER_N_CUBES,
            overlap=MAPPER_OVERLAP,
        )

        n_nodes = len(result["graph"]["nodes"])
        n_edges = sum(len(v) for v in result["graph"]["links"].values())
        logger.info(f"  {filter_name}: {n_nodes} nodes, {n_edges} edges")

        if n_nodes == 0:
            logger.warning(f"  {filter_name}: 0 nodes, skipping save")
            continue

        # Save Mapper graph
        with open(mapper_dir / f"mapper_{filter_name}.pkl", "wb") as f:
            pickle.dump(result["graph"], f)

        # Save HTML visualization
        result["mapper"].visualize(
            result["graph"],
            path_html=str(mapper_dir / f"mapper_{filter_name}.html"),
            title=f"Mapper ({filter_name})",
        )

        mapper_results[filter_name] = result

    # ── Mapper parameter sensitivity ────────────────────────
    logger.info("\n--- Mapper Parameter Sensitivity ---")
    sensitivity = []
    best_filter = "emt_score" if "emt_score" in filters else "pc1"

    for n_cubes in [10, 15, 20, 25, 30]:
        for overlap in [0.2, 0.3, 0.4, 0.5]:
            result = run_mapper(
                X_pca, filters[best_filter], n_cubes=n_cubes, overlap=overlap
            )
            n_nodes = len(result["graph"]["nodes"])
            n_edges = sum(len(v) for v in result["graph"]["links"].values())
            sensitivity.append({
                "n_cubes": n_cubes,
                "overlap": overlap,
                "n_nodes": n_nodes,
                "n_edges": n_edges,
            })

    pd.DataFrame(sensitivity).to_csv(
        TABLES_DIR / "mapper_sensitivity.csv", index=False
    )

    logger.info("\n=== Topological Analysis Complete ===")
    logger.info(f"Results saved to {RESULTS_DIR}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Topological data analysis")
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    main(seed=args.seed)
