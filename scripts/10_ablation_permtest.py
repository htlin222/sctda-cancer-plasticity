#!/usr/bin/env python
"""
Cell cycle ablation permutation tests.

For each timepoint, after removing CC genes and re-running PCA,
run 100 permutations to test H1 significance.

Run: python scripts/10_ablation_permtest.py --seed 42
"""

import argparse
import logging
import sys

sys.path.insert(0, "src")

import numpy as np
import pandas as pd
import scanpy as sc

from sctda_plasticity.config import (
    DATA_PROCESSED,
    PH_MAX_DIM,
    PH_N_PCS,
    RESULTS_DIR,
    SEED,
    TABLES_DIR,
)
from sctda_plasticity.data import remove_cell_cycle_genes
from sctda_plasticity.topology import (
    compute_persistent_homology,
    permutation_test_h1,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def main(seed: int = SEED, n_perm: int = 100, max_cells: int = 1000):
    logger.info(f"=== CC Ablation Permutation Tests (seed={seed}) ===")

    adata = sc.read_h5ad(DATA_PROCESSED / "pc9_erlotinib_analyzed.h5ad")
    logger.info(f"Loaded: {adata.n_obs} cells")

    # Remove CC genes and re-do PCA
    adata_nocc = remove_cell_cycle_genes(adata)
    sc.pp.scale(adata_nocc, max_value=10)
    sc.tl.pca(
        adata_nocc, n_comps=min(50, adata_nocc.n_vars - 1), random_state=seed,
    )

    ablation_dir = RESULTS_DIR / "ablation"
    ablation_dir.mkdir(parents=True, exist_ok=True)

    timepoints = sorted(adata.obs["timepoint"].unique())
    rng = np.random.RandomState(seed)

    results = []
    for tp in timepoints:
        cells = adata_nocc[adata_nocc.obs["timepoint"] == tp]
        if cells.n_obs < 30:
            continue

        if cells.n_obs > max_cells:
            idx = rng.choice(cells.n_obs, max_cells, replace=False)
            cells = cells[idx]

        # Observed
        ph = compute_persistent_homology(
            cells.obsm["X_pca"], maxdim=PH_MAX_DIM, n_pcs=PH_N_PCS,
        )
        obs_h1 = ph["max_persistence"].get(1, 0)
        if obs_h1 == 0:
            continue

        logger.info(f"\n--- {tp} (n={cells.n_obs}, obs H1={obs_h1:.3f}) ---")

        X = cells.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = np.array(X, dtype=np.float64)

        perm = permutation_test_h1(
            X, obs_h1, n_permutations=n_perm, n_pcs=PH_N_PCS, seed=seed,
        )

        logger.info(f"  p = {perm['p_value']:.4f}")
        results.append({
            "timepoint": tp,
            "n_cells": cells.n_obs,
            "observed_h1_nocc": obs_h1,
            "p_value_nocc": perm["p_value"],
            "null_mean": float(np.mean(perm["null_distribution"])),
            "null_95th": float(np.percentile(perm["null_distribution"], 95)),
        })
        np.save(ablation_dir / f"null_dist_h1_{tp}_nocc.npy", perm["null_distribution"])

    df = pd.DataFrame(results)
    df.to_csv(TABLES_DIR / "permutation_test_h1_ablation.csv", index=False)

    logger.info("\n=== FINAL: CC ablation permutation tests ===")
    logger.info(df.to_string(index=False))

    sig = df[df["p_value_nocc"] < 0.05]
    logger.info(f"\nSignificant after CC ablation (p<0.05): {len(sig)}/{len(df)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--n-perm", type=int, default=100)
    parser.add_argument("--max-cells", type=int, default=1000)
    args = parser.parse_args()
    main(seed=args.seed, n_perm=args.n_perm, max_cells=args.max_cells)
