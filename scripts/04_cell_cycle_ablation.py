#!/usr/bin/env python
"""
Step 5: Cell cycle ablation — the critical control experiment.

Remove cell cycle genes → re-run PCA → re-run PH.
Compare H₁ with vs. without cell cycle genes.

Run: python scripts/04_cell_cycle_ablation.py --seed 42
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
    FIGURES_DIR,
    N_PERMUTATIONS,
    PH_MAX_DIM,
    PH_N_PCS,
    RESULTS_DIR,
    SEED,
    TABLES_DIR,
)
from sctda_plasticity.data import G2M_GENES, S_GENES, remove_cell_cycle_genes
from sctda_plasticity.topology import compute_persistent_homology, permutation_test_h1
from sctda_plasticity.visualize import (
    plot_permutation_test,
    plot_persistence_barcodes,
    save_figure,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def main(seed: int = SEED):
    logger.info(f"=== Cell Cycle Ablation (seed={seed}) ===")

    h5ad_path = DATA_PROCESSED / "pc9_erlotinib_analyzed.h5ad"
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Run standard analysis first: {h5ad_path}")

    adata = sc.read_h5ad(h5ad_path)
    logger.info(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    ablation_dir = RESULTS_DIR / "ablation"
    ablation_dir.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # ── Remove cell cycle genes ─────────────────────────────
    cc_genes_in_data = [g for g in S_GENES + G2M_GENES if g in adata.var_names]
    logger.info(
        f"Cell cycle genes in HVG set: {len(cc_genes_in_data)} / {len(S_GENES) + len(G2M_GENES)}"
    )

    adata_nocc = remove_cell_cycle_genes(adata)

    # Re-run PCA on cell-cycle-free data
    logger.info("Re-running PCA without cell cycle genes...")
    sc.pp.scale(adata_nocc, max_value=10)
    sc.tl.pca(adata_nocc, n_comps=min(50, adata_nocc.n_vars - 1), random_state=seed)

    # ── Compare PH with and without CC genes ────────────────
    timepoints = (
        sorted(adata.obs["timepoint"].unique())
        if "timepoint" in adata.obs.columns
        else ["all"]
    )

    comparison = []
    dgms_with_cc = {}
    dgms_without_cc = {}

    for tp in timepoints:
        logger.info(f"\n--- Timepoint: {tp} ---")

        if tp == "all":
            cells_with = adata
            cells_without = adata_nocc
        else:
            cells_with = adata[adata.obs["timepoint"] == tp]
            cells_without = adata_nocc[adata_nocc.obs["timepoint"] == tp]

        ph_with = compute_persistent_homology(
            cells_with.obsm["X_pca"], maxdim=PH_MAX_DIM, n_pcs=PH_N_PCS
        )
        ph_without = compute_persistent_homology(
            cells_without.obsm["X_pca"], maxdim=PH_MAX_DIM, n_pcs=PH_N_PCS
        )

        dgms_with_cc[tp] = ph_with["dgms"]
        dgms_without_cc[tp] = ph_without["dgms"]

        max_h1_with = ph_with["max_persistence"].get(1, 0.0)
        max_h1_without = ph_without["max_persistence"].get(1, 0.0)

        logger.info(f"  With CC genes:    max H₁ = {max_h1_with:.4f}")
        logger.info(f"  Without CC genes: max H₁ = {max_h1_without:.4f}")

        ratio = max_h1_without / max_h1_with if max_h1_with > 0 else float("nan")
        logger.info(f"  Ratio (without/with): {ratio:.4f}")

        comparison.append({
            "timepoint": tp,
            "max_h1_with_cc": max_h1_with,
            "max_h1_without_cc": max_h1_without,
            "ratio": ratio,
            "n_h1_with_cc": len(ph_with["dgms"][1]) if PH_MAX_DIM >= 1 else 0,
            "n_h1_without_cc": len(ph_without["dgms"][1]) if PH_MAX_DIM >= 1 else 0,
        })

        for dim in range(PH_MAX_DIM + 1):
            np.save(ablation_dir / f"dgm{dim}_{tp}_nocc.npy", ph_without["dgms"][dim])

    # ── Permutation test on CC-ablated data ─────────────────
    logger.info("\n--- Permutation Tests (CC-ablated) ---")

    for tp in timepoints:
        if tp == "all":
            cells = adata_nocc
        else:
            cells = adata_nocc[adata_nocc.obs["timepoint"] == tp]

        max_h1 = [r for r in comparison if r["timepoint"] == tp][0]["max_h1_without_cc"]
        if max_h1 == 0:
            logger.info(f"  {tp}: no H₁ after ablation, skipping permutation test")
            continue

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

        for r in comparison:
            if r["timepoint"] == tp:
                r["p_value_nocc"] = perm["p_value"]

        fig = plot_permutation_test(
            perm["null_distribution"], perm["observed"], perm["p_value"]
        )
        save_figure(fig, f"fig4_permtest_nocc_{tp}")

    # ── Save comparison table ───────────────────────────────
    comp_df = pd.DataFrame(comparison)
    comp_df.to_csv(TABLES_DIR / "cell_cycle_ablation_comparison.csv", index=False)

    # ── Plot side-by-side barcodes ──────────────────────────
    fig_with = plot_persistence_barcodes(
        dgms_with_cc, max_dim=PH_MAX_DIM, title="With Cell Cycle Genes"
    )
    save_figure(fig_with, "fig4_a_barcodes_with_cc")

    fig_without = plot_persistence_barcodes(
        dgms_without_cc, max_dim=PH_MAX_DIM, title="Without Cell Cycle Genes"
    )
    save_figure(fig_without, "fig4_b_barcodes_without_cc")

    # ── Interpretation ──────────────────────────────────────
    logger.info("\n=== Cell Cycle Ablation Results ===")
    for r in comparison:
        tp = r["timepoint"]
        ratio = r["ratio"]

        if ratio > 0.7:
            verdict = "H₁ PERSISTS — loop is NOT cell cycle"
        elif ratio > 0.3:
            verdict = "H₁ PARTIALLY REDUCED — mixed cell cycle + other"
        elif r["max_h1_without_cc"] > 0:
            verdict = "H₁ MOSTLY REDUCED — likely cell cycle driven"
        else:
            verdict = "H₁ ELIMINATED — loop was cell cycle"

        logger.info(f"  {tp}: {verdict} (ratio={ratio:.2f})")

    logger.info(f"\nFull results: {TABLES_DIR / 'cell_cycle_ablation_comparison.csv'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cell cycle ablation experiment")
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    main(seed=args.seed)
